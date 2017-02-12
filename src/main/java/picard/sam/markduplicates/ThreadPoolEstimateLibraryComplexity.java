package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 10.02.2017
 * Time: 23:44
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static java.lang.Math.pow;

/*
 * DEFAULT COMMENT
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ThreadPoolEstimateLibraryComplexity extends ThreadedEstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ThreadPoolEstimateLibraryComplexity.class);

    public ThreadPoolEstimateLibraryComplexity(){
        super();
    }

    protected int doWork() {

        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                || null != READ_ONE_BARCODE_TAG
                || null != READ_TWO_BARCODE_TAG);

        final ELCSortResponse response = doSmartSort(useBarcodes);

        long startTime = System.nanoTime();

        final SortingCollection<PairedReadSequence> sorter = response.getSorter();
        final ProgressLogger progress = response.getProgress();
        final List<SAMReadGroupRecord> readGroups = response.getReadGroup();

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<PairedReadSequence>(sorter.iterator());

        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new HashMap<String, Histogram<Integer>>();
        final Map<String, Histogram<Integer>> opticalHistosByLibrary = new HashMap<String, Histogram<Integer>>();

        int groupsProcessed = 0;
        long lastLogTime = System.currentTimeMillis();
        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));

        ElcDuplicatesFinderResolver algorithmResolver = new ElcDuplicatesFinderResolver(
                MAX_DIFF_RATE,
                MAX_READ_LENGTH,
                MIN_IDENTICAL_BASES,
                useBarcodes,
                opticalDuplicateFinder
        );

        final int CORE_COUNT = 4 * 2;
        BlockingQueue iterateQueue = new ArrayBlockingQueue(CORE_COUNT);

        long startSortIterateTime = System.nanoTime();
        int streamReady = (int)(progress.getCount() / meanGroupSize) / 2;
        int streamable = streamReady / 10000;

        ForkJoinPool pool = new ForkJoinPool();

        final List<List<PairedReadSequence>> temporaryGroups = new ArrayList<List<PairedReadSequence>>();
        ConcurrentLinkedQueue<List<List<PairedReadSequence>>> groupQueue
                = new ConcurrentLinkedQueue<List<List<PairedReadSequence>>>();

        pool.execute(() -> {
            final List<List<PairedReadSequence>> groupList = groupQueue.poll();
            pool.submit(() ->
            {
                for(List<PairedReadSequence> groupPairs : groupList) {
                    long startWork = System.nanoTime();
                    // Get the next group and split it apart by library
                    final List<PairedReadSequence> group = getNextGroup(iterator);

                    if (group.size() > meanGroupSize * MAX_GROUP_RATIO) {
                        final PairedReadSequence prs = group.get(0);
                        log.warn("Omitting group with over " + MAX_GROUP_RATIO + " times the expected mean number of read pairs. " +
                                "Mean=" + meanGroupSize + ", Actual=" + group.size() + ". Prefixes: " +
                                StringUtil.bytesToString(prs.read1, 0, MIN_IDENTICAL_BASES) +
                                " / " +
                                StringUtil.bytesToString(prs.read2, 0, MIN_IDENTICAL_BASES));
                    } else {
                        final Map<String, List<PairedReadSequence>> sequencesByLibrary = splitByLibrary(group, readGroups);

                        // Now process the reads by library
                        for (final Map.Entry<String, List<PairedReadSequence>> entry : sequencesByLibrary.entrySet()) {
                            final String library = entry.getKey();
                            final List<PairedReadSequence> seqs = entry.getValue();

                            Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                            Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);

                            if (duplicationHisto == null) {
                                duplicationHisto = new Histogram<>("duplication_group_count", library);
                                opticalHisto = new Histogram<>("duplication_group_count", "optical_duplicates");
                                duplicationHistosByLibrary.put(library, duplicationHisto);
                                opticalHistosByLibrary.put(library, opticalHisto);
                            }
                            algorithmResolver.resolveAndSearch(seqs, duplicationHisto, opticalHisto);
                        }
                    }
                }
            });
        });

        while (iterator.hasNext()) {
            final List<PairedReadSequence> group = getNextGroup(iterator);

            if (group.size() > meanGroupSize * MAX_GROUP_RATIO)
            {
                final PairedReadSequence prs = group.get(0);
                log.warn("Omitting group with over " + MAX_GROUP_RATIO
                        + " times the expected mean number of read pairs. Mean=" + meanGroupSize
                        + ", Actual="    + group.size()
                        + ". Prefixes: " + StringUtil.bytesToString(prs.read1, 0, MIN_IDENTICAL_BASES)
                        + " / "          + StringUtil.bytesToString(prs.read2, 0, MIN_IDENTICAL_BASES));
            }
            else
            {
                groupsProcessed++;
                temporaryGroups.add(group);

                if(temporaryGroups.size() >= streamable || !iterator.hasNext())
                {
                    groupQueue.add(new ArrayList<>(temporaryGroups));
                    temporaryGroups.clear();
                }
            }
        }

        log.info("SORTER PROCESS - (ms) : "
                + (double)(System.nanoTime() - startSortIterateTime)/ 1000000);

        iterator.close();
        sorter.cleanup();

        long startMetricFile = System.nanoTime();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        for (final String library : duplicationHistosByLibrary.keySet())
        {
            long startHistogram = System.nanoTime();

            final Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
            final Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);
            final DuplicationMetrics metrics = new DuplicationMetrics();
            metrics.LIBRARY = library;

            // Filter out any bins that have fewer than MIN_GROUP_COUNT entries in them and calculate derived metrics
            for (final Integer bin : duplicationHisto.keySet())
            {
                final double duplicateGroups = duplicationHisto.get(bin).getValue();
                final double opticalDuplicates = opticalHisto.get(bin) == null ? 0 : opticalHisto.get(bin).getValue();

                if (duplicateGroups >= MIN_GROUP_COUNT)
                {
                    metrics.READ_PAIRS_EXAMINED += (bin * duplicateGroups);
                    metrics.READ_PAIR_DUPLICATES += ((bin - 1) * duplicateGroups);
                    metrics.READ_PAIR_OPTICAL_DUPLICATES += opticalDuplicates;
                }
            }
            metrics.calculateDerivedFields();
            file.addMetric(metrics);
            file.addHistogram(duplicationHisto);
/*
            log.info("HISTOGRAM INFO - "
                    + (double)(System.nanoTime() - startHistogram)
                    + " LIBRARY : "
                    + library
                    + " | duplicationHisto SIZE : "
                    + duplicationHisto.size()
                    + " | opticalHisto SIZE : "
                    + opticalHisto.size());
*/
        }

        double elcTime;
        log.info("METRIC - THREAD POOL ELC (ms) : "
                + ((double)((System.nanoTime() - startMetricFile)/1000000)));

        log.info("doWork - THREAD POOL ELC (ms) : "
                + (elcTime = (double)((System.nanoTime() - startTime)/1000000)));
        log.info("TOTAL TIME - THREAD POOL ELC (ms) : " + (sortTime + elcTime));

        file.write(OUTPUT);
        return 0;
    }
}
