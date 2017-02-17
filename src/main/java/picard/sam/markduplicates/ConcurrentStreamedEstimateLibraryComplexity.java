package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:29
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.ConcurrentSortingCollection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static java.lang.Math.decrementExact;
import static java.lang.Math.pow;

/*
 * Streamed version of ELC
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ConcurrentStreamedEstimateLibraryComplexity extends ConcurrentExecutorEstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ConcurrentStreamedEstimateLibraryComplexity.class);

    public ConcurrentStreamedEstimateLibraryComplexity() {
        super();
    }

    public void Initiate(String[] args) {
        instanceMain(args);
    }

    protected int doWork() {

        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                                    || null != READ_ONE_BARCODE_TAG
                                    || null != READ_TWO_BARCODE_TAG);

        final ElcSmartSortResponse response = doSmartSort(useBarcodes);

        long startTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter = response.getSorter();
        final ProgressLogger progress = response.getProgress();
        final List<SAMReadGroupRecord> readGroups = response.getReadGroup();

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<>(sorter.iterator());

        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new ConcurrentHashMap<>();
        final Map<String, Histogram<Integer>> opticalHistosByLibrary = new ConcurrentHashMap<>();

        int groupsProcessed = 0;

        //long lastLogTime = System.currentTimeMillis();
        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));

        ElcDuplicatesFinderResolver algorithmResolver = new ElcDuplicatesFinderResolver(
                MAX_DIFF_RATE,
                MAX_READ_LENGTH,
                MIN_IDENTICAL_BASES,
                useBarcodes,
                opticalDuplicateFinder
        );

        //Predicate<List<PairedReadSequence>> pairFilter = (grp) -> grp.size() > meanGroupSize * MAX_GROUP_RATIO;

        int streamReady = (int)(progress.getCount() / meanGroupSize) / 2;
        int streamable = streamReady / 1000;

        final ForkJoinPool pool = new ForkJoinPool();
        final BlockingQueue<List<List<PairedReadSequence>>> groupQueue = new LinkedBlockingQueue<>();
        final AtomicInteger locker = new AtomicInteger(0);

        List<List<PairedReadSequence>> temporaryGroups = new ArrayList<>(streamable);

        long startSortIterateTime = System.nanoTime();
        // Pool process sorted groups
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                try {
                    final List<List<PairedReadSequence>> groupList = groupQueue.take();
                    pool.submit(() -> {
                        groupList.stream()
                                .parallel()
                                .unordered()
                                .map(grps -> splitByLibrary(grps, readGroups))
                                .flatMap(entry -> entry.entrySet().stream())
                                .forEach(entry -> {
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
                                });
                        locker.decrementAndGet();
                    });
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
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
                ++groupsProcessed;
                temporaryGroups.add(group);
                if(temporaryGroups.size() >= GROUP_PROCESS_STACK_SIZE)
                {
                    try {
                        groupQueue.put(temporaryGroups);
                        locker.incrementAndGet();
                        temporaryGroups = new ArrayList<>(streamable);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            }
        }

        if (!temporaryGroups.isEmpty()) {
            try {
                groupQueue.put(temporaryGroups);
                locker.incrementAndGet();
                temporaryGroups = new ArrayList<>();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        while (locker.get() != 0) {
          if(locker.get() == 0)
              break;
        }

        log.info("SORTER - STREAMED (ms) : "
                + (double)(System.nanoTime() - startSortIterateTime)/ 1000000);

        iterator.close();
        sorter.cleanup();

        final long startMetricFile = System.nanoTime();
        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        final Object sync = new Object();

        for (final String library : duplicationHistosByLibrary.keySet()) {
            locker.incrementAndGet();
            pool.execute(() ->
            {
                final Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                final Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);
                final DuplicationMetrics metrics = new DuplicationMetrics();

                metrics.LIBRARY = library;
                // Filter out any bins that have fewer than MIN_GROUP_COUNT entries in them and calculate derived metrics

                for (final Integer bin : duplicationHisto.keySet()) {
                    final double duplicateGroups = duplicationHisto.get(bin).getValue();
                    final double opticalDuplicates = (opticalHisto.get(bin) == null)
                                                   ? 0
                                                   : opticalHisto.get(bin).getValue();

                    if (duplicateGroups >= MIN_GROUP_COUNT) {
                        metrics.READ_PAIRS_EXAMINED += (bin * duplicateGroups);
                        metrics.READ_PAIR_DUPLICATES += ((bin - 1) * duplicateGroups);
                        metrics.READ_PAIR_OPTICAL_DUPLICATES += opticalDuplicates;
                    }
                }
                metrics.calculateDerivedFields();
                synchronized (sync) {
                    file.addMetric(metrics);
                    file.addHistogram(duplicationHisto);
                }
                locker.decrementAndGet();
            });
        }

        while (locker.get() != 0) {
          if(locker.get() == 0)
              break;
        }

        pool.shutdown();
        try {
            pool.awaitTermination(1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        long elcTime = System.nanoTime() / 1000000;
        log.info("METRIC - THREAD POOL ELC (ms) : " + ((double)(elcTime - startMetricFile/1000000)));
        log.info("TOTAL - THREAD POOL ELC (ms) : " + (sortTime + (elcTime - startTime/1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
