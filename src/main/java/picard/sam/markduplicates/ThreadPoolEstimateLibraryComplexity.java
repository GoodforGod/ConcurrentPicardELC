package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 10.02.2017
 * Time: 23:44
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.sam.DuplicationMetrics;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import static java.lang.Math.pow;

/*
 * DEFAULT COMMENT
 */
public class ThreadPoolEstimateLibraryComplexity extends ThreadedEstimateLibraryComplexity
{
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

        while (iterator.hasNext())
        {
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
                for (final Map.Entry<String, List<PairedReadSequence>> entry : sequencesByLibrary.entrySet())
                {
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
                ++groupsProcessed;

                if (lastLogTime < System.currentTimeMillis() - 60000) {
                    log.info("Processed " + groupsProcessed + " groups.");
                    lastLogTime = System.currentTimeMillis();
                }
            }
        }

        iterator.close();
        sorter.cleanup();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        for (final String library : duplicationHistosByLibrary.keySet())
        {
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
        }

        double elcTime;
        log.info("DUPLICATE - THREADED ELC (ms) : "
                + (elcTime = (double)((System.nanoTime() - startTime)/1000000)));
        log.info("TOTAL TIME - THREADED ELC (ms) : " + (sortTime + elcTime));

        file.write(OUTPUT);
        return 0;
    }
}
