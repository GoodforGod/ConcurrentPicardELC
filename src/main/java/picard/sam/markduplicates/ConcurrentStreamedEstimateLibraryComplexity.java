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

import java.util.List;
import java.util.concurrent.*;
import java.util.function.Predicate;

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

    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                                    || null != READ_ONE_BARCODE_TAG
                                    || null != READ_TWO_BARCODE_TAG);

        final ElcSmartSortResponse response = doSmartSort(useBarcodes);

        final long startTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter = response.getSorter();
        final ProgressLogger progress                                = response.getProgress();
        final List<SAMReadGroupRecord> readGroups                    = response.getReadGroup();

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<>(sorter.iterator());

        //long lastLogTime = System.currentTimeMillis();
        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));
        final Predicate<List<PairedReadSequence>> isGroupValid = (grp) -> grp.size() > meanGroupSize * MAX_GROUP_RATIO;

        final ConcurrentHistoCollection histoCollection = new ConcurrentHistoCollection(useBarcodes);
        final ConcurrentSupplier<List<PairedReadSequence>> groupSupplier
                = new ConcurrentSupplier<>(GROUP_PROCESS_STACK_SIZE, USED_THREADS);

        final ForkJoinPool pool = new ForkJoinPool();
        final Object sync = new Object();

        long startSortIterateTime = System.nanoTime();
        // Pool process sorted groups
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                final List<List<PairedReadSequence>> groupList = groupSupplier.getJob();

                if (groupList.isEmpty())
                    return;

                pool.submit(() ->
                        groupList.stream()
                                .parallel()
                                .unordered()
                                .map(grps -> splitByLibrary(grps, readGroups))
                                .flatMap(entry -> entry.entrySet().stream())
                                .forEach(histoCollection::processGroup));
            }
        });

        while (iterator.hasNext()) {
            final List<PairedReadSequence> group = getNextGroup(iterator);

            if (isGroupValid.test(group))
                logInvalidGroup(group, meanGroupSize);
            else
                groupSupplier.add(group);
        }

        groupSupplier.tryAddRest();
        groupSupplier.awaitConfirmation();
        groupSupplier.finish();

        log.info("SORTER - STREAMED (ms) : "
                + (double)(System.nanoTime() - startSortIterateTime)/ 1000000);

        iterator.close();
        sorter.cleanup();

        final long startMetricFile = System.nanoTime();
        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        for (final String library : histoCollection.getLibraries()) {
            pool.execute(() ->
            {
                final Histogram<Integer> duplicationHisto = histoCollection.getDuplicationHisto(library);
                final Histogram<Integer> opticalHisto = histoCollection.getOpticalHisto(library);
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
            });
        }

        pool.shutdown();
        try                             { pool.awaitTermination(1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        double elcTime = System.nanoTime() / 1000000;
        log.info("METRIC - THREAD POOL ELC (ms) : " + (elcTime - startMetricFile / 1000000));
        log.info("TOTAL - THREAD POOL ELC (ms) : " + (sortTime + (elcTime - startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
