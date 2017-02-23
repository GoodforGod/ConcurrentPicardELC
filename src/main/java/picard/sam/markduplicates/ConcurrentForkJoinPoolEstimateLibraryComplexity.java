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
import picard.sam.markduplicates.util.ConcurrentSortingCollection;

import java.util.*;
import java.util.concurrent.*;

import static java.lang.Math.pow;

/*
 * FJPOOL Implementation of ELC
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ConcurrentForkJoinPoolEstimateLibraryComplexity extends ConcurrentExecutorEstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ConcurrentForkJoinPoolEstimateLibraryComplexity.class);

    public ConcurrentForkJoinPoolEstimateLibraryComplexity(){
        super();
    }

    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                || null != READ_ONE_BARCODE_TAG
                || null != READ_TWO_BARCODE_TAG);

        // Results from the doSort
        final ElcSmartSortResponse response = doSmartSort(useBarcodes);
        final long startTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter     = response.getSorter();
        final ProgressLogger                                  progress   = response.getProgress();
        final List<SAMReadGroupRecord>                        readGroups = response.getReadGroup();

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<>(sorter.iterator());

        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));
        final ConcurrentHistoCollection histoCollection = new ConcurrentHistoCollection(useBarcodes);
        final ConcurrentSupplier<List<PairedReadSequence>> groupSupplier
                = new ConcurrentSupplier<>(GROUP_PROCESS_STACK_SIZE, USED_THREADS);

        final ForkJoinPool pool = new ForkJoinPool();
        final Object sync = new Object();

        // pool manager, receives stack of groups, and make worker to process them
        final long startSortIterateTime = System.nanoTime();
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                final List<List<PairedReadSequence>> groupList = groupSupplier.getJob();

                // Poison pill check
                if (groupList.isEmpty())
                    return;

                pool.submit(() ->
                {
                    // Get the next group and split it apart by library
                    for (List<PairedReadSequence> group : groupList) {
                        if (group.size() > meanGroupSize * MAX_GROUP_RATIO)
                            logInvalidGroup(group, meanGroupSize);
                        else
                            splitByLibrary(group, readGroups).forEach(histoCollection::processGroup);
                    }
                    groupSupplier.confirm();
                });
            }
        });

        // Iterating through sorted groups, and making stack to process them
        while (iterator.hasNext()) {
            try                              { groupSupplier.add(getNextGroup(iterator)); }
            catch (NoSuchElementException e) { e.printStackTrace(); }
        }

        groupSupplier.tryAddRest();
        groupSupplier.awaitConfirmation();
        groupSupplier.finish();

        log.info("SORTER - EXECUTOR (ms) : " + (double)(System.nanoTime() - startSortIterateTime) / 1000000);

        iterator.close();
        sorter.cleanup();

        final long startMetricFile = System.nanoTime();
        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        for (final String library : histoCollection.getLibraries())
        {
            pool.execute(() ->
            {
                final Histogram<Integer> duplicationHisto = histoCollection.getDuplicationHisto(library);
                final Histogram<Integer> opticalHisto = histoCollection.getOpticalHisto(library);
                final DuplicationMetrics metrics = new DuplicationMetrics();

                metrics.LIBRARY = library;

                // Filter out any bins that have fewer than MIN_GROUP_COUNT entries in them and calculate derived metrics
                for (final Integer bin : duplicationHisto.keySet())
                {
                    final double duplicateGroups = duplicationHisto.get(bin).getValue();
                    final double opticalDuplicates = opticalHisto.get(bin) == null
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
        try                             { pool.awaitTermination( 1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        double elcTime = System.nanoTime() / 1000000;
        log.info("METRIC - EXECUTOR ELC (ms) : " + (elcTime - startMetricFile / 1000000));
        log.info("TOTAL - EXECUTOR ELC (ms) : " + (sortTime + (elcTime - startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
