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
import picard.sam.markduplicates.util.QueueIteratorProducer;

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
public class ConcurrentPoolEstimateLibraryComplexity extends ConcurrentExecutorEstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ConcurrentPoolEstimateLibraryComplexity.class);

    public ConcurrentPoolEstimateLibraryComplexity(){
        super();
    }

    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                || null != READ_ONE_BARCODE_TAG
                || null != READ_TWO_BARCODE_TAG);

        // Results from the doSort
        final ElcSmartSortResponse response = doSmartSort(useBarcodes);
        final long doWorkStartTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter     = response.getSorter();
        final ProgressLogger                                  progress   = response.getProgress();
        final List<SAMReadGroupRecord>                        readGroups = response.getReadGroup();

        final int meanGroupSize = (int)(Math.max(1, (progress.getCount() / 2) / (int)pow(4, MIN_IDENTICAL_BASES * 2)));

        final QueueIteratorProducer<PairedReadSequence, List<PairedReadSequence>> pairProducer
                = new QueueIteratorProducer<>(new PeekableIterator<>(sorter.iterator()), pairHandler);

        final ConcurrentSupplier<List<PairedReadSequence>> groupSupplier
                = new ConcurrentSupplier<>(GROUP_PROCESS_STACK_SIZE, USED_THREADS);

        final ConcurrentHistoCollection histoCollection = new ConcurrentHistoCollection(useBarcodes);

        //final ExecutorService pool = Executors.newFixedThreadPool(USED_THREADS);
        //final ForkJoinPool pool = new ForkJoinPool();
        final ExecutorService pool = Executors.newCachedThreadPool();
        final Object sync = new Object();

        // pool manager, receives job of groups, and make worker to process them
        // Now go through the sorted reads and attempt to find duplicates
        final long groupStartTime = System.nanoTime();
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                final List<List<PairedReadSequence>> groupList = groupSupplier.getJob();

                // Poison pill check
                if (isPairsPoisonPill.test(groupList))
                    return;

                pool.submit(() ->
                {
                    // Get the next group and split it apart by library
                    for (List<PairedReadSequence> group : groupList) {

                        if (group.size() > meanGroupSize * MAX_GROUP_RATIO)
                            logInvalidGroup(group, meanGroupSize);
                        else
                            splitByLibrary(group, readGroups).entrySet().forEach(histoCollection::processGroup);
                    }
                    groupSupplier.confirm();
                });
            }
        });

        // Iterating through sorted groups, and making job to process them
        while (pairProducer.hasNext()) {
            try                              { groupSupplier.add(pairProducer.next()); }
            catch (NoSuchElementException e) { e.printStackTrace(); }
        }

        groupSupplier.tryAddRest();
        groupSupplier.awaitConfirmation();
        groupSupplier.finish();

        final double groupEndTime = System.nanoTime();

        pairProducer.finish();
        sorter.cleanup();

        final long metricStartTime = System.nanoTime();
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

        double doWorkEndTime = System.nanoTime();
        double doWorkTotal = (doWorkEndTime - doWorkStartTime) / 1000000;
        log.info("GROUPS - FJ POOL (ms) : " + ((groupEndTime - groupStartTime) / 1000000));
        log.info("METRIC - FJ POOL (ms) : " + ((doWorkEndTime - metricStartTime) / 1000000));
        log.info("doWork - FJ POOL (ms) : " + doWorkTotal);
        log.info("----------------------------------------");
        log.info("TOTAL  - FJ POOL (ms) : " + (sortTime + doWorkTotal));

        file.write(OUTPUT);
        return 0;
    }
}
