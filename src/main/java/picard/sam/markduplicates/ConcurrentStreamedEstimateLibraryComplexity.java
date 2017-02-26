package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:29
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.sam.markduplicates.util.ConcurrentSortingCollection;
import picard.sam.markduplicates.util.QueueProducer;

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

        final long doWorkStartTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter = response.getSorter();
        final ProgressLogger progress                                = response.getProgress();
        final List<SAMReadGroupRecord> readGroups                    = response.getReadGroup();

        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));
        final Predicate<List<PairedReadSequence>> isGroupValid = (grp) -> grp.size() > meanGroupSize * MAX_GROUP_RATIO;

        final ConcurrentHistoCollection histoCollection = new ConcurrentHistoCollection(useBarcodes);
        final ConcurrentSupplier<List<PairedReadSequence>> groupSupplier
                = new ConcurrentSupplier<>(GROUP_PROCESS_STACK_SIZE, USED_THREADS);
        final QueueProducer<PairedReadSequence, List<PairedReadSequence>> pairProducer
                = new QueueProducer<>(new PeekableIterator<>(sorter.iterator()), pairHandler);

        final ForkJoinPool pool = new ForkJoinPool();
        final long groupStartTime = System.nanoTime();

        // Pool process sorted groups
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                final List<List<PairedReadSequence>> groupList = groupSupplier.getJob();

                if (isPairsPoisonPill.test(groupList))
                    return;

                pool.submit(() -> {
                    groupList.stream()
                            .parallel()
                            .unordered()
                            .flatMap(grp -> splitByLibrary(grp, readGroups).entrySet().stream())
                            .forEach(histoCollection::processGroup);

                    groupSupplier.confirm();
                });
            }
        });

        // Now go through the sorted reads and attempt to find duplicates
        while (pairProducer.hasNext()) {
            final List<PairedReadSequence> group = pairProducer.next();

            if (isGroupValid.test(group))
                logInvalidGroup(group, meanGroupSize);
            else
                groupSupplier.add(group);
        }

        groupSupplier.tryAddRest();
        groupSupplier.awaitConfirmation();
        groupSupplier.finish();

        double groupEndTime = System.nanoTime();

        pairProducer.finish();
        sorter.cleanup();

        final long metricStartTime = System.nanoTime();
        final ConcurrentMetrics concurrentMetrics = new ConcurrentMetrics(histoCollection);

        histoCollection.getLibraries().stream()
                .parallel()
                .unordered()
                .forEach(concurrentMetrics::add);

        concurrentMetrics.fillFile();

        pool.shutdown();
        try                             { pool.awaitTermination(1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        double doWorkEndTime = System.nanoTime();
        double doWorkTotal = (doWorkEndTime - doWorkStartTime) / 1000000;
        log.info("GROUPS - STREAMED (ms) : " + ((groupEndTime - groupStartTime) / 1000000));
        log.info("METRIC - STREAMED (ms) : " + ((doWorkEndTime - metricStartTime) / 1000000));
        log.info("doWork - STREAMED (ms) : " + doWorkTotal);
        log.info("----------------------------------------");
        log.info("TOTAL  - STREAMED (ms) : " + (sortTime + doWorkTotal));

        concurrentMetrics.writeFile();
        return 0;
    }
}
