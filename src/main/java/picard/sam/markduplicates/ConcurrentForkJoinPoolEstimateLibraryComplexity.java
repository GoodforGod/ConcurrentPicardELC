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
import java.util.concurrent.atomic.AtomicInteger;

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

        long startTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter     = response.getSorter();
        final ProgressLogger                                  progress   = response.getProgress();
        final List<SAMReadGroupRecord>                        readGroups = response.getReadGroup();

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<>(sorter.iterator());

        final Map<String, Histogram<Integer>> opticalHistosByLibrary     = new ConcurrentHashMap<>();
        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new ConcurrentHashMap<>();

        final AtomicInteger groupsProcessed = new AtomicInteger(0);
        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));

        final ElcDuplicatesFinderResolver algorithmResolver = new ElcDuplicatesFinderResolver(
                MAX_DIFF_RATE,
                MAX_READ_LENGTH,
                MIN_IDENTICAL_BASES,
                useBarcodes,
                opticalDuplicateFinder
        );

        final int streamReady = (int)(progress.getCount() / meanGroupSize) / 2;
        final int streamable = streamReady / 10000;

        final ForkJoinPool pool = new ForkJoinPool();

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);
        final Object sync = new Object();
        final BlockingQueue<List<List<PairedReadSequence>>> groupQueue = new LinkedBlockingQueue<>(2);

        List<List<PairedReadSequence>> groupStack = new ArrayList<>();

        final long startSortIterateTime = System.nanoTime();
        // pool manager, receives stack of groups, and make worker to process them
        pool.execute(() -> {
            while (!Thread.interrupted()) {
                try {
                    final List<List<PairedReadSequence>> groupList = groupQueue.take();
                    pool.submit(() ->
                    {
                        for (List<PairedReadSequence> group : groupList) {

                            // Get the next group and split it apart by library
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
                                groupsProcessed.incrementAndGet();
                            }
                        }
                        locker.decrementAndGet();
                    });
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
        });

        // Iterating through sorted groups, and making stack to process them
        while (iterator.hasNext()) {
            try {
                groupStack.add(getNextGroup(iterator));

                if(groupStack.size() >= GROUP_PROCESS_STACK_SIZE)
                {
                    try {
                        groupQueue.put(groupStack);
                        locker.incrementAndGet();
                        groupStack = new ArrayList<>(streamable);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            }
            catch (NoSuchElementException e) {
                e.printStackTrace();
            }
        }

        if (!groupStack.isEmpty()) {
            try {
                groupQueue.put(groupStack);
                locker.incrementAndGet();
                groupStack = new ArrayList<>();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // Waiting for all threads finish their job and check for missed group stacks
        while (locker.get() != 0) {
           if(locker.get() == 0)
               break;
        }

        log.info("SORTER PROCESSED FJPOOL - (ms) : " + (double)(System.nanoTime() - startSortIterateTime)/ 1000000);

        iterator.close();
        sorter.cleanup();

        long startMetricFile = System.nanoTime();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();

        for (final String library : duplicationHistosByLibrary.keySet())
        {
            locker.incrementAndGet();
            pool.execute(() ->
            {
                final Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                final Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);
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
                locker.decrementAndGet();
            });
        }

        // Waiting for all threads finish their job and check for missed group stacks
        while (locker.get() != 0) {
            if(locker.get() == 0)
                break;
        }

        pool.shutdown();
        try {
            pool.awaitTermination( 1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        double elcTime = System.nanoTime() / 1000000;
        log.info("METRIC - FJPOOL ELC (ms) : " + ((elcTime - (double)startMetricFile / 1000000)));
        log.info("TOTAL - FJPOOL ELC (ms) : " + (sortTime + (elcTime - (double)startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
