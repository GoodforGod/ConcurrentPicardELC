package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:28
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.ConcurrentSortingCollection;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.io.File;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;

import static java.lang.Math.nextUp;
import static java.lang.Math.pow;

/*
 * Thread Executor version of the ELC
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ConcurrentExecutorEstimateLibraryComplexity extends EstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ConcurrentExecutorEstimateLibraryComplexity.class);

    public ConcurrentExecutorEstimateLibraryComplexity() {
        super();
    }

    //  Temporary structure, used in sorting via StreamAPI (doStreamSort)
    private class PrsAndRec {
        public final SAMRecord rec;
        public final PairedReadSequence prs;
        //public final PairedReadSequenceWithBarcodes prsWithBarcodes;

        public PrsAndRec(PairedReadSequence prs, SAMRecord rec)
        {
            this.rec = rec;
            this.prs = prs;
            //this.prsWithBarcodes = prsWithBarcodes;
        }
    }

    // Temporary structure, used while processing histograms via StreamAPI
    protected class HistoAndMetric {
        final Histogram<Integer> duplicationHisto;
        final DuplicationMetrics metrics;

        public HistoAndMetric(DuplicationMetrics metrics,
                              Histogram<Integer> duplicationHisto)
        {
            this.duplicationHisto = duplicationHisto;
            this.metrics = metrics;
        }
    }

    // Takes temporary object (PrsAndRec), process it, fill reads and returns PairedReadSequence
    private PairedReadSequence FillPairedSequence(final PairedReadSequence prs, final SAMRecord record) {
/*        final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes)
                ? (PairedReadSequenceWithBarcodes) prs
                : null;*/

        // Get the bases and restore them to their original orientation if necessary
        final byte[] bases = record.getReadBases();

        if (record.getReadNegativeStrandFlag())
            SequenceUtil.reverseComplement(bases);

        if (record.getFirstOfPairFlag()) {
            prs.read1 = bases;

/*            if (useBarcodes) {
                prsWithBarcodes.barcode = getBarcodeValue(record);
                prsWithBarcodes.readOneBarcode = getReadOneBarcodeValue(record);
            }*/
        } else {
            prs.read2 = bases;

/*            if (useBarcodes)
                prsWithBarcodes.readTwoBarcode = getReadTwoBarcodeValue(record);*/
        }
        return prs;
    }

    // Replacement of the pendingByName (Concurrent)HashMap
    protected class ConcurrentPendingByNameCollection {
        private final int MAP_INIT_CAPACITY = 100;

        private final boolean useBarcodes;
        private final OpticalDuplicateFinder opticalDuplicateFinder;
        private final List<SAMReadGroupRecord> groupRecords;

        private final ConcurrentMap<String, PairedReadSequence> PairedRSMap = new ConcurrentHashMap<>(MAP_INIT_CAPACITY);

        public ConcurrentPendingByNameCollection(boolean useBarcodes,
                                                 OpticalDuplicateFinder opticalDuplicateFinder,
                                                 List<SAMReadGroupRecord> groupRecords) {
            this.useBarcodes = useBarcodes;
            this.opticalDuplicateFinder = opticalDuplicateFinder;
            this.groupRecords = groupRecords;
        }

        public PairedReadSequence remove(SAMRecord record) {
            PairedReadSequence prs = PairedRSMap.remove(record.getReadName());

            if (prs == null) {
                // Make a new paired read object and add RG and physical location information to it
                prs = useBarcodes
                        ? new PairedReadSequenceWithBarcodes()
                        : new PairedReadSequence();

                if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), prs)) {
                    final SAMReadGroupRecord rg = record.getReadGroup();
                    if (rg != null)
                        prs.setReadGroup((short) groupRecords.indexOf(rg));
                }
                PairedRSMap.put(record.getReadName(), prs);
            }
            return prs;
        }

        // REDUNDANT due to loosing reference on collection after cycle
        public void clear() {
            PairedRSMap.clear();
        }
    }

    //Check if the pair os OK to add to sorter
    private final Predicate<PairedReadSequence> isPairValid = (pair) -> (pair.read1 != null
                                                                      && pair.read2 != null
                                                                      && pair.qualityOk);
    // Filter is Record valid or not
    private final Predicate<SAMRecord> isRecValid = (rec) -> (rec.getReadPairedFlag()
                                                          && (rec.getFirstOfPairFlag() || rec.getSecondOfPairFlag())
                                                          && !rec.isSecondaryOrSupplementary());

    //
    protected final int RECORD_PROCESS_STACK_SIZE = MAX_RECORDS_IN_RAM / 1000;

    // Sorting via StreamAPI
    // IS NOT WORKING PROPERLY, some bugs in parallel process pairs, when fill bases or wahtever
    protected ElcSmartSortResponse doStreamSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        final long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        final ConcurrentSortingCollection<PairedReadSequence>
                sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                    !useBarcodes
                            ? new PairedReadCodec()
                            : new PairedReadWithBarcodesCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);

        // Check for valid pairs & records
        final Predicate<PrsAndRec> isPrsAndRecValid = (prsAndRec) -> (prsAndRec.prs != null
                                                                    && (prsAndRec.prs.qualityOk = prsAndRec.prs.qualityOk
                                                                    && passesQualityCheck(prsAndRec.rec.getReadBases(),
                                                                                          prsAndRec.rec.getBaseQualities(),
                                                                                          MIN_IDENTICAL_BASES,
                                                                                          MIN_MEAN_QUALITY)));
        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());
            //final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final ConcurrentPendingByNameCollection pendingByName
                    = new ConcurrentPendingByNameCollection(useBarcodes, opticalDuplicateFinder, readGroups);

            in.iterator().stream()
                    .parallel()
                    .unordered()
                    .filter(isRecValid)
                    .map(rec -> {
                        progress.record(rec);
                        PairedReadSequence prs = pendingByName.remove(rec);
                        // prs.qualityOk check - partly skipped!!!
                        if(isPrsAndRecValid.test(new PrsAndRec(prs, rec)))
                            return FillPairedSequence(prs, rec);
                        else
                            return new PairedReadSequence();
                    })
                    //.filter(isPrsAndRecValid)
                    //.map(pairAndRec -> FillPairedSequence(pairAndRec.prs, pairAndRec.rec))
                    .filter(isPairValid)
                    .forEachOrdered(sorter::add);

            CloserUtil.close(in);
        }

        log.info("FILL SORTING - STREAM-SORT ELC (ms) : "
                + (sortTime = (double) ((System.nanoTime() - startTime) / 1000000)));
        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates." ,
                progress.getCount()));

        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Sorting with ExecutorService
    protected ElcSmartSortResponse doSmartSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        final ConcurrentSortingCollection<PairedReadSequence>
                sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                    !useBarcodes
                            ? new PairedReadCodec()
                            : new PairedReadWithBarcodesCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);

        final int recordListSize = RECORD_PROCESS_STACK_SIZE;
        final ExecutorService pool = Executors.newFixedThreadPool(2);
        final BlockingQueue<List<SAMRecord>> queueRecords = new ArrayBlockingQueue<>(10);
        List<SAMRecord> records = new ArrayList<>(recordListSize);

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);

        //final Object sync = new Object();

        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());
            //final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final ConcurrentPendingByNameCollection pendingByName
                    = new ConcurrentPendingByNameCollection(useBarcodes, opticalDuplicateFinder, readGroups);

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    try {
                        final List<SAMRecord> rec = queueRecords.take();
                        pool.submit(() -> {
                            for (SAMRecord record : rec) {
                                PairedReadSequence prs = pendingByName.remove(record);

                                // Read passes quality check if both ends meet the mean quality criteria
                                prs.qualityOk = prs.qualityOk && passesQualityCheck(record.getReadBases(),
                                        record.getBaseQualities(),
                                        MIN_IDENTICAL_BASES,
                                        MIN_MEAN_QUALITY);

                                final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes)
                                        ? (PairedReadSequenceWithBarcodes) prs
                                        : null;

                                // Get the bases and restore them to their original orientation if necessary
                                final byte[] bases = record.getReadBases();

                                if (record.getReadNegativeStrandFlag())
                                    SequenceUtil.reverseComplement(bases);

                                if (record.getFirstOfPairFlag()) {
                                    prs.read1 = bases;

                                    if (useBarcodes) {
                                        prsWithBarcodes.barcode = getBarcodeValue(record);
                                        prsWithBarcodes.readOneBarcode = getReadOneBarcodeValue(record);
                                    }
                                } else {
                                    prs.read2 = bases;

                                    if (useBarcodes)
                                        prsWithBarcodes.readTwoBarcode = getReadTwoBarcodeValue(record);
                                }

                                progress.record(record);
                                if (isPairValid.test(prs))
                                    sorter.add(prs);
                            }
                            locker.decrementAndGet();
                        });
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            });

            for (SAMRecord rec : in) {
                if (!isRecValid.test(rec))
                    continue;

                records.add(rec);

                if (records.size() >= recordListSize) {
                    try {
                        queueRecords.put(records);
                        locker.incrementAndGet();
                        records = new ArrayList<>(recordListSize);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            }

            if (!records.isEmpty()) {
                try {
                    queueRecords.put(records);
                    locker.incrementAndGet();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            // Some useless work to wait for locks & sorter to spill
            while (locker.get() != 0 || sorter.isSpillingToDisk()){
                if(locker.get() == 0 && !sorter.isSpillingToDisk())
                    break;
            }

            CloserUtil.close(in);
        }

        pool.shutdown();
        try {
            pool.awaitTermination(1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        log.info("SORTING - SMART-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime) / 1000000)));

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));
        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Sorting with ForkJoinPool difference is in ForkJoinPool instead of doSmartSort
    protected ElcSmartSortResponse doPoolSort(boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        final ConcurrentSortingCollection<PairedReadSequence>
                sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                    !useBarcodes
                            ? new PairedReadCodec()
                            : new PairedReadWithBarcodesCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);

        final int recordListSize = RECORD_PROCESS_STACK_SIZE;
        final ForkJoinPool pool = new ForkJoinPool();
        final BlockingQueue<List<SAMRecord>> queueRecords = new ArrayBlockingQueue<>(10);
        List<SAMRecord> records = new ArrayList<>(recordListSize);

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);

        //final Object sync = new Object();

        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());
            //final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final ConcurrentPendingByNameCollection pendingByName
                    = new ConcurrentPendingByNameCollection(useBarcodes, opticalDuplicateFinder, readGroups);

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    try {
                        final List<SAMRecord> rec = queueRecords.take();
                        pool.submit(() -> {
                            for (SAMRecord record : rec) {
                                PairedReadSequence prs = pendingByName.remove(record);

                                // Read passes quality check if both ends meet the mean quality criteria
                                prs.qualityOk = prs.qualityOk && passesQualityCheck(record.getReadBases(),
                                        record.getBaseQualities(),
                                        MIN_IDENTICAL_BASES,
                                        MIN_MEAN_QUALITY);

                                final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes)
                                        ? (PairedReadSequenceWithBarcodes) prs
                                        : null;

                                // Get the bases and restore them to their original orientation if necessary
                                final byte[] bases = record.getReadBases();

                                if (record.getReadNegativeStrandFlag())
                                    SequenceUtil.reverseComplement(bases);

                                if (record.getFirstOfPairFlag()) {
                                    prs.read1 = bases;

                                    if (useBarcodes) {
                                        prsWithBarcodes.barcode = getBarcodeValue(record);
                                        prsWithBarcodes.readOneBarcode = getReadOneBarcodeValue(record);
                                    }
                                } else {
                                    prs.read2 = bases;

                                    if (useBarcodes)
                                        prsWithBarcodes.readTwoBarcode = getReadTwoBarcodeValue(record);
                                }

                                progress.record(record);
                                if (isPairValid.test(prs))
                                    sorter.add(prs);
                            }
                            locker.decrementAndGet();
                        });
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            });

            for (SAMRecord rec : in) {
                if (!isRecValid.test(rec))
                    continue;

                records.add(rec);

                if (records.size() >= recordListSize) {
                    try {
                        queueRecords.put(records);
                        locker.incrementAndGet();
                        records = new ArrayList<>(recordListSize);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            }

            if (!records.isEmpty()) {
                try {
                    queueRecords.put(records);
                    locker.incrementAndGet();
                    records = new ArrayList<>();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            // Some useless work to wait for locks
            while (locker.get() != 0){
                if(locker.get() == 0)
                    break;
            }

            CloserUtil.close(in);
        }

        pool.shutdown();
        try {
            pool.awaitTermination(1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        log.info("SORTING - POOL-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime)/1000000)));

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));
        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Group size to process stack in doWork
    protected final int GROUP_PROCESS_STACK_SIZE = 50000;

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

        final ExecutorService pool = Executors.newCachedThreadPool();

        final Object sync = new Object();

        List<List<PairedReadSequence>> groupStack = new ArrayList<>();
        final BlockingQueue<List<List<PairedReadSequence>>> groupQueue = new LinkedBlockingQueue<>(2);

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);

        final long startSortIterateTime = System.nanoTime();

        // pool manager, receives stack of groups, and make worker to process them
        // Somehow we process less groups then in sequential doWork 241000 vs 715000,
        // concurrentSorter has bugs, but where... who knows
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
                                // Now process the reads by library
                                for (final Map.Entry<String, List<PairedReadSequence>> entry : splitByLibrary(group, readGroups).entrySet()) {
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

                if (groupStack.size() >= GROUP_PROCESS_STACK_SIZE) {
                    groupQueue.put(groupStack);
                    locker.incrementAndGet();
                    groupStack = new ArrayList<>(streamable);
                }
            } catch (NoSuchElementException | InterruptedException e) {
                e.printStackTrace();
            }
        }

        // Waiting for all threads finish their job and check for missed group stacks
        while (locker.get() != 0) {
            if (!groupStack.isEmpty()) {
                try {
                    groupQueue.put(groupStack);
                    locker.incrementAndGet();
                    groupStack = new ArrayList<>();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

        // Shutting pool down after all work is finished
      /*  pool.shutdown();
        try {
            pool.awaitTermination(1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }*/

        log.info("SORTER - EXECUTOR (ms) : " + (double)(System.nanoTime() - startSortIterateTime) / 1000000);

        iterator.close();
        sorter.cleanup();

        long startMetricFile = System.nanoTime();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();
        //final ExecutorService histoPool = Executors.newCachedThreadPool();

        for (final String library : duplicationHistosByLibrary.keySet())
        {
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
            });
        }

        pool.shutdown();
        try {
            pool.awaitTermination( 1, TimeUnit.MICROSECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        double elcTime = System.nanoTime() / 1000000;
        log.info("METRIC - EXECUTOR ELC (ms) : " + ((elcTime - (double)startMetricFile / 1000000)));
        log.info("TOTAL - EXECUTOR ELC (ms) : " + (sortTime + (elcTime - (double)startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
