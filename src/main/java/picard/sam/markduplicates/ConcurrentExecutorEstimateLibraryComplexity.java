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
import picard.sam.markduplicates.util.QueueIteratorProducer;
import picard.sam.markduplicates.util.QueuePeekableIterator;

import java.io.File;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
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

    /**
     * Classes
     */
    // Temporary structure, used in sorting via StreamAPI (doStreamSort)
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

    // Temporary structure, used while processing histograms via StreamAPI
    protected class HistoTuple {
        public final Histogram<Integer> duplicate;
        public final Histogram<Integer> optical;

        public HistoTuple(Histogram<Integer> duplicate, Histogram<Integer> optical) {
            this.duplicate = duplicate;
            this.optical = optical;
        }
    }

    // Replacement of the pendingByName (Concurrent)HashMap
    protected class ConcurrentPendingByNameCollection {
        private final int MAP_INIT_CAPACITY = 100;

        private final boolean useBarcodes;
        private final OpticalDuplicateFinder opticalDuplicateFinder;
        private final List<SAMReadGroupRecord> groupRecords;

        private final Map<String, PairedReadSequence> PairedRSMap = new ConcurrentHashMap<>(MAP_INIT_CAPACITY);

        public ConcurrentPendingByNameCollection(boolean useBarcodes,
                                                 OpticalDuplicateFinder opticalDuplicateFinder,
                                                 List<SAMReadGroupRecord> groupRecords) {
            this.useBarcodes = useBarcodes;
            this.opticalDuplicateFinder = opticalDuplicateFinder;
            this.groupRecords = groupRecords;
        }

        public synchronized PairedReadSequence getAsync(SAMRecord record) {
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

            prs.qualityOk = prs.qualityOk && passesQualityCheck(record.getReadBases(),
                    record.getBaseQualities(),
                    MIN_IDENTICAL_BASES,
                    MIN_MEAN_QUALITY);

            return prs;
        }

        public PairedReadSequence get(SAMRecord record) {
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

            prs.qualityOk = prs.qualityOk && passesQualityCheck(record.getReadBases(),
                    record.getBaseQualities(),
                    MIN_IDENTICAL_BASES,
                    MIN_MEAN_QUALITY);

            return prs;
        }

        // REDUNDANT due to loosing reference on collection after cycle
        public void clear() {
            PairedRSMap.clear();
        }
    }

    //
    protected class ConcurrentHistoCollection{
        private final ElcDuplicatesFinderResolver algorithmResolver;

        private final Map<String, Histogram<Integer>> opticalHistosByLibrary     = new ConcurrentHashMap<>();
        private final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new ConcurrentHashMap<>();

        public ConcurrentHistoCollection(boolean useBarcodes) {
            this.algorithmResolver = new ElcDuplicatesFinderResolver(
                    MAX_DIFF_RATE,
                    MAX_READ_LENGTH,
                    MIN_IDENTICAL_BASES,
                    useBarcodes,
                    opticalDuplicateFinder
            );
        }

        public Set<String> getLibraries() {
            return duplicationHistosByLibrary.keySet();
        }

        public Histogram<Integer> getDuplicationHisto(String library) {
            return duplicationHistosByLibrary.get(library);
        }

        public Histogram<Integer> getOpticalHisto(String library) {
            return opticalHistosByLibrary.get(library);
        }

        public void processGroup(final String library, final List<PairedReadSequence> seqs) {
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

        public void processGroup(final Map.Entry<String, List<PairedReadSequence>> entry) {
            processGroup(entry.getKey(), entry.getValue());
        }
    }

    //
    protected class ConcurrentSupplier<T> {
        private final int jobCapacity;
        private final int DEFAULT_JOB_QUEUE_CAPACITY = 2;

        private final List<T> poisonPill = Collections.emptyList();
        private final AtomicInteger jobsInProgress = new AtomicInteger(0);
        private final BlockingQueue<List<T>> jobQueue;

        private long counter = 0;
        private List<T> job;

        public ConcurrentSupplier(int jobCapacity) {
            this.jobCapacity = jobCapacity;
            this.jobQueue = new LinkedBlockingDeque<>(DEFAULT_JOB_QUEUE_CAPACITY);
            this.job = new ArrayList<>();
        }

        public ConcurrentSupplier(int jobCapacity, int queueCapacity) {
            this.jobCapacity = jobCapacity;
            this.jobQueue = new LinkedBlockingDeque<>(queueCapacity);
            this.job = new ArrayList<>(jobCapacity);
        }

        public void add(T item) {
            job.add(item);

            counter++;

            if(job.size() == jobCapacity) {
                putToQueue();
                job = new ArrayList<>(jobCapacity);
            }
        }

        private void putToQueue() {
            try                             {
                jobQueue.put(job);
                jobsInProgress.incrementAndGet();
            }
            catch (InterruptedException e)  { e.printStackTrace(); }
        }

        //
        public void tryAddRest() {
            if(!job.isEmpty()) {
                putToQueue();
                job = new ArrayList<T>();
            }
        }

        public List<T> getJob() {
            try                             { return jobQueue.take(); }
            catch (InterruptedException e)  { e.printStackTrace(); }
            return poisonPill;
        }

        public void confirm() {
            jobsInProgress.decrementAndGet();
        }

        // Waiting for all threads finish their job and check for missed group stacks
        public void awaitConfirmation() {
            while (jobsInProgress.get() != 0) {
                if(jobsInProgress.get() == 0)
                    break;
            }
        }

        public void finish() {
            try                             { jobQueue.put(poisonPill); }
            catch (InterruptedException e)  { e.printStackTrace(); }
        }
    }

    //
    protected class ConcurrentMetrics {

        private Map<Histogram<Integer>, DuplicationMetrics> duplicationMetrics = new ConcurrentHashMap<>();

        public void add(final String library,
                        final Histogram<Integer> duplicationHisto,
                        final Histogram<Integer> opticalHisto) {
            final DuplicationMetrics metrics = new DuplicationMetrics();

            metrics.LIBRARY = library;
            // Filter out any bins that have fewer than MIN_GROUP_COUNT entries in them and calculate derived metrics

            //duplicationHisto.keySet().stream()

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
            duplicationMetrics.put(duplicationHisto, metrics);
        }

        public Set<Map.Entry<Histogram<Integer>, DuplicationMetrics>> getMetrics() {
            return duplicationMetrics.entrySet();
        }
    }

    /**
     * Predicates
     */
    // Check for poison pill into doWork method
    protected final Predicate<List<List<PairedReadSequence>>> isPairsPoisonPill = (pairs -> (pairs.isEmpty()));

    // Check for poison pill into doSort method
    private final Predicate<List<SAMRecord>> isRecordPoisonPill = (samRecords -> (samRecords.isEmpty()));

    // Check if the pair is OK to add to sorter
    private final Predicate<PairedReadSequence> isPairValid = (pair) -> (pair.read1 != null
                                                                      && pair.read2 != null
                                                                      && pair.qualityOk);
    // Filter is Record valid or not
    private final Predicate<SAMRecord> isRecValid = (rec) -> (rec.getReadPairedFlag()
                                                          && (rec.getFirstOfPairFlag() || rec.getSecondOfPairFlag())
                                                          && !rec.isSecondaryOrSupplementary());
    /**
     * Functions
     */
    //
    protected final Function<QueuePeekableIterator<PairedReadSequence>,
            List<PairedReadSequence>> pairHandler = (iterator) ->
    {
        final List<PairedReadSequence> group = new ArrayList<>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        outer:
        while (iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();

            if(next == null)
                return group;

            for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
                if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i])
                    break outer;
            }
            group.add(iterator.next());
        }
        return group;
    };

    /**
     * Different sort algorithms
     */
    //Sorting with ExecutorService
    protected ElcSmartSortResponse doSmartSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        final long startTime = System.nanoTime();
        final List<SAMRecord> recordsPoisonPill = Collections.emptyList();

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

                        // If poison pill
                        if(rec.equals(recordsPoisonPill))
                            return;

                        pool.submit(() -> {
                            for (SAMRecord record : rec) {
                                // Read passes quality check if both ends meet the mean quality criteria
                                PairedReadSequence prs = processPairs(pendingByName.get(record),
                                        record,
                                        useBarcodes);
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
                    queueRecords.put(recordsPoisonPill);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            // Some useless work to wait for locks
            while (locker.get() != 0 || sorter.isSpillingToDisk()){
                if(locker.get() == 0 && !sorter.isSpillingToDisk())
                    break;
            }

            CloserUtil.close(in);
        }

        pool.shutdown();
        try                             { pool.awaitTermination(1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        log.info("SORTING - POOL-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime)/1000000)));
        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));

        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Sorting via Streams
    //IS NOT WORKING PROPERLY, some bugs in parallel process pairs, when fill bases or what ever
    protected ElcSmartSortResponse doStreamSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        final long startTime = System.nanoTime();
        final List<SAMRecord> recordsPoisonPill = Collections.emptyList();

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
                        PairedReadSequence prs = pendingByName.get(rec);
                        // prs.qualityOk check - partly skipped!!!
                        if(isPrsAndRecValid.test(new PrsAndRec(prs, rec)))
                            return processPairs(prs, rec, useBarcodes);
                        else
                            return new PairedReadSequence();
                    })
                    //.filter(isPrsAndRecValid)
                    //.map(pairAndRec -> processPairs(pairAndRec.prs, pairAndRec.rec))
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

    //Sorting with ForkJoinPool difference is in ForkJoinPool instead of doSmartSort
    protected ElcSmartSortResponse doPoolSort(final boolean useBarcodes) {
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

        final ForkJoinPool pool = new ForkJoinPool();
        final ConcurrentSupplier<SAMRecord> recordSupplier
                = new ConcurrentSupplier<>(RECORD_PROCESS_STACK_SIZE, 8);

        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());
            final ConcurrentPendingByNameCollection pendingByName
                    = new ConcurrentPendingByNameCollection(useBarcodes, opticalDuplicateFinder, readGroups);

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    final List<SAMRecord> rec = recordSupplier.getJob();

                    // If poison pill
                    if (rec.isEmpty())
                        return;

                    pool.submit(() -> {
                        for (SAMRecord record : rec) {
                            // Read passes quality check if both ends meet the mean quality criteria
                            PairedReadSequence prs = processPairs(pendingByName.get(record),
                                                                  record,
                                                                  useBarcodes);
                            progress.record(record);

                            if (isPairValid.test(prs))
                                sorter.add(prs);
                        }
                        recordSupplier.confirm();
                    });
                }
            });

            for (SAMRecord rec : in) {
                if (isRecValid.test(rec))
                    recordSupplier.add(rec);
            }

            recordSupplier.tryAddRest();
            recordSupplier.awaitConfirmation();
            sorter.awaitSpillingToDisk();

            CloserUtil.close(in);
        }

        recordSupplier.finish();
        pool.shutdown();
        try                             { pool.awaitTermination(1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        log.info("SORTING - SMART-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime) / 1000000)));
        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));

        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    //Sorting with ForkJoinPool and stream integration
    protected ElcSmartSortResponse doStreamPoolSort(final boolean useBarcodes){
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

        final ForkJoinPool pool = new ForkJoinPool();
        final ConcurrentSupplier<SAMRecord> recordSupplier
                = new ConcurrentSupplier<>(RECORD_PROCESS_STACK_SIZE, USED_THREADS);

        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());
            final ConcurrentPendingByNameCollection pendingByName
                    = new ConcurrentPendingByNameCollection(useBarcodes, opticalDuplicateFinder, readGroups);

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    final List<SAMRecord> rec = recordSupplier.getJob();
                    // If poison pill
                    if (isRecordPoisonPill.test(rec))
                        return;

                    pool.submit(() -> {
                        rec.stream()
                                .parallel()
                                .unordered()
                                .map(record -> {
                                    progress.record(record);
                                    return processPairs(pendingByName.get(record),
                                            record,
                                            useBarcodes);
                                })
                                .filter(isPairValid)
                                .forEach(sorter::add);
                        recordSupplier.confirm();
                    });
                }
            });

            for (SAMRecord rec : in) {
                if (isRecValid.test(rec))
                    recordSupplier.add(rec);
            }

/*            in.iterator().stream()
                    .parallel()
                    .unordered()
                    .filter(isRecValid)
                    .forEach(recordSupplier::add);*/

            recordSupplier.tryAddRest();
            recordSupplier.awaitConfirmation();
            sorter.awaitSpillingToDisk();

            CloserUtil.close(in);
        }

        recordSupplier.finish();
        pool.shutdown();
        try                             { pool.awaitTermination(1000, TimeUnit.SECONDS); }
        catch (InterruptedException e)  { e.printStackTrace(); }

        log.info("SORTING - POOL-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime)/1000000)));
        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));

        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    /**
     * ENCAPSULATED LOGIC
     */
    //Log invalid histogram prs
    protected void logInvalidGroup(List<PairedReadSequence> group, int meanGroupSize) {
        final PairedReadSequence prs = group.get(0);
        log.warn("Omitting group with over " + MAX_GROUP_RATIO + " times the expected mean number of read pairs. " +
                "Mean=" + meanGroupSize + ", Actual=" + group.size() + ". Prefixes: " +
                StringUtil.bytesToString(prs.read1, 0, MIN_IDENTICAL_BASES) +
                " / " +
                StringUtil.bytesToString(prs.read2, 0, MIN_IDENTICAL_BASES));
    }

    //Takes temporary object (PrsAndRec), process it, fill reads and returns PairedReadSequence
    private PairedReadSequence processPairs(final PairedReadSequence prs,
                                            final SAMRecord record,
                                            final boolean useBarcodes) {
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


    /**
     * CONSTANTS
     */
    // Stack size to process in sort method
    protected final int RECORD_PROCESS_STACK_SIZE = MAX_RECORDS_IN_RAM / 1000;

    //Group size to process job in doWork
    protected final int GROUP_PROCESS_STACK_SIZE = 50000;

    //
    protected final int QUEUE_CAPACITY = 4;

    //
    protected final int USED_THREADS = 8;

    /**
     * doWork main calculation method
     * Used to process histograms
     */
    protected int doWork() {
        IOUtil.assertFilesAreReadable(INPUT);

        final boolean useBarcodes   = (null != BARCODE_TAG
                                    || null != READ_ONE_BARCODE_TAG
                                    || null != READ_TWO_BARCODE_TAG);

        // Results from the doSort
        final ElcSmartSortResponse response = doPoolSort(useBarcodes);
        final long startTime = System.nanoTime();

        final ConcurrentSortingCollection<PairedReadSequence> sorter     = response.getSorter();
        final ProgressLogger                                  progress   = response.getProgress();
        final List<SAMReadGroupRecord>                        readGroups = response.getReadGroup();

        final QueuePeekableIterator<PairedReadSequence> iterator
                = new QueuePeekableIterator<>(sorter.iterator());

        final QueueIteratorProducer<PairedReadSequence, List<PairedReadSequence>> wrapper
                = new QueueIteratorProducer<>(iterator, pairHandler);

        final int meanGroupSize = (int) (Math.max(1, (progress.getCount() / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2)));
        final ConcurrentHistoCollection histoCollection = new ConcurrentHistoCollection(useBarcodes);
        final ConcurrentSupplier<List<PairedReadSequence>> groupSupplier
                = new ConcurrentSupplier<>(GROUP_PROCESS_STACK_SIZE, USED_THREADS);

        //final ExecutorService pool = Executors.newFixedThreadPool(USED_THREADS);
        final ExecutorService pool = Executors.newCachedThreadPool();
        final Object sync = new Object();

        // pool manager, receives job of groups, and make worker to process them
        // Now go through the sorted reads and attempt to find duplicates
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

        // Iterating through sorted groups, and making job to process them
        while (wrapper.hasNext()) {
            try                              { groupSupplier.add(wrapper.next()); }
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
        log.info("METRIC - EXECUTOR ELC (ms) : " + ((elcTime - startMetricFile / 1000000)));
        log.info("TOTAL - EXECUTOR ELC (ms) : " + (sortTime + (elcTime - startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }

    protected List<PairedReadSequence> tester(final QueuePeekableIterator<PairedReadSequence> iterator){
        Function<PairedReadSequence, List<PairedReadSequence>> caller = (PairedReadSequence first) -> {
            final List<PairedReadSequence> group = new ArrayList<>();
            //final PairedReadSequence first = iterator.next();
            group.add(first);

            outer:
            while (iterator.hasNext()) {
                final PairedReadSequence next = iterator.peek();
                for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
                    if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i])
                        break outer;
                }
                group.add(iterator.next());
            }
            return group;
        };
        return caller.apply(iterator.next());
    }

    protected List<PairedReadSequence> getNextGroup(final QueuePeekableIterator<PairedReadSequence> iterator) {
        final List<PairedReadSequence> group = new ArrayList<>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        outer:
        while (iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();
            for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
                if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i])
                    break outer;
            }
            group.add(iterator.next());
        }
        return group;

        /////////

    /*    final List<PairedReadSequence> group = new ArrayList<PairedReadSequence>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        outer:
        while (iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();
            for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
                if (first.read1[i] != next.read1[i] || first.read2[i] != next.read2[i])
                    break outer;
            }
            group.add(iterator.next());
        }
        return group;*/
    }
}
