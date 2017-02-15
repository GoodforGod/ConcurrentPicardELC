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

import java.io.File;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;

import static java.lang.Math.pow;

/*
 * Thread Executor version of the ELC
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ThreadExecutorEstimateLibraryComplexity extends EstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ThreadExecutorEstimateLibraryComplexity.class);

    public ThreadExecutorEstimateLibraryComplexity() {
        super();
    }

    //  Temporary structure, used in sorting via StreamAPI (doStreamSort)
    protected class PairAndRec {
        public SAMRecord rec;
        public PairedReadSequence prs;
        public PairedReadSequenceWithBarcodes prsWithBarcodes;

        public PairAndRec(SAMRecord rec, PairedReadSequence prs)
        {
            this.rec = rec;
            this.prs = prs;
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

    // Takes temporary object (PairAndRec), process it, fill reads and returns PairedReadSequence
    protected PairedReadSequence FillPairedSequence(PairAndRec prsAndRec) {
        final byte[] bases = prsAndRec.rec.getReadBases();

        if (prsAndRec.rec.getReadNegativeStrandFlag())
            SequenceUtil.reverseComplement(bases);

        if (prsAndRec.rec.getFirstOfPairFlag())
            prsAndRec.prs.read1 = bases;
        else
            prsAndRec.prs.read2 = bases;

        return prsAndRec.prs;
    }

    // Sorting via StreamAPI
    protected ElcSmartSortResponse doStreamSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final Object sync = new Object();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final ConcurrentSortingCollection<PairedReadSequence> sorter;
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                !useBarcodes
                        ? new PairedReadCodec()
                        : new PairedReadWithBarcodesCodec(),
                new PairedReadComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);

        AtomicInteger count = new AtomicInteger(0);

        // Filter is Record valid or not
        Predicate<SAMRecord> isRecValid = (rec) -> (rec.getReadPairedFlag()
                                                && (rec.getFirstOfPairFlag() || rec.getSecondOfPairFlag())
                                                && !rec.isSecondaryOrSupplementary());

        Predicate<PairAndRec> isPairAndRecValid = (prsAndRec) -> (prsAndRec.rec != null
                                                                && prsAndRec.prs != null
                                                                && prsAndRec.prs.qualityOk
                                                                && prsAndRec.rec.getReadBases() != null
                                                                && passesQualityCheck(prsAndRec.rec.getReadBases(),
                                                                                      prsAndRec.rec.getBaseQualities(),
                                                                                      MIN_IDENTICAL_BASES,
                                                                                      MIN_MEAN_QUALITY));
        // Loop through the input files and pick out the read sequences etc.
        for (final File f : INPUT) {
            final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            long startRecLoop = System.nanoTime();

            BlockingQueue<PairedReadSequence> parallelListPRS = new LinkedBlockingQueue<>();

             in.iterator().stream()
                        .parallel()
                        .unordered()
                        .filter(isRecValid)
                        .map(rec -> new PairAndRec(rec, pendingByName.remove(rec.getReadName())))
                        .map(prsAndRec -> {
                            if (prsAndRec.prs == null) {
                                // Make a new paired read object and add RG and physical location information to it
                                prsAndRec.prs = useBarcodes
                                              ? new PairedReadSequenceWithBarcodes()
                                              : new PairedReadSequence();

                                if (opticalDuplicateFinder.addLocationInformation(prsAndRec.rec.getReadName(), prsAndRec.prs)) {
                                    final SAMReadGroupRecord rg = prsAndRec.rec.getReadGroup();
                                    if (rg != null)
                                        prsAndRec.prs.setReadGroup((short) readGroups.indexOf(rg));
                                }
                                prsAndRec.prs = pendingByName.put(prsAndRec.rec.getReadName(), prsAndRec.prs);
                            }
                            prsAndRec.prsWithBarcodes = useBarcodes
                                                      ? (PairedReadSequenceWithBarcodes) prsAndRec.prs
                                                      : null;
                            progress.record(prsAndRec.rec);
                            return prsAndRec;
                        })
                        .filter(isPairAndRecValid)
                        .map(this::FillPairedSequence)
                        .forEachOrdered(parallelListPRS::add);

            CloserUtil.close(in);
        }

        log.info("FILL SORTING - STREAM-SORT ELC (ms) : "
                + (sortTime = (double) ((System.nanoTime() - startTime) / 1000000)));

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));

        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Sorting with ExecutorService
    protected ElcSmartSortResponse doSmartSort(final boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final ConcurrentSortingCollection<PairedReadSequence> sorter;
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                !useBarcodes
                        ? new PairedReadCodec()
                        : new PairedReadWithBarcodesCodec(),
                new PairedReadComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);


        // Loop through the input files and pick out the read sequences etc.

        final int recordListSize = MAX_RECORDS_IN_RAM / 1000;
        final ExecutorService pool = Executors.newCachedThreadPool();
        final BlockingQueue<List<SAMRecord>> queueRecords = new ArrayBlockingQueue<>(10);
        List<SAMRecord> records = new ArrayList<>(recordListSize);

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);

        final Object sync = new Object();

        for (final File f : INPUT) {
            final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    try {
                        final List<SAMRecord> rec = queueRecords.take();

                        if (rec != null) {
                            pool.submit(() -> {
                                for (SAMRecord record : rec) {
                                    PairedReadSequence prs = pendingByName.remove(record.getReadName());

                                    if (prs == null) {
                                        // Make a new paired read object and add RG and physical location information to it
                                        prs = useBarcodes
                                                ? new PairedReadSequenceWithBarcodes()
                                                : new PairedReadSequence();

                                        if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), prs)) {
                                            final SAMReadGroupRecord rg = record.getReadGroup();
                                            if (rg != null)
                                                prs.setReadGroup((short) readGroups.indexOf(rg));
                                        }
                                        pendingByName.put(record.getReadName(), prs);
                                    }

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

                                    if (prs.read1 != null && prs.read2 != null && prs.qualityOk) {
                                        sorter.add(prs);
                                    }
                                    progress.record(record);
                                }
                                locker.decrementAndGet();
                            });
                        }
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            });

            for (SAMRecord rec : in) {
                if (!rec.getReadPairedFlag()
                        || !rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag()
                        || rec.isSecondaryOrSupplementary())
                    continue;

                records.add(rec);

                if (records.size() <= recordListSize)
                    continue;

                try {
                    queueRecords.put(records);
                    locker.incrementAndGet();
                    records = new ArrayList<>(recordListSize);
                } catch (InterruptedException e) {
                    e.printStackTrace();
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

        log.info("SORTING - SMART-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime)/1000000)));

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));
        return new ElcSmartSortResponse(sorter, readGroups, progress);
    }

    // Sorting with ForkJoinPool difference is in ForkJoinPool instead of doSmartSort
    protected ElcSmartSortResponse doPoolSort(boolean useBarcodes) {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final ConcurrentSortingCollection<PairedReadSequence> sorter;
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        sorter = ConcurrentSortingCollection.newInstance(PairedReadSequence.class,
                !useBarcodes
                        ? new PairedReadCodec()
                        : new PairedReadWithBarcodesCodec(),
                new PairedReadComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);


        // Loop through the input files and pick out the read sequences etc.

        final int recordListSize = MAX_RECORDS_IN_RAM / 1000;
        final ForkJoinPool pool = new ForkJoinPool();
        final BlockingQueue<List<SAMRecord>> queueRecords = new ArrayBlockingQueue<>(10);
        List<SAMRecord> records = new ArrayList<>(recordListSize);

        // Countdown for active threads (cause Executor is kinda stuck in awaitTermination, don't know why)
        final AtomicInteger locker = new AtomicInteger(0);

        final Object sync = new Object();

        for (final File f : INPUT) {
            final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            pool.execute(() -> {
                while (!Thread.interrupted()) {
                    try {
                        final List<SAMRecord> rec = queueRecords.take();

                        if (rec != null) {
                            pool.submit(() -> {
                                for (SAMRecord record : rec) {
                                    PairedReadSequence prs = pendingByName.remove(record.getReadName());

                                    if (prs == null) {
                                        // Make a new paired read object and add RG and physical location information to it
                                        prs = useBarcodes
                                                ? new PairedReadSequenceWithBarcodes()
                                                : new PairedReadSequence();

                                        if (opticalDuplicateFinder.addLocationInformation(record.getReadName(), prs)) {
                                            final SAMReadGroupRecord rg = record.getReadGroup();
                                            if (rg != null)
                                                prs.setReadGroup((short) readGroups.indexOf(rg));
                                        }
                                        pendingByName.put(record.getReadName(), prs);
                                    }

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

                                    if (prs.read1 != null && prs.read2 != null && prs.qualityOk) {
                                        sorter.add(prs);
                                    }
                                    progress.record(record);
                                }
                                locker.decrementAndGet();
                            });
                        }
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            });

            for (SAMRecord rec : in) {
                if (!rec.getReadPairedFlag()
                        || !rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag()
                        || rec.isSecondaryOrSupplementary())
                    continue;

                records.add(rec);

                if (records.size() < recordListSize)
                    continue;

                try {
                    queueRecords.put(records);
                    locker.incrementAndGet();
                    records = new ArrayList<>(recordListSize);
                } catch (InterruptedException e) {
                    e.printStackTrace();
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
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<PairedReadSequence>(sorter.iterator());

        final Map<String, Histogram<Integer>> opticalHistosByLibrary     = new ConcurrentHashMap<String, Histogram<Integer>>();
        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new ConcurrentHashMap<String, Histogram<Integer>>();

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

                if(groupStack.size() >= 50000)
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

        log.info("SORTER PROCESS - (ms) : " + (double)(System.nanoTime() - startSortIterateTime)/ 1000000);

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
        log.info("METRIC - THREAD POOL ELC (ms) : " + ((elcTime - (double)startMetricFile / 1000000)));
        log.info("TOTAL - THREAD POOL ELC (ms) : " + (sortTime + (elcTime - (double)startTime / 1000000)));

        file.write(OUTPUT);
        return 0;
    }
}
