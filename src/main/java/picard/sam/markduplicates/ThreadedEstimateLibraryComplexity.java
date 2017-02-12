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

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static java.lang.Math.pow;

/*
 * Threaded version of the ELC
 */
@CommandLineProgramProperties(
        usage = EstimateLibraryComplexity.USAGE_SUMMARY + EstimateLibraryComplexity.USAGE_DETAILS,
        usageShort = EstimateLibraryComplexity.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ThreadedEstimateLibraryComplexity extends EstimateLibraryComplexity
{
    protected final Log log = Log.getInstance(ThreadedEstimateLibraryComplexity.class);

    public class PairAndRec
    {
        public SAMRecord rec;
        public PairedReadSequence prs;
        public PairedReadSequenceWithBarcodes prsWithBarcodes;

        public PairAndRec(SAMRecord rec, PairedReadSequence prs)
        {
            this.rec = rec;
            this.prs = prs;
        }
    }

    public ThreadedEstimateLibraryComplexity() {
        super();
    }

    protected ELCSortResponse doSmartSort(final boolean useBarcodes)
    {
        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        long startTime = System.nanoTime();

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final SortingCollection<PairedReadSequence> sorter;
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");

        if (!useBarcodes)
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);
        else
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadWithBarcodesCodec(),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);

        // Loop through the input files and pick out the read sequences etc.
        final ExecutorService executorService = Executors.newCachedThreadPool();

        for (final File f : INPUT)
        {
            final Map<String, PairedReadSequence> pendingByName = new ConcurrentHashMap<>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(f);
            readGroups.addAll(in.getFileHeader().getReadGroups());

            Optional<PairedReadSequence> prsOptional = Optional.ofNullable(null);

            Supplier<PairedReadSequence> prsSupplier = () -> {
                return useBarcodes ? new PairedReadSequenceWithBarcodes() : new PairedReadSequence();
            };

            /*
            List<SAMRecord> pairAndRecs = in.iterator().stream()
                    .unordered()
                    .filter(rec -> (rec.getReadPairedFlag()
                                && rec.getFirstOfPairFlag()
                                && rec.getSecondOfPairFlag()
                                && !rec.isSecondaryOrSupplementary()))
                    .map(rec -> new PairAndRec(rec, pendingByName.remove(rec.getReadName())))
                    .map(prsAndRec -> {
                        if (prsAndRec.prs == null)
                        {
                            // Make a new paired read object and add RG and physical location information to it
                            prsAndRec.prs = useBarcodes ? new PairedReadSequenceWithBarcodes() : new PairedReadSequence();

                            if (opticalDuplicateFinder.addLocationInformation(prsAndRec.rec.getReadName(), prsAndRec.prs))
                            {
                                final SAMReadGroupRecord rg = prsAndRec.rec.getReadGroup();
                                if (rg != null)
                                    prsAndRec.prs.setReadGroup((short) readGroups.indexOf(rg));
                            }
                            pendingByName.put(prsAndRec.rec.getReadName(), prsAndRec.prs);
                        }
                        prsAndRec.prsWithBarcodes = (useBarcodes)
                                ? (PairedReadSequenceWithBarcodes) prsAndRec.prs
                                : null;

                        //progress.record(prsAndRec.rec);

                        return prsAndRec;
                    })
                    .filter(prsAndRec -> prsAndRec.prs.qualityOk
                            && prsAndRec.rec.getReadBases() != null
                            && passesQualityCheck(prsAndRec.rec.getReadBases(),
                                                  prsAndRec.rec.getBaseQualities(),
                                                  MIN_IDENTICAL_BASES,
                                                  MIN_MEAN_QUALITY))
                    .map(prsAndRec -> {
                        final byte[] bases = prsAndRec.rec.getReadBases();

                        if (prsAndRec.rec.getReadNegativeStrandFlag())
                            SequenceUtil.reverseComplement(bases);

                        if (prsAndRec.rec.getFirstOfPairFlag())
                            prsAndRec.prs.read1 = bases;
                        else
                            prsAndRec.prs.read2 = bases;

                        return prsAndRec;
                    })
                    .forEach(prsAndRec -> sorter.add(prsAndRec.prs));
*/

            for (final SAMRecord rec : in)
            {
                // TO LOG TimeInLoop
                long startRecLoop = System.nanoTime();

                if (!rec.getReadPairedFlag()
                    || (!rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag())
                    || rec.isSecondaryOrSupplementary())
                        continue;

                PairedReadSequence prs = pendingByName.remove(rec.getReadName());

                if (prs == null)
                {
                    // Make a new paired read object and add RG and physical location information to it
                    prs = useBarcodes ? new PairedReadSequenceWithBarcodes() : new PairedReadSequence();

                    if (opticalDuplicateFinder.addLocationInformation(rec.getReadName(), prs))
                    {
                        final SAMReadGroupRecord rg = rec.getReadGroup();
                        if (rg != null)
                            prs.setReadGroup((short) readGroups.indexOf(rg));
                    }
                    pendingByName.put(rec.getReadName(), prs);
                }

                // Read passes quality check if both ends meet the mean quality criteria
                prs.qualityOk = prs.qualityOk && passesQualityCheck(rec.getReadBases(),
                        rec.getBaseQualities(),
                        MIN_IDENTICAL_BASES,
                        MIN_MEAN_QUALITY);

                final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes)
                        ? (PairedReadSequenceWithBarcodes) prs
                        : null;

                // Get the bases and restore them to their original orientation if necessary
                final byte[] bases = rec.getReadBases();

                if (rec.getReadNegativeStrandFlag())
                    SequenceUtil.reverseComplement(bases);

                if (rec.getFirstOfPairFlag())
                {
                    prs.read1 = bases;

                    if (useBarcodes) {
                        prsWithBarcodes.barcode = getBarcodeValue(rec);
                        prsWithBarcodes.readOneBarcode = getReadOneBarcodeValue(rec);
                    }
                } else {
                    prs.read2 = bases;

                    if (useBarcodes)
                        prsWithBarcodes.readTwoBarcode = getReadTwoBarcodeValue(rec);
                }

                if (prs.read1 != null && prs.read2 != null && prs.qualityOk)
                    sorter.add(prs);

                log.info("REC - SMART-SORT ELC (microsecond) : "
                        + ((double)((System.nanoTime() - startRecLoop)/1000)));

                progress.record(rec);
            }
            CloserUtil.close(in);
        }

        log.info("SORTING - SMART-SORT ELC (ms) : "
                + (sortTime = (double)((System.nanoTime() - startTime)/1000000)));

        log.info(String.format("Finished reading - read %d records - moving on to scanning for duplicates.", progress.getCount()));
        return new ELCSortResponse(sorter, readGroups, progress);
    }


    /*
     * WORK
     */
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

        final Map<String, Histogram<Integer>> opticalHistosByLibrary = new HashMap<String, Histogram<Integer>>();
        final Map<String, Histogram<Integer>> duplicationHistosByLibrary = new HashMap<String, Histogram<Integer>>();

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

        while (iterator.hasNext())
        {
            // Get the next group and split it apart by library
            final List<PairedReadSequence> group = getNextGroup(iterator);

            if (group.size() > meanGroupSize * MAX_GROUP_RATIO)
            {
                final PairedReadSequence prs = group.get(0);

                log.warn("Omitting group with over "
                        + MAX_GROUP_RATIO + " times the expected mean number of read pairs. "
                        + "Mean=" + meanGroupSize
                        + ", Actual=" + group.size()
                        + ". Prefixes: " + StringUtil.bytesToString(prs.read1, 0, MIN_IDENTICAL_BASES)
                        + " / "
                        + StringUtil.bytesToString(prs.read2, 0, MIN_IDENTICAL_BASES));
            }
            else
            {
                final Map<String, List<PairedReadSequence>> sequencesByLibrary = splitByLibrary(group, readGroups);

                // Now process the reads by library
                for (final Map.Entry<String, List<PairedReadSequence>> entry : sequencesByLibrary.entrySet())
                {
                    final String library = entry.getKey();
                    final List<PairedReadSequence> seqs = entry.getValue();

                    Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);

                    Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);

                    if (duplicationHisto == null)
                    {

                        duplicationHisto = new Histogram<>("duplication_group_count", library);
                        opticalHisto = new Histogram<>("duplication_group_count", "optical_duplicates");
                        duplicationHistosByLibrary.put(library, duplicationHisto);
                        opticalHistosByLibrary.put(library, opticalHisto);
                    }

                    algorithmResolver.resolveAndSearch(seqs, duplicationHisto, opticalHisto);
                }
                ++groupsProcessed;

                if (lastLogTime < System.currentTimeMillis() - 60000)
                {
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
