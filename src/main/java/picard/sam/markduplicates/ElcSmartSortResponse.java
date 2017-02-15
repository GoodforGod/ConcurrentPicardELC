package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 15.02.2017
 * Time: 17:41
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.ProgressLogger;
import picard.sam.markduplicates.util.ConcurrentSortingCollection;

import java.util.List;

/*
 * Result of the doSmartSort, doStreamSort & etc methods (ELC)
 */
public class ElcSmartSortResponse
{
    private ConcurrentSortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter;
    private List<SAMReadGroupRecord> readGroup;
    private ProgressLogger progress;

    public ElcSmartSortResponse(ConcurrentSortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter,
                                List<SAMReadGroupRecord> readGroup,
                                ProgressLogger progress) {
        this.sorter = sorter;
        this.progress = progress;
        this.readGroup = readGroup;
    }

    public ConcurrentSortingCollection<EstimateLibraryComplexity.PairedReadSequence> getSorter() {
        return sorter;
    }

    public List<SAMReadGroupRecord> getReadGroup() {
        return readGroup;
    }

    public ProgressLogger getProgress() {
        return progress;
    }
}
