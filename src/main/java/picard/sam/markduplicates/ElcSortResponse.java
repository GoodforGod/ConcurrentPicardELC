package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 18:28
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import picard.sam.markduplicates.util.ConcurrentSortingCollection;

import java.util.List;

/*
 * Result of the doSort method (ELC)
 */
public class ElcSortResponse
{
    private SortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter;
    private List<SAMReadGroupRecord> readGroup;
    private ProgressLogger progress;

    public ElcSortResponse(SortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter,
                           List<SAMReadGroupRecord> readGroup,
                           ProgressLogger progress) {
        this.sorter = sorter;
        this.progress = progress;
        this.readGroup = readGroup;
    }

    public SortingCollection<EstimateLibraryComplexity.PairedReadSequence> getSorter() {
        return sorter;
    }

    public List<SAMReadGroupRecord> getReadGroup() {
        return readGroup;
    }

    public ProgressLogger getProgress() {
        return progress;
    }
}
