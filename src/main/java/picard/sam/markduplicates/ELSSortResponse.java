package picard.sam.markduplicates;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 18:28
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.util.List;

/*
 * DEFAULT COMMENT
 */
public class ELSSortResponse
{
    private SortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter;
    private List<SAMReadGroupRecord> readGroup;
    private ProgressLogger progress;

    public ELSSortResponse(SortingCollection<EstimateLibraryComplexity.PairedReadSequence> sorter,
                            List<SAMReadGroupRecord> readGroup,
                           ProgressLogger progress)
    {
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
