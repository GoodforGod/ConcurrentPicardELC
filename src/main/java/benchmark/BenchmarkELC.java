package benchmark;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:17
 */

import picard.sam.markduplicates.ConcurrentPoolEstimateLibraryComplexity;
import picard.sam.markduplicates.ConcurrentStreamedEstimateLibraryComplexity;
import picard.sam.markduplicates.ConcurrentExecutorEstimateLibraryComplexity;
import picard.sam.markduplicates.EstimateLibraryComplexity;

/*
 * Entrance point for ELC benchmarks & tests
 */
public class BenchmarkELC
{
    public static void main(String[] args) {
        new ConcurrentStreamedEstimateLibraryComplexity().instanceMain(BenchmarkVariables.small_bam_args);
    }
}
