package benchmark;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:17
 */

import picard.sam.markduplicates.ConcurrentExecutorEstimateLibraryComplexity;
import picard.sam.markduplicates.ConcurrentStreamedEstimateLibraryComplexity;
import picard.sam.markduplicates.ConcurrentForkJoinPoolEstimateLibraryComplexity;
import picard.sam.markduplicates.EstimateLibraryComplexity;

/*
 * Entrance point for ELC benchmarks & tests
 */
public class BenchmarkELC
{
    private static void Initialize() {

    }

    public static void main(String[] args)
    {
        //new EstimateLibraryComplexity().instanceMain(BenchmarkVariables.small_bam_args);
        new ConcurrentStreamedEstimateLibraryComplexity().instanceMain(BenchmarkVariables.small_bam_args);
    }

    private static void Filalize() {

    }
}
