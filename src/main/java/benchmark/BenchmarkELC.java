package benchmark;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:17
 */

import picard.sam.markduplicates.ConcurrentExecutorEstimateLibraryComplexity;

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
        new ConcurrentExecutorEstimateLibraryComplexity().instanceMain(BenchmarkVariables.small_bam_args);
    }

    private static void Filalize() {

    }
}
