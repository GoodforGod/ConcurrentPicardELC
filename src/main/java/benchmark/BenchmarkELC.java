package benchmark;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:17
 */


import picard.sam.markduplicates.EstimateLibraryComplexity;

/*
 * Entrance point for ELC benchmarks & tests
 */
public class BenchmarkELC
{
    private static void Initialize()
    {

    }

    public static void main(String[] args)
    {
        String[] arguments = new String[]
        {
            "EstimateLibraryComplexity",
            "\\",
            "I=" + BenchmarkVariables.IN_FILE_PATH//.replaceAll("\\\\", "/")
            ,"\\",
            "O=" + BenchmarkVariables.OUT_FILE_PATH//.replaceAll("\\\\", "/")
        };

        new EstimateLibraryComplexity().instanceMain(arguments);
    }

    private static void Filalize()
    {

    }
}
