package benchmark;

/*
 * Created by GoodforGod 
 * Date: 09.02.2017
 * Time: 17:21
 */

/*
 * Variables & consts used in benchmarking
 */
public class BenchmarkVariables
{
    public static final int CORES = 4;

    public static final String IN_FILE_PATH = "C:\\data\\bam\\na_s.bam";
    public static final String OUT_FILE_PATH = "C:\\data\\bam\\out_na_s.txt";

    public static final String[] small_bam_args = new String[]
                        {
                                "I=" + BenchmarkVariables.IN_FILE_PATH,
                                "O=" + BenchmarkVariables.OUT_FILE_PATH
                        };

}
