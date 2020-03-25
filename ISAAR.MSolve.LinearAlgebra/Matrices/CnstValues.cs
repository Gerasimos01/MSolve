using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public class CnstValues
    {
        public static bool run_plate = false;

        public static int[,] diakritopoihsh_examples = new int[,]
        {{2,6,},{3,6},{4,6},{5,6},{5,6},{6,6},
            {4,2 },{4,3},{4,4 },{4,5},
            {3,3 },{5,5} };

        public static int current_example_no;

        public static int exampleNo { get; set; }

        bool runCluster = false;

        public static bool runOnlyHexaModel { get; set; } = false;
        //public string exampleOutputPathGen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example1\input_matlab";

        public string exampleOutputPathGen
        {
            get
            {
                if (runCluster) { return @"C:\Users\cluster\Documents\Large_rves\examples\example" + exampleNo + @"\input_matlab"; }
                else { return @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\examples\example" + exampleNo + @"\input_matlab"; }
            }
        }
        public string solverPath { get { return exampleOutputPathGen + @"\model_overwrite\subdomain_data_solver"; } }

        public string debugString = @"C:\Users\acivi\Documents\notes_elegxoi_2\develp_3D";

        public string exampleDiscrInputPathGen { get { return exampleOutputPathGen + @"\Msolve_input"; } }

        public string interfaceSolverStatsPath { get { return exampleOutputPathGen + @"\Msolve_solution"; } }

        public string rand_data_vec_path
        {
            get
            {
                if (runCluster) { return @"C:\Users\cluster\Documents\Large_rves\files_from_old_paths\rand_data.txt"; }
                else { return @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_random_geometry_Msolve\REF2_50_000_renu_new_multiple_algorithms_check_develop_copy_for_progr_random_direct_in_C\rand_data.txt"; }
            }
        }
        public CnstValues()
        {

        }

        public static bool printInterfaceSolutionStats { get; set; } = true;



        public static bool printSolver_run_overwrite_data_region { get; set; } = true;

        public static bool printGlobalSolutionStats { get; set; } = true;

        public static bool printPcgMatRhsEtc_AndInterfaceProblemStats { get; set; } = false; // 

        public static bool printPreconditoner { get; set; } = false;//true; // 

        public void WriteToFileStringArray(string[] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {

                writer.Write(array[i]);
                //writer.Write(' ');

                writer.WriteLine();
            }
            writer.Flush();

            writer.Dispose();
        }

        public static void preventOutputFileWrite()
        {
            SaveDefaultValuesForRestore();
            printInterfaceSolutionStats = false;
            printSolver_run_overwrite_data_region = false;
            printGlobalSolutionStats = false;
            printPcgMatRhsEtc_AndInterfaceProblemStats = false;
            printPreconditoner = false;
        }

        private static bool[] defaultBooleanValues;
        private static void SaveDefaultValuesForRestore()
        {
            defaultBooleanValues = new bool[] {
                printInterfaceSolutionStats ,
                printSolver_run_overwrite_data_region,
                printGlobalSolutionStats,
                printPcgMatRhsEtc_AndInterfaceProblemStats,
                printPreconditoner};
        }

        public static void RestoreDefaultBoolValues()
        {
            printInterfaceSolutionStats = defaultBooleanValues[0];
            printSolver_run_overwrite_data_region = defaultBooleanValues[1];
            printGlobalSolutionStats = defaultBooleanValues[2];
            printPcgMatRhsEtc_AndInterfaceProblemStats = defaultBooleanValues[3];
            printPreconditoner = defaultBooleanValues[4];
        }
    }
}
