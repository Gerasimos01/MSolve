using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Elements;

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    class SeparateCodeCheckingClass
    {
        private static void check01()
        {
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double sigma_f= 0.2; // apo to arxeio create_random_data_for_geom_programing_in_C tou fakelou tou parapanw rand data vec path
            int n_graphene_sheets = 10; // omoiws

            string rand_data_vec_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria\REF2_50_000_renu_new_multiple_algorithms_check_develop_copy_for_progr_random_direct_in_C\rand_data.txt";
            savedRandomDataClass a = new savedRandomDataClass(PrintUtilities.ReadVector(rand_data_vec_path));
            Tuple<double[], double[], double[][]> RandomDataForGeomGiaSugkekrimenoRand = RandomOrientations.CreateRandomDataForGeom(n_graphene_sheets, gp, mp, sigma_f, a);
        }
    }
}
