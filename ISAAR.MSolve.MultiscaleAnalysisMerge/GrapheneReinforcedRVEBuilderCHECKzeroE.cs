using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;


namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class GrapheneReinforcedRVEBuilderCHECKzeroE : IRVEbuilder
    {
        public GrapheneReinforcedRVEBuilderCHECKzeroE()
        { }

        public Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes()
        {
            return Reference2RVEExample1000withRenumberingwithInput_forMS();
        }

        private static Tuple<Model, Dictionary<int, Node>,double> Reference2RVEExample1000withRenumberingwithInput_forMS()
        {
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();

            //PROELEFSI public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
            double[,] Dq;
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integrationZeroE\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integrationZeroE\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integrationZeroE\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 6;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 6;
            int subdiscr1_shell = 3;
            int discr1_shell = 1;
            mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1; mp.hexa1 = 6; mp.hexa2 = 6; mp.hexa3 = 6;
            gp = mpgp.Item2;
            gp.E_shell = 0.000001;


            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
            double volume = mp.L01 * mp.L02 * mp.L03;

            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads         
            //RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            // commented out MS

            // model: add constraints
            //RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            // commented out MS

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);

            return new Tuple<Model, Dictionary<int, Node>,double>(model, boundaryNodes,volume);

        }

               
    }
}
