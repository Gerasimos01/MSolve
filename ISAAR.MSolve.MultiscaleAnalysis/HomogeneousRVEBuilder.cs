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
    public class HomogeneousRVEBuilder : IRVEbuilder
    {
        //TODOGerasimos gia na ta krataei mesa kai na kanei build model oses fores tou zhththei
        // omoiws na ginei kai to RVE me graphene sheets 
        // string renumbering_vector_path; 
        // int subdiscr1;
       

        public HomogeneousRVEBuilder()
        {
            //TODOGerasimos
            // this.renumbering_vector_path=renumbering_vector_path,
            // this.subdiscr1=subdiscr1
        }

        public Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes()
        {
           return Reference2RVEExample10_000withRenumbering_mono_hexa();
        }

        public static Tuple<Model, Dictionary<int, Node>,double> Reference2RVEExample10_000withRenumbering_mono_hexa()
        {
            Model model= new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            Dictionary<int, Node> boundaryNodes= new Dictionary<int, Node>();
            // COPY APO: Reference2RVEExample100_000withRenumbering_mono_hexa
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            //C:\Users\cluster 5\Desktop\Gerasimos\REFERENCE_Examples_me_develop\10_000_mono_hexa\REF_new_total_numbering.txt einai link sto PC LAB
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_mono_hexa\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_mono_hexa\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;

            mpgp = FEMMeshBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];


            FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
            double volume = mp.L01 * mp.L02 * mp.L03;

            // MS: oi epomenes 6 grammes aforoun embedding commented out 
            //int hexaElementsNumber = model.ElementsDictionary.Count();
            //IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            //List<int> EmbeddedElementsIDs = new List<int>();
            //int element_counter_after_Adding_sheet;
            //element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            //int shellElementsNumber;


            //  MS: to loop commented out afou  graphene_sheets_number=0
            //for (int j = 0; j < graphene_sheets_number; j++)
            //{
            //    RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
            //    shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
            //    //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
            //    //embdeddedGroup.Concat(embdeddedGroup_adittion);
            //    for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
            //    {
            //        EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
            //    }
            //    element_counter_after_Adding_sheet = model.ElementsDictionary.Count();
            //}

            // model: add loads
            //RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            // commented out MS

            // model: add constraints
            //RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            // commented out MS

            // MS: den uparxei pia embedding oi 3 epomenes grammes commented out 
            //int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            //IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            //var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);

            return new Tuple<Model, Dictionary<int, Node>,double>(model, boundaryNodes,volume);
        }
    }
}
