﻿using System;
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
using ISAAR.MSolve.Discretization.Interfaces;


namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataGPgrSh_v2 : IRVEbuilder // IdegenerateRVEbuilder
    {
        //PROELEFSI GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostDataGPgrSh
        //allages:i)undo allages sxetikes me to linear kai me to peripheral
        //compared Proelefsi to GrapheneReinforcedRVEBuilderExample5GrSh1RVEstif and determined which changes will be kept


        Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
        rveMatrixParameters mp;
        grapheneSheetParameters gp;
        grapheneSheetParameters gpNeo;
        string renumbering_vector_path;
        int RVE_id;
        int graphene_sheets_number;


        public GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataGPgrSh_v2(int RVE_id, grapheneSheetParameters gpNeo, int graphene_sheets_number)
        {
            this.RVE_id = RVE_id;
            this.gpNeo = gpNeo;
            this.graphene_sheets_number = graphene_sheets_number;
        }

        public IRVEbuilder Clone(int a) => new GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataGPgrSh_v2(a,gpNeo,graphene_sheets_number);
    
        public Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes()
        {
            return Reference2RVEExample10000withRenumberingwithInput_forMS();
        }

        private Tuple<Model, Dictionary<int, Node>,double> Reference2RVEExample10000withRenumberingwithInput_forMS()
        {
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();

            //PROELEFSI public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
            double[,] Dq;
            //Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            //rveMatrixParameters mp;
            //grapheneSheetParameters gp;
            var rve_id_data = RVE_id.ToString();
            string nGraphenes = graphene_sheets_number.ToString();
            renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\RVE_database"+nGraphenes+ @"\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
            renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\RVE_database" + nGraphenes + @"\RVE_database\rve_no_{0}\Fxk_p_komvoi_rve.txt";
            Fxk_p_komvoi_rve_path = string.Format(Fxk_p_komvoi_rve_path, rve_id_data);
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\RVE_database" + nGraphenes + @"\RVE_database\rve_no_{0}\o_xsunol_gs_";
            o_xsunol_input_path_gen = string.Format(o_xsunol_input_path_gen, rve_id_data);
            o_xsunol_input_path_gen = o_xsunol_input_path_gen + "{0}.txt";
            int subdiscr1 = 2;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 8;
            int subdiscr1_shell = 6;
            int discr1_shell = 1;
            mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
            gp = mpgp.Item2;
            if (!(gpNeo == null))
            { gp.D_f_1 = gpNeo.D_f_1; gp.D_f_3 = gpNeo.D_f_3; gp.D_o_1 = gpNeo.D_o_1; gp.D_o_3 = gpNeo.D_o_3; gp.T_o_1 = gpNeo.T_o_1; gp.T_o_3 = gpNeo.T_o_3; gp.n_curve = gpNeo.n_curve; }

            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];
            //double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };


            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
            double volume = mp.L01 * mp.L02 * mp.L03;

            int hexaElementsNumber = model.ElementsDictionary.Count();

            //IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_withRenumberingBondSlip(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
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
            //var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);

            //var CohesiveGroupings = new EmbeddedCohesiveGrouping[EmbElementsIds.GetLength(0)];

            var hostSubGroups = new Dictionary<int, IEnumerable<Element>>();
            for (int i1 = 0; i1 < EmbElementsIds.GetLength(0); i1++)
            {
                hostSubGroups.Add(EmbElementsIds[i1], FEMMeshBuilder.GetHostGroupForCohesiveElement(model.ElementsDictionary[EmbElementsIds[i1]], mp, model, renumbering_vector_path));
                //var embeddedGroup_i1 = new List<Element>(1) { model.ElementsDictionary[EmbElementsIds[i1]] };
                //CohesiveGroupings[i1] = new EmbeddedCohesiveGrouping(model, hostGroup_i1, embeddedGroup_i1);
            }

            var CohesiveGroupping = new EmbeddedCohesiveSubGrouping(model, hostSubGroups, embdeddedGroup);

            return new Tuple<Model, Dictionary<int, Node>,double>(model, boundaryNodes,volume);

        }

        // PROSOXH DEN ARKEI MONO TO PARAKATW NA GINEI UNCOMMENT WSTE NA GINEI IMPLEMENT TO IDegenerateRVEBuilder 
        //xreiazetai kai na xrhsimopoithei h katallhlh methodos tou femmeshbuilder gia to model and boundary nodes na dinei mono ta peripheral
        //public Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model model)
        //{
        //    return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
        //    //TODO:  Pithanws na epistrefetai apo GetModelAndBoundaryNodes ... AndConstraints.
        //}

    }
}
