﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.MSAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Creates a elastic matrix rve with embedded graphene sheets, with the asumption of damage behaviour for the material interface.
    /// Use of model separation methods is made.
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte : IRVEbuilder //IdegenerateRVEbuilder
    {
        //origin: RveGrShMultipleSeparatedDevelop
        //changes: discretization 

        //GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataDdm
        //Origin branch: example/ms_development_nl_elements_merge (xwris sto telos )
        // modifications update se v2

        public int[] hexaPrint { get; private set; }
        public int[] cohePrint { get; private set; }
        public int[] shellPrint { get; private set; }
        public Dictionary<int, double[]> CornerNodesIds { get; private set; }
        public Dictionary<int, int[]> CornerNodesIdAndsubdomains { get; private set; }
        public Dictionary<int, int[]> subdFreeBRNodes { get; private set; }
        public string subdomainOutputPath { get; private set; }
        public IList<Node> EmbeddedNodes { get; private set; }
        public Dictionary<ISubdomain, List<Node>> RveMatrixSubdomainInnerNodes { get; private set; }

        private bool decomposeModel;
        public Dictionary<int, HashSet<INode>> cornerNodes;
        public ISolver GetAppropriateSolver(Model model)
        {
            if (decomposeModel)
            {
                //Setup solver
                var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
                interfaceSolverBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1);
                interfaceSolverBuilder.PcgConvergenceTolerance = 1E-10;
                var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
                //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
                var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
                var fetiSolverBuilder = new FetiDPSolver.Builder(cornerNodeSelection, fetiMatrices);
                fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
                fetiSolverBuilder.ProblemIsHomogeneous = false;
                fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();
                FetiDPSolver fetiSolver = fetiSolverBuilder.BuildSolver(model);
                return fetiSolver;
            }
            else
            {
                return (new SkylineSolver.Builder()).BuildSolver(model);
            }
        }

        public Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
        public rveMatrixParameters mp;
        public grapheneSheetParameters gp;
        public string renumbering_vector_path;
        int RVE_id;
        public int[][] CornerNodesData;

        public RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(int RVE_id, bool decomposeModel)
        {
            this.RVE_id = RVE_id;
            this.decomposeModel = decomposeModel;
        }

        public IRVEbuilder Clone(int a) => new RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(a, decomposeModel);

        public Tuple<Model, Dictionary<int, Node>, double> GetModelAndBoundaryNodes()
        {
            return Reference2RVEExample10000withRenumberingwithInput_forMS();
        }

        private Tuple<Model, Dictionary<int, Node>, double> Reference2RVEExample10000withRenumberingwithInput_forMS()
        {
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();

            //Origin public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
            double[,] Dq;
            //Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            //rveMatrixParameters mp;
            //grapheneSheetParameters gp;
            var rve_id_data = RVE_id.ToString();

            int path = 2;

            //renumbering_vector_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\REF_new_total_numbering.txt";
            if (path == 1) renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
            if (path == 2) renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2c_alte\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
            renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);

            //string Fxk_p_komvoi_rve_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\Fxk_p_komvoi_rve.txt";
            string Fxk_p_komvoi_rve_path = null;
            if (path == 1) Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\Fxk_p_komvoi_rve.txt";
            if (path == 2) Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2c_alte\RVE_database\rve_no_{0}\Fxk_p_komvoi_rve.txt";
            Fxk_p_komvoi_rve_path = string.Format(Fxk_p_komvoi_rve_path, rve_id_data);


            //string o_xsunol_input_path_gen = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\o_xsunol_gs_";
            string o_xsunol_input_path_gen = null;
            if (path == 1) o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}\o_xsunol_gs_";
            if (path == 2) o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2c_alte\RVE_database\rve_no_{0}\o_xsunol_gs_";
            o_xsunol_input_path_gen = string.Format(o_xsunol_input_path_gen, rve_id_data);
            o_xsunol_input_path_gen = o_xsunol_input_path_gen + "{0}.txt";
            string subdomainOutputPath_gen = null;
            if (path == 1) subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2c\RVE_database\rve_no_{0}";
            if (path == 2) subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2c_alte\RVE_database\rve_no_{0}";
            subdomainOutputPath = string.Format(subdomainOutputPath_gen, rve_id_data);
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 16;
            int subdiscr1_shell = 6;
            int discr1_shell = 1;
            mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
            gp = mpgp.Item2;


            int graphene_sheets_number = 3;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];
            double[][] ekk_xyz = new double[graphene_sheets_number][];

            int[][] CornerNodesData = new int[27][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            int thesi = 0;
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int i3 = 0; i3 < 3; i3++)
                    {
                        CornerNodesData[thesi] = new int[3] { (5 - 1) + i1 * 4, (5 - 1) + i2 * 4, (5 - 1) + i3 * 4 };
                        thesi++;
                    }
                }
            }

            Dq = new double[9, 3 * (((mp.hexa1 + 1) * (mp.hexa2 + 1) * (mp.hexa3 + 1)) - ((mp.hexa1 - 1) * (mp.hexa2 - 1) * (mp.hexa3 - 1)))];
            FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
            //domain separation ds1
            int totalSubdomains = 64;
            if (decomposeModel) DdmCalculationsGeneral.BuildModelInterconnectionData(model);
            var decomposer = new AutomaticDomainDecomposer2(model, totalSubdomains);
            if (decomposeModel) decomposer.UpdateModel();
            var subdHexaIds = DdmCalculationsGeneral.DetermineHexaElementsSubdomainsFromModel(model);
            RveMatrixSubdomainInnerNodes = DdmCalculationsGeneral.DetermineRveSubdomainsInnerNodesFromModel(model);

            double volume = mp.L01 * mp.L02 * mp.L03;
            int hexaElementsNumber = model.ElementsDictionary.Count();

            //IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            //ds2
            int[] lowerCohesiveBound = new int[graphene_sheets_number]; int[] upperCohesiveBound = new int[graphene_sheets_number]; int[] grShElementssnumber = new int[graphene_sheets_number];
            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                FEMMeshBuilder.AddGrapheneSheet_with_o_x_Input_withRenumberingBondSlip(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                lowerCohesiveBound[j] = shellElementsNumber + element_counter_after_Adding_sheet + 1; //ds3
                upperCohesiveBound[j] = 2 * shellElementsNumber + element_counter_after_Adding_sheet;
                grShElementssnumber[j] = shellElementsNumber;
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }



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

            #region corner node data generation
                                         

            if (!decomposeModel)
            {
                (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model);
            }


            #endregion


            //ds4
            if (decomposeModel)
            {
                //int[][] subdCohElementIds = DdmCalculationsGeneral.DetermineCoheiveELementsSubdomainsSimple(model, totalSubdomains);
                (int[][] subdCohElementIds, Dictionary<int, List<int>> reassignedHexas, Dictionary<int, int> hexaOriginalSubdomains,
                    Dictionary<int, List<int>> SubdomainNeedsHexas) =
                    DdmCalculationsGeneral.DetermineCoheiveELementsSubdomainsSimple_Alte4(model, totalSubdomains, lowerCohesiveBound,
                    upperCohesiveBound, grShElementssnumber);

                (Dictionary<int, List<int>> hexaAndTheirSharingSubdomains, Dictionary<int, int> DuplicateHexaOriginalSubdomain) =
                   DdmCalculationsAlterna2.GetHexaSharingSubdomains(SubdomainNeedsHexas, subdHexaIds, false);
                int[][] SubdomainNeedsHexasIds = DdmCalculationsPartb.ConvertIntListToArray(SubdomainNeedsHexas, totalSubdomains);
                
                (int[][] subdAdditionalHexaIds, Dictionary<int, int> NewHexaIdsAndTheirOriginals) = DdmCalculationsAlterna2.ModelAddDuplicateHexas(model, mp, hexaAndTheirSharingSubdomains, totalSubdomains);

                //int[][] subdHexaIdsNew = DdmCalculationsGeneral.ReassignHexas(subdHexaIds, reassignedHexas, hexaOriginalSubdomains);
                int[][] subdShellElementIds = DdmCalculationsGeneral.DetermineShellELementsSubdomains(model, totalSubdomains, subdCohElementIds,
                lowerCohesiveBound, upperCohesiveBound, grShElementssnumber);
                int[][] subdElementIds1 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdHexaIds, subdCohElementIds);
                int[][] subdElementIds2 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdElementIds1, subdShellElementIds);
                int[][] subdElementIds3 = DdmCalculationsGeneral.CombineSubdomainElementsIdsArraysIntoOne(subdElementIds2, subdAdditionalHexaIds);
                DdmCalculationsGeneral.UndoModelInterconnectionDataBuild(model);
                DdmCalculations.SeparateSubdomains(model, subdElementIds3);

                model.ConnectDataStructures();
                bool isTrue = DdmCalculationsGeneral.CheckSubdomainsEmbeddingHostNodes(model, RveMatrixSubdomainInnerNodes);

                
                //PerfmormNesessaryChecksSubdomains(model);

                #region print extra data 

                (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes) = DefineCornerNodesFromCornerNodeData(CornerNodesData, model);

                #region find embedded
                EmbeddedNodes = new List<Node>();
                foreach (Element element in model.Elements)
                {
                    if (element.ID == 572)
                    {
                        string breakpoint = "here";
                    }

                    if (element.ElementType is IEmbeddedElement)
                    {
                        var element_I = element.ElementType as IEmbeddedElement;
                        foreach (EmbeddedNode embeddedNode in element_I.EmbeddedNodes)
                        {
                            if (!EmbeddedNodes.Contains(embeddedNode.Node))
                            { EmbeddedNodes.Add(embeddedNode.Node); }
                        }
                    }
                }
                #endregion

                DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);

                subdFreeBRNodes = new Dictionary<int, int[]>(); //nodeID, subdomainIDs  
                foreach (Node node in model.Nodes)
                {
                    if ((node.SubdomainsDictionary.Keys.Count() > 1 && !EmbeddedNodes.Contains(node)) && (!CornerNodesIds.Keys.Contains(node.ID) && (node.Constraints.Count() == 0)))
                    {
                        int[] subdIDs = node.SubdomainsDictionary.Keys.ToArray();
                        subdFreeBRNodes.Add(node.ID, subdIDs);
                    }
                }

                #endregion

                DdmCalculationsGeneral.UndoModelInterconnectionDataBuild(model);


                bool print_subdomain_data = true;
                if (print_subdomain_data)
                {
                    //DdmCalculationsGeneral.PrintSubdomainDataForPostPro(subdHexaIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro(subdHexaIds, SubdomainNeedsHexasIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdFreeBRNodes, subdomainOutputPath, @"\subdomainBRNodesAndSubd.txt");
                    DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndsubdomains, subdomainOutputPath, @"\CornerNodesAndSubdIds.txt");
                }

                bool get_subdomain_data = true;
                if (get_subdomain_data)
                {
                    (hexaPrint, cohePrint, shellPrint) = DdmCalculationsGeneral.GetSubdomainDataForPostPro(subdHexaIds, subdCohElementIds, subdShellElementIds, subdomainOutputPath);
                }

            }

            return new Tuple<Model, Dictionary<int, Node>, double>(model, boundaryNodes, volume);

        }

        private void PerfmormNesessaryChecksSubdomains(Model model)
        {
            int subdToCheck = 41;

            List<Node> TotalNodes = new List<Node>();
            List<Node> TotalNodes2 = new List<Node>();
            List<Node> Hexa8Nodes = new List<Node>();
            List<Node> Hexa8Nodes2 = new List<Node>();

            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements)
            {
                foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                {
                    TotalNodes.Add(node);
                }
            }
            TotalNodes = TotalNodes.Union(TotalNodes2).ToList();
            int ch01 = TotalNodes.Count();


            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements)
            {
                if (element.ElementType is Hexa8NonLinear)
                {
                    foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                    {
                        Hexa8Nodes.Add(node);
                    }
                }
            }
            foreach (Element element in model.SubdomainsDictionary[subdToCheck].Elements)
            {
                if (element.ElementType is Shell8NonLinear)
                {
                    foreach (Node node in element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element))
                    {
                        Hexa8Nodes2.Add(node);
                    }
                }
            }
            Hexa8Nodes = Hexa8Nodes.Union(Hexa8Nodes2).ToList();
            int ch020 = Hexa8Nodes.Count();

            var unconnected = TotalNodes.RemoveAll(x => Hexa8Nodes.Contains(x));

            var elements1 = TotalNodes[0].ElementsDictionary.Values.ToList();
            elements1.RemoveAll(x => (!(x.Subdomain.ID == 41)));
            elements1.RemoveAll(x => (!(x.ElementType is Hexa8NonLinear)));
            var elements2 = TotalNodes[1].ElementsDictionary.Values.ToList();
            elements2.RemoveAll(x => (!(x.Subdomain.ID == 41)));
            elements2.RemoveAll(x => (!(x.ElementType is Hexa8NonLinear)));

            var elemtns = model.ElementsDictionary.Values.ToList();
            elemtns.RemoveAll(x => (!CheckElementSameNodes(x, model.ElementsDictionary[4182])));

            var embElement2 = elemtns[2].ElementType as IEmbeddedElement;
            int[] subdIds2 = new int[embElement2.EmbeddedNodes.Count()];
            int[] nodesIDs2 = new int[embElement2.EmbeddedNodes.Count()];
            int[] elmntsIDs2 = new int[embElement2.EmbeddedNodes.Count()];
            for (int i1 = 0; i1 < embElement2.EmbeddedNodes.Count(); i1++)
            {
                subdIds2[i1] = model.ElementsDictionary[embElement2.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID].Subdomain.ID;
                //embElement2.EmbeddedNodes.ElementAt(0).EmbeddedInElement.Subdomain.ID dinei null apotelesma ara bug
                nodesIDs2[i1] = embElement2.EmbeddedNodes.ElementAt(i1).Node.ID;
                elmntsIDs2[i1] = embElement2.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID;
            }

            var embElement1 = elemtns[1].ElementType as IEmbeddedElement;
            int[] subdIds1 = new int[embElement1.EmbeddedNodes.Count()];
            int[] nodesIDs1 = new int[embElement1.EmbeddedNodes.Count()];
            int[] elmntsIDs1 = new int[embElement1.EmbeddedNodes.Count()];
            for (int i1 = 0; i1 < embElement1.EmbeddedNodes.Count(); i1++)
            {
                subdIds1[i1] = model.ElementsDictionary[embElement1.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID].Subdomain.ID;
                nodesIDs1[i1] = embElement1.EmbeddedNodes.ElementAt(i1).Node.ID;
                elmntsIDs1[i1] = embElement1.EmbeddedNodes.ElementAt(i1).EmbeddedInElement.ID;
            }

            List<int> elmnt1AssemblyNodes = elemtns[1].ElementType.DofEnumerator.GetNodesForMatrixAssembly(elemtns[1]).ToList().ConvertAll(x => x.ID);
            List<int> elmnt2AssemblyNodes = elemtns[2].ElementType.DofEnumerator.GetNodesForMatrixAssembly(elemtns[2]).ToList().ConvertAll(x => x.ID);

            bool ch03 = elmnt1AssemblyNodes.Contains(2643);
            //false
            bool ch04 = elmnt1AssemblyNodes.Contains(2666);
            //false
            bool ch05 = elmnt2AssemblyNodes.Contains(2666);
            //true
            bool ch06 = elmnt2AssemblyNodes.Contains(2643);
            //true


        }

        public bool CheckElementSameNodes(Element element1, Element element2 )
        {
            bool have8FirstSameNodes = true;
            for(int i1=0; i1<8;i1++)
            {
                if (!(element1.Nodes.ElementAt(i1).ID == element2.Nodes.ElementAt(i1).ID))
                {
                    have8FirstSameNodes = false;
                }
            }
            return have8FirstSameNodes;
        }

        // PROSOXH DEN ARKEI MONO TO PARAKATW NA GINEI UNCOMMENT WSTE NA GINEI IMPLEMENT TO IDegenerateRVEBuilder 
        //xreiazetai kai na xrhsimopoithei h katallhlh methodos tou femmeshbuilder gia to model and boundary nodes na dinei mono ta peripheral
        //public Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model model)
        //{
        //    return FEMMeshBuilder.GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
        //    //TODO:  Pithanws na epistrefetai apo GetModelAndBoundaryNodes ... AndConstraints.
        //}
        private void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        private Dictionary<int, HashSet<INode>> DefineCornerNodesPerSubdomainAndOtherwise(Dictionary<int, int[]> CornerNodesIdAndsubdomains, Model model)
        {
            Dictionary<int, HashSet<INode>> cornerNodesList = new Dictionary<int, HashSet<INode>>(model.Subdomains.Count());
            Dictionary<int, HashSet<INode>> cornerNodes = new Dictionary<int, HashSet<INode>>(model.Subdomains.Count());

            foreach (Subdomain subdomain in model.Subdomains)
            {
                cornerNodesList.Add(subdomain.ID, new HashSet<INode>());
            }

            foreach (int CornerNodeID in CornerNodesIdAndsubdomains.Keys)
            {
                Node node1 = model.NodesDictionary[CornerNodeID];
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            //foreach (Subdomain subdomain in model.Subdomains)
            //{
            //    cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID]);
            //}

            return cornerNodesList;
        }

        private (Dictionary<int, double[]> CornerNodesIds, Dictionary<int, int[]> CornerNodesIdAndsubdomains, Dictionary<int, HashSet<INode>> cornerNodes)
            DefineCornerNodesFromCornerNodeData(int[][] CornerNodesData, Model model)
        {
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            CornerNodesIds = new Dictionary<int, double[]>(CornerNodesData.Length); //nodeID, coordinate data
            CornerNodesIdAndsubdomains = new Dictionary<int, int[]>(CornerNodesData.Length);//nodeID, subdIds            
                                                                                            //var CornerNodesSubdomains = new Dictionary<int, int[]> (CornerNodesData.Length); //nodeID, subdomains opou anhkei
            for (int i1 = 0; i1 < CornerNodesData.Length; i1++)
            {
                int h1 = CornerNodesData[i1][0]; int h2 = CornerNodesData[i1][1]; int h3 = CornerNodesData[i1][2];
                int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                double nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                double nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                double nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                double[] coordinates = new double[3] { nodeCoordX, nodeCoordY, nodeCoordZ };
                CornerNodesIds.Add(nodeID, coordinates);
                CornerNodesIdAndsubdomains.Add(nodeID, model.NodesDictionary[nodeID].SubdomainsDictionary.Keys.ToArray());
            }

            cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(CornerNodesIdAndsubdomains, model);

            return (CornerNodesIds, CornerNodesIdAndsubdomains, cornerNodes);
        }

    }
}
