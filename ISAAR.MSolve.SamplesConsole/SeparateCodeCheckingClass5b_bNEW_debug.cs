using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.SamplesConsole.SupportiveClasses;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass5b_bNEW_debug
    {

        //prosthiki model.ConnectDataStructures entos rve gia na vrei to output node.Subdomains =/=0
        public static (Model, double[]) RunExample()
        {
            // EPILOGH RVE
            int subdiscr1 = 2;// 4;// 6;
            int discr1 = 2;// 3;//4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = discr1 * subdiscr1;// 23;
            int subdiscr1_shell = 6;//14;
            int discr1_shell = 1;
            int graphene_sheets_number = 0; //periektikothta 0.525% 

            double scale_factor = 1; //PROSOXH
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = SeperateIntegrationClassCheck.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90;
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;


            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3D(1, true, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);

            // EPILOGH RVE
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateDevelop(1, true); //edw ginetai develop h feti dp gia provlhmata 3d
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateLARGE(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop(1, true);

            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelopHSTAM(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopb(1, true);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbLARGE(1, true);
            //var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elem(1, true);

            bool WRITESTIFFNESSES = true;


            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            //model.ConnectDataStructures();

            double load_value = 1; //A.7
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);
            var cornerNodesAndSubds = rveBuilder.CornerNodesIdAndsubdomains;

            Dictionary<int, HashSet<INode>> cornerNodes = rveBuilder.cornerNodes;




            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1);
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-10;
            var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();
            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes);
            var fetiSolverBuilder = new FetiDPSolverPrint.Builder(cornerNodeSelection, fetiMatrices);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();
            FetiDPSolverPrint fetiSolver = fetiSolverBuilder.BuildSolver(model);



            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            // @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}Iter{1}Stiffness.txt";
            string print_path_gen = rveBuilder.subdomainOutputPath + @"\mat\Subdomain{0}Iter{1}Stiffness.txt";
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                //var subdMatrix= provider.CalculateMatrix(subdomain);

                string subdID = subdomain.ID.ToString();
                var subdMatrix = fetiSolver.LinearSystems[subdomain.ID].Matrix; // subdomainMatrixes[subdomain.ID];

                string counter_data = 1.ToString();
                string print_path = string.Format(print_path_gen, subdID, counter_data);

                var writer = new MatlabWriter();
                if (WRITESTIFFNESSES) writer.WriteToFile((ISparseMatrix)subdMatrix, print_path, false);
            }
            //string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}GlobalDofs.txt";
            string print_path_gen2 = rveBuilder.subdomainOutputPath + @"\mat\Subdomain{0}GlobalDofs.txt";
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                double[] subdomainGlobalDofs = new double[subdomain.FreeDofOrdering.NumFreeDofs];

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };

                foreach (Node node in subdomain.Nodes)
                {
                    foreach (StructuralDof dof in dofs)
                    {
                        bool check = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(node, dof, out int dofValue);
                        if (check)
                        {
                            subdomainGlobalDofs[dofValue] = model.GlobalDofOrdering.GlobalFreeDofs[node, dof];

                            if (model.GlobalDofOrdering.GlobalFreeDofs[node, dof] == 0)
                            {
                                string breakpoint = "here";
                            }
                        }
                    }
                }

                string subdID = subdomain.ID.ToString();
                string print_path = string.Format(print_path_gen2, subdID);
                var writer = new MatlabWriter();
                if (WRITESTIFFNESSES) writer.WriteToFile(Vector.CreateFromArray(subdomainGlobalDofs), print_path, false);
            }
            #region print 


            Dictionary<int, int[]> CornerNodesIdAndGlobalDofs = new Dictionary<int, int[]>(rveBuilder.CornerNodesIds.Keys.Count());//nodeID, globalDofs
            Dictionary<int, int[]> subdBRNodesAndGlobalDOfs = new Dictionary<int, int[]>(rveBuilder.subdFreeBRNodes.Keys.Count());//nodeID, globalDofs
            foreach (int corrnerNodeID in rveBuilder.CornerNodesIds.Keys)
            {
                Node CornerNode = model.NodesDictionary[corrnerNodeID];
                CornerNodesIdAndGlobalDofs.Add(corrnerNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationZ]});

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(CornerNode, StructuralDof.RotationX, out int globalDofId3);
                if (check)
                {
                    string breakpoint = "here";
                }
            }
            foreach (int boundaryNodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node boundaryNode = model.NodesDictionary[boundaryNodeID];

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(boundaryNode, StructuralDof.RotationX, out int globalDofId4);
                if (!check)
                {
                    subdBRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ]});
                }
                else
                {
                    subdBRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[5] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationY]});
                }
            }
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndGlobalDofs, rveBuilder.subdomainOutputPath, @"\CornerNodesAndGlobalDofsIds.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdBRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\SubdBRNodesAndGlobalDofsIds.txt");

            List<List<int>> extraConstraintsNoeds = rveBuilder.extraConstraintsNoeds;
            bool changeOrder = true;
            if (changeOrder)
            {
                extraConstraintsNoeds[0] = new List<int>() { 38 };
                extraConstraintsNoeds[1] = new List<int>() { 58 };
                extraConstraintsNoeds[2] = new List<int>() { 62 };
                extraConstraintsNoeds[3] = new List<int>() { 64 };
                extraConstraintsNoeds[4] = new List<int>() { 70 };
                extraConstraintsNoeds[5] = new List<int>() { 100 };
            }
            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = GetExtraConstrNodesPositions(subdBRNodesAndGlobalDOfs, extraConstraintsNoeds, model);

            #region  overwrite data model region

            PrintHexaModelData(model, rveBuilder.subdomainOutputPath);

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesIds = GetExtraConstrNodesIds(subdBRNodesAndGlobalDOfs, extraConstraintsNoeds, model);
            int[] brNodesMsolveWise = ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(rveBuilder.subdomainOutputPath + @"\RB_Nodes_IDs_MSOLVE_wise" + ".txt");
            Dictionary<int, int[]> subdBRNodesMsolveWiseAndGlobalDOfs = new Dictionary<int, int[]>();
            for (int i1 = 0; i1 < brNodesMsolveWise.Length; i1++) subdBRNodesMsolveWiseAndGlobalDOfs.Add(brNodesMsolveWise[i1], new int[1]);
            Dictionary<int, int[]> ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis = GetExtraConstrNodesPositions(subdBRNodesMsolveWiseAndGlobalDOfs, extraConstraintsNoeds, model);


            Dictionary<int, int[]> subd_msolve_BRNodesAndGlobalDOfs = new Dictionary<int, int[]>(rveBuilder.subdFreeBRNodes.Keys.Count());
            foreach (int boundaryNodeID in subdBRNodesMsolveWiseAndGlobalDOfs.Keys)
            {
                Node boundaryNode = model.NodesDictionary[boundaryNodeID];

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(boundaryNode, StructuralDof.RotationX, out int globalDofId4);
                if (!check)
                {
                    subd_msolve_BRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ]});
                }
                else
                {
                    subd_msolve_BRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[5] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.RotationY]});
                }
            }
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subd_msolve_BRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\subd_msolve_BRNodesAndGlobalDOfs.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBRNodesIds, rveBuilder.subdomainOutputPath, @"\ExtraConstrIdAndTheirBRNodesIds.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis, rveBuilder.subdomainOutputPath, @"\ExtraConstrIdAndTheirBR_msolveWise_NodesTheseis.txt");

            #endregion



            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(ExtraConstrIdAndTheirBRNodesTheseis, rveBuilder.subdomainOutputPath, @"\ExtraConstrIdAndTheirBRNodesTheseis.txt");
            
            #region coupled data arrays
            Dictionary<int, int[]> GlobalDofCoupledDataSubdIds = new Dictionary<int, int[]>(3 * (CornerNodesIdAndGlobalDofs.Count() + subdBRNodesAndGlobalDOfs.Count));
            Dictionary<int, int[]> GlobalDofCoupledDataLocalDofsInSubdIds = new Dictionary<int, int[]>(3 * (CornerNodesIdAndGlobalDofs.Count() + subdBRNodesAndGlobalDOfs.Count));
            foreach (int nodeID in rveBuilder.CornerNodesIds.Keys)
            {
                Node node = model.NodesDictionary[nodeID];
                int[] subdIds = node.SubdomainsDictionary.Keys.ToArray();

                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                foreach (StructuralDof doftype in dofs)
                {
                    int globalDofId1 = model.GlobalDofOrdering.GlobalFreeDofs[node, doftype];
                    int[] localIds = new int[subdIds.Length];

                    for (int i1 = 0; i1 < subdIds.Length; i1++)
                    {
                        localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                    }
                    GlobalDofCoupledDataSubdIds.Add(globalDofId1, subdIds);
                    GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId1, localIds);
                }
            }
            foreach (int nodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node node = model.NodesDictionary[nodeID];
                int[] subdIds = node.SubdomainsDictionary.Keys.ToArray();

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };

                foreach (StructuralDof doftype in dofs)
                {
                    bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, doftype, out int globalDofId2);
                    if (check)
                    {
                        int[] localIds = new int[subdIds.Length];

                        for (int i1 = 0; i1 < subdIds.Length; i1++)
                        {
                            localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                        }
                        GlobalDofCoupledDataSubdIds.Add(globalDofId2, subdIds);
                        GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId2, localIds);
                    }
                }

            }

            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataSubdIds.txt");
            if (WRITESTIFFNESSES) DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataLocalDofsInSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataLocalDofsInSubdIds.txt");
            #endregion
            #endregion

            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);

            Node monitoredNode = model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key];
            int globalDofId = model.GlobalDofOrdering.GlobalFreeDofs[monitoredNode, StructuralDof.TranslationZ];
            double solution = globalU[globalDofId];


            double[] uc = new double[3 * cornerNodesAndSubds.Count()];

            int node_counter = 0;
            foreach (int nodeId in cornerNodesAndSubds.Keys)
            {
                //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                Node node = model.NodesDictionary[nodeId];
                for (int i1 = 0; i1 < 3; i1++)
                {
                    int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
                    uc[3 * node_counter + i1] = globalU[globalDof];

                }
                node_counter++;
            }

            return (model, uc);

        }

        private static void PrintHexaModelData(Model model, string subdomainOutputPath)
        {
            Dictionary<int, List<int>> ElementIdsAndModelIds = new Dictionary<int, List<int>>();
            foreach (var element in model.Elements)
            {
                List<int> ElementNodesIds = element.Nodes.Select(x => x.ID).ToList();
                ElementIdsAndModelIds.Add(element.ID, ElementNodesIds);
            }

            int[] ElementIds = model.Elements.Select(x => x.ID).ToArray();
            int[] subdomainIds = model.Subdomains.Select(x => x.ID).ToArray();
            int[] NodeIds = model.Nodes.Select(x => x.ID).ToArray();

            int[,] ElementNodes = new int[ElementIds.GetLength(0), 8];
            int thesi = 0;
            foreach (var elementID in ElementIds)
            {
                int thesi2 = 0;
                foreach(Node node in model.ElementsDictionary[elementID].Nodes)
                {
                    ElementNodes[thesi, thesi2] = node.ID;
                    thesi2++;
                }
                thesi++;
            }

            double[,] NodeCoordinates = new double[model.Nodes.Count, 3];
            thesi = 0;
            foreach (var nodeID in NodeIds)
            {
                NodeCoordinates[thesi, 0] = model.NodesDictionary[nodeID].X1;
                NodeCoordinates[thesi, 1] = model.NodesDictionary[nodeID].X2;
                NodeCoordinates[thesi, 2] = model.NodesDictionary[nodeID].X3;
                thesi++;
            }

            int[,] SubdElements = new int[model.Subdomains.Count, model.Subdomains.ElementAt(0).Elements.Count];
            thesi = 0;
            foreach (var SubdId in subdomainIds)
            {
                int thesi2 = 0;
                foreach (var element in model.SubdomainsDictionary[SubdId].Elements)
                {
                    var elementID = element.ID;
                    SubdElements[thesi,thesi2]=elementID;
                    thesi2++;
                }
                thesi++;

            }

            List<int> constrainedNodes = new List<int>();
            foreach (var constraint in model.Constraints)
            {
                if (!(constrainedNodes.Contains(constraint.row.ID)))
                { constrainedNodes.Add(constraint.row.ID); }
            }
            var constraintIds = constrainedNodes.ToArray();

            //print model reconstruction data 
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileVectorMsolveInput(ElementIds,  subdomainOutputPath +@"\MsolveModel\"+ @"\ElementIds.txt");
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileVectorMsolveInput(subdomainIds, subdomainOutputPath + @"\MsolveModel\" + @"\subdomainIds.txt");
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileVectorMsolveInput(NodeIds, subdomainOutputPath + @"\MsolveModel\" + @"\NodeIds.txt");
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileVectorMsolveInput(constraintIds, subdomainOutputPath + @"\MsolveModel\" + @"\constraintIds.txt");



            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileMsolveInput(ElementNodes, subdomainOutputPath + @"\MsolveModel\" + @"\ElementNodes.txt");
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileMsolveInput(NodeCoordinates, subdomainOutputPath + @"\MsolveModel\" + @"\NodeCoordinates.txt");
            ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.WriteToFileMsolveInput(SubdElements, subdomainOutputPath + @"\MsolveModel\" + @"\SubdElements.txt");

        }

        
    

        private static Dictionary<int, int[]> GetExtraConstrNodesPositions(Dictionary<int, int[]> subdBRNodesAndGlobalDOfs, List<List<int>> extraConstraintsNoedsIds, Model model)
        {
            int[] positionsOfBRNodes = new int[subdBRNodesAndGlobalDOfs.Keys.Max()];

            for (int i1 = 0; i1 < subdBRNodesAndGlobalDOfs.Keys.Count(); i1++)
            {
                int BRnodeID = subdBRNodesAndGlobalDOfs.ElementAt(i1).Key;
                positionsOfBRNodes[BRnodeID - 1] = i1 + 1;
            }

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = new Dictionary<int, int[]>();

            for (int i1 = 0; i1 < extraConstraintsNoedsIds.Count(); i1++)
            {
                int[] BRnodesTheseis = new int[extraConstraintsNoedsIds.ElementAt(i1).Count()];
                for (int i2 = 0; i2 < extraConstraintsNoedsIds.ElementAt(i1).Count(); i2++)
                {
                    BRnodesTheseis[i2] = positionsOfBRNodes[extraConstraintsNoedsIds.ElementAt(i1).ElementAt(i2) - 1];
                }

                ExtraConstrIdAndTheirBRNodesTheseis.Add(i1 + 1, BRnodesTheseis);
            }

            return ExtraConstrIdAndTheirBRNodesTheseis;
        }

        private static Dictionary<int, int[]> GetExtraConstrNodesIds(Dictionary<int, int[]> subdBRNodesAndGlobalDOfs, List<List<int>> extraConstraintsNoedsIds, Model model)
        {
            int[] positionsOfBRNodes = new int[subdBRNodesAndGlobalDOfs.Keys.Max()];

            for (int i1 = 0; i1 < subdBRNodesAndGlobalDOfs.Keys.Count(); i1++)
            {
                int BRnodeID = subdBRNodesAndGlobalDOfs.ElementAt(i1).Key;
                positionsOfBRNodes[BRnodeID - 1] = i1 + 1;
            }

            Dictionary<int, int[]> ExtraConstrIdAndTheirBRNodesTheseis = new Dictionary<int, int[]>();

            for (int i1 = 0; i1 < extraConstraintsNoedsIds.Count(); i1++)
            {
                int[] BRnodesTheseis = new int[extraConstraintsNoedsIds.ElementAt(i1).Count()];
                for (int i2 = 0; i2 < extraConstraintsNoedsIds.ElementAt(i1).Count(); i2++)
                {
                    BRnodesTheseis[i2] = extraConstraintsNoedsIds.ElementAt(i1).ElementAt(i2);
                }

                ExtraConstrIdAndTheirBRNodesTheseis.Add(i1 + 1, BRnodesTheseis);
            }

            return ExtraConstrIdAndTheirBRNodesTheseis;
        }
        //needs to be corrected rve_multiple -> b kai to path kai ta stoixeia diakritopoihshs pou einai afhmena exwterika (Genika elegxoume connectDataStructures kai defineAppropriateConstraintsForBoundaryNodes)
        public static (Model, double[]) RunExampleSerial()
        {
            // EPILOGH RVE
            int subdiscr1 = 2;// 4;// 6;
            int discr1 = 2;// 3;//4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = discr1 * subdiscr1;// 23;
            int subdiscr1_shell = 6;//14;
            int discr1_shell = 1;
            int graphene_sheets_number = 0; //periektikothta 0.525% 

            double scale_factor = 1;
            //tvra ginontai scale input tou mpgp = getRe... methodou
            graphene_sheets_number = (int)Math.Floor(scale_factor * scale_factor * scale_factor * graphene_sheets_number);
            subdiscr1 = (int)Math.Floor(scale_factor * subdiscr1);


            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp = SeperateIntegrationClassCheck.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //mpgp.Item2.E_shell = 0.0000001;
            mpgp.Item1.L01 = scale_factor * 90; mpgp.Item1.L02 = scale_factor * 90; mpgp.Item1.L03 = scale_factor * 90;
            mpgp.Item1.L01 = scale_factor * mpgp.Item1.L01; mpgp.Item1.L02 = scale_factor * mpgp.Item1.L02; mpgp.Item1.L03 = scale_factor * mpgp.Item1.L03;


            var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop3D(1, false, mpgp,
            subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell, graphene_sheets_number);

            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicateLARGE(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2c_alte(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbDuplicate_2d_alteDevelop(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopb(1, false);
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbLARGE(1, false); // diorthose kai to parakatw path apla gia na mhn xtupaei.
            //var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elem(1, false); //A.1

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            var boundaryNodes = ModelAndNodes.Item2;

            var renumbering_vector_path = rveBuilder.renumbering_vector_path;

            
            //var mpgp = rveBuilder.mpgp; 
            var mp = mpgp.Item1; 
            var gp = mpgp.Item2;
            renumbering renumbering = new renumbering(ISAAR.MSolve.SamplesConsole.SupportiveClasses.PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            Dictionary<int, double[]> CornerNodesIds;
            Dictionary<int, int[]> CornerNodesIdAndsubdomains;
            int[][] CornerNodesData = rveBuilder.CornerNodesData;
            CornerNodesIds = rveBuilder.CornerNodesIds;
            CornerNodesIdAndsubdomains = rveBuilder.CornerNodesIdAndsubdomains;

            //A.6
            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);

            //A.7
            bool boundaryNodesOrPaktwsh=true;
            if (boundaryNodesOrPaktwsh)
            {
                DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);
            }
            else
            {
                int[][] paktwsiNodesData = new int[4][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
                int thesi = 0;
                int j1 = 0;
                for (int i2 = 0; i2 < 2; i2++)
                {
                    for (int i3 = 0; i3 < 2; i3++)
                    {
                        paktwsiNodesData[thesi] = new int[3] { j1, i2 * 1, i3 * 1 };
                        thesi++;
                    }
                }
                Dictionary<int, Node> paktwshAristeraNodes = new Dictionary<int, Node>();
                for (int i1 = 0; i1 < paktwsiNodesData.Length; i1++)
                {
                    int h1 = paktwsiNodesData[i1][0]; int h2 = paktwsiNodesData[i1][1]; int h3 = paktwsiNodesData[i1][2];
                    int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                    paktwshAristeraNodes.Add(nodeID, model.NodesDictionary[nodeID]);


                }
                DefineAppropriateConstraintsForBoundaryNodes(model, paktwshAristeraNodes);
            }


            #region define solver from Quad4LinearCantilever %81 and IntegrationElastic... %37 tests.
            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            // Request output
            //int[] outputPositions = new int[model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs];
            //for( int i1 = 0;  i1 < model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs ; i1++ )
            //{

            //}
            //childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { 0 });

            // Run the anlaysis - part a
            parentAnalyzer.Initialize();

            //A.8    
            bool print_data = false;
            if (print_data)
            {
                printElementStiffnessAndData(model);
            }

            // Run the anlaysis - part b
            parentAnalyzer.Solve();
            #endregion

            var globalU = solver.LinearSystems[1].Solution.CopyToArray();// fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
            double[] uc = new double[3 * CornerNodesIds.Count()];

            int node_counter = 0;
            foreach (int nodeId in CornerNodesIds.Keys)
            {
                //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                Node node = model.NodesDictionary[nodeId];
                for (int i1 = 0; i1 < 3; i1++)
                {
                    int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
                    uc[3 * node_counter + i1] = globalU[globalDof];

                }
                node_counter++;
            }

            return (model, uc);
        }

        public static void printElementStiffnessAndData(Model model)
        {
            var dofOrdering = model.SubdomainsDictionary[0].FreeDofOrdering;
            var matrixProvider = new ElementStructuralStiffnessProvider();
            
            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}Stiffness.txt";
            string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}elementDofIndices.txt";
            string print_path_gen3 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}subdomainDofIndices.txt";
            foreach (Element element in model.ElementsDictionary.Values)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                var elementMatrix = matrixProvider.Matrix(element).CopytoArray2D();
                string counterData = element.ID.ToString();
                string print_path = string.Format(print_path_gen, counterData);
                var writer = new Array2DWriter();
                writer.WriteToFile(elementMatrix, print_path, false);

                //Aray1Dwriter
                string print_path2 = string.Format(print_path_gen2, counterData);
                var writer2 = new MatlabWriter();
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(elementDofIndices)), print_path2, false);

                string print_path3 = string.Format(print_path_gen3, counterData);
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(subdomainDofIndices)), print_path3, false);
            }
        }

        private static void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        public static double[] Cnvrt(int[] array)
        {
            double[] array2 = new double[array.Length];
            for(int i1=0; i1<array.Length; i1++)
            {
                array2[i1] = (double)array[i1];
            }
            return array2;
        }

    }
}
