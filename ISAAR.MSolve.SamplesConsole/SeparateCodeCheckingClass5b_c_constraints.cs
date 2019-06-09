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

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass5b_c_constraints
    {
        public static double /*(double[], double[], double[], double[], IVector, IVector)*/ StiffnessMatrixOutputWrite()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelopc_constraints(1);

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            int[] hexaPrint = rveBuilder.hexaPrint;
            int[] cohePrint = rveBuilder.cohePrint;
            int[] shellPrint = rveBuilder.shellPrint;


            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
            // Solver
            var solverBuilder = new Feti1Solver.Builder(1e-4); //factorizationTolerance
            solverBuilder.ProblemIsHomogeneous = false;
            solverBuilder.PreconditionerFactory = new LumpedPreconditioner.Factory();
            solverBuilder.ProblemIsHomogeneous = true;
            Feti1Solver solver = solverBuilder.BuildSolver(model);

            //model.ConnectDataStructures(); // tha ginei entos tou rve builder             

            solver.OrderDofs(false);

            Dictionary<int, int[]> CornerNodesIdAndGlobalDofs = new Dictionary<int, int[]>( rveBuilder.CornerNodesIds.Keys.Count());//nodeID, globalDofs
            Dictionary<int, int[]> subdBRNodesAndGlobalDOfs = new Dictionary<int, int[]>(rveBuilder.subdFreeBRNodes.Keys.Count());//nodeID, globalDofs
            foreach(int corrnerNodeID in rveBuilder.CornerNodesIds.Keys )
            {
                Node CornerNode = model.NodesDictionary[corrnerNodeID];
                CornerNodesIdAndGlobalDofs.Add(corrnerNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[CornerNode, StructuralDof.TranslationZ]});

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(CornerNode, StructuralDof.RotationX, out int globalDofId);
                if(check)
                {
                    string breakpoint = "here";
                }
            }
            foreach (int boundaryNodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node boundaryNode = model.NodesDictionary[boundaryNodeID];

                bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(boundaryNode, StructuralDof.RotationX, out int globalDofId);
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
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndGlobalDofs, rveBuilder.subdomainOutputPath, @"\CornerNodesAndGlobalDofsIds.txt");
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdBRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\SubdBRNodesAndGlobalDofsIds.txt");

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
                    int globalDofId = model.GlobalDofOrdering.GlobalFreeDofs[node, doftype];
                    int[] localIds = new int[subdIds.Length];

                    for (int i1 = 0; i1 < subdIds.Length; i1++)
                    {
                        localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                    }
                    GlobalDofCoupledDataSubdIds.Add(globalDofId, subdIds);
                    GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId, localIds);
                }                
            }
            foreach (int nodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node node = model.NodesDictionary[nodeID];
                int[] subdIds = node.SubdomainsDictionary.Keys.ToArray();

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ,StructuralDof.RotationX,StructuralDof.RotationY };

                foreach (StructuralDof doftype in dofs)
                {
                    bool check = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, doftype, out int globalDofId);
                    if (check)
                    {
                        int[] localIds = new int[subdIds.Length];

                        for (int i1 = 0; i1 < subdIds.Length; i1++)
                        {
                            localIds[i1] = model.SubdomainsDictionary[subdIds[i1]].FreeDofOrdering.FreeDofs[node, doftype];
                        }
                        GlobalDofCoupledDataSubdIds.Add(globalDofId, subdIds);
                        GlobalDofCoupledDataLocalDofsInSubdIds.Add(globalDofId, localIds);
                    }
                }

            }

            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataSubdIds.txt");
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataLocalDofsInSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataLocalDofsInSubdIds.txt");
            #endregion


            foreach (ILinearSystem linearSystem in solver.LinearSystems.Values)
            {
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
            }

            solver.Initialize();

            Dictionary<int,IMatrix> subdomainMatrixes= solver.BuildGlobalMatrices(elementProvider); // to kanei o provide mesa
            

            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}Iter{1}Stiffness.txt";

            //var provider = new ProblemStructural(model, solver); //apo paradeigma Papagiannakis test 
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                //var subdMatrix= provider.CalculateMatrix(subdomain);

                string subdID = subdomain.ID.ToString();
                var subdMatrix = subdomainMatrixes[subdomain.ID];

                string counter_data = 1.ToString();
                string print_path = string.Format(print_path_gen, subdID, counter_data);

                var writer = new MatlabWriter();
                writer.WriteToFile((ISparseMatrix)subdMatrix,print_path,false);
                

            }

            string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}GlobalDofs.txt";
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                double[] subdomainGlobalDofs = new double[subdomain.FreeDofOrdering.NumFreeDofs];

                StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };

                foreach (Node node in subdomain.Nodes)
                {
                    foreach (StructuralDof dof in dofs)
                    {
                        bool check = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(node, dof, out int dofValue);
                        if(check)
                        {
                            subdomainGlobalDofs[dofValue] = model.GlobalDofOrdering.GlobalFreeDofs[node, dof];

                            if (model.GlobalDofOrdering.GlobalFreeDofs[node, dof]==0)
                            {
                                string breakpoint = "here";
                            }
                        }
                    }
                }

                string subdID = subdomain.ID.ToString();
                string print_path = string.Format(print_path_gen2, subdID);
                var writer = new MatlabWriter();
                writer.WriteToFile(Vector.CreateFromArray(subdomainGlobalDofs) , print_path, false);
            }



            return new double() ;
        }

        //change rve etc.
        public static void CheckSolution()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelop(1);

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            int[] hexaPrint = rveBuilder.hexaPrint;
            int[] cohePrint = rveBuilder.cohePrint;
            int[] shellPrint = rveBuilder.shellPrint;
            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();

            Dictionary<int, List<INode>> cornerNodesList = new Dictionary<int, List<INode>>(model.Subdomains.Count());
            Dictionary<int, INode[]> cornerNodes = new Dictionary<int, INode[]>(model.Subdomains.Count());

            foreach (Subdomain subdomain in model.Subdomains)
            {
                cornerNodesList.Add(subdomain.ID, new List<INode>());
            }

            foreach(int CornerNodeID in rveBuilder.CornerNodesIdAndsubdomains.Keys)
            {
                Node node1 = model.NodesDictionary[CornerNodeID];
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            foreach(Subdomain subdomain in model.Subdomains)
            {
                cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID].ToArray());
            }

            double load_value = 1;
            Load load1;            
            load1 = new Load()
            {
                Node =model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);
            

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-7;
            var fetiSolverBuilder = new FetiDPSolver.Builder(cornerNodes);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();
            FetiDPSolver fetiSolver = fetiSolverBuilder.BuildSolver(model);

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in fetiSolver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
        }

        public static void RunExample()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelopc_constraints(1);
            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[rveBuilder.CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);
            var cornerNodesAndSubds = rveBuilder.CornerNodesIdAndsubdomains;

            Dictionary<int, List<INode>> CornerNodes = new Dictionary<int, List<INode>>();

            foreach (int nodeId in cornerNodesAndSubds.Keys)
            {
                var subdIDs = cornerNodesAndSubds[nodeId];
                foreach (int subdId in subdIDs)
                {
                    if (CornerNodes.ContainsKey(subdId))
                    {
                        if (!CornerNodes[subdId].Contains(model.NodesDictionary[nodeId]))
                        {
                            CornerNodes[subdId].Add(model.NodesDictionary[nodeId]);
                        }
                    }
                    else
                    {
                        CornerNodes.Add(subdId, new List<INode>());
                        CornerNodes[subdId].Add(model.NodesDictionary[nodeId]);
                    }
                }
            }

            Dictionary<int, INode[]> cornerNodes = new Dictionary<int, INode[]>();
            foreach (int nodeId in CornerNodes.Keys)
            {
                cornerNodes[nodeId] = CornerNodes[nodeId].ToArray();
            }

            // Setup solver
            var interfaceSolverBuilder = new FetiDPInterfaceProblemSolver.Builder();
            interfaceSolverBuilder.MaxIterationsProvider = new PercentageMaxIterationsProvider(1);
            interfaceSolverBuilder.PcgConvergenceTolerance = 1E-10;
            var fetiSolverBuilder = new FetiDPSolver.Builder(cornerNodes);
            fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.ProblemIsHomogeneous = false;
            fetiSolverBuilder.PreconditionerFactory = new DirichletPreconditioner.Factory();
            FetiDPSolver fetiSolver = fetiSolverBuilder.BuildSolver(model);

            

            // Run the analysis
            var problem = new ProblemStructural(model, fetiSolver);
            var linearAnalyzer = new LinearAnalyzer(model, fetiSolver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, fetiSolver, problem, linearAnalyzer);
            staticAnalyzer.Initialize();
            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}Iter{1}Stiffness.txt";
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                //var subdMatrix= provider.CalculateMatrix(subdomain);

                string subdID = subdomain.ID.ToString();
                var subdMatrix = fetiSolver.LinearSystems[subdomain.ID].Matrix; // subdomainMatrixes[subdomain.ID];

                string counter_data = 1.ToString();
                string print_path = string.Format(print_path_gen, subdID, counter_data);

                var writer = new MatlabWriter();
                writer.WriteToFile((ISparseMatrix)subdMatrix, print_path, false);
            }
            string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}GlobalDofs.txt";
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
                writer.WriteToFile(Vector.CreateFromArray(subdomainGlobalDofs), print_path, false);
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
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(CornerNodesIdAndGlobalDofs, rveBuilder.subdomainOutputPath, @"\CornerNodesAndGlobalDofsIds.txt");
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(subdBRNodesAndGlobalDOfs, rveBuilder.subdomainOutputPath, @"\SubdBRNodesAndGlobalDofsIds.txt");

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

            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataSubdIds.txt");
            DdmCalculationsGeneral.PrintSubdomainDataForPostPro2(GlobalDofCoupledDataLocalDofsInSubdIds, rveBuilder.subdomainOutputPath, @"\GlobalDofCoupledDataLocalDofsInSubdIds.txt");
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

        }

        //needs to be corrected rve_multiple -> b
        public static void RunExampleSerial()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelopc_constraints(1);

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            var boundaryNodes = ModelAndNodes.Item2;

            #region corner nodes Data
            var RVE_id = 1;
            var rve_id_data = RVE_id.ToString();

            //renumbering_vector_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\REF_new_total_numbering.txt";
            var renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
            renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);
            int subdiscr1 = 2;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 8;
            int subdiscr1_shell = 6;
            int discr1_shell = 1;
            var mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            var mp = mpgp.Item1; //mp.hexa1 = 9; mp.hexa2 = 9; mp.hexa3 = 9;
            var gp = mpgp.Item2;
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));
            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            Dictionary<int, double[]> CornerNodesIds;
            Dictionary<int, int[]> CornerNodesIdAndsubdomains;
            int[][] CornerNodesData = new int[7][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            CornerNodesData[0] = new int[3] { 5 - 1, 5 - 1, 5 - 1 };

            // kai extra corner nodes
            CornerNodesData[1] = new int[3] { 2 - 1, 5 - 1, 5 - 1 };
            CornerNodesData[2] = new int[3] { 8 - 1, 5 - 1, 5 - 1 };
            CornerNodesData[3] = new int[3] { 5 - 1, 2 - 1, 5 - 1 };
            CornerNodesData[4] = new int[3] { 5 - 1, 8 - 1, 5 - 1 };
            CornerNodesData[5] = new int[3] { 5 - 1, 5 - 1, 2 - 1 };
            CornerNodesData[6] = new int[3] { 5 - 1, 5 - 1, 8 - 1 };


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
            #endregion 

            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);

            DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);


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
            int[] outputPositions = new int[model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs];
            for( int i1 = 0;  i1 < model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs ; i1++ )
            {

            }
            childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { 0 });

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion

            var globalU = solver.LinearSystems[0].Solution.CopyToArray();// fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
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
        }


        private static void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

    }
}
