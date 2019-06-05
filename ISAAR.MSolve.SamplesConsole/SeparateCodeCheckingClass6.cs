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

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass6
    {
        public static double /*(double[], double[], double[], double[], IVector, IVector)*/ StiffnessMatrixOutputWrite()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelopHexaOnly(1);

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
            }
            foreach (int boundaryNodeID in rveBuilder.subdFreeBRNodes.Keys)
            {
                Node boundaryNode = model.NodesDictionary[boundaryNodeID];
                subdBRNodesAndGlobalDOfs.Add(boundaryNodeID, new int[3] { model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationX],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationY],
                                                                           model.GlobalDofOrdering.GlobalFreeDofs[boundaryNode, StructuralDof.TranslationZ]});
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
    }
}
