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
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;
using ISAAR.MSolve.Analyzers;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass_c_alte_5_elem
    {
        public static (Model, double[]) RunExampleSerial()
        {
            var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elem(1, false); // diorthose kai to parakatw path apla gia na mhn xtupaei.

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            var boundaryNodes = ModelAndNodes.Item2;

            
            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[model.NodesDictionary.ElementAt(10).Value.ID],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);


            
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

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion

            //var globalU = solver.LinearSystems[1].Solution.CopyToArray();// fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
            //double[] uc = new double[3 * CornerNodesIds.Count()];

            //int node_counter = 0;
            //foreach (int nodeId in CornerNodesIds.Keys)
            //{
            //    //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
            //    StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

            //    Node node = model.NodesDictionary[nodeId];
            //    for (int i1 = 0; i1 < 3; i1++)
            //    {
            //        int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
            //        uc[3 * node_counter + i1] = globalU[globalDof];

            //    }
            //    node_counter++;
            //}



            double[] uc = new double[2];
            return (model, uc);
        }

    }

}

