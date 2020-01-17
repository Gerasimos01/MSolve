﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Logging;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP.IntegrationTests
{
    public static class Cantilever3DFetiDPTests
    {
        private const string plotPath = @"C:\Users\Serafeim\Desktop\FETI-DP\Plots";

        private const int singleSubdomainID = 0;

        public enum Crosspoints { Minimum, FullyRedundant }
        public enum Corners { C4848, C0404, C2626, C0426, C1526 }
        public enum Augmented { None, Nodes21 }


        [Theory]
        [InlineData(Corners.C4848, Crosspoints.Minimum, Augmented.None)]
        [InlineData(Corners.C2626, Crosspoints.Minimum, Augmented.None)]
        [InlineData(Corners.C2626, Crosspoints.FullyRedundant, Augmented.None)]
        [InlineData(Corners.C2626, Crosspoints.Minimum, Augmented.Nodes21)]
        [InlineData(Corners.C2626, Crosspoints.FullyRedundant, Augmented.Nodes21)]
        //[InlineData(Corners.C1526, Crosspoints.Minimum, Augmented.None)]
        //[InlineData(Corners.C1526, Crosspoints.FullyRedundant, Augmented.None)]
        // Each subdomain must have at least 1 corner node. Otherwise some matrices degenerate
        //[InlineData(Corners.C0426, Crosspoints.Minimum, Augmented.None)]          
        //[InlineData(Corners.C0426, Crosspoints.FullyRedundant, Augmented.None)]
        //[InlineData(Corners.C0404, Crosspoints.Minimum)]
        //[InlineData(Corners.C0404, Crosspoints.FullyRedundant)]
        public static void Run(Corners corners, Crosspoints crosspoints, Augmented augmented)
        {
            double pcgConvergenceTol = 1E-5;
            IVectorView directDisplacements = SolveModelWithoutSubdomains();
            (IVectorView ddDisplacements, ISolverLogger logger) =
                SolveModelWithSubdomains(corners, crosspoints, augmented, pcgConvergenceTol);
            double normalizedError = directDisplacements.Subtract(ddDisplacements).Norm2() / directDisplacements.Norm2();

            // The error is provided in the reference solution the, but it is almost impossible for two different codes run on 
            // different machines to achieve the exact same accuracy.
            Assert.Equal(0.0, normalizedError, 6);
            //Assert.True(directDisplacements.Equals(ddDisplacements, 1E-5));
        }

        internal static Model CreateModel()
        {
            // Subdomains:

            double E0 = 2.1E7;
            //double E1 = stiffnessRatio * E0;

            var builder = new Uniform3DModelBuilder();
            builder.DomainLengthX = 4.0;
            builder.DomainLengthY = 2.0;
            builder.DomainLengthZ = 1.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumSubdomainsZ = 1;
            builder.NumTotalElementsX = 4;
            builder.NumTotalElementsY = 2;
            builder.NumTotalElementsZ = 1;
            builder.YoungModulus = E0;
            //builder.YoungModuliOfSubdomains = new double[,] { { E1, E0, E0, E0 }, { E1, E0, E0, E0 } };

            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationY, 0.0);
            builder.PrescribeDisplacement(Uniform3DModelBuilder.BoundaryRegion.MinX, StructuralDof.TranslationZ, 0.0);
            builder.DistributeLoadAtNodes(Uniform3DModelBuilder.BoundaryRegion.MaxX, StructuralDof.TranslationY, -600.0);

            return builder.BuildModel();
        }

        internal static IVectorView SolveModelWithoutSubdomains()
        {
            Model model = CreateSingleSubdomainModel();

            // Solver
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[singleSubdomainID].Solution;
        }

        private static Model CreateSingleSubdomainModel()
        {
            // Replace the existing subdomains with a single one 
            Model model = CreateModel();
            model.SubdomainsDictionary.Clear();
            var subdomain = new Subdomain(singleSubdomainID);
            model.SubdomainsDictionary.Add(singleSubdomainID, subdomain);
            foreach (Element element in model.ElementsDictionary.Values) subdomain.Elements.Add(element.ID, element);
            return model;
        }

        private static (IVectorView globalDisplacements, ISolverLogger logger) SolveModelWithSubdomains(
            Corners corners, Crosspoints crosspoints, Augmented augmented, double pcgConvergenceTolerance)
        {
            // Model
            Model model = CreateModel();
            model.ConnectDataStructures();

            // Corner, midside nodes
            ICornerNodeSelection cornerNodeSelection = DefineCornerNodes(model, corners);
            IMidsideNodesSelection midsideNodes = DefineMidsideNodes(model, augmented);

            // Crosspoint strategy
            ICrosspointStrategy crosspointStrategy = null;
            if (crosspoints == Crosspoints.FullyRedundant) crosspointStrategy = new FullyRedundantConstraints();
            else if (crosspoints == Crosspoints.Minimum) crosspointStrategy = new MinimumConstraints();
            else throw new ArgumentException();

            //// Plot for debugging
            //var logger = new DomainDecompositionLoggerFetiDP(cornerNodeSelection, plotPath, true);
            //logger.PlotSubdomains(model);

            // Specify PCG settings
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = pcgConvergenceTolerance,
                MaxIterationsProvider = new FixedMaxIterationsProvider(100)
            };

            // Solver
            if (augmented == Augmented.None)
            {
                var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
                var solverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);
                solverBuilder.ProblemIsHomogeneous = true;
                solverBuilder.Preconditioning = new DirichletPreconditioning();
                solverBuilder.CrosspointStrategy = crosspointStrategy;
                solverBuilder.PcgSettings = pcgSettings;
                FetiDPSolverSerial fetiSolver = solverBuilder.Build(model, cornerNodeSelection);

                // Run the analysis
                RunAnalysis(model, fetiSolver); //check dof separator

                // Gather the global displacements
                Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();
                return (globalDisplacements, fetiSolver.Logger);
            }
            else if (augmented == Augmented.Nodes21)
            {

                var fetiMatrices = new FetiDP3dMatrixManagerFactoryDense();
                var solverBuilder = new FetiDP3dSolverSerial.Builder(fetiMatrices);
                solverBuilder.ProblemIsHomogeneous = true;
                solverBuilder.Preconditioning = new DirichletPreconditioning();
                solverBuilder.CrosspointStrategy = crosspointStrategy;
                solverBuilder.PcgSettings = pcgSettings;
                FetiDP3dSolverSerial fetiSolver = solverBuilder.Build(model, cornerNodeSelection, midsideNodes);

                // Run the analysis
                RunAnalysis(model, fetiSolver); //check dof separator

                // Gather the global displacements
                Vector globalDisplacements = fetiSolver.GatherGlobalDisplacements();
                return (globalDisplacements, fetiSolver.Logger);
            }
            else throw new ArgumentException();
        }

        private static ICornerNodeSelection DefineCornerNodes(Model model, Corners corners)
        {
            // Important nodes
            double meshTol = 1E-6;
            Node n200 = FindNode(2, 0, 0, model, meshTol);
            Node n210 = FindNode(2, 1, 0, model, meshTol);
            Node n220 = FindNode(2, 2, 0, model, meshTol);
            Node n201 = FindNode(2, 0, 1, model, meshTol);
            Node n211 = FindNode(2, 1, 1, model, meshTol);
            Node n221 = FindNode(2, 2, 1, model, meshTol);
            IEnumerable<Node> nodes2 = FindNodesWithX(2, model, meshTol);

            Node n400 = FindNode(4, 0, 0, model, meshTol);
            Node n410 = FindNode(4, 1, 0, model, meshTol);
            Node n420 = FindNode(4, 2, 0, model, meshTol);
            Node n401 = FindNode(4, 0, 1, model, meshTol);
            Node n411 = FindNode(4, 1, 1, model, meshTol);
            Node n421 = FindNode(4, 2, 1, model, meshTol);
            IEnumerable<Node> nodes4 = FindNodesWithX(4, model, meshTol);

            // Corner nodes
            var cornerNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                cornerNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            if (corners == Corners.C4848)
            {
                foreach (Node node in nodes4.Union(nodes2))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C0404)
            {
                foreach (Node node in nodes4)
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C2626)
            {
                var extraCorners = new Node[] { n200, n201, n220, n221 };
                foreach (Node node in nodes4.Union(extraCorners))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C0426)
            {
                var extraCorners = new Node[] { n220, n221 };
                foreach (Node node in nodes4.Union(extraCorners))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else if (corners == Corners.C1526)
            {
                var extraCorners = new Node[] { n201, n220, n221 };
                foreach (Node node in nodes4.Union(extraCorners))
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        cornerNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else throw new ArgumentException();
            return new UsedDefinedCornerNodes(cornerNodesOfEachSubdomain);
        }

        private static IMidsideNodesSelection DefineMidsideNodes(Model model, Augmented augmented)
        {
            // Important nodes
            double meshTol = 1E-6;
            Node n200 = FindNode(2, 0, 0, model, meshTol);
            Node n210 = FindNode(2, 1, 0, model, meshTol);
            Node n220 = FindNode(2, 2, 0, model, meshTol);
            Node n201 = FindNode(2, 0, 1, model, meshTol);
            Node n211 = FindNode(2, 1, 1, model, meshTol);
            Node n221 = FindNode(2, 2, 1, model, meshTol);
            IEnumerable<Node> nodes2 = FindNodesWithX(2, model, meshTol);

            Node n400 = FindNode(4, 0, 0, model, meshTol);
            Node n410 = FindNode(4, 1, 0, model, meshTol);
            Node n420 = FindNode(4, 2, 0, model, meshTol);
            Node n401 = FindNode(4, 0, 1, model, meshTol);
            Node n411 = FindNode(4, 1, 1, model, meshTol);
            Node n421 = FindNode(4, 2, 1, model, meshTol);
            IEnumerable<Node> nodes4 = FindNodesWithX(4, model, meshTol);

            // Midside nodes
            var midsideNodesOfEachSubdomain = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                midsideNodesOfEachSubdomain[subdomain] = new HashSet<INode>();
            }
            if (augmented == Augmented.None) return null;
            else if (augmented == Augmented.Nodes21)
            {
                var midsideNodes = new Node[] { n210, n211 };
                foreach (Node node in midsideNodes)
                {
                    foreach (Subdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        midsideNodesOfEachSubdomain[subdomain].Add(node);
                    }
                }
            }
            else throw new ArgumentException();
            return new UserDefinedMidsideNodes(midsideNodesOfEachSubdomain,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
        }

        private static Node FindNode(double x, double y, double z, Model model, double tol)
        {
            return model.NodesDictionary.Values.Where(
                n => (Math.Abs(n.X - x) <= tol) && (Math.Abs(n.Y - y) <= tol) && (Math.Abs(n.Z - z) <= tol)).First();
        }

        private static IEnumerable<Node> FindNodesWithX(double x, Model model, double tol)
        {
            return model.NodesDictionary.Values.Where(n => (Math.Abs(n.X - x) <= tol));
        }

        private static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis
            solver.OrderDofs(false);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ILinearSystem linearSystem = solver.GetLinearSystem(subdomain);
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
            solver.BuildGlobalMatrix(new ElementStructuralStiffnessProvider());
            model.ApplyLoads();
            LoadingUtilities.ApplyNodalLoads(model, solver);
            solver.Solve();
        }
    }
}
