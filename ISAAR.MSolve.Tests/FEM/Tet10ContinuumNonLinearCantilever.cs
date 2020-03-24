﻿using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using System.Linq;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class Tet10ContinuumNonLinearCantilever
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
            TotalDisplacementsPerIterationLog computedDisplacements = SolveModel();
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
        }

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, 
            TotalDisplacementsPerIterationLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[11]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

            expectedDisplacements[0] = new Dictionary<int, double> {
                { 0, 0.039075524153873623}, {11, -0.032541895181220408}, {23, -0.057387148941853101}, {35, -0.071994381984550326}, {47, -0.077053554770404833}
            };

            expectedDisplacements[0] = new Dictionary<int, double> {
    { 0,3.907552415387362300e-02 }, {11,-3.254189518122040800e-02 }, {23,-5.738714894185310100e-02 }, {35,-7.199438198455032600e-02 }, {47,-7.705355477040483300e-02 }};
            expectedDisplacements[1] = new Dictionary<int, double> {
    { 0,4.061313406968563400e-02 }, {11,-3.418876666892714500e-02 }, {23,-6.682708262609965400e-02 }, {35,-9.647418428408424700e-02 }, {47,-1.214556593711370000e-01 }};
            expectedDisplacements[2] = new Dictionary<int, double> {
    { 0,4.036171804663909300e-02 }, {11,-3.396515033613205900e-02 }, {23,-6.665084050819490600e-02 }, {35,-9.713633946904017000e-02 }, {47,-1.236631490430697600e-01 }};
            expectedDisplacements[3] = new Dictionary<int, double> {
    { 0,4.032905162001462800e-02 }, {11,-3.393260905426281900e-02 }, {23,-6.657423779424630200e-02 }, {35,-9.701032579889114200e-02 }, {47,-1.234941821043235900e-01 }};
            expectedDisplacements[4] = new Dictionary<int, double> {
    { 0,4.032900093364350700e-02 }, {11,-3.393255831972321500e-02 }, {23,-6.657411965268195100e-02 }, {35,-9.701012513482368300e-02 }, {47,-1.234939001150344400e-01 }};
            expectedDisplacements[5] = new Dictionary<int, double> {
    { 0,8.095088461395548400e-02 }, {11,-6.826589092291023000e-02 }, {23,-1.393261307096994000e-01 }, {35,-2.129883579558797000e-01 }, {47,-2.840192458274605800e-01 }};
            expectedDisplacements[6] = new Dictionary<int, double> {
    { 0,8.179065808895391600e-02 }, {11,-6.914910025670165100e-02 }, {23,-1.449912527358244700e-01 }, {35,-2.283048858573358000e-01 }, {47,-3.126785624370127000e-01 }};
            expectedDisplacements[7] = new Dictionary<int, double> {
    { 0,8.008398180684392400e-02 }, {11,-6.747544383562544000e-02 }, {23,-1.408463169597064000e-01 }, {35,-2.210877012127209200e-01 }, {47,-3.022981704019522300e-01 }};
            expectedDisplacements[8] = new Dictionary<int, double> {
    { 0,7.976397887674688300e-02 }, {11,-6.715673915988762400e-02 }, {23,-1.400151566610138300e-01 }, {35,-2.195056794855129700e-01 }, {47,-2.998365539162924900e-01 }};
            expectedDisplacements[9] = new Dictionary<int, double> {
    { 0,7.975945236918889600e-02 }, {11,-6.715223199537226400e-02 }, {23,-1.400036710136937400e-01 }, {35,-2.194845023343510200e-01 }, {47,-2.998046100841828000e-01 }};
            expectedDisplacements[10] = new Dictionary<int, double> {
    { 0,7.975944951878896600e-02 }, {11,-6.715222916021290600e-02 }, {23,-1.400036636464831200e-01 }, {35,-2.194844883932760600e-01 }, {47,-2.998045884933974200e-01 }};


            return expectedDisplacements;
        }

        private static TotalDisplacementsPerIterationLog SolveModel()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            BuildCantileverModel(model, 850);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            int increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-8;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Output
            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
            var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
            childAnalyzer.TotalDisplacementsPerIterationLog = log1;

            // Run the anlaysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return log1;
        }

        private static void BuildCantileverModel(Model model, double load_value)
        {
            //xrhsimopoiithike to  ParadeigmataElegxwnBuilder.HexaCantileverBuilder(Model model, double load_value)
            // allagh tou element kai tou material

            //ElasticMaterial3DTemp material1 = new ElasticMaterial3DTemp()
            //{
            //    YoungModulus = 1353000,
            //    PoissonRatio = 0.3,
            //};


            //VonMisesMaterial3D material1 = new VonMisesMaterial3D(1353000, 0.30, 1353000, 0.15);
            var material1 = new ElasticMaterial3D() { PoissonRatio = 0.3, YoungModulus = 1353000 };

            double[,] nodeData = new double[,] { {-0.2500000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,0.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,0.0000000000000000},
                {0.0000000000000000,0.2500000000000000,-1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,-1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,-1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,-1.0000000000000000},
                {0.0000000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,0.2500000000000000,0.0000000000000000},
                {0.2500000000000000,0.0000000000000000,0.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,0.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,0.0000000000000000},
                {-0.2500000000000000,0.2500000000000000,-0.5000000000000000},
                {0.2500000000000000,0.2500000000000000,-0.5000000000000000},
                {0.2500000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.2500000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.2500000000000000,0.2500000000000000,0.5000000000000000},
                {0.2500000000000000,0.2500000000000000,0.5000000000000000},
                {0.2500000000000000,-0.2500000000000000,0.5000000000000000},
                {-0.2500000000000000,-0.2500000000000000,0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,-1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,-1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,-1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,-1.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,-1.0000000000000000},
                {0.0000000000000000,0.0000000000000000,1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,1.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,1.0000000000000000},
                {0.0000000000000000,0.0000000000000000,0.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,0.0000000000000000},
                {0.1250000000000000,0.1250000000000000,0.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,0.0000000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.0000000000000000},
                {0.0000000000000000,0.2500000000000000,-0.5000000000000000},
                {0.1250000000000000,0.2500000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.2500000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.2500000000000000,-0.2500000000000000},
                {0.1250000000000000,0.2500000000000000,-0.2500000000000000},
                {0.2500000000000000,0.0000000000000000,-0.5000000000000000},
                {0.2500000000000000,-0.1250000000000000,-0.7500000000000000},
                {0.2500000000000000,0.1250000000000000,-0.7500000000000000},
                {0.2500000000000000,0.1250000000000000,-0.2500000000000000},
                {0.2500000000000000,-0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,-0.2500000000000000,-0.5000000000000000},
                {-0.1250000000000000,-0.2500000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.2500000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.2500000000000000,-0.2500000000000000},
                {-0.1250000000000000,-0.2500000000000000,-0.2500000000000000},
                {-0.2500000000000000,0.0000000000000000,-0.5000000000000000},
                {-0.2500000000000000,0.1250000000000000,-0.7500000000000000},
                {-0.2500000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.2500000000000000,-0.1250000000000000,-0.2500000000000000},
                {-0.2500000000000000,0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,0.2500000000000000,0.5000000000000000},
                {0.1250000000000000,0.2500000000000000,0.2500000000000000},
                {-0.1250000000000000,0.2500000000000000,0.2500000000000000},
                {-0.1250000000000000,0.2500000000000000,0.7500000000000000},
                {0.1250000000000000,0.2500000000000000,0.7500000000000000},
                {0.2500000000000000,0.0000000000000000,0.5000000000000000},
                {0.2500000000000000,-0.1250000000000000,0.2500000000000000},
                {0.2500000000000000,0.1250000000000000,0.2500000000000000},
                {0.2500000000000000,0.1250000000000000,0.7500000000000000},
                {0.2500000000000000,-0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,-0.2500000000000000,0.5000000000000000},
                {-0.1250000000000000,-0.2500000000000000,0.2500000000000000},
                {0.1250000000000000,-0.2500000000000000,0.2500000000000000},
                {0.1250000000000000,-0.2500000000000000,0.7500000000000000},
                {-0.1250000000000000,-0.2500000000000000,0.7500000000000000},
                {-0.2500000000000000,0.0000000000000000,0.5000000000000000},
                {-0.2500000000000000,0.1250000000000000,0.2500000000000000},
                {-0.2500000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.2500000000000000,-0.1250000000000000,0.7500000000000000},
                {-0.2500000000000000,0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,0.0000000000000000,-0.5000000000000000},
                {-0.1250000000000000,0.1250000000000000,-0.7500000000000000},
                {0.1250000000000000,0.1250000000000000,-0.7500000000000000},
                {0.1250000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.1250000000000000,-0.1250000000000000,-0.7500000000000000},
                {-0.1250000000000000,0.1250000000000000,-0.2500000000000000},
                {0.1250000000000000,0.1250000000000000,-0.2500000000000000},
                {0.1250000000000000,-0.1250000000000000,-0.2500000000000000},
                {-0.1250000000000000,-0.1250000000000000,-0.2500000000000000},
                {0.0000000000000000,0.0000000000000000,-0.7500000000000000},
                {0.0000000000000000,0.1250000000000000,-0.5000000000000000},
                {0.1250000000000000,0.0000000000000000,-0.5000000000000000},
                {0.0000000000000000,-0.1250000000000000,-0.5000000000000000},
                {-0.1250000000000000,0.0000000000000000,-0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,-0.2500000000000000},
                {0.0000000000000000,0.0000000000000000,0.5000000000000000},
                {-0.1250000000000000,0.1250000000000000,0.2500000000000000},
                {0.1250000000000000,0.1250000000000000,0.2500000000000000},
                {0.1250000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.2500000000000000},
                {-0.1250000000000000,0.1250000000000000,0.7500000000000000},
                {0.1250000000000000,0.1250000000000000,0.7500000000000000},
                {0.1250000000000000,-0.1250000000000000,0.7500000000000000},
                {-0.1250000000000000,-0.1250000000000000,0.7500000000000000},
                {0.0000000000000000,0.0000000000000000,0.2500000000000000},
                {0.0000000000000000,0.1250000000000000,0.5000000000000000},
                {0.1250000000000000,0.0000000000000000,0.5000000000000000},
                {0.0000000000000000,-0.1250000000000000,0.5000000000000000},
                {-0.1250000000000000,0.0000000000000000,0.5000000000000000},
                {0.0000000000000000,0.0000000000000000,0.7500000000000000},};

            int[] test = { 33, 2, 1, 88, 35, 13, 34, 90, 89, 97 };

            int[,] elementData = new int[,] {{33,2,1,88,35,13,34,90,89,97},
                {33,3,2,88,36,14,35,91,90,97},
                {33,4,3,88,37,15,36,92,91,97},
                {33,1,4,88,34,16,37,89,92,97},
                {48,1,2,88,50,13,49,89,90,98},
                {48,2,10,88,49,26,52,90,94,98},
                {48,10,9,88,52,21,51,94,93,98},
                {48,9,1,88,51,25,50,93,89,98},
                {53,2,3,88,55,14,54,90,91,99},
                {53,3,11,88,54,27,57,91,95,99},
                {53,11,10,88,57,22,56,95,94,99},
                {53,10,2,88,56,26,55,94,90,99},
                {58,3,4,88,60,15,59,91,92,100},
                {58,11,3,88,61,27,60,95,91,100},
                {58,12,11,88,62,23,61,96,95,100},
                {58,4,12,88,59,28,62,92,96,100},
                {63,4,1,88,65,16,64,92,89,101},
                {63,12,4,88,66,28,65,96,92,101},
                {63,9,12,88,67,24,66,93,96,101},
                {63,1,9,88,64,25,67,89,93,101},
                {43,9,10,88,44,21,45,93,94,102},
                {43,10,11,88,45,22,46,94,95,102},
                {43,11,12,88,46,23,47,95,96,102},
                {43,12,9,88,47,24,44,96,93,102},
                {43,10,9,103,45,21,44,105,104,112},
                {43,11,10,103,46,22,45,106,105,112},
                {43,12,11,103,47,23,46,107,106,112},
                {43,9,12,103,44,24,47,104,107,112},
                {68,9,10,103,70,21,69,104,105,113},
                {68,10,6,103,69,30,72,105,109,113},
                {68,6,5,103,72,17,71,109,108,113},
                {68,5,9,103,71,29,70,108,104,113},
                {73,10,11,103,75,22,74,105,106,114},
                {73,11,7,103,74,31,77,106,110,114},
                {73,7,6,103,77,18,76,110,109,114},
                {73,6,10,103,76,30,75,109,105,114},
                {78,11,12,103,80,23,79,106,107,115},
                {78,7,11,103,81,31,80,110,106,115},
                {78,8,7,103,82,19,81,111,110,115},
                {78,12,8,103,79,32,82,107,111,115},
                {83,12,9,103,85,24,84,107,104,116},
                {83,8,12,103,86,32,85,111,107,116},
                {83,5,8,103,87,20,86,108,111,116},
                {83,9,5,103,84,29,87,104,108,116},
                {38,5,6,103,39,17,40,108,109,117},
                {38,6,7,103,40,18,41,109,110,117},
                {38,7,8,103,41,19,42,110,111,117},
                {38,8,5,103,42,20,39,111,108,117},
                 };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y:  nodeData[nNode, 1], z: nodeData[nNode, 2] ));

            }

            // orismos elements 
            Element e1;
            int subdomainID = Tet10ContinuumNonLinearCantilever.subdomainID;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);
                //Dictionary<int,Node3D >
                List<Node> nodeSet = new List<Node>(8);
                for (int j = 0; j < 8; j++)
                {
                    int nodeID = elementData[nElement, j];
                    nodeSet.Add((Node)model.NodesDictionary[nodeID]);
                }

                var factory = new ContinuumElement3DFactory(material1, DynamicMaterial);

                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinear_v2(nodeSet, material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3), InterpolationHexa8Reverse_v2.UniqueInstance)// dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1.ID, e1);
            }

            int[] constrainedIds = new int[] { 1,2,3,4,13, 14, 15, 16, 33, 34, 35, 36, 37 };

            
            // constraint vashh opou z=-1
            for (int k1 = 1; k1 < constrainedIds.Length; k1++)
            {
                int k = constrainedIds[k1];
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
            }

            // fortish korufhs
            Load load1;
            for (int k = 5; k < 9; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = StructuralDof.TranslationX,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
    }

}