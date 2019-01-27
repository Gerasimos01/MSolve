//using ISAAR.MSolve.PreProcessor.Elements;
//using ISAAR.MSolve.PreProcessor.Materials;
using System.Collections.Generic;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;

namespace ISAAR.MSolve.SamplesConsole.DdmBenchmarks1
{
    public static class GrapheneCantileverBenchmarkFe2 
    {
        //PROELEFSI: IntegrationElasticCantileverBenchmark
        //allages: homogeneousRveBuilder-->Graphene

        public static void RunExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            int subdomainID = 1; model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });
            HexaCantileverBuilder_copyMS_222(model, 0.00219881744271988174427);

            model.ConnectDataStructures();

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);

            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            var solver = new SolverSkyline(linearSystems[subdomainID]);
            var linearSystemsArray = new[] { linearSystems[subdomainID] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);

            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 47 });
            var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
            childAnalyzer.IncrementalDisplacementsLog = log1;


            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;

            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            //return log1;
        }

        public static void HexaCantileverBuilder_copyMS_222(Model model, double load_value)
        {
            //PROELEFSI: ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222(Model model, double load_value)
            //H UPDATED PROELEFSI VRISKETAI PANW PANW STO ARXEIO

            var gpNeo = new grapheneSheetParameters()
            {
                //parametroi cohesive epifaneias
                T_o_3 = 0.12,//0.20,                   //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.15,//0.25,                   //0.5, // nm
                D_f_3 = 4,                      // nm
                T_o_1 = 0.12,//0.20,                   //0.05,// Gpa
                D_o_1 = 0.15,//0.25,                   //0.5, // nm
                D_f_1 = 4,                      // nm
                n_curve = 1.4
            };

            //VectorExtensions.AssignTotalAffinityCount();            
            IRVEbuilder RveBuilder1 = new GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostData(1);
            var material1 = new Microstructure3DevelopMultipleSubdomainsUseBaseSimuRand(RveBuilder1, false, 1);

            double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
            {0.250000,-0.250000,-1.000000},
            {-0.250000,0.250000,-1.000000},
            {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000},
            {0.250000,-0.250000,-0.500000},
            {-0.250000,0.250000,-0.500000},
            {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000},
            {0.250000,-0.250000,0.000000},
            {-0.250000,0.250000,0.000000},
            {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000},
            {0.250000,-0.250000,0.500000},
            {-0.250000,0.250000,0.500000},
            {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000},
            {0.250000,-0.250000,1.000000},
            {-0.250000,0.250000,1.000000},
            {0.250000,0.250000,1.000000}};

            int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14}, };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            }

            // constraint vashh opou z=-1
            for (int k = 1; k < 5; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // fortish korufhs
            Load load1;
            for (int k = 17; k < 21; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = DOFType.X,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
    }
}

