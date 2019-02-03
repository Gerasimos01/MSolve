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
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Materials;

namespace ISAAR.MSolve.SamplesConsole.IntegrationTests2
{
    public static class NRNLAnalyzerDevelopTest
    {
        //PROELEFSI: programElegxoiDdm opou eixan ginei comment out kai den htan updated apo ekdosh feat/prosthiki_allagwn 

        public static void SolveDisplLoadsExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice = 1;
            if (model__builder_choice == 1) // 
            { HexaCantileverBuilderDispControl(model, 850); }            
            
            model.ConnectDataStructures();            

            Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, Vector>();
            uInitialFreeDOFDisplacementsPerSubdomain.Add(model.SubdomainsDictionary[1].ID, new Vector(model.SubdomainsDictionary[1].TotalDOFs));// prosoxh sto Id twn subdomain
            Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();
            for (int k = 17; k < 21; k++)
            {
                boundaryNodes.Add(model.NodesDictionary[k].ID, model.NodesDictionary[k]);
            }
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            Dictionary<DOFType, double> initialConvergedBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
            initialConvergedBoundaryNodalDisplacements.Add(DOFType.X, 0);
            for (int k = 17; k < 21; k++)
            {
                initialConvergedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, initialConvergedBoundaryNodalDisplacements);
            }
            Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            double[] prescribedDisplacmentXValues = new double[4] { 7.81614E-01, 7.07355E-01, 7.81614E-01, 7.07355E-01 };
            for (int k = 17; k < 21; k++)
            {
                Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
                totalBoundaryNodalDisplacements.Add(DOFType.X, prescribedDisplacmentXValues[k - 17]);
                totalBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
            }

            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
            Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));

            //comment section 1 palaia version
            //SolverSkyline solver = new SolverSkyline(model);
            //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces); // elegxos me model.subdomainsDictionary[1]

            ProblemStructural provider = new ProblemStructural(model, linearSystems);


            // PARADEIGMA A: LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //SolverSkyline2 solver = new SolverSkyline2(linearSystems[1]); //H MARIA XRHSIMOPOIEI TON sklinesolver 
            //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            //---------------------------------------------------------------------------------------------------------------------------------

            // PARADEIGMA B: Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //PALIA DIATUPWSH: NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48); 
            // NEA DIATUPWSH:
            var solver = new SolverSkyline(linearSystems[1]);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelop(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, equivalentContributionsAssemblers);
            //h epomenhgrammh einai gia paradeigma ws pros to access
            //IAnalyzer childAnalyzer2 = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);


            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            //---------------------------------------------------------------------------------------------------------------------------------


            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);           

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


        }

        public static void HexaCantileverBuilderDispControl(Model model, double load_value)
        {
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 1353000,
                PoissonRatio = 0.3,
            };

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
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
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

            // thetoume constraint tous prescribed
            //Load load1;
            for (int k = 17; k < 21; k++)
            {
                //load1 = new Load()
                //{
                //    Node = model.NodesDictionary[k],
                //    DOF = DOFType.X,
                //    Amount = 1 * load_value
                //};
                //model.Loads.Add(load1);
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
            }
        }
    }
}
