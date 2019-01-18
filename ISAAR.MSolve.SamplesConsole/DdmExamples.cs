using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.SamplesConsole.SupportiveClasses;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.SamplesConsole;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;

namespace ISAAR.MSolve.SamplesConsole
{
    class DdmExamples
    {
        public static void SolveRVEExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice =208;   // 9 einai to megalo me to renumbering pou tsekaretai
            
            if (model__builder_choice == 208) // 
            { //RVEkanoninkhsGewmetriasBuilder.Reference2RVEExample50_000withRenumberingwithInput(model);

                var Rvebuilder = new GrapheneReinforcedRVEBuilderCHECKddmExample();
                model = Rvebuilder.GetModelAndBoundaryNodes().Item1;
                //RVEkanoninkhsGewmetriasBuilder.Reference2RVEExample1000ddm(model);
                int[][] subdElementIds = DdmCalculations.CalculateSubdElementIds(6, 6, 6, 3, 3, model);
                DdmCalculations.SeparateSubdomains(model, subdElementIds);
                foreach (int subdomainID in model.SubdomainsDictionary.Keys)
                {
                    Subdomain subd = model.SubdomainsDictionary[subdomainID];
                    DdmCalculations.PrintDictionary(subd.GlobalNodalDOFsDictionary, subd.TotalDOFs, subd.ID);
                }
            }

            


            bool use_domain_decomposer = false;
            if (use_domain_decomposer)
            {
                //i)
                DdmCalculationsPartb.MakeModelDictionariesZeroBasedForDecomposer(model);

                model.ConnectDataStructures();
                // ii)
                AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 8); //2o orisma arithmoos subdomains
                domainDecomposer.UpdateModel();
            }
            else
            {
                model.ConnectDataStructures();
                foreach (int subdomainID in model.SubdomainsDictionary.Keys)
                {
                    Subdomain subd = model.SubdomainsDictionary[subdomainID];
                    DdmCalculations.PrintDictionary(subd.GlobalNodalDOFsDictionary, subd.TotalDOFs, subd.ID);
                }
            }

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);

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
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 1;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            //h epomenhgrammh einai gia paradeigma ws pros to access
            //IAnalyzer childAnalyzer2 = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);


            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            //---------------------------------------------------------------------------------------------------------------------------------


            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });

            //comment section 2 palaia version
            //int increments = 1;
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            //analyzer.SetMaxIterations = 100;
            //analyzer.SetIterationsForMatrixRebuild = 1;



            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            #region

            //model.ConnectDataStructures();

            //SolverSkyline solver = new SolverSkyline(model);
            //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            ////LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            ////gia 2CZM
            ////Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            ////gia 3CZM
            //int increments = 1;
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            //analyzer.SetMaxIterations = 100;
            //analyzer.SetIterationsForMatrixRebuild = 1;

            //if (model__builder_choice==1)
            //{
            //    analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[12][DOFType.X],
            //model.NodalDOFsDictionary[12][DOFType.Y],
            //model.NodalDOFsDictionary[12][DOFType.Z]});
            //}

            //if (model__builder_choice == 2)
            //{
            //    analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[12][DOFType.X],
            //model.NodalDOFsDictionary[12][DOFType.Y],
            //model.NodalDOFsDictionary[12][DOFType.Z]});
            //}

            //if (model__builder_choice == 3)
            //{
            //    analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[12][DOFType.X],
            //model.NodalDOFsDictionary[12][DOFType.Y],
            //model.NodalDOFsDictionary[12][DOFType.Z]});
            //}

            //parentAnalyzer.BuildMatrices();
            //parentAnalyzer.Initialize();
            //parentAnalyzer.Solve();

            ////Console.WriteLine("checkPoint1 reached");
            //Console.WriteLine("Writing results for node 5");
            //Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            //Console.WriteLine(analyzer.Logs[1][0]);

            #endregion
        }

        static void Main(string[] args)
        {
            SolveRVEExample(); //|
        }

    }
}
