using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program2NL
    {
        private static void SolveHexaBuilding()
        { 
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
    
            ArchBuilder.MakeArchBuilding(model);

            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);
            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            Analyzers.NewtonRaphsonNonLinearAnalyzer2 analyzer = new NewtonRaphsonNonLinearAnalyzer2(solver, solver.SubdomainsDictionary, provider, 10, model.TotalDOFs);//1. increments einai to 10
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            analyzer.SetMaxIterations = 100;
            analyzer.SetIterationsForMatrixRebuild = 1;
            //analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            // apo theofilo
            analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[4][DOFType.X],
            model.NodalDOFsDictionary[4][DOFType.Y],
            model.NodalDOFsDictionary[4][DOFType.Z] });
            //model.ElementsDictionary[1][]
            //model.NodalDOFsDictionary[17][DOFType.Z]});
            //ews edw

            //// apo theofilo
            //analyzer.LogFactories[2] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[5][DOFType.X],
            //model.NodalDOFsDictionary[5][DOFType.Y],
            //model.NodalDOFsDictionary[5][DOFType.Z] });
            ////ews edw

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Console.WriteLine("checkPoint1 reached");
            Console.WriteLine("Writing results for node 5");
            Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            Console.WriteLine(analyzer.Logs[1][0]);

        }

        static void Main(string[] args)
        {
            SolveHexaBuilding();
        }
    }
}
