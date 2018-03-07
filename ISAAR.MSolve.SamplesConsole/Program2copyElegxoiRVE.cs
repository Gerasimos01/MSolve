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
    class Program2copyElegxoiRVE
    {
        private static void SolveRVEExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice = 7;

            if (model__builder_choice == 1) // 
            { RVEExamplesBuilder.OneElementRVECheckExampleConstrained(model); }
            if (model__builder_choice == 2) // 
            { RVEExamplesBuilder.FewElementsRVECheckExample(model); }
            if (model__builder_choice == 3) // Beam3d me NL analyzer //<- mallon lathos sxolio exei apomeinei
            { RVEExamplesBuilder.OriginalRVECholExample(model); }
            if (model__builder_choice == 4) // 
            { RVEExamplesBuilder.FewElementsRVECheckExample2GrapheneSheets(model); }
            if (model__builder_choice == 5) // 
            { RVEExamplesBuilder.Reference1RVEExample10000(model); }
            if (model__builder_choice == 6) // 
            { RVEExamplesBuilder.Reference2RVEExample50000(model); }
            if (model__builder_choice == 7) // 
            { RVEExamplesBuilder.Reference2RVEExample50000withRenumbering(model); }
            if (model__builder_choice == 8) // 
            { RVEExamplesBuilder.Reference2RVEExample10000withRenumbering(model); }


            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);
            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //gia 2CZM
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //gia 3CZM
            int increments = 1;
            Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            analyzer.SetMaxIterations = 100;
            analyzer.SetIterationsForMatrixRebuild = 1;

            if (model__builder_choice==1)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[12][DOFType.X],
            model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[12][DOFType.Z]});
            }

            if (model__builder_choice == 2)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[12][DOFType.X],
            model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[12][DOFType.Z]});
            }

            if (model__builder_choice == 3)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[12][DOFType.X],
            model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[12][DOFType.Z]});
            }

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Console.WriteLine("checkPoint1 reached");
            Console.WriteLine("Writing results for node 5");
            Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            Console.WriteLine(analyzer.Logs[1][0]);

        }

        //static void Main(string[] args)
        //{
        //    SolveRVEExample(); //|
        //}

    }
}
