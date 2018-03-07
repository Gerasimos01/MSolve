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
    class Program2copyElegxoiEmbedded
    {
        private static void SolveExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice = 4;

            if (model__builder_choice == 1) // Hexa8 kanoniko me NL analyzer
            { EmbeddedExamplesBuilder.HexaElementsOnly(model); }
            if (model__builder_choice == 2) // Beam3d me NL analyzer
            { EmbeddedExamplesBuilder.BeamElementOnly(model); }
            if (model__builder_choice == 3) // Beam3d me NL analyzer
            { EmbeddedExamplesBuilder.ExampleWithEmbedded(model); }

            if (model__builder_choice == 4) // Hexa8 kanoniko me NL analyzer paradeigma me Vasili Von mises
            { EmbeddedExamplesBuilder.HexaElementsOnlyVonMises(model); }


            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);
            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //gia 2CZM
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //gia 3CZM
            int increments = 2;
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
            model.NodalDOFsDictionary[14][DOFType.X],
            model.NodalDOFsDictionary[14][DOFType.Y],
            model.NodalDOFsDictionary[14][DOFType.Z]});
            }

            if (model__builder_choice == 3)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[12][DOFType.X],
            model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[12][DOFType.Z]});
            }


            if (model__builder_choice == 4)
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
        //    SolveExample();
        //}

    }
}
