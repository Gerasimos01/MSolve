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
    class Program2copyElegxoiBenc
    {

        private static void SolveBencExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice = 3;

            if (model__builder_choice == 1) // 
            { ParadeigmataElegxwnBuilder.Example_cohesive_hexa_orthi_constr_anw_benc1(model); }
            if (model__builder_choice == 2) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveShellPaktwsh(model); }
            if (model__builder_choice == 3) // diatmhsh elastic
            { RVEExamplesBuilder.FewElementsRVECheckExample(model); }
            if (model__builder_choice == 4) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder(model,850); }
            if (model__builder_choice == 5) // 
            { RVEExamplesBuilder.OriginalRVECholExample(model); }

            if (model__builder_choice == 6) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveShellPaktwshRAM(model); }
            if (model__builder_choice == 7) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilderRAM(model, 850); }
            if (model__builder_choice == 8) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilderRAM2(model, 850); }
            if (model__builder_choice == 9) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveShellPaktwshRAM2(model); }
            if (model__builder_choice == 10) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveRAMShellPaktwsh(model); }
            if (model__builder_choice == 11) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveShellRAM_1Paktwsh(model); }
            if (model__builder_choice == 12) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilderRAM_1(model, 850); }
            if (model__builder_choice == 13) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveRAM_1ShellPaktwsh(model); }
            if (model__builder_choice == 14) // 
            { ParadeigmataElegxwnBuilder.ShellAndCohesiveRAM_11ShellPaktwsh(model); }
            if (model__builder_choice == 15) // diatmhsh elastic
            { RVEExamplesBuilder.FewElementsRVECheckExample_embedder_2(model); }


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

            if (model__builder_choice == 1)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[6][DOFType.X],
            //model.NodalDOFsDictionary[6][DOFType.Y],
            model.NodalDOFsDictionary[6][DOFType.Z]});
            }

            if (model__builder_choice == 2)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            //model.NodalDOFsDictionary[12][DOFType.X],
            //model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
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
            //model.NodalDOFsDictionary[12][DOFType.X],
            //model.NodalDOFsDictionary[12][DOFType.Y],
            model.NodalDOFsDictionary[20][DOFType.X]});
            }

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //Console.WriteLine("checkPoint1 reached");
            Console.WriteLine("Writing results for node 6");
            Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            Console.WriteLine(analyzer.Logs[1][0]);
        }

        //static void Main(string[] args)
        //{
        //    SolveBencExample(); //|
        //}



    }
}
