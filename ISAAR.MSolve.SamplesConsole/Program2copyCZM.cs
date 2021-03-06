﻿using ISAAR.MSolve.Analyzers;
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
    class Program2copyCZM
    {
        private static void SolveHexaBuilding()
        { 
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int HexaBuilder__CZM_choice;
            HexaBuilder__CZM_choice = 7;
            if (HexaBuilder__CZM_choice == 2)
            { HexaBuilder2CZM.MakeHexaBuilding(model); }
            if (HexaBuilder__CZM_choice == 3)
            { // shell kai cohesive shell
                HexaBuilder3CZM.MakeHexaBuilding(model); }
            if (HexaBuilder__CZM_choice == 4)
            { // shell kai cohesive shell anestrameno
                HexaBuilder3CZM.MakeHexaBuilding2(model); }
            if (HexaBuilder__CZM_choice == 5)
            { // shell mono tou
                HexaBuilder3CZM.MakeHexaBuilding3(model);}
            if (HexaBuilder__CZM_choice == 6)
            { // shell mono tou anestrameno
                HexaBuilder3CZM.MakeHexaBuilding4(model);}
            if (HexaBuilder__CZM_choice == 7)
            { // shell kai cohesive shell anestrameno isio
                HexaBuilder3CZM.MakeHexaBuilding5(model);}
            if (HexaBuilder__CZM_choice == 8)
            { // hexa kai cohesive building
                HexaBuilder4CZM.MakeHexaBuilding(model);
            }
            model.ConnectDataStructures();

            SolverSkyline solver = new SolverSkyline(model);
            ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //gia 2CZM
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //gia 3CZM
            Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 104, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            analyzer.SetMaxIterations = 100;
            analyzer.SetIterationsForMatrixRebuild = 1;
            //analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 420 });

            // apo theofilo
            if (HexaBuilder__CZM_choice == 2)
            { analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[4][DOFType.X],
            model.NodalDOFsDictionary[4][DOFType.Y],
            model.NodalDOFsDictionary[4][DOFType.Z]}); }
            if (HexaBuilder__CZM_choice == 3)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[1][DOFType.X],
            model.NodalDOFsDictionary[1][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
            }
            if (HexaBuilder__CZM_choice == 4)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[1][DOFType.X],
            model.NodalDOFsDictionary[1][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
            }
            if (HexaBuilder__CZM_choice == 5)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[1][DOFType.X],
            model.NodalDOFsDictionary[1][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
            }
            if (HexaBuilder__CZM_choice == 6)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[1][DOFType.X],
            model.NodalDOFsDictionary[1][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
            }
            if (HexaBuilder__CZM_choice == 7)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[1][DOFType.X],
            model.NodalDOFsDictionary[1][DOFType.Y],
            model.NodalDOFsDictionary[1][DOFType.Z]});
            }
            if (HexaBuilder__CZM_choice == 8)
            {
                analyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] {
            model.NodalDOFsDictionary[8][DOFType.X],
            model.NodalDOFsDictionary[8][DOFType.Y],
            model.NodalDOFsDictionary[8][DOFType.Z]});
            }


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
            Console.WriteLine("Writing results for node 4");
            Console.WriteLine("Dof and Values for Displacement X, Y, Z");
            Console.WriteLine(analyzer.Logs[1][0]);

        }

        //static void Main(string[] args)
        //{
        //    SolveHexaBuilding();
        //}
    }
}
