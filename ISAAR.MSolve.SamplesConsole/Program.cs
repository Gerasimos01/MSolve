using System;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MGroup.Stochastic;
using MGroup.Stochastic.Structural;
using MGroup.Stochastic.Structural.Example;

namespace ISAAR.MSolve.SamplesConsole
{
    class Program
    {
        private const int subdomainID = 0;

        static void Main(string[] args)
        {
            //SolveBuildingInNoSoilSmall();
            //TrussExample.Run();
            //FEM.Cantilever2D.Run();
            //FEM.Cantilever2DPreprocessor.Run();
            //FEM.WallWithOpenings.Run();
            //SeparateCodeCheckingClass.Check06();
            //SolveBuildingInNoSoilSmall();
            //SolveBuildingInNoSoilSmallDynamic();
            //SolveStochasticMaterialBeam2DWithBruteForceMonteCarlo();
            //CNTExamples.CNT_4_4_DisplacementControl();
            //CNTExamples.CNT_4_4_NewtonRaphson();
            //Tests.FEM.Shell8andCohesiveNonLinear.RunTest();
            //AppliedDisplacementExample.Run();

            //Logging.PrintForceDisplacementCurve.CantileverBeam2DCorotationalLoadControl();

            //SuiteSparseBenchmarks.MemoryConsumptionDebugging();
            //SolverBenchmarks.SuiteSparseMemoryConsumptionDebugging();
            //NRNLAnalyzerDevelopTest.SolveDisplLoadsExample();
            //SeparateCodeCheckingClass4.Check05bStressIntegrationObje_Integration();
            //SeparateCodeCheckingClass4.Check_Graphene_rve_Obje_Integration();
            //IntegrationElasticCantileverBenchmark.RunExample();
            //OneRveExample.Check_Graphene_rve_serial();
            //BondSlipTest.CheckStressStrainBonSlipMaterial();
            //OneRveExample.Check_Graphene_rve_parallel();
            //LinearRves.CheckShellScaleTransitionsAndMicrostructure();
            //SolveCantileverWithStochasticMaterial();
            //SeparateCodeCheckingClass5.StiffnessMatrixOutputWrite();
            //SeparateCodeCheckingClass6.StiffnessMatrixOutputWrite();
            //SeparateCodeCheckingClass5b.CheckSolution();
            //SeparateCodeCheckingClass6b.RunExample();
            //SeparateCodeCheckingClass5b.RunExample();
            //ISAAR.MSolve.Tests.SeparateCodeCheckingClass5b.RunExampleSerial();
            //SeparateCodeCheckingClass5b_b.StiffnessMatrixOutputWrite();
            // //(Model model1, double[] uc1) = SeparateCodeCheckingClass5b_b.RunExample();
            //(Model model2, double[] uc2) = SeparateCodeCheckingClass5b_b.RunExampleSerial();
            //SeparateCodeCheckingClass5b_c.StiffnessMatrixOutputWrite();
            SeparateCodeCheckingClass5b_c.RunExample();
            ///SeparateCodeCheckingClass5b_c.RunExampleSerial();
            (Model model, double[] uc) = SeparateCodeCheckingClass5b_c.RunExampleSerial();
            //SeparateCodeCheckingClass5b_c1.StiffnessMatrixOutputWrite();
            //SeparateCodeCheckingClass5b_c1.RunExample();

            //SeparateCodeCheckingClass7_b_b.Check(); //nonlinear strains example etc.

            //SeparateCodeCheckingClass5b_c_constraints.RunExample();
        }

        private static void SolveBuildingInNoSoilSmall()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);
            model.Loads.Add(new Load() { Amount = -100, Node = model.Nodes[21], DOF = StructuralDof.TranslationX });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");
        }

        private static void SolveBuildingInNoSoilSmallDynamic()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));
            BeamBuildingBuilder.MakeBeamBuilding(model, 20, 20, 20, 5, 4, model.NodesDictionary.Count + 1,
                model.ElementsDictionary.Count + 1, subdomainID, 4, false, false);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.01, 0.1);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            int monitorDof = 420;
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monitorDof}, u = {log.DOFValues[monitorDof]}");

            //TODO: No loads have been defined so the result is bound to be 0.
        }

        private static void SolveCantileverWithStochasticMaterial()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new[] { 0d, 0d, 0d });
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }
    }
}
