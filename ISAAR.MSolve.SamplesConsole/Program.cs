using System;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
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

            //Random geometry and integration of constitutive
            //SeperateIntegrationClassCheck.RunExample();

            //Random geometry and solution of one rhs
            //(Model model3, double[] uc3) = SeparateCodeCheckingClass5b_bNEWhstam.RunExample();
            //(Model model4, double[] uc4) = SeparateCodeCheckingClass5b_bNEW.RunExampleSerial();

            //Check model data output
            ///SeparateCodeCheckingClass8Output.CheckOutputWriteFile();

            //dokimi olwn px develop extraConstraints(devvelopbDuplicateDevelop) klp.
            //(Model model3, double[] uc3) = SeparateCodeCheckingClass5b_bNEW.RunExample();
            //(Model model4, double[] uc4) = SeparateCodeCheckingClass5b_bNEW.RunExampleSerial();
            (Model model3, double[] uc3) = SeparateCodeCheckingClass5b_bNEW_debug.RunExample();
            (Model model4, double[] uc4) = SeparateCodeCheckingClass5b_bNEW_debug.RunExampleSerial();

            (Model model1, double[] uc1) = SeparateCodeCheckingClass5b_b.RunExample();
            (Model model2, double[] uc2) = SeparateCodeCheckingClass5b_b.RunExampleSerial();
            //(Model model2, double[] uc2) = SeparateCodeCheckingClass5b_b_cantilever.RunExampleSerial();
            ////PrintReorderingModel1ToModel2(model1, model2);
            //(int num1, int num2) = CountElements(model1, model2);

            //SeparateCodeCheckingClass5b_c.StiffnessMatrixOutputWrite();
            SeparateCodeCheckingClass5b_c.RunExample();  // // 
            ///SeparateCodeCheckingClass5b_c.RunExampleSerial();
            (Model model, double[] uc) = SeparateCodeCheckingClass5b_c.RunExampleSerial(); // // 

            SeparateCodeCheckingClass5b_c1.StiffnessMatrixOutputWrite();
            SeparateCodeCheckingClass5b_c1.RunExample();

            //SeparateCodeCheckingClass7_b_b.Check(); //nonlinear strains example etc.

            SeparateCodeCheckingClass5b_c_constraints.RunExample();

            SeparateCodeCheckingClass_c_alte_develop.NLRVEStrainParralelSolution();
            SeparateCodeCheckingClass_c_alte_develop.NLRVEStrainParralelSolution_FETI1(); //TWRINO




        }

        private static (int num1, int num2) CountElements(Model model1, Model model2)
        {
            int num1 = 0;
            foreach(Subdomain subdomain in model1.Subdomains)
            {
                num1 += subdomain.Elements.Count;
            }

            int num2 = 0;
            foreach (Subdomain subdomain in model2.Subdomains)
            {
                num2 += subdomain.Elements.Count;
            }
            return (num1, num2);
        }

        private static void PrintReorderingModel1ToModel2(Model model1, Model model2)
        {
            var destinationNumbers = new int[model1.GlobalDofOrdering.NumGlobalFreeDofs];
            foreach (Node node in model1.GlobalDofOrdering.GlobalFreeDofs.GetRows())
            {
                foreach(IDofType doftype in model1.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node))
                {
                    destinationNumbers[model1.GlobalDofOrdering.GlobalFreeDofs[node, doftype]] = model2.GlobalDofOrdering.GlobalFreeDofs[model2.NodesDictionary[node.ID], doftype];
                }
            }

            var path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\destinationNumbering.txt";
            DdmCalculationsGeneral.WriteToFileVector(destinationNumbers, path);

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
