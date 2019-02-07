//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using ISAAR.MSolve.Analyzers;
//using ISAAR.MSolve.Discretization;
//using ISAAR.MSolve.Discretization.Interfaces;
//using ISAAR.MSolve.IGA.Elements;
//using ISAAR.MSolve.IGA.Entities;
//using ISAAR.MSolve.IGA.Postprocessing;
//using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
//using ISAAR.MSolve.IGA.Readers;
//using ISAAR.MSolve.LinearAlgebra;
//using ISAAR.MSolve.Materials;
//using ISAAR.MSolve.MultiscaleAnalysis;
//using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
//using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
//using ISAAR.MSolve.MultiscaleAnalysisMerge;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
//using ISAAR.MSolve.Problems;
//using ISAAR.MSolve.Solvers;
//using ISAAR.MSolve.Solvers.Direct;
//using ISAAR.MSolve.Solvers.Interfaces;
//using ISAAR.MSolve.Solvers.Ordering;
//using ISAAR.MSolve.Solvers.Ordering.Reordering;
//using ISAAR.MSolve.Solvers.Skyline;
//using Xunit;
//namespace ISAAR.MSolve.IGA.Tests
//{
//    public class TSplineKirchhoffLoveShellsElasticStresses
//    {
//        [Fact] //commented out: requires mkl and suitesparse can't be test
//        public void SimpleHoodBenchmarkMKLPrintStresses()
//        {
//            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
//            VectorExtensions.AssignTotalAffinityCount();
//            Model model = new Model();
//            var filename = "attempt2";
//            string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
//            IGAFileReader modelReader = new IGAFileReader(model, filepath);

//            var runMs = false;
//            var transformationA = false;

//            if (runMs)
//            {
//                IdegenerateRVEbuilder RveBuilder3 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostData(1);
//                var BasicMaterial = new Shell2dRVEMaterialHost(2, 2, 0, RveBuilder3);
//                var thickness = 0.015;
//                modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, BasicMaterial, thickness);
//            }
//            else
//            {
//                if (transformationA)
//                {
//                    var thickness = 0.015;
//                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D()
//                    {
//                        PoissonRatio = 0.4,
//                        YoungModulus = 3.5
//                    }, thickness);
//                }
//                else
//                {
//                    var thickness = 0.015;
//                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2Dtransformationb()
//                    {
//                        PoissonRatio = 0.4,
//                        YoungModulus = 3.5
//                    }, thickness);
//                }
//            }

//            for (int i = 0; i < 100; i++)
//            {
//                var id = model.ControlPoints[i].ID;
//                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.X });
//                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
//                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
//            }

//            for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
//            {
//                var id = model.ControlPoints[i].ID;
//                model.Loads.Add(new Load()
//                {
//                    Amount = 100* 3.2e-7, //wste max u =0.0896m gia graphene
//                    ControlPoint = model.ControlPointsDictionary[id],
//                    DOF = DOFType.Z
//                });
//            }


//            var solverBuilder = new DenseMatrixSolver.Builder();
//            solverBuilder.DofOrderer = new DofOrderer(
//                new NodeMajorDofOrderingStrategy(), new NullReordering());
//            ISolver_v2 solver = solverBuilder.BuildSolver(model);

//            var provider = new ProblemStructural_v2(model, solver);

//            var childAnalyzer = new LoadControlAnalyzer_v2.Builder(model, solver, provider, 2).Build();

//            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

//            parentAnalyzer.Initialize();
//            parentAnalyzer.Solve();

//            foreach (Element element in model.Patches[0].Elements)
//            {
//                var e1 = element.IElementType as TSplineKirchhoffLoveShellElementMaterial;
//                Dictionary<GaussLegendrePoint3D, double[]> stressesTodisplay1 = e1.GetStressesAtLowerSurfaceForLogging();
//                Dictionary<GaussLegendrePoint3D, double[]> stressesTodisplay2 = e1.GetStressesAtLowerSurfaceForLogging();
//            }
//            var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, filename);
//            paraview.CreateParaviewFile();

//            double[] solutiondata = solver.LinearSystems[0].Solution.CopyToArray();
//            PrintUtilities.WriteToFileVector(new double[1] { new Vector(solutiondata).Norm }, $"..\\..\\..\\OutputFiles\\{filename}SolutionNorm.txt");
//        }
//    }
//}
