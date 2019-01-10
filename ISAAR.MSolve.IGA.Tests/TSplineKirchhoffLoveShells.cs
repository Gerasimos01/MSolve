using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public class TSplineKirchhoffLoveShells
	{
		[Fact]
		public void CantileverShellBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "CantileverShell";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);
			modelReader.CreateTSplineShellsModelFromFile();

			model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStress)
			{
				PoissonRatio = 0.0,
				YoungModulus = 100
			};
			model.PatchesDictionary[0].Thickness = 1;

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp=>cp.X<3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){ DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X >49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
			}

			var solverBuilder = new DenseMatrixSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, filename);
			paraview.CreateParaviewFile();

			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-8));
			}
			
		}


		[Fact]
		public void CantileverShellMaterialBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\CantileverShell.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);
			
			var thickness = 1.0;

			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D
			{
				PoissonRatio = 0.0,
				YoungModulus = 100,
			}, thickness);
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X > 49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
			}

			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-6));
			}

		}

		[Fact]
		public void CantileverShellMaterialBenchmark_v2()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\CantileverShell.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);

			var thickness = 1.0;

			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D
			{
				PoissonRatio = 0.0,
				YoungModulus = 100,
			}, thickness);
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 3))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X > 49.8))
			{
				model.Loads.Add(new Load()
				{
					Amount = -0.5,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
			}

			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();


			var expectedSolutionVector = new Vector(new double[]
			{
				0, 0, -306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924153, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0,
				-3454.810407, 0, 0, -5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0,
				-15000.00025, 0, 0, -306.1224493, 0, 0, -1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0,
				-8702.623837, 0, 0, -11785.71389, 0, 0, -13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224494, 0, 0,
				-1552.478088, 0, 0, -3454.810449, 0, 0, -5881.924038, 0, 0, -8702.623837, 0, 0, -11785.71389, 0, 0,
				-13928.57157, 0, 0, -14999.99909, 0, 0, -306.1224369, 0, 0, -1552.47811, 0, 0, -3454.810407, 0, 0,
				-5881.924117, 0, 0, -8702.623683, 0, 0, -11785.71423, 0, 0, -13928.57093, 0, 0, -15000.00025, 0, 0,
				-306.122431, 0, 0, -1552.478121, 0, 0, -3454.810388, 0, 0, -5881.924154, 0, 0, -8702.62361, 0, 0,
				-11785.71439, 0, 0, -13928.57064, 0, 0, -15000.0008
			});
			for (int i = 0; i < expectedSolutionVector.Length; i++)
			{
				Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], solver.LinearSystems[0].Solution[i],
					1e-6));
			}

		}

		[Fact]
		public void SquareShellMaterialMultiscaleBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			string filename = "..\\..\\..\\InputFiles\\square_unstructured.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filename);

			var thickness = 1.0;

			//VectorExtensions.AssignTotalAffinityCount();
			IdegenerateRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHost();
			var material1 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimu(homogeneousRveBuilder1, true);


			modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, material1, thickness);
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X < 1e-6))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint {DOF=DOFType.X});
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Z });
			}
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.Y < 1e-6))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Z });
			}
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.X > 1 - 1e-6))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Z });
			}
			foreach (var controlPoint in model.ControlPointsDictionary.Values.Where(cp => cp.Y > 1 - 1e-6))
			{
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.X });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Y });
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint { DOF = DOFType.Z });
			}

			foreach (var controlPoint in model.ControlPointsDictionary.Values)
			{
				model.Loads.Add(new Load()
				{
					Amount = -10,
					ControlPoint = model.ControlPointsDictionary[controlPoint.ID],
					DOF = DOFType.Z
				});
			}
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			var provider = new ProblemStructural_v2(model, solver);

			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var solutionData=solver.LinearSystems[0].Solution.CopyToArray();
			PrintUtilities.WriteToFileVector(solutionData, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\U_sunol_1.txt");
		}


		//[Fact]
		public void SimpleHoodBenchmark()
		{
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "attempt2";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);

            //var thickness = 1.0;
            //modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.LinearMaterial,new ShellElasticMaterial2D
            //{
            //	PoissonRatio = 0.3,
            //	YoungModulus = 1e5,
            //}, thickness);

            modelReader.CreateTSplineShellsModelFromFile();
            model.PatchesDictionary[0].Material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                PoissonRatio = 0.3,
                YoungModulus = 10000 //PROSOXH: na antikatastathei me to kanoniko material
            };
            model.PatchesDictionary[0].Thickness = 1;

            for (int i = 0; i < 100; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			for (int i = model.ControlPoints.Count-100; i < model.ControlPoints.Count; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.Loads.Add(new Load()
				{
					Amount = 100,
					ControlPoint = model.ControlPointsDictionary[id],
					DOF = DOFType.Z
				});
			}

			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), new NullReordering());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var paraview= new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution,filename);
			paraview.CreateParaviewFile();
		}

		[Fact] //commented out: requires mkl and suitesparse can't be test
		public void SimpleHoodBenchmarkMKL()
		{
            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
			VectorExtensions.AssignTotalAffinityCount();
			Model model = new Model();
			var filename = "attempt2";
			string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
			IGAFileReader modelReader = new IGAFileReader(model, filepath);

            var runMs = true;
            var transformationA = false;

            if (runMs)
            {
                IdegenerateRVEbuilder RveBuilder3 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostData(1);
                var BasicMaterial = new Shell2dRVEMaterialHost(2, 2, 0, RveBuilder3);
                var thickness = 1.0;
                modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, BasicMaterial, thickness);
            }
            else
            {
                if (transformationA)
                {
                    var thickness = 1.0;
                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D()
                    {
                        PoissonRatio = 0.4,
                        YoungModulus = 3.5
                    }, thickness);
                }
                else
                {
                    var thickness = 1.0;
                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2Dtransformationb()
                    {
                        PoissonRatio = 0.4,
                        YoungModulus = 3.5
                    }, thickness);
                }
            }

			for (int i = 0; i < 100; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint(){DOF = DOFType.X});
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
				model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
			}

			for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
			{
				var id = model.ControlPoints[i].ID;
				model.Loads.Add(new Load()
				{
					Amount = 100,
					ControlPoint = model.ControlPointsDictionary[id],
					DOF = DOFType.Z
				});
			}
			var solverBuilder = new SuiteSparseSolver.Builder();
			solverBuilder.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
			ISolver_v2 solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural_v2(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer_v2(solver);
			var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

            var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, filename);
            paraview.CreateParaviewFile();

            double[] solutiondata =solver.LinearSystems[0].Solution.CopyToArray();
            PrintUtilities.WriteToFileVector(new double[1] { new Vector(solutiondata).Norm }, $"..\\..\\..\\OutputFiles\\{filename}SolutionNorm.txt");
        }

        [Fact] //commented out: requires mkl and suitesparse can't be test
        public void SimpleHoodBenchmarkMKLStochastic()
        {
            var runMs = true;
            IdegenerateRVEbuilder RveBuilder3 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostData(1);
            var BasicMaterial = new Shell2dRVEMaterialHost(50, 2, 0, RveBuilder3);
            int totalsimulations = 3;

            #region Genika settings
            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
            VectorExtensions.AssignTotalAffinityCount();
            #endregion

            for (int simulation_id = 1; simulation_id < totalsimulations + 1; simulation_id++)
            {
                Model model = new Model();
                var filename = "attempt2";
                string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
                IGAFileReader modelReader = new IGAFileReader(model, filepath);

                var transformationA = false;

                if (runMs)
                {
                    var thickness = 0.015;
                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, BasicMaterial, thickness);
                }
                else
                {
                    if (transformationA)
                    {
                        var thickness = 0.015;
                        modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D()
                        {
                            PoissonRatio = 0.4,
                            YoungModulus = 3.5
                        }, thickness);
                    }
                    else
                    {
                        var thickness = 0.015;
                        modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2Dtransformationb()
                        {
                            PoissonRatio = 0.4,
                            YoungModulus = 3.5
                        }, thickness);
                    }
                }

                for (int i = 0; i < 100; i++)
                {
                    var id = model.ControlPoints[i].ID;
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.X });
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
                }

                for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
                {
                    var id = model.ControlPoints[i].ID;
                    model.Loads.Add(new Load()
                    {
                        Amount = 1,
                        ControlPoint = model.ControlPointsDictionary[id],
                        DOF = DOFType.Z
                    });
                }
                var solverBuilder = new SuiteSparseSolver.Builder();
                solverBuilder.DofOrderer = new DofOrderer(
                    new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
                SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

                //var solverBuilder = new SkylineSolver.Builder();
                //SkylineSolver solver = solverBuilder2.BuildSolver(model);

                // Structural problem provider
                var provider = new ProblemStructural_v2(model, solver);

                // Linear static analysis
                var childAnalyzer = new LinearAnalyzer_v2(solver);
                var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

                // Run the analysis
                parentAnalyzer.Initialize();
                parentAnalyzer.Solve();

                string outputString = "Analysis_no_" + simulation_id.ToString() + "_output";

                var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, outputString);
                paraview.CreateParaviewFile();

                double[] solutiondata = solver.LinearSystems[0].Solution.CopyToArray();
                PrintUtilities.WriteToFileVector(new double[1] { new Vector(solutiondata).Norm }, $"..\\..\\..\\OutputFiles\\{outputString}SolutionNorm.txt");

                #region empty data
                model.Clear();
                model = null;
                modelReader = null;
                solverBuilder = null;
                solver.Dispose();
                provider = null; childAnalyzer = null; parentAnalyzer = null; paraview = null; solutiondata = null;
                #endregion
            }
        }

        [Fact] //commented out: requires mkl and suitesparse can't be test
        public void SimpleHoodBenchmarkMKLStochasticElastic()
        {                        
            int totalsimulations = 3;

            #region Genika settings
            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
            VectorExtensions.AssignTotalAffinityCount();
            #endregion

            for (int simulation_id = 1; simulation_id < totalsimulations + 1; simulation_id++)
            {
                Model model = new Model();
                var filename = "attempt2";
                string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
                IGAFileReader modelReader = new IGAFileReader(model, filepath);

                var runMs = false;
                var transformationA = false;

                if (runMs)
                {
                   

                }
                else
                {
                    if (transformationA)
                    {
                        var thickness = 0.015;
                        modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D()
                        {
                            PoissonRatio = 0.4,
                            YoungModulus = 3.5
                        }, thickness);
                    }
                    else
                    {
                        var thickness = 0.015;
                        modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2Dtransformationb()
                        {
                            PoissonRatio = 0.4,
                            YoungModulus = 3.5
                        }, thickness);
                    }
                }

                for (int i = 0; i < 100; i++)
                {
                    var id = model.ControlPoints[i].ID;
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.X });
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
                    model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
                }

                for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
                {
                    var id = model.ControlPoints[i].ID;
                    model.Loads.Add(new Load()
                    {
                        Amount = 1,
                        ControlPoint = model.ControlPointsDictionary[id],
                        DOF = DOFType.Z
                    });
                }
                var solverBuilder = new SuiteSparseSolver.Builder();
                solverBuilder.DofOrderer = new DofOrderer(
                    new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
                SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

                //var solverBuilder = new SkylineSolver.Builder();
                //SkylineSolver solver = solverBuilder2.BuildSolver(model);

                // Structural problem provider
                var provider = new ProblemStructural_v2(model, solver);

                // Linear static analysis
                var childAnalyzer = new LinearAnalyzer_v2(solver);
                var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

                // Run the analysis
                parentAnalyzer.Initialize();
                parentAnalyzer.Solve();

                string outputString = "Analysis_no_" + simulation_id.ToString() + "_output";

                var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, outputString);
                paraview.CreateParaviewFile();

                double[] solutiondata = solver.LinearSystems[0].Solution.CopyToArray();
                PrintUtilities.WriteToFileVector(new double[1] { new Vector(solutiondata).Norm }, $"..\\..\\..\\OutputFiles\\{outputString}SolutionNorm.txt");

                #region empty data
                model.Clear();
                model = null;
                modelReader = null;
                solverBuilder = null;
                solver.Dispose();
                provider = null; childAnalyzer = null; parentAnalyzer = null; paraview = null; solutiondata = null;
                #endregion
            }
        }

        [Fact]
        public void SimpleHoodBenchmarkMKLStochasticElastic2()
        {

            int totalsimulations = 3;

            #region Genika settings
            //LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
            VectorExtensions.AssignTotalAffinityCount();
            #endregion

            for (int simulation_id = 1; simulation_id < totalsimulations + 1; simulation_id++)
            {
                Analyze(simulation_id);
            }
        }

        private void Analyze(int analysis_No)
        {
            Model model = new Model();
            var filename = "attempt2";
            string filepath = $"..\\..\\..\\InputFiles\\{filename}.iga";
            IGAFileReader modelReader = new IGAFileReader(model, filepath);

            var runMs = false;
            var transformationA = false;

            if (runMs)
            {


            }
            else
            {
                if (transformationA)
                {
                    var thickness = 0.015;
                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2D()
                    {
                        PoissonRatio = 0.4,
                        YoungModulus = 3.5
                    }, thickness);
                }
                else
                {
                    var thickness = 0.015;
                    modelReader.CreateTSplineShellsModelFromFile(IGAFileReader.TSplineShellTypes.ThicknessMaterial, new ShellElasticMaterial2Dtransformationb()
                    {
                        PoissonRatio = 0.4,
                        YoungModulus = 3.5
                    }, thickness);
                }
            }

            for (int i = 0; i < 100; i++)
            {
                var id = model.ControlPoints[i].ID;
                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.X });
                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Y });
                model.ControlPointsDictionary[id].Constrains.Add(new Constraint() { DOF = DOFType.Z });
            }

            for (int i = model.ControlPoints.Count - 100; i < model.ControlPoints.Count; i++)
            {
                var id = model.ControlPoints[i].ID;
                model.Loads.Add(new Load()
                {
                    Amount = 1,
                    ControlPoint = model.ControlPointsDictionary[id],
                    DOF = DOFType.Z
                });
            }
            var solverBuilder = new SuiteSparseSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(
                new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
            SuiteSparseSolver solver = solverBuilder.BuildSolver(model);

            //var solverBuilder = new SkylineSolver.Builder();
            //SkylineSolver solver = solverBuilder2.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural_v2(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer_v2(solver);
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            string outputString = "Analysis_no_" + analysis_No.ToString() + "_output";

            var paraview = new ParaviewTsplineShells(model, solver.LinearSystems[0].Solution, outputString);
            paraview.CreateParaviewFile();

            double[] solutiondata = solver.LinearSystems[0].Solution.CopyToArray();
            PrintUtilities.WriteToFileVector(new double[1] { new Vector(solutiondata).Norm }, $"..\\..\\..\\OutputFiles\\{outputString}SolutionNorm.txt");

            #region empty data
            model.Clear();
            model = null;
            modelReader = null;
            solverBuilder = null;
            solver.Dispose();
            provider = null; childAnalyzer = null; parentAnalyzer = null; paraview = null; solutiondata = null;
            #endregion
        }
    }
}
