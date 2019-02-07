using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.Solvers.Skyline;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
	public static class TSplineKirchhoffLoveShellsIncremental
	{
        private const int pathcID = 0;

        [Fact]
        private static void RunTest()
        {
            int increments = 2;
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements(increments);
            TotalDisplacementsPerIterationLog computedDisplacements = SolveModel(increments);
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
        }

        private static TotalDisplacementsPerIterationLog SolveModel(int increments)
        {
            Model model = GetCantileverShellMaterialBenchmarkModel(); //exei assignAffinity kai model.connectDataStructures

            var solverBuilder = new DenseMatrixSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(
	            new NodeMajorDofOrderingStrategy(), new NullReordering());
            ISolver_v2 solver = solverBuilder.BuildSolver(model);

            var provider = new ProblemStructural_v2(model, solver);
			
			var childAnalyzer = new LoadControlAnalyzer_v2.Builder(model,solver,provider,2).Build();

			var watchDofs = new Dictionary<int, int[]>(); //TODO
			var subdomainDofsIDs = new int[model.GlobalDofOrdering.NumGlobalFreeDofs]; for (int i1 = 0; i1 < model.GlobalDofOrdering.NumGlobalFreeDofs; i1++) { subdomainDofsIDs[i1] = i1; }
			watchDofs.Add(pathcID, subdomainDofsIDs);
			var log1= new TotalDisplacementsPerIterationLog(watchDofs);
			//childAnalyzer.TotalDisplacementsPerIterationLog = log1;
            var parentAnalyzer = new StaticAnalyzer_v2(model, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return log1;

			//var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
			//linearSystems[pathcID] = new SkylineLinearSystem(pathcID, model.Patches[0].Forces);

			//ProblemStructural provider = new ProblemStructural(model, linearSystems);

			//var solver = new SolverSkyline(linearSystems[pathcID]);
			//var linearSystemsArray = new[] { linearSystems[pathcID] };
			//var subdomainUpdaters = new[] { new NonLinearPatchUpdater(model.Patches[0]) };
			//var subdomainMappers = new[] { new PatchGlobalMapping(model.Patches[0]) };

			//var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);


			//var watchDofs = new Dictionary<int, int[]>(); //TODO
			//var subdomainDofsIDs = new int[model.TotalDOFs]; for (int i1=0; i1 < model.TotalDOFs; i1++) { subdomainDofsIDs[i1] = i1; }
			//watchDofs.Add(pathcID, subdomainDofsIDs);
			//var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
			//childAnalyzer.IncrementalDisplacementsLog = log1;

			//childAnalyzer.SetMaxIterations = 100;
			//childAnalyzer.SetIterationsForMatrixRebuild = 1;

			//StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

			//parentAnalyzer.BuildMatrices();
			//parentAnalyzer.Initialize();
			//parentAnalyzer.Solve();



		}

        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements(int increments )
        {
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

            int totalDofs = expectedSolutionVector.Length;

            var expectedDisplacements = new Dictionary<int, double>[increments];

            for(int i1=0; i1<increments; i1++ )
            {
                var incrementalDisplacements = new Dictionary<int, double>(totalDofs);
                for (int i2 = 0; i2 < totalDofs; i2++)
                {
                    double dofValue = ((double)i1 + 1) / ((double)increments)*expectedSolutionVector[i2];
                    incrementalDisplacements.Add(i2, dofValue);
                }
                expectedDisplacements[i1] = incrementalDisplacements;
            }

            return expectedDisplacements;

        }

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, TotalDisplacementsPerIterationLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-6);

            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, pathcID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        
		public static Model GetCantileverShellMaterialBenchmarkModel()
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
				model.ControlPointsDictionary[controlPoint.ID].Constrains.Add(new Constraint() {DOF=DOFType.X});
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
			model.ConnectDataStructures();

            //var linearSystems = new Dictionary<int, ILinearSystem>();
            //linearSystems[0] = new SkylineLinearSystem(0, model.PatchesDictionary[0].Forces);
            //SolverSkyline solver = new SolverSkyline(linearSystems[0]);
            //ProblemStructural provider = new ProblemStructural(model, linearSystems);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, linearSystems);

            //parentAnalyzer.BuildMatrices();
            //parentAnalyzer.Initialize();
            //parentAnalyzer.Solve();


            //for (int i = 0; i < expectedSolutionVector.Length; i++)
            //{
            //	Assert.True(Utilities.AreValuesEqual(expectedSolutionVector[i], linearSystems[0].Solution[i],
            //		1e-6));
            //}

            return model;
		}

		
	}
}
