using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.LineLoads;
using ISAAR.MSolve.IGA.Loading.LoadElementFactories;
using ISAAR.MSolve.IGA.Loading.NodalLoads;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Newtonsoft.Json;
using Xunit;
using MatlabWriter = ISAAR.MSolve.LinearAlgebra.Output.MatlabWriter;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsKirchhoffLoveShellsDevelopDefGrad
	{


        [Fact]
        public void IsogeometricCantileverShell()
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "CantileverShell.txt");
            //var material = new ShellElasticMaterial2Dtransformationb()
            //{
            //    YoungModulus = 100,
            //    PoissonRatio = 0
            //};
            //var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filename, material: material);

            var material = new ShellElasticMaterial2DtransformationbDefGrad()
            {
                YoungModulus = 100,
                PoissonRatio = 0
            };
            var modelReader = new IsogeometricShellReader(GeometricalFormulation.DefGrad, filename, defGradMaterial: material);


            var model = modelReader.GenerateModelFromFile();

            //Value verticalDistributedLoad = delegate (double x, double y, double z)
            //{
            //    return new double[] { 0, 0, 4 };
            //};
            //model.Patches[0].cr.EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));
            //var rightEdgeLoads = model.Patches[0].
            //    .CreateLoadForEdge(model, NurbsSurfaceEdges.Right, new DistributedLineLoad2D(-0.3, 0));

            

            for (int i = 0; i < 6; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 20); //1000);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var logger = new TotalLoadsDisplacementsPerIncrementLogIGA(model.PatchesDictionary[0], 20, //1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "CantileverBenchmarkLog16x1.txt");
            childAnalyzer.IncrementalLogs.Add(0, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            model.Patches[0].Forces[11] = 1.3333333333333344;//..Forces[2] = 1.33333333334;
            model.Patches[0].Forces[14] = 1.3333333333333335;//..Forces[5] = 1.33333333334;
            model.Patches[0].Forces[17] = 1.3333333333333344;//..Forces[8] = 1.33333333334;

            parentAnalyzer.Solve();
        }

        [Fact]
		public void IsogeometricCantileverShellDimitriRefactoring1()
		{
			string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "CantileverShell.txt");
            var material = new ShellElasticSectionMaterial2D()
            {
                YoungModulus = 100,
                PoissonRatio = 0
            };
			var modelReader = new IsogeometricShellReader(GeometricalFormulation.Linear,filename,sectionMaterial: material);
			var model=modelReader.GenerateModelFromFile();

			model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[9],StructuralDof.TranslationZ,-1));
            model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[10],StructuralDof.TranslationZ,-1));
            model.Loads.Add(new NodalLoad(model.ControlPoints.ToList()[11],StructuralDof.TranslationZ,-1));

			for (int i = 0; i < 6; i++)
			{
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationX});
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationY});
				model.ControlPointsDictionary[i].Constraints.Add(new Constraint() {DOF = StructuralDof.TranslationZ});
			}

			// Solvers
			var solverBuilder = new SkylineSolver.Builder();
			ISolver solver = solverBuilder.BuildSolver(model);

			// Structural problem provider
			var provider = new ProblemStructural(model, solver);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer(model, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			var expectedSolution = new double[18]
			{
				0.0, 0.0, -7499.999986865148, 0.0, 0.0, -7499.99998660616, 0.0, 0.0, -7499.999986347174, 0.0, 0.0,
				-14999.999980230163, 0.0, 0.0, -14999.999980050825, 0.0, 0.0, -14999.999979871487
			};
			for (int i = 0; i < expectedSolution.Length; i++)
				Utilities.AreValuesEqual(expectedSolution[i], solver.LinearSystems[0].Solution[i], 7);
		}

        [Fact]
        public static void ScordelisLoShellDimitriRefactoring1()
        {
            var filename = "ScordelisLoShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles", $"{filename}.txt")
                .ToString(CultureInfo.InvariantCulture);
			var material = new ShellElasticMaterial2Dtransformationb()
            {
                YoungModulus = 4.3210e8,
                PoissonRatio = 0.0
            };
            var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filepath, material);
            var model = modelReader.GenerateModelFromFile();

			//model.SurfaceLoads.Add(new SurfaceDistributedLoad(-90, StructuralDof.TranslationY));

            // Rigid diaphragm for AB
            for (var i = 0; i < 19; i++)
            {
                model.ControlPointsDictionary[i * 19].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i * 19].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
            }

            // Symmetry for CD
            for (var i = 0; i < 19; i++)
            {
                model.ControlPointsDictionary[i * 19 + 18].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });

                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationX),
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationX)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 18], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[i * 19 + 17], StructuralDof.TranslationY)));
            }

            // Symmetry for AD
            for (var j = 0; j < 19; j++)
            {
                model.ControlPointsDictionary[j].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationY)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[j], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[j + 19], StructuralDof.TranslationZ)));
            }

            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            var cp = model.ControlPointsDictionary.Values.Last();
            var dofA = model.GlobalDofOrdering.GlobalFreeDofs[cp, StructuralDof.TranslationY];

            var solution = solver.LinearSystems[0].Solution[dofA];
        }

        


    }
}