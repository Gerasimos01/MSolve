﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Running;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using MathNet.Numerics.Data.Matlab;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace ISAAR.MSolve.IGA.Tests
{
    public class NurbsNonLinearShellsSemicylinder
    {
        private const double Tolerance = 1e-9;

        #region Fields
        private List<ControlPoint> ElementControlPoints()
        {
            return new List<ControlPoint>
            {
                new ControlPoint {ID = 0 , X = 1.016       , Y =  0		     , Z = 0		  , WeightFactor = 1.0},
                new ControlPoint {ID = 1 , X = 1.016002552 , Y =  0		     , Z = 0.018343973, WeightFactor = 1.0},
                new ControlPoint {ID = 2 , X = 1.015007727 , Y =  0		     , Z = 0.055031929, WeightFactor = 1.0},
                new ControlPoint {ID = 3 , X = 1.010535881 , Y =  0		     , Z = 0.109902059, WeightFactor = 1.0},
                new ControlPoint {ID = 32, X = 1.016       , Y =  0.035034483, Z = 0	      , WeightFactor = 1.0},
                new ControlPoint {ID = 33, X = 1.016002552 , Y =  0.035034483, Z = 0.018343973, WeightFactor = 1.0},
                new ControlPoint {ID = 34, X = 1.015007727 , Y =  0.035034483, Z = 0.055031929, WeightFactor = 1.0},
                new ControlPoint {ID = 35, X = 1.010535881 , Y =  0.035034483, Z = 0.109902059, WeightFactor = 1.0},
                new ControlPoint {ID = 64, X = 1.016	   , Y =  0.105103448, Z = 0	      , WeightFactor = 1.0},
                new ControlPoint {ID = 65, X = 1.016002552 , Y =  0.105103448, Z = 0.018343973, WeightFactor = 1.0},
                new ControlPoint {ID = 66, X = 1.015007727 , Y =  0.105103448, Z = 0.055031929, WeightFactor = 1.0},
                new ControlPoint {ID = 67, X = 1.010535881 , Y =  0.105103448, Z = 0.109902059, WeightFactor = 1.0},
                new ControlPoint {ID = 96, X = 1.016	   , Y =  0.210206897, Z = 0		  , WeightFactor = 1.0},
                new ControlPoint {ID = 97, X = 1.016002552 , Y =  0.210206897, Z = 0.018343973, WeightFactor = 1.0},
                new ControlPoint {ID = 98, X = 1.015007727 , Y =  0.210206897, Z = 0.055031929, WeightFactor = 1.0},
                new ControlPoint {ID = 99, X = 1.010535881 , Y =  0.210206897, Z = 0.109902059, WeightFactor = 1.0},
            };
        }

        private List<Knot> ElementKnots()
        {
            return new List<Knot>()
            {
                new Knot(){ID=0 ,Ksi=0.000,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=1 ,Ksi=0.000,Heta=1.0,Zeta =0.0 },
                new Knot(){ID=30,Ksi=1.000,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=31,Ksi=1.000,Heta=1.0,Zeta =0.0 }
            };
        }

        private Vector KnotValueVectorKsi()
        {
            return Vector.CreateFromArray(new double[]
            {
                0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29, 29, 29, 
            });
        }

        private Vector KnotValueVectorHeta()
        {
            return Vector.CreateFromArray(new double[]
            {
                0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29, 29, 29,
            });
        }

        private ShellElasticMaterial2D Material => new ShellElasticMaterial2D()
        {
            YoungModulus = 20685000,
            PoissonRatio = 0.3,
            //ta vectors einai tuxaia giati den tha xrhsimopoiithoun
            TangentVectorV1 = new double[] { 1.000000000000000000000000000000000000, 0.000000000000000053342746886286800000, 0.000000000000000000000000000000000000 },
            TangentVectorV2 = new double[] { 3.90312782094781000000000000000000E-18, 9.99999999999999000000000000000000E-01, 0.00000000000000000000000000000000E+00, },
            NormalVectorV3 = new double[] { 0, 0, 1 }
        };
        
        private NurbsKirchhoffLoveShellElementNLDevelop Element
        {
            get
            {
                var patch = new Patch();
                patch.DegreeKsi = 3;
                patch.DegreeHeta = 3;
                patch.NumberOfControlPointsHeta = 32;
                patch.KnotValueVectorKsi = KnotValueVectorKsi();
                patch.KnotValueVectorHeta = KnotValueVectorHeta();
                var element =
                    new NurbsKirchhoffLoveShellElementNLDevelop(Material, ElementKnots(), ElementControlPoints(), patch, 0.03);
                element._solution= localSolution;
                return element;
            }
        }

        private double[] localSolution => new double[]
        {
            0.8098531, 0.8075662, 0, 0.8100202, 0.8086463, 0, 8.129682, 0.01817363, 4.367361, 0.5137918,
            2.247061, -0.8942891, 0.4346529, -0.1461345, 0, 0.4353087, -0.1463096, 0, 6.092152, 4.127695, 
           -5.340579, 11.74476, -2.574487, 1.302691, 1.265022, -0.009938874, 0, 1.266052, -0.008793021, 
            0, 8.890049, -2.830236, 0.6063107, -2.602719, -0.07121834, 0.76033, -1.600903, 0.1575092,
            0, -1.601199, 0.1575786, 0, -6.643739, 0.8328065, 1.065195, 3.57901, 0.9146189, 0.09372011,

        };

        private Forces MembraneForces => new Forces()
        {
            v0 = -0.0327209647710278000000000000000000000000000000000,
            v1 = -0.0000000000000000000000000000193609885449360000000,
            v2 = 0.0000000000001690940892323080000000000000000000000
        };

        private Forces BendingMoments => new Forces()
        {
            v0 = -0.42618363401210900000000000000000000000000000000000,
            v1 = -0.00000000000000000075669834331405100000000000000000,
            v2 = 0.00000000000000683985998006887000000000000000000000,
        };

            #endregion

        [Fact]
        public void IsogeometricCantileverShell()
        {
            Model model = new Model();
            var filename = "CantileverShellBenchmark16x1";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(),"InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { 0, 0, 4 };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

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
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 1000);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var logger = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "CantileverBenchmarkLog16x1.txt");
            childAnalyzer.IncrementalLogs.Add(0, logger);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void SlitAnnularPlate()
        {
            Model model = new Model();
            var filename = "SplitAnnularPlate";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { 0, 0, 0.8 };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));

            for (int i = 0; i < 20; i++)
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
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 500);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 500,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "SplitAnnularPlateWa.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void HemisphericalShell()
        {
            Model model = new Model();
            var filename = "PinchedHemisphere";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            model.Loads.Add(new Load()
            {
                Amount = -200,
                Node = model.ControlPoints.ToList()[0],
                DOF = StructuralDof.TranslationX
            });

            model.Loads.Add(new Load()
            {
                Amount = 200,
                Node = model.ControlPoints.ToList()[240],
                DOF = StructuralDof.TranslationY
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            for (int i = 0; i < 16; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });

                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i], StructuralDof.TranslationX),
                    new NodalDof(model.ControlPointsDictionary[i + 16], StructuralDof.TranslationX)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[i + 16], StructuralDof.TranslationZ)));
            }

            for (int i = 256-16; i < 256; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });

                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[i - 16], StructuralDof.TranslationY)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[i], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[i - 16], StructuralDof.TranslationZ)));
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 20);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 20,
                model.ControlPoints.ToList()[0], StructuralDof.TranslationX, "PinchedHemisphereNegativeLoadNode.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 20,
                model.ControlPoints.ToList()[240], StructuralDof.TranslationY, "PinchedHemispherePositiveLoadNode.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void PulloutCylinderShell()
        {
            Model model = new Model();
            var filename = "PulloutCylinderShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            model.Loads.Add(new Load()
            {
                Amount = -40000,
                Node = model.ControlPoints.ToList().Last(),
                DOF = StructuralDof.TranslationZ
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            //TODO:Check boundary conditions
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[1].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationY});
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[2].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationZ});
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[3].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationX});
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 500);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            //var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 500,
            //    model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "SplitAnnularPlateWa.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            //childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void PinchedCylinderShell()
        {
            Model model = new Model();
            var filename = "PinchedCylinderShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            model.Loads.Add(new Load()
            {
                Amount = -12000,
                Node = model.ControlPoints.ToList().Last(),
                DOF = StructuralDof.TranslationZ
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            //TODO:Check boundary conditions
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[1].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationY});
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[2].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationZ});
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[3].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint(){DOF = StructuralDof.TranslationX});
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 500);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            //var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 500,
            //    model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "SplitAnnularPlateWa.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            //childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }


        [Fact]
        public void PinchedSemiCylinderShell()
        {
            Model model = new Model();
            var filename = "PinchedSemiCylindricalShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            model.Loads.Add(new Load()
            {
                Amount = -1000,
                Node = model.ControlPoints.ToList()[31],
                DOF = StructuralDof.TranslationZ
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            //TODO:Check boundary conditions
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[1].ControlPointsDictionary)
            {
                var id = controlPoint.Value.ID;
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });

                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            //TODO: constrain rotation
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[2].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id + 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationX),
                    new NodalDof(model.ControlPointsDictionary[id + 1], StructuralDof.TranslationX)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[id + 1], StructuralDof.TranslationY)));
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[3].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });

                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id - 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[id - 1], StructuralDof.TranslationZ)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[id - 1], StructuralDof.TranslationY)));
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 100);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 100,
                model.ControlPoints.ToList()[31], StructuralDof.TranslationZ, "PinchedSemiCylinderShell.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            var paraview = new ParaviewNurbsShells(model, childAnalyzer.uPlusdu[0], filename);
            //paraview.CreateParaviewFile();
        }

        [Fact]
        public void PinchedSemiCylinderShellCoarseAndCluster5()
        {
            Model model = new Model();
            var filename = "PinchedSemiCylindricalShell8x8Scaled";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            //TODO:Find load from previous papers
            model.Loads.Add(new Load()
            {
                Amount = -1000,
                Node = model.ControlPoints.ToList()[7],
                DOF = StructuralDof.TranslationZ
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            //TODO:Check boundary conditions
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[1].ControlPointsDictionary)
            {
                var id = controlPoint.Value.ID;
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });

                model.ControlPointsDictionary[id - 8].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[id - 8].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[id - 8].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            //TODO: constrain rotation
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[2].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id + 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationX),
                    new NodalDof(model.ControlPointsDictionary[id + 1], StructuralDof.TranslationX)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[id + 1], StructuralDof.TranslationY)));
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[3].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });

                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id - 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationZ),
                    new NodalDof(model.ControlPointsDictionary[id - 1], StructuralDof.TranslationZ)));
                model.AddPenaltyConstrainedDofPair(new PenaltyDofPair(
                    new NodalDof(model.ControlPointsDictionary[id], StructuralDof.TranslationY),
                    new NodalDof(model.ControlPointsDictionary[id - 1], StructuralDof.TranslationY)));
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 100);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 100,
                model.ControlPoints.ToList()[7], StructuralDof.TranslationZ, "PinchedSemiCylindricalShell8x8Scaled.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            var paraview = new ParaviewNurbsShells(model, childAnalyzer.uPlusdu[0], filename);
            //paraview.CreateParaviewFile();
        }


        [Fact]
        public void PinchedSemiCylinderShellNoPenalty()
        {
            Model model = new Model();
            var filename = "PinchedSemiCylindricalShell";
            var filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            //TODO:Find load from previous papers
            model.Loads.Add(new Load()
            {
                Amount = -2000,
                Node = model.ControlPoints.ToList()[31],
                DOF = StructuralDof.TranslationZ
            });

            //TODO: Possibly the tangent should also be fixes due to symmetry
            //TODO:Check boundary conditions
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[1].ControlPointsDictionary)
            {
                var id = controlPoint.Value.ID;
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });

                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[id - 32].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            //TODO: constrain rotation
            foreach (var controlPoint in model.Patches[0].EdgesDictionary[2].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id + 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            foreach (var controlPoint in model.Patches[0].EdgesDictionary[3].ControlPointsDictionary)
            {
                controlPoint.Value.Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });

                var id = controlPoint.Value.ID;
                model.ControlPointsDictionary[id - 1].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
            }

            // Solvers
            var solverBuilder = new SuiteSparseSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 100);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 100,
                model.ControlPoints.ToList()[31], StructuralDof.TranslationZ, "PinchedSemiCylinderShellNoPenaltyPointA.txt");
            //var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
            //    model.ControlPointsDictionary[790], StructuralDof.TranslationZ, "SplitAnnularPlateWb.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            //childAnalyzer.IncrementalLogs.Add(1, loggerB);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        [Fact]
        public void IsogeometricSquareShell10x10Straight()
        {
            Model model = new Model();
            var filename = "SquareShell10x10Straight";
            string filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

            for (int i = 0; i < 20; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            double load_factor = 300;
            int increments = 30;
            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { 0, 0, load_factor };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));


            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], increments,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "squarePlateZ.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationX, "squareplateX.txt");
            var loggerC = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationY, "squareplateY.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            childAnalyzer.IncrementalLogs.Add(1, loggerB);
            childAnalyzer.IncrementalLogs.Add(2, loggerC);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //var paraview = new ParaviewNurbsShells(model, solver.LinearSystems[0].Solution, filename);
            //paraview.CreateParaview2DFile();

            //var a = solver.LinearSystems[0].Solution;
        }

        [Fact]// exei ginei develop kai efelkusmos
        public void IsogeometricSquareShell10x10Straight30Degrees()
        {
            Model model = new Model();
            var filename = "SquareShell10x10Straight30Degrees";
            string filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            for (int i = 0; i < 20; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            double load_factor = 300;
            int increments = 30;

            double[] loadValues = new double[] { load_factor, 0, 0, 0, 0, 0 };
            double rot_phi_1 = 0;
            double rot_phi_2 = Math.PI * 2 * 30 / 360;
            double[] ekk_xyz = new double[] { 0, 0, 0 };
            loadValues = SupportiveMethods.modify_ox_sunol_forRotationAndTranslation(loadValues, rot_phi_1, rot_phi_2, ekk_xyz);

            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                //return new double[] { 0, 0, load_factor };
                return new double[] { loadValues[0], loadValues[1], loadValues[2] };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));


            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], increments,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "squarePlate30DegreesZ.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationX, "squareplate30DegreesX.txt");
            var loggerC = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationY, "squareplate30DegreesY.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            childAnalyzer.IncrementalLogs.Add(1, loggerB);
            childAnalyzer.IncrementalLogs.Add(2, loggerC);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            var paraview = new ParaviewNurbsShells(model, childAnalyzer.uPlusdu[0], filename);
            paraview.CreateParaview2DFile();

            //var a = solver.LinearSystems[0].Solution;
        }

        [Fact]
        public void IsogeometricSquareShell10x10StraightDevelop()
        {
            Model model = new Model();
            var filename = "SquareShell10x10Straight";
            string filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            for (int i = 0; i < 20; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            double load_factor = 300;
            int increments = 30;
            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { load_factor, 0, 0};
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));


            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], increments,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "squarePlateZDevelopTension.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationX, "squareplateXDevelopTension.txt");
            var loggerC = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationY, "squareplateYDevelopTension.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            childAnalyzer.IncrementalLogs.Add(1, loggerB);
            childAnalyzer.IncrementalLogs.Add(2, loggerC);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //var paraview = new ParaviewNurbsShells(model, solver.LinearSystems[0].Solution, filename);
            //paraview.CreateParaview2DFile();

            //var a = solver.LinearSystems[0].Solution;
        }

        [Fact]
        public void IsogeometricSquareShell10x10DeformedDevelop()
        {
            Model model = new Model();
            var filename = "SquareShell10x10Deformed";
            string filepath = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", $"{filename}.txt");
            IsogeometricShellReader modelReader = new IsogeometricShellReader(model, filepath);
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinearDevelop);

            for (int i = 0; i < 20; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            double load_factor = 300;
            int increments = 30;
            Value verticalDistributedLoad = delegate (double x, double y, double z)
            {
                return new double[] { load_factor, 0, 0 };
            };
            model.Patches[0].EdgesDictionary[1].LoadingConditions.Add(new NeumannBoundaryCondition(verticalDistributedLoad));


            // Solvers
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], increments,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "squarePlateZDevelopDeformedTension.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationX, "squareplateXDevelopDeformedTension.txt");
            var loggerC = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationY, "squareplateYDevelopDeformedTension.txt");
            childAnalyzer.IncrementalLogs.Add(0, loggerA);
            childAnalyzer.IncrementalLogs.Add(1, loggerB);
            childAnalyzer.IncrementalLogs.Add(2, loggerC);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //var paraview = new ParaviewNurbsShells(model, solver.LinearSystems[0].Solution, filename);
            //paraview.CreateParaview2DFile();

            //var a = solver.LinearSystems[0].Solution;
        }

        [Fact]
        public void JacobianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = new double[2, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var expectedJacobian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "jacobian");
            
            for (var i = 0; i < jacobianMatrix.GetLength(0); i++)
            {
                for (var j = 0; j < jacobianMatrix.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedJacobian[i, j], jacobianMatrix[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void HessianTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var expectedHessian = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "hessian");

            for (var i = 0; i < hessian.GetLength(0); i++)
            {
                for (var j = 0; j < hessian.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedHessian[i, j], hessian[i, j],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void MembraneDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);

            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;
            var Bmembrane = new double[3, controlPoints.Length * 3];
            shellElement.CalculateMembraneDeformationMatrix(elementControlPoints.Length, nurbs, 0, surfaceBasisVector1, surfaceBasisVector2, Bmembrane);

            var expectedBmembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bmembrane");

            for (var i = 0; i < Bmembrane.GetLength(0); i++)
            {
                for (var j = 0; j < Bmembrane.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBmembrane[i, j], Bmembrane[i, j],
                        Tolerance));
                }
            }
        }
        
        [Fact]
        public void BendingDeformationMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);

            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector1(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector1(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector1(hessian, 2);
            var Bbending = new double[3, controlPoints.Length * 3];
            shellElement.CalculateBendingDeformationMatrix(elementControlPoints.Length, surfaceBasisVector3, nurbs, 0, surfaceBasisVector2,
                surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, Bbending);

            var expectedBbending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Bbending");

            for (var i = 0; i < Bbending.GetLength(0); i++)
            {
                for (var j = 0; j < Bbending.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedBbending[i, j], Bbending[i, j],
                        Tolerance));
                }
            }
        }
        
        
        [Fact]
        public void StiffnessMatrixMembraneNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var KmembraneNL = new double[elementControlPoints.Length * 3, elementControlPoints.Length * 3];

            var membraneForces = new Forces()
            {
                v0 = -0.0327209647710278000000000000000000000000000000000,
                v1 = -0.0000000000000000000000000000193609885449360000000,
                v2 = 0.0000000000001690940892323080000000000000000000000
            };

            shellElement.CalculateKmembraneNL(elementControlPoints, ref membraneForces, nurbs, 0, KmembraneNL);

            var expectedKmembraneNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KmembraneNL");

            for (var i = 0; i < KmembraneNL.GetLength(0); i++)
            {
                for (var j = 0; j < KmembraneNL.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKmembraneNL[i, j], KmembraneNL[i, j], Tolerance));
                }
            }
        }

        [Fact]
        public void StiffnessMatrixBendingNLTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            shellElement._solution = localSolution;
            var elementControlPoints = shellElement.CurrentControlPoint(controlPoints);
            var jacobianMatrix = new double[3, 3];
            shellElement.CalculateJacobian(elementControlPoints, nurbs, 0, jacobianMatrix);
            var hessian = shellElement.CalculateHessian(elementControlPoints, nurbs, 0);

            var surfaceBasisVector1 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 0);
            var surfaceBasisVector2 = shellElement.CalculateSurfaceBasisVector1(jacobianMatrix, 1);
            var surfaceBasisVector3 = new[]
            {
                surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
            };
            double norm = surfaceBasisVector3.Sum(t => t * t);
            var J1 = Math.Sqrt(norm);

            for (int i = 0; i < surfaceBasisVector3.Length; i++)
                surfaceBasisVector3[i] = surfaceBasisVector3[i] / J1;

            var surfaceBasisVectorDerivative1 = shellElement.CalculateSurfaceBasisVector1(hessian, 0);
            var surfaceBasisVectorDerivative2 = shellElement.CalculateSurfaceBasisVector1(hessian, 1);
            var surfaceBasisVectorDerivative12 = shellElement.CalculateSurfaceBasisVector1(hessian, 2);

            var KbendingNL = new double[elementControlPoints.Length * 3, elementControlPoints.Length * 3];

            var BendingMoments = new Forces()
            {
                v0 = -0.42618363401210900000000000000000000000000000000000,
                v1 = -0.00000000000000000075669834331405100000000000000000,
                v2 = 0.00000000000000683985998006887000000000000000000000,
            };

            shellElement.CalculateKbendingNL(elementControlPoints, ref BendingMoments, nurbs,
                surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3,
                surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, J1, 0, KbendingNL);

            var expectedKbendingNL = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "KbendingNL");


            for (var i = 0; i < KbendingNL.GetLength(0); i++)
            {
                for (var j = 0; j < KbendingNL.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKbendingNL[i, j], KbendingNL[i, j],1e-04));
                }
            }
        }
        
        [Fact]
        public void ConstitutiveMatrixThicknessIntegration()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
                shellElement.IntegratedConstitutiveOverThickness(gaussPoints[0]);

            var expectedConstitutiveMembrane = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "MembraneConstitutiveMatrix");
            var expectedConstitutiveBending = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "BendingConstitutiveMatrix");

            for (int i = 0; i < MembraneConstitutiveMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < MembraneConstitutiveMatrix.GetLength(1); j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveMembrane[i, j], MembraneConstitutiveMatrix[i, j], Tolerance));
                    Assert.True(Utilities.AreValuesEqual(expectedConstitutiveBending[i, j], BendingConstitutiveMatrix[i, j], Tolerance));
                }
            }
        }

        [Fact]
        public void StressesThicknessIntegration()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);

            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var MembraneForces = new Forces();
            var BendingMoments = new Forces();
            shellElement.IntegratedStressesOverThickness(gaussPoints[0], ref MembraneForces, ref BendingMoments);

            var expectedMembraneForces = new double[3]
            {
                -0.03272096477102780000000000000000000000000000000,
                -0.00000000000000000000000000001936098854493600000,
                0.00000000000016909408923230800000000000000000000,
            };

            var expectedBendingMoments = new double[3]
            {
                -0.426183634012109000000000000000000000000,
                -0.000000000000000000756698343314051000000,
                0.000000000000006839859980068870000000000,

            };

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[0], MembraneForces.v0, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[0], BendingMoments.v0, Tolerance));

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[1], MembraneForces.v1, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[1], BendingMoments.v1, Tolerance));

            Assert.True(Utilities.AreValuesEqual(expectedMembraneForces[2], MembraneForces.v2, Tolerance));
            Assert.True(Utilities.AreValuesEqual(expectedBendingMoments[2], BendingMoments.v2, Tolerance));
        }

        [Fact]
        public void TotalElementStiffnessMatrixTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;
            
            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var Ktotal=shellElement.StiffnessMatrix(shellElement);

            var expectedKtotal = MatlabReader.Read<double>(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "NurbsNonLinearThicknessShell.mat"), "Ktotal");
            
            for (var i = 0; i < Ktotal.NumRows; i++)
            {
                for (var j = 0; j < Ktotal.NumColumns; j++)
                {
                    Assert.True(Utilities.AreValuesEqual(expectedKtotal[i, j], Ktotal[i, j],
                        1e-6));
                }
            }
        }

        [Fact]
        public void InternalForcesTestSermiCylinder()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var f = shellElement.CalculateForces(shellElement,localSolution, new double[27]);

            var expectedForces = new double[48]
            {
                -4475874511.6703205, -460505704.67915344, -521958924.80092776, -1636389571.34632, -454686747.3297115, 2596964237.0245924, 6212880155.002533, 
                30746666.69804227, 4359693512.940589, 1713850335.7118824, 28763849.307121132, 928636487.32193434, -6728006086.4731855, -1111253448.2442992,
                916806505.492374, -4552631231.625246, 138153058.40934575, -2499355211.8332977, 7970041652.2430191, 1760091337.3699565, -4318875872.3011389,
                2461338167.8724332, 398795166.24829733, -881449694.7845453, -3382002623.9320226, -430550944.3491326, 829847975.84043252, -2273984471.4064822,
                -233602885.72271252, 37058342.897073761, 3802334834.4444332, 339683015.47409087, -1232819191.862458, 1166886598.116899, 88098778.670311511, 
                -293671857.34619135, -464657651.00243622, -38888505.0231959, 118532858.48272517, -375335435.19537985, -52073310.575925149, 79684644.195540756,
                423734126.86863887, -4290467.937362656, -92428023.980797991, 137815712.39155602, 1520141.6843277283, -26665787.285910059,

            };

            for (var j = 0; j < expectedForces.Length; j++)
            {
                Assert.True(Utilities.AreValuesEqual(expectedForces[j], f[j],
                    Tolerance));
            }
        }

        [Fact]
        public void StiffnessMatrixTestSermiCylinderForNonlinearity()// TODO: metakinhseis tha perasoun sto calculateStresses kai tha xreiazomaste ekeino pia.
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);


            var forces = shellElement.CalculateForces(shellElement, localSolution, new double[27]);

            var Stiffness1 = shellElement.StiffnessMatrix(shellElement);

            var expectedStiffnessMatr1 = GetExpectedStiffnessMatrixForNonlinearityCheck();  //.

            for (var j = 0; j < expectedStiffnessMatr1.GetLength(0); j++)
            {
                for (int k = 0; k < expectedStiffnessMatr1.GetLength(0); k++)
                {


                    Assert.True(Utilities.AreValuesEqual(Stiffness1[j,k], expectedStiffnessMatr1[j,k],
                        Tolerance));
                }
            }
        }

        [Fact]
        public void InternalForcesTestSermiCylinderDevelop()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStressesDevelop(shellElement, localSolution, new double[27]);
            var f = shellElement.CalculateForcesDevelop1(shellElement, localSolution, new double[27]);

            var expectedForces = new double[48]
            {
                -4475874511.6703205, -460505704.67915344, -521958924.80092776, -1636389571.34632, -454686747.3297115, 2596964237.0245924, 6212880155.002533,
                30746666.69804227, 4359693512.940589, 1713850335.7118824, 28763849.307121132, 928636487.32193434, -6728006086.4731855, -1111253448.2442992,
                916806505.492374, -4552631231.625246, 138153058.40934575, -2499355211.8332977, 7970041652.2430191, 1760091337.3699565, -4318875872.3011389,
                2461338167.8724332, 398795166.24829733, -881449694.7845453, -3382002623.9320226, -430550944.3491326, 829847975.84043252, -2273984471.4064822,
                -233602885.72271252, 37058342.897073761, 3802334834.4444332, 339683015.47409087, -1232819191.862458, 1166886598.116899, 88098778.670311511,
                -293671857.34619135, -464657651.00243622, -38888505.0231959, 118532858.48272517, -375335435.19537985, -52073310.575925149, 79684644.195540756,
                423734126.86863887, -4290467.937362656, -92428023.980797991, 137815712.39155602, 1520141.6843277283, -26665787.285910059,

            };

            for (var j = 0; j < expectedForces.Length; j++)
            {
                Assert.True(Utilities.AreValuesEqual(expectedForces[j], f[j],
                    Tolerance));
            }
        }

        [Fact]
        public void StiffnessMatrixTestSermiCylinderForNonlinearityDevelop()// TODO: metakinhseis tha perasoun sto calculateStresses kai tha xreiazomaste ekeino pia.
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStressesDevelop(shellElement, localSolution, new double[27]);


            

            var Stiffness1 = shellElement.StiffnessMatrixDevelop(shellElement);

            var expectedStiffnessMatr1 = GetExpectedStiffnessMatrixForNonlinearityCheck();  //.

            for (var j = 0; j < expectedStiffnessMatr1.GetLength(0); j++)
            {
                for (int k = 0; k < expectedStiffnessMatr1.GetLength(0); k++)
                {


                    Assert.True(Utilities.AreValuesEqual(Stiffness1[j, k], expectedStiffnessMatr1[j, k],
                        Tolerance));
                }
            }
        }


        [Fact]
        public void TestA3r()
        {
            var surfaceBasisVector1 = new double[] { 0.49281343102090647, -0.0027907619686642761, -2.0647413184286493E-07 };
            var surfaceBasisVector2 = new double[] {0.0032374497982368571, 0.57141938776238976, 2.3218349204021492E-09};
            var surfaceBasisVector3 = new double[] { 4.1893372881671206E-07, -6.4367991618857595E-09, 0.99999999999991229};

            double dksi_r=-2.0934480533670992;
            double dheta_r=-2.0934480533670996;
            double J1=0.28161218398684618;
            var a3r= new a3r();

            NurbsKirchhoffLoveShellElementNL.CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r, dheta_r, J1, ref a3r);
             
            Assert.Equal(1.7882446740031873E-06,a3r.a3r00);
            Assert.Equal(-2.7475877512613432E-08,a3r.a3r01);
            Assert.Equal(4.2685621877569231,a3r.a3r02);
                         
            Assert.Equal(1.5246711353875095E-06,a3r.a3r10);
            Assert.Equal(-2.3426144068498867E-08,a3r.a3r11);
            Assert.Equal(3.639408886207,a3r.a3r12);

            Assert.Equal(-7.38802532863679E-13,a3r.a3r20);
            Assert.Equal(1.1827148338265939E-14,a3r.a3r21);
            Assert.Equal(-1.7648185299346879E-06,a3r.a3r22);
        }

        [Fact]
        public void TestCalculateCrossProduct()
        {
            var v1 = new double[] {1, 20, 50};
            var v2 = new double[] {4, 30, 45};
            var c = new double[3];
            NurbsKirchhoffLoveShellElementNL.CalculateCrossProduct(v1, v2, c);

            var vector1 = Vector.CreateFromArray(v1);
            var vector2 = Vector.CreateFromArray(v2);

            var cross = vector1.CrossProduct(vector2);

            Assert.Equal(cross[0],c[0]);
            Assert.Equal(cross[1],c[1]);
            Assert.Equal(cross[2],c[2]);
        }

        [Fact]
        public void TestBab_rs()
        {
            #region Initialization
            var surfaceBasisVectorDerivative1 = new double[3]
            {
                -0.0002215697315227748,
                -0.040236627065225024,
                -3.1283539134440684E-06
            };
            var surfaceBasisVectorDerivative2 = new double[3]
            {
                4.3550534925664095E-11,
                1.0323258744575857E-08,
                -1.0202863489210098E-09
            };
            var surfaceBasisVectorDerivative12 = new double[3]
            {
                0.046627259953136886,
                -0.00026442650364636821,
                6.5202798291920747E-08
            };

            var d2ksi_dr2 = 4.4992901171737891;
            var d2ksi_ds2 = 4.4992901171737891;

            var d2heta_dr2 = 4.49929011717379;
            var d2heta_ds2 = 4.49929011717379;

            var d2ksiheta_dr2 = 6.7489351757606837;
            var d2ksiheta_ds2 = 6.7489351757606837;
            var a3rs= new a3rs();
            a3rs.a3rs00_0=	1.5266467195814444E-05	;
            a3rs.a3rs00_1=	1.3016307114560258E-05	;
            a3rs.a3rs00_2=	-1.1837641977763269E-11	;
            a3rs.a3rs01_0=	6.390871065454354E-06	;
            a3rs.a3rs01_1=	5.448905725897356E-06	;
            a3rs.a3rs01_2=	-2.55440113505756E-12	;
            a3rs.a3rs02_0=	18.220623150741094	;
            a3rs.a3rs02_1=	15.535043157448497	;
            a3rs.a3rs02_2=	-2.0715372921711805E-05	;
            a3rs.a3rs10_0=	6.390871065454354E-06	;
            a3rs.a3rs10_1=	5.448905725897356E-06	;
            a3rs.a3rs10_2=	-2.55440113505756E-12	;
            a3rs.a3rs11_0=	-1.9999190555149591E-07	;
            a3rs.a3rs11_1=	-1.7051463378493538E-07	;
            a3rs.a3rs11_2=	8.8817841970012523E-14	;
            a3rs.a3rs12_0=	15.53504315745124	;
            a3rs.a3rs12_1=	13.245297041003683	;
            a3rs.a3rs12_2=	-6.22035643166942E-06	;
            a3rs.a3rs20_0=	18.220623150741094	;
            a3rs.a3rs20_1=	15.535043157448497	;
            a3rs.a3rs20_2=	-2.0715372921711805E-05	;
            a3rs.a3rs21_0=	15.53504315745124	;
            a3rs.a3rs21_1=	13.245297041003683	;
            a3rs.a3rs21_2=	-6.22035643166942E-06	;
            a3rs.a3rs22_0=	-2.8248610566845741E-05	;
            a3rs.a3rs22_1=	-1.2643252672057042E-05	;
            a3rs.a3rs22_2=	-31.465920191744779	;

            var a3r= new a3r();
            a3r.a3r00=	1.7882446740031873E-06  ;
            a3r.a3r01=	-2.7475877512613432E-08	;
            a3r.a3r02=	4.2685621877569231	;
            a3r.a3r10=	1.5246711353875095E-06	;
            a3r.a3r11=	-2.3426144068498867E-08	;
            a3r.a3r12=	3.639408886207	;
            a3r.a3r20=	-7.38802532863679E-13	;
            a3r.a3r21=	1.1827148338265939E-14	;
            a3r.a3r22=	-1.7648185299346879E-06	;

            var a3s= new a3r();
            a3s.a3r00=1.7882446740031873E-06	;
            a3s.a3r01=-2.7475877512613432E-08	;
            a3s.a3r02=4.2685621877569231	;
            a3s.a3r10=1.5246711353875095E-06	;
            a3s.a3r11=-2.3426144068498867E-08	;
            a3s.a3r12=3.639408886207	;
            a3s.a3r20=-7.38802532863679E-13	;
            a3s.a3r21=1.1827148338265939E-14	;
            a3s.a3r22=	-1.7648185299346879E-06	;
            
            var da3_drds = new double[3,3][];
            da3_drds[0, 0] = new double[]
            {
                1.5266467195814444E-05,
                1.3016307114560258E-05,
                -1.1837641977763269E-11,
            };
            da3_drds[0,1] = new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            };
            da3_drds[0, 2] = new double[]
            {
                18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05
            };

            da3_drds[1, 0] = new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            };
            da3_drds[1, 1] = new double[]
            {
                -1.9999190555149591E-07,
                -1.7051463378493538E-07,
                8.8817841970012523E-14
            };
            da3_drds[1, 2] = new double[]
            {
                15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06
            };

            da3_drds[2, 0] = new double[]
            {
                18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05
            };
            da3_drds[2, 1] = new double[]
            {
                15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06
            };
            da3_drds[2, 2] = new double[]
            {
                -2.8248610566845741E-05,
                -1.2643252672057042E-05,
                -31.465920191744779
            };

            Bab_rs Bab_rsExpected= new Bab_rs();
            Bab_rsExpected.Bab_rs00_0=	1.5564548295526568E-05	;
            Bab_rsExpected.Bab_rs00_1=	1.6091663312697994E-05	;
            Bab_rsExpected.Bab_rs00_2=	4.9691772888834864E-05	;
            Bab_rsExpected.Bab_rs01_0=	6.5156542160513032E-06	;
            Bab_rsExpected.Bab_rs01_1=	6.7363158837647762E-06	;
            Bab_rsExpected.Bab_rs01_2=	2.0802043424519865E-05	;
            Bab_rsExpected.Bab_rs02_0=	18.576384789429813	;
            Bab_rsExpected.Bab_rs02_1=	19.205499827078938	;
            Bab_rsExpected.Bab_rs02_2=	59.307438707759943	;
            Bab_rsExpected.Bab_rs10_0=	6.5156542160513032E-06	;
            Bab_rsExpected.Bab_rs10_1=	6.7363158837647762E-06	;
            Bab_rsExpected.Bab_rs10_2=	2.0802043424519861E-05	;
            Bab_rsExpected.Bab_rs11_0=	-2.0389679110046287E-07	;
            Bab_rsExpected.Bab_rs11_1=	-2.1080203875074925E-07	;
            Bab_rsExpected.Bab_rs11_2=	-6.5096608290578748E-07	;
            Bab_rsExpected.Bab_rs12_0=	15.838368261336552	;
            Bab_rsExpected.Bab_rs12_1=	16.374756571476876	;
            Bab_rsExpected.Bab_rs12_2=	50.565977478394949	;
            Bab_rsExpected.Bab_rs20_0=	18.576384789429813	;
            Bab_rsExpected.Bab_rs20_1=	19.205499827078938	;
            Bab_rsExpected.Bab_rs20_2=	59.307438707759943	;
            Bab_rsExpected.Bab_rs21_0=	15.838368261336552	;
            Bab_rsExpected.Bab_rs21_1=	16.374756571476876	;
            Bab_rsExpected.Bab_rs21_2=	50.565977478394956	;
            Bab_rsExpected.Bab_rs22_0=	8.3070654310999034E-05	;
            Bab_rsExpected.Bab_rs22_1=	-1.5848757023602572E-05	;
            Bab_rsExpected.Bab_rs22_2=	-5.4373539710938839E-05	;

            #endregion
            
            Bab_rs Bab_rsAlternative = NurbsKirchhoffLoveShellElementNL.CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12,d2ksi_dr2,d2ksi_ds2, d2heta_dr2,d2heta_ds2, d2ksiheta_dr2,d2ksiheta_ds2,
                a3rs, a3r, a3s, da3_drds);


            Assert.Equal(Bab_rsExpected.Bab_rs00_0,Bab_rsAlternative.Bab_rs00_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs00_1,Bab_rsAlternative.Bab_rs00_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs00_2,Bab_rsAlternative.Bab_rs00_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_0,Bab_rsAlternative.Bab_rs01_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_1,Bab_rsAlternative.Bab_rs01_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs01_2,Bab_rsAlternative.Bab_rs01_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_0,Bab_rsAlternative.Bab_rs02_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_1,Bab_rsAlternative.Bab_rs02_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs02_2,Bab_rsAlternative.Bab_rs02_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_0,Bab_rsAlternative.Bab_rs10_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_1,Bab_rsAlternative.Bab_rs10_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs10_2,Bab_rsAlternative.Bab_rs10_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_0,Bab_rsAlternative.Bab_rs11_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_1,Bab_rsAlternative.Bab_rs11_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs11_2,Bab_rsAlternative.Bab_rs11_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_0,Bab_rsAlternative.Bab_rs12_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_1,Bab_rsAlternative.Bab_rs12_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs12_2,Bab_rsAlternative.Bab_rs12_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_0,Bab_rsAlternative.Bab_rs20_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_1,Bab_rsAlternative.Bab_rs20_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs20_2,Bab_rsAlternative.Bab_rs20_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_0,Bab_rsAlternative.Bab_rs21_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_1,Bab_rsAlternative.Bab_rs21_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs21_2,Bab_rsAlternative.Bab_rs21_2  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_0,Bab_rsAlternative.Bab_rs22_0  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_1,Bab_rsAlternative.Bab_rs22_1  ,9);
            Assert.Equal(Bab_rsExpected.Bab_rs22_2,Bab_rsAlternative.Bab_rs22_2  ,9);
        }

        [Fact]
        public void Test_a3rs()
        {
            var surfaceBasisVector1 = Vector.CreateFromArray(new double[]{0.49281343102090647,
                                                                          -0.0027907619686642761,
                                                                          -2.0647413184286493E-07});
            var surfaceBasisVector2 =Vector.CreateFromArray(new double[]{0.0032374497982368571,
                0.57141938776238976,
                2.3218349204021492E-09});
            var surfaceBasisVector3 = Vector.CreateFromArray(new double[]{4.1893372881671206E-07,
                                                                          -6.4367991618857595E-09,
                0.99999999999991229});
            double J1 = 0.28161218398684618;
            var dksi_r = -2.0934480533670992;
            var dheta_r = -2.0934480533670996;

            var dksi_s = -2.0934480533670992;
            var dheta_s= -2.0934480533670996;

            #region Expected Values
            var a3rsExpected= new a3rs();
            a3rsExpected.a3rs00_0=	1.5266467195814444E-05	;
            a3rsExpected.a3rs00_1=	1.3016307114560258E-05	;
            a3rsExpected.a3rs00_2=	-1.1837641977763269E-11	;
            a3rsExpected.a3rs01_0=	6.390871065454354E-06	;
            a3rsExpected.a3rs01_1=	5.448905725897356E-06	;
            a3rsExpected.a3rs01_2=	-2.55440113505756E-12	;
            a3rsExpected.a3rs02_0=	18.220623150741094	;
            a3rsExpected.a3rs02_1=	15.535043157448497	;
            a3rsExpected.a3rs02_2=	-2.0715372921711805E-05	;
            a3rsExpected.a3rs10_0=	6.390871065454354E-06	;
            a3rsExpected.a3rs10_1=	5.448905725897356E-06	;
            a3rsExpected.a3rs10_2=	-2.55440113505756E-12	;
            a3rsExpected.a3rs11_0=	-1.9999190555149591E-07	;
            a3rsExpected.a3rs11_1=	-1.7051463378493538E-07	;
            a3rsExpected.a3rs11_2=	8.8817841970012523E-14	;
            a3rsExpected.a3rs12_0=	15.53504315745124	;
            a3rsExpected.a3rs12_1=	13.245297041003683	;
            a3rsExpected.a3rs12_2=	-6.22035643166942E-06	;
            a3rsExpected.a3rs20_0=	18.220623150741094	; 
            a3rsExpected.a3rs20_1=	15.535043157448497	;
            a3rsExpected.a3rs20_2=	-2.0715372921711805E-05	;
            a3rsExpected.a3rs21_0=	15.53504315745124	;
            a3rsExpected.a3rs21_1=	13.245297041003683	;
            a3rsExpected.a3rs21_2=	-6.22035643166942E-06	;
            a3rsExpected.a3rs22_0=	-2.8248610566845741E-05	;
            a3rsExpected.a3rs22_1=	-1.2643252672057042E-05	;
            a3rsExpected.a3rs22_2=	-31.465920191744779	;

            var da3_drdsExpected = new Vector[3, 3];
            da3_drdsExpected[0, 0] = Vector.CreateFromArray(new double[] { 1.5266467195814444E-05,
                1.3016307114560258E-05,
                -1.1837641977763269E-11 });
            da3_drdsExpected[0, 1] = Vector.CreateFromArray(new double[] { 6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12 });
            da3_drdsExpected[0, 2] = Vector.CreateFromArray(new double[] {18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05 });

            da3_drdsExpected[1, 0] = Vector.CreateFromArray(new double[]
            {
                6.390871065454354E-06,
                5.448905725897356E-06,
                -2.55440113505756E-12
            });
            da3_drdsExpected[1, 1] = Vector.CreateFromArray(new double[] { -1.9999190555149591E-07,
                                                                           -1.7051463378493538E-07,
                8.8817841970012523E-14 });
            da3_drdsExpected[1, 2] = Vector.CreateFromArray(new double[] { 15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06 });

            da3_drdsExpected[2, 0] = Vector.CreateFromArray(new double[] { 18.220623150741094,
                15.535043157448497,
                -2.0715372921711805E-05 });
            da3_drdsExpected[2, 1] = Vector.CreateFromArray(new double[] { 15.53504315745124,
                13.245297041003683,
                -6.22035643166942E-06 });
            da3_drdsExpected[2, 2] = Vector.CreateFromArray(new double[] { -2.8248610566845741E-05,
                                                                           -1.2643252672057042E-05,
                                                                           -31.465920191744779 });
            #endregion

            (a3rs a3rsAlternative, var da3tilde_drds, var da3tilde_dr, var da3tilde_ds,
                    double[] dnorma3_dr, double[] dnorma3_ds, double[,] dnorma3_drds, var a3_tilde, var da3_drds) =
                NurbsKirchhoffLoveShellElementNL.Calculate_a3rs(surfaceBasisVector1, surfaceBasisVector2,
                    surfaceBasisVector3, J1, dksi_r, dksi_s, dheta_r,dheta_s);
            
            #region Assertions
            Assert.Equal(a3rsExpected.a3rs00_0 ,a3rsAlternative.a3rs00_0,9);
            Assert.Equal(a3rsExpected.a3rs00_1 ,a3rsAlternative.a3rs00_1,9);
            Assert.Equal(a3rsExpected.a3rs00_2 ,a3rsAlternative.a3rs00_2,9);
            Assert.Equal(a3rsExpected.a3rs01_0 ,a3rsAlternative.a3rs01_0,9);
            Assert.Equal(a3rsExpected.a3rs01_1 ,a3rsAlternative.a3rs01_1,9);
            Assert.Equal(a3rsExpected.a3rs01_2 ,a3rsAlternative.a3rs01_2,9);
            Assert.Equal(a3rsExpected.a3rs02_0 ,a3rsAlternative.a3rs02_0,9);
            Assert.Equal(a3rsExpected.a3rs02_1 ,a3rsAlternative.a3rs02_1,9);
            Assert.Equal(a3rsExpected.a3rs02_2 ,a3rsAlternative.a3rs02_2,9);
            Assert.Equal(a3rsExpected.a3rs10_0 ,a3rsAlternative.a3rs10_0,9);
            Assert.Equal(a3rsExpected.a3rs10_1 ,a3rsAlternative.a3rs10_1,9);
            Assert.Equal(a3rsExpected.a3rs10_2 ,a3rsAlternative.a3rs10_2,9);
            Assert.Equal(a3rsExpected.a3rs11_0 ,a3rsAlternative.a3rs11_0,9);
            Assert.Equal(a3rsExpected.a3rs11_1 ,a3rsAlternative.a3rs11_1,9);
            Assert.Equal(a3rsExpected.a3rs11_2 ,a3rsAlternative.a3rs11_2,9);
            Assert.Equal(a3rsExpected.a3rs12_0 ,a3rsAlternative.a3rs12_0,9);
            Assert.Equal(a3rsExpected.a3rs12_1 ,a3rsAlternative.a3rs12_1,9);
            Assert.Equal(a3rsExpected.a3rs12_2 ,a3rsAlternative.a3rs12_2,9);
            Assert.Equal(a3rsExpected.a3rs20_0 ,a3rsAlternative.a3rs20_0,9);
            Assert.Equal(a3rsExpected.a3rs20_1 ,a3rsAlternative.a3rs20_1,9);
            Assert.Equal(a3rsExpected.a3rs20_2 ,a3rsAlternative.a3rs20_2,9);
            Assert.Equal(a3rsExpected.a3rs21_0 ,a3rsAlternative.a3rs21_0,9);
            Assert.Equal(a3rsExpected.a3rs21_1 ,a3rsAlternative.a3rs21_1,9);
            Assert.Equal(a3rsExpected.a3rs21_2 ,a3rsAlternative.a3rs21_2,9);
            Assert.Equal(a3rsExpected.a3rs22_0 ,a3rsAlternative.a3rs22_0,9);
            Assert.Equal(a3rsExpected.a3rs22_1 ,a3rsAlternative.a3rs22_1,9);
            Assert.Equal(a3rsExpected.a3rs22_2 ,a3rsAlternative.a3rs22_2,9);


            Assert.Equal(da3_drdsExpected[0, 0][0], da3_drds[0, 0][0],9);
            Assert.Equal(da3_drdsExpected[0, 0][1], da3_drds[0, 0][1],9);
            Assert.Equal(da3_drdsExpected[0, 0][2], da3_drds[0, 0][2],9);
            Assert.Equal(da3_drdsExpected[0, 1][0], da3_drds[0, 1][0],9);
            Assert.Equal(da3_drdsExpected[0, 1][1], da3_drds[0, 1][1],9);
            Assert.Equal(da3_drdsExpected[0, 1][2], da3_drds[0, 1][2],9);
            Assert.Equal(da3_drdsExpected[0, 2][0], da3_drds[0, 2][0],9);
            Assert.Equal(da3_drdsExpected[0, 2][1], da3_drds[0, 2][1],9);
            Assert.Equal(da3_drdsExpected[0, 2][2], da3_drds[0, 2][2],9);
                                                                     
            Assert.Equal(da3_drdsExpected[1, 0][0], da3_drds[1, 0][0],9);
            Assert.Equal(da3_drdsExpected[1, 0][1], da3_drds[1, 0][1],9);
            Assert.Equal(da3_drdsExpected[1, 0][2], da3_drds[1, 0][2],9);
            Assert.Equal(da3_drdsExpected[1, 1][0], da3_drds[1, 1][0],9);
            Assert.Equal(da3_drdsExpected[1, 1][1], da3_drds[1, 1][1],9);
            Assert.Equal(da3_drdsExpected[1, 1][2], da3_drds[1, 1][2],9);
            Assert.Equal(da3_drdsExpected[1, 2][0], da3_drds[1, 2][0],9);
            Assert.Equal(da3_drdsExpected[1, 2][1], da3_drds[1, 2][1],9);
            Assert.Equal(da3_drdsExpected[1, 2][2], da3_drds[1, 2][2],9);
                                                                     
            Assert.Equal(da3_drdsExpected[2, 0][0], da3_drds[2, 0][0],9);
            Assert.Equal(da3_drdsExpected[2, 0][1], da3_drds[2, 0][1],9);
            Assert.Equal(da3_drdsExpected[2, 0][2], da3_drds[2, 0][2],9);
            Assert.Equal(da3_drdsExpected[2, 1][0], da3_drds[2, 1][0],9);
            Assert.Equal(da3_drdsExpected[2, 1][1], da3_drds[2, 1][1],9);
            Assert.Equal(da3_drdsExpected[2, 1][2], da3_drds[2, 1][2],9);
            Assert.Equal(da3_drdsExpected[2, 2][0], da3_drds[2, 2][0],9);
            Assert.Equal(da3_drdsExpected[2, 2][1], da3_drds[2, 2][1],9);
            Assert.Equal(da3_drdsExpected[2, 2][2], da3_drds[2, 2][2],9);
            #endregion
        }

        private double[,] GetExpectedStiffnessMatrixForNonlinearityCheck()
        {

            return new double[48, 48]
            {
                {1716444723.664000,48908752.316760,213089244.821000,-448002390.873400,-12380544.874560,-165550004.526600,-966407266.488500,-37771219.743610,-172792633.144400,-137982232.076600,-5778441.871523,-21805477.189570,928376324.158000,105269438.332400,-17713640.040720,-337150051.807000,-25148377.157940,72168439.883600,-596524338.979500,-61681926.982530,56679557.301050,-85313796.885140,-8674070.198144,7114083.895979,172337495.184300,23414047.376740,-22770197.776290,-91546859.810120,-8908094.293958,24102679.618870,-128743897.474400,-14807967.725500,22792782.919870,-18063722.009260,-1966619.946085,2872857.970792,9579449.707148,1095658.160488,-1873379.793432,-7412799.726635,-635193.539235,1856495.253439,-8444136.840209,-830298.026321,1628909.845061,-1146499.742688,-105141.826984,200280.961398},
{48908752.316760,698077083.535200,-8773559.772368,-17473808.648510,-196480327.093900,-9972372.453455,-38502949.844640,-370619071.976200,1982542.102834,-5761599.796508,-53607523.352450,1409636.597223,111116296.810600,348545272.653300,-10509825.399080,-25395606.762060,-133695349.659400,11712140.651800,-62860666.194480,-227236504.833900,14362244.594330,-8853120.314119,-32553289.716040,2035629.077692,24697920.354930,62752169.663180,-9040106.241809,-8666162.773001,-35457522.942690,3458067.962965,-14865179.154960,-49169992.241780,3906194.901479,-1983490.020482,-6845439.654393,480023.550587,1177387.877662,2761647.001018,-1526718.694218,-608617.512885,-2863007.421718,194627.045608,-824080.024888,-3179296.243139,252150.438110,-105076.313441,-428847.717194,29325.638307},
{213089244.821000,-8773559.772368,778590775.819300,-182277586.311400,-8419550.885810,-203705375.676000,-183731945.866700,1426411.086863,-466722610.038500,-23137695.761760,1279715.191553,-71182311.100070,-10981139.738420,-12236940.718520,401165836.573100,84724263.879390,12020061.172190,-114264156.120300,61767900.784380,14379092.424770,-246007504.120400,7608848.039342,2035715.995364,-36653848.816980,-21698743.933650,-9146801.371653,62745553.990000,26066137.336600,3969932.966259,-38077846.089730,23619701.601640,3921430.582137,-53194772.537520,2956619.484618,476438.197608,-7448905.919630,-1832494.388443,-1494032.128482,1928864.691948,1952955.907904,278072.510935,-3177845.101799,1669369.367870,255084.865091,-3523076.878087,204564.777644,28929.884065,-472778.675372},
{-448002390.873400,-17473808.648510,-182277586.311400,2058416088.680000,175445192.549600,64541270.220650,305666772.394700,19849391.444740,74812385.461170,-157074448.982800,-13386285.904480,14012744.775030,-527759692.043500,-87270854.468620,-86652610.299570,74397079.443510,23022376.184660,-95575043.860650,-657054101.317400,-38213644.255790,119937376.935200,-231167753.909200,-15468582.195540,39976013.359510,-114466666.396900,-24295430.590650,4834703.166486,-44240888.936240,-7689520.748165,-650134.573982,-171629325.924800,-8535750.868367,31981791.774200,-52180038.384280,-2236744.603651,9100271.310895,-8647574.160560,-1748719.699067,1321396.980309,-8934893.059884,-1361451.890614,1665600.731800,-13795125.998250,-568852.500040,2421713.049512,-3527040.530731,-67313.805491,550107.280843},
{-12380544.874560,-196480327.093900,-8419550.885810,175445192.549600,962561264.733500,-86999713.981290,25766252.969530,246625689.922100,-52157534.394220,-11579262.083050,-46089630.520600,-213888.452694,-87518070.797380,-244261409.029500,-8176995.199610,17313010.932540,-67713552.576730,36654788.866670,-43210282.483890,-355850688.948500,84226467.669280,-16317685.075780,-108822134.886200,19600218.115980,-24024707.761000,-52133592.142860,1948573.883489,-7910983.179951,-27729161.172360,2965700.726961,-9369196.446418,-74213840.485560,7847437.074308,-2463885.747021,-21317561.397190,1838373.548065,-1724702.058698,-3894122.589645,236123.643760,-1329124.003283,-3857526.223617,274719.029317,-611471.869302,-5456343.458962,310223.180735,-84540.071311,-1367064.130021,65057.175057},
{-165550004.526600,-9972372.453455,-203705375.676000,64541270.220650,-86999713.981290,1501077371.632000,31067807.514370,-57603460.292580,580765674.003900,2387493.741746,-1812756.932541,-30580636.140950,-80468396.070610,-6297441.869654,-283078277.000700,-74795007.024200,41552918.493040,-386175360.005200,126826178.173500,84978829.487180,-709758597.508800,39662199.289020,19410150.735720,-182790891.060200,5452221.287132,2164183.775798,-51768848.090360,2208999.485632,3970399.979942,-66002588.860830,33362306.272980,7869761.547328,-119740323.720000,9168881.209442,1770663.726887,-30437104.185760,1310006.242899,241568.167802,-3784613.010429,1759377.416756,356933.880860,-5198477.140684,2504456.654063,310815.047397,-7127443.304868,562210.113196,59520.687559,-1694509.931406},
{-966407266.488500,-38502949.844640,-183731945.866700,305666772.394700,25766252.969530,31067807.514370,1853815603.758000,41790847.380210,447272722.393300,423191872.408800,4584936.359525,119324112.039800,-687075729.531100,-74825096.039540,-42777518.385660,-863505505.377900,-19239723.809830,-253367169.660800,173653881.422800,73644726.422470,-95979784.781700,100459802.604300,17070332.465980,5205433.627653,-140565064.268600,-20638370.581680,10452012.318880,-206598327.575900,-20706650.102960,-17301833.189400,15535644.495980,10348307.496680,-19719666.034890,19189004.114400,3984334.997136,-1922119.126476,-9199324.835989,-1365344.233040,1280402.548882,-16847152.262980,-2206898.048733,931686.737105,-2108093.604720,87219.420844,-579387.196203,793882.746982,208075.148048,-154752.938147},
{-37771219.743610,-370619071.976200,1426411.086863,19849391.444740,246625689.922100,-57603460.292580,41790847.380210,941105962.177800,-118259334.710400,5665559.878569,217786271.962700,-28239898.900620,-74581793.732270,-280403705.689800,8396802.220426,-17260753.835860,-483906308.030800,72850268.694740,75791413.818570,-106934057.528400,88555936.032010,17403925.082760,4732504.192706,18176737.185100,-20505425.366850,-56583237.159710,3210665.381662,-20508694.612950,-92398568.624340,5313374.637186,9768073.233813,-11149111.615680,4491305.918524,3746655.270669,3568029.605760,1007673.549151,-1353572.265297,-3660963.930574,281652.827787,-2204972.174367,-6887112.637398,335837.926301,-5293.526083,-1425981.742924,54440.461926,175859.147947,149661.074757,1587.981951},
{-172792633.144400,1982542.102834,-466722610.038500,74812385.461170,-52157534.394220,580765674.003900,447272722.393300,-118259334.710400,1754561671.189000,113301804.581300,-29127076.962950,417228095.068800,-42611358.498900,8837570.433150,-340035433.668900,-264878383.360800,72518957.143720,-994024660.318800,-120256040.990200,85343236.089310,-588479288.995800,-1178583.160808,17236325.283460,-78382183.025410,10356477.721510,3219883.841453,-59856928.234590,-19275632.549810,5138602.611173,-139867355.039200,-23168684.921550,3836633.230307,-64013620.015990,-2788644.988566,826237.099901,-6442666.041524,1265973.919209,277201.398790,-3782380.830926,826205.142838,324303.298184,-8090340.484320,-704359.966755,12819.608890,-2754466.593690,-181247.637494,-10366.073558,-103506.974432},
{-137982232.076600,-5761599.796508,-23137695.761760,-157074448.982800,-11579262.083050,2387493.741746,423191872.408800,5665559.878569,113301804.581300,148321860.979900,2908832.634881,33858762.746940,-94229693.603020,-8968819.037903,-4628934.702248,-274148097.870600,-7248423.451332,-57951221.792940,84312263.764780,22222964.009770,-45013617.754360,63623749.378170,6876419.856473,-5090849.723091,-19316616.331910,-2472272.456948,1390831.024252,-60685573.178110,-4998844.860563,-4150655.847740,15138522.259480,2676565.556025,-8921878.548189,13371578.407100,1274411.646772,-1847037.780204,-1236219.250727,-160856.811145,157779.956758,-4361920.358794,-483617.790898,139298.241057,342249.338646,-3973.077457,-378606.931776,732705.115681,52915.783313,-115471.449797},
{-5778441.871523,-53607523.352450,1279715.191553,-13386285.904480,-46089630.520600,-1812756.932541,4584936.359525,217786271.962700,-29127076.962950,2908832.634881,74899439.374400,-10375668.174510,-8884560.081059,-37877262.328840,1600980.143446,-6004205.464473,-136463402.550000,16740247.035510,23591883.452200,-6312516.784454,16502494.379820,7132786.474177,16843739.363880,2737896.816575,-2451999.141638,-7568979.115796,398444.601130,-4909765.540610,-25433028.323710,1176733.899038,2598346.143281,1751255.833549,673009.808715,1229186.285453,4065025.966709,113804.934091,-159498.946221,-479036.767810,31433.583892,-486970.883775,-1721647.034314,65676.768376,-28781.353729,-24604.120500,-600.371967,44537.837987,231898.397256,-4334.720175},
{-21805477.189570,1409636.597223,-71182311.100070,14012744.775030,-213888.452694,-30580636.140950,119324112.039800,-28239898.900620,417228095.068800,33858762.746940,-10375668.174510,141920656.324000,-4741886.817244,1613108.124541,-46988454.076140,-63131704.326970,16104575.578760,-250711400.434200,-53863266.944790,15262070.393430,-109350114.867900,-7294063.075466,2420732.452420,2301594.449685,1361532.467731,393569.257970,-8156817.109472,-4989885.722183,1046337.472254,-35400126.320400,-10256987.220950,453515.862688,-8846356.950818,-2169444.929799,59557.879949,2273738.285814,155540.486292,30688.536794,-501605.251642,99484.709801,57309.970162,-1957273.037467,-432211.347354,-13959.254852,-246793.527412,-127249.651297,-7687.343516,197804.688182},
{928376324.158000,111116296.810600,-10981139.738420,-527759692.043500,-87518070.797380,-80468396.070610,-687075729.531100,-74581793.732270,-42611358.498900,-94229693.603020,-8884560.081059,-4741886.817244,2255802955.655000,290980533.447600,-279112661.926200,-482833654.433700,-91273005.170520,168176838.682300,-1259892775.994000,-150825718.658700,205483411.791200,-186305723.793900,-19747977.984260,26544733.424880,1106962215.555000,106156226.025000,-187896725.086200,-307035655.201600,-23163068.335530,75992096.135210,-653751836.780500,-49804790.534870,115490641.615700,-94677299.117630,-6519272.560106,15274389.828920,138292911.025100,8745646.553818,-24415934.575990,-41696963.396020,-869080.915804,7875899.655807,-82398317.143170,-3376524.539920,13581987.372560,-11777065.354750,-434839.526662,1808104.206964},
{105269438.332400,348545272.653300,-12236940.718520,-87270854.468620,-244261409.029500,-6297441.869654,-74825096.039540,-280403705.689800,8837570.433150,-8968819.037903,-37877262.328840,1613108.124541,290980533.447600,923598815.181000,-46997516.816290,-91011840.391150,-188334738.632900,29280679.288110,-151251636.090300,-480166923.511400,33750455.875600,-19851468.193150,-69997681.204120,4060215.995050,108612677.868300,421057646.227800,-35781815.261080,-21340121.619870,-114555901.612000,9532750.887928,-49140700.057080,-243074145.705400,15544051.280520,-6459364.528269,-34919211.162460,1946573.335264,9246026.268208,49643994.006680,-5197386.104423,-393177.734187,-15083872.736610,329868.807989,-3181234.926675,-29907001.887510,1426556.874272,-414362.829803,-4263874.568292,189269.867533},
{-17713640.040720,-10509825.399080,401165836.573100,-86652610.299570,-8176995.199610,-283078277.000700,-42777518.385660,8396802.220426,-340035433.668900,-4628934.702248,1600980.143446,-46988454.076140,-279112661.926200,-46997516.816290,964619026.642400,172904836.102500,30420663.517400,-149808591.883100,209451105.334600,34044044.949360,-488885510.800200,27054087.891260,4087697.579862,-72658265.544040,-187132387.830600,-37452550.815160,445945359.345500,77570783.666320,10615021.556250,-127724320.945400,116695975.915600,15634840.434630,-264042041.193300,15436764.231420,1945393.886226,-37677855.043890,-24340787.197640,-5617509.992776,53198619.155040,7840790.957314,417624.752484,-16685275.926670,13591178.687010,1405153.288935,-32721912.956230,1813017.596513,186175.893907,-4622902.677496},
{-337150051.807000,-25395606.762060,84724263.879390,74397079.443510,17313010.932540,-74795007.024200,-863505505.377900,-17260753.835860,-264878383.360800,-274148097.870600,-6004205.464473,-63131704.326970,-482833654.433700,-91011840.391150,172904836.102500,2295498083.436000,141034053.234300,-24038164.327610,-149618273.298400,-20983952.535200,80146490.005240,-333459183.475300,-19314670.674180,28014915.475460,-355550596.782800,-44281476.700190,80317533.906050,855881347.457200,54522555.338700,-79047627.451290,-215573530.721600,10891706.177610,33034472.325870,-178957004.896700,-2365426.982283,18829934.691820,-56445029.279380,-4684705.128863,11346809.639560,80246092.824800,4146037.233691,-8786705.169506,-36292709.336580,2990000.331596,3246147.733998,-22488965.881720,405275.225853,2112187.900466},
{-25148377.157940,-133695349.659400,12020061.172190,23022376.184660,-67713552.576730,41552918.493040,-19239723.809830,-483906308.030800,72518957.143720,-7248423.451332,-136463402.550000,16104575.578760,-91273005.170520,-188334738.632900,30420663.517400,141034053.234300,992880368.355400,-87197180.669260,-24337411.851190,76800368.143470,-70362567.979640,-20360706.999020,-99797439.077460,-9666882.185976,-44364308.225660,-142872989.031200,11284345.466290,56035741.081160,326708987.721800,-16208892.867030,10935649.707400,-65133397.091650,-3043100.442725,-2566584.078075,-63250722.646020,1419489.164051,-4735261.680471,-22161585.329160,717892.522966,4590174.757577,28821233.388950,-855280.105701,3237840.999354,-13539449.450680,860974.289716,417966.459625,-8342023.533753,434026.902207},
{72168439.883600,11712140.651800,-114264156.120300,-95575043.860650,36654788.866670,-386175360.005200,-253367169.660800,72850268.694740,-994024660.318800,-57951221.792940,16740247.035510,-250711400.434200,168176838.682300,29280679.288110,-149808591.883100,-24038164.327610,-87197180.669260,1425457013.489000,94761358.255490,-69114398.411770,406872081.104900,32040397.195360,-9348694.805370,-54484045.482860,78937559.884910,11434233.185340,-154404684.189300,-79175922.144170,-14169028.610270,380901189.255900,36754736.868030,-2296361.444790,-34883674.244490,19891246.455890,1452789.810961,-60571691.674760,11112988.143770,928462.558117,-24791480.280960,-9194450.615818,-365295.915551,31983724.488860,3289942.358756,1000334.766914,-12677003.779740,2168464.673920,437014.998841,-8417259.924734},
{-596524338.979500,-62860666.194480,61767900.784380,-657054101.317400,-43210282.483890,126826178.173500,173653881.422800,75791413.818570,-120256040.990200,84312263.764780,23591883.452200,-53863266.944790,-1259892775.994000,-151251636.090300,209451105.334600,-149618273.298400,-24337411.851190,94761358.255490,1984498039.934000,180883677.605200,-258390855.457200,502050639.038600,42596136.156400,-72290402.800500,-683378774.162900,-61811954.738540,118714397.142200,-335822299.720300,-32510327.631740,43124960.526630,784758676.369100,50527546.988950,-129543722.113900,218262059.992200,13527573.619140,-33398067.893630,-91275363.574810,-5894894.589281,15853979.885430,-71759779.944020,-7011945.118445,11670178.684150,74246258.927720,1212677.575692,-11061627.476560,23543887.542480,758209.481661,-3366075.109544},
{-61681926.982530,-227236504.833900,14379092.424770,-38213644.255790,-355850688.948500,84978829.487180,73644726.422470,-106934057.528400,85343236.089310,22222964.009770,-6312516.784454,15262070.393430,-150825718.658700,-480166923.511400,34044044.949360,-20983952.535200,76800368.143470,-69114398.411770,180883677.605200,900897821.008100,-135253831.441200,42052199.025690,217787955.886000,-28692497.572880,-61892462.639110,-257952899.867900,17553180.556860,-33345241.934490,-115476782.571200,2273337.517970,47730577.582910,300084827.813300,-20497364.659800,12677969.791180,81999968.822480,-4761795.372160,-5940084.786358,-33845767.656730,1985396.401441,-7408389.211598,-27533798.939010,2480045.212946,504483.907336,25515666.199280,115929.296742,574822.659147,8223332.768709,-95274.872196},
{56679557.301050,14362244.594330,-246007504.120400,119937376.935200,84226467.669280,-709758597.508800,-95979784.781700,88555936.032010,-588479288.995800,-45013617.754360,16502494.379820,-109350114.867900,205483411.791200,33750455.875600,-488885510.800200,80146490.005240,-70362567.979640,406872081.104900,-258390855.457200,-135253831.441200,1350480443.264000,-70216982.607800,-28514586.683780,309895665.540800,117557319.888800,17279666.378430,-277700313.605100,39173042.723950,1511327.920324,-85098727.125990,-129864911.366100,-21284157.733620,368621485.664300,-32936516.810930,-4928543.727287,96830594.139590,15762041.340040,1950163.642895,-37080396.996490,11621375.659890,2369160.156373,-28500488.487890,-10713297.088150,-33879.871580,29056897.205830,-3244649.778939,-130349.211964,9103775.589286},
{-85313796.885140,-8853120.314119,7608848.039342,-231167753.909200,-16317685.075780,39662199.289020,100459802.604300,17403925.082760,-1178583.160808,63623749.378170,7132786.474177,-7294063.075466,-186305723.793900,-19851468.193150,27054087.891260,-333459183.475300,-20360706.999020,32040397.195360,502050639.038600,42052199.025690,-70216982.607800,204552421.622900,13173436.037320,-24040429.168260,-98201965.790940,-7885302.321522,15670612.629640,-211899944.054800,-13714165.516800,21350666.895900,200415702.646200,7660164.687331,-32272240.855610,91849300.792700,3145493.857903,-11428640.015560,-12825834.518520,-738935.892960,2080856.860091,-32163041.895080,-2245697.703312,4317287.576366,18324009.874710,-617834.985464,-2252648.949894,10061618.365270,16911.836948,-1101368.543586},
{-8674070.198144,-32553289.716040,2035715.995364,-15468582.195540,-108822134.886200,19410150.735720,17070332.465980,4732504.192706,17236325.283460,6876419.856473,16843739.363880,2420732.452420,-19747977.984260,-69997681.204120,4087697.579862,-19314670.674180,-99797439.077460,-9348694.805370,42596136.156400,217787955.886000,-28514586.683780,13173436.037320,81861261.464840,-7064535.610458,-7897888.723472,-36581429.758000,2160696.743556,-14128077.088300,-76170028.438510,2497722.886527,6733419.020780,75485667.307590,-4378822.566152,2883014.787589,33920651.319170,-1499180.871443,-747048.522210,-4716355.790654,252896.689432,-2413551.676275,-11975320.826740,769679.470101,-890310.322875,6376335.018678,22727.636983,-50580.939294,3605565.144919,-88524.936220},
{7114083.895979,2035629.077692,-36653848.816980,39976013.359510,19600218.115980,-182790891.060200,5205433.627653,18176737.185100,-78382183.025410,-5090849.723091,2737896.816575,2301594.449685,26544733.424880,4060215.995050,-72658265.544040,28014915.475460,-9666882.185976,-54484045.482860,-72290402.800500,-28692497.572880,309895665.540800,-24040429.168260,-7064535.610458,104949330.902800,15513697.371310,2124954.383424,-39255354.998790,20212241.326040,2184167.497151,-73836772.447640,-32908466.015800,-4731806.652328,89945541.969650,-11451814.772030,-1570328.036707,37839737.114260,2070081.903508,246269.299765,-5128440.835189,4335807.976772,710595.108998,-12552252.947990,-2139152.862993,-47378.883395,6990305.644296,-1065893.018436,-103254.537990,3819879.537561},
{172337495.184300,24697920.354930,-21698743.933650,-114466666.396900,-24024707.761000,5452221.287132,-140565064.268600,-20505425.366850,10356477.721510,-19316616.331910,-2451999.141638,1361532.467731,1106962215.555000,108612677.868300,-187132387.830600,-355550596.782800,-44364308.225660,78937559.884910,-683378774.162900,-61892462.639110,117557319.888800,-98201965.790940,-7897888.723472,15513697.371310,958827176.371600,58394209.280500,-158163032.297100,-241093484.186200,-9056410.077552,44183521.763510,-543238136.919300,-23979250.007310,89054304.378780,-78400991.739470,-3113728.996488,12055035.512450,165338235.167500,5806238.890732,-25384881.677910,-30759763.945130,1000123.221615,3712233.524821,-85988546.244530,-1080759.138167,12464669.343750,-12504515.510220,-144229.538810,1730472.594505},
{23414047.376740,62752169.663180,-9146801.371653,-24295430.590650,-52133592.142860,2164183.775798,-20638370.581680,-56583237.159710,3219883.841453,-2472272.456948,-7568979.115796,393569.257970,106156226.025000,421057646.227800,-37452550.815160,-44281476.700190,-142872989.031200,11434233.185340,-61811954.738540,-257952899.867900,17279666.378430,-7885302.321522,-36581429.758000,2124954.383424,58394209.280500,365562770.025600,-6102481.781487,-7195224.114423,-86775207.737170,3433813.110982,-23063687.722410,-197841041.897900,9025628.891111,-3001896.861650,-28452858.388780,1218785.816592,6029100.685142,62631837.972110,1513975.341143,1570927.478594,-10161443.228470,-23053.221523,-806918.484188,-30618235.827940,792240.240686,-111976.273746,-4462509.733002,123952.966899},
{-22770197.776290,-9040106.241809,62745553.990000,4834703.166486,1948573.883489,-51768848.090360,10452012.318880,3210665.381662,-59856928.234590,1390831.024252,398444.601130,-8156817.109472,-187896725.086200,-35781815.261080,445945359.345500,80317533.906050,11284345.466290,-154404684.189300,118714397.142200,17553180.556860,-277700313.605100,15670612.629640,2160696.743556,-39255354.998790,-158163032.297100,-6102481.781487,410668162.673500,44099018.362220,2307165.115358,-101926331.084200,89125231.721910,9057429.228083,-216904148.912400,12078127.313510,1233481.135963,-30860795.410670,-25366908.571310,1223991.447162,71426397.425220,3463309.763592,-348564.854041,-11895014.586230,12333188.151570,770266.055049,-33256204.437930,1717898.230540,124728.523822,-4800032.775178},
{-91546859.810120,-8666162.773001,26066137.336600,-44240888.936240,-7910983.179951,2208999.485632,-206598327.575900,-20508694.612950,-19275632.549810,-60685573.178110,-4909765.540610,-4989885.722183,-307035655.201600,-21340121.619870,77570783.666320,855881347.457200,56035741.081160,-79175922.144170,-335822299.720300,-33345241.934490,39173042.723950,-211899944.054800,-14128077.088300,20212241.326040,-241093484.186200,-7195224.114423,44099018.362220,850431053.387100,51095660.690000,-118409377.272500,-144585391.424500,3066555.088328,13307239.273150,-149514008.957400,-3684390.588349,16606130.855660,-40237405.501060,-39499.498269,5994347.206194,155558573.235000,7865835.831414,-23462815.887160,-7254653.119259,3423217.329642,-1936476.354374,-21356482.413900,241150.929666,2012169.694460},
{-8908094.293958,-35457522.942690,3969932.966259,-7689520.748165,-27729161.172360,3970399.979942,-20706650.102960,-92398568.624340,5138602.611173,-4998844.860563,-25433028.323710,1046337.472254,-23163068.335530,-114555901.612000,10615021.556250,54522555.338700,326708987.721800,-14169028.610270,-32510327.631740,-115476782.571200,1511327.920324,-13714165.516800,-76170028.438510,2184167.497151,-9056410.077552,-86775207.737170,2307165.115358,51095660.690000,316198834.994500,-15355809.365810,5745325.119478,-48796536.225620,286813.147454,-2930120.800173,-53931721.338520,2072993.353522,-374513.142289,-13989862.087180,-195108.735465,8102719.109361,57386393.793900,-2972326.615777,4165235.204269,-1913053.891105,-636757.348295,420220.047917,-7666841.545754,226269.055921},
{24102679.618870,3458067.962965,-38077846.089730,-650134.573982,2965700.726961,-66002588.860830,-17301833.189400,5313374.637186,-139867355.039200,-4150655.847740,1176733.899038,-35400126.320400,75992096.135210,9532750.887928,-127724320.945400,-79047627.451290,-16208892.867030,380901189.255900,43124960.526630,2273337.517970,-85098727.125990,21350666.895900,2497722.886527,-73836772.447640,44183521.763510,3433813.110982,-101926331.084200,-118409377.272500,-15355809.365810,347268873.236200,12666679.595690,1239606.533492,-42053258.644120,16472071.230590,2305271.342795,-55131252.997010,6050315.957028,231817.543181,-16667972.267930,-23753570.454950,-2737018.004887,62404560.943580,-2517013.950713,-397283.813299,-900920.612860,1887221.017161,270807.002004,-7887151.000433},
{-128743897.474400,-14865179.154960,23619701.601640,-171629325.924800,-9369196.446418,33362306.272980,15535644.495980,9768073.233813,-23168684.921550,15138522.259480,2598346.143281,-10256987.220950,-653751836.780500,-49140700.057080,116695975.915600,-215573530.721600,10935649.707400,36754736.868030,784758676.369100,47730577.582910,-129864911.366100,200415702.646200,6733419.020780,-32908466.015800,-543238136.919300,-23063687.722410,89125231.721910,-144585391.424500,5745325.119478,12666679.595690,664105972.691900,15364736.203490,-96069467.638380,164955724.879400,739856.789389,-21216052.294520,-91481653.549180,-1941796.844622,13931211.143940,-28157661.597300,-360491.429028,3351226.095744,105807453.776400,-488219.779281,-13120068.657450,26443737.272920,-386712.366741,-2902431.100749},
{-14807967.725500,-49169992.241780,3921430.582137,-8535750.868367,-74213840.485560,7869761.547328,10348307.496680,-11149111.615680,3836633.230307,2676565.556025,1751255.833549,453515.862688,-49804790.534870,-243074145.705400,15634840.434630,10891706.177610,-65133397.091650,-2296361.444790,50527546.988950,300084827.813300,-21284157.733620,7660164.687331,75485667.307590,-4731806.652328,-23979250.007310,-197841041.897900,9057429.228083,3066555.088328,-48796536.225620,1239606.533492,15364736.203490,244723964.806200,-10233559.294150,1117956.374293,60861386.361860,-2546903.460088,-2141678.190693,-32782251.807510,1009993.255632,-1165364.304662,-9591527.005800,150697.948320,-832119.043703,39024031.081900,-1628692.259342,-386617.897597,9820710.872437,-452427.778293},
{22792782.919870,3906194.901479,-53194772.537520,31981791.774200,7847437.074308,-119740323.720000,-19719666.034890,4491305.918524,-64013620.015990,-8921878.548189,673009.808715,-8846356.950818,115490641.615700,15544051.280520,-264042041.193300,33034472.325870,-3043100.442725,-34883674.244490,-129543722.113900,-20497364.659800,368621485.664300,-32272240.855610,-4378822.566152,89945541.969650,89054304.378780,9025628.891111,-216904148.912400,13307239.273150,286813.147454,-42053258.644120,-96069467.638380,-10233559.294150,272865972.790900,-21328616.286180,-2417035.606719,66411451.886790,13997217.854710,997882.558655,-35878973.816800,3940864.335096,-66090.474248,-10126578.309400,-12843047.041200,-1692102.441514,41519763.787920,-2900675.959081,-444248.095458,10319532.245390},
{-18063722.009260,-1983490.020482,2956619.484618,-52180038.384280,-2463885.747021,9168881.209442,19189004.114400,3746655.270669,-2788644.988566,13371578.407100,1229186.285453,-2169444.929799,-94677299.117630,-6459364.528269,15436764.231420,-178957004.896700,-2566584.078075,19891246.455890,218262059.992200,12677969.791180,-32936516.810930,91849300.792700,2883014.787589,-11451814.772030,-78400991.739470,-3001896.861650,12078127.313510,-149514008.957400,-2930120.800173,16472071.230590,164955724.879400,1117956.374293,-21328616.286180,70217535.686630,-69293.599260,-7480571.367998,-13140195.697100,-253448.954165,1906153.207341,-26976672.273750,-757686.027762,3390594.923517,23451787.768520,-903099.737582,-2213004.818660,10612941.434630,-265912.154749,-931844.082159},
{-1966619.946085,-6845439.654393,476438.197608,-2236744.603651,-21317561.397190,1770663.726887,3984334.997136,3568029.605760,826237.099901,1274411.646772,4065025.966709,59557.879949,-6519272.560106,-34919211.162460,1945393.886226,-2365426.982283,-63250722.646020,1452789.810961,13527573.619140,81999968.822480,-4928543.727287,3145493.857903,33920651.319170,-1570328.036707,-3113728.996488,-28452858.388780,1233481.135963,-3684390.588349,-53931721.338520,2305271.342795,739856.789389,60861386.361860,-2417035.606719,-69293.599260,25894825.235990,-1093840.692601,-279783.159578,-4711479.226500,151825.712557,-1012515.201936,-9659435.854940,385365.041592,-1125182.164677,8808294.114570,-397199.863459,-298713.107926,3970248.242256,-200075.907665},
{2872857.970792,480023.550587,-7448905.919630,9100271.310895,1838373.548065,-30437104.185760,-1922119.126476,1007673.549151,-6442666.041524,-1847037.780204,113804.934091,2273738.285814,15274389.828920,1946573.335264,-37677855.043890,18829934.691820,1419489.164051,-60571691.674760,-33398067.893630,-4761795.372160,96830594.139590,-11428640.015560,-1499180.871443,37839737.114260,12055035.512450,1218785.816592,-30860795.410670,16606130.855660,2072993.353522,-55131252.997010,-21216052.294520,-2546903.460088,66411451.886790,-7480571.367998,-1093840.692601,27377654.331590,1913767.741007,145965.969146,-5101121.744646,3575233.261711,316611.927382,-10269199.651790,-2030838.003498,-451873.104003,9118256.301172,-904294.691361,-206701.647556,4089160.610467},
{9579449.707148,1177387.877662,-1832494.388443,-8647574.160560,-1724702.058698,1310006.242899,-9199324.835989,-1353572.265297,1265973.919209,-1236219.250727,-159498.946221,155540.486292,138292911.025100,9246026.268208,-24340787.197640,-56445029.279380,-4735261.680471,11112988.143770,-91275363.574810,-5940084.786358,15762041.340040,-12825834.518520,-747048.522210,2070081.903508,165338235.167500,6029100.685142,-25366908.571310,-40237405.501060,-374513.142289,6050315.957028,-91481653.549180,-2141678.190693,13997217.854710,-13140195.697100,-279783.159578,1913767.741007,33770361.853070,569833.855842,-4770783.830951,-4073691.070847,368347.111151,228305.511133,-16055393.973700,59412.882578,2139069.283674,-2363272.340931,6034.071231,305665.605089},
{1095658.160488,2761647.001018,-1494032.128482,-1748719.699067,-3894122.589645,241568.167802,-1365344.233040,-3660963.930574,277201.398790,-160856.811145,-479036.767810,30688.536794,8745646.553818,49643994.006680,-5617509.992776,-4684705.128863,-22161585.329160,928462.558117,-5894894.589281,-33845767.656730,1950163.642895,-738935.892960,-4716355.790654,246269.299765,5806238.890732,62631837.972110,1223991.447162,-39499.498269,-13989862.087180,231817.543181,-1941796.844622,-32782251.807510,997882.558655,-253448.954165,-4711479.226500,145965.969146,569833.855842,12852575.222340,776505.513365,476911.213052,-1175468.580372,-15657.023217,120308.596032,-5637898.016186,63560.350112,13604.381448,-835262.419825,13122.158693},
{-1873379.793432,-1526718.694218,1928864.691948,1321396.980309,236123.643760,-3784613.010429,1280402.548882,281652.827787,-3782380.830926,157779.956758,31433.583892,-501605.251642,-24415934.575990,-5197386.104423,53198619.155040,11346809.639560,717892.522966,-24791480.280960,15853979.885430,1985396.401441,-37080396.996490,2080856.860091,252896.689432,-5128440.835189,-25384881.677910,1513975.341143,71426397.425220,5994347.206194,-195108.735465,-16667972.267930,13931211.143940,1009993.255632,-35878973.816800,1906153.207341,151825.712557,-5101121.744646,-4770783.830951,776505.513365,14680878.680420,170530.053201,-116824.794092,-1524969.551531,2100297.029551,64123.646023,-6098645.681070,301215.367023,14219.190200,-894159.685003},
{-7412799.726635,-608617.512885,1952955.907904,-8934893.059884,-1329124.003283,1759377.416756,-16847152.262980,-2204972.174367,826205.142838,-4361920.358794,-486970.883775,99484.709801,-41696963.396020,-393177.734187,7840790.957314,80246092.824800,4590174.757577,-9194450.615818,-71759779.944020,-7408389.211598,11621375.659890,-32163041.895080,-2413551.676275,4335807.976772,-30759763.945130,1570927.478594,3463309.763592,155558573.235000,8102719.109361,-23753570.454950,-28157661.597300,-1165364.304662,3940864.335096,-26976672.273750,-1012515.201936,3575233.261711,-4073691.070847,476911.213052,170530.053201,39190316.024070,1728623.151699,-6268098.018340,2068811.686421,552463.058460,-789702.429206,-3919454.239810,863.934226,419886.333439},
{-635193.539235,-2863007.421718,278072.510935,-1361451.890614,-3857526.223617,356933.880860,-2206898.048733,-6887112.637398,324303.298184,-483617.790898,-1721647.034314,57309.970162,-869080.915804,-15083872.736610,417624.752484,4146037.233691,28821233.388950,-365295.915551,-7011945.118445,-27533798.939010,2369160.156373,-2245697.703312,-11975320.826740,710595.108998,1000123.221615,-10161443.228470,-348564.854041,7865835.831414,57386393.793900,-2737018.004887,-360491.429028,-9591527.005800,-66090.474248,-757686.027762,-9659435.854940,316611.927382,368347.111151,-1175468.580372,-116824.794092,1728623.151699,14526858.070410,-783794.891099,763376.801544,1118604.954427,-399706.241690,59719.112717,-1342929.718709,-13316.429769},
{1856495.253439,194627.045608,-3177845.101799,1665600.731800,274719.029317,-5198477.140684,931686.737105,335837.926301,-8090340.484320,139298.241057,65676.768376,-1957273.037467,7875899.655807,329868.807989,-16685275.926670,-8786705.169506,-855280.105701,31983724.488860,11670178.684150,2480045.212946,-28500488.487890,4317287.576366,769679.470101,-12552252.947990,3712233.524821,-23053.221523,-11895014.586230,-23462815.887160,-2972326.615777,62404560.943580,3351226.095744,150697.948320,-10126578.309400,3390594.923517,385365.041592,-10269199.651790,228305.511133,-15657.023217,-1524969.551531,-6268098.018340,-783794.891099,15712298.331360,-987601.776848,-338140.791454,1291189.233593,366413.916915,1735.398220,-1414057.771627},
{-8444136.840209,-824080.024888,1669369.367870,-13795125.998250,-611471.869302,2504456.654063,-2108093.604720,-5293.526083,-704359.966755,342249.338646,-28781.353729,-432211.347354,-82398317.143170,-3181234.926675,13591178.687010,-36292709.336580,3237840.999354,3289942.358756,74246258.927720,504483.907336,-10713297.088150,18324009.874710,-890310.322875,-2139152.862993,-85988546.244530,-806918.484188,12333188.151570,-7254653.119259,4165235.204269,-2517013.950713,105807453.776400,-832119.043703,-12843047.041200,23451787.768520,-1125182.164677,-2030838.003498,-16055393.973700,120308.596032,2100297.029551,2068811.686421,763376.801544,-987601.776848,23155577.754290,-278557.166285,-2699294.812035,4940827.133680,-207296.626129,-421615.399273},
{-830298.026321,-3179296.243139,255084.865091,-568852.500040,-5456343.458962,310815.047397,87219.420844,-1425981.742924,12819.608890,-3973.077457,-24604.120500,-13959.254852,-3376524.539920,-29907001.887510,1405153.288935,2990000.331596,-13539449.450680,1000334.766914,1212677.575692,25515666.199280,-33879.871580,-617834.985464,6376335.018678,-47378.883395,-1080759.138167,-30618235.827940,770266.055049,3423217.329642,-1913053.891105,-397283.813299,-488219.779281,39024031.081900,-1692102.441514,-903099.737582,8808294.114570,-451873.104003,59412.882578,-5637898.016186,64123.646023,552463.058460,1118604.954427,-338140.791454,-278557.166285,8911700.963311,-674892.418393,-176871.648296,1947232.306776,-169086.699809},
{1628909.845061,252150.438110,-3523076.878087,2421713.049512,310223.180735,-7127443.304868,-579387.196203,54440.461926,-2754466.593690,-378606.931776,-600.371967,-246793.527412,13581987.372560,1426556.874272,-32721912.956230,3246147.733998,860974.289716,-12677003.779740,-11061627.476560,115929.296742,29056897.205830,-2252648.949894,22727.636983,6990305.644296,12464669.343750,792240.240686,-33256204.437930,-1936476.354374,-636757.348295,-900920.612860,-13120068.657450,-1628692.259342,41519763.787920,-2213004.818660,-397199.863459,9118256.301172,2139069.283674,63560.350112,-6098645.681070,-789702.429206,-399706.241690,1291189.233593,-2699294.812035,-674892.418393,9330587.098848,-451679.002395,-160954.266137,1999468.500232},
{-1146499.742688,-105076.313441,204564.777644,-3527040.530731,-84540.071311,562210.113196,793882.746982,175859.147947,-181247.637494,732705.115681,44537.837987,-127249.651297,-11777065.354750,-414362.829803,1813017.596513,-22488965.881720,417966.459625,2168464.673920,23543887.542480,574822.659147,-3244649.778939,10061618.365270,-50580.939294,-1065893.018436,-12504515.510220,-111976.273746,1717898.230540,-21356482.413900,420220.047917,1887221.017161,26443737.272920,-386617.897597,-2900675.959081,10612941.434630,-298713.107926,-904294.691361,-2363272.340931,13604.381448,301215.367023,-3919454.239810,59719.112717,366413.916915,4940827.133680,-176871.648296,-451679.002395,1953696.403104,-77990.565378,-145315.953910},
{-105141.826984,-428847.717194,28929.884065,-67313.805491,-1367064.130021,59520.687559,208075.148048,149661.074757,-10366.073558,52915.783313,231898.397256,-7687.343516,-434839.526662,-4263874.568292,186175.893907,405275.225853,-8342023.533753,437014.998841,758209.481661,8223332.768709,-130349.211964,16911.836948,3605565.144919,-103254.537990,-144229.538810,-4462509.733002,124728.523822,241150.929666,-7666841.545754,270807.002004,-386712.366741,9820710.872437,-444248.095458,-265912.154749,3970248.242256,-206701.647556,6034.071231,-835262.419825,14219.190200,863.934226,-1342929.718709,1735.398220,-207296.626129,1947232.306776,-160954.266137,-77990.565378,760704.559440,-59570.402439},
{200280.961398,29325.638307,-472778.675372,550107.280843,65057.175057,-1694509.931406,-154752.938147,1587.981951,-103506.974432,-115471.449797,-4334.720175,197804.688182,1808104.206964,189269.867533,-4622902.677496,2112187.900466,434026.902207,-8417259.924734,-3366075.109544,-95274.872196,9103775.589286,-1101368.543586,-88524.936220,3819879.537561,1730472.594505,123952.966899,-4800032.775178,2012169.694460,226269.055921,-7887151.000433,-2902431.100749,-452427.778293,10319532.245390,-931844.082159,-200075.907665,4089160.610467,305665.605089,13122.158693,-894159.685003,419886.333439,-13316.429769,-1414057.771627,-421615.399273,-169086.699809,1999468.500232,-145315.953910,-59570.402439,776738.244560},

            };
        }
    }
}
