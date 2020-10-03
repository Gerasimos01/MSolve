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
    public class NurbsNonLinearShellsDevelop
    {
        private const double Tolerance = 1e-9;

        #region Fields
        private List<ControlPoint> ElementControlPoints()
        {
            return new List<ControlPoint>
            {
                new ControlPoint {ID = 0 , X = 10          , Y =  0,           Z = 0, WeightFactor = 1.0},
                new ControlPoint {ID = 1 , X = 10.00000572 , Y =  0,           Z = 0.36249144, WeightFactor = 1.0},
                new ControlPoint {ID = 2 , X = 9.96053846  , Y =  0,           Z = 1.087473387, WeightFactor = 1.0},
                new ControlPoint {ID = 3 , X = 9.783676177 , Y =  0,           Z = 2.162099139, WeightFactor = 1.0},
                new ControlPoint {ID = 16, X = 10.00000917 , Y =  0.4027681,   Z = 0, WeightFactor = 1.0},
                new ControlPoint {ID = 17, X = 10.0000149  , Y =  0.40276833,  Z = 0.36249144, WeightFactor = 1.0},
                new ControlPoint {ID = 18, X = 9.960547598 , Y =  0.401178715, Z = 1.087473387, WeightFactor = 1.0},
                new ControlPoint {ID = 19, X = 9.783685153 , Y =  0.394055266, Z = 2.162099139, WeightFactor = 1.0},
                new ControlPoint {ID = 32, X = 9.951271191 , Y =  1.208302888, Z = 0, WeightFactor = 1.0},
                new ControlPoint {ID = 33, X = 9.951276884 , Y =  1.208303579, Z = 0.36249144, WeightFactor = 1.0},
                new ControlPoint {ID = 34, X = 9.912001942 , Y =  1.203534738, Z = 1.087473387, WeightFactor = 1.0},
                new ControlPoint {ID = 35, X = 9.736001488 , Y =  1.182164418, Z = 2.162099139, WeightFactor = 1.0},
                new ControlPoint {ID = 48, X = 9.733071058 , Y =  2.398986433, Z = 0, WeightFactor = 1.0},
                new ControlPoint {ID = 49, X = 9.733076627 , Y =  2.398987806, Z = 0.36249144, WeightFactor = 1.0},
                new ControlPoint {ID = 50, X = 9.694662861 , Y =  2.389519663, Z = 1.087473387, WeightFactor = 1.0},
                new ControlPoint {ID = 51, X = 9.522521544 , Y =  2.347090642, Z = 2.162099139, WeightFactor = 1.0},
            };
        }

        private List<Knot> ElementKnots()
        {
            return new List<Knot>()
            {
                new Knot(){ID=0 ,Ksi=0.000,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=1 ,Ksi=0.000,Heta=1.0,Zeta =0.0 },
                new Knot(){ID=14,Ksi=1.000,Heta=0.0,Zeta =0.0 },
                new Knot(){ID=15,Ksi=1.000,Heta=1.0,Zeta =0.0 }
            };
        }

        private Vector KnotValueVectorKsi()
        {
            return Vector.CreateFromArray(new double[]
            {
                0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 13,             });
        }

        private Vector KnotValueVectorHeta()
        {
            return Vector.CreateFromArray(new double[]
            {
                0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 13            });
        }

        private ShellElasticMaterial2D Material => new ShellElasticMaterial2D()
        {
            YoungModulus = 68250000,
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
                patch.NumberOfControlPointsHeta = 16;
                patch.KnotValueVectorKsi = KnotValueVectorKsi();
                patch.KnotValueVectorHeta = KnotValueVectorHeta();
                var element =
                    new NurbsKirchhoffLoveShellElementNLDevelop(Material, ElementKnots(), ElementControlPoints(), patch, 0.04);
                element._solution= localSolution;
                return element;
            }
        }

        private double[] localSolution => new double[]
        {
            -0.073810726680524588, 0, -0.0048017449915675444, -0.068280013293049779, 0, -0.0048066252245659684, -0.057841468788129791,
            0, -0.0041856398772149731, -0.0451468058143807, 0, -0.0020263926209961023, -0.073810725359380713, 5.409657482823016E-06, -0.0048017454378473283,
            -0.068280011925985373, -8.3550104875203951E-06, -0.0048066261867408056, -0.057841467343111866, -1.1739995839693453E-06, -0.0041856405549517486, -0.045146804399930685, 7.2981249203213578E-07,
            -0.0020263928646768238, -0.071670829076570269, 6.971822182063464E-05, -0.0039875076812788731, -0.066595662362529576, 6.7553000990231045E-06, -0.0039920290411740063, -0.056722841520517792, 1.4477121129733133E-06,
            -0.00341855470621827, -0.044372004315413965, 5.1246110226694326E-06, -0.0013318231144180889, -0.0654240849246012, 0.0011587487507479593, -0.00062590319250741856, -0.060956634696546413, 0.00095412179398749125, -0.00062916764286217,
            -0.05230357168163733, 0.00073111204034573426, -0.00014611305505038129, -0.041166807382690719, 0.00053554065379350516, 0.0016763898708845439
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
                Amount = -400,
                Node = model.ControlPoints.ToList()[0],
                DOF = StructuralDof.TranslationX
            });

            model.Loads.Add(new Load()
            {
                Amount = 400,
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
            var newtonRaphsonBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 500);
            var childAnalyzer = newtonRaphsonBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 500,
                model.ControlPoints.ToList()[0], StructuralDof.TranslationX, "PinchedHemisphereNegativeLoadNode.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 500,
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
            modelReader.CreateShellModelFromFile(GeometricalFormulation.NonLinear);

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

        [Fact]
        public void IsogeometricSquareShell10x10Straight30Degrees()
        {
            Model model = new Model();
            var filename = "SquareShell10x10Straight30Degrees";
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

            double[] loadValues = new double[] { 0, 0, load_factor, 0, 0, 0 };
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
                return new double[] { 0,  0, load_factor };
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
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "squarePlateZDevelop.txt");
            var loggerB = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationX, "squareplateXDevelop.txt");
            var loggerC = new TotalLoadsDisplacementsPerIncrementLog(model.PatchesDictionary[0], 1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationY, "squareplateYDevelop.txt");
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
        public void InternalForcesTest()
        {
            var controlPoints = ElementControlPoints().ToArray();
            var shellElement = Element;

            var nurbs = new Nurbs2D(shellElement, controlPoints);
            var gaussPoints = shellElement.materialsAtThicknessGP.Keys.ToArray();
            shellElement.CalculateInitialConfigurationData(controlPoints, nurbs, gaussPoints);
            shellElement.CalculateStresses(shellElement, localSolution, new double[27]);
            var f = shellElement.CalculateForces(shellElement,localSolution, new double[27]);

            var expectedForces = new double[]
            {
                2.61204053011917,
                -39.584849368062144,
                -103.03746965256518,
                4.9982594203934383    ,
                -40.111015182213791    ,
                23.360111410592239    ,
                -0.62805478656865732    ,
                -19.413739122907394    ,
                48.727164158213192    ,
                -0.74904012830686639    ,
                -2.5744170929770003    ,
                13.268740078566783    ,
                2.3413840865438709    ,
                -8.1345865170999332    ,
                -167.79611381137644    ,
                3.1014721308458424    ,
                16.295734951861554    ,
                48.607258542185086    ,
                -4.4037945882815768    ,
                13.077361401481163    ,
                88.0044692825881    ,
                -1.8260883603937252    ,
                2.7292007261331395    ,
                23.345513379564196    ,
                1.298264842541768    ,
                7.3909901792263284    ,
                -93.8527141477876    ,
                0.0094726305208331328    ,
                30.370674862028892    ,
                38.844910319781562    ,
                -4.0220911303029929    ,
                19.515645895531911    ,
                58.751956360217626    ,
                -1.3553024649973651    ,
                3.4691183046398004    ,
                14.947130025132434    ,
                0.08735789234790528    ,
                3.3728097835975759    ,
                -13.040761864928726    ,
                -0.3607396504444938    ,
                7.9596527225602651    ,
                7.5738362204668706    ,
                -0.85478639873529116    ,
                4.8280325879039321    ,
                9.8833413070180534    ,
                -0.24835402528185749    ,
                0.80938586829566861    ,
                2.412628392331706    ,
            };

            for (var j = 0; j < expectedForces.Length; j++)
            {
                Assert.True(Utilities.AreValuesEqual(expectedForces[j], f[j],
                    Tolerance));
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
    }
}