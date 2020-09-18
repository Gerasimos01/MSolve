﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using System;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ConvectionDiffusionStructuralMulti1D9Hexa
    {
        private const int subdomainID = 0;
        private static DenseMatrixSolver.Builder builder = new DenseMatrixSolver.Builder();
        private static Dictionary<int, IVector>[] Solutions;
        private static Dictionary<int, IVector> Accelerations;
        private static Dictionary<int, IVector> Velocities;
        private static Dictionary<int, IVector> Displacements;

        [Fact]
        private static void RunTest()
        {
            //var models = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1 };
            //var modelReaders = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2 };

            //var modelTuples = new[] { CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0), CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0) };
            var modelTuples = new[] { CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0), CreateStructuralModel(10e3, 0,
                new DynamicMaterial(1,25,0,true), 0, 1) };
            var models = new[] { modelTuples[0].Item1, modelTuples[1].Item1 };
            var modelReaders = new[] { modelTuples[0].Item2, modelTuples[1].Item2 };
            //IVectorView[] solutions = SolveModels(models, modelReaders);
            IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);
            Assert.True(CompareResults(solutions[0]));
        }

        private static void UpdateModels(Dictionary<int, IVector>[] solutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
            IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            Solutions = solutions;
            double[] sol0 = solutions[0][0].CopyToArray();
            double[] sol1 = solutions[1][0].CopyToArray();
            modelsToReplace[0] = CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1;
            modelsToReplace[1] = CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1;

            for (int i = 0; i < modelsToReplace.Length; i++)
            {
                solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
                providersToReplace[i] = new ProblemConvectionDiffusion2((Model)modelsToReplace[i], solversToReplace[i]);
                childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
            }
        }

        private static void UpdateNewmarkModel(Dictionary<int, IVector> accelerations, Dictionary<int, IVector> velocities, Dictionary<int, IVector> displacements, IStructuralModel[] modelsToReplace,
            ISolver[] solversToReplace, IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            double[] sol0 = Solutions[0][0].CopyToArray();
            double[] sol1 = Solutions[1][0].CopyToArray();
            IDynamicMaterial commonDynamicMaterialProperties = new DynamicMaterial(1, 25, 0, true);
            Accelerations = accelerations;
            Velocities = velocities;
            Displacements = displacements;
            modelsToReplace[0] = CreateStructuralModel(10e3, 0, new DynamicMaterial(1, 25, 0, true), 0, sol1[10]).Item1;
            solversToReplace[0] = builder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
        }

        private static bool CompareResults(IVectorView solution)
        {
            var comparer = new ValueComparer(1E-5);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
            int numFreeDofs = 6;
            if (solution.Length != 6) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }

        private static Tuple<Model, IModelReader> CreateConvectionDiffusionModel(double k, double[] U, double L, double b1, double b2, double f1, double f2)
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            var modelReader = new ComsolMeshReader2(filename, k, U, L);
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var flux1 = new FluxLoad(f1);
            var flux2 = new FluxLoad(f2);
            var dir1 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, b1);
            });
            var dir2 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, b2);
            });
            var weakDirichlet1 = new WeakDirichlet(dir1, k);
            var weakDirichlet2 = new WeakDirichlet(dir2, k);

            var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var fluxFactory1 = new SurfaceLoadElementFactory(flux1);
            var fluxFactory2 = new SurfaceLoadElementFactory(flux2);
            var boundaryFactory3D = new SurfaceBoundaryFactory3D(0,
                new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));


            int[] boundaryIDs = new int[] { 0, };
            int QuadID = model.ElementsDictionary.Count + 1;
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(fluxElement1);
                    var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(dirichletElement1);
                    var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    var element = new Element();
                    element.ID = QuadID;
                    element.ElementType = SurfaceBoundaryElement;
                    model.SubdomainsDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(QuadID, element);
                    foreach (Node node in nodes)
                    {
                        element.AddNode(node);
                    }
                    QuadID += 1;
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(fluxElement1);
                    var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(dirichletElement2);
                    var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    var element = new Element();
                    element.ID = QuadID;
                    element.ElementType = SurfaceBoundaryElement;
                    model.SubdomainsDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(QuadID, element);
                    foreach (Node node in nodes)
                    {
                        element.AddNode(node);
                    }
                    QuadID += 1;
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateStructuralModel(double C1, double C2, IDynamicMaterial commonDynamicMaterialProperties, double b, double p)
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            var modelReader = new ComsolMeshReader1(filename, 10e3, 0, commonDynamicMaterialProperties);
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var pressure = new PressureLoad(p);
            //var flux2 = new FluxLoad(f2);
            //var dir1 = new DirichletDistribution(list => {
            //    return Vector.CreateWithValue(list.Count, b1);
            //});
            //var dir2 = new DirichletDistribution(list => {
            //    return Vector.CreateWithValue(list.Count, b2);
            //});
            //var weakDirichlet1 = new WeakDirichlet(dir1, k);
            //var weakDirichlet2 = new WeakDirichlet(dir2, k);

            //var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            //var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var pressureFactory = new SurfaceLoadElementFactory(pressure);
            //var fluxFactory2 = new SurfaceLoadElementFactory(flux2);
            //var boundaryFactory3D = new SurfaceBoundaryFactory3D(0,
            //    new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));


            int[] boundaryIDs = new int[] { 0, };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationX
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationY
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationZ
                        });
                    }
                }
            }
            boundaryIDs = new int[] { 5 };
            int QuadID = model.ElementsDictionary.Count + 1;
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    var pressureElement = pressureFactory.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(pressureElement);
                    //var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(dirichletElement2);
                    //var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    var element = new Element();
                    element.ID = QuadID;
                    //element.ElementType = SurfaceBoundaryElement;
                    model.SubdomainsDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(QuadID, element);
                    foreach (Node node in nodes)
                    {
                        element.AddNode(node);
                    }
                    QuadID += 1;
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }


        private static IVectorView[] SolveModels(Model[] models, IModelReader[] modelReaders)
        {
            Vector[] initialTemps = new Vector[models.Length];
            var temp0 = new Dictionary<int, double[]>();
            for (int i = 0; i < models.Length; i++)
            {
                double[] t0 = new double[models[i].Nodes.Count];
                temp0.Add(i, t0);
            }
            int[] boundaryIDs = new int[] { 0 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReaders[0].quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[0][node.ID] = 1;
                        temp0[1][node.ID] = 1;
                    }
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReaders[0].quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[0][node.ID] = 0;
                        temp0[1][node.ID] = 0;
                    }
                }
            }

            DenseMatrixSolver[] solvers = new DenseMatrixSolver[models.Length];
            IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                initialTemps[i] = Vector.CreateFromArray(temp0[i]);
                //var builder = new DenseMatrixSolver.Builder();
                builder.IsMatrixPositiveDefinite = false;
                solvers[i] = builder.BuildSolver(models[i]);
                providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers, providers, childAnalyzers, 2.5e-3, 10, initialTemperature: initialTemps);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }

        private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
        {
            Vector[] initialTemps = new Vector[models.Length];
            var temp0 = new Dictionary<int, double[]>();
            for (int i = 0; i < models.Length; i++)
            {
                double[] t0 = new double[models[i].Nodes.Count];
                temp0.Add(i, t0);
            }
            int[] boundaryIDs = new int[] { 0 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReaders[0].quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[0][node.ID] = 1;
                        temp0[1][node.ID] = 1;
                    }
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReaders[0].quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[0][node.ID] = 0;
                        temp0[1][node.ID] = 0;
                    }
                }
            }

            DenseMatrixSolver[] solvers = new DenseMatrixSolver[models.Length];
            IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                initialTemps[i] = Vector.CreateFromArray(temp0[i]);
                //var builder = new DenseMatrixSolver.Builder();
                builder.IsMatrixPositiveDefinite = false;
                solvers[i] = builder.BuildSolver(models[i]);
                providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            const double timestep = 2.5e-3;
            const double time = 10;
            var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
                providers, childAnalyzers, timestep, time, initialTemperature: initialTemps);
            parentAnalyzer.Initialize();

            var structuralModel = CreateStructuralModel(10e3, 0, new DynamicMaterial(1, 25, 0, true), 0, 1).Item1; // new Model();
            var structuralSolver = builder.BuildSolver(structuralModel);
            var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
            //var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
            var increments = 2;
            var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
            structuralChildAnalyzerBuilder.ResidualTolerance = 1E-8;
            structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
            var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
                structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);

            for (int i = 0; i < time / timestep; i++)
            {
                parentAnalyzer.SolveTimestep(i);
                structuralParentAnalyzer.SolveTimestep(i);
            }

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }
    }
}