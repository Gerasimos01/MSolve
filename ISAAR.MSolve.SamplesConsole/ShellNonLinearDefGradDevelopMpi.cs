using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Analyzers.ObjectManagers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Postprocessing;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Distributed;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MSAnalysis.remoteMatImplementations;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class ShellNonLinearDefGradDevelopMpi
    {
        //origin: Hexa8NonLinearCantileverDefGradDevelop4Multiscale
        //changes: use of multiscale shell formulation
        private const int subdomainID = 0;
                    
        public static void ParallelNonLinearCantilever(int numProcesses)
        {
                   
            int numSubdomains = numProcesses;
            var procs = ProcessDistribution.CreateDistribution(numProcesses, numSubdomains); // FetiDPDofSeparatorMpiTests .CreateModelAndDofSeparator

            Console.Write("thread sleeping for sychronization");
            //Console.Write($"waiting time = " + procs.OwnRank*20000);
            System.Threading.Thread.Sleep(/*(procs.OwnRank+1)**/20000);


            #region material choice
            //var material = new ShellElasticMaterial2Dtransformationb()
            //{
            //    YoungModulus = 100,
            //    PoissonRatio = 0
            //};
            //var modelReader = new IsogeometricShellReader(GeometricalFormulation.NonLinear, filename, material: material);

            //var material = new ShellElasticMaterial2DtransformationbDefGrad()
            //{
            //    YoungModulus = 100,
            //    PoissonRatio = 0
            //};

            IdegenerateRVEbuilder homogenousRve = new HomogeneousRVEBuilderNonLinearAndDegenerate()
            {
                Young_s_Modulus = 100,
                Poisson_s_Ration = 0
            };
            var material1 = new MicrostructureDefGrad2D(homogenousRve,
                model1 => (new SkylineSolver.Builder()).BuildSolver(model1), false, 1);
            #endregion
            IMaterialManager materialManager = new MaterialManagerMpi2(material1, procs);

            Model model = new Model();//'

            if(procs.IsMasterProcess)
            {
                BuildIsogeometricCantileverShell(ref model, materialManager);
            }


            var increments = 2;//.
            var resTol = 1E-3;// sto multiscale xrhsimopoihsame thn default.
            int maxIters = 100;
            int itersRebuild = 1;


            (StaticAnalyzerDevelopMpi parentAnalyzer, LoadControlAnalyzerDevelop4Mpi childAnalyzer) = GetAnalyzers(procs, model, increments, maxIters,
                itersRebuild, resTol, materialManager);
            


            materialManager.Initialize();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            if (procs.IsMasterProcess)
            {
                var log1 = childAnalyzer.TotalDisplacementsPerIterationLog;//
                IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
                bool isProblemSolvedCorrectly = AreDisplacementsSame(expectedDisplacements, log1);
                if (isProblemSolvedCorrectly)
                {
                    Console.WriteLine($"Problem is solved correctly ");
                }
                else
                {
                    Console.WriteLine($"the problem has not been solved correctly");
                }

            }


        }

        private static (StaticAnalyzerDevelopMpi parentAnalyzer, LoadControlAnalyzerDevelop4Mpi childAnalyzer) GetAnalyzers(ProcessDistribution procs,
            Model model, int increments, int maxIters, int itersRebuild, double resTol, IMaterialManager materialManager)
        {
            if (procs.IsMasterProcess)
            {
                // Solver
                var solverBuilder = new SkylineSolver.Builder();
                ISolver solver = solverBuilder.BuildSolver(model);
                // Problem type
                var provider = new ProblemStructural(model, solver);
                var childAnalyzer = new LoadControlAnalyzerDevelop4Mpi(model, solver, provider, increments, maxIters, itersRebuild, resTol, materialManager, procs);
                var parentAnalyzer = new StaticAnalyzerDevelopMpi(model, solver, provider, childAnalyzer, procs);
                var watchDofs = new Dictionary<int, int[]>();
                watchDofs.Add(subdomainID, new int[1] { 17 });
                var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
                childAnalyzer.TotalDisplacementsPerIterationLog = log1; //.

                return (parentAnalyzer, childAnalyzer);
            }
            else
            {
                var childAnalyzer = new LoadControlAnalyzerDevelop4Mpi(materialManager, procs, increments, maxIters, itersRebuild);
                var parentAnalyzer = new StaticAnalyzerDevelopMpi(childAnalyzer, procs);

                return (parentAnalyzer, childAnalyzer);
            }
        }

        

        public static void BuildIsogeometricCantileverShell(ref Model model , IMaterialManager materialManager)
        {
            IContinuumMaterial3DDefGrad material = new RemoteMaterial(materialManager);
            #region choose model
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "CantileverShell.txt");

            var modelReader = new IsogeometricShellReader(GeometricalFormulation.DefGrad, filename, defGradMaterial: material);
            model = modelReader.GenerateModelFromFile();

            
            for (int i = 0; i < 6; i++)
            {
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationX });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationY });
                model.ControlPointsDictionary[i].Constraints.Add(new Constraint() { DOF = StructuralDof.TranslationZ });
            }

            //TODO: add logs to loadcontrol analyzer

            var logger = new TotalLoadsDisplacementsPerIncrementLogIGA(model.PatchesDictionary[0], 20, //1000,
                model.ControlPointsDictionary.Values.Last(), StructuralDof.TranslationZ, "CantileverBenchmarkLog16x1.txt");
            //childAnalyzer.IncrementalLogs.Add(0, logger);

            //model.Patches[0].Forces[11] = 1.3333333333333344;//..Forces[2] = 1.33333333334;
            //model.Patches[0].Forces[14] = 1.3333333333333335;//..Forces[5] = 1.33333333334;
            //model.Patches[0].Forces[17] = 1.3333333333333344;//..Forces[8] = 1.33333333334;            
            model.Loads.Add(new Load() { Node = model.GetNode(9), DOF = StructuralDof.TranslationZ, Amount = 1.3333333333333344 });
            model.Loads.Add(new Load() { Node = model.GetNode(10), DOF = StructuralDof.TranslationZ, Amount = 1.3333333333333335 });
            model.Loads.Add(new Load() { Node = model.GetNode(11), DOF = StructuralDof.TranslationZ, Amount = 1.3333333333333344 });

            #endregion
        }


        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, TotalDisplacementsPerIterationLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }



        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[9]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES            

            expectedDisplacements[0]  = new Dictionary<int, double> {{ 17,999.999996814643}, };
            expectedDisplacements[1]  = new Dictionary<int, double> {{ 17,667.045140505124}, };
            expectedDisplacements[2]  = new Dictionary<int, double> {{ 17,445.26607807662066}, };
            expectedDisplacements[3]  = new Dictionary<int, double> {{ 17,297.69977429659536}, };
            expectedDisplacements[4]  = new Dictionary<int, double> {{ 17,199.75236602177949}, };
            expectedDisplacements[5]  = new Dictionary<int, double> {{ 17,135.10107181694747}, };
            expectedDisplacements[6]  = new Dictionary<int, double> {{ 17,92.975524369045274}, };
            expectedDisplacements[7]  = new Dictionary<int, double> {{ 17,66.359166004254561}, };
            expectedDisplacements[8]  = new Dictionary<int, double> {{ 17,50.739928431320266}, };
            expectedDisplacements[9]  = new Dictionary<int, double> {{ 17,42.651539056979033}, };
            expectedDisplacements[10] = new Dictionary<int, double> {{ 17,37.5409811863373}, };
            expectedDisplacements[11] = new Dictionary<int, double> {{ 17,36.27396142265728}, };
            expectedDisplacements[12] = new Dictionary<int, double> {{ 17,36.568414624233107}, };
            expectedDisplacements[13] = new Dictionary<int, double> {{ 17,36.68237139312906}, };
            expectedDisplacements[14] = new Dictionary<int, double> {{ 17,36.725607667922112}, };
            expectedDisplacements[15] = new Dictionary<int, double> {{ 17,37.854378438146036}, };
            expectedDisplacements[16] = new Dictionary<int, double> {{ 17,37.505658524669116}, };
            expectedDisplacements[17] = new Dictionary<int, double> {{ 17,37.399100491285296}, };
            expectedDisplacements[18] = new Dictionary<int, double> {{ 17,37.376274802449359}, };
            expectedDisplacements[19] = new Dictionary<int, double> {{ 17,37.747726005757194}, };
            expectedDisplacements[20] = new Dictionary<int, double> {{ 17,37.694915696859866}, };
            expectedDisplacements[21] = new Dictionary<int, double> {{ 17,37.93214994905761}, };
            expectedDisplacements[22] = new Dictionary<int, double> {{ 17,37.922430949221038}, };
            expectedDisplacements[23] = new Dictionary<int, double> {{ 17,38.130105761482021}, };
            expectedDisplacements[24] = new Dictionary<int, double> {{ 17,38.11987712590426}, };
            expectedDisplacements[25] = new Dictionary<int, double> {{ 17,38.304979201650404}, };
            expectedDisplacements[26] = new Dictionary<int, double> {{ 17,38.300116410916885}, };
            expectedDisplacements[27] = new Dictionary<int, double> {{ 17,38.473881011594933}, };
            expectedDisplacements[28] = new Dictionary<int, double> {{ 17,38.647645612272981}, };
            expectedDisplacements[29] = new Dictionary<int, double> {{ 17,38.634267016110783}, };
            expectedDisplacements[30] = new Dictionary<int, double> {{ 17,38.793904106013088}, };
            expectedDisplacements[31] = new Dictionary<int, double> {{ 17,38.953541195915392}, };
            expectedDisplacements[32] = new Dictionary<int, double> {{ 17,38.947944723937113}, };
            expectedDisplacements[33] = new Dictionary<int, double> {{ 17,39.101044345246848}, };
            expectedDisplacements[34] = new Dictionary<int, double> {{ 17,39.254143966556583}, };
            expectedDisplacements[35] = new Dictionary<int, double> {{ 17,39.249273570902687}, };
            expectedDisplacements[36] = new Dictionary<int, double> {{ 17,39.397206280584371}, };
            expectedDisplacements[37] = new Dictionary<int, double> {{ 17,39.545138990266054}, };
            expectedDisplacements[38] = new Dictionary<int, double> {{ 17,39.54119857973707}, };
            expectedDisplacements[39] = new Dictionary<int, double> {{ 17,39.684998533568816}, };
            expectedDisplacements[40] = new Dictionary<int, double> {{ 17,39.828798487400562}, };
            expectedDisplacements[41] = new Dictionary<int, double> {{ 17,39.8253350771238}, };
            expectedDisplacements[42] = new Dictionary<int, double> {{ 17,39.96554957357349}, };
            expectedDisplacements[43] = new Dictionary<int, double> {{ 17,40.105764070023177}, };
            expectedDisplacements[44] = new Dictionary<int, double> {{ 17,40.102636823288712}, };
            expectedDisplacements[45] = new Dictionary<int, double> {{ 17,40.239631781858641}, };
            expectedDisplacements[46] = new Dictionary<int, double> {{ 17,40.37662674042857}, };
            expectedDisplacements[47] = new Dictionary<int, double> {{ 17, 40.37374531856247 }, };
            
            return expectedDisplacements;
        }



    }

}
