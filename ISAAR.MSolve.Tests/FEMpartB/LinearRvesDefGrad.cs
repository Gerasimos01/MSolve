using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public static class LinearRvesDefGrad
    {

        public static double lerance = 1e-5;

        [Fact]
        public static void CheckDEfGradScaleTransitions2DAndMicrostructure()
        {
            //Origin: linearRves 
            //alllages: strain history


            double E_disp = 20685000; /*Gpa*/ double ni_disp = 0.3; // stather Poisson

            IdegenerateRVEbuilder homogeneousRveBuilder2 = new HomogeneousRVEBuilderNonLinearAndDegenerate()
            {
                Young_s_Modulus = E_disp,
                Poisson_s_Ration = ni_disp
            };


            var microstructure2DdefGrad = new MicrostructureDefGrad2D(homogeneousRveBuilder2,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);

            //microstructure3.UpdateMaterial(new double[4] { 1.01,1, 0.01,0 });
            //.var consCHECK = microstructure2DdefGrad.ConstitutiveMatrix;
            //materialDevelop.UpdateMaterial(new double[] { F_rve[0, 0], F_rve[1, 1], F_rve[0, 1], F_rve[1, 0] });
            microstructure2DdefGrad.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748, 0 });

            var Aijkl = microstructure2DdefGrad.ConstitutiveMatrix.CopytoArray2D();

            double[,] FPKrve = new double[2, 2] { { microstructure2DdefGrad.Stresses[0], microstructure2DdefGrad.Stresses[2] }, { microstructure2DdefGrad.Stresses[3], microstructure2DdefGrad.Stresses[1] } };

            var material = new ShellElasticMaterial2DtransformationbDefGrad()
            {
                YoungModulus = E_disp,
                PoissonRatio= ni_disp
            };

            material.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748, 0 });

            var Aijkl_target = material.ConstitutiveMatrix.CopytoArray2D();

            double[,] FPKrve_target = new double[2, 2] { { material.Stresses[0], material.Stresses[2] }, { material.Stresses[3], material.Stresses[1] } };

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesCheck3, stressesCheck4));
            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(new double[3] { 2 * stressesCheck3[0], 2 * stressesCheck3[1], 2 * stressesCheck3[2] },
            //                                                                stressesCheck5));
            //Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), consCheck1));
            Assert.True(AreDisplacementsSame(Aijkl, Aijkl_target));
            Assert.True(AreDisplacementsSame(FPKrve, FPKrve_target));

        }

        [Fact]
        public static void CheckDEfGradScaleTransitions2DAndMicrostructure_strain_load_case_2()
        {
            //Origin: linearRves 
            //alllages: strain history


            double E_disp = 20685000; /*Gpa*/ double ni_disp = 0.3; // stather Poisson

            IdegenerateRVEbuilder homogeneousRveBuilder2 = new HomogeneousRVEBuilderNonLinearAndDegenerate()
            {
                Young_s_Modulus = E_disp,
                Poisson_s_Ration = ni_disp
            };


            var microstructure2DdefGrad = new MicrostructureDefGrad2D(homogeneousRveBuilder2,
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);

            //microstructure3.UpdateMaterial(new double[4] { 1.01,1, 0.01,0 });
            //.var consCHECK = microstructure2DdefGrad.ConstitutiveMatrix;
            //materialDevelop.UpdateMaterial(new double[] { F_rve[0, 0], F_rve[1, 1], F_rve[0, 1], F_rve[1, 0] });
            microstructure2DdefGrad.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748, 0 });
            microstructure2DdefGrad.SaveState();
            microstructure2DdefGrad.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748*2, 0 });


            var Aijkl = microstructure2DdefGrad.ConstitutiveMatrix.CopytoArray2D();

            double[,] FPKrve = new double[2, 2] { { microstructure2DdefGrad.Stresses[0], microstructure2DdefGrad.Stresses[2] }, { microstructure2DdefGrad.Stresses[3], microstructure2DdefGrad.Stresses[1] } };

            var material = new ShellElasticMaterial2DtransformationbDefGrad()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp
            };

            material.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748, 0 });
            material.SaveState();
            material.UpdateMaterial(new double[] { 1.0005174997415096, 0.99948276378867351, 0.101916564938067748 * 2, 0 });

            var Aijkl_target = material.ConstitutiveMatrix.CopytoArray2D();

            double[,] FPKrve_target = new double[2, 2] { { material.Stresses[0], material.Stresses[2] }, { material.Stresses[3], material.Stresses[1] } };

            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(stressesCheck3, stressesCheck4));
            //Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(new double[3] { 2 * stressesCheck3[0], 2 * stressesCheck3[1], 2 * stressesCheck3[2] },
            //                                                                stressesCheck5));
            //Assert.True(BondSlipTest.AreDisplacementsSame(Matrix1.CopyToArray2D(), consCheck1));
            Assert.True(AreDisplacementsSame(Aijkl, Aijkl_target));
            Assert.True(AreDisplacementsSame(FPKrve, FPKrve_target));

        }


        public static bool AreDisplacementsSame(double[,] expectedValues,
            IMatrixView computedValues)
        {
            var comparer = new ValueComparer(lerance);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public static bool AreDisplacementsSame(double[,] expectedValues,
            double[,] computedValues)
        {
            var comparer = new ValueComparer(lerance);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(1); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1, i2], computedValues[i1, i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

    }
}
