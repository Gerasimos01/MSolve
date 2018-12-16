using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.FEM.Materials;

namespace ISAAR.MSolve.SamplesConsole
{
    class SeparateCodeCheckingClass
    {
        public static void Check01()
        {
            double[,] Aijkl = new double[9, 9];
            for (int k = 0; k < 9; k++)
            {
                for (int l = 0; l < 9; l++)
                {
                    Aijkl[k, l] = 100 * (k+1) + (l+1);
                }
            }

            double[,] SPK = new double[3, 3] { { 0.5, -0.1, -0.0 },{ -0.1, 0.2, 0.1 },{ 0, 0.1, 0.15 } } ;
            double[,] F = new double[3, 3] { { 1.2, -0.1, -0.0 }, { -0.1,1.1, 0.1 }, { 0, 0.1, 1.15 } };

            double[,] Cinpk = Transform_d2Wdfdf_to_Cijrs(Aijkl, SPK, F);
        }

        public static void Check02()
        {
            double[,] F__F__ = new double[9, 9];
            

            for (int i1 = 0; i1 < 9; i1++)
            {
                F__F__[i1, i1] = 1;
            }

            Vector solution = new Vector(new double[9]);

            SkylineMatrix2D F__F__Mat = new SkylineMatrix2D(F__F__);
            int linearsystemID = 1;
            SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[9]);
            var solver = new SolverSkyline(linearSystem);
            // BuildMatrices();
            linearSystem.Matrix = F__F__Mat;
            //solver.Initialize();
            solver.Initialize(); // dld factorize

            double[] rhs = new double[9];
            for (int l = 0; l < 9; l++)
            {
                rhs[l] = 10;
            }
            Vector RHS = new Vector(rhs);
            SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
            k.Solve(RHS, solution);

        }

        public static void Check03()
        {
            Dictionary<int, double> dictionary1 = new Dictionary<int, double>();
            dictionary1.Add(4, 0.111112);
            dictionary1.Add(8, 0.222223);
            var dok1 = dictionary1.Keys;
            var dok2 = dictionary1.GetEnumerator();
        }

        public static void Check04integration()
        {
            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilder();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheck(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3copyConsCheck.UpdateMaterial(new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 });

        }

        public static void Check05aStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp,PoissonRatio = ni_disp,};
            double[,] DGtr = new double[3, 3] { { 1.01, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D( GLVec));
            double[] stressesCheck1 = material1.Stresses.Data;
            DGtr = new double[3, 3] { { 1.02, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilder();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.01, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
        }

        public static void Check05bStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27Hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });            
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void Check05cStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial2D( StressState2D.PlaneStress ) 
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck1 = new double[3] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2]};
            material1.SaveState();
            GLVec = new double[3] { 0.02, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck2 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27HexaDegenerateAndLinear();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial2D microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3,3];
            for (int i1 = 0; i1 <3 ; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D( new double[3] { 0.010, 0, 0}));
            double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.020, 0, 0 }));
            double[] stressesCheck4 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };

        }

        public static void Check05c2StressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };
            //material1.SaveState();
            GLVec = new double[3] { 0.02, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck2 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27HexaDegenerateAndLinearPeripheral();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial2D microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.010, 0, 0 }));
            double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.020, 0, 0 }));
            double[] stressesCheck4 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };

            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.030, 0, 0 }));
            double[] stressesCheck5 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };


        }

        public static void Check05dStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };
            //material1.SaveState();
            //GLVec = new double[3] { 0.02, 0, 0 };
            //material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            //double[] stressesCheck2 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVEstifDegenAndLinear();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial2D microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.010, 0, 0 }));
            double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.020, 0, 0 }));
            double[] stressesCheck4 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };

        }

        public static void Check05d2StressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };
            //material1.SaveState();
            //GLVec = new double[3] { 0.02, 0, 0 };
            //material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            //double[] stressesCheck2 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            VectorExtensions.AssignTotalAffinityCount();
            IdegenerateRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVEstifDegenAndLinearPeripheral();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial2D microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.010, 0, 0 }));
            double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new StressStrainVectorContinuum2D(new double[3] { 0.020, 0, 0 }));
            double[] stressesCheck4 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };

        }

        public static void Check05eStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial2D(StressState2D.PlaneStress)
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[] GLVec = new double[3] { 0.01, 0, 0 };
            material1.UpdateMaterial(new StressStrainVectorContinuum2D(GLVec));
            double[] stressesCheck1 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };

            var Vec1 = new Vector(new double[3] { 2, 0, 1 });
            var Vec2 = new Vector(new double[3] { 1, 2, 1 });
            var strain = new double[3] { 0.01, 0, 0 };

            var material2 = new ShellElasticMaterial2D() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1=new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            material2.UpdateMaterial(strain);
            double[] stressesCheck2 = new double[3] { material2.Stresses[0], material2.Stresses[1], material2.Stresses[2] };

            var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } };
            var Matrix1 = new Matrix2D(3,3); for (int i1 = 0; i1 < 3; i1++) { for (int i2 = 0; i2 < 3; i2++) { Matrix1[i1, i2] = material3.ConstitutiveMatrix[i1, i2]; } }
            material3.UpdateMaterial(strain);
            double[] stressesCheck3 = new double[3] { material3.Stresses[0], material3.Stresses[1], material3.Stresses[2] };


        }



        public static void Check06bStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECK();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) {consCheck1[i1,i2]= microstructure3.ConstitutiveMatrix[i1,i2]; } }
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void Check07bStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroE();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void Check08bStressIntegration()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroE8hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void Check09bStressIntegrationStrainFE2()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroEstrainfe2();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1,0.995548518728336,1.03039263387409,-5.18077793066249e-17,0.00608083650156188,-0.0222489521716709,0.0667140176791647,-0.00445148127166406,6.68129711107579e-16 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            microstructure3.SaveState();
            //microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            //double[] stressesCheck4 = microstructure3.Stresses;
        }

        public static void Check10bStressIntegration5GrSh1RveForDemonstration()
        {

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVE();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            int cycles = 4;
            double[] totalStrain = new double[9] { 1.1, 1, 1, -5.18077793066249e-17, 0, 0, 0, 0, 0 };
            double[][] strainHistory = new double[cycles][];
            double[] initialDefGrad = new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 };
            double[] strainIncr = new double[totalStrain.Length]; for (int l = 0; l < totalStrain.Length; l++) { strainIncr[l] = (totalStrain[l] - initialDefGrad[l]) / ((double)cycles); }
            for (int l = 0; l < cycles; l++)
            {
                strainHistory[l] = new double[totalStrain.Length];
                for (int i1 = 0; i1 < totalStrain.Length; i1++) { strainHistory[l][i1] = initialDefGrad[i1] + (strainIncr[i1] * ((double)l + 1)); }
            }

            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(strainHistory, microstructure3);



            //microstructure3.UpdateMaterial(new double[9] { 1, 0.995548518728336, 1.03039263387409, -5.18077793066249e-17, 0.00608083650156188, -0.0222489521716709, 0.0667140176791647, -0.00445148127166406, 6.68129711107579e-16 });
            //double[] stressesCheck3 = microstructure3.Stresses;
            //double[,] consCheck1 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            //microstructure3.SaveState();

            //microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            //double[] stressesCheck4 = microstructure3.Stresses;
        }

        public static void Check10cStressIntegration5GrSh1RveForDemonstration()
        {
            //proelefsi Check10bStressIntegration5GrSh1RveForDemonstration
            //allages to 1)strainhistory, 2)print to output kai 3)parametroi cohesive
            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVEstif();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            int cycles = 10;
            double[] totalStrain = new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 };
            double[] initialDefGrad = new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 };
            double[][] strainHistory1 = CreateStrainHistoryInterpolation(totalStrain, initialDefGrad, cycles);
            double[][] strainHistory2 = CreateStrainHistoryInterpolation(initialDefGrad, totalStrain, cycles);
            double[][] strainHistory = CombineStrainHistoriesIntoOne(strainHistory1, strainHistory2);



            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(strainHistory, microstructure3);

            double[,] strainHistoryArray = ConvertStressHistoryTodoubleArray(strainHistory);
            double[,] stressHistoryArray = ConvertStressHistoryTodoubleArray(stressHistory);

            PrintUtilities.WriteToFile(strainHistoryArray, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\strainHistory.txt");
            PrintUtilities.WriteToFile(stressHistoryArray, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressHistory.txt");


            //microstructure3.UpdateMaterial(new double[9] { 1, 0.995548518728336, 1.03039263387409, -5.18077793066249e-17, 0.00608083650156188, -0.0222489521716709, 0.0667140176791647, -0.00445148127166406, 6.68129711107579e-16 });
            //double[] stressesCheck3 = microstructure3.Stresses;
            //double[,] consCheck1 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            //microstructure3.SaveState();

            //microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            //double[] stressesCheck4 = microstructure3.Stresses;
        }

        public static void CheckProxeiroStressElastic()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            double[,] DGtr = new double[3, 3] { { 1.05, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
        }

        public static void Check10bStressIntegration0GrSh1RveForDemonstration()
        {

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroE8hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            int cycles = 4;
            double[] totalStrain = new double[9] { 1, 0.995548518728336, 1.03039263387409, -5.18077793066249e-17, 0.00608083650156188, -0.0222489521716709, 0.0667140176791647, -0.00445148127166406, 6.68129711107579e-16 };
            double[][] strainHistory = new double[cycles][];
            double[] initialDefGrad = new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 };
            double[] strainIncr = new double[totalStrain.Length]; for (int l = 0; l < totalStrain.Length; l++) { strainIncr[l] = (totalStrain[l] - initialDefGrad[l]) / ((double)cycles); }
            for (int l = 0; l < cycles; l++)
            {
                strainHistory[l] = new double[totalStrain.Length];
                for (int i1 = 0; i1 < totalStrain.Length; i1++) { strainHistory[l][i1] = initialDefGrad[i1] + (strainIncr[i1] * ((double)l + 1)); }
            }

            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(strainHistory, microstructure3);



            //microstructure3.UpdateMaterial(new double[9] { 1, 0.995548518728336, 1.03039263387409, -5.18077793066249e-17, 0.00608083650156188, -0.0222489521716709, 0.0667140176791647, -0.00445148127166406, 6.68129711107579e-16 });
            //double[] stressesCheck3 = microstructure3.Stresses;
            //double[,] consCheck1 = new double[6, 6];
            //for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            //microstructure3.SaveState();

            //microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            //double[] stressesCheck4 = microstructure3.Stresses;
        }

        public static void Check11bStressIntegrationMultipleSubdomainMethods()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.05, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroE8hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3DevelopMultipleSubdomains(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void Check12bStressIntegrationMultipleSubdomainMethodsUseBase()
        {
            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.05, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            //IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample1GrSh();
            IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderCHECKzeroE8hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBase(homogeneousRveBuilder1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);

            microstructure3.UpdateMaterial(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }

        public static void CheckPartitionLimitsArray()
        {
            int[] ch01 = new int[10];

            var ch02 = ch01.PartitionLimits(7);
        }

        public static void CheckListUnionInDictionary()
        {
            List<int> list1 = new List<int>() { 1, 1, 1, 2, 3, 4, 5, 5, 7 };
            List<int> list2 = new List<int>() { 4, 5, 5, 7, 8, 9, 10, 10, 11 };
            Dictionary<int, List<int>> lists = new Dictionary<int, List<int>>();
            lists.Add(1,list1);
            lists.Add(2, list2);

            lists[1]= lists[1].Union(lists[2]).ToList();
            lists.Remove(2);

            List<int> ch01 = new List<int>(10);
            int ch02 = ch01.Capacity;
        }

        public static void CheckStressStrainBonSlipMaterial()
        {
            VectorExtensions.AssignTotalAffinityCount();
            BondSlipCohMat material1 = new BondSlipCohMat(100, 10, 100, 100, new double[2], new double[2], 1e-10);
            int loadsteps_2 = 120;
            double[][] DeltaEhist = new double[2*loadsteps_2][];
            double phi_metakinhshs = ((double)30 /(double) 360) * 2 * Math.PI;
            double[] epsilon_max = new double[3] { 2.8 * Math.Cos(phi_metakinhshs), 2.8 * Math.Sin(phi_metakinhshs), 2.8 };
            for (int i1 = 0; i1 < loadsteps_2; i1++)
            { DeltaEhist[i1] = new double[3] { (1 / (double)loadsteps_2) * epsilon_max[0], (1 / (double)loadsteps_2) * epsilon_max[1], (1 / (double)loadsteps_2) * epsilon_max[2] }; }
            for (int i1 = loadsteps_2; i1 <2* loadsteps_2; i1++)
            { DeltaEhist[i1] = new double[3] { -1.5*(1 / (double)loadsteps_2) * epsilon_max[0], -1.5*(1 / (double)loadsteps_2) * epsilon_max[1], -1.5*(1 / (double)loadsteps_2) * epsilon_max[2] }; }
            double[][] Ehist = new double[2 * loadsteps_2][];
            Ehist[0] = new double[3] { DeltaEhist[0][0], DeltaEhist[0][1], DeltaEhist[0][2] };
            for (int i1 = 1; i1 <2*loadsteps_2; i1++)
            { Ehist[i1] = new double[3] { Ehist[i1-1][0]+DeltaEhist[i1][0], Ehist[i1 - 1][1] + DeltaEhist[i1][1], Ehist[i1 - 1][2] + DeltaEhist[i1][2] }; }


            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(Ehist, material1);

        }

        public static void CheckConstitutiveMatrixTransformation()
        {
            double[,] d2W_dfdf = new double[9, 9]
                {{7.502343750000,5.000000000000,5.000000000000,0.187500000000,-0.000000000000,-0.000000000000,-0.000000000000,0.031250000000,0.000000000000},
                {5.000000000000,7.502343750000,5.000000000000,0.187500000000,-0.000000000000,-0.000000000000,0.000000000000,0.031250000000,0.000000000000},
                {5.000000000000,5.000000000000,7.501562500228,0.125000000000,-0.000000000000,0.000000000000,-0.000000000000,-0.000000000000,0.000000000000},
                {0.187500000000,0.187500000000,0.125000000000,1.257031250000,-0.000000000000,0.000000000000,-0.000000000000,1.250000000000,0.000000000000},
                {0.000000000000,-0.000000000000,0.000000000000,0.000000000000,1.251562500000,-0.000000000000,0.031250000000,-0.000000000000,1.250000000000},
                {0.000000000000,-0.000000000000,0.000000000000,0.000000000000,-0.000000000000,1.251562500000,1.250000000000,0.000000000000,0.031250000000},
                {-0.000000000000,-0.000000000000,0.000000000000,-0.000000000000,0.031250000000,1.250000000000,1.252343750000,-0.000000000000,0.031250000000},
                {0.031250000000,0.031250000000,-0.000000000000,1.250000000000,-0.000000000000,-0.000000000000,0.000000000000,1.251562500000,-0.000000000000},
                {0.000000000000,0.000000000000,0.000000000000,0.000000000000,1.250000000000,0.031250000000,0.031250000000,0.000000000000,1.252343750000},

                };

            double[,] DefGradMat = new double[3, 3]
            {
                {1.000000000000,0.025000000000,0.000000000000},
                {0.000000000000,1.000000000000,0.000000000000},
                {0.000000000000,0.000000000000,1.000000000000},
            };

            double[,] SPK_mat = new double[3, 3]
            {
                {0.001562500000,0.031250000000,0.000000000000},
                {0.031250000000,0.002343750000,0.000000000000},
                {-0.000000000000,0.000000000000,0.001562500000},
            };

            double[,] Cinpk = Transform_d2Wdfdf_to_Cijrs(d2W_dfdf, SPK_mat, DefGradMat); // to onomazoume Cinpk epeidh einai to 9x9 kai to diakrinoume etsi apo to Cijrs 6x6
            // transformation se 6x6 se 2 vhmata
            double[,] Cijrs_columns = new double[9, 6];
            for (int i1 = 0; i1 < 9; i1++)
            {
                Cijrs_columns[i1, 0] = Cinpk[i1, 0];
                Cijrs_columns[i1, 1] = Cinpk[i1, 1];
                Cijrs_columns[i1, 2] = Cinpk[i1, 2];
                Cijrs_columns[i1, 3] = 0.5 * (Cinpk[i1, 3] + Cinpk[i1, 7]);
                Cijrs_columns[i1, 4] = 0.5 * (Cinpk[i1, 4] + Cinpk[i1, 8]);
                Cijrs_columns[i1, 5] = 0.5 * (Cinpk[i1, 5] + Cinpk[i1, 6]);
            }

            double[,] Cijrs = new double[6, 6];

            for (int j1 = 0; j1 < 6; j1++)
            {
                Cijrs[0, j1] = Cijrs_columns[0, j1];
                Cijrs[1, j1] = Cijrs_columns[1, j1];
                Cijrs[2, j1] = Cijrs_columns[2, j1];
                Cijrs[3, j1] = 0.5 * (Cijrs_columns[3, j1] + Cijrs_columns[7, j1]);
                Cijrs[4, j1] = 0.5 * (Cijrs_columns[4, j1] + Cijrs_columns[8, j1]);
                Cijrs[5, j1] = 0.5 * (Cijrs_columns[5, j1] + Cijrs_columns[6, j1]);
            }


        }

        private static double[] Transform_DGtr_to_GLvec(double[,] DGtr)
        {
            double[,] GL = new double[3, 3];

            //
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0;
                    for (int p = 0; p < 3; p++)
                    {
                        GL[m, n] += DGtr[m, p] * DGtr[n, p];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                GL[m, m] += -1;
            }
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0.5 * GL[m, n];
                }
            }

            double[] GLvec = new double[6];
            //
            for (int m = 0; m < 3; m++)
            {
                GLvec[m] = GL[m, m];
            }
            GLvec[3] = 2 * GL[0, 1];
            GLvec[4] = 2 * GL[1, 2];
            GLvec[5] = 2 * GL[2, 0];

            return GLvec;
        }

        public static void Check06CohMatBKstraining()
        {
            //thumizetai 
            //for (int j = 0; j < 3; j++)
            //{
            //    D_tan[j, j] = (1 - d_prev_step) * E;
            //}
            //if (Delta[2] < 0)
            //{
            //    D_tan[2, 2] += d_prev_step * E;
            //}
            //dld sth tlipsi den ephreazei to damage
            //paramaetroi apo 
            //mpgp = RVEkanoninkhsGewmetriasBuilder.GetReferenceKanonikhGewmetriaRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //parametroi cohesive epifaneias
            double T_o_3 = 0.05;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = 0.5; // nm
            double D_f_3 = 4; // nm
            double T_o_1 = 0.05;// Gpa
            double D_o_1 = 0.5; // nm
            double D_f_1 = 4; // nm
            double n_curve = 1.4;

            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };


        }

        private static double[,] Transform_d2Wdfdf_to_Cijrs(double[,] Aijkl, double[,] SPK, double[,] F)
        {
            int[,] i_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };
            int[,] k_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };

            double[,] Cinpk = new double[9, 9];

            double[,] F__F__ = new double[9, 9];
            //[F(:, 1) * F(1,:), F(:, 2) * F(1,:), F(:, 3) * F(1,:);
            //F(:, 1) * F(2,:),F(:, 2) * F(2,:),F(:, 3) * F(2,:);
            //F(:, 1) * F(3,:),F(:, 2) * F(3,:),F(:, 3) * F(3,:)];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            F__F__[3 * i1 + k, 3 * j1 + l] = F[k, j1] * F[i1, l];
                        }
                    }

                }
            }

            //TODO upologismos F_inv (F__F__^(-1))

            double[,] F_inv = new double[9, 9];
            double[,] eye9 = new double[9, 9];
            for (int i1 = 0; i1 < 9; i1++)
            { eye9[i1, i1] = 1; }
            Vector solution = new Vector(new double[9]);

            Matrix2D F__F__Mat = new Matrix2D(F__F__);
            //int linearsystemID = 1;
            //SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[9]);
            //var solver = new SolverSkyline(linearSystem);
            //// BuildMatrices();
            //linearSystem.Matrix = F__F__Mat;
            ////solver.Initialize();
            //solver.Initialize(); // dld factorize

            //for (int j1 = 0; j1 < 9; j1++)
            //{
            //    Vector RHS = new Vector(new double[9] { eye9[0, j1], eye9[1, j1], eye9[2, j1], eye9[3, j1], eye9[4, j1], eye9[5, j1], eye9[6, j1], eye9[7, j1], eye9[8, j1] });
            //    SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
            //    k.Solve(RHS, solution);
            //    for (int i1 = 0; i1 < 9; i1++)
            //    {
            //        F_inv[i1, j1] = solution[i1];
            //    }
            //}

            double[,] multipleRHSs = new double[9, 9];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    double[] A_j_l = new double[9] { Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[k1, 2] - 1) + k1]};

                    double[] sec_term = new double[9] { -SPK[i1, k1], 0, 0, 0, -SPK[i1, k1], 0, 0, 0, -SPK[i1, k1] };

                    //Matrix2D F_invMat = new Matrix2D(F_inv);
                    //Vector A_j_lVec = new Vector(A_j_l);
                    //Vector sec_termVec = new Vector(sec_term); //comment out eginan afta afou den xreiazontai gia ta parakatw

                    int RHScolumn = i1 * 3 + k1;
                    for (int a1 = 0; a1 < 9; a1++)
                    {
                        multipleRHSs[a1, RHScolumn] = A_j_l[a1] + sec_term[a1];
                    }

                    //Vector C_np_ = F_invMat * (new Vector(A_j_lVec + sec_termVec)); //comment out afta

                    //Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[0];
                    //Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[1];
                    //Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[2];
                    //Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[3];
                    //Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[4];
                    //Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[5];
                    //Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[6];
                    //Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[7];
                    //Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[8];

                }
            }

            //TODO right here inversion:
            Matrix2D MultipleSolutions = F__F__Mat.SolveLU(new Matrix2D(multipleRHSs), true);
            //MultipleSolutions=F__F__Mat\multipleRHSs;


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {


                    //Matrix2D F_invMat = new Matrix2D(F_inv);
                    //Vector A_j_lVec = new Vector(A_j_l);
                    //Vector sec_termVec = new Vector(sec_term); //comment out eginan afta afou den xreiazontai gia ta parakatw

                    int RHScolumn = i1 * 3 + k1;
                    //for (int a1 = 0; a1 < 9; a1++)
                    //{
                    //    multipleRHSs[a1, RHScolumn] = A_j_l[a1] + sec_term[a1];
                    //}

                    //Vector C_np_ = F_invMat * (new Vector(A_j_lVec + sec_termVec)); //comment out afta

                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[0, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[1, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[2, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[3, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[4, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[5, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 0] = MultipleSolutions[6, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 1] = MultipleSolutions[7, RHScolumn];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 2] = MultipleSolutions[8, RHScolumn];

                }
            }


            return Cinpk;
        }

        //private static Tuple< double [][], double[] [,]> StressStrainHistoryCopy(double [][] strainHistory, IFiniteElementMaterial3D testedMaterial)
        //{
        //    double[][] stressHistory = new double[strainHistory.GetLength(0)][];
        //    double[][,] constitutiveMatrixHistory = new double[strainHistory.GetLength(0)][,];

        //    for (int l = 0; l < strainHistory.GetLength(0); l++)
        //    {
        //        testedMaterial.UpdateMaterial(strainHistory[l]);
        //        testedMaterial.SaveState();
        //        testedMaterial.Stresses.CopyTo(stressHistory[l],0);
        //        constitutiveMatrixHistory[l] = new double[testedMaterial.ConstitutiveMatrix.Columns, testedMaterial.ConstitutiveMatrix.Rows];

        //        for (int m = 0; m < testedMaterial.ConstitutiveMatrix.Columns; m++)
        //        {
        //            for (int n = 0; n < testedMaterial.ConstitutiveMatrix.Rows; n++)
        //            { constitutiveMatrixHistory[l][m,n] = testedMaterial.ConstitutiveMatrix[m, n]; }
        //        }
        //    }

        //    return new Tuple<double[][], double[][,]>(stressHistory, constitutiveMatrixHistory);
        //}

        private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, ICohesiveZoneMaterial3D testedMaterial)
        {
            double[][] stressHistory = new double[strainHistory.GetLength(0)][];
            double[][,] constitutiveMatrixHistory = new double[strainHistory.GetLength(0)][,];

            for (int l = 0; l < strainHistory.GetLength(0); l++)
            {
                if (l==42)
                {
                    Console.Write("breakPointIsHere");
                }
                testedMaterial.UpdateMaterial(strainHistory[l]);
                testedMaterial.SaveState();
                stressHistory[l] = new double[testedMaterial.Tractions.Length];
                testedMaterial.Tractions.CopyTo(stressHistory[l], 0);
                constitutiveMatrixHistory[l] = new double[testedMaterial.ConstitutiveMatrix.Columns, testedMaterial.ConstitutiveMatrix.Rows];

                for (int m = 0; m < testedMaterial.ConstitutiveMatrix.Columns; m++)
                {
                    for (int n = 0; n < testedMaterial.ConstitutiveMatrix.Rows; n++)
                    { constitutiveMatrixHistory[l][m, n] = testedMaterial.ConstitutiveMatrix[m, n]; }
                }
            }

            return (stressHistory, constitutiveMatrixHistory);
        }

        private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IContinuumMaterial3DDefGrad testedMaterial)
        {
            double[][] stressHistory = new double[strainHistory.GetLength(0)][];
            double[][,] constitutiveMatrixHistory = new double[strainHistory.GetLength(0)][,];

            for (int l = 0; l < strainHistory.GetLength(0); l++)
            {
                if (l == 42)
                {
                    Console.Write("breakPointIsHere");
                }
                testedMaterial.UpdateMaterial(strainHistory[l]);
                testedMaterial.SaveState();
                stressHistory[l] = new double[testedMaterial.Stresses.Length];
                testedMaterial.Stresses.CopyTo(stressHistory[l], 0);
                constitutiveMatrixHistory[l] = new double[testedMaterial.ConstitutiveMatrix.Columns, testedMaterial.ConstitutiveMatrix.Rows];

                for (int m = 0; m < testedMaterial.ConstitutiveMatrix.Columns; m++)
                {
                    for (int n = 0; n < testedMaterial.ConstitutiveMatrix.Rows; n++)
                    { constitutiveMatrixHistory[l][m, n] = testedMaterial.ConstitutiveMatrix[m, n]; }
                }
            }

            return (stressHistory, constitutiveMatrixHistory);
        }

        private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IContinuumMaterial3D testedMaterial)
        {
            double[][] stressHistory = new double[strainHistory.GetLength(0)][];
            double[][,] constitutiveMatrixHistory = new double[strainHistory.GetLength(0)][,];

            for (int l = 0; l < strainHistory.GetLength(0); l++)
            {
                if (l == 42)
                {
                    Console.Write("breakPointIsHere");
                }
                testedMaterial.UpdateMaterial(new StressStrainVectorContinuum3D(strainHistory[l]));
                testedMaterial.SaveState();
                stressHistory[l] = new double[testedMaterial.Stresses.Length];
                testedMaterial.Stresses.CopyTo(stressHistory[l], 0);
                constitutiveMatrixHistory[l] = new double[testedMaterial.ConstitutiveMatrix.Columns, testedMaterial.ConstitutiveMatrix.Rows];

                for (int m = 0; m < testedMaterial.ConstitutiveMatrix.Columns; m++)
                {
                    for (int n = 0; n < testedMaterial.ConstitutiveMatrix.Rows; n++)
                    { constitutiveMatrixHistory[l][m, n] = testedMaterial.ConstitutiveMatrix[m, n]; }
                }
            }

            return (stressHistory, constitutiveMatrixHistory);
        }

        private static double [,] ConvertStressHistoryTodoubleArray(double[][] stressHistory)
        {
            double[,] stressHistoryArray = new double[stressHistory.Length, stressHistory[0].Length];

            for (int i1 = 0; i1 < stressHistory.Length; i1++)
            {
                for (int i2 = 0; i2 < stressHistory[0].Length; i2++)
                {
                    stressHistoryArray[i1, i2] = stressHistory[i1][i2];
                }
            }

            return stressHistoryArray;
        }

        private static double [][] CreateStrainHistoryInterpolation(double[] targetStrain, double [] initialStrain, int cycles)
        {
            double[][] strainHistory = new double[cycles][];            
            double[] strainIncr = new double[targetStrain.Length]; for (int l = 0; l < targetStrain.Length; l++) { strainIncr[l] = (targetStrain[l] - initialStrain[l]) / ((double)cycles); }
            for (int l = 0; l < cycles; l++)
            {
                strainHistory[l] = new double[targetStrain.Length];
                for (int i1 = 0; i1 < targetStrain.Length; i1++) { strainHistory[l][i1] = initialStrain[i1] + (strainIncr[i1] * ((double)l + 1)); }
            }

            return strainHistory;
        }

        private static double [][] CombineStrainHistoriesIntoOne(double[][] firstStrainHistory, double[][] secondStrainHistory)
        {
            double[][] strainHistory = new double[firstStrainHistory.Length + secondStrainHistory.Length][];

            for (int l = 0; l < firstStrainHistory.Length; l++)
            {
                strainHistory[l] = new double[firstStrainHistory[0].Length];
                for (int i1 = 0; i1 < firstStrainHistory[0].Length; i1++) { strainHistory[l][i1] = firstStrainHistory[l][i1]; }
            }

            for (int l = 0; l < secondStrainHistory.Length; l++)
            {
                strainHistory[l + firstStrainHistory.Length]= new double[firstStrainHistory[0].Length];
                for (int i1 = 0; i1 < secondStrainHistory[0].Length; i1++) { strainHistory[l+ firstStrainHistory.Length][i1] = secondStrainHistory[l][i1]; }
            }

            return strainHistory;
        }
    }
}
