using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.IGA.Tests;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ISAAR.MSolve.SamplesConsole
{
    class SeparateCodeCheckingClass2
    {
        public static void ExampleParametricStudy2CyclicConferenceEarlySlip()
        {
            //proelefsi ExampleParametricStudy2CyclicConferenceSlip
            //allages gp dosmeno exwterika
            var gpNeo = new grapheneSheetParameters()
            {
                //parametroi cohesive epifaneias
                T_o_3 = 0.12,//0.20,                   //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.15,//0.25,                   //0.5, // nm
                D_f_3 = 4,                      // nm
                T_o_1 = 0.12,//0.20,                   //0.05,// Gpa
                D_o_1 = 0.15,//0.25,                   //0.5, // nm
                D_f_1 = 4,                      // nm
                n_curve = 1.4
            };

            VectorExtensions.AssignTotalAffinityCount();
            //arxiko paradeigma finite strains:
            //IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVEstif();           
            //IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //neo paradeigma:
            IdegenerateRVEbuilder RveBuilder1 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostDataGP(9, gpNeo);
            var Vec1 = new Vector(new double[3] { 1, 0, 0 });
            var Vec2 = new Vector(new double[3] { 0, 1, 0 });
            var material4 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimu(RveBuilder1, false) { TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } }; ;


            int cycles = 5;
            double[] totalStrain = new double[3] { 0.07, 0, 0 };
            double[] initialstrain = new double[3] { 0, 0, 0 };
            double[][] strainHistory1 = CreateStrainHistoryInterpolation(totalStrain, initialstrain, cycles);
            double[][] strainHistory2 = CreateStrainHistoryInterpolation(initialstrain, totalStrain, cycles);
            double[][] strainHistory = CombineStrainHistoriesIntoOne(strainHistory1, strainHistory2);



            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(strainHistory, material4);

            double[,] strainHistoryArray = ConvertStressHistoryTodoubleArray(strainHistory);
            double[,] stressHistoryArray = ConvertStressHistoryTodoubleArray(stressHistory);

            PrintUtilities.WriteToFile(strainHistoryArray, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\strainHistoryShell.txt");
            PrintUtilities.WriteToFile(stressHistoryArray, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressHistoryShell.txt");

        }

        public static void ExampleParametricStudy2CyclicConferenceEarlySlipGrShNum(int cycles, double[] totalStrain, int graphene_sheets_number)
        {
            //proelefsi ExampleParametricStudy2CyclicConferenceSlip
            //allages gp dosmeno exwterika
            var gpNeo = new grapheneSheetParameters()
            {
                //parametroi cohesive epifaneias
                T_o_3 = 0.12,//0.20,                   //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.15,//0.25,                   //0.5, // nm
                D_f_3 = 4,                      // nm
                T_o_1 = 0.12,//0.20,                   //0.05,// Gpa
                D_o_1 = 0.15,//0.25,                   //0.5, // nm
                D_f_1 = 4,                      // nm
                n_curve = 1.4
            };

            VectorExtensions.AssignTotalAffinityCount();
            //arxiko paradeigma finite strains:
            //IRVEbuilder homogeneousRveBuilder1 = new GrapheneReinforcedRVEBuilderExample5GrSh1RVEstif();           
            //IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3Develop(homogeneousRveBuilder1);
            //neo paradeigma:                                                                                                                //EDW DIORTHOSI TOU RVE_id
            IdegenerateRVEbuilder RveBuilder1 = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostTestPostDataGPgrSh(7, gpNeo, graphene_sheets_number);
            var Vec1 = new Vector(new double[3] { 1, 0, 0 });
            var Vec2 = new Vector(new double[3] { 0, 1, 0 });
            var material4 = new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformationSimu(RveBuilder1, false) { TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] }, TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] } }; ;


            //int cycles = 5;
            //double[] totalStrain = new double[3] { 0.07, 0, 0 };
            double[] initialstrain = new double[3] { 0, 0, 0 };
            double[][] strainHistory1 = CreateStrainHistoryInterpolation(totalStrain, initialstrain, cycles);
            double[][] strainHistory2 = CreateStrainHistoryInterpolation(initialstrain, totalStrain, cycles);
            double[][] strainHistory = CombineStrainHistoriesIntoOne(strainHistory1, strainHistory2);



            (double[][] stressHistory, double[][,] constitutiveMatrixHistory) = StressStrainHistory(strainHistory, material4);

            double[,] strainHistoryArray = ConvertStressHistoryTodoubleArray(strainHistory);
            double[,] stressHistoryArray = ConvertStressHistoryTodoubleArray(stressHistory);

            string strainString = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\strainHistoryShell" + graphene_sheets_number.ToString() + @"Graphenes.txt";
            string stressString = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressHistoryShell" + graphene_sheets_number.ToString() + @"Graphenes.txt";
            PrintUtilities.WriteToFile(strainHistoryArray, strainString);
            PrintUtilities.WriteToFile(stressHistoryArray, stressString);

        }


















        private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, ICohesiveZoneMaterial3D testedMaterial)
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

        private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IShellMaterial testedMaterial)
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

        private static double[,] ConvertStressHistoryTodoubleArray(double[][] stressHistory)
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

        private static double[][] CreateStrainHistoryInterpolation(double[] targetStrain, double[] initialStrain, int cycles)
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

        private static double[][] CombineStrainHistoriesIntoOne(double[][] firstStrainHistory, double[][] secondStrainHistory)
        {
            double[][] strainHistory = new double[firstStrainHistory.Length + secondStrainHistory.Length][];

            for (int l = 0; l < firstStrainHistory.Length; l++)
            {
                strainHistory[l] = new double[firstStrainHistory[0].Length];
                for (int i1 = 0; i1 < firstStrainHistory[0].Length; i1++) { strainHistory[l][i1] = firstStrainHistory[l][i1]; }
            }

            for (int l = 0; l < secondStrainHistory.Length; l++)
            {
                strainHistory[l + firstStrainHistory.Length] = new double[firstStrainHistory[0].Length];
                for (int i1 = 0; i1 < secondStrainHistory[0].Length; i1++) { strainHistory[l + firstStrainHistory.Length][i1] = secondStrainHistory[l][i1]; }
            }

            return strainHistory;
        }
    }
}
