using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Readers;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public static class ElementStiffnesses
    {

        public static bool saveStiffnessMatrixState { get; set; } = false;
        public static bool saveForcesState1 { get; set; } = false;
        public static bool saveForcesState2s { get; set; } = false;
        public static bool saveForcesState0 { get; set; } = false;
        public static bool saveOriginalState { get; set; } = false;

        public static bool saveVariationStates { get; set; } = false;

        public static bool performCalculations { get; set; } = false;


        public static Dictionary<int, double[,]> OriginalLinearStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalMembraneNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalBendingNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalMembraneLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalBendingLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrixpublic static Dictionary<int, double[,]> OriginalNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix

        public static Dictionary<int, double[,]> Stiffnesses1FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> Stiffnesses2FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> MembraneStiffnesses1FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> MembraneStiffnesses2FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> BendingStiffnesses1FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> BendingStiffnesses2FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> variatedForcesStep1 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedForcesStep2 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[]> originalForces = new Dictionary<int, double[]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedMembraneForcesStep1 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedMembraneForcesStep2 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[]> originalMemebraneForces = new Dictionary<int, double[]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedBendingForcesStep1 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedBendingForcesStep2 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[]> originalBendingForces = new Dictionary<int, double[]>(); // element, Forces


        public static Dictionary<int, Dictionary<int, double[,]>> derivativesCalculated = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId,  variableId, variableDerivatives

        public static Dictionary<int, Dictionary<int, double[,]>> derivativesExpected = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId,  variableId, variableDerivatives

        public static Dictionary<int, Dictionary<int, double[]>> variableOriginalValues = new Dictionary<int, Dictionary<int, double[]>>(); // elementId, variableId, variableOriginalValues

        /// <summary>
        /// elementId, variableId, variableVariatedValues
        /// </summary>
        public static Dictionary<int, Dictionary<int, double[,]>> variableVariatedValues = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId, variableId, variableVariatedValues

        public static int variationDofOrdef;

        public static int ForceVariationDofOrdef;

        public static double variationSize = 0.0000001;

        public static int[] variableIds = new int[] { 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,30,31,32 }; //a1,1,ur   //a2,2,ur  // a1,2,ur  // a1,r
        public static int[] variableSizes = new int[] { 3, 3, 3, 3, 3, 3, 3,3,3,3,1,3,3,1,1,3,3,9,3,3,1,1,3,3,9,4,4,9,3,3  };

        public static int elementId { get; set; }

        public static int gpNumber { get; set; }
        public static int gpNumberToCheck { get; set; } = 0;

        public static void ProcessElement(double[] localSolution, double[] localdSolution, Element element)
        {
            ElementStiffnesses.elementId = element.ID;

            #region initialize variables
            derivativesCalculated[ElementStiffnesses.elementId] = new Dictionary<int, double[,]>();
            derivativesExpected[ElementStiffnesses.elementId] = new Dictionary<int, double[,]>();
            variableOriginalValues[ElementStiffnesses.elementId] = new Dictionary<int, double[]>();
            variableVariatedValues[ElementStiffnesses.elementId] = new Dictionary<int, double[,]>();
            for (int i1 = 0; i1 < variableIds.Length; i1++)
            {
                derivativesCalculated[ElementStiffnesses.elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
                derivativesExpected[ElementStiffnesses.elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
                variableOriginalValues[ElementStiffnesses.elementId][variableIds[i1]] = new double[variableSizes[i1]];
                variableVariatedValues[ElementStiffnesses.elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
            }
            variatedForcesStep1[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];
            variatedForcesStep2[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];

            originalMemebraneForces[ElementStiffnesses.elementId] = new double[localSolution.Length];
            originalBendingForces[ElementStiffnesses.elementId] = new double[localSolution.Length];

            variatedMembraneForcesStep1[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];
            variatedBendingForcesStep1[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];

            variatedMembraneForcesStep2[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];
            variatedBendingForcesStep2[ElementStiffnesses.elementId] = new double[localSolution.Length, localSolution.Length];
            #endregion


            #region original state
            performCalculations=true;
            saveOriginalState = true;
            element.ElementType.CalculateStresses(element, localSolution, localdSolution);
            saveForcesState0 = true;
            originalForces[ElementStiffnesses.elementId] = element.ElementType.CalculateForces(element, localSolution, localdSolution);
            saveForcesState0 = false;
            saveOriginalState = false;
            saveStiffnessMatrixState = true;
            OriginalStiffnesses[ElementStiffnesses.elementId] = element.ElementType.StiffnessMatrix(element).CopytoArray2D();
            saveStiffnessMatrixState = false;

            #endregion

            #region variation states
            int nDofs = localSolution.Length;
            Stiffnesses1FromVariations[element.ID] = new double[nDofs, nDofs];
            Stiffnesses2FromVariations[element.ID] = new double[nDofs, nDofs];
            MembraneStiffnesses1FromVariations[element.ID] = new double[nDofs, nDofs];
            MembraneStiffnesses2FromVariations[element.ID] = new double[nDofs, nDofs];
            BendingStiffnesses1FromVariations[element.ID] = new double[nDofs, nDofs];
            BendingStiffnesses2FromVariations[element.ID] = new double[nDofs, nDofs];

            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                variationDofOrdef = i1;
                localSolution[i1] += variationSize;

                saveForcesState1 = true; ForceVariationDofOrdef = i1;
                var forces1 = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                saveForcesState1 = false;
                for (int i2 = 0; i2 < localSolution.Length; i2++) variatedForcesStep1[ElementStiffnesses.elementId][i2, i1] = forces1[i2];

                localSolution[i1] += -variationSize;

            }
            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                variationDofOrdef = i1;
                localSolution[i1] += variationSize;
                saveVariationStates = true;
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                element.ElementType.StiffnessMatrix(element);
                saveVariationStates = false;
                saveForcesState2s = true; ForceVariationDofOrdef = i1;
                var forces2 = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                saveForcesState2s = false;
                for (int i2 = 0; i2 < localSolution.Length; i2++) variatedForcesStep2[ElementStiffnesses.elementId][i2, i1] = forces2[i2];
                localSolution[i1] += -variationSize;

                element.ElementType.CalculateStresses(element, localSolution, localdSolution);

            }
            #endregion

            #region calculations
            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                for (int i3 = 0; i3 < localSolution.Length; i3++)
                {
                    Stiffnesses1FromVariations[element.ID][i3, i1] = (variatedForcesStep1[element.ID][i3, i1] - originalForces[element.ID][i3]) / variationSize;
                    Stiffnesses2FromVariations[element.ID][i3, i1] = (variatedForcesStep2[element.ID][i3, i1] - variatedForcesStep1[element.ID][i3, i1]) / variationSize;
                    MembraneStiffnesses1FromVariations[element.ID][i3, i1] = (variatedMembraneForcesStep1[element.ID][i3, i1] - originalMemebraneForces[element.ID][i3]) / variationSize;
                    MembraneStiffnesses2FromVariations[element.ID][i3, i1] = (variatedMembraneForcesStep2[element.ID][i3, i1] - variatedMembraneForcesStep1[element.ID][i3, i1]) / variationSize;
                    BendingStiffnesses1FromVariations[element.ID][i3, i1] = (variatedBendingForcesStep1[element.ID][i3, i1] - originalBendingForces[element.ID][i3]) / variationSize;
                    BendingStiffnesses2FromVariations[element.ID][i3, i1] = (variatedBendingForcesStep2[element.ID][i3, i1] - variatedBendingForcesStep1[element.ID][i3, i1]) / variationSize;

                }
            }

            for (int i1 = 0; i1 < variableIds.Length; i1++)
            {
                for (int i2 = 0; i2 < localSolution.Length; i2++)
                {
                    for (int i3 = 0; i3 < variableSizes[i1]; i3++)
                    {
                        derivativesCalculated[ElementStiffnesses.elementId][variableIds[i1]][i3, i2] = (variableVariatedValues[ElementStiffnesses.elementId][variableIds[i1]][i3, i2] - variableOriginalValues[ElementStiffnesses.elementId][variableIds[i1]][i3]) / variationSize;
                    }

                }

            }
            #endregion
            performCalculations = false;

            ///int elementId = 0;
            #region print matrices
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(OriginalMembraneNLStiffnesses[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\OriginalMembraneNLStiffnesses.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(OriginalBendingNLStiffnesses[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\OriginalBendingNLStiffnesses.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(OriginalMembraneLStiffnesses[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\OriginalMembraneLStiffnesses.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(OriginalBendingLStiffnesses[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\OriginalBendingLStiffnesses.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(MembraneStiffnesses1FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\MembraneStiffnesses1FromVariations.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(MembraneStiffnesses2FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\MembraneStiffnesses2FromVariations.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(BendingStiffnesses1FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\BendingStiffnesses1FromVariations.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(BendingStiffnesses2FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\BendingStiffnesses2FromVariations.txt");

            #region print variables 
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][8], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated8.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][8], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected8.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][9], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated9.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][9], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected9.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][1], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated1.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][1], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected1.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][4], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected4.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][4], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated4.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][5], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected5.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][5], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated5.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][7], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated7.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][7], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected7.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][6], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated6.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][6], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected6.txt");

            // new variables 
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][11], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated11.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][11], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected11.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][12], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated12.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][12], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected12.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][13], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated13.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][13], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected13.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][14], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected14.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][14], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated14.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][15], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected15.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][15], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated15.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][16], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated16.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][16], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected16.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][17], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated17.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][17], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected17.txt");

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][18], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated18.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][18], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected18.txt");

            for (int i1 = 21; i1 < 33; i1++) //mporei na ginei me variableids
            {
                (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesCalculated[elementId][i1], $@"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesCalculated{i1}.txt");
                (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(derivativesExpected[elementId][i1], $@"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\derivativesExpected{i1}.txt");
            }

            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(Stiffnesses2FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\Stiffnesses2FromVariations.txt");
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(Stiffnesses1FromVariations[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\Stiffnesses1FromVariations.txt");
            
            (new ISAAR.MSolve.LinearAlgebra.Output.Array2DWriter()).WriteToFile(OriginalStiffnesses[elementId], @"C:\Users\acivi\Documents\notes_elegxoi_2\develop_nl_iga_shell\MSolve_output\OriginalStiffnesses.txt");
            #endregion


            #endregion

        }

        public static void AccessChecks()
        {
            var breakPoint = "here";
        }

        public static void ProccessVariable(int variableId, double[] variableValue, bool isExpectedValue, int variatedDofFor = 0)
        {
            if (saveOriginalState)
            {
                if (isExpectedValue)
                {
                    for (int i1 = 0; i1 < variableValue.Length; i1++)
                    {
                        derivativesExpected[elementId][variableId][i1, variatedDofFor] = variableValue[i1];

                    }
                }
                else
                {
                    for (int i1 = 0; i1 < variableValue.Length; i1++)
                    {
                        variableOriginalValues[elementId][variableId][i1] = variableValue[i1];
                    }

                }
            }

            if (saveVariationStates)
            {
                for (int i1 = 0; i1 < variableValue.Length; i1++)
                {
                    variableVariatedValues[elementId][variableId][i1, variationDofOrdef] = variableValue[i1];
                }
            }
        }

        public static void ProccessVariable(int variableId, double[,] variableValue, bool isExpectedValue)
        {
            if (saveOriginalState)
            {
                if (isExpectedValue)
                {
                    for (int i1 = 0; i1 < variableValue.GetLength(0); i1++)
                    {
                        for (int i2 = 0; i2 < variableValue.GetLength(1); i2++)
                        {
                            derivativesExpected[elementId][variableId][i1, i2] = variableValue[i1, i2];
                        }
                    }
                }
            }
        }

        internal static void SaveStiffnessMatrixes(double[,] kL_total, double[,] kmembraneNL_total, double[,] kbendingNL_total, Matrix KNLTotal,
             double[,] KmembraneL, double[,] KbendingL, IElement element)
        {
            OriginalLinearStiffnesses[element.ID] = kL_total;
            OriginalMembraneNLStiffnesses[element.ID] = kmembraneNL_total;
            OriginalBendingNLStiffnesses[element.ID] = kbendingNL_total;
            OriginalMembraneLStiffnesses[element.ID] = KmembraneL;
            OriginalBendingLStiffnesses[element.ID] = KbendingL;
            OriginalNLStiffnesses[element.ID] = KNLTotal.CopytoArray2D();
        }

        internal static void SaveNodalForces(double[] elementNodalMembraneForces, double[] elementNodalBendingForces, IElement element)
        {
            if (ElementStiffnesses.saveForcesState0)
            {
                originalMemebraneForces[element.ID]=elementNodalMembraneForces; // element, Forces
                originalBendingForces[element.ID] = elementNodalBendingForces;
            }

            if (ElementStiffnesses.saveForcesState1)
            {
                for (int i = 0; i < elementNodalMembraneForces.Length; i++)
                {
                    variatedMembraneForcesStep1[element.ID][i, ForceVariationDofOrdef] = elementNodalMembraneForces[i];
                    variatedBendingForcesStep1[element.ID][i, ForceVariationDofOrdef] = elementNodalBendingForces[i];

                }
            }

            if (ElementStiffnesses.saveForcesState2s)
            {
                for (int i = 0; i < elementNodalMembraneForces.Length; i++)
                {
                    variatedMembraneForcesStep2[element.ID][i, ForceVariationDofOrdef] = elementNodalMembraneForces[i];
                    variatedBendingForcesStep2[element.ID][i, ForceVariationDofOrdef] = elementNodalBendingForces[i];

                }
            }
        }
    }
}

