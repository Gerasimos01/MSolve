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
        public static bool saveOriginalState { get; set; } = false;

        public static bool saveVariationStates { get; set; } = false;

        
        public static Dictionary<int, double[,]> OriginalLinearStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalMembraneNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalBendingNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalNLStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix
        public static Dictionary<int, double[,]> OriginalStiffnesses = new Dictionary<int, double[,]>(); // element, StiffnesMatrix

        public static Dictionary<int, double[,]> Stiffnesses1FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> Stiffnesses2FromVariations = new Dictionary<int, double[,]>(); //element, StiffnessMatrix

        public static Dictionary<int, double[,]> variatedForcesStep1 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[,]> variatedForcesStep2 = new Dictionary<int, double[,]>(); // element, Forces

        public static Dictionary<int, double[]> originalForces = new Dictionary<int, double[]>(); // element, Forces

        public static Dictionary<int, Dictionary<int, double[,]>> derivativesCalculated = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId,  variableId, variableDerivatives

        public static Dictionary<int, Dictionary<int, double[,]>> derivativesExpected = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId,  variableId, variableDerivatives

        public static Dictionary<int, Dictionary<int, double[]>> variableOriginalValues = new Dictionary<int, Dictionary<int, double[]>>(); // elementId, variableId, variableOriginalValues

        /// <summary>
        /// elementId, variableId, variableVariatedValues
        /// </summary>
        public static Dictionary<int, Dictionary<int, double[,]>> variableVariatedValues = new Dictionary<int, Dictionary<int, double[,]>>(); // elementId, variableId, variableVariatedValues

        public static int variationDofOrdef;

        public static double variationSize = 0.000001;

        public static int[] variableIds = new int[] { 1, 2, 3, 4,5 }; //a1,1,ur   //a2,2,ur  // a1,2,ur  // a1,r
        public static int[] variableSizes = new int[] { 3, 3, 3, 3,3 };

        public static int elementId { get; set; }

        public static int gpNumber { get; set; }
        public static int gpNumberToCheck { get; set; } = 0;

        public static void ProcessElement(double[] localSolution, double[] localdSolution, Element element)
        {
            elementId = element.ID;

            #region initialize variables
            derivativesCalculated[elementId] = new Dictionary<int, double[,]>();
            derivativesExpected[elementId] = new Dictionary<int, double[,]>();
            variableOriginalValues[elementId] = new Dictionary<int, double[]>();
            variableVariatedValues[elementId] = new Dictionary<int, double[,]>();
            for (int i1 = 0; i1 < variableIds.Length; i1++)
            {
                derivativesCalculated[elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
                derivativesExpected[elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
                variableOriginalValues[elementId][variableIds[i1]] = new double[variableSizes[i1]];
                variableVariatedValues[elementId][variableIds[i1]] = new double[variableSizes[i1], localSolution.Length];
            }
            variatedForcesStep1[elementId] = new double[localSolution.Length, localSolution.Length];
            variatedForcesStep2[elementId] = new double[localSolution.Length, localSolution.Length];
            #endregion


            #region original state
            saveOriginalState = true;
            element.ElementType.CalculateStresses(element, localSolution, localdSolution);
            originalForces[elementId] = element.ElementType.CalculateForces(element, localSolution, localdSolution);
            OriginalStiffnesses[elementId] = element.ElementType.StiffnessMatrix(element).CopytoArray2D();
            saveOriginalState = false;
            #endregion

            #region variation states
            int nDofs = localSolution.Length;
            Stiffnesses1FromVariations[element.ID] = new double[nDofs, nDofs];
            Stiffnesses2FromVariations[element.ID] = new double[nDofs, nDofs];
            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                variationDofOrdef = i1;
                localSolution[i1] += variationSize;

                var forces1 = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                for (int i2 = 0; i2 < localSolution.Length; i2++) variatedForcesStep1[elementId][i2, i1] = forces1[i2];

                localSolution[i1] += -variationSize;

            }
            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                variationDofOrdef = i1;
                localSolution[i1] += variationSize;
                saveVariationStates = true;
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                saveVariationStates = false;
                var forces2 = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                for (int i2 = 0; i2 < localSolution.Length; i2++) variatedForcesStep2[elementId][i2, i1] = forces2[i2];
                localSolution[i1] += -variationSize;
                
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);

            }
            #endregion

            #region calculations
            for (int i1 = 0; i1 < localSolution.Length; i1++)
            {
                for (int i3 = 0; i3 < localSolution.Length; i3++)
                {
                    Stiffnesses1FromVariations[element.ID][i3, i1] = (variatedForcesStep1[element.ID][i3, i1] - originalForces[element.ID][i3])/variationSize;
                    Stiffnesses2FromVariations[element.ID][i3, i1] = (variatedForcesStep2[element.ID][i3, i1] - variatedForcesStep1[element.ID][i3, i1]) / variationSize;
                }
            }

            for (int i1 = 0; i1 < variableIds.Length; i1++)
            {
               for (int i2 = 0; i2 < localSolution.Length; i2++)
                {
                    for (int i3 = 0; i3 < variableSizes[i1]; i3++)
                    {
                        derivativesCalculated[elementId][variableIds[i1]][i3, i2] = (variableVariatedValues[elementId][variableIds[i1]][i3, i2] - variableOriginalValues[elementId][variableIds[i1]][i3]) / variationSize;
                    }

                }

            }
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

        internal static void SaveStiffnessMatrixes(double[,] kL_total, double[,] kmembraneNL_total, double[,] kbendingNL_total, Matrix KNLTotal, IElement element)
        {
            OriginalLinearStiffnesses[element.ID] = kL_total;
            OriginalMembraneNLStiffnesses[element.ID] = kmembraneNL_total;
            OriginalBendingNLStiffnesses[element.ID] = kbendingNL_total;
            OriginalNLStiffnesses[element.ID] = KNLTotal.CopytoArray2D();
        }
    }
}

