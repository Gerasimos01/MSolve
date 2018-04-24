﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

/// <summary>
/// The matrix that will be "inverted" is a unified DOK matrix, where enriched dofs are numbered after all 
/// standard dofs. 
/// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
/// needed. Also the global constrained submatrix must be sparse.
/// </summary>
namespace ISAAR.MSolve.XFEM.Assemblers
{
    static class SingleGlobalDOKAssembler
    {
        public static (DOKSymmetricColMajor Kuu, Matrix Kuc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int constrainedDofsCount = dofEnumerator.ConstrainedDofsCount;
            int allFreeDofsCount = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            var Kuu = DOKSymmetricColMajor.CreateEmpty(allFreeDofsCount);

            // TODO: this should be in a sparse format. Only used for SpMV and perhaps transpose SpMV.
            // Row = standard free dofs + enriched dofs. Columns = standard constrained dofs. 
            Matrix Kuc = Matrix.CreateZero(dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount,
                dofEnumerator.ConstrainedDofsCount);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                Matrix kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                // 2) The Elements have two versions of BuildStiffness(). 
                // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                Matrix kse = kes.Transpose(); 

                // Element to global dofs mappings
                // TODO: perhaps that could be done during the assembly to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element, 
                    out IReadOnlyDictionary<int, int> mapFree, 
                    out IReadOnlyDictionary<int, int> mapConstrained);
                IReadOnlyDictionary<int, int> mapEnriched = dofEnumerator.MatchElementToGlobalEnrichedDofsOf(element);

                // Add the element contributions to the global matrices
                Kuu.AddSubmatrixSymmetric(kss, mapFree);
                Kuu.AddSubmatrix(kse, mapFree, mapEnriched);
                Kuu.AddSubmatrixSymmetric(kee, mapEnriched);

                AddElementToGlobalMatrix(Kuc, kss, mapFree, mapConstrained);
                AddElementToGlobalMatrix(Kuc, kes, mapEnriched, mapConstrained);
            }

            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofEnumerator);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            return (Kuu, Kuc);
        }

        /// <summary>
        /// This method is for the non symmetric part of the matrix (free - constrained).
        /// TODO: The global matrix should be sparse. Probably in CSR/DOK format.
        /// </summary>
        /// <param name="globalMatrix"></param>
        /// <param name="elementMatrix"></param>
        /// <param name="elementRowsToGlobalRows"></param>
        /// <param name="elementColsToGlobalCols"></param>
        private static void AddElementToGlobalMatrix(Matrix globalMatrix, Matrix elementMatrix,
            IReadOnlyDictionary<int, int> elementRowsToGlobalRows, IReadOnlyDictionary<int, int> elementColsToGlobalCols)
        {
            foreach (var rowPair in elementRowsToGlobalRows)
            {
                int elementRow = rowPair.Key;
                int globalRow = rowPair.Value;
                foreach (var colPair in elementColsToGlobalCols)
                {
                    int elementCol = colPair.Key;
                    int globalCol = colPair.Value;

                    globalMatrix[globalRow, globalCol] += elementMatrix[elementRow, elementCol];
                }
            }
        }

        private static void CheckMatrix(IIndexable2D expected, IIndexable2D computed)
        {
            bool isCorrect = true;
            ValueComparer comparer = new ValueComparer(1e-6);
            if ((computed.NumRows != expected.NumRows) || (computed.NumColumns != expected.NumColumns))
            {
                Console.WriteLine("Invalid dimensions");
            }
            for (int i = 0; i < computed.NumRows; ++i)
            {
                for (int j = 0; j < computed.NumColumns; ++j)
                {
                    if (!comparer.AreEqual(computed[i, j], expected[i, j]))
                    {
                        Console.WriteLine($"Computed[{i}, {j}] = {computed[i, j]}   -   Expected[{i}, {j}] = {expected[i, j]}");
                        isCorrect = false;
                    }
                }
            }
            if (isCorrect) Console.WriteLine("Correct");
        }
    }
}
