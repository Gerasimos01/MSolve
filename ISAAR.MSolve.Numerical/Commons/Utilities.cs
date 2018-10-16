using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.Commons
{
    public static class Utilities
    {
        public static bool AreEqual(Matrix2D matrix1, Matrix2D matrix2, double tolerance)
        {
            if ((matrix1.Rows != matrix2.Rows) || (matrix1.Columns != matrix2.Columns)) return false;

            var comparer = new ValueComparer(tolerance);
            for (int i = 0; i < matrix1.Rows; ++i)
            {
                for (int j = 0; j < matrix1.Columns; ++j)
                {
                    if (!comparer.AreEqual(matrix1[i, j], matrix2[i, j])) return false;
                }
            }
            return true;
        }
    }
}
