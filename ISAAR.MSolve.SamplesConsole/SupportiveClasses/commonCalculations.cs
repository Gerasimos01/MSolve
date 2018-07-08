using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    public class commonCalculations
    {
        public static void cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }
        public static double dot_product(double[] vec1, double[] vec2)
        {
            double dot_product = 0;
            for (int i1 = 0; i1 < vec1.GetLength(0); i1++)
            {
                dot_product += vec1[i1] * vec2[i1];
            }
            return dot_product;
        }
        public static double[] MatVecMult(double[,] matrix, double[] vec)
        {
            double[] product = new double[matrix.GetLength(0)];
            for (int q1 = 0; q1 < matrix.GetLength(0); q1++)
            {
                for (int q2 = 0; q2 < matrix.GetLength(1); q2++)
                {
                    product[q1] += matrix[q1, q2] * vec[q2];
                }
            }
            return product;
        }
    }
}
