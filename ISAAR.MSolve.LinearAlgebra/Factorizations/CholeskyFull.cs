﻿using System;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.MKL;
using ISAAR.MSolve.LinearAlgebra.Vectors;


// TODO: check if the last minor is non-negative, during factorization. Is it possible that it isn't. Does it affect system 
// solution or inversion?
namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    /// <summary>
    /// Cholesky factorization of a symmetric positive definite matrix, stored in full column major format. Only the upper
    /// triangle part of the matrix is factorized. The subdiagonal part of the original matrix is still stored in the array but 
    /// it is ignored. Uses Intel MKL.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskyFull: ITriangulation
    {
        private readonly double[] data;

        private CholeskyFull(int order, double[] data)
        {
            this.Order = order;
            this.data = data;
        }

        /// <summary>
        /// If true, the internal data of this object are overwritten and used by another object. No property or method of
        /// this object must be called as it would throw exceptions or lead to data corruption. If false, this object can be 
        /// used normally.
        /// </summary>
        public bool IsOverwritten { get; private set; }

        /// <summary>
        /// The number of rows/columns of the original square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Calculates the cholesky factorization of a symmetric positive definite matrix, such that A = transpose(U) * U. 
        /// </summary>
        /// <param name="order">The number of rows/columns of the square matrix.</param>
        /// <param name="matrix">The entries of the original matrix in full column major layout.</param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not symmetric positive definite.</exception>
        public static CholeskyFull Factorize(int order, double[] matrix)
        {
            // Call MKL
            int info = MklUtilities.DefaultInfo;
            Lapack.Dpotrf("U", ref order, ref matrix[0], ref order, ref info);

            // Check MKL execution
            if (info == 0) return new CholeskyFull(order, matrix);
            else if (info > 0)
            {
                string msg = $"The leading minor of order {info-1} (and therefore the matrix itself) is not"
                + " positive-definite, and the factorization could not be completed.";
                throw new IndefiniteMatrixException(msg);
            }
            else throw MklUtilities.ProcessNegativeInfo(info); // info < 0
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        /// <remarks>
        /// A = U^T*U => det(A) = det(U^T)* det(U) => det(A) = (det(U))^2, where det(U) = U[0,0] * U[1,1] * ... * U[n,n]
        /// </remarks>
        public double CalcDeterminant()
        {
            CheckOverwritten();
            double det = 1.0;
            for (int i = 0; i < Order; ++i)
            {
                det *= data[i * Order + i];
            }
            return det * det;
        }

        /// <summary>
        /// Explicitly creates the upper triangular matrix U that resulted from the Cholesky factorization: A = transpose(U) * U,
        /// where A and U are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public Matrix GetFactorU() //TODO: Return TriangleUpper instead.
        {
            CheckOverwritten();
            double[] u = Conversions.FullColMajorToFullUpperColMajor(data, false);
            return Matrix.CreateFromArray(u, Order, Order, false);
        }

        /// <summary>
        /// Calculates the inverse of the original matrix and returns it in a new <see cref="Matrix"/> instance. 
        /// WARNING: If <paramref name="inPlace"/> is set to true, this object must not be used again, otherwise a 
        /// <see cref="InvalidOperationException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal factorization data before inversion. True, to overwrite it with
        ///     the inverse matrix, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        public Matrix Invert(bool inPlace)
        {
            CheckOverwritten();

            // Call MKL
            int info = MklUtilities.DefaultInfo;
            double[] inverse;
            if (inPlace)
            {
                inverse = data;
                IsOverwritten = true;
            }
            else
            {
                inverse = new double[data.Length];
                Array.Copy(data, inverse, data.Length);
            }
            info = LAPACKE.Dpotri(LAPACKE.LAPACK_COL_MAJOR, LAPACKE.LAPACK_UPPER, Order, inverse, Order);
            Conversions.CopyUpperToLowerColMajor(inverse, Order); //So far the lower triangle was the same as the original matrix

            // Check MKL execution
            if (info == 0) return Matrix.CreateFromArray(inverse, Order, Order, false);
            else if (info > 0) // this should not have happened
            {
                throw new IndefiniteMatrixException($"The leading minor of order {info - 1} (and therefore the matrix itself)"
                + "is not positive-definite, and the factorization could not be completed.");
            }
            else throw MklUtilities.ProcessNegativeInfo(info); // info < 0
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector)"/>.
        /// </summary>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid arguments.</exception>
        public Vector SolveLinearSystem(Vector rhs)
        {
            CheckOverwritten();
            Preconditions.CheckSystemSolutionDimensions(this.Order, this.Order, rhs.Length);

            // Back & forward substitution using MKL
            int n = Order;
            double[] b = rhs.CopyToArray();
            int info = MklUtilities.DefaultInfo;
            int nRhs = 1; // rhs is a n x nRhs matrix, stored in b
            int ldb = n; // column major ordering: leading dimension of b is n 
            Lapack.Dpotrs("U", ref n, ref nRhs, ref data[0], ref n, ref b[0], ref ldb, ref info);

            // Check MKL execution
            if (info == 0) return Vector.CreateFromArray(b, false);
            else throw MklUtilities.ProcessNegativeInfo(info); // info < 0. This function does not return info > 0
        }

        private void CheckOverwritten()
        {
            if (IsOverwritten) throw new InvalidOperationException(
                "The internal buffer of this factorization has been overwritten and thus cannot be used anymore.");
        }
    }
}