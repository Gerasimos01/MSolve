﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Factorizations
{
    public static class TriangulationExtensions
    {
        /// <summary>
        /// Solves the linear system A * x = b, where A is the original matrix (before the factorization), 
        /// b = <paramref name="rhsVector"/> and x is the solution vector, which will be returned.
        /// </summary>
        /// <param name="rhsVector">
        /// The right hand side vector. Its <see cref="IIndexable1D.Length"/> must be equal to 
        /// <see cref="Matrices.IIndexable2D.NumRows"/> of the original matrix A.
        /// </param>
        /// Thrown if the length of <paramref name="rhsVector"/> is different than <see cref="Matrices.IIndexable2D.NumRows"/> 
        /// of the original matrix A.
        /// </exception>
        public static Vector SolveLinearSystem(this ITriangulation triangulation, Vector rhsVector)
        {
            var solution = Vector.CreateZero(rhsVector.Length);
            triangulation.SolveLinearSystem(rhsVector, solution);
            return solution;
        }
    }
}
