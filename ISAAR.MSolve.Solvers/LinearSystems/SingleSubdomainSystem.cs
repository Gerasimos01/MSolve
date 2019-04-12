﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// Implementation of <see cref="ILinearSystem"/> that can be used with solvers withour domain decomposition.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    /// <typeparam name="TMatrix"></typeparam>
    public class SingleSubdomainSystem<TMatrix> : LinearSystemBase<TMatrix, Vector>
        where TMatrix : class, IMatrix
    {
        internal SingleSubdomainSystem(ISubdomain subdomain) : base(subdomain) { }
        internal override Vector CreateZeroVector() => Vector.CreateZero(this.Size);
    }
}