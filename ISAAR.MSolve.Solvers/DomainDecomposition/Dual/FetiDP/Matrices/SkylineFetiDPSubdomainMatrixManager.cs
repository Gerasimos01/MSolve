﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.SchurComplements;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    /// <summary>
    /// Dense format for Kbb, Kcc, KccStar, skyline for Kff, Krr, Kii and CSC for Kib, Krc.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineFetiDPSubdomainMatrixManager : IFetiDPSubdomainMatrixManager
    {
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly SingleSubdomainSystem<SkylineMatrix> linearSystem;

        private DiagonalMatrix inverseKiiDiagonal;
        private LdlSkyline inverseKii;
        private LdlSkyline inverseKrr;
        private Matrix Kbb;
        private CscMatrix Kib;
        private Matrix Kcc;
        private Matrix KccStar;
        private CscMatrix Krc;
        private SkylineMatrix Krr;

        public SkylineFetiDPSubdomainMatrixManager(ISubdomain subdomain)
        {
            this.linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
        }

        public ISingleSubdomainLinearSystem LinearSystem => linearSystem;

        public Matrix SchurComplementOfRemainderDofs => KccStar;

        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalMatrix(dofOrdering, elements, matrixProvider);


        public (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalSubmatrices(freeDofOrdering, constrainedDofOrdering, elements, matrixProvider);

        public void Clear()
        {
            inverseKii = null;
            inverseKiiDiagonal = null;
            inverseKrr = null;
            Kbb = null;
            Kib = null;
            Kcc = null;
            Krc = null;
            Krr = null;
            //linearSystem.Matrix = null; // DO NOT DO THAT!!! The analyzer manages that.
        }

        public void CalcSchurComplementOfRemainderDofs() //TODO: This should be done in a dedicated class
        {
            // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s]
            KccStar = SchurComplementCsc.CalcSchurComplementFull(Kcc, Krc, inverseKrr);
        }

        public void ExtractAndInvertKii(int[] internalDofs)
        {
            try
            {
                SkylineMatrix Kii = Krr.GetSubmatrixSymmetricSkyline(internalDofs);
                inverseKii = Kii.FactorLdl(true);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractAndInvertKiiDiagonal(int[] internalDofs)
        {
            try
            {
                var diagonal = new double[internalDofs.Length];
                for (int i = 0; i < diagonal.Length; ++i)
                {
                    int idx = internalDofs[i];
                    diagonal[i] = 1.0 / Krr[idx, idx];
                    //diagonal[i] = Krr[idx, idx];
                }
                inverseKiiDiagonal = DiagonalMatrix.CreateFromArray(diagonal, false);
                //inverseKiiDiagonal.Invert();
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractKbb(int[] boundaryDofs)
        {
            try
            {
                Kbb = Krr.GetSubmatrixSymmetricFull(boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractKbiKib(int[] boundaryDofs, int[] internalDofs)
        {
            try
            {
                Kib = Krr.GetSubmatrixCsc(internalDofs, boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractKcc(int[] cornerDofs)
        {
            Kcc = linearSystem.Matrix.GetSubmatrixFull(cornerDofs, cornerDofs);
        }

        public void ExtractKcrKrc(int[] cornerDofs, int[] remainderDofs)
        {
            Krc = linearSystem.Matrix.GetSubmatrixCsc(remainderDofs, cornerDofs);
        }

        public void ExtractKrr(int[] remainderDofs)
        {
            Krr = linearSystem.Matrix.GetSubmatrixSymmetricSkyline(remainderDofs);
        }

        public void HandleDofOrderingWillBeModified() => assembler.HandleDofOrderingWillBeModified();

        public void InvertKrr(bool inPlace)
        {
            inverseKrr = Krr.FactorLdl(inPlace);
        }

        public Vector MultiplyInverseKiiDiagonalTimes(Vector vector)
        {
            if (inverseKiiDiagonal == null)
            {
                throw new InvalidOperationException("The inverse of the diagonal of the internal-internal stiffness submatrix"
                    + " 'inv(diag(Kii))' of this subdomain must be calculated first.");
            }
            return inverseKiiDiagonal * vector;
        }

        public Matrix MultiplyInverseKiiDiagonalTimes(Matrix matrix)
        {
            if (inverseKiiDiagonal == null)
            {
                throw new InvalidOperationException("The inverse of the diagonal of the internal-internal stiffness submatrix"
                    + " 'inv(diag(Kii))' of this subdomain must be calculated first.");
            }
            return inverseKiiDiagonal * matrix;
        }

        public Vector MultiplyInverseKiiTimes(Vector vector)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal_remainder - internal_remainder stiffness"
                    + " submatrix 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystem(vector);
        }

        public Matrix MultiplyInverseKiiTimes(Matrix matrix)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal_remainder - internal_remainder stiffness"
                    + " submatrix 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystems(matrix);
        }

        public Vector MultiplyInverseKrrTimes(Vector vector)
        {
            if (inverseKrr == null) throw new InvalidOperationException("The inverse of the remainder-remainder stiffness "
                + " submatrix 'Krr' of this subdomain must be calculated first.");
            return inverseKrr.SolveLinearSystem(vector);
        }

        public Vector MultiplyKbbTimes(Vector vector)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException("The boundary_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * vector;
        }

        public Matrix MultiplyKbbTimes(Matrix matrix)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException("The boundary_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * matrix;
        }

        public Vector MultiplyKbiTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector, true);
        }

        public Matrix MultiplyKbiTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix, true);
        }

        public Vector MultiplyKccTimes(Vector vector)
        {
            if (Kcc == null)
            {
                throw new InvalidOperationException("The corner-corner stiffness submatrix"
                    + " 'Kcc' of this subdomain must be calculated first.");
            }
            return Kcc * vector;
        }

        public Vector MultiplyKcrTimes(Vector vector)
        {
            if (Krc == null)
            {
                throw new InvalidOperationException("The remainder-corner stiffness submatrix"
                    + " 'Kcr' of this subdomain must be calculated first.");
            }
            return Krc.Multiply(vector, true);
        }

        public Vector MultiplyKibTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector);
        }

        public Matrix MultiplyKibTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix);
        }

        public Vector MultiplyKrcTimes(Vector vector)
        {
            if (Krc == null)
            {
                throw new InvalidOperationException("The remainder-corner stiffness submatrix"
                    + " 'Krc' of this subdomain must be calculated first.");
            }
            return Krc.Multiply(vector);
        }

        public class Factory : IFetiDPSubdomainMatrixManagerFactory
        {
            public IFetiDPSubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain)
                => new SkylineFetiDPSubdomainMatrixManager(subdomain);
        }
    }
}