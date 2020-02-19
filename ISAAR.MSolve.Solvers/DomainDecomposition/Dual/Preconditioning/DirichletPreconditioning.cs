﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public class DirichletPreconditioning : IFetiPreconditioningOperations
    {
        public bool ReorderInternalDofsForFactorization => true;

        public void PrepareSubdomainSubmatrices(IFetiSubdomainMatrixManager matrixManager)
        {
            ISubdomain subdomain = matrixManager.LinearSystem.Subdomain;
            if (subdomain.StiffnessModified)
            {
                Debug.WriteLine($"{this.GetType().Name}:"
                    + $" Extracting boundary/internal submatrices of subdomain {subdomain.ID} for preconditioning");
                matrixManager.ExtractKbb();
                matrixManager.ExtractKbiKib();
                matrixManager.CalcInverseKii(false);
            }
        }

        public Matrix PreconditionSubdomainMatrix(Matrix rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb)
        {
            // inv(F) * Y = Bpb * S * Bpb^T * Y
            // S = Kbb - Kbi * inv(Kii) * Kib
            Matrix BY = Bpb.MultiplyRight(rhs, true);
            Matrix temp = matrixManager.MultiplyKibTimes(BY);
            temp = matrixManager.MultiplyInverseKiiTimes(temp, false);
            temp = matrixManager.MultiplyKbiTimes(temp);
            temp = matrixManager.MultiplyKbbTimes(BY) - temp;
            Matrix subdomainContribution = Bpb.MultiplyRight(temp);
            return subdomainContribution;
        }

        public Vector PreconditionSubdomainVector(Vector rhs, IFetiSubdomainMatrixManager matrixManager, IMappingMatrix Bpb)
        {
            // inv(F) * y = Bpb * S * Bpb^T * y
            // S = Kbb - Kbi * inv(Kii) * Kib
            Vector By = Bpb.Multiply(rhs, true);
            Vector temp = matrixManager.MultiplyKibTimes(By);
            temp = matrixManager.MultiplyInverseKiiTimes(temp, false);
            temp = matrixManager.MultiplyKbiTimes(temp);
            temp = matrixManager.MultiplyKbbTimes(By) - temp;
            Vector subdomainContribution = Bpb.Multiply(temp);
            return subdomainContribution;
        }
    }
}