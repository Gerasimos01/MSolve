﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public class LumpedPreconditioner : IFetiPreconditioner
    {
        private readonly Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryBoundary;
        private readonly int[] subdomainIDs;

        private LumpedPreconditioner(int[] subdomainIDs, Dictionary<int, Matrix> stiffnessesBoundaryBoundary,
            Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices)
        {
            this.subdomainIDs = subdomainIDs;
            this.preconditioningBoundarySignedBooleanMatrices = preconditioningBoundarySignedBooleanMatrices;
            this.stiffnessesBoundaryBoundary = stiffnessesBoundaryBoundary;
        }

        //TODO: This can be moved to a base class. Only the S matrix is different for these preconditioners. 
        //      Other ones might be different though.
        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            // inv(F) * y = Bpb * Kbb * Bpb^T * y
            lhs.Clear(); //TODO: this should be avoided
            foreach (int id in subdomainIDs)
            {
                Matrix Bpb = preconditioningBoundarySignedBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Vector subdomainContribution = Bpb.Multiply(Kbb.Multiply(Bpb.Multiply(rhs, true)));
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs)
        {
            // inv(F) * y: Bpb * Kbb * Bpb^T * Y
            lhs.Clear(); //TODO: this should be avoided
            foreach (int id in subdomainIDs)
            {
                Matrix Bpb = preconditioningBoundarySignedBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Matrix subdomainContribution = Bpb.MultiplyRight(Kbb.MultiplyRight(Bpb.MultiplyRight(rhs, true)));
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public class Factory : FetiPreconditionerFactoryBase
        {
            public override IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution,
                IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, IMatrixView> stiffnessMatrices)
            {
                int[] subdomainIDs = dofSeparator.BoundaryDofIndices.Keys.ToArray();
                Dictionary<int, Matrix> boundaryBooleans = CalcBoundaryPreconditioningBooleanMatrices(stiffnessDistribution, 
                    dofSeparator, lagrangeEnumerator);
                Dictionary<int, Matrix> stiffnessesBoundaryBoundary = 
                    ExtractStiffnessesBoundaryBoundary(dofSeparator, stiffnessMatrices);
                return new LumpedPreconditioner(subdomainIDs, stiffnessesBoundaryBoundary, boundaryBooleans);
            }
        }
    }
}