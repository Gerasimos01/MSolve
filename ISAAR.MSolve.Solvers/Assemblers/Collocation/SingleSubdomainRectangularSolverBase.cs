﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;

namespace ISAAR.MSolve.Solvers.Assemblers.Collocation
{
	public abstract class SingleSubdomainRectangularSolverBase<TMatrix>:ISolver
	where TMatrix: class,IMatrix
	{
        protected readonly IGlobalMatrixRectangularAssembler<TMatrix> assembler;
        protected readonly IStructuralAsymmetricModel model;
        protected readonly IAsymmetricSubdomain subdomain;
        protected readonly SingleSubdomainSystem<TMatrix> linearSystem;
        protected readonly IAsymmetricDofOrderer dofRowOrderer;
        protected readonly IDofOrderer dofColOrderer;

		protected SingleSubdomainRectangularSolverBase(IStructuralAsymmetricModel model, IAsymmetricDofOrderer dofRowOrderer,
			IDofOrderer dofColOrderer, IGlobalMatrixRectangularAssembler<TMatrix> assembler, string name)
		{
			if (model.Subdomains.Count != 1) throw new InvalidSolverException(
				$"{name} can be used if there is only 1 subdomain");
			this.model = model;
			subdomain= model.Subdomains[0];
			linearSystem= new SingleSubdomainSystem<TMatrix>(subdomain);
			LinearSystems= new Dictionary<int, ILinearSystem>(){{subdomain.ID, linearSystem}};
			linearSystem.MatrixObservers.Add(this);

			this.dofRowOrderer = dofRowOrderer;
			this.dofColOrderer = dofColOrderer;
			this.assembler = assembler;
		}
        
		public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }

        public virtual Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            return new Dictionary<int, IMatrix>
            {
                {
                    subdomain.ID,
                    assembler.BuildGlobalMatrix(subdomain.FreeDofRowOrdering, subdomain.FreeDofColOrdering,
                        subdomain.Elements, elementMatrixProvider)
                }
            };
        }
        

        public virtual Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView
            matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            if (subdomain.ConstrainedDofOrdering == null)
            {
                throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs,"
                                                    + " they must have been ordered first.");
            }

            return new Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView
                matrixConstrFree, IMatrixView matrixConstrConstr)>
            {
                {
                    subdomain.ID,
                    assembler.BuildGlobalSubmatrices(subdomain.FreeDofRowOrdering, subdomain.FreeDofColOrdering,
                        subdomain.ConstrainedDofRowOrdering, subdomain.ConstrainedDofColOrdering, subdomain.Elements,
                        elementMatrixProvider)
                }
            };
        }

        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, IDofType, double> globalNodalLoads)
        {
            return new Dictionary<int, SparseVector>();
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            if (otherMatrix.Count!=1) throw new InvalidSolverException("There can only be 1 subdomain when using this solver");
            KeyValuePair<int, IMatrixView> idMatrixPair = otherMatrix.First();
            int id = idMatrixPair.Key;
            Debug.Assert(id == subdomain.ID,
                "The matrix that will be multiplied with the inverse system matrix belongs to a different subdomain.");
            Matrix result = InverseSystemMatrixTimesOtherMatrix(idMatrixPair.Value);
            return new Dictionary<int, Matrix>{{id, result}};
        }


        public void OrderDofs(bool alsoOrderConstrainedDofs)
		{
			IGlobalFreeDofOrdering globalRowOrdering = dofRowOrderer.OrderFreeDofs(model);
			IGlobalFreeDofOrdering globalColOrdering = dofColOrderer.OrderFreeDofs(model);
			assembler.HandleDofOrderingWillBeModified();

			model.GlobalRowDofOrdering = globalRowOrdering;
			model.GlobalColDofOrdering = globalColOrdering;

			foreach (var subdomain in model.Subdomains)
			{
				subdomain.FreeDofRowOrdering = globalRowOrdering.SubdomainDofOrderings[subdomain];
				subdomain.FreeDofColOrdering = globalColOrdering.SubdomainDofOrderings[subdomain];

                if (alsoOrderConstrainedDofs)
                {
                    subdomain.ConstrainedDofRowOrdering= dofRowOrderer.OrderConstrainedDofs(subdomain);
                    subdomain.ConstrainedDofColOrdering = dofColOrderer.OrderConstrainedDofs(subdomain);
                }
			}
		}

        public abstract void Initialize();
        public abstract void HandleMatrixWillBeSet();
        public abstract void PreventFromOverwrittingSystemMatrices();
        public abstract void Solve();
        protected abstract Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix);
    }
}
