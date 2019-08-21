﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MPI;

//TODO: The solver should decide which subdomains will be reused. This class only provides functionality.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Orders the unconstrained freedom degrees of each subdomain and the shole model. Also applies any reordering and other 
    /// optimizations.
    /// </summary>
    public class DofOrdererMpi : DofOrdererBase
    {
        private readonly ProcessDistribution procs;

        public DofOrdererMpi(ProcessDistribution processDistribution, IFreeDofOrderingStrategy freeOrderingStrategy,
            IDofReorderingStrategy reorderingStrategy, bool cacheElementToSubdomainDofMaps = true):
            base(freeOrderingStrategy, reorderingStrategy, cacheElementToSubdomainDofMaps)
        {
            this.procs = processDistribution;
        }

        public override void OrderFreeDofs(IModel model)
        {
            // Each process orders its subdomain dofs
            ISubdomain subdomain = model.GetSubdomain(procs.OwnSubdomainID);
            ISubdomainFreeDofOrdering subdomainOrdering = OrderFreeDofs(subdomain);
            subdomain.FreeDofOrdering = subdomainOrdering;

            // Order global dofs
            int numGlobalFreeDofs = -1;
            DofTable globalFreeDofs = null;
            if (procs.IsMasterProcess) (numGlobalFreeDofs, globalFreeDofs) = freeOrderingStrategy.OrderGlobalDofs(model);
            var globalOrdering = new GlobalFreeDofOrderingMpi(procs, numGlobalFreeDofs, globalFreeDofs, model);
            globalOrdering.CreateSubdomainGlobalMaps(model);
            model.GlobalDofOrdering = globalOrdering;
        }
    }
}
