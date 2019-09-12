﻿using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation
{
    public class FetiDPSubdomainGlobalMappingMpi
    {
        private readonly IStiffnessDistribution distribution;
        private readonly IFetiDPDofSeparator dofSeparator;
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        public FetiDPSubdomainGlobalMappingMpi(ProcessDistribution processDistribution, IModel model,
            IFetiDPDofSeparator dofSeparator, IStiffnessDistribution distribution)
        {
            this.procs = processDistribution;
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.distribution = distribution;
        }

        public double CalcGlobalForcesNorm(Func<ISubdomain, Vector> getSubdomainForces)
        {
            //TODO: This can be optimized: calculate the dot product f*f for the internal dofs of each subdomain separately,
            //      only assemble global vector for the boundary dofs, find its dot product with itself, add the contributions
            //      for the internal dofs and finally apply SQRT(). This would greatly reduce the communication requirements.
            //TODO: this should be used for non linear analyzers as well (instead of building the global RHS)
            //TODO: Is this correct? For the residual, it would be wrong to find f-K*u for each subdomain and then call this.

            Vector globalForces = AssembleSubdomainVectors(getSubdomainForces);
            double norm = double.NaN;
            if (procs.IsMasterProcess) norm = globalForces.Norm2();
            procs.Communicator.Broadcast(ref norm, procs.MasterProcess); //TODO: Not sure if this is needed.
            return norm;
        }

        public Vector GatherGlobalDisplacements(Func<ISubdomain, Vector> getSubdomainFreeDisplacements)
        {
            return AssembleSubdomainVectors(sub =>
            {
                Vector u = getSubdomainFreeDisplacements(sub);
                FetiDPSubdomainGlobalMappingUtilities.ScaleSubdomainFreeDisplacements(dofSeparator, distribution, sub, u);
                return u;
            });
        }

        private Vector AssembleSubdomainVectors(Func<ISubdomain, Vector> getSubdomainVector)
        {
            ISubdomain subdomain = model.GetSubdomain(procs.OwnSubdomainID);
            Vector subdomainVector = getSubdomainVector(subdomain);
            Vector globalVector = null;
            if (procs.IsMasterProcess) globalVector = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
            model.GlobalDofOrdering.AddVectorSubdomainToGlobal(subdomain, subdomainVector, globalVector);
            return globalVector;
        }
    }
}