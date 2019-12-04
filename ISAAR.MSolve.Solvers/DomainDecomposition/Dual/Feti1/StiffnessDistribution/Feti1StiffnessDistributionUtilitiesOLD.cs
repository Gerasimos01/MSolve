﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.DofSeparation;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    internal static class Feti1StiffnessDistributionUtilitiesOLD
    {
        //TODO: This was part of a now obsolete design. However it was more efficient than the current one for serial code.
        public static Dictionary<int, SparseVector> DistributeNodalLoadsOLD(Feti1DofSeparator dofSeparator,
            Dictionary<int, ISubdomain> subdomains, Table<INode, IDofType, double> globalNodalLoads,
           Func<INode, IDofType, Dictionary<int, double>> calcBoundaryDofCoefficients)
        {
            //TODO: Should I implement this as fb(s) = Lpb(s) * fb, with a) Lpb(s) = Lb(s) * inv(Mb) for homogeneous and 
            //      b) Lpb(s) = Db(s)*Lb(s) * inv(Lb^T*Db*Lb) for heterogeneous?

            var subdomainLoads = new Dictionary<int, SortedDictionary<int, double>>();
            foreach (var subdomainID in subdomains.Keys) subdomainLoads[subdomainID] = new SortedDictionary<int, double>();

            foreach ((INode node, IDofType dofType, double loadAmount) in globalNodalLoads)
            {
                if (node.Multiplicity == 1) // optimization for internal dof
                {
                    ISubdomain subdomain = node.SubdomainsDictionary.First().Value;
                    int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                    subdomainLoads[subdomain.ID][subdomainDofIdx] = loadAmount;
                }
                else // boundary dof
                {
                    Dictionary<int, double> boundaryDofCoeffs = calcBoundaryDofCoefficients(node, dofType);
                    foreach (var idSubdomain in node.SubdomainsDictionary)
                    {
                        int id = idSubdomain.Key;
                        ISubdomain subdomain = idSubdomain.Value;
                        int subdomainDofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        subdomainLoads[id][subdomainDofIdx] = loadAmount * boundaryDofCoeffs[id];
                    }
                }
            }

            var vectors = new Dictionary<int, SparseVector>();
            foreach (var idSubdomains in subdomains)
            {
                int id = idSubdomains.Key;
                int numSubdomainDofs = idSubdomains.Value.FreeDofOrdering.NumFreeDofs;
                vectors[id] = SparseVector.CreateFromDictionary(numSubdomainDofs, subdomainLoads[id]);
            }
            return vectors;
        }
    }
}