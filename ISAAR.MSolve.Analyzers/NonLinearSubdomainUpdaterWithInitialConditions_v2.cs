﻿using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearSubdomainUpdaterWithInitialConditions_v2 : INonLinearSubdomainUpdater_v2
	{
		private readonly Subdomain_v2 subdomain;

		public NonLinearSubdomainUpdaterWithInitialConditions_v2(Subdomain_v2 subdomain)
		{
			this.subdomain = subdomain;
		}

		public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, Node> boundaryNodes,
			Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements) //TODO leave 
		{
			return this.subdomain.GetRHSFromSolutionWithInitialDisplacemntsEffect(solution, dSolution, boundaryNodes,
			 initialConvergedBoundaryDisplacements, totalBoundaryDisplacements,
			 nIncrement, totalIncrements);
		}

		public void ResetState()
		{
			this.subdomain.ClearMaterialStresses();
		}

		public void UpdateState()
		{
			this.subdomain.SaveMaterialState();
		}

		public void ScaleConstraints(double scalingFactor)
		{
			throw new NotSupportedException();
		}

		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution) //TODO leave 
		{
			throw new NotSupportedException();
			return this.subdomain.GetRhsFromSolution(solution, dSolution);
		}
	}
}