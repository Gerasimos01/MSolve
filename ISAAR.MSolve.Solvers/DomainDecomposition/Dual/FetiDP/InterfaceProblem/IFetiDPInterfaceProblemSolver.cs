﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;

//TODO: This could be split into an interface with the same name and an IFtiDPCoarseProblemSolver.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem
{
    public interface IFetiDPInterfaceProblemSolver
    {
        (Vector lagrangeMultipliers, Vector cornerDisplacements) SolveInterfaceProblem(FetiDPFlexibilityMatrixOLD flexibility, 
            IFetiPreconditioner preconditioner, IFetiDPCoarseProblemSolverOLD coarseProblemSolver, Vector globalFcStar, Vector dr,
            double globalForcesNorm, SolverLogger logger);
    }
}
