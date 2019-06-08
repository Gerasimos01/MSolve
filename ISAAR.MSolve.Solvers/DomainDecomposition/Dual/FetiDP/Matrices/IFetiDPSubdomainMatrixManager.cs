﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    /// <summary>
    /// Implements the linear algebra operations needed by <see cref="Feti1Solver"/> depending on the underlying matrix storage
    /// format. All the matrices represented by this interface belong to a single subdomain.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IFetiDPSubdomainMatrixManager : IFetiSubdomainMatrixManager
    {
        Matrix SchurComplementOfRemainderDofs { get; } //TODO: Perhaps static condensations should be handled by a different interface

        void CalcSchurComplementOfRemainderDofs();

        //TODO: E.g. Once Kcc* is calculated Kcc and Krc can be cleared. There are 2 options:
        //      a) Each matrix must be able to be cleared independently, if the FETI-DP solver and its strategies decide when.
        //      b) Otherwise this matrix decides when to clear what and these methods are optional/risky.
        //void ClearKcc();
        //void ClearKcrKrc();

        void ExtractKcc(int[] cornerDofs);  //TODO: perhaps SymmetricMatrix
        void ExtractKcrKrc(int[] cornerDofs, int[] remainderDofs); //TODO: perhaps CSR or CSC
        void ExtractKrr(int[] remainderDofs); //TODO: perhaps SkylineMatrix or SymmetricCSC

        void InvertKrr(bool inPlace);

        Vector MultiplyInverseKrrTimes(Vector vector);
        Vector MultiplyKccTimes(Vector vector);
        Vector MultiplyKcrTimes(Vector vector);
        Vector MultiplyKrcTimes(Vector vector);
    }
}