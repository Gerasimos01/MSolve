﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests.Mocks;
using ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.Example4x4x4Quads;
using ISAAR.MSolve.Solvers.Tests.Utilities;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP3d.UnitTests
{
    public static class FetiDP3dFlexibilityMatrixSerialTests
    {
        [Fact]
        public static void TestFlexibilityMatrices()
        {
            (IModel model, FetiDPDofSeparatorSerial dofSeparator, LagrangeMultipliersEnumeratorSerial lagrangesEnumerator) =
                FetiDP3dLagrangesEnumeratorSerialTests.CreateModelDofSeparatorLagrangesEnumerator();
            IAugmentationConstraints augmentationConstraints =
                FetiDP3dAugmentedConstraintsTests.CalcAugmentationConstraintsSimple(model, lagrangesEnumerator);

            // Setup matrix manager
            IFetiDP3dMatrixManager matrixManager = new MockMatrixManager(model);

            // Create explicit matrices that can be checked
            var flexibility = new FetiDP3dFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator,
                augmentationConstraints, matrixManager);
            int coarseProblemSize = dofSeparator.NumGlobalCornerDofs + augmentationConstraints.NumGlobalAugmentationConstraints;
            int numLagranges = lagrangesEnumerator.NumLagrangeMultipliers;
            Matrix FIrr = ImplicitMatrixUtilities.MultiplyWithIdentity(
                numLagranges, numLagranges, flexibility.MultiplyGlobalFIrr);
            Matrix FIrc = ImplicitMatrixUtilities.MultiplyWithIdentity(
                numLagranges, coarseProblemSize, flexibility.MultiplyGlobalFIrc);
            Matrix FIrcTranspose = ImplicitMatrixUtilities.MultiplyWithIdentity(
                coarseProblemSize, numLagranges, flexibility.MultiplyGlobalFIrcTransposed); //TODO: Add this test in FETI-DP 2D (serial, MPI) as well

            // Check
            double tol = 1E-3;
            Assert.True(ExpectedGlobalMatrices.MatrixFIrr.Equals(FIrr, tol));
            Assert.True(ExpectedGlobalMatrices.MatrixFIrcTildeSimple.Equals(FIrc, tol));
            Assert.True(ExpectedGlobalMatrices.MatrixFIrcTildeSimple.Transpose().Equals(FIrcTranspose, tol));
        }
    }
}
