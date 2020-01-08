﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Exceptions;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.DofSeparation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using MPI;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.FlexibilityMatrix;
using ISAAR.MSolve.Solvers.Logging;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Displacements;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.FlexibilityMatrix;

//TODO: Add time logging
//TODO: Use a base class for the code that is identical between FETI-1 and FETI-DP.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d
{
    public class FetiDP3dSolverSerial : ISolverMpi
    {
        internal const string name = "FETI-DP Solver"; // for error messages and logging
        private readonly IAugmentationConstraints augmentationConstraints;
        private readonly IFreeDofDisplacementsCalculator displacementsCalculator;
        private readonly DofOrderer dofOrderer;
        private readonly FetiDPDofSeparatorSerial dofSeparator;
        private readonly IFetiDPInterfaceProblemSolver interfaceProblemSolver;
        private readonly LagrangeMultipliersEnumeratorSerial lagrangesEnumerator;
        private readonly FetiDP3dMatrixManagerSerial matrixManager;
        private readonly IModel model;
        private readonly string msgHeader;
        private readonly bool problemIsHomogeneous;
        private readonly IFetiPreconditionerFactory precondFactory;
        private readonly IFetiPreconditioningOperations preconditioning; //TODO: perhaps this should be hidden inside IFetiPreconditionerFactory
        private readonly IStiffnessDistribution stiffnessDistribution;
        private readonly FetiDPSubdomainGlobalMappingSerial subdomainGlobalMapping;

        private bool factorizeInPlace = true;
        private IFetiDPFlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private IFetiPreconditioner preconditioner;

        public FetiDP3dSolverSerial(IModel model, ICornerNodeSelection cornerNodeSelection, IMidsideNodesSelection midsideNodesSelection, IAugmentationConstraintsFactory augmentationConstraintsFactory,
            IFetiDP3dMatrixManagerFactory matrixManagerFactory, IFetiPreconditioningOperations preconditioning,
            ICrosspointStrategy crosspointStrategy, PcgSettings pcgSettings, bool problemIsHomogeneous)
        {
            this.msgHeader = $"{this.GetType().Name}: ";

            if (model.NumSubdomains == 1) throw new InvalidSolverException(msgHeader 
                + $"This solver cannot be used if there is only 1 subdomain");
            this.model = model;

            this.Logger = new SolverLoggerSerial(name);

            // Connectivity
            this.CornerNodes = cornerNodeSelection;
            this.dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            this.dofSeparator = new FetiDPDofSeparatorSerial(model, cornerNodeSelection);
            this.lagrangesEnumerator = new LagrangeMultipliersEnumeratorSerial(model, crosspointStrategy, dofSeparator);
            this.augmentationConstraints = augmentationConstraintsFactory.CreateAugmentationConstraints(model, 
                midsideNodesSelection, dofSeparator, lagrangesEnumerator);

            // Matrix managers and linear systems
            this.matrixManager = new FetiDP3dMatrixManagerSerial(model, this.dofSeparator, lagrangesEnumerator, 
                augmentationConstraints ,matrixManagerFactory);
            //TODO: This will call HandleMatrixWillBeSet() once for each subdomain. For now I will clear the data when 
            //      BuildMatrices() is called. Redesign this.
            //matrixManager.LinearSystem.MatrixObservers.Add(this); 

            // Preconditioning
            this.preconditioning = preconditioning;
            this.precondFactory = new FetiPreconditionerSerial.Factory();

            // Interface problem
            this.interfaceProblemSolver = new FetiDP3dInterfaceProblemSolverSerial(model, pcgSettings,augmentationConstraints);
            this.displacementsCalculator = new FetiDP3dFreeDofDisplacementsCalculatorSerial(model, dofSeparator, lagrangesEnumerator,
                augmentationConstraints, matrixManager);

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
            if (problemIsHomogeneous)
            {
                this.stiffnessDistribution = new HomogeneousStiffnessDistributionSerial(model, dofSeparator,
                    new FetiDPHomogeneousDistributionLoadScaling(dofSeparator));
            }
            else throw new NotImplementedException();

            this.subdomainGlobalMapping = new FetiDPSubdomainGlobalMappingSerial(model, dofSeparator, stiffnessDistribution);
        }

        public ICornerNodeSelection CornerNodes { get; }
        public ISolverLogger Logger { get; }
        public string Name => name;
        public INodalLoadDistributor NodalLoadDistributor => stiffnessDistribution;

        /// <summary>
        ///  builds Kff of each subdomain
        /// </summary>
        /// <param name="elementMatrixProvider"></param>
        public void BuildGlobalMatrix(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            Logger.StartMeasuringTime();

            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader
                        + $" Assembling the free-free stiffness matrix of subdomain {subdomain.ID}");
                    IFetiDPSubdomainMatrixManager subdomainMatrices = matrixManager.GetFetiDPSubdomainMatrixManager(subdomain);
                    subdomainMatrices.BuildFreeDofsMatrix(subdomain.FreeDofOrdering, elementMatrixProvider); 
                }
            }

            Logger.LogCurrentTaskDuration("Matrix assembly");

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
        }

        public Vector GatherGlobalDisplacements()
        {
            return subdomainGlobalMapping.GatherGlobalDisplacements(
                sub => matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.SolutionConcrete);
        }

        public ILinearSystem GetLinearSystem(ISubdomain subdomain)
            => matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).LinearSystem;
        
        public void HandleMatrixWillBeSet()
        {
            isStiffnessModified = true;
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine(msgHeader + $"Clearing saved matrices of subdomain {subdomain.ID}.");
                    matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ClearMatrices();
                }
            }

            flexibility = null;
            preconditioner = null;
            matrixManager.ClearInverseCoarseProblemMatrix();
        }

        public void Initialize()
        {
            Logger.StartMeasuringTime();

            CornerNodes.Update();

            // Define the various dof groups
            dofSeparator.SeparateDofs(matrixManager);
            //FetiDPDofSeparationLogging.PrintDofSeparationSerial(model, dofSeparator);

            //TODO: B matrices could also be reused in some cases
            // Define lagrange multipliers and boolean matrices. 
            lagrangesEnumerator.CalcBooleanMatrices(dofSeparator.GetRemainderDofOrdering);

            //Define augmentation constraints and boolean matrices.
            augmentationConstraints.CalcAugmentationMappingMatrices();


            // Log dof statistics
            Logger.LogCurrentTaskDuration("Dof ordering");
            Logger.LogNumDofs("Lagrange multipliers", lagrangesEnumerator.NumLagrangeMultipliers);
            Logger.LogNumDofs("Corner dofs", dofSeparator.NumGlobalCornerDofs);
            Logger.LogNumDofs("AugmentationConstraints", augmentationConstraints.NumGlobalAugmentationConstraints);

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            stiffnessDistribution.Update();
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs)
        {
            Logger.StartMeasuringTime();

            // Order dofs
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                if (subdomain.ConnectivityModified)
                {
                    matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).HandleDofOrderingWillBeModified(); //TODO: Not sure about this
                }
            }

            // This should not create subdomain-global mappings which require MPI communication
            //TODO: What about subdomain-global mappings, especially for boundary dofs? Who should create them? 
            dofOrderer.OrderFreeDofs(model);

            if (alsoOrderConstrainedDofs)
            {
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);
                }
            }

            // Log dof statistics
            Logger.LogCurrentTaskDuration("Dof ordering");
            Logger.LogNumDofs("Global dofs", model.GlobalDofOrdering.NumGlobalFreeDofs);
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            if (isStiffnessModified)
            {
                // Separate the stiffness matrix
                Logger.StartMeasuringTime();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ExtractCornerRemainderSubmatrices();
                    }
                }
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                // Calculate the preconditioner before factorizing each subdomain's Krr.
                // The inter-subdomain stiffness distribution may have changed even if a subdomain's stiffness is the same.
                Logger.StartMeasuringTime();
                if (preconditioning.ReorderInternalDofsForFactorization) dofSeparator.ReorderInternalDofs(matrixManager);
                preconditioner = precondFactory.CreatePreconditioner(preconditioning, model, dofSeparator, lagrangesEnumerator,
                    matrixManager, stiffnessDistribution);
                Logger.LogCurrentTaskDuration("Calculating preconditioner");



                // Factorize each subdomain's Krr
                Logger.StartMeasuringTime();
                foreach (ISubdomain subdomain in model.EnumerateSubdomains())
                {
                    if (subdomain.StiffnessModified)
                    {
                        //TODO: If I can reuse Krr, I can also reuse its factorization. Therefore this must be inPlace. In contrast, FETI-1 needs Kff intact for Stiffness distribution, in the current design).
                        Debug.WriteLine(msgHeader
                            + $"Inverting the remainder-remainder stiffness matrix of subdomain {subdomain.ID} in place.");
                        matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).InvertKrr(true);
                    }
                }

                // Calculate FETI-DP coarse problem matrix
                matrixManager.CalcInverseCoarseProblemMatrix(CornerNodes);
                flexibility = new FetiDP3dFlexibilityMatrixSerial(model, dofSeparator, lagrangesEnumerator, 
                    augmentationConstraints, matrixManager);
                Logger.LogCurrentTaskDuration("Calculating coarse problem matrix");

                isStiffnessModified = false;
            }

            // Calculate FETI-DP coarse problem rhs 
            //TODO: rename this it is not coarse problem rhs only
            Logger.StartMeasuringTime();
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                matrixManager.GetFetiDPSubdomainMatrixManager(subdomain).ExtractCornerRemainderRhsSubvectors();
            }
            matrixManager.CalcCoarseProblemRhs();
            Logger.LogCurrentTaskDuration("Calculating coarse problem rhs");

            Logger.StartMeasuringTime();
            // Calculate the norm of the forces vector |f| = |K*u|. It is needed to check the convergence of PCG.
            double globalForcesNorm = globalForcesNorm = subdomainGlobalMapping.CalcGlobalForcesNorm(
                    sub => matrixManager.GetFetiDPSubdomainMatrixManager(sub).LinearSystem.RhsConcrete);

            // Solve interface problem
            Vector lagranges = interfaceProblemSolver.SolveInterfaceProblem(matrixManager, lagrangesEnumerator,
                flexibility, preconditioner, globalForcesNorm, Logger);
            Logger.LogCurrentTaskDuration("Solving interface problem");

            // Calculate the displacements of each subdomain
            Logger.StartMeasuringTime();
            displacementsCalculator.CalculateSubdomainDisplacements(lagranges, flexibility);
            Logger.LogCurrentTaskDuration("Calculate displacements from lagrange multipliers");

            Logger.IncrementAnalysisStep();
        }

        public class Builder
        {
            private readonly IFetiDP3dMatrixManagerFactory matrixManagerFactory;

            public Builder(IFetiDP3dMatrixManagerFactory matrixManagerFactory)
            {
                this.matrixManagerFactory = matrixManagerFactory;
            }

            public IAugmentationConstraintsFactory AugmentationConstraintsFactory { get; set; } = new AugmentationConstraints.Factory();
            
            public ICrosspointStrategy CrosspointStrategy { get; set; } = new FullyRedundantConstraints();

            public IDofOrderer DofOrderer { get; set; } =
                new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public PcgSettings PcgSettings { get; set; } = new PcgSettings();

            public IFetiPreconditioningOperations Preconditioning { get; set; } = new DirichletPreconditioning();

            public bool ProblemIsHomogeneous { get; set; } = true;

            public FetiDP3dSolverSerial Build(IModel model, ICornerNodeSelection cornerNodeSelection,
                IMidsideNodesSelection midsideNodesSelection)
            {
                return new FetiDP3dSolverSerial(model, cornerNodeSelection, midsideNodesSelection, 
                    AugmentationConstraintsFactory, matrixManagerFactory, Preconditioning, CrosspointStrategy,
                    PcgSettings, ProblemIsHomogeneous);
            }
        }
    }
}
