﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.DomainDecomposition;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;

// TODO: fix a bug that happens when the crack has almost reached the boundary, is inside but no tip can be found
namespace ISAAR.MSolve.XFEM.Analyzers
{
    /// <summary>
    /// Implements crack propagation under static loading with linear material behavior. Based on Linear Elastic Fracture 
    /// Mechanics. Appropriate for brittle materials or fatigue crack propagation analysis. For now, it only works with XFEM.
    /// </summary>
    public class QuasiStaticCrackPropagationAnalyzer //: IAnalyzer
    {
        private readonly ICrackDescription crack;
        private readonly double fractureToughness;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly int maxIterations;
        private readonly XModel model;
        private readonly IStaticProvider problem;
        private readonly ISolver solver;

        public QuasiStaticCrackPropagationAnalyzer(XModel model, ISolver solver, IStaticProvider problem,
            ICrackDescription crack, double fractureToughness, int maxIterations)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
            this.problem = problem;
            this.crack = crack;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
        }

        public IDomainDecompositionLogger DDLogger { get; set; }
        public CrackPropagationTermination Termination { get; private set;}

        public void Initialize(bool isFirstAnalysis = true)
        {
            // The order in which the next initializations happen is very important.
            if (isFirstAnalysis) model.ConnectDataStructures();

            solver.Initialize(); //TODO: not sure about this one.
        }

        /// <summary>
        /// Returns the crack path after repeatedly executing: XFEM analysis, SIF calculation, crack propagation
        /// </summary>
        /// <returns></returns>
        public void Analyze()
        {
            int iteration;
            for (iteration = 0; iteration < maxIterations; ++iteration)
            {
                Debug.WriteLine($"Iteration {iteration}");

                // Apply the updated enrichements.
                crack.UpdateEnrichments();

                // Plot domain decomposition data, if necessary
                if (DDLogger != null) DDLogger.PlotSubdomains(model);

                // Order and count dofs
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }

                // Create the stiffness matrix and then the forces vector
                problem.ClearMatrices();
                BuildMatrices();
                model.AssignLoads(solver.DistributeNodalLoads);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                }
                AddEquivalentNodalLoadsToRhs();

                // Solve the linear system
                solver.Solve();
                //Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DofOrderer);
                var freeDisplacements = new Dictionary<int, Vector>();
                foreach (int s in linearSystems.Keys) freeDisplacements[s] = (Vector)(linearSystems[s].Solution); //TODO: avoid this cast.

                //// Output field data
                //if (fieldOutput != null)
                //{
                //    fieldOutput.WriteOutputData(solver.DofOrderer, freeDisplacements, constrainedDisplacements, iteration);
                //}

                // Let the crack propagate
                crack.Propagate(freeDisplacements);

                // Check convergence 
                //TODO: Perhaps this should be done by the crack geometry or the Propagator itself and handled via exceptions 

                foreach (var tipPropagator in crack.CrackTipPropagators)
                {
                    double sifEffective = CalculateEquivalentSIF(tipPropagator.Value.Logger.SIFsMode1[iteration],
                    tipPropagator.Value.Logger.SIFsMode2[iteration]);
                    if (sifEffective >= fractureToughness)
                    {
                        Termination = CrackPropagationTermination.FractureToughnessIsExceeded;
                        return;
                    }
                    if (!model.Boundary.IsInside(tipPropagator.Key))
                    {
                        Termination = CrackPropagationTermination.CrackExitsDomainBoundary;
                        return;
                    }
                }
            }
            Termination = CrackPropagationTermination.RequiredIterationsWereCompleted;
        }

        private void AddEquivalentNodalLoadsToRhs()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                try
                {
                    // Make sure there is at least one non zero prescribed displacement.
                    (INode node, IDofType dof, double displacement) = linearSystem.Subdomain.Constraints.Find(du => du != 0.0);

                    //TODO: the following 2 lines are meaningless outside diplacement control (and even then, they are not so clear).
                    double scalingFactor = 1;
                    IVector initialFreeSolution = linearSystem.CreateZeroVector();

                    IVector equivalentNodalLoads = problem.DirichletLoadsAssembler.GetEquivalentNodalLoads(
                        linearSystem.Subdomain, initialFreeSolution, scalingFactor);
                    linearSystem.RhsVector.SubtractIntoThis(equivalentNodalLoads);
                }
                catch (KeyNotFoundException)
                {
                    // There aren't any non zero prescribed displacements, therefore we do not have to calculate the equivalent 
                    // nodal loads, which is an expensive operation (all elements are accessed, their stiffness is built, etc..)
                }
            }
        }

        private void BuildMatrices()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = problem.CalculateMatrix(linearSystem.Subdomain);
            }
        }

        // TODO: Abstract this and add Tanaka_1974 approach
        private double CalculateEquivalentSIF(double sifMode1, double sifMode2)
        {
            return Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
        }
    }
}
