﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public abstract class NonLinearAnalyzerBuilderBaseDevelop
    {
        protected int maxIterationsPerIncrement = 1000;
        protected readonly IModel model;
        protected readonly int numIncrements;
        protected int numIterationsForMatrixRebuild = 1;
        protected readonly INonLinearProvider provider;
        protected double residualTolerance = 1e-8;
        protected readonly ISolver solver;

        protected NonLinearAnalyzerBuilderBaseDevelop(IModel model, ISolver solver, INonLinearProvider provider, 
            int numIncrements)
        {
            //TODO: this should belong to all (child) analyzer builders
            this.model = model;
            this.solver = solver;
            this.provider = provider;
            this.numIncrements = numIncrements;
            SubdomainUpdaters = CreateDefaultSubdomainUpdaters();
        }

        public int MaxIterationsPerIncrement
        {
            get => maxIterationsPerIncrement;
            set
            {
                if (value < 1) throw new ArgumentException("Max iterations per increment must be >= 1");
                maxIterationsPerIncrement = value;
            }
        }

        public int NumIterationsForMatrixRebuild
        {
            get => numIterationsForMatrixRebuild;
            set
            {
                if (value < 1) throw new ArgumentException("Iterations number for matrix rebuild must be >= 1");
                numIterationsForMatrixRebuild = value;
            }
        }

        public double ResidualTolerance
        {
            get => residualTolerance;
            set
            {
                if (value <= 0.0) throw new ArgumentException("Residual tolerance must be positive");
                residualTolerance = value;
            }
        }

        public IReadOnlyDictionary<int, INonLinearSubdomainUpdater> SubdomainUpdaters { get; set; }

        private IReadOnlyDictionary<int, INonLinearSubdomainUpdater> CreateDefaultSubdomainUpdaters()
        {
            int numSubdomains = model.NumSubdomains;
            var subdomainUpdaters = new Dictionary<int, INonLinearSubdomainUpdater>(numSubdomains);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                subdomainUpdaters[subdomain.ID] = new NonLinearSubdomainUpdaterDevelop(subdomain);

            }
            return subdomainUpdaters;

        }
    }
}