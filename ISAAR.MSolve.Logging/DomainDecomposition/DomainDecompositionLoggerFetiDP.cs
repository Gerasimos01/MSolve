﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;

namespace ISAAR.MSolve.Logging.DomainDecomposition
{
    public class DomainDecompositionLoggerFetiDP : IDomainDecompositionLogger
    {
        private readonly string plotDirectoryPath;
        private readonly FetiDPSolver solver;
        private int analysisStep;

        public DomainDecompositionLoggerFetiDP(string plotDirectoryPath, FetiDPSolver solver) //TODO: make sure the path does not end in "\"
        {
            this.plotDirectoryPath = plotDirectoryPath;
            this.solver = solver;
            analysisStep = 0;
        }

        public void PlotSubdomains(IStructuralModel model)
        {
            var writer = new MeshPartitionWriter();
            writer.WriteSubdomainElements($"{plotDirectoryPath}\\subdomains_{analysisStep}.vtk", model);
            writer.WriteBoundaryNodes($"{plotDirectoryPath}\\boundary_nodes_{analysisStep}.vtk", model);

            var allCornerNodes = new HashSet<INode>();
            foreach (INode[] cornerNodes in solver.CornerNodesOfSubdomains.Values) allCornerNodes.UnionWith(cornerNodes);
            writer.WriteSpecialNodes($"{plotDirectoryPath}\\corner_nodes_{analysisStep}.vtk", "corner_nodes", allCornerNodes);

            ++analysisStep;
        }
    }
}
