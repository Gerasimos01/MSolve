using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates the nesessary methods that should be implemented by builders of models intended to be used as rves in Multiscale Problems
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IRVEbuilderDevelop
    {
        Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes();
        IRVEbuilderDevelop Clone(int a);
        ISolver GetAppropriateSolver(Model model);
        void ReassignModel(Model model,bool decomposeModel);
        ISolver GetParallelSolver(Model model);
        ISolver GetSerialSolver(Model model);

    }
}
