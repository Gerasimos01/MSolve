using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IScaleTransitions
    {
        double[] MacroToMicroTransition(Node boundaryNode, double[] MacroScaleVariable);
        double[] MicroToMacroTransition(Node boundaryNode, double[] MicroScaleVariable);
        int PrescribedDofsPerNode();
    }
}
