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
        void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(Node boundaryNode,
            double[] MacroScaleVariable, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements);
        void ImposeAppropriateConstraintsPerBoundaryNode(Model model, Node boundaryNode);
        int PrescribedDofsPerNode();
        int MacroscaleVariableDimension();
    }
}
