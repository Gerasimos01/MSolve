﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    public interface IScaleTransitions_v2
    {
        double[] MacroToMicroTransition(Node boundaryNode, double[] MacroScaleVariable);
        double[] MicroToMacroTransition(INode boundaryNode, double[] MicroScaleVariable);
        void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(Node boundaryNode,
            double[] MacroScaleVariable, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements);
        void ImposeAppropriateConstraintsPerBoundaryNode(Model_v2 model, Node boundaryNode);
        void ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(Model_v2 model_v2, Node boundaryNode, Dictionary<Node, IList<DOFType>> RigidBodyNodeConstraints); //TODO: enopoihsh twn duo duplicate
        int PrescribedDofsPerNode(); // TODO: pithanws epistrofh kai to poioi einai me input sugkekrimeno node
        int MacroscaleVariableDimension();
    }
}
