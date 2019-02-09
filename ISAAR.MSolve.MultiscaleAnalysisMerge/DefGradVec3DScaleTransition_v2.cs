﻿using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class DefGradVec3DScaleTransition_v2 : IScaleTransitions_v2
    {
        public DefGradVec3DScaleTransition_v2()
        { }

        public double[] MacroToMicroTransition(Node boundaryNode, double[] MacroScaleVariable)
        {
            double[,] Dq_nodal = new double[9, 3];
            Dq_nodal[0, +0] = boundaryNode.X; // h kai katedtheian boundaryNode.X 
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +2] = boundaryNode.Z;
            Dq_nodal[3, +0] = boundaryNode.Y;
            Dq_nodal[4, +1] = boundaryNode.Z;
            Dq_nodal[5, +2] = boundaryNode.X;
            Dq_nodal[6, +0] = boundaryNode.Z;
            Dq_nodal[7, +1] = boundaryNode.X;
            Dq_nodal[8, +2] = boundaryNode.Y;

            double[] microVariable = new double[3];            

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 9; j1++)
                {
                    microVariable[i1] += Dq_nodal[j1, i1] * MacroScaleVariable[j1]; //einai sunolikh 
                }
            }

            return microVariable;
        }

        public double[] MicroToMacroTransition(INode boundaryNode, double[] MicroScaleVariable)
        {
            double[,] Dq_nodal = new double[9, 3];
            Dq_nodal[0, +0] = boundaryNode.X; // h kai katedtheian boundaryNode.X 
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +2] = boundaryNode.Z;
            Dq_nodal[3, +0] = boundaryNode.Y;
            Dq_nodal[4, +1] = boundaryNode.Z;
            Dq_nodal[5, +2] = boundaryNode.X;
            Dq_nodal[6, +0] = boundaryNode.Z;
            Dq_nodal[7, +1] = boundaryNode.X;
            Dq_nodal[8, +2] = boundaryNode.Y;

            double[] macroVariable = new double[9];

            for (int i1 = 0; i1 < 9; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    macroVariable[i1] += Dq_nodal[ i1, j1] * MicroScaleVariable[j1]; //einai sunolikh 
                }
            }

            return macroVariable;
        }

        public int PrescribedDofsPerNode()
        {
            return 3;
        }

        public int MacroscaleVariableDimension()
        {
            return 9;
        }

        public void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(Node boundaryNode,
            double[] DefGradVec, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements)
        {
            double[,] Dq_nodal = new double[9, 3];
            Dq_nodal[0, +0] = boundaryNode.X; // h kai katedtheian boundaryNode.X 
            Dq_nodal[1, +1] = boundaryNode.Y;
            Dq_nodal[2, +2] = boundaryNode.Z;
            Dq_nodal[3, +0] = boundaryNode.Y;
            Dq_nodal[4, +1] = boundaryNode.Z;
            Dq_nodal[5, +2] = boundaryNode.X;
            Dq_nodal[6, +0] = boundaryNode.Z;
            Dq_nodal[7, +1] = boundaryNode.X;
            Dq_nodal[8, +2] = boundaryNode.Y;

            double[] thesi_prescr_xyz = new double[3];
            double[] u_prescr_xyz_sunol = new double[3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 9; j1++)
                {
                    thesi_prescr_xyz[i1] += Dq_nodal[j1, i1] * DefGradVec[j1]; //einai sunolikh 
                }
            }
            u_prescr_xyz_sunol = new double[3] { thesi_prescr_xyz[0] - boundaryNode.X,
                                                     thesi_prescr_xyz[1] - boundaryNode.Y,
                                                     thesi_prescr_xyz[2] - boundaryNode.Z };

            Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
            totalBoundaryNodalDisplacements.Add(DOFType.X, u_prescr_xyz_sunol[0]);
            totalBoundaryNodalDisplacements.Add(DOFType.Y, u_prescr_xyz_sunol[1]);
            totalBoundaryNodalDisplacements.Add(DOFType.Z, u_prescr_xyz_sunol[2]);

            totalPrescribedBoundaryDisplacements.Add(boundaryNode.ID, totalBoundaryNodalDisplacements);
        }

        public void ImposeAppropriateConstraintsPerBoundaryNode(Model_v2 model, Node boundaryNode)
        {
            model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint { DOF = DOFType.X });
            model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint { DOF = DOFType.Y });
            model.NodesDictionary[boundaryNode.ID].Constraints.Add(new Constraint { DOF = DOFType.Z });
        }

        public void ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(Model_v2 model, Node boundaryNode, Dictionary<Node, IList<DOFType>> RigidBodyNodeConstraints)
        {
            throw new System.NotSupportedException();
        }
    }
}
