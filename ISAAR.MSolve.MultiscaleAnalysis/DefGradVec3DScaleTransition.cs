using ISAAR.MSolve.PreProcessor;
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


namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class DefGradVec3DScaleTransition : IScaleTransitions
    {
        public DefGradVec3DScaleTransition()
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

        public double[] MicroToMacroTransition(Node boundaryNode, double[] MicroScaleVariable)
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
    }
}
