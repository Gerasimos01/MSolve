using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class DebugMethods
    {
        [Fact]
        public static void WriteInputFilesForTetRve()
        {
            var outterMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal
            var innerMaterial = new ElasticMaterial3DtotalStrain() { PoissonRatio = 0.4, YoungModulus = 3.5e9 }; //pascal
            IdegenerateRVEbuilder rveBuilder = new CompositeMaterialModeluilderTet(outterMaterial, innerMaterial, 100, 100, 100);

            Tuple<Model, Dictionary<int, Node>, double> modelAndBoundaryNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = modelAndBoundaryNodes.Item1;
            Dictionary<int, Node> boundaryNodes = modelAndBoundaryNodes.Item2;
            Dictionary<Node, IList<IDofType>> RigidBodyNodeConstraints = rveBuilder.GetModelRigidBodyNodeConstraints(model);

            // Create constraints
            Dictionary<Node, IList<IDofType>> RigidBodyConstraintsToPrint = new Dictionary<Node, IList<IDofType>>();
            foreach (Node node in RigidBodyNodeConstraints.Keys)
            { RigidBodyConstraintsToPrint.Add(node, new List<IDofType>() { StructuralDof.TranslationZ });}

            //Write Constraints
            string OutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi_2\develop_nl_iga_shell\Adina_Input\";
            string constraints_output = OutputPath_gen + "constraints_for_adina_input.txt";
            WriteConstraints(model, RigidBodyConstraintsToPrint, constraints_output);

            //Create displacement LOADs
            double[] smallStrainVec = new double[3] { 0.01, 0, 0 };
            IScaleTransitions scaleTransitions = new SmallStrain3Dto2DplaneStressScaleTransition();
            Dictionary<int, Dictionary<IDofType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<IDofType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, totalPrescribedBoundaryDisplacements);
            }

            //write displacement loads
            string loads_output = OutputPath_gen + "loads_for_adina_input.txt";
            WriteDisplacementLoads(totalPrescribedBoundaryDisplacements, 1, loads_output);

        }

        public static void WriteConstraints(Model model, Dictionary<Node, IList<IDofType>> RigidBodyConstraintsToPrint, string path)
        {
            StringBuilder s = new StringBuilder();
            foreach(Node node in RigidBodyConstraintsToPrint.Keys)
            {
                foreach (var doftype in RigidBodyConstraintsToPrint[node])
                {
                    if (doftype == StructuralDof.TranslationX)
                    {
                        s.Append(node.ID + "\t" + "Fixed" + "\t" + "Free" + "\t" + "Free");
                        s.AppendLine();
                    }
                    if (doftype == StructuralDof.TranslationY)
                    {
                        s.Append(node.ID + "\t" + "Free" + "\t" + "Fixed" + "\t" + "Free");
                        s.AppendLine();
                    }
                    
                    if (doftype == StructuralDof.TranslationZ)
                    { s.Append(node.ID + "\t" + "Free" + "\t" + "Free" + "\t" + "Fixed");
                        s.AppendLine();
                    }
                    
                }
            }


            var writer = new StreamWriter(path);
            writer.Write(s.ToString());
            writer.Flush();
            writer.Dispose();
            

            //sb.AppendLine($"Iterations required = {IterationsRequired}");
            //text = text.Replace("@", "@" + System.Environment.NewLine);

        }

        public static void WriteDisplacementLoads(Dictionary<int, Dictionary<IDofType, double>> totalPrescribedBoundaryDisplacements, int TimeFunctionID, string path)
        {
            StringBuilder s = new StringBuilder();
            
            foreach(int nodeId in totalPrescribedBoundaryDisplacements.Keys)
            {
                foreach (var doftype in totalPrescribedBoundaryDisplacements[nodeId].Keys)
                {
                    if (doftype == StructuralDof.TranslationX)
                    {
                        s.Append(nodeId + "\t" + "X-Translation" + "\t" + totalPrescribedBoundaryDisplacements[nodeId][doftype]  + "\t" + TimeFunctionID);
                        s.AppendLine();
                    }
                    if (doftype == StructuralDof.TranslationY)
                    {
                        s.Append(nodeId + "\t" + "Y-Translation" + "\t" + totalPrescribedBoundaryDisplacements[nodeId][doftype] + "\t" + TimeFunctionID);
                        s.AppendLine();
                    }
                    if (doftype == StructuralDof.TranslationZ)
                    {
                        s.Append(nodeId + "\t" + "Z-Translation" + "\t" + totalPrescribedBoundaryDisplacements[nodeId][doftype] + "\t" + TimeFunctionID);
                        s.AppendLine();
                    }

                }
            }

            var writer = new StreamWriter(path);
            writer.Write(s.ToString());
            writer.Flush();
            writer.Dispose();


        }

    }
}
