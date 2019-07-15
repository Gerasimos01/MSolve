using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Multiscale;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass5b_b_cantilever
    {
        //original:SeparateCodeCheckingClass5b_b
        //changes: cantilever geometry

        //needs to be corrected rve_multiple -> b kai to path kai ta stoixeia diakritopoihshs pou einai afhmena exwterika (Genika elegxoume connectDataStructures kai defineAppropriateConstraintsForBoundaryNodes)
        public static (Model, double[]) RunExampleSerial()
        {
            //var rveBuilder = new RveGrShMultipleSeparatedDevelopbLARGE(1, false); // diorthose kai to parakatw path apla gia na mhn xtupaei.
            var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop5elemCantilever(1, false); //A.1

            var ModelAndNodes = rveBuilder.GetModelAndBoundaryNodes();
            Model model = ModelAndNodes.Item1;
            var boundaryNodes = ModelAndNodes.Item2;

            #region corner nodes Data
            var RVE_id = 1;
            var rve_id_data = RVE_id.ToString();

            //renumbering_vector_path = "..\\..\\..\\RveTemplates\\Input\\RveGrShMultiple\\rve_no_{0}\\REF_new_total_numbering.txt";


            //var renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2bLARGE\RVE_database\rve_no_{0}\REF_new_total_numbering.txt";
            var renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\2c_alte_5elemCantilever\RVE_database\rve_no_{0}\REF_new_total_numbering.txt"; //A.2
            renumbering_vector_path = string.Format(renumbering_vector_path, rve_id_data);

            //A.3
            //int subdiscr1 = 3;
            //int discr1 = 5;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 15;
            //int subdiscr1_shell = 6;
            //int discr1_shell = 1;
            int subdiscr1 = 2;
            int discr1 = 1;
            // int discr2 dn xrhsimopoieitai 
            int discr3 = 2; //A.9 prosoxh kai alles allages 
            int subdiscr1_shell = 1;
            int discr1_shell = 1;
            var mpgp = FEMMeshBuilder.GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            //A.9
            var mp = mpgp.Item1; mp.hexa1 = 4; mp.hexa2 = 1; mp.hexa3 = 1;
            var gp = mpgp.Item2;
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumbering_vector_path));

            //A.4
            //mp.L01 = 120; mp.L02 = 120; mp.L03 = 120;
            mp.L01 = 100; mp.L02 = 25; mp.L03 = 25;
            gp = mpgp.Item2;
            gp.elem1 = 4; gp.elem2 = 1;


            double L01 = mp.L01; double L02 = mp.L02; double L03 = mp.L03;
            int hexa1 = mp.hexa1; int hexa2 = mp.hexa2; int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            Dictionary<int, double[]> CornerNodesIds;
            Dictionary<int, int[]> CornerNodesIdAndsubdomains;

            //A.5
            //int[][] CornerNodesData = new int[8][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            //CornerNodesData[0] = new int[3] { 5 - 1, 5 - 1, 5 - 1 };
            //CornerNodesData[1] = new int[3] { 9 - 1, 5 - 1, 5 - 1 };
            //CornerNodesData[2] = new int[3] { 5 - 1, 9 - 1, 5 - 1 };
            //CornerNodesData[3] = new int[3] { 9 - 1, 9 - 1, 5 - 1 };
            //CornerNodesData[4] = new int[3] { 5 - 1, 5 - 1, 9 - 1 };
            //CornerNodesData[5] = new int[3] { 9 - 1, 5 - 1, 9 - 1 };
            //CornerNodesData[6] = new int[3] { 5 - 1, 9 - 1, 9 - 1 };
            //CornerNodesData[7] = new int[3] { 9 - 1, 9 - 1, 9 - 1 };
            int[][] CornerNodesData = new int[3][];// new int[27][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            CornerNodesData[0] = new int[3] { 1, 0, 0 };
            CornerNodesData[1] = new int[3] { 1, 1, 0 };
            CornerNodesData[2] = new int[3] { 1, 0, 1 };




            CornerNodesIds = new Dictionary<int, double[]>(CornerNodesData.Length); //nodeID, coordinate data
            CornerNodesIdAndsubdomains = new Dictionary<int, int[]>(CornerNodesData.Length);//nodeID, subdIds            
                                                                                            //var CornerNodesSubdomains = new Dictionary<int, int[]> (CornerNodesData.Length); //nodeID, subdomains opou anhkei
            for (int i1 = 0; i1 < CornerNodesData.Length; i1++)
            {
                int h1 = CornerNodesData[i1][0]; int h2 = CornerNodesData[i1][1]; int h3 = CornerNodesData[i1][2];
                int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                double nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                double nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                double nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                double[] coordinates = new double[3] { nodeCoordX, nodeCoordY, nodeCoordZ };
                CornerNodesIds.Add(nodeID, coordinates);
                CornerNodesIdAndsubdomains.Add(nodeID, model.NodesDictionary[nodeID].SubdomainsDictionary.Keys.ToArray());
            }
            #endregion 

            //A.6
            double load_value = 1;
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[CornerNodesIds.ElementAt(0).Key],
                DOF = StructuralDof.TranslationZ,
                Amount = 1 * load_value
            };
            model.Loads.Add(load1);

            //A.7
            //DefineAppropriateConstraintsForBoundaryNodes(model, boundaryNodes);
            int[][] paktwsiNodesData = new int[4][]; //arithmos corner nodes,  h1 h2 h3 data (afairoume 1 apo ta pragmatika)
            int thesi = 0;
            int j1 = 0;
            for (int i2 = 0; i2 < 2; i2++)
            {
                for (int i3 = 0; i3 < 2; i3++)
                {
                    paktwsiNodesData[thesi] = new int[3] { j1, i2 * 1, i3 * 1 };
                    thesi++;
                }
            }
            Dictionary<int, Node> paktwshAristeraNodes = new Dictionary<int, Node>();
            for (int i1 = 0; i1 < paktwsiNodesData.Length; i1++)
            {
                int h1 = paktwsiNodesData[i1][0]; int h2 = paktwsiNodesData[i1][1]; int h3 = paktwsiNodesData[i1][2];
                int nodeID = renumbering.GetNewNodeNumbering(FEMMeshBuilder.Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                paktwshAristeraNodes.Add(nodeID, model.NodesDictionary[nodeID]);


            }
            DefineAppropriateConstraintsForBoundaryNodes(model, paktwshAristeraNodes);



            #region define solver from Quad4LinearCantilever %81 and IntegrationElastic... %37 tests.
            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            // Request output
            //int[] outputPositions = new int[model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs];
            //for( int i1 = 0;  i1 < model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs ; i1++ )
            //{

            //}
            //childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { 0 });

            // Run the anlaysis - part a
            parentAnalyzer.Initialize();

            //A.8            
            printElementStiffnessAndData(model);            


            // Run the anlaysis - part b
            parentAnalyzer.Solve();
            #endregion

            var globalU = solver.LinearSystems[1].Solution.CopyToArray();// fetiSolver.GatherGlobalDisplacements(sudomainDisplacements);
            double[] uc = new double[3 * CornerNodesIds.Count()];

            int node_counter = 0;
            foreach (int nodeId in CornerNodesIds.Keys)
            {
                //StructuralDof[] dofs = new StructuralDof[5] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ, StructuralDof.RotationX, StructuralDof.RotationY };
                StructuralDof[] dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

                Node node = model.NodesDictionary[nodeId];
                for (int i1 = 0; i1 < 3; i1++)
                {
                    int globalDof = model.GlobalDofOrdering.GlobalFreeDofs[node, dofs[i1]];
                    uc[3 * node_counter + i1] = globalU[globalDof];

                }
                node_counter++;
            }

            return (model, uc);
        }

        public static void printElementStiffnessAndData(Model model)
        {
            var dofOrdering = model.SubdomainsDictionary[0].FreeDofOrdering;
            var matrixProvider = new ElementStructuralStiffnessProvider();
            
            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}Stiffness.txt";
            string print_path_gen2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}elementDofIndices.txt";
            string print_path_gen3 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Element{0}subdomainDofIndices.txt";
            foreach (Element element in model.ElementsDictionary.Values)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                var elementMatrix = matrixProvider.Matrix(element).CopytoArray2D();
                string counterData = element.ID.ToString();
                string print_path = string.Format(print_path_gen, counterData);
                var writer = new Array2DWriter();
                writer.WriteToFile(elementMatrix, print_path, false);

                //Aray1Dwriter
                string print_path2 = string.Format(print_path_gen2, counterData);
                var writer2 = new MatlabWriter();
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(elementDofIndices)), print_path2, false);

                string print_path3 = string.Format(print_path_gen3, counterData);
                writer2.WriteToFile(Vector.CreateFromArray(Cnvrt(subdomainDofIndices)), print_path3, false);
            }
        }

        private static void DefineAppropriateConstraintsForBoundaryNodes(Model model, Dictionary<int, Node> boundaryNodes)
        {
            IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ImposeAppropriateConstraintsPerBoundaryNode(model, boundaryNode);
            }
        }

        public static double[] Cnvrt(int[] array)
        {
            double[] array2 = new double[array.Length];
            for(int i1=0; i1<array.Length; i1++)
            {
                array2[i1] = (double)array[i1];
            }
            return array2;
        }

    }
}
