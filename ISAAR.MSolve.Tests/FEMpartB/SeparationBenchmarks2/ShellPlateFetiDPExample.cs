using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Solvers.Direct;
using System;
using System.IO;
using System.Linq;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Tests
{
    public static class ShellPlateFetiDPExample
    {
        [Fact]
        public static void Solve()
        {
            var model = CreateModel();
        }
        public static Model CreateModel()
        {
            int subdiscrShell = 1;
            int discrShell = 2;

            
            int elem1 = subdiscrShell * discrShell;
            int elem2 = subdiscrShell * discrShell;
            var mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(2, 2, 2, subdiscrShell, discrShell);
            var gp = mpgp.Item2;
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;

            var coefficientsProvider = new SpectralRepresentation2DRandomField(0, 0, 0, 0.01, 2 * Math.PI / ((double)20), 20);
            coefficientsProvider.RandomVariables = null;
            double[] o_xsunol = FEMMeshBuilder.ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, gp.L1, gp.L2, gp.elem1, gp.elem2, gp.a1_shell, new double[] { 0, 0, 0 },
               new o_x_parameters() , coefficientsProvider);


            //model and node coordinates
            Model model = new Model();
            int totalNodes = o_xsunol.Length / 6;
            int[] newNumbering = new int[totalNodes];
            for(int i1 = 0; i1 < totalNodes; i1++) { newNumbering[i1] = i1 + 1; }
           
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int eswterikosNodeCounter = 0;
            var renumbering = new renumbering(newNumbering);
            for(int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter +  1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;

            //material
            var material2 = new ShellElasticMaterial3D()
            {
                YoungModulus = gp.E_shell,
                PoissonRatio = gp.ni_shell,
                ShearCorrectionCoefficientK = 5 / 6,
            };


            //perioxh orismou twn elements
            #region
            int elements = elem1 * elem2;
            int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
            int komvoi_8 = fdof_8 / 5;
            int[,] t_shell;
            t_shell = topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

            double[] Tk_vec = new double[8];
            double[][] VH = new double[8][];
            int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
            Element e2;
            int ElementID;

            int[] subdIDs = new int[discrShell * discrShell];
            for (int i1 = 0; i1 < discrShell*discrShell; i1++)
            {
                model.SubdomainsDictionary.Add(i1+1, new Subdomain(i1+1));
            }
            

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = gp.tk;
            }

            int eswterikosElementCounter = 0;
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter +  1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    //
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] ), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] )]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                int subdomainID = GetNodeSubdomainID(discrShell, discrShell, e2, gp.L1, gp.L2);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            #endregion

            model.ConnectDataStructures();

            List<Node> cornerNodeList = model.NodesDictionary.Values.Where(x => x.SubdomainsDictionary.Keys.Count == 4).ToList();
            


            return model;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 95, //150, // diastaseis
                L02 = 95, //150,
                L03 = 95, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr1 * subdiscr1,
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.20, //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.25, //0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.20, //0.05,// Gpa
                D_o_1 = 0.25, //0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
        {
            int elem;
            int[,] t_shell = new int[elements, 8];
            for (int nrow = 0; nrow < elem1; nrow++)
            {
                for (int nline = 0; nline < elem2; nline++)
                {
                    elem = (nrow + 1 - 1) * elem2 + nline + 1;//nrow+ 1 nline+1 einai zero based 
                    t_shell[elem - 1, -1 + 1] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 8] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 4] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;

                    t_shell[elem - 1, -1 + 5] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 2;
                    t_shell[elem - 1, -1 + 7] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 1;

                    t_shell[elem - 1, -1 + 2] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 6] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 3] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;
                }
            }
            return t_shell;

        }

        public static int GetNodeSubdomainID(int n1, int n2, Element element, double L01, double L02)
        {
            var x_meso = element.NodesDictionary.Values.Select(x => x.X).ToArray().Average()/8 +0.5*L01;
            var y_meso = element.NodesDictionary.Values.Select(x => x.Y).ToArray().Average()/8 +0.5*L02;

            int x_gemata = (int)Math.Truncate(x_meso/(L01 / n1));
            int y_gemata = (int)Math.Truncate(y_meso / (L02 / n2));

            int subd_ID = n2 * x_gemata + y_gemata + 1;
            return subd_ID;

        }

        private  static Dictionary<int, HashSet<INode>> DefineCornerNodesPerSubdomainAndOtherwise(List<Node>  CornerNodesList, Model model)
        {
            // a copy of this is used in modelCreatorInput class
            Dictionary<int, HashSet<INode>> cornerNodesList = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());
            
            foreach (Subdomain subdomain in model.EnumerateSubdomains())
            {
                cornerNodesList.Add(subdomain.ID, new HashSet<INode>());
            }

            foreach (Node node1 in CornerNodesList)
            {
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            //foreach (Subdomain subdomain in model.Subdomains)
            //{
            //    cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID]);
            //}

            return cornerNodesList;
        }
    }
}
