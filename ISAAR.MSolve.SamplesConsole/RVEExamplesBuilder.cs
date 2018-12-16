using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Materials;//using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.FEM.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;//using ISAAR.MSolve.Matrices.Interfaces;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.SamplesConsole
{
    class RVEExamplesBuilder
    {       
        public static void HexaElementsOnlyRVE(Model model,rveMatrixParameters mp,double [,] Dq)
        {
            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1+1; h1++)
            {
                for (int h2 = 0; h2 < hexa2+1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3+1; h3++)
                    {
                        nodeID = Topol_rve(h1+1, h2+1, h3+1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1+1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY = -0.5 * L02 + (h2+1 - 1) * (L02 / hexa2);  
                        nodeCoordZ = -0.5 * L03 + (h3+1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int [] globalNodeIDforlocalNode_i=new int [8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1+1 + (h2+1 - 1) * hexa1 + (h3+1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[0] = Topol_rve(h1+1 + 1, h2+1 + 1, h3+1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[1] = Topol_rve(h1+1, h2+1 + 1, h3+1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[2] = Topol_rve(h1+1, h2+1, h3+1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[3] = Topol_rve(h1+1 + 1, h2+1, h3+1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[4] = Topol_rve(h1+1 + 1, h2+1 + 1, h3+1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[5] = Topol_rve(h1+1, h2+1 + 1, h3+1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[6] = Topol_rve(h1+1, h2+1, h3+1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
                        globalNodeIDforlocalNode_i[7] = Topol_rve(h1+1 + 1, h2+1, h3+1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;
            Dq = new double[9, 3 * p_komvoi_rve];
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = f_komvoi_rve + j + 1;
                Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
                Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
                Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
                Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
                Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
                Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
                Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
                Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
                Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
            }


        }

        public static void HexaElementsOnlyRVEwithRenumbering(Model model, rveMatrixParameters mp, double[,] Dq,string renumberingVectorPath)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1 + 1; h1++)
            {
                for (int h2 = 0; h2 < hexa2 + 1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3 + 1; h3++)
                    {
                        nodeID = renumbering.GetNewNodeNumbering( Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                        nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int[] globalNodeIDforlocalNode_i = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;
            Dq = new double[9, 3 * p_komvoi_rve];
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
                Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
                Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
                Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
                Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
                Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
                Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
                Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
                Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
                Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
            }


        }

        public static void AddGrapheneSheet(Model model,grapheneSheetParameters gp,double [] ekk_xyz)
        {
            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm
    
            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            o_xsunol = ox_sunol_Builder_ekk(new_rows,new_lines, L1, L2,elem1,elem2, a1_shell,ekk_xyz);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0)/6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK=1,
            };

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
            int subdomainID = 1;
            
            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }
            
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
                    {
                        oVn_i = VH,
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalShellNode_i[j1]+PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1 ;
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = VH,
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue+arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue+arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1 ;
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = VH,
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2*arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2*arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        public static void AddGrapheneSheet2(Model model, grapheneSheetParameters gp, double[] ekk_xyz)
        {
            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            o_xsunol = ox_sunol_Builder_ekk(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK=1,
            };

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
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
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
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        public static void AddGrapheneSheet_with_o_x_parameters(Model model, grapheneSheetParameters gp, double[] ekk_xyz,o_x_parameters o_x_parameters)
        {
            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz,o_x_parameters);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK = 1,
            };

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
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
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
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = eswterikosNodeCounter + PreviousNodesNumberValue + 1;
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue, model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface,
                        model.NodesDictionary[midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        public static void AddGrapheneSheet_with_o_x_parameters_withRenumbering(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, string renumberingVectorPath)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK=1,
            };

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
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3))//ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
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
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3))
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        public static void AddGrapheneSheet_with_o_x_Input_withRenumbering(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, string renumberingVectorPath, string o_xsunol_input_path)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            //o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
            o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK=1,
            };

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
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        //TOUSE Gerasimos
        public static void AddGrapheneSheet_with_o_x_Input_from_MSOLVE_withRenumbering_from_MSOLVE(Model model, grapheneSheetParameters gp, renumbering renumbering,  double[] o_xsunol)
        {
            // Perioxh renumbering initialization 
            //renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            //double[] o_xsunol; //o_xsunol apo Msolve pia
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            //o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
            //o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path); //o_xsunol apo Msolve pia

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK=1,
            };

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
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
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
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node() { ID = NodeID, X = nodeCoordX, Y = nodeCoordY, Z = nodeCoordZ });
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }


        public static void AddConstraintsForNonSingularStiffnessMatrix(Model model, int hexa1, int hexa2, int hexa3)
        {
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            int nodeID;

            nodeID = Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka); 
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = Topol_rve(hexa1+1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = Topol_rve(1,hexa2+ 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);

        }

        public static void AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, string renumberingVectorPath)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            int nodeID;

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
            model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);

        }

        public static void AddConstraintsForFe2RveBounds(Model model, Dictionary<int, Node> boundaryNodes)
        {
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                int nodeID = boundaryNode.ID;
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
            }
        }


        public static void AddLoadsOnRveFromFile(Model model, int hexa1, int hexa2, int hexa3,string vectorpath)
        {
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            double[] Fxk_p_komvoi_rve;
            //Fxk_p_komvoi_rve = PrintUtilities.ReadVector(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
            Fxk_p_komvoi_rve = PrintUtilities.ReadVector(vectorpath);
            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;

            Load load_i;
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = f_komvoi_rve + j + 1;
                load_i = new Load()
                {
                    Node = model.NodesDictionary[komvos],
                    DOF = DOFType.X,
                    Amount = Fxk_p_komvoi_rve[3*(j)+0]
                };
                model.Loads.Add(load_i);

                load_i = new Load()
                {
                    Node = model.NodesDictionary[komvos],
                    DOF = DOFType.Y,
                    Amount = Fxk_p_komvoi_rve[3 * (j) + 1]
                };
                model.Loads.Add(load_i);

                load_i = new Load()
                {
                    Node = model.NodesDictionary[komvos],
                    DOF = DOFType.Z,
                    Amount = Fxk_p_komvoi_rve[3 * (j) + 2]
                };
                model.Loads.Add(load_i);
            }

            // Afairesh fortiwn apo tous desmevmenous vathmous eleftherias 
            int nodeID;
            int[] supportedDOFs=new int[9];
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            nodeID = Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[0] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[1] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
            supportedDOFs[2] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[3] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
            supportedDOFs[4] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[5] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[6] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[7] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[8] = 3 * (nodeID - f_komvoi_rve - 1) + 1;

            for (int j = 0; j < 9; j++)
            {
                model.Loads.RemoveAt(supportedDOFs[8 - j]); // afairoume apo pisw pros ta mpros gia na mh xalaei h thesh twn epomenwn pou tha afairethoun
            }

        }

        public static void AddLoadsOnRveFromFile_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, string vectorpath, string renumberingVectorPath)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            double[] Fxk_p_komvoi_rve;
            //Fxk_p_komvoi_rve = PrintUtilities.ReadVector(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
            Fxk_p_komvoi_rve = PrintUtilities.ReadVector(vectorpath);
            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;

            Load load_i;
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = f_komvoi_rve + j + 1;
                load_i = new Load()
                {
                    Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
                    DOF = DOFType.X,
                    Amount = Fxk_p_komvoi_rve[3 * (j) + 0]
                };
                model.Loads.Add(load_i);

                load_i = new Load()
                {
                    Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
                    DOF = DOFType.Y,
                    Amount = Fxk_p_komvoi_rve[3 * (j) + 1]
                };
                model.Loads.Add(load_i);

                load_i = new Load()
                {
                    Node = model.NodesDictionary[renumbering.GetNewNodeNumbering(komvos)],
                    DOF = DOFType.Z,
                    Amount = Fxk_p_komvoi_rve[3 * (j) + 2]
                };
                model.Loads.Add(load_i);
            }

            // Afairesh fortiwn apo tous desmevmenous vathmous eleftherias 
            int nodeID;
            int[] supportedDOFs = new int[9];
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            nodeID = Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[0] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[1] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
            supportedDOFs[2] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[3] = 3 * (nodeID - f_komvoi_rve - 1) + 1;
            supportedDOFs[4] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[5] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[6] = 3 * (nodeID - f_komvoi_rve - 1) + 2;

            nodeID = Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka);
            supportedDOFs[7] = 3 * (nodeID - f_komvoi_rve - 1) + 0;
            supportedDOFs[8] = 3 * (nodeID - f_komvoi_rve - 1) + 1;

            for (int j = 0; j < 9; j++)
            {
                model.Loads.RemoveAt(supportedDOFs[8 - j]); // afairoume apo pisw pros ta mpros gia na mh xalaei h thesh twn epomenwn pou tha afairethoun
            }

        }

        public static void AddXLoadsOnYZConstrainedOneElementRVE(Model model)
        {

            Load load_i;
            int[] LoadedNodes = new int[] { 2, 4, 6, 8 };
            for (int j = 0; j < 4; j++)
            {
                load_i = new Load()
                {
                    Node = model.NodesDictionary[LoadedNodes[j]],
                    DOF = DOFType.X,
                    Amount = 5.4218569154
                };
                model.Loads.Add(load_i);
            }
        }

        public static void AddConstraintsForYZConstraindeOneElementRVE(Model model)
        {
            int[] paktwmenaNodes = new int[] { 1, 3, 5, 7 };
            int[] NodesKulhshX = new int[] { 2, 4, 6, 8  };

            for (int j = 0; j < 4; j++)
            {
                model.NodesDictionary[paktwmenaNodes[j]].Constraints.Add(DOFType.X);
                model.NodesDictionary[paktwmenaNodes[j]].Constraints.Add(DOFType.Y);
                model.NodesDictionary[paktwmenaNodes[j]].Constraints.Add(DOFType.Z);
                model.NodesDictionary[NodesKulhshX[j]].Constraints.Add(DOFType.Y);
                model.NodesDictionary[NodesKulhshX[j]].Constraints.Add(DOFType.Z);
            }
                
        }



        private static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
        {
            int elem;
            int[,] t_shell = new int[elements, 8];
            for (int nrow = 0; nrow < elem1; nrow++)
            {
                for (int nline = 0; nline < elem2; nline++)
                {
                    elem = (nrow+1 - 1) * elem2 + nline+1;//nrow+ 1 nline+1 einai zero based 
                    t_shell[elem-1,-1+ 1]= (nrow+1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 3;
                    t_shell[elem-1,-1+ 8]= (nrow+1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 2;
                    t_shell[elem-1,-1+ 4]= (nrow+1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 1;

                    t_shell[elem-1,-1+ 5]= (nrow+1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline+1 - 1) * 1 + 2;
                    t_shell[elem-1,-1+ 7]= (nrow+1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline+1 - 1) * 1 + 1;

                    t_shell[elem-1,-1+ 2]= (nrow+1 - 1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 3;
                    t_shell[elem-1,-1+ 6]= (nrow+1 - 1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 2;
                    t_shell[elem-1,-1+ 3]= (nrow+1 - 1) * (3 * elem2 + 2) + (nline+1 - 1) * 2 + 1;
                }
            }
            return t_shell;

        }

        public static int Topol_rve(int h1, int h2, int h3, int hexa1, int hexa2, int hexa3, int kuvos, int endiam_plaka, int katw_plaka)
        {
            int arith;
            if (h3==1)
            { arith = h1 + (h2 - 1) * (hexa1 + 1) + kuvos; }
            else
            {
                if (h3== hexa3+1)
                { arith = hexa3 * (hexa1 + 1) * (hexa2 + 1) + h1 + (h2 - 1) * (hexa1 + 1); }
                else
                {
                    if (h2==1)
                    { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + h1; }
                    else
                    {
                        if (h2==hexa2+1)
                        { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + (hexa1 + 1) + 2 * (hexa2 - 1) + h1; }
                        else
                        {
                            if (h1==1)
                            { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 1; }
                            else
                            {
                                if(h1==hexa1+1)
                                { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 2; }
                                else
                                { arith = (h1 - 1) + (h2 - 2) * (hexa1 - 1) + (h3 - 2) * (hexa1 - 1) * (hexa2 - 1); }
                            }
                        }
                    }

                }
            }
            return arith;
        }

        public static double[] ox_sunol_Builder(int new_rows, int new_lines, double L1, double L2, int elem1, int elem2, double a1_shell)
        {
            double[] o_xsunol = new double[6 * new_lines * new_rows];
            int npoint;
            double ksi_1;
            double heta_1;
            double k1;
            double l1;
            double[] e_ksi_1 = new double[3];
            double[] e_heta_1 = new double[3];
            double[] e_3 = new double[3];
            double e_3norm;

            k1 = 2 * Math.PI / L1;
            l1 = 2 * Math.PI / L2;
            for (int nrow = 0; nrow < new_rows; nrow++)
            {
                for (int nline = 0; nline < new_lines; nline++)
                {
                    npoint = (nrow+1 - 1) * new_lines + nline+1; // nrow+1 kai nline +1 dioti einai zero based
                    ksi_1 = (nrow+1 - 1) * (L1 / (2 * elem1));
                    heta_1 = (nline+1 - 1) * (L2 / (2 * elem2));
                    o_xsunol[6 * (npoint - 1) + 1 - 1] = ksi_1 - 0.5 * L1;
                    o_xsunol[6 * (npoint - 1) + 2 - 1] = heta_1 - 0.5 * L2;
                    o_xsunol[6 * (npoint - 1) + 3 - 1] = a1_shell * Math.Sin(k1 * ksi_1) * Math.Sin(l1 * heta_1);

                    e_ksi_1 = new double[] { 1, 0, a1_shell * k1 * Math.Cos(k1 * ksi_1) * Math.Sin(l1 * heta_1) };
                    e_heta_1 = new double[] { 0, 1, a1_shell * l1 * Math.Sin(k1 * ksi_1) * Math.Cos(l1 * heta_1) };

                    cross(e_ksi_1, e_heta_1, e_3);
                    e_3norm = Math.Sqrt(Math.Pow(e_3[0], 2) + Math.Pow(e_3[1], 2) + Math.Pow(e_3[2], 2));

                    o_xsunol[(6 * (npoint - 1) + 3)] = e_3[0] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 4)] = e_3[1] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 5)] = e_3[2] / e_3norm;
                }
            }

            //grammes 71-72
            double[,] V_epil = new double[6 * (2 * elem1 + 1) * (2 * elem2 + 1), 6 * (2 * elem1 + 1) * (2 * elem2 + 1)];
            double[,] upomhtrwoV = new double[6 * (2 * elem2 + 1) + 6 * elem2 + 6, 6 * (2 * elem2 + 1) + 12 * elem2 + 6];
            // morfwsi upomhtrwouV           
            for (int j = 0; j < 6 * (2 * elem2 + 1); j++)
            {
                upomhtrwoV[j, j] = 1;
            }
            for (int i7 = 0; i7 < elem2; i7++)
            {
                for (int j = 0; j < 6; j++)
                {
                    upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * i7 + j), (6 * (2 * elem2 + 1) + 12 * i7 + j)] = 1;
                }
            }
            for (int j = 0; j < 6; j++)
            {
                upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * elem2 + j), (6 * (2 * elem2 + 1) + 12 * elem2 + j)] = 1;
            }

            // morfwsi mhtrwou V_epil
            for (int i7 = 0; i7 < elem1; i7++)
            {
                for (int j1 = 0; j1 < upomhtrwoV.GetLength(0); j1++)
                {
                    for (int j2 = 0; j2 < upomhtrwoV.GetLength(1); j2++)
                    {
                        V_epil[i7 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, i7 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j2] = upomhtrwoV[j1, j2];
                    }
                }
            }
            for (int j1 = 0; j1 < 6 * (2 * elem2 + 1); j1++)
            {
                V_epil[elem1 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, elem1 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j1] = 1;
            }

            // grammes 89-90
            double[] o_x_endiam = new double[o_xsunol.GetLength(0)];
            for (int j1 = 0; j1 < o_xsunol.GetLength(0); j1++)
            {
                o_x_endiam[j1] = o_xsunol[j1];
            }
            double[] o_x_endiam2;
            Matrix2D V_epil_mat = new Matrix2D(V_epil);
            o_x_endiam2 = (V_epil_mat * new Vector(o_x_endiam)).Data;

            // grammh 91
            o_xsunol = new double[elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1)];
            for (int j1 = 0; j1 < elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1); j1++)
            {
                o_xsunol[j1] = o_x_endiam2[j1];
            }

            return o_xsunol;

        }

        public static double[] ox_sunol_Builder_ekk(int new_rows, int new_lines, double L1, double L2, int elem1, int elem2, double a1_shell,double[] ekk_xyz)
        {
            double[] o_xsunol = new double[6 * new_lines * new_rows];
            int npoint;
            double ksi_1;
            double heta_1;
            double k1;
            double l1;
            double[] e_ksi_1 = new double[3];
            double[] e_heta_1 = new double[3];
            double[] e_3 = new double[3];
            double e_3norm;

            k1 = 2 * Math.PI / L1;
            l1 = 2 * Math.PI / L2;
            for (int nrow = 0; nrow < new_rows; nrow++)
            {
                for (int nline = 0; nline < new_lines; nline++)
                {
                    npoint = (nrow + 1 - 1) * new_lines + nline + 1; // nrow+1 kai nline +1 dioti einai zero based
                    ksi_1 = (nrow + 1 - 1) * (L1 / (2 * elem1));
                    heta_1 = (nline + 1 - 1) * (L2 / (2 * elem2));
                    o_xsunol[6 * (npoint - 1) + 1 - 1] = ksi_1 - 0.5 * L1+ekk_xyz[0];
                    o_xsunol[6 * (npoint - 1) + 2 - 1] = heta_1 - 0.5 * L2+ekk_xyz[1];
                    o_xsunol[6 * (npoint - 1) + 3 - 1] = a1_shell * Math.Sin(k1 * ksi_1) * Math.Sin(l1 * heta_1)+ekk_xyz[2];

                    e_ksi_1 = new double[] { 1, 0, a1_shell * k1 * Math.Cos(k1 * ksi_1) * Math.Sin(l1 * heta_1) };
                    e_heta_1 = new double[] { 0, 1, a1_shell * l1 * Math.Sin(k1 * ksi_1) * Math.Cos(l1 * heta_1) };

                    cross(e_ksi_1, e_heta_1, e_3);
                    e_3norm = Math.Sqrt(Math.Pow(e_3[0], 2) + Math.Pow(e_3[1], 2) + Math.Pow(e_3[2], 2));

                    o_xsunol[(6 * (npoint - 1) + 3)] = e_3[0] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 4)] = e_3[1] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 5)] = e_3[2] / e_3norm;
                }
            }

            //grammes 71-72
            double[,] V_epil = new double[6 * (2 * elem1 + 1) * (2 * elem2 + 1), 6 * (2 * elem1 + 1) * (2 * elem2 + 1)];
            double[,] upomhtrwoV = new double[6 * (2 * elem2 + 1) + 6 * elem2 + 6, 6 * (2 * elem2 + 1) + 12 * elem2 + 6];
            // morfwsi upomhtrwouV           
            for (int j = 0; j < 6 * (2 * elem2 + 1); j++)
            {
                upomhtrwoV[j, j] = 1;
            }
            for (int i7 = 0; i7 < elem2; i7++)
            {
                for (int j = 0; j < 6; j++)
                {
                    upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * i7 + j), (6 * (2 * elem2 + 1) + 12 * i7 + j)] = 1;
                }
            }
            for (int j = 0; j < 6; j++)
            {
                upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * elem2 + j), (6 * (2 * elem2 + 1) + 12 * elem2 + j)] = 1;
            }

            // morfwsi mhtrwou V_epil
            for (int i7 = 0; i7 < elem1; i7++)
            {
                for (int j1 = 0; j1 < upomhtrwoV.GetLength(0); j1++)
                {
                    for (int j2 = 0; j2 < upomhtrwoV.GetLength(1); j2++)
                    {
                        V_epil[i7 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, i7 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j2] = upomhtrwoV[j1, j2];
                    }
                }
            }
            for (int j1 = 0; j1 < 6 * (2 * elem2 + 1); j1++)
            {
                V_epil[elem1 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, elem1 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j1] = 1;
            }

            // grammes 89-90
            double[] o_x_endiam = new double[o_xsunol.GetLength(0)];
            for (int j1 = 0; j1 < o_xsunol.GetLength(0); j1++)
            {
                o_x_endiam[j1] = o_xsunol[j1];
            }
            double[] o_x_endiam2;
            Matrix2D V_epil_mat = new Matrix2D(V_epil);
            o_x_endiam2 = (V_epil_mat * new Vector(o_x_endiam)).Data;

            // grammh 91
            o_xsunol = new double[elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1)];
            for (int j1 = 0; j1 < elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1); j1++)
            {
                o_xsunol[j1] = o_x_endiam2[j1];
            }

            return o_xsunol;

        }

        // prosthiki print
        static int ox_sunol_counter = 0;
        static string string2 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\ox_sunol_MOSLVE_{0}.txt";
        // prosthiki print ews edw

        public static double[] ox_sunol_Builder_ekk_with_o_x_parameters(int new_rows, int new_lines, double L1, double L2, int elem1, int elem2, double a1_shell, double[] ekk_xyz,o_x_parameters o_x_paramters)
        {
            double[] o_xsunol = new double[6 * new_lines * new_rows];
            int npoint;
            double ksi_1;
            double heta_1;
            double k1;
            double l1;
            double[] e_ksi_1 = new double[3];
            double[] e_heta_1 = new double[3];
            double[] e_3 = new double[3];
            double e_3norm;

            k1 = 2 * Math.PI / L1;
            l1 = 2 * Math.PI / L2;
            for (int nrow = 0; nrow < new_rows; nrow++)
            {
                for (int nline = 0; nline < new_lines; nline++)
                {
                    npoint = (nrow + 1 - 1) * new_lines + nline + 1; // nrow+1 kai nline +1 dioti einai zero based
                    ksi_1 = (nrow + 1 - 1) * (L1 / (2 * elem1));
                    heta_1 = (nline + 1 - 1) * (L2 / (2 * elem2));
                    o_xsunol[6 * (npoint - 1) + 1 - 1] = ksi_1 - 0.5 * L1 + ekk_xyz[0];
                    o_xsunol[6 * (npoint - 1) + 2 - 1] = heta_1 - 0.5 * L2 + ekk_xyz[1];
                    o_xsunol[6 * (npoint - 1) + 3 - 1] = a1_shell * Math.Sin(k1 * ksi_1) * Math.Sin(l1 * heta_1) + ekk_xyz[2];

                    e_ksi_1 = new double[] { 1, 0, a1_shell * k1 * Math.Cos(k1 * ksi_1) * Math.Sin(l1 * heta_1) };
                    e_heta_1 = new double[] { 0, 1, a1_shell * l1 * Math.Sin(k1 * ksi_1) * Math.Cos(l1 * heta_1) };

                    cross(e_ksi_1, e_heta_1, e_3);
                    e_3norm = Math.Sqrt(Math.Pow(e_3[0], 2) + Math.Pow(e_3[1], 2) + Math.Pow(e_3[2], 2));

                    o_xsunol[(6 * (npoint - 1) + 3)] = e_3[0] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 4)] = e_3[1] / e_3norm;
                    o_xsunol[(6 * (npoint - 1) + 5)] = e_3[2] / e_3norm;
                }
            }

            //grammes 71-72
            double[,] V_epil = new double[6 * (2 * elem1 + 1) * (2 * elem2 + 1), 6 * (2 * elem1 + 1) * (2 * elem2 + 1)];
            double[,] upomhtrwoV = new double[6 * (2 * elem2 + 1) + 6 * elem2 + 6, 6 * (2 * elem2 + 1) + 12 * elem2 + 6];
            // morfwsi upomhtrwouV           
            for (int j = 0; j < 6 * (2 * elem2 + 1); j++)
            {
                upomhtrwoV[j, j] = 1;
            }
            for (int i7 = 0; i7 < elem2; i7++)
            {
                for (int j = 0; j < 6; j++)
                {
                    upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * i7 + j), (6 * (2 * elem2 + 1) + 12 * i7 + j)] = 1;
                }
            }
            for (int j = 0; j < 6; j++)
            {
                upomhtrwoV[(6 * (2 * elem2 + 1) + 6 * elem2 + j), (6 * (2 * elem2 + 1) + 12 * elem2 + j)] = 1;
            }

            // morfwsi mhtrwou V_epil
            for (int i7 = 0; i7 < elem1; i7++)
            {
                for (int j1 = 0; j1 < upomhtrwoV.GetLength(0); j1++)
                {
                    for (int j2 = 0; j2 < upomhtrwoV.GetLength(1); j2++)
                    {
                        V_epil[i7 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, i7 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j2] = upomhtrwoV[j1, j2];
                    }
                }
            }
            for (int j1 = 0; j1 < 6 * (2 * elem2 + 1); j1++)
            {
                V_epil[elem1 * (6 * (2 * elem2 + 1) + 6 * elem2 + 6) + j1, elem1 * (6 * (2 * elem2 + 1) + 12 * elem2 + 6) + j1] = 1;
            }

            // grammes 89-90
            double[] o_x_endiam = new double[o_xsunol.GetLength(0)];
            for (int j1 = 0; j1 < o_xsunol.GetLength(0); j1++)
            {
                o_x_endiam[j1] = o_xsunol[j1];
            }
            double[] o_x_endiam2;
            Matrix2D V_epil_mat = new Matrix2D(V_epil);
            o_x_endiam2 = (V_epil_mat * new Vector(o_x_endiam)).Data;

            // grammh 91
            o_xsunol = new double[elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1)];
            for (int j1 = 0; j1 < elem1 * (6 * (2 * elem2 + 1) + 6 * (elem2 + 1)) + 6 * (2 * elem2 + 1); j1++)
            {
                o_xsunol[j1] = o_x_endiam2[j1];
            }

            // gia print 
            ox_sunol_counter += 1;
            string counter_data = ox_sunol_counter.ToString();
            string path = string.Format(string2, counter_data);
            Vector o_xsunol_print= new Vector(o_xsunol);
            o_xsunol_print.WriteToFile(path);
            //gia print ews edw

            return o_xsunol;

        }


        public static void cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }
        public static double [,] eye(int dimension)
        {
            double[,] eyeMatrix = new double[dimension, dimension];
            for (int l = 0; l < dimension; l++)
            {
                eyeMatrix[l, l] = 1;
            }
            return eyeMatrix;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetOriginalRveCholExampleParameters()
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 150, // diastaseis
                L02 = 150,
                L03 = 40,
                hexa1 = 8,// diakritopoihsh
                hexa2 = 8,
                hexa3 = 5
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem2 = 8,
                elem1 = 8,
                L1 = 105,// nm
                L2 = 110,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 3, //6, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);

            return gpmp;

        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetOneElementRveCheckExampleParameters()
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 150, // diastaseis
                L02 = 150,
                L03 = 40,
                hexa1 = 2,// diakritopoihsh
                hexa2 = 2,
                hexa3 = 2
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = 2,
                elem2 = 2,
                L1 = 105,// nm
                L2 = 110,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk =  0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> Get2GrapheneSheetsRveExampleParameters()
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 150, // diastaseis
                L02 = 150,
                L03 = 40,
                hexa1 = 2,// diakritopoihsh
                hexa2 = 2,
                hexa3 = 2
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = 1, // DIORTHOSI 2 graohene sheets
                elem2 = 2,
                L1 = 0.5*105,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 110,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceRveExampleParameters(int subdiscr1,int discr1,int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 100, //150, // diastaseis
                L02 = 100, //150,
                L03 = 100, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr3
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
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }


        public static void OriginalRVECholExample(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            mpgp =RVEExamplesBuilder.GetOriginalRveCholExampleParameters();
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[] ekk_xyz = new double[] { 0, 0, 0 };
            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();
            RVEExamplesBuilder.AddGrapheneSheet2(model, gp,ekk_xyz);
            int shellElementsNumber = (model.ElementsDictionary.Count() - hexaElementsNumber) / 3;
            // model: add loads NA ALLAXTHEI TO PATH
            //RVEExamplesBuilder.AddLoadsOnRveFromFile(model,mp.hexa1,mp.hexa2,mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol\Fxk_p_komvoi_rve.txt");
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\Fxk_p_komvoi_rve.txt");
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);

            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value), model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + hexaElementsNumber + 1)).Select(kv => kv.Value));

        }

        public static void OneElementRVECheckExampleConstrained(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            mpgp = RVEExamplesBuilder.GetOneElementRveCheckExampleParameters();
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[] ekk_xyz = new double[] { 0, 0, 0 };
            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp,Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();
            RVEExamplesBuilder.AddGrapheneSheet(model, gp,ekk_xyz); // kai edw thelei graphene sheet2 an kai den paizei rolo 
            int shellElementsNumber = (model.ElementsDictionary.Count()-hexaElementsNumber)/3;
            // model: add loads
            //RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
            RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            //RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa2);
            RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber+1) ).Select(kv => kv.Value), model.ElementsDictionary.Where(x => (x.Key>=shellElementsNumber+hexaElementsNumber+1)).Select(kv => kv.Value));

        }

        public static void FewElementsRVECheckExample(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            mpgp = RVEExamplesBuilder.GetOneElementRveCheckExampleParameters();
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[] ekk_xyz = new double[] { 0, 0, 0.25*40 };
            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();
            RVEExamplesBuilder.AddGrapheneSheet2(model, gp,ekk_xyz);
            int shellElementsNumber = (model.ElementsDictionary.Count() - hexaElementsNumber) / 3;
            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value), model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + hexaElementsNumber + 1)).Select(kv => kv.Value));
        }

        //public static void FewElementsRVECheckExample_embedder_2(Model model)
        //{
        //    double[,] Dq = new double[1, 1];
        //    Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
        //    rveMatrixParameters mp;
        //    grapheneSheetParameters gp;
        //    mpgp = RVEExamplesBuilder.GetOneElementRveCheckExampleParameters();
        //    mp = mpgp.Item1;
        //    gp = mpgp.Item2;
        //    double[] ekk_xyz = new double[] { 0, 0, 0.25 * 40 };
        //    RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
        //    int hexaElementsNumber = model.ElementsDictionary.Count();
        //    RVEExamplesBuilder.AddGrapheneSheet2(model, gp, ekk_xyz);
        //    int shellElementsNumber = (model.ElementsDictionary.Count() - hexaElementsNumber) / 3;
        //    // model: add loads
        //    RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
        //    //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
        //    // model: add constraints
        //    RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
        //    //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);
        //    var embeddedGrouping = new EmbeddedCohesiveGrouping_v2(model, model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value), model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + hexaElementsNumber + 1)).Select(kv => kv.Value));
        //}

        public static void FewElementsRVECheckExample2GrapheneSheets(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            mpgp = RVEExamplesBuilder.Get2GrapheneSheetsRveExampleParameters();
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { -0.25*105, 0, 0.25 * 40 }, new double[] { 0.25*105, 0, 0.25 * 40 } };
            
            int graphene_sheets_number = 2;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count()+1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();
                
            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\elegxos_alalgwn_fe2_tax_me1_arxiko_chol_dixws_me1_OneElementRVECheckExample\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference1RVEExample10000(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF1_10000\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference1RVEExample100_000(Model model) // copy apo 10.000 kai dedomena apo to 100_000 me renumbering
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 8;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 28;
            int subdiscr1_shell = 18;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample50000(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 6;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 23;
            int subdiscr1_shell = 14;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_50000\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample50000withRenumbering(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_50000_renu_new\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_50000_renu_new\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 6;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 23;
            int subdiscr1_shell = 14;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j],renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample100_000withRenumbering(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 8;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 28;
            int subdiscr1_shell = 18;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample100_000withRenumbering_me_dev(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 8;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 28;
            int subdiscr1_shell = 18;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample500_000withRenumbering_me_dev(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_500_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_500_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 14;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 50;
            int subdiscr1_shell = 26;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample_1_000_000withRenumbering_me_dev(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_1_000_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_me_develop\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_1_000_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 17;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 68;
            int subdiscr1_shell = 36;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumbering(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10000_renu\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10000_renu\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumbering_2GrSh(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi_multple_graphene\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_2_graphene_renu_new_multiple_algorithms\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi_multple_graphene\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_2_graphene_renu_new_multiple_algorithms\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 4;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            gp.elem2 = 8; gp.L1 = 0.5 * 50;
            double[][] ekk_xyz = new double[2][] { new double[] { -12.5, 0, 0 }, new double[] { 12.5, 0, 0 } };

            int graphene_sheets_number = 2;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumbering_1GrSh(Model model)
        {
            // dhmiourghthike kata to Reference2RVEExample100_000withRenumbering
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi_multple_graphene\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_1_graphene_antistoixo_twn_duo_renu_new_multiple_algorithms\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi_multple_graphene\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_1_graphene_antistoixo_twn_duo_renu_new_multiple_algorithms\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 8;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample100_000withRenumbering_mono_hexa(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_mono_hexa\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_mono_hexa\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 8;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 28;
            int subdiscr1_shell = 18;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample500_000withRenumbering_mono_hexa(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_mono_hexa\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_mono_hexa\Fxk_p_komvoi_rve.txt";
            //int subdiscr1 = 4;
            //int discr1 = 4;
            //// int discr2 dn xrhsimopoieitai
            //int discr3 = 10;
            //int subdiscr1_shell = 7;
            //int discr1_shell = 1;

            int subdiscr1 = 14;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 50;
            int subdiscr1_shell = 26;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10_000withRenumbering_mono_hexa(Model model)
        {
            // COPY APO: Reference2RVEExample100_000withRenumbering_mono_hexa
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            //C:\Users\cluster 5\Desktop\Gerasimos\REFERENCE_Examples_me_develop\10_000_mono_hexa\REF_new_total_numbering.txt einai link sto PC LAB
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_mono_hexa\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_mono_hexa\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;

            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumberingwithInput(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_develop\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_develop\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen= @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_develop\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j+1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference2RVEExample10000withRenumberingwithInput_develop2(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            string renumbering_vector_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_twn_10_000\REF_new_total_numbering.txt";
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_twn_10_000\Fxk_p_komvoi_rve.txt";
            string o_xsunol_input_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_100_000_renu_new_multiple_algorithms_check_twn_10_000\o_xsunol_gs_{0}.txt";
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVEwithRenumbering(model, mp, Dq, renumbering_vector_path);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                string file_no = (j + 1).ToString();
                string ox_sunol_input_path = string.Format(o_xsunol_input_path_gen, file_no);
                //RVEExamplesBuilder.AddGrapheneSheet_with_o_x_Input_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path, ox_sunol_input_path);
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters_withRenumbering(model, gp, ekk_xyz[j], model_o_x_parameteroi[j], renumbering_vector_path);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, Fxk_p_komvoi_rve_path, renumbering_vector_path);
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering(model, mp.hexa1, mp.hexa2, mp.hexa3, renumbering_vector_path);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        

    }

    public class rveMatrixParameters
    {
        public double E_disp { get; set; }
        public double ni_disp { get; set; }
        public double L01 { get; set; }
        public double L02 { get; set; }
        public double L03 { get; set; }
        public int hexa1 { get; set; }
        public int hexa2 { get; set; }
        public int hexa3 { get; set; }

        public rveMatrixParameters()
        {

        }
        public rveMatrixParameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
        {
            this.E_disp = E_disp;
            this.ni_disp = ni_disp;
            this.L01 = L01;
            this.L02 = L02;
            this.L03 = L03;
            this.hexa1 = hexa1;
            this.hexa2 = hexa2;
            this.hexa3 = hexa3;
        }
    }

    public class o_x_parameters
    {
        //public double E_disp { get; set; }
        //public double ni_disp { get; set; }
        //public double L01 { get; set; }
        //public double L02 { get; set; }
        //public double L03 { get; set; }
        //public int hexa1 { get; set; }
        //public int hexa2 { get; set; }
        //public int hexa3 { get; set; }

        public o_x_parameters()
        {

        }
        public o_x_parameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
        {
            //this.E_disp = E_disp;
            //this.ni_disp = ni_disp;
            //this.L01 = L01;
            //this.L02 = L02;
            //this.L03 = L03;
            //this.hexa1 = hexa1;
            //this.hexa2 = hexa2;
            //this.hexa3 = hexa3;
        }
    }

    public class grapheneSheetParameters
    {
        // parametroi shell
        public double E_shell; // GPa = 1000Mpa = 1000N / mm2
        public double ni_shell; // stathera poisson
        public int elem1;
        public int elem2;
        public double L1;// nm
        public double L2;// nm
        public double L3; // nm
        public double a1_shell; // nm
        public double tk;  // 0.0125016478913782nm
                           //parametroi cohesive epifaneias
        public double T_o_3;// Gpa = 1000Mpa = 1000N / mm2
        public double D_o_3; // nm
        public double D_f_3; // nm
        public double T_o_1;// Gpa
        public double D_o_1; // nm
        public double D_f_1; // nm
        public double n_curve = 1.4;

        public grapheneSheetParameters()
        {

        }
        public grapheneSheetParameters(double E_shell, double ni_shell, int elem1, int elem2, double L1, double L2, double L3, double a1_shell, double tk,
            double T_o_3, double D_o_3, double D_f_3, double T_o_1, double D_o_1, double D_f_1, double n_curve)
        {
            this.E_shell = E_shell; // GPa = 1000Mpa = 1000N / mm2
            this.ni_shell = ni_shell; // stathera poisson
            this.elem1 = elem1;
            this.elem2 = elem2;
            this.L1 = L1;// nm
            this.L2 = L2;// nm
            this.L3 = L3; // nm
            this.a1_shell = a1_shell; // nm
            this.tk = tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            this.T_o_3 = T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            this.D_o_3 = D_o_3; // nm
            this.D_f_3 = D_f_3; // nm

            this.T_o_1 = T_o_1;// Gpa
            this.D_o_1 = D_o_1; // nm
            this.D_f_1 = D_f_1; // nm

            this.n_curve = n_curve;
        }
    }

    public class renumbering
    {
        public int [] sunol_nodes_numbering { get; set; }

        public renumbering()
        {

        }
        public renumbering(int[] sunol_nodes_numbering)
        {
            this.sunol_nodes_numbering = sunol_nodes_numbering;
        }

        public int GetNewNodeNumbering (int initial_node_number)
        {
            return sunol_nodes_numbering[initial_node_number-1];
        }

    }
    
}
