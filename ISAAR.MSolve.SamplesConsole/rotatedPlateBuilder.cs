using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.SamplesConsole
{
    class rotatedPlateBuilder
    {
        public static double [] VectorTransformation(double[,] Qij, double [] vec)
        {
            double[] VecOut = new double[Qij.GetLength(0)];
            for (int i1 = 0; i1 < Qij.GetLength(0); i1++)
            {
                VecOut[i1] = 0;
                for (int j1 = 0; j1 < Qij.GetLength(1); j1++)
                {
                    VecOut[i1] += Qij[i1, j1] * vec[j1];
                }
            }
            return VecOut;
        }


        public static void ShellAndCohesiveShellPaktwsh(Model model) // proelefsi
        {
            //mhtrwo peristrofhs
            double rot_phi_2 = 0.7853981634; // (45 / 360) * 2 * Math.PI;
            double rot_phi_1 = 0;
            double e1_new_z = Math.Sin(rot_phi_2);
            double e1_new_y = Math.Sin(rot_phi_1) * Math.Cos(rot_phi_2);  //e1_new_xy = Math.Cos(rot_phi_2);
            double e1_new_x = Math.Cos(rot_phi_1) * Math.Cos(rot_phi_2);

            double e2_new_y = Math.Cos(rot_phi_1);
            double e2_new_x = -Math.Sin(rot_phi_1);
            double e2_new_z = 0;

            double[,] e_new = new double[3, 3] { { e1_new_x, e2_new_x, 0 }, { e1_new_y, e2_new_y, 0 }, { e1_new_z, e2_new_z, 0 } };
            e_new[0,2] = e_new[1, 0] * e_new[2, 1] - e_new[2, 0] * e_new[1, 1];
            e_new[1,2] = e_new[2, 0] * e_new[0, 1] - e_new[0, 0] * e_new[2, 1];
            e_new[2, 2] = e_new[0, 0] * e_new[1, 1] - e_new[1, 0] * e_new[0, 1];

            double[,] e_old = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

            double[,] Qij = new double[3, 3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    Qij[i1, j1] = e_old[0, i1] * e_new[0, j1] + e_old[1, i1] * e_new[1, j1] + e_old[2, i1] * e_new[2, j1];
                }
            }


            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                double[] orig_vec = new double[3] { startX, startY + l * 0.25, startZ };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);

                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                double [] orig_vec = new double [3] { startX, startY + l * 0.5, startZ };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);

                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                double[] orig_vec = new double[3] { startX, startY + l * 0.25, startZ };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            // katw strwsh pou tha paktwthei

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                double[] orig_vec = new double[3] { startX, startY + l * 0.25, startZ - 0.5 * Tk };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                double[] orig_vec = new double[3] { startX, startY + l * 0.5, startZ - 0.5 * Tk };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                double[] orig_vec = new double[3] { startX, startY + l * 0.25, startZ - 0.5 * Tk };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = upd_vec[0], Y = upd_vec[1], Z = upd_vec[2] });
                nodeID++;
            }

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                double[] orig_vec = new double[3] { 0,0,1 };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                VH[j] = new double[3];
                VH[j][0] = upd_vec[0];
                VH[j][1] = upd_vec[1];
                VH[j][2] = upd_vec[2];
            }
            // perioxh gewmetrias ews edw

            // constraints

            nodeID = 9;
            for (int j = 0; j < 8; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                nodeID++;
            }
            //perioxh constraints ews edw

            // perioxh materials 
            BenzeggaghKenaneCohMat material1 = new PreProcessor.Materials.BenzeggaghKenaneCohMat()
            {
                T_o_3 = 57, // New load case argurhs NR_shell_coh.m
                D_o_3 = 5.7e-5,
                D_f_3 = 0.0098245610,
                T_o_1 = 57,
                D_o_1 = 5.7e-5,
                D_f_1 = 0.0098245610,
                n_curve = 1.4,
            };

            ElasticMaterial3D material2 = new ElasticMaterial3D()
            {
                YoungModulus = 1353000,
                PoissonRatio = 0.3,
            };
            // perioxh materials ews edw




            //eisagwgh tou shell element
            double[] Tk_vec = new double[8];
            for (int j = 0; j < 8; j++)
            {
                Tk_vec[j] = Tk;
            }

            Element e1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Shell8dispCopyGet(material2, 3, 3, 3)
                {
                    oVn_i = VH, // mono epeidh einai ena element einai ok o orismos etsi kai de tha allaxei meta to VH
                    tk = Tk_vec,
                }
            };
            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh shell ews edw

            // eisagwgh tou cohesive element
            int[] coh_global_nodes;
            coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

            Element e2;
            e2 = new Element()
            {
                ID = 2,
                ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                    endeixi_element_2 = 0,
                }
            };

            for (int j = 0; j < 16; j++)
            {
                e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
            }

            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
            // eisagwgh cohesive ews edw

            // perioxh loads
            double value_ext;
            value_ext = 2 * 2.5 * 0.5;

            int[] points_with_negative_load;
            points_with_negative_load = new int[] { 1, 3, 6, 8 };
            int[] points_with_positive_load;
            points_with_positive_load = new int[] { 2, 4, 5, 7 };

            Load load1;
            Load load2;

            // LOADCASE '' orthi ''
            //for (int j = 0; j < 4; j++)
            //{
            //    load1 = new Load()
            //    {
            //        Node = model.NodesDictionary[points_with_negative_load[j]],
            //        DOF = DOFType.Z,
            //        Amount = -0.3333333 * value_ext,
            //    };
            //    model.Loads.Add(load1);

            //    load2 = new Load()
            //    {
            //        Node = model.NodesDictionary[points_with_positive_load[j]],
            //        DOF = DOFType.Z,
            //        Amount = 1.3333333 * value_ext,
            //    };
            //    model.Loads.Add(load2);
            //}

            // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
            for (int j = 0; j < 3; j++)
            {
                double[] orig_vec = new double[3] { 0, 0, -0.3333333 * value_ext };
                double[] upd_vec = VectorTransformation(Qij, orig_vec);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.X,
                    Amount = upd_vec[0],
                };
                model.Loads.Add(load1);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.Y,
                    Amount = upd_vec[1],
                };
                model.Loads.Add(load1);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = upd_vec[2],
                };
                model.Loads.Add(load1);

                orig_vec = new double[3] { 0, 0, 1.3333333 * value_ext };
                upd_vec = VectorTransformation(Qij, orig_vec);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.X,
                    Amount = upd_vec[0],
                };
                model.Loads.Add(load1);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.Y,
                    Amount = upd_vec[1],
                };
                model.Loads.Add(load1);
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = upd_vec[2],
                };
                model.Loads.Add(load1);


                //      //
                //load1 = new Load()
                //{
                //    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                //    DOF = DOFType.Z,
                //    Amount = -0.3333333 * value_ext,
                //};
                //model.Loads.Add(load1);

                //load2 = new Load()
                //{
                //    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                //    DOF = DOFType.Z,
                //    Amount = 1.3333333 * value_ext,
                //};
                //model.Loads.Add(load2);
            }


            // perioxh loads ews edw
        }
    }
}
