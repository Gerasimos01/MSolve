using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class HexaBuilder3CZM
    {
        public static void MakeHexaBuilding(Model model)
        {
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX , Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            // katw strwsh pou tha paktwthei

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = 1;
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
                T_o_3=57, D_o_3=5.7e-5, D_f_3=0.0098245610,
                T_o_1 = 57, D_o_1 = 5.7e-5, D_f_1 = 0.0098245610,
                n_curve=1.4,
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
                ElementType = new Shell8dispCopyPrint(material2, 3, 3, 3)
                {
                    oVn_i = VH,
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
                    endeixi_element_2=0,
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
            value_ext = 2.5 * 0.5 ;
            
            int[] points_with_negative_load;
            points_with_negative_load= new int [] { 1,3,6,8 };
            int[] points_with_positive_load;
            points_with_positive_load = new int[] { 2,4,5,7 };

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
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j+1]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j+1]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }


            // perioxh loads ews edw
        }

        public static void MakeHexaBuilding2(Model model)
        {
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            // katw strwsh pou tha paktwthei

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = -1;
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
                T_o_3 = 100*57,
                D_o_3 = 100*5.7e-5,
                D_f_3 = 100*0.0098245610,
                T_o_1 = 100*57,
                D_o_1 = 100*5.7e-5,
                D_f_1 = 100*0.0098245610,
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
                ElementType = new Shell8dispCopyPrint(material2, 3, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                }
            };
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);

            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh shell ews edw

            // eisagwgh tou cohesive element
            int[] coh_global_nodes;
            coh_global_nodes = new int[] { 3,8,6,1,5,7,4,2,11,16,14,9,13,15,12,10 };

            Element e2;
            e2 = new Element()
            {
                ID = 2,
                ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                    endeixi_element_2=1,
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
            value_ext = 100*2.5 * 0.5;

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
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }
            // perioxh loads ews edw
        }

        public static void MakeHexaBuilding3(Model model)
        {
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;
            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            double[][] VH = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = 1;
            }
            // perioxh gewmetrias ews edw

            // constraints
            nodeID = 3;
            for (int j = 0; j < 6; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                nodeID++;
            }
            //perioxh constraints ews edw

            // perioxh materials 
            ElasticMaterial3D material1 = new ElasticMaterial3D()
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
                ElementType = new Shell8dispCopyPrint(material1, 3, 3, 3)
                {
                    oVn_i = VH,
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

            // perioxh loads
            double value_ext;
            value_ext = 25 * 0.5;

            int[] points_with_negative_load;
            points_with_negative_load = new int[] { 1 };
            int[] points_with_positive_load;
            points_with_positive_load = new int[] { 2  };

            Load load1;
            Load load2;

            for (int j = 0; j < 1 ; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j ]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j ]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }
            // perioxh loads ews edw
        }

        public static void MakeHexaBuilding4(Model model)
        {
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }
            
            double[][] VH = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = -1;
            }
            // perioxh gewmetrias ews edw

            // constraints

            nodeID = 3;
            for (int j = 0; j < 6; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                nodeID++;
            }
            //perioxh constraints ews edw

            // perioxh materials 
            ElasticMaterial3D material1 = new ElasticMaterial3D()
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
                ElementType = new Shell8dispCopyPrint(material1, 3, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                }
            };
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);

            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh shell ews edw

            // perioxh loads
            double value_ext;
            value_ext = 25 * 0.5;

            int[] points_with_negative_load;
            points_with_negative_load = new int[] { 1 };
            int[] points_with_positive_load;
            points_with_positive_load = new int[] { 2 };

            Load load1;
            Load load2;

            for (int j = 0; j < 1; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j ]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }
            // perioxh loads ews edw
        }

        public static void MakeHexaBuilding5(Model model)
        {
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0.5;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            // katw strwsh pou tha paktwthei

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ + 0.5 * Tk });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ + 0.5 * Tk });
                nodeID++;
            }

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ + 0.5 * Tk });
                nodeID++;
            }

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = 1;
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
                T_o_3 = 100*57,
                D_o_3 = 100*5.7e-5,
                D_f_3 = 100*0.0098245610,
                T_o_1 = 100*57,
                D_o_1 = 100*5.7e-5,
                D_f_1 = 100*0.0098245610,
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
                ElementType = new Shell8dispCopyPrint(material2, 3, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                }
            };
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);

            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh shell ews edw

            // eisagwgh tou cohesive element
            int[] coh_global_nodes;
            coh_global_nodes = new int[] { 3, 8, 6, 1, 5, 7, 4, 2, 11, 16, 14, 9, 13, 15, 12, 10 };

            Element e2;
            e2 = new Element()
            {
                ID = 2,
                ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                    endeixi_element_2 = 1,
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
            value_ext = -100*2.5 * 0.5;

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
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }
            // perioxh loads ews edw
        }
    }
}
