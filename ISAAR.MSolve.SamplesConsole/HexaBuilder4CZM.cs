using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class HexaBuilder4CZM
    {
        public static void MakeHexaBuilding(Model model)
        {
            // gewmetria
            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 1;

            for (int n = 0; n < 2; n++)
            {
                for (int m = 0; m < 2; m++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + m * 0.5, Y = startY + l * 0.5, Z = startZ-n*0.5 });
                        nodeID++;
                    }
                }
            }
            // katw strwsh pou tha paktwthei
            startX = 0;
            startY = 0;
            startZ = 0.5;
            for (int m = 0; m < 2; m++)
            {
                for (int l = 0; l < 2; l++)
                {
                    model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + m * 0.5, Y = startY + l * 0.5, Z = startZ });
                    nodeID++;
                }
            }
            // perioxh gewmetrias ews edw

            // constraints
            nodeID = 1;
            for (int j = 0; j < 6; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                nodeID++;
            }
            nodeID = 9;
            for (int j = 0; j < 4; j++)
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
                YoungModulus = 135300,
                PoissonRatio = 0.3,
            };
            // perioxh materials ews edw




            //eisagwgh tou hexa element
            Element e1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Hexa8NL(material2,3,3,3) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
            };
            int[] hexa_global_nodes;
            hexa_global_nodes = new int[] { 4, 2, 1, 3, 8, 6, 5, 7 };
            for (int j = 0; j < 8; j++)
            {
                e1.NodesDictionary.Add(hexa_global_nodes[j], model.NodesDictionary[hexa_global_nodes[j]]);
            }
            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh hexa ews edw

            // eisagwgh tou cohesive element
            int[] coh_global_nodes;
            coh_global_nodes = new int[] { 8,6,5,7,12,10,9,11 };

            Element e2;
            e2 = new Element()
            {
                ID = 2,
                ElementType = new cohesive8node(material1, 3, 3)
            };

            for (int j = 0; j < 8; j++)
            {
                e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
            }
            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
            // eisagwgh cohesive ews edw



            // perioxh loads
            double value_ext;
            value_ext = 0.0807692307692308*104 ;
            
            int[] points_with_load;
            points_with_load= new int [] { 7,8 };



            Load load1;
            Load load2;
            for (int j = 0; j < 2; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_load[j]],
                    DOF = DOFType.X,
                    Amount = value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_load[j]],
                    DOF = DOFType.Z,
                    Amount = 0.2 * value_ext,
                };
                model.Loads.Add(load2);
            }
            // perioxh loads ews edw
        }

        public static void MakeHexaBuilding2(Model model)
        {
            //  var embeddedGrouping = new EmbeddedGrouping(Hexa8MeshGenerator.model, Hexa8MeshGenerator.model.ElementsDictionary.Where(x => x.Key < 5).Select(kv => kv.Value), Hexa8MeshGenerator.model.ElementsDictionary.Where(x => x.Key >= 5).Select(kv => kv.Value));


            //int[] nodesPerDirection = new int[] { 5, 2, 2 };
            //double[] distancePerDirection = new double[] { 1, 1, 1 };
            //int[] elementsPerSubdomainSide = new int[] { 4, 1, 1 };

            //Surface pinSurface = Surface.XYFront;
            //Surface poreSurface = Surface.None;
            //Surface loadSurface = Surface.None;
            //Surface massLoadSurface = Surface.XZDown;
            //DOFType loadDirection = DOFType.Y;
            //MassAccelerationHistoryLoad massLoadAmount = new EmptyMassAccelerationHistoryLoad();
            //double loadAmount = 0;
            //Hexa8MeshGenerator.BuildEmbeddedMesh(nodesPerDirection, distancePerDirection, elementsPerSubdomainSide,
            //    pinSurface, poreSurface, loadSurface, massLoadSurface, loadDirection, loadAmount, massLoadAmount);

            //Hexa8MeshGenerator.model.NodesDictionary.Add(21, new Node() { ID = 21, X = 0, Y = 0.6, Z = 0.5 });
            //Hexa8MeshGenerator.model.NodesDictionary.Add(22, new Node() { ID = 22, X = 1, Y = 0.6, Z = 0.5 });
            //Hexa8MeshGenerator.model.NodesDictionary.Add(23, new Node() { ID = 23, X = 2, Y = 0.6, Z = 0.5 });
            //Hexa8MeshGenerator.model.NodesDictionary.Add(24, new Node() { ID = 24, X = 3, Y = 0.6, Z = 0.5 });
            //Hexa8MeshGenerator.model.NodesDictionary.Add(25, new Node() { ID = 25, X = 4, Y = 0.6, Z = 0.5 });

            //var e = new Element() { ID = 5, ElementType = new Beam3D(new ElasticMaterial3D() { YoungModulus = 2.1E+11, PoissonRatio = 0.32 }) { SectionArea = 3.141519E-4, MomentOfInertiaPolar = 1.56947E-9, MomentOfInertiaY = 7.85398E-9, MomentOfInertiaZ = 7.85398E-9 } };
            //e.NodesDictionary.Add(21, Hexa8MeshGenerator.model.NodesDictionary[21]);
            //e.NodesDictionary.Add(22, Hexa8MeshGenerator.model.NodesDictionary[22]);
            //Hexa8MeshGenerator.model.ElementsDictionary.Add(e.ID, e);
            //Hexa8MeshGenerator.model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            //e = new Element() { ID = 6, ElementType = new Beam3D(new ElasticMaterial3D() { YoungModulus = 2.1E+11, PoissonRatio = 0.32 }) { SectionArea = 3.141519E-4, MomentOfInertiaPolar = 1.56947E-9, MomentOfInertiaY = 7.85398E-9, MomentOfInertiaZ = 7.85398E-9 } };
            //e.NodesDictionary.Add(22, Hexa8MeshGenerator.model.NodesDictionary[22]);
            //e.NodesDictionary.Add(23, Hexa8MeshGenerator.model.NodesDictionary[23]);
            //Hexa8MeshGenerator.model.ElementsDictionary.Add(e.ID, e);
            //Hexa8MeshGenerator.model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            //e = new Element() { ID = 7, ElementType = new Beam3D(new ElasticMaterial3D() { YoungModulus = 2.1E+11, PoissonRatio = 0.32 }) { SectionArea = 3.141519E-4, MomentOfInertiaPolar = 1.56947E-9, MomentOfInertiaY = 7.85398E-9, MomentOfInertiaZ = 7.85398E-9 } };
            //e.NodesDictionary.Add(23, Hexa8MeshGenerator.model.NodesDictionary[23]);
            //e.NodesDictionary.Add(24, Hexa8MeshGenerator.model.NodesDictionary[24]);
            //Hexa8MeshGenerator.model.ElementsDictionary.Add(e.ID, e);
            //Hexa8MeshGenerator.model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            //e = new Element() { ID = 8, ElementType = new Beam3D(new ElasticMaterial3D() { YoungModulus = 2.1E+11, PoissonRatio = 0.32 }) { SectionArea = 3.141519E-4, MomentOfInertiaPolar = 1.56947E-9, MomentOfInertiaY = 7.85398E-9, MomentOfInertiaZ = 7.85398E-9 } };
            //e.NodesDictionary.Add(24, Hexa8MeshGenerator.model.NodesDictionary[24]);
            //e.NodesDictionary.Add(25, Hexa8MeshGenerator.model.NodesDictionary[25]);
            //Hexa8MeshGenerator.model.ElementsDictionary.Add(e.ID, e);
            //Hexa8MeshGenerator.model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            //var embeddedGrouping = new EmbeddedGrouping(Hexa8MeshGenerator.model, Hexa8MeshGenerator.model.ElementsDictionary.Where(x => x.Key < 5).Select(kv => kv.Value), Hexa8MeshGenerator.model.ElementsDictionary.Where(x => x.Key >= 5).Select(kv => kv.Value));

            //Hexa8MeshGenerator.model.Loads.Add(new Load() { Node = Hexa8MeshGenerator.model.NodesDictionary[5], DOF = DOFType.Y, Amount = -25000 });
            //Hexa8MeshGenerator.model.Loads.Add(new Load() { Node = Hexa8MeshGenerator.model.NodesDictionary[10], DOF = DOFType.Y, Amount = -25000 });
            //Hexa8MeshGenerator.model.Loads.Add(new Load() { Node = Hexa8MeshGenerator.model.NodesDictionary[15], DOF = DOFType.Y, Amount = -25000 });
            //Hexa8MeshGenerator.model.Loads.Add(new Load() { Node = Hexa8MeshGenerator.model.NodesDictionary[20], DOF = DOFType.Y, Amount = -25000 });

            //Hexa8MeshGenerator.model.ConnectDataStructures();

            //SolverSkyline solver = new SolverSkyline(Hexa8MeshGenerator.model);
            //ProblemStructural provider = new ProblemStructural(Hexa8MeshGenerator.model, solver.SubdomainsDictionary);
            //LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);

            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            //parentAnalyzer.BuildMatrices();
            //parentAnalyzer.Initialize();
            //parentAnalyzer.Solve();

        }

        //public static void MakeHexaBuilding2(Model model)
        //{
        //    // gewmetria
        //    double Tk = 0.5;

        //    int nodeID = 1;

        //    double startX = 0;
        //    double startY = 0;
        //    double startZ = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
        //        nodeID++;
        //    }

        //    // katw strwsh pou tha paktwthei

        //    startX = 0;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.25;
        //    for (int l = 0; l < 2; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    startX = 0.5;
        //    for (int l = 0; l < 3; l++)
        //    {
        //        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
        //        nodeID++;
        //    }

        //    double[][] VH = new double[8][];

        //    for (int j = 0; j < 8; j++)
        //    {
        //        VH[j] = new double[3];
        //        VH[j][0] = 0;
        //        VH[j][1] = 0;
        //        VH[j][2] = -1;
        //    }
        //    // perioxh gewmetrias ews edw

        //    // constraints

        //    nodeID = 9;
        //    for (int j = 0; j < 8; j++)
        //    {
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
        //        model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
        //        nodeID++;
        //    }
        //    //perioxh constraints ews edw

        //    // perioxh materials 
        //    BenzeggaghKenaneCohMat material1 = new PreProcessor.Materials.BenzeggaghKenaneCohMat()
        //    {
        //        T_o_3 = 57,
        //        D_o_3 = 5.7e-5,
        //        D_f_3 = 0.0098245610,
        //        T_o_1 = 57,
        //        D_o_1 = 5.7e-5,
        //        D_f_1 = 0.0098245610,
        //        n_curve = 1.4,
        //    };

        //    ElasticMaterial3D material2 = new ElasticMaterial3D()
        //    {
        //        YoungModulus = 1353000,
        //        PoissonRatio = 0.3,
        //    };
        //    // perioxh materials ews edw




        //    //eisagwgh tou shell element
        //    double[] Tk_vec = new double[8];
        //    for (int j = 0; j < 8; j++)
        //    {
        //        Tk_vec[j] = Tk;
        //    }

        //    Element e1;
        //    e1 = new Element()
        //    {
        //        ID = 1,
        //        ElementType = new Shell8dispCopyPrint(material2, 3, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //        }
        //    };
        //    e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
        //    e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
        //    e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
        //    e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
        //    e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
        //    e1.NodesDictionary.Add(7, model.NodesDictionary[7]);
        //    e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
        //    e1.NodesDictionary.Add(2, model.NodesDictionary[2]);

        //    int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
        //    model.ElementsDictionary.Add(e1.ID, e1);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
        //    //eisagwgh shell ews edw

        //    // eisagwgh tou cohesive element
        //    int[] coh_global_nodes;
        //    coh_global_nodes = new int[] { 3,8,6,1,5,7,4,2,11,16,14,9,13,15,12,10 };

        //    Element e2;
        //    e2 = new Element()
        //    {
        //        ID = 2,
        //        ElementType = new cohesive_shell_to_hexaCopyGet(material1, 3, 3)
        //        {
        //            oVn_i = VH,
        //            tk = Tk_vec,
        //            endeixi_element_2=1,
        //        }
        //    };


        //    for (int j = 0; j < 16; j++)
        //    {
        //        e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
        //    }

        //    model.ElementsDictionary.Add(e2.ID, e2);
        //    model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
        //    // eisagwgh cohesive ews edw

        //    // perioxh loads
        //    double value_ext;
        //    value_ext = 2.5 * 0.5;

        //    int[] points_with_negative_load;
        //    points_with_negative_load = new int[] { 1, 3, 6, 8 };
        //    int[] points_with_positive_load;
        //    points_with_positive_load = new int[] { 2, 4, 5, 7 };

        //    Load load1;
        //    Load load2;

        //    // LOADCASE '' orthi ''
        //    //for (int j = 0; j < 4; j++)
        //    //{
        //    //    load1 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_negative_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = -0.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load1);

        //    //    load2 = new Load()
        //    //    {
        //    //        Node = model.NodesDictionary[points_with_positive_load[j]],
        //    //        DOF = DOFType.Z,
        //    //        Amount = 1.3333333 * value_ext,
        //    //    };
        //    //    model.Loads.Add(load2);
        //    //}

        //    // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
        //    for (int j = 0; j < 3; j++)
        //    {
        //        load1 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_negative_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = -0.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load1);

        //        load2 = new Load()
        //        {
        //            Node = model.NodesDictionary[points_with_positive_load[j + 1]],
        //            DOF = DOFType.Z,
        //            Amount = 1.3333333 * value_ext,
        //        };
        //        model.Loads.Add(load2);
        //    }


        //    // perioxh loads ews edw
        //}
    }
}
