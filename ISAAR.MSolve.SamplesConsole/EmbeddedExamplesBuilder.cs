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
    public static class EmbeddedExamplesBuilder
    {
        public static void HexaElementsOnly(Model model)
        {
            int startX = 0;
            int startY = 0;
            int startZ = 0;

            int nodeID = 1;
            for (int l = 0; l < 3; l++)
            {
                for (int k = 0; k < 2; k++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + j * 1, Y = startY + k * 1, Z = startZ + l * 1 });

                        nodeID++;
                    }
                }
            }
            nodeID = 1;
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                    model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                    nodeID++;
                }
            }
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };

            //eisagwgh enos element
            Element e1;
            int ID2 = 1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Hexa8(material1) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
            };
            ID2 = 1;
            for (int j = 0; j < 2; j++)
            {
                e1.NodesDictionary.Add(4 * (ID2-1) + 1, model.NodesDictionary[4 * (ID2-1) + 1]); // na allaxthei h arithmisi swsth seira
                e1.NodesDictionary.Add(4 * (ID2-1) + 2, model.NodesDictionary[4 * (ID2-1) + 2]);
                e1.NodesDictionary.Add(4 * (ID2-1) + 4, model.NodesDictionary[4 * (ID2-1) + 4]);
                e1.NodesDictionary.Add(4 * (ID2-1) + 3, model.NodesDictionary[4 * (ID2-1) + 3]);
                ID2++;
            }


            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            // ews edw

            //eisagwgh defterou element
            Element e2 = new Element()
            {
                ID = 2,
                ElementType = new Hexa8(material1) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
            };
            ID2 = 1;
            for (int j = 0; j < 2; j++)
            {
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 5, model.NodesDictionary[4 * (ID2 - 1) + 5]); // na allaxthei h arithmisi swsth seira
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 6, model.NodesDictionary[4 * (ID2 - 1) + 6]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 8, model.NodesDictionary[4 * (ID2 - 1) + 8]);
                e2.NodesDictionary.Add(4 * (ID2 - 1) + 7, model.NodesDictionary[4 * (ID2 - 1) + 7]);
                ID2++;
            }
            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
            // ews edw




            //model.Loads.Add()
            ID2 = 9;
            //apait 1
            DOFType doftype1 ;
            doftype1 = new DOFType();


            // apait 2
            Load load1;

            for (int j = 0; j < 4; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[ID2],
                    //DOF = doftype1,
                    DOF = DOFType.Z,
                    Amount = 500
                    
                };
                model.Loads.Add(load1);
                ID2++;
            }
        }

        public static void BeamElementOnly(Model model)
        {
            //gewmetria
            model.NodesDictionary.Add(13, new Node() { ID = 13, X = 0.5, Y = 0.5, Z = 0.5 });
            model.NodesDictionary.Add(14, new Node() { ID = 14, X = 0.5, Y = 0.5, Z = 1.5 });

            // constraints
            model.NodesDictionary[13].Constraints.Add(DOFType.X);
            model.NodesDictionary[13].Constraints.Add(DOFType.Y);
            model.NodesDictionary[13].Constraints.Add(DOFType.Z);
            model.NodesDictionary[13].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[13].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[13].Constraints.Add(DOFType.RotZ);
            model.NodesDictionary[14].Constraints.Add(DOFType.RotX);
            model.NodesDictionary[14].Constraints.Add(DOFType.RotY);
            model.NodesDictionary[14].Constraints.Add(DOFType.RotZ);

            //uliko
            ElasticMaterial3D material2 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };

            //diatomh
            double b = 0.3;
            double h = 0.1;

            // orismos enos element
            int subdomainID = 1;
            Element e;
            e = new Element()
            {
                ID = 3,
                ElementType = new Beam3D(material2, null, null)// kaleitai kai me tria orismata kai me ligotera
                {
                    Density = 7.85,
                    SectionArea = b * h,
                    MomentOfInertiaY = b * b * b * h,
                    MomentOfInertiaZ = b * h * h * h,
                }
            };
            e.NodesDictionary.Add(13, model.NodesDictionary[13]);
            e.NodesDictionary.Add(14, model.NodesDictionary[14]);

            model.ElementsDictionary.Add(e.ID, e);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);

            //Loads
            Load load1;
            load1 = new Load()
            {
                Node = model.NodesDictionary[14],
                //DOF = doftype1,
                DOF = DOFType.Z,
                Amount = 50

            };
            model.Loads.Add(load1);
            
        }

        public static void ExampleWithEmbedded(Model model)
        {
            HexaElementsOnly(model);
            BeamElementOnly(model);
            model.Loads.RemoveAt(4);
            model.NodesDictionary[13].Constraints.Remove(DOFType.X);
            model.NodesDictionary[13].Constraints.Remove(DOFType.Y);
            model.NodesDictionary[13].Constraints.Remove(DOFType.Z);
            var embeddedGrouping = new EmbeddedGrouping(model, model.ElementsDictionary.Where(x => x.Key < 3).Select(kv => kv.Value), model.ElementsDictionary.Where(x => x.Key >= 3).Select(kv => kv.Value));
        }
    }
}
