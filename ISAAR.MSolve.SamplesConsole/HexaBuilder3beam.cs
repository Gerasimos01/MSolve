using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class HexaBuilder3beam
    {
        public static void MakeHexaBuilding(Model model)
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
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 1, model.NodesDictionary[4 * (ID2 - 1) + 1]); // na allaxthei h arithmisi swsth seira
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 2, model.NodesDictionary[4 * (ID2 - 1) + 2]);
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 4, model.NodesDictionary[4 * (ID2 - 1) + 4]);
                e1.NodesDictionary.Add(4 * (ID2 - 1) + 3, model.NodesDictionary[4 * (ID2 - 1) + 3]);
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
            DOFType doftype1;
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
                    Amount = 5

                };
                model.Loads.Add(load1);
                ID2++;
            }

            //prosthiki embedded
            //prosthiki shmeiwn 
            double startX1 = 0.25;
            double startY1 = 0.25;
            double startZ1 = 0.5;
            nodeID = 13;
            for (int l = 0; l < 1; l++)
            {                                                                      
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX1 + l * 0.5, Y = startY1 + l * 0.5, Z = startZ1 + l * 1 });
                        nodeID++;                                  
            }
            //prosthiki ulikou gia to beam 3d 
            ElasticMaterial material2 = new ElasticMaterial()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };
            //prosthiki element
            double sctn = 0.05 ;
            e2 = new Element()
            {
                ID = 3,
                ElementType = new Beam3D(material2) {SectionArea=sctn,MomentOfInertiaY=0.001*sctn,MomentOfInertiaZ=0.001*sctn },                
            };
            
            e2.NodesDictionary.Add(13, model.NodesDictionary[13]); 
            e2.NodesDictionary.Add(14, model.NodesDictionary[14]);
  
            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);

            //embeding shmeiwn

            //EmbeddedNode node1 = e2.BuildHostElementEmbeddedNode();
            //IList<DOFType> doftypelist1 = new IList<DOFType>(); <--den ginetai

            ////IList<DOFType> doftypelist1 = new List<DOFType>();
            //////List<DOFType> doftypelist1 = new List<DOFType>();    <--ginetai        
            ////doftypelist1.Add(DOFType.X);
            ////doftypelist1.Add(DOFType.Y);
            ////doftypelist1.Add(DOFType.Z);
            ////EmbeddedNode n4;
            ////n4 = new EmbeddedNode(model.NodesDictionary[13], model.ElementsDictionary[1], doftypelist1);

            ////ID2 = 13;
            ////for (int j = 0; j < 4; j++)
            ////{
            ////    n4 = new EmbeddedNode(model.NodesDictionary[ID2], model.ElementsDictionary[1], doftypelist1);
            ////    model.ElementsDictionary[1].EmbeddedNodes.Add(n4);
            ////    ID2++;
            ////}

            //ID2 = 13;
            //for (int j = 0; j < 4; j++)
            //{
            //    model.ElementsDictionary[1].EmbeddedNodes.Add(model.NodesDictionary[ID2]);
            //    ID2++;
            //}

            //ID2 = 17;
            //for (int j = 0; j < 4; j++)
            //{
            //    model.ElementsDictionary[2].EmbeddedNodes.Add(model.NodesDictionary[ID2]);
            //    ID2++;
            //}

            
            

        }
    }
}
