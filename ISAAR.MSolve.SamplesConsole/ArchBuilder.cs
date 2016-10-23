using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class ArchBuilder
    {
        public static void MakeArchBuilding(Model model)
        {
            int startX = 0;
            int startY = 0;
            int startZ = 0;

            int nodeID = 1;
            for (int l = 0; l < 3; l++)
            {
               model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ});
               nodeID++;                                   
            }
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX+0.5, Y = startY + l * 1, Z = startZ+0.5 });
                nodeID++;
            }
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX+1.0, Y = startY + l * 0.5, Z = startZ+1.0 });
                nodeID++;
            }

            nodeID = 1;
            for (int j = 0; j < 3; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                nodeID++;
                
            }

            nodeID = 6;
            for (int j = 0; j < 3; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.RotY);
                nodeID++;

            }
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e8,
                PoissonRatio = 0.3,
            };

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = -0.7071067812;
                VH[j][1] = 0;
                VH[j][2] = 0.7071067812;
            }

            double[] Tk = new double[8];
            for (int j = 0; j < 8; j++)
            {
                Tk[j] = 0.1;
            }
            //eisagwgh enos element
            Element e1;
            int ID2 = 1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Shell8disp(material1, 3, 3, 2)
                {
                    oVn_i = VH,
                    tk=Tk,
                }// dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                                
            };

            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);
            
            //e1.NodesDictionary.Add(4 * (ID2-1) + 1, model.NodesDictionary[4 * (ID2-1) + 1]); // na allaxthei h arithmisi swsth seira
                


            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            // ews edw
            


            //model.Loads.Add()
            ID2 = 4;           
            Load load1;            
            load1 = new Load()
            {
                    Node = model.NodesDictionary[ID2],
                    //DOF = doftype1,
                    DOF = DOFType.Z,
                    Amount = -50
                    
             };
             model.Loads.Add(load1);
             ID2++;
        }
        
    }
}
