using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class HexaBuilder
    {
        public static void MakeHexaBuilding(Model model)
        {
            int startX = 0;
            int startY = 0;
            int startZ = 0;

            int nodeID = 1;
            for (int l = 0; l < 3; l++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
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
            ElasticMaterial3D material = new ElasticMaterial3D()
            {
                YoungModulus = 2.1e5,
                PoissonRatio = 0.35,
            };
            Element e;
            int ID = 1;
            e = new Element();
            {
                ID = 1;
                e.ElementType = new Hexa8(material); // dixws to e. exoume sfalma enw sto beambuilding oxi
            }
            ID = 1;
            for (int j = 0; j < 8; j++)
            {
                e.NodesDictionary.Add(ID, model.NodesDictionary[ID]);
                ID++;
            }
            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e.ID, e);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);

        }
    }
}
