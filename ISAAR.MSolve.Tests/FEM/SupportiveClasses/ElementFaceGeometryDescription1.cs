using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.Tests.FEM.SupportiveClasses
{
    public class ElementFaceGeometryDescription
    {
        public static List<Node>[] GetElementFaceNodesTet10(List<Node> msolveOrderedElementNodes)
        {
            var nodes = msolveOrderedElementNodes;
            List<Node>[] facesNodes = new List<Node>[4];

            facesNodes[0] = new List<Node>() { nodes[1], nodes[2], nodes[3], nodes[5], nodes[8], nodes[9] };
            facesNodes[1] = new List<Node>() { nodes[1], nodes[0], nodes[2], nodes[4], nodes[6], nodes[5] };
            facesNodes[2] = new List<Node>() { nodes[0], nodes[3], nodes[2], nodes[7], nodes[8], nodes[6] };
            facesNodes[3] = new List<Node>() { nodes[0], nodes[1], nodes[3], nodes[4], nodes[9], nodes[7] };

            return facesNodes;
        }

        public static List<List<Node>> GetElementFaceNodesTet10( Model model, double zCoordinate)
        {
            List<Node> interfaceNodes = model.NodesDictionary.Values.Where(x => (Math.Abs(x.Z - zCoordinate) < 0.00001)).ToList();

            List<Element> interfaceElements = new List<Element>();

            foreach(var interfaceNode in interfaceNodes)
            {
                interfaceElements = interfaceElements.Union(interfaceNode.ElementsDictionary.Values.ToList()).ToList();
            }

            bool[] isUpper = new bool[interfaceElements.Count()];

            for (int i1 = 0; i1 < interfaceElements.Count(); i1++)
            {
               foreach(Node node in interfaceElements[i1].Nodes)
                {
                    if(node.Z-zCoordinate>0.00001)
                    {
                        isUpper[i1] = true;
                    }
                }
            }

            int upperCounter = 0;

            for (int i1 = 0; i1 < isUpper.Length; i1++)
            {
                if (isUpper[i1])
                {
                    upperCounter++;
                }
            }



            if(!(upperCounter==(isUpper.Length/2)))
            {
                //throw new Exception("there is a mistake in element search");
                var breakpoint = "here";
            }

            List<Element> downerInterfaceElements = new List<Element>();

            for (int i1 = 0; i1 < isUpper.Length; i1++)
            {
                if(!isUpper[i1])
                {
                    downerInterfaceElements.Add(interfaceElements[i1]);
                }
            }

            List<List<Node>> interfaceNodeSetTri6 = new List<List<Node>>();

            foreach (var elemnent in downerInterfaceElements)
            {
                var ElemnentNodesSets = GetElementFaceNodesTet10(elemnent.NodesDictionary.Values.ToList());
                foreach (var elementNOdeSet in ElemnentNodesSets)
                {
                    bool isSetInterfaceSet = true;
                    foreach (var setNode in elementNOdeSet)
                    {
                        if (Math.Abs(setNode.Z - zCoordinate) > 0.00001)
                        {
                            isSetInterfaceSet = false;
                        }
                    }
                    if(isSetInterfaceSet)
                    {
                        interfaceNodeSetTri6.Add(elementNOdeSet);
                    }
                }
                
            }
            return interfaceNodeSetTri6;
        }
    }
}
