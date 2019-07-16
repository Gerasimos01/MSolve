﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;

namespace ISAAR.MSolve.MSAnalysis.SupportiveClasses
{
    public static class DdmCalculationsAlterna2
    {
        internal static (Dictionary<int, List<int>> AssignedSubdomainsFirstLevelOfCohesive, Dictionary<int, List<int>> reassignedHexas,
            Dictionary<int, int> hexaOriginalSubdomains, Dictionary<int, List<int>> SubdomainNeedsHexas) FindEmbeddedElementsSubdomainsCorrectedSimpleFirstLevel2(Model model, int totalSubdomains,
            int[] lowerCohesiveBound, int[] upperCohesiveBound, int[] grShElementssnumber)
        {
            // origin: DdmCalculationsAlterna2
            // changes: implement subdNeedsHexasShared

            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            Dictionary<int, List<int>> SubdomainNeedsHexas = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            Dictionary<int, List<int>> reassignedHexas = new Dictionary<int, List<int>>();
            Dictionary<int, int> hexaOriginalSubdomains = new Dictionary<int, int>();

            Dictionary<Element, Dictionary<int, IList<int>>> ElementStored_HostSubdomains = new Dictionary<Element, Dictionary<int, IList<int>>>();
            Dictionary<Element, Dictionary<Subdomain, List<EmbeddedNode>>> ElementStored_elementHostSubdomainsAndNodesInThem = new Dictionary<Element, Dictionary<Subdomain, List<EmbeddedNode>>>();
            Dictionary<Element, int> ElementStored_chosenSubdomainID = new Dictionary<Element, int>();

            //1
            //Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            //Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                bool isAmbiguous = false;
                if (element.ElementType is IEmbeddedElement)
                {
                    int cohID = element.ID;
                    bool isFirstLevelCoheive = false;
                    for (int i3 = 0; i3 < lowerCohesiveBound.Length; i3++)
                    {
                        if ((lowerCohesiveBound[i3] <= cohID) & (upperCohesiveBound[i3] >= cohID))
                        {
                            isFirstLevelCoheive = true;
                        }
                    }
                    if (isFirstLevelCoheive)
                    {
                        Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                        var e1 = element.ElementType as IEmbeddedElement;
                        Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                        Dictionary<Subdomain, List<EmbeddedNode>> elementHostSubdomainsAndNodesInThem = new Dictionary<Subdomain, List<EmbeddedNode>>();//alte
                        foreach (var embeddedNode in (e1).EmbeddedNodes)
                        {
                            //1
                            Element hostELement = embeddedNode.EmbeddedInElement;
                            if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                            {
                                if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                                {
                                    HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                                    //alte
                                    elementHostSubdomainsAndNodesInThem[hostELement.Subdomain].Add(embeddedNode);
                                }
                            }
                            else
                            {
                                List<int> specificElementsIDs = new List<int>();
                                specificElementsIDs.Add(hostELement.ID);
                                HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);

                                //alte
                                List<EmbeddedNode> specificNodesIds = new List<EmbeddedNode>();
                                specificNodesIds.Add(embeddedNode);
                                elementHostSubdomainsAndNodesInThem.Add(hostELement.Subdomain, specificNodesIds);
                            }
                            //2
                            //if (hexaConnectsShellsLocal.ContainsKey(hostELement.ID))
                            //{
                            //    if (!hexaConnectsShellsLocal[hostELement.ID].Contains(element.ID))
                            //    {
                            //        hexaConnectsShellsLocal[hostELement.ID].Add(element.ID);
                            //    }
                            //}
                            //else
                            //{
                            //    List<int> connectionElementsData1 = new List<int>();
                            //    connectionElementsData1.Add(element.ID);
                            //    hexaConnectsShellsLocal.Add(hostELement.ID, connectionElementsData1);
                            //}
                        }
                        if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                        {
                            int chosenSubdomainId = 0;
                            int hexaListlength = 0;
                            foreach (int subdId in HostSubdomains.Keys)
                            {
                                if (HostSubdomains[subdId].Count > hexaListlength)
                                {
                                    chosenSubdomainId = subdId;
                                    hexaListlength = HostSubdomains[subdId].Count;
                                }
                            }
                            if (AssignedSubdomains.ContainsKey(chosenSubdomainId))
                            {
                                AssignedSubdomains[chosenSubdomainId].Add(element.ID);
                            }
                            else
                            {
                                List<int> subdElementsIds = new List<int>();
                                subdElementsIds.Add(element.ID);
                                AssignedSubdomains.Add(chosenSubdomainId, subdElementsIds);
                            }

                            foreach (int subdomainID in HostSubdomains.Keys)
                            {
                                if (!(subdomainID == chosenSubdomainId))
                                {
                                    foreach (var remoteHexa in HostSubdomains[subdomainID])
                                    {
                                        if (SubdomainNeedsHexas.ContainsKey(chosenSubdomainId))
                                        { SubdomainNeedsHexas[chosenSubdomainId].Add(remoteHexa); }
                                        else
                                        {
                                            List<int> neededHexasIds = new List<int>();
                                            neededHexasIds.Add(remoteHexa);
                                            SubdomainNeedsHexas.Add(chosenSubdomainId, neededHexasIds);
                                        }

                                        //optionally: mporei kai na mhn xreiazetai 
                                        if (hexaOriginalSubdomains.ContainsKey(remoteHexa))
                                        { }
                                        else
                                        {
                                            int originalSubdID = model.ElementsDictionary[remoteHexa].Subdomain.ID;
                                            hexaOriginalSubdomains.Add(remoteHexa, originalSubdID);
                                        }
                                    }
                                }
                            }

                            ElementStored_HostSubdomains.Add(element, HostSubdomains);
                            ElementStored_elementHostSubdomainsAndNodesInThem.Add(element, elementHostSubdomainsAndNodesInThem);
                            ElementStored_chosenSubdomainID.Add(element, chosenSubdomainId);

                        }
                        if (HostSubdomains.Count == 1)
                        {
                            if (AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key))
                            {
                                AssignedSubdomains[HostSubdomains.ElementAt(0).Key].Add(element.ID);
                            }
                            else
                            {
                                List<int> subdElementsIds = new List<int>();
                                subdElementsIds.Add(element.ID);
                                AssignedSubdomains.Add(HostSubdomains.ElementAt(0).Key, subdElementsIds);
                            }
                        }
                    }

                }

            }

            foreach (Element element in ElementStored_HostSubdomains.Keys)
            {
                Dictionary<int, IList<int>> HostSubdomains = ElementStored_HostSubdomains[element];
                Dictionary<Subdomain, List<EmbeddedNode>> elementHostSubdomainsAndNodesInThem = ElementStored_elementHostSubdomainsAndNodesInThem[element];
                int chosenSubdomainId = ElementStored_chosenSubdomainID[element];

                foreach (int subdomainID in HostSubdomains.Keys)
                {
                    if (!(subdomainID == chosenSubdomainId))
                    {
                        foreach (var remoteHexa in HostSubdomains[subdomainID])
                        {
                            int[] adjacentHexaIds = GetAdjacentHexaElementsIds(model, new int[1] { remoteHexa });
                            bool isRemoteHexaConnected = false;
                            for (int i1 = 0; i1 < adjacentHexaIds.Length; i1++)
                            {
                                if (model.ElementsDictionary[adjacentHexaIds[i1]].Subdomain.ID == chosenSubdomainId)
                                {
                                    isRemoteHexaConnected = true;
                                }
                            }

                            //int ConnectorHexaID;
                            //bool ConnectorHexaFound = false;                            
                            if (!isRemoteHexaConnected)
                            {
                                int[] chosenSubdomainAdjacentHexas = GetAdjacentHexaElementsIds(model, HostSubdomains[chosenSubdomainId].ToArray());
                                var possibleSolutions = adjacentHexaIds.Intersect(chosenSubdomainAdjacentHexas);
                                if (possibleSolutions.Count() > 0)
                                {
                                    bool ConnectionAlreadyExists = false;
                                    foreach (int possibleSolutionID in possibleSolutions)
                                    {
                                        if (!(SubdomainNeedsHexas[chosenSubdomainId] == null))
                                        {
                                            if (SubdomainNeedsHexas[chosenSubdomainId].Contains(possibleSolutionID))
                                            {
                                                ConnectionAlreadyExists = true;
                                            }
                                        }
                                    }
                                    if (!ConnectionAlreadyExists)
                                    {
                                        if (!(SubdomainNeedsHexas[chosenSubdomainId] == null))
                                        {
                                            SubdomainNeedsHexas[chosenSubdomainId].Add(possibleSolutions.ElementAt(0));
                                        }
                                        else
                                        {
                                            SubdomainNeedsHexas[chosenSubdomainId] = new List<int>() { possibleSolutions.ElementAt(0) };
                                        }
                                    }

                                }
                                else
                                {
                                    bool wereAdjacentElementsFound = false;
                                    foreach(int adjacentHexaId in adjacentHexaIds)
                                    {
                                        int[] adjacentHexaIds2 = GetAdjacentHexaElementsIds(model, new int[1] { adjacentHexaId });
                                        possibleSolutions = adjacentHexaIds2.Intersect(chosenSubdomainAdjacentHexas);
                                        if(possibleSolutions.Count()>0)
                                        {
                                            
                                            wereAdjacentElementsFound = true;
                                            //ADD adjacent to chosenSubdomain elements
                                            bool ConnectionAlreadyExists = false;
                                            foreach (int possibleSolutionID in possibleSolutions)
                                            {
                                                if (!(SubdomainNeedsHexas[chosenSubdomainId] == null))
                                                {
                                                    if (SubdomainNeedsHexas[chosenSubdomainId].Contains(possibleSolutionID))
                                                    {
                                                        ConnectionAlreadyExists = true;
                                                    }
                                                }
                                            }
                                            if (!ConnectionAlreadyExists)
                                            {
                                                if (!(SubdomainNeedsHexas[chosenSubdomainId] == null))
                                                {
                                                    SubdomainNeedsHexas[chosenSubdomainId].Add(possibleSolutions.ElementAt(0));
                                                }
                                                else
                                                {
                                                    SubdomainNeedsHexas[chosenSubdomainId] = new List<int>() { possibleSolutions.ElementAt(0) };
                                                }
                                            }
                                            //ADD adjacent to remote needed hexa element as well
                                            if (!(SubdomainNeedsHexas[chosenSubdomainId] == null))
                                            {
                                                SubdomainNeedsHexas[chosenSubdomainId].Add(adjacentHexaId);
                                            }
                                            else
                                            {
                                                SubdomainNeedsHexas[chosenSubdomainId] = new List<int>() { adjacentHexaId };
                                            }
                                            break;

                                        }
                                        
                                    }
                                    if (!wereAdjacentElementsFound)
                                    {
                                        throw new NotImplementedException();
                                    }
                                    
                                }
                            }


                        }
                    }
                }

            }

            return (AssignedSubdomains, reassignedHexas, hexaOriginalSubdomains, SubdomainNeedsHexas);
        }

        public static int[] GetAdjacentHexaElementsIds(Model model, int[] hexaIDs)
        {
            int[][] hexaNodePositions = new int[6][]
                {new int[2] {0,2},
                 new int[2] {4,6},
                 new int[2] {4,3},
                 new int[2] {2,5},
                 new int[2] {0,5},
                 new int[2] {3,6} };
            List<int> adjacentHexas = new List<int>();
            foreach (int HexaId in hexaIDs)
            {
                Element hexaElement = model.ElementsDictionary[HexaId];
                for (int nodePairPosition = 0; nodePairPosition < hexaNodePositions.Length; nodePairPosition++)
                {
                    int[] nodePair = hexaNodePositions[nodePairPosition];
                    Node node1 = model.ElementsDictionary[HexaId].NodesDictionary.ElementAt(nodePair[0]).Value;
                    Node node2 = model.ElementsDictionary[HexaId].NodesDictionary.ElementAt(nodePair[1]).Value;
                    List<Element> node1ELements = new List<Element>();
                    List<Element> node2ELements = new List<Element>();

                    foreach (Element element in node1.ElementsDictionary.Values)
                    {
                        if (!(element.ElementType is IEmbeddedElement))
                        {
                            node1ELements.Add(element);
                        }
                    }
                    foreach (Element element in node2.ElementsDictionary.Values)
                    {
                        if (!(element.ElementType is IEmbeddedElement))
                        {
                            node2ELements.Add(element);
                        }
                    }

                    foreach (Element element in node1ELements)
                    {
                        if (node2ELements.Contains(element))
                        {
                            if (!(hexaIDs.Contains(element.ID)))
                            {
                                if (!(adjacentHexas.Contains(element.ID)))
                                {
                                    adjacentHexas.Add(element.ID);
                                }
                            }
                        }
                    }

                }
            }

            return adjacentHexas.ToArray();

        }

        public static (Dictionary<int, List<int>>, Dictionary<int, int>) GetHexaSharingSubdomains(Dictionary<int, List<int>> SubdomainNeedsHexas, int[][] hexaIds,
            bool withOriginal)
        {
            Dictionary<int, List<int>> hexaAndTheirSharingSubdomains = new Dictionary<int, List<int>>();
            Dictionary<int, int> DuplicateHexaOriginalSubdomain= new Dictionary<int, int>();

            foreach (int subdomainID in SubdomainNeedsHexas.Keys)
            {
                foreach (int hexaID in SubdomainNeedsHexas[subdomainID])
                {
                    if (!(hexaAndTheirSharingSubdomains.ContainsKey(hexaID)))
                    {
                        List<int> hexaSharingSubdomainsofElement = new List<int>() { subdomainID };
                        hexaAndTheirSharingSubdomains[hexaID] = hexaSharingSubdomainsofElement;
                    }
                    else
                    {
                        if (!(hexaAndTheirSharingSubdomains[hexaID].Contains(subdomainID)))
                        {
                            hexaAndTheirSharingSubdomains[hexaID].Add(subdomainID);
                        }
                    }
                }
                
            }


            for(int i1=0; i1<hexaIds.Length; i1++)
            {
                if(!(hexaIds[i1]==null))
                {
                    for (int i2 = 0; i2 < hexaIds[i1].Length; i2++)
                    {
                        if (withOriginal)
                        {
                            hexaAndTheirSharingSubdomains[hexaIds[i1][i2]].Add(i1);
                        }
                        DuplicateHexaOriginalSubdomain[hexaIds[i1][i2]] = i1;
                    }
                }
            }

            return (hexaAndTheirSharingSubdomains,DuplicateHexaOriginalSubdomain);
        }



        public static int[][] ModelAddDuplicateHexas(Model model, rveMatrixParameters  mp, Dictionary<int, List<int>> hexaAndTheirSharingSubdomains,int totalSubdomains)
        {

            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            int ElementID = model.Elements.Count() + 1;
            Dictionary<int, List<int>> subdAdditionalHexas = new Dictionary<int, List<int>>(totalSubdomains);


            foreach (int hexaID in hexaAndTheirSharingSubdomains.Keys)
            {
                //ORIGINAL HEXA
                double totalDuplicates = (double)(hexaAndTheirSharingSubdomains[hexaID].Count()+1);
                ElasticMaterial3D material1 = new ElasticMaterial3D()
                {
                    YoungModulus = E_disp/totalDuplicates,
                    PoissonRatio = ni_disp,
                };

                Element e1 = new Element()
                {
                    ID = hexaID,
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };

                foreach(Node node in model.ElementsDictionary[hexaID].NodesDictionary.Values)
                {
                    e1.NodesDictionary.Add(node.ID,node);
                }
                int originalSubdomainID = model.ElementsDictionary[hexaID].Subdomain.ID;
                model.ElementsDictionary[hexaID] = e1;
                model.SubdomainsDictionary[originalSubdomainID].Elements.Add(e1);

                //ADDITIONAL HEXAS
                foreach (int subdomainID in hexaAndTheirSharingSubdomains[hexaID])
                {
                    e1 = new Element()
                    {
                        ID = ElementID,
                        ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                    };

                    foreach (Node node in model.ElementsDictionary[hexaID].NodesDictionary.Values)
                    {
                        e1.NodesDictionary.Add(node.ID, node);
                    }

                    model.ElementsDictionary[ElementID] = e1;
                    model.SubdomainsDictionary[subdomainID].Elements.Add(e1);

                    if (!subdAdditionalHexas.ContainsKey(subdomainID))
                    {
                        List<int> additionalHexas1 = new List<int>() { ElementID };
                        subdAdditionalHexas[subdomainID] = additionalHexas1;
                    }
                    else
                    {
                        subdAdditionalHexas[subdomainID].Add(ElementID);
                    }

                    ElementID++;
                }
            }

            int[][] subdAdditionalHexasArray = DdmCalculationsPartb.ConvertIntListToArray(subdAdditionalHexas, totalSubdomains);
            return subdAdditionalHexasArray;
        }
    }
}