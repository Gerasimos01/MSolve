using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.MSAnalysis.SupportiveClasses
{
    public static class DdmCalculationsAlterna
    {
        public static Dictionary<int, List<int>> FindEmbeddedElementsSubdomainsCorrectedSimple(Model model, int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            Dictionary<int, Subdomain> elementOriginalSubdomains = new Dictionary<int, Subdomain>(); //element kai h subdomain sthn opoia anhkei
            Dictionary<int, Dictionary<Subdomain, List<int>>> revisitedElementsIds = new Dictionary<int, Dictionary<Subdomain, List<int>>>(); // <elementId <subdomain kai poia embedded nodes se afthn>>


            //1
            //Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            //Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                    var e1 = element.ElementType as IEmbeddedElement;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    Dictionary<Subdomain, List<int>> elementHostSubdomainsAndNodesInThem = new Dictionary<Subdomain, List<int>>();//alte
                    int embeddedNodesNesessary = 3;
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
                                elementHostSubdomainsAndNodesInThem[hostELement.Subdomain].Add(embeddedNode.Node.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);

                            //alte
                            List<int> specificNodesIds = new List<int>();
                            specificNodesIds.Add(embeddedNode.Node.ID);
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

                        //alte
                        if (elementHostSubdomainsAndNodesInThem[model.SubdomainsDictionary[chosenSubdomainId]].Count() < embeddedNodesNesessary)
                        {
                            revisitedElementsIds.Add(element.ID, elementHostSubdomainsAndNodesInThem);
                            elementOriginalSubdomains.Add(element.ID, model.SubdomainsDictionary[chosenSubdomainId]);
                        }

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

            foreach (int elementId in revisitedElementsIds.Keys)
            {
                Subdomain assignedSubdomain;
                var element = model.ElementsDictionary[elementId];
                List<Element> elementAdkjacents = new List<Element>();
                foreach (Node node in element.Nodes.Skip(8))
                {
                    if (node.ElementsDictionary.Count() == 2)
                    {
                        foreach (Element possibleAdjacent in node.ElementsDictionary.Values)
                        {
                            if (!(possibleAdjacent.ID == element.ID))
                            {
                                elementAdkjacents.Add(possibleAdjacent);
                            }
                        }
                    }
                }
                if (elementAdkjacents.Count > 4)
                {
                    throw new NotImplementedException();
                }
                else
                {
                    //gather adjacent in subdomains 
                    Dictionary<int, List<Element>> subdomainsAndTheirAdjacents = new Dictionary<int, List<Element>>();
                    foreach (Element adjacent in elementAdkjacents)
                    {
                        Subdomain hostSubdomain = adjacent.Subdomain;
                        if (subdomainsAndTheirAdjacents.Keys.Contains(hostSubdomain.ID))
                        {
                            if (!subdomainsAndTheirAdjacents[hostSubdomain.ID].Contains(adjacent))
                            {
                                subdomainsAndTheirAdjacents[hostSubdomain.ID].Add(adjacent);
                            }
                        }
                        else
                        {
                            List<Element> adjacentElements = new List<Element>() { adjacent };
                            subdomainsAndTheirAdjacents.Add(hostSubdomain.ID, adjacentElements);
                        }
                    }

                    //choose subdomain 
                    bool foundSubdomain = false;
                    foreach (int hostSubdomainID in subdomainsAndTheirAdjacents.Keys)
                    {
                        if (subdomainsAndTheirAdjacents[hostSubdomainID].Count() > 1)
                        {
                            foundSubdomain = true;
                            assignedSubdomain = model.SubdomainsDictionary[hostSubdomainID];
                        }
                    }

                    if (!foundSubdomain)
                    {
                        Subdomain chosenSubdomainId = null;
                        int hexaListlength = 0;
                        foreach(Subdomain subdId in revisitedElementsIds[elementId].Keys)
                        {
                            if (revisitedElementsIds[elementId][subdId].Count > hexaListlength)
                            {
                                chosenSubdomainId = subdId;
                                hexaListlength = revisitedElementsIds[elementId][subdId].Count;
                            }
                        }
                        assignedSubdomain = chosenSubdomainId;

                        var originalSubdomainID = elementOriginalSubdomains[element.ID].ID;
                        AssignedSubdomains[originalSubdomainID].Remove(element.ID);
                        AssignedSubdomains[assignedSubdomain.ID].Add(element.ID);
                    }

                }
                
            }

            return AssignedSubdomains;
        }

        internal static Dictionary<int, List<int>> FindEmbeddedElementsSubdomainsCorrectedSimpleFirstLevel(Model model, int totalSubdomains,
            int[] lowerCohesiveBound, int[] upperCohesiveBound, int[] grShElementssnumber)
        {
            // origin: DdmCalculations.FindEmbeddedElementsSubdomainsCorrectedSimple()
            // changes: efarmogh mono sta cohesive tou prwtou level 

            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            //1
            //Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            //Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
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
                    { Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                        var e1 = element.ElementType as IEmbeddedElement;
                        Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                        foreach (var embeddedNode in (e1).EmbeddedNodes)
                        {
                            //1
                            Element hostELement = embeddedNode.EmbeddedInElement;
                            if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                            {
                                if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                                {
                                    HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                                }
                            }
                            else
                            {
                                List<int> specificElementsIDs = new List<int>();
                                specificElementsIDs.Add(hostELement.ID);
                                HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);
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
            return AssignedSubdomains;
        }

        internal static Dictionary<int, List<int>> FindEmbeddedElementsSubdomainsCorrectedSimpleSecondLevel(Model model, int totalSubdomains,
            int[] lowerCohesiveBound, int[] upperCohesiveBound, int[] grShElementssnumber, Dictionary<int, List<int>> assignedSubdomainsFirstLevelOfCohesive)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(assignedSubdomainsFirstLevelOfCohesive.Keys.Count());
            foreach (int subdomainID in assignedSubdomainsFirstLevelOfCohesive.Keys)
            {
                AssignedSubdomains.Add(subdomainID, new List<int>(2 * assignedSubdomainsFirstLevelOfCohesive[subdomainID].Count()));
                foreach (int firsLevelElemId in assignedSubdomainsFirstLevelOfCohesive[subdomainID])
                {
                    AssignedSubdomains[subdomainID].Add(firsLevelElemId);
                }
            }

            foreach (int subdomainID in assignedSubdomainsFirstLevelOfCohesive.Keys)
            {
                foreach(int firsLevelElemId in assignedSubdomainsFirstLevelOfCohesive[subdomainID])
                {
                    int graoheneSheetId = 0;
                    for (int i3 = 0; i3 < lowerCohesiveBound.Length; i3++)
                    {
                        if ((lowerCohesiveBound[i3] <= firsLevelElemId) & (upperCohesiveBound[i3] >= firsLevelElemId))
                        {
                            graoheneSheetId = i3;
                        }
                    }

                    AssignedSubdomains[subdomainID].Add(firsLevelElemId + grShElementssnumber[graoheneSheetId]);
                }

            }

            return AssignedSubdomains;
        }

        internal static (Dictionary<int, List<int>> AssignedSubdomainsFirstLevelOfCohesive, Dictionary<int, List<int>> reassignedHexas,
            Dictionary<int, int> hexaOriginalSubdomains) FindEmbeddedElementsSubdomainsCorrectedSimpleFirstLevel2(Model model, int totalSubdomains,
            int[] lowerCohesiveBound, int[] upperCohesiveBound, int[] grShElementssnumber)
        {
            // origin: DdmCalculations.FindEmbeddedElementsSubdomainsCorrectedSimple()
            // changes: efarmogh mono sta cohesive tou prwtou level 

            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            Dictionary<int, List<int>> reassignedHexas = new Dictionary<int, List<int>>();
            Dictionary<int, int> hexaOriginalSubdomains = new Dictionary<int, int>();

            //1
            //Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            //Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
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

                            ////alte TODO: adjust this 
                            //foreach(Subdomain subdomain in elementHostSubdomainsAndNodesInThem.Keys)
                            //{
                            //    if (elementHostSubdomainsAndNodesInThem[subdomain].Count == 1)
                            //    {

                            //        var newSubdID = chosenSubdomainId;
                            //        EmbeddedNode node = elementHostSubdomainsAndNodesInThem[subdomain].First();
                            //        int reassignedHexaID = node.EmbeddedInElement.ID;
                            //        int originalSubdIID = node.EmbeddedInElement.Subdomain.ID;

                            //        if (hexaOriginalSubdomains.ContainsKey(reassignedHexaID))
                            //        { }
                            //        else
                            //        {
                            //            if (reassignedHexas.Keys.Contains(newSubdID))
                            //            {
                            //                reassignedHexas[newSubdID].Add(reassignedHexaID);
                            //            }
                            //            else
                            //            {
                            //                var subdHexas = new List<int>() { reassignedHexaID };
                            //                reassignedHexas.Add(newSubdID, subdHexas);
                            //            }
                            //            hexaOriginalSubdomains.Add(reassignedHexaID, originalSubdIID);
                            //        }

                            //    }
                            //    if (elementHostSubdomainsAndNodesInThem[subdomain].Count == 2)
                            //    {
                            //        {
                            //            var newSubdID = chosenSubdomainId;
                            //            EmbeddedNode node = elementHostSubdomainsAndNodesInThem[subdomain].ElementAt(0);
                            //            int reassignedHexaID = node.EmbeddedInElement.ID;
                            //            int originalSubdIID = node.EmbeddedInElement.Subdomain.ID;

                            //            if (hexaOriginalSubdomains.ContainsKey(reassignedHexaID))
                            //            { }
                            //            else
                            //            {
                            //                if (reassignedHexas.Keys.Contains(newSubdID))
                            //                {
                            //                    reassignedHexas[newSubdID].Add(reassignedHexaID);
                            //                }
                            //                else
                            //                {
                            //                    var subdHexas = new List<int>() { reassignedHexaID };
                            //                    reassignedHexas.Add(newSubdID, subdHexas);
                            //                }
                            //                hexaOriginalSubdomains.Add(reassignedHexaID, originalSubdIID);
                            //            }
                            //        }
                            //        {
                            //            var newSubdID = chosenSubdomainId;
                            //            EmbeddedNode node = elementHostSubdomainsAndNodesInThem[subdomain].ElementAt(1);
                            //            int reassignedHexaID = node.EmbeddedInElement.ID;
                            //            int originalSubdIID = node.EmbeddedInElement.Subdomain.ID;

                            //            if (hexaOriginalSubdomains.ContainsKey(reassignedHexaID))
                            //            { }
                            //            else
                            //            {
                            //                if (reassignedHexas.Keys.Contains(newSubdID))
                            //                {
                            //                    reassignedHexas[newSubdID].Add(reassignedHexaID);
                            //                }
                            //                else
                            //                {
                            //                    var subdHexas = new List<int>() { reassignedHexaID };
                            //                    reassignedHexas.Add(newSubdID, subdHexas);
                            //                }
                            //                hexaOriginalSubdomains.Add(reassignedHexaID, originalSubdIID);
                            //            }
                            //        }
                            //    }
                            //}

                            //alte TODO: adjust this 
                            foreach (Subdomain subdomain in elementHostSubdomainsAndNodesInThem.Keys)
                            {
                                if(!(subdomain.ID==chosenSubdomainId))
                                {
                                    IList<int> hexaIds = HostSubdomains[subdomain.ID];

                                    foreach(int reassignedHexaId in hexaIds)
                                    {
                                        var newSubdID = chosenSubdomainId;
                                        //EmbeddedNode node = elementHostSubdomainsAndNodesInThem[subdomain].First();
                                        //int reassignedHexaId = node.EmbeddedInElement.ID;
                                        int originalSubdIID = model.ElementsDictionary[reassignedHexaId].Subdomain.ID;

                                        if (hexaOriginalSubdomains.ContainsKey(reassignedHexaId))
                                        { }
                                        else
                                        {
                                            if (reassignedHexas.Keys.Contains(newSubdID))
                                            {
                                                reassignedHexas[newSubdID].Add(reassignedHexaId);
                                            }
                                            else
                                            {
                                                var subdHexas = new List<int>() { reassignedHexaId };
                                                reassignedHexas.Add(newSubdID, subdHexas);
                                            }
                                            hexaOriginalSubdomains.Add(reassignedHexaId, originalSubdIID);
                                        }
                                    }
                                }
                            }


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
            return (AssignedSubdomains, reassignedHexas, hexaOriginalSubdomains);
        }
    }
}
