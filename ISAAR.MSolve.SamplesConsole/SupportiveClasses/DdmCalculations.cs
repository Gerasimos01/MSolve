using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    public static class  DdmCalculations
    {
        public static int[][] CalculateSubdElementIds(int hexa1, int hexa2, int hexa3, int elem1, int elem2, Model model)
        {
            int[][] subdElementIds = new int[8][];
            for (int i1 = 0; i1 < 7; i1++)
            {subdElementIds[i1] = new int[hexa1 * hexa2 * hexa3 / 8];}
            subdElementIds[7] = new int[(hexa1 * hexa2 * hexa3 / 8) + 3 * elem1 * elem2];

            int [] subdElementCounters = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        int ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based

                        int s1; int s2; int s3;
                        if (h1 <= 0.5 * hexa1-1) { s1 = 1; } else { s1 = 2; };
                        if (h2 <= 0.5 * hexa2-1) { s2 = 1; } else { s2 = 2; };
                        if (h3 <= 0.5 * hexa3-1) { s3 = 1; } else { s3 = 2; };

                        int subdID = s1 + (s2 - 1) * 2 + (s3 - 1) * 4;

                        subdElementIds[subdID - 1][subdElementCounters[subdID - 1]] = ElementID;
                        subdElementCounters[subdID - 1] += 1;
                        //model.ElementsDictionary.Add(e1.ID, e1);
                        //model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);

                    }
                }
            }

            
            //int subdID = 8;

            for (int ElementID = hexa1 * hexa2 * hexa3 + 1; ElementID < hexa1 * hexa2 * hexa3+ 3*(elem1*elem2)+1; ElementID++)
            {
                subdElementIds[7][subdElementCounters[7]] = ElementID;
                subdElementCounters[7] += 1;
            }
            return subdElementIds;
        }

        public static int[][] CalculateSubdElementIdsVerticalHexaOnly(int hexa1, int hexa2, int hexa3, int elem1, int elem2, Model model,int nSubdomains)
        {
            //int[][] subdElementIds = new int[8][];
            int[][] subdElementIds = new int[nSubdomains][];
            for (int i1 = 0; i1 < nSubdomains; i1++)
            { subdElementIds[i1] = new int[hexa1 * hexa2 * hexa3 / nSubdomains]; }
           


            int[] subdElementCounters = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        int ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based

                        int s1; int s2; int s3;
                        if (h1 <= 0.5 * hexa1 - 1) { s1 = 1; } else { s1 = 2; };
                        if (h2 <= 0.5 * hexa2 - 1) { s2 = 1; } else { s2 = 2; };
                        if (h3 <= 0.5 * hexa3 - 1) { s3 = 1; } else { s3 = 2; };

                        int subdID = s1 + (s2 - 1) * 2 + (s3 - 1) * 4;

                        subdElementIds[subdID - 1][subdElementCounters[subdID - 1]] = ElementID;
                        subdElementCounters[subdID - 1] += 1;
                        //model.ElementsDictionary.Add(e1.ID, e1);
                        //model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);

                    }
                }
            }


            //int subdID = 8;

            for (int ElementID = hexa1 * hexa2 * hexa3 + 1; ElementID < hexa1 * hexa2 * hexa3 + 3 * (elem1 * elem2) + 1; ElementID++)
            {
                subdElementIds[7][subdElementCounters[7]] = ElementID;
                subdElementCounters[7] += 1;
            }
            return subdElementIds;
        }

        public static void SeparateSubdomains(Model model, int[][] subdElementIds)
        {
            model.SubdomainsDictionary.Clear();

            for (int subdID = 0; subdID < subdElementIds.GetLength(0); subdID++)
            {
                model.SubdomainsDictionary.Add(subdID, new Subdomain() { ID = subdID });
                for (int i1 = 0; i1 < subdElementIds[subdID].GetLength(0); i1++)
                {
                    model.SubdomainsDictionary[subdID].ElementsDictionary.Add(subdElementIds[subdID][i1], model.ElementsDictionary[subdElementIds[subdID][i1]]);
                }
            }

        }

        public static void PrintDictionary (Dictionary<int, Dictionary<DOFType, int>> globalNodalDOFsDictionary,int TotalDOFs, int subdomainID)
        {
            double[] globalDOFs = new double[TotalDOFs];
            int counter = 0;

            foreach (int nodeID in globalNodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = globalNodalDOFsDictionary[nodeID];
                //Dictionary<DOFType, int> globalDOFTypes = new Dictionary<DOFType, int>(dofTypes.Count);
                foreach (DOFType dofType in dofTypes.Keys)
                {
                    if (dofTypes[dofType]!=-1)
                    {
                        globalDOFs[counter] = dofTypes[dofType];
                        counter += 1;
                    }
                }

            }

            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}globalDOFs.txt";
            string file_no = subdomainID.ToString();
            string print_path = string.Format(print_path_gen, file_no);
            Vector globalDOFS = new Vector(globalDOFs);
            globalDOFS.WriteToFile(print_path);
        }

        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomains(Model model)
        {
            //1
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            //2
            Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values)
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in ((IEmbeddedElement)element).EmbeddedNodes)
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
                        if (hexaConnectsShells.ContainsKey(hostELement.ID))
                        {
                            if (!hexaConnectsShells[hostELement.ID].Contains(element.ID))
                            {
                                hexaConnectsShells[hostELement.ID].Add(element.ID);
                            }
                        }
                        else
                        {
                            List<int> connectionElementsData1 = new List<int>();
                            connectionElementsData1.Add(element.ID);
                            hexaConnectsShells.Add(hostELement.ID, connectionElementsData1);
                        }
                    }
                    if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                    { EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem.Add(element.ID, HostSubdomains); }
                }
            }
            return (EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,hexaConnectsShells);
        }

        public static Dictionary<int, Dictionary<int, IList<int>>> FindAmbiguousEmbeddedElementsSubdomains(Model model)
        {
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndElements = new Dictionary<int, Dictionary<int, IList<int>>>();

            foreach (Element element in model.ElementsDictionary.Values)
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in ((IEmbeddedElement)element).EmbeddedNodes)
                    {
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
                    }

                    if (HostSubdomains.Count > 1)
                    { EmbeddedElementsHostSubdomainsAndElements.Add(element.ID, HostSubdomains); }
                }
            }
            return EmbeddedElementsHostSubdomainsAndElements;
        }

        public static void DetermineOnlyNeededCombinations( 
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
            Dictionary<int, List<int>> hexaConnectsShells)
        {
            //List<List<>> CombinationElements
            Dictionary<int, List<int>> connectedShellElementsLists = new Dictionary<int, List<int>>();

            //List<int> connectedShells1 = new List<int>();

            foreach (int hexaID in hexaConnectsShells.Keys)
            {
                //foreach (int connectedShell in hexaConnectsShells[hexaID])
                //{
                List<int> foundInLists = new List<int>();

                foreach (int ListID in connectedShellElementsLists.Keys)
                {
                    bool isShellFoundInList = false;
                    
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        if (connectedShellElementsLists[ListID].Contains(shellELement))
                        {
                            isShellFoundInList = true;

                            if (!foundInLists.Contains(ListID))
                            { foundInLists.Add(ListID); }
                            //prosthese kai ta upoloipa shell sth lista
                            //kai oles tis upoloipes listes pou ta periexoun 
                        }
                    }
                }

                int foundInListsNum = foundInLists.Count();
                if (foundInListsNum==0)
                {
                    List<int> newConnectedShellsList = new List<int>();
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        newConnectedShellsList.Add(shellELement);
                    }

                    connectedShellElementsLists.Add(connectedShellElementsLists.Count()+1, newConnectedShellsList);
                }
                if (foundInListsNum==1)
                {
                    var updatedList = connectedShellElementsLists[foundInLists.ElementAt(0)];
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        if(!updatedList.Contains(shellELement))
                        { updatedList.Add(shellELement); }

                    }
                }
                if (foundInListsNum>1)
                {
                    //lists[1] = lists[1].Union(lists[2]).ToList();
                    //lists.Remove(2);

                    for (int i1=1; i1<foundInListsNum; i1++)
                    {
                        connectedShellElementsLists[foundInLists.ElementAt(0)] = 
                            connectedShellElementsLists[foundInLists.ElementAt(0)].Union(connectedShellElementsLists[foundInLists.ElementAt(i1)]).ToList();
                        //todo: concat can be used as well if it is known that there are not duplicates
                        connectedShellElementsLists.Remove(foundInLists.ElementAt(i1));
                    }
                }



                //connectedShellElementsLists[0].Union(connectedShellElementsLists[1]);

                //}
            }
            


        }

        public static void CalculateCombinationSolution(List<int> connectedShellElementsLists, Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndElements)
        {
            int solutionVectorSize = connectedShellElementsLists.Count();
            int possibleSolutions = 1;            
            foreach(int shellId in connectedShellElementsLists)
            {
                possibleSolutions *= EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
            }
            List<int>[] possibleSolutionHexas = new List<int>[possibleSolutions];
            //element ids twn hexas pou tha einai boundary

            int previousDivider = 1;
            int choices = 2;

            foreach (int shellId in connectedShellElementsLists)
            {
                for (int i1 = 0; i1 < solutionVectorSize; i1++)
                {
                    int SubregionsSize = solutionVectorSize / previousDivider;
                    int subregionsSeparationSize = SubregionsSize / EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
                }
            }

            //var solutionhexaElementsIds = 
            int[] hexaIds = new int[2] { 0, 1 };
            var unique = hexaIds.Distinct().Count();

            int maxSubdElements = 0;
            foreach (var subdmNhexas in EmbeddedElementsHostSubdomainsAndElements.Values)
            {
                foreach(var hexaList in subdmNhexas.Values)
                {

                }
            }

        }
    }


}
