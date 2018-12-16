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

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    public static class DdmCalculationsPartb
    {
        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomains(Model model,int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            //1
            Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
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
                    { AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem.Add(element.ID, HostSubdomains); }
                    if (HostSubdomains.Count == 1) 
                    {
                        if ( AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key) )
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
            return (AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, hexaConnectsShells,AssignedSubdomains);
        }

        public static int[][] DetermineAmbiguousSimple(Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, Dictionary<int, List<int>> hexaConnectsShells, Dictionary<int, List<int>> AssignedSubdomains)
        {
            Dictionary<int, Dictionary<int, IList<int>>> AmbigElements = AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem;

            foreach (int ambElementID in AmbigElements.Keys)
            {
                int[] numSpecificElements = new int[AmbigElements[ambElementID].Keys.Count];

                int subdCounter = 0;
                foreach(int subdID in AmbigElements[ambElementID].Keys)
                {
                    //numSpecificElements[subdCounter] = AmbigElements[ambElementID][subdID].Count; // Count dld ta element ths inner listas tou AmbigElements...
                    numSpecificElements[subdCounter] = AmbigElements[ambElementID].ElementAt(subdCounter).Value.Count; // Value einai h inner lista tou AmbigElements...
                    subdCounter++;
                }

                int thesiChosenSubd = 0;
                int nElements = 0;
                for (int i1=0; i1<numSpecificElements.GetLength(0); i1++ )
                {
                    if (numSpecificElements[i1]>nElements)
                    {
                        thesiChosenSubd = i1;
                        nElements = numSpecificElements[i1];
                    }
                }

                int chosenSubdomainID = AmbigElements[ambElementID].ElementAt(thesiChosenSubd).Key;

                if (AssignedSubdomains.Keys.Contains(chosenSubdomainID))
                {
                    AssignedSubdomains[chosenSubdomainID].Add(ambElementID);
                }
                else
                {
                    List<int> subdElementsIds = new List<int>();
                    subdElementsIds.Add(ambElementID);
                    AssignedSubdomains.Add(chosenSubdomainID, subdElementsIds);
                }

            }

            return ConvertIntListToArray(AssignedSubdomains);
        }

        public static int[][] ConvertIntListToArray(Dictionary<int, List<int>> AssignedSubdomains)
        {
            int maxSubdId = AssignedSubdomains.Keys.Max();
            int[][] subdIdsAndElements = new int[maxSubdId][];

            foreach (int subdID in AssignedSubdomains.Keys)
            {
                subdIdsAndElements[subdID] = AssignedSubdomains[subdID].ToArray();
            }

            return subdIdsAndElements;
        }
    }
}
