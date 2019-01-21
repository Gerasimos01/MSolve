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
    public static class DdmCalculationsGeneral
    {
        public static void BuildModelInterconnectionData(Model model)
        {
            //private void BuildSubdomainOfEachElement()
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.ElementsDictionary.Values)
                { element.Subdomain = subdomain; }
            }

            //private void BuildElementDictionaryOfEachNode()            
            foreach (Element element in model.ElementsDictionary.Values)
            {
                foreach (Node node in element.Nodes)
                { node.ElementsDictionary.Add(element.ID, element); }
            }

            foreach (Node node in model.NodesDictionary.Values)
            { node.BuildSubdomainDictionary(); }

            

            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            { subdomain.BuildNodesDictionary(); }
        }

        public static void UndoModelInterconnectionDataBuild(Model model)
        {
            //private void BuildSubdomainOfEachElement()
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.ElementsDictionary.Values)
                { element.Subdomain = null; }  // subdomain; }
            }

            //private void BuildElementDictionaryOfEachNode()            
            //foreach (Element element in model.ElementsDictionary.Values)
            //{
            //    foreach (Node node in element.Nodes)
            //    { node.ElementsDictionary.Add(element.ID, element); }
            //}
            foreach (Node node in model.NodesDictionary.Values)
            {
                node.ElementsDictionary.Clear();
                node.SubdomainsDictionary.Clear();
            }
            //foreach (Node node in model.NodesDictionary.Values)
            //{ node.BuildSubdomainDictionary(); }



            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                //subdomain.BuildNodesDictionary();
                subdomain.NodesDictionary.Clear();
            }
        }

        public static int[][] DetermineHexaElementsSubdomainsFromModel(Model model)
        {
            int[][] subdomainsAndHexas = new int[model.Subdomains.Count()][];

            for(int subdomainId=0; subdomainId < model.Subdomains.Count(); subdomainId++ )
            {
                var subdomain = model.Subdomains[subdomainId]; //ZERo based model.subdomainsDictionary access == model.Subdomains access
                subdomainsAndHexas[subdomainId] = new int[subdomain.ElementsDictionary.Count()];
                int hexaPositionInArray = 0;
                foreach(Element element in subdomain.ElementsDictionary.Values)
                {
                    subdomainsAndHexas[subdomainId][hexaPositionInArray] = element.ID;
                    hexaPositionInArray++;
                }
            }

            return subdomainsAndHexas;
        }

        public static int[][] CombineSubdomainElementsIdsArraysIntoOne(int[][] subdomainsAndElements1, int[][] subdomainsAndElements2)
        {
            // input example: int[][] subdomainsAndHexas, int[][] subdomainsAndEmbedded
            int[][] subdomainsAndElements = new int[subdomainsAndElements1.GetLength(0)][];

            for(int i1=0; i1< subdomainsAndElements1.GetLength(0);i1++)
            {
                subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
            }
            return subdomainsAndElements;
        }

        public static int[][] DetermineCoheiveELementsSubdomains(Model model, int totalSubdomains)
        {
            (Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
               Dictionary<int, List<int>> hexaConnectsShells,
               Dictionary<int, List<int>> AssignedSubdomains) =
               DdmCalculationsPartb.FindEmbeddedElementsSubdomains(model, totalSubdomains);

            Dictionary<int, List<int>> connectedShellElementsLists = 
                DdmCalculations.DetermineOnlyNeededCombinations(
            AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
             hexaConnectsShells);

            int numlistsForCombinations = connectedShellElementsLists.Keys.Count();
            int[][] combinationSolutions = new int[numlistsForCombinations][];
            for (int i1=0; i1<numlistsForCombinations;i1++)
            {
                combinationSolutions[i1] =
                    DdmCalculations.CalculateCombinationSolution2(connectedShellElementsLists.ElementAt(i1).Value,
                    AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem);
            }

            //gather solutions 
            var subdElementsNums = new int[totalSubdomains];
            for (int i1 = 0; i1 < numlistsForCombinations; i1++)
            {
                for (int i2 = 0; i2 < combinationSolutions[i1].Length; i2++)
                {
                    subdElementsNums[combinationSolutions[i1][i2]] += 1;
                }
            }

            int[][] subdAmbCohElementIds = new int[totalSubdomains][];
            for(int i1 = 0; i1 < totalSubdomains; i1++)
            {
                subdAmbCohElementIds[i1] = new int[subdElementsNums[i1]];
            }

            int[] subdCohElemCounters = new int[totalSubdomains];
            for(int i1 = 0; i1 < numlistsForCombinations; i1++)
            {
                var cohIDsList = connectedShellElementsLists.ElementAt(i1).Value;
                var solutionSubdomainIDs = combinationSolutions[i1];

                for (int i2 = 0; i2 < numlistsForCombinations; i2++)
                {
                    int shellID = cohIDsList[i2];
                    int subdId = solutionSubdomainIDs[i2];
                    subdAmbCohElementIds[subdId][subdCohElemCounters[subdId]] = shellID;
                    subdCohElemCounters[subdId] += 1;
                }
            }

            int[][] subdCohElementIdsDirect = DdmCalculationsPartb.ConvertIntListToArray(AssignedSubdomains);
            int[][] subdCohElementIds = CombineSubdomainElementsIdsArraysIntoOne(subdAmbCohElementIds, subdCohElementIdsDirect);

            return subdCohElementIds;
        }

        public static int[][] DetermineShellELementsSubdomains(Model model, int totalSubdomains, int[][] subdCohElementIds,
            int[] lowerCohesiveBound, int[] upperCohesiveBound,int[] grShElementssnumber)
        {
            List<int>[] subdShellElementIds = new List<int>[totalSubdomains];
            for (int i1 = 0; i1 < totalSubdomains; i1++)
            {
                subdShellElementIds[i1] = new List<int>(subdCohElementIds[i1].Length);
            }

            for (int i1 = 0; i1 < totalSubdomains; i1++)
            {
                for (int i2 = 0; i2 < subdCohElementIds[i1].Length; i2++)
                {
                    int cohID = subdCohElementIds[i1][i2];
                    for (int i3 = 0; i3 < lowerCohesiveBound.Length; i3++)
                    {
                        if ((lowerCohesiveBound[i3]<=cohID)&(upperCohesiveBound[i3]>=cohID))
                        {
                            //subdID=i1;
                            if (!subdShellElementIds[i1].Contains(cohID - grShElementssnumber[i3]))
                            { subdShellElementIds[i1].Add(cohID - grShElementssnumber[i3]); }
                            break;
                        }

                    }
                }
            }

            int[][] subdShellElementIdsArrays = new int[subdShellElementIds.Length][];
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                subdShellElementIdsArrays[i1] = subdShellElementIds[i1].ToArray();
            }

            return subdShellElementIdsArrays;
        }


    }
}
