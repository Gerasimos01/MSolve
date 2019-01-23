﻿using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;
using System;
using System.IO;

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
            int[][] subdomainsAndElements = new int[Math.Max(subdomainsAndElements1.GetLength(0), subdomainsAndElements2.GetLength(0))][];

            if (subdomainsAndElements1.GetLength(0) == subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }
            }

            if (subdomainsAndElements1.GetLength(0) < subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }

                for (int i1 = subdomainsAndElements1.GetLength(0); i1 < subdomainsAndElements2.GetLength(0); i1++)
                {
                    subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                    for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                }
            }

            if (subdomainsAndElements1.GetLength(0) < subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements2.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }

                for (int i1 = subdomainsAndElements2.GetLength(0); i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                    for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                }
            }


            return subdomainsAndElements;
        }

        public static int[][] CombineSubdomainElementsIdsArraysIntoOneCopy(int[][] subdomainsAndElements1, int[][] subdomainsAndElements2)
        {
            // input example: int[][] subdomainsAndHexas, int[][] subdomainsAndEmbedded
            int[][] subdomainsAndElements = new int[subdomainsAndElements1.GetLength(0)][];

            for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
            {
                if (subdomainsAndElements1[i1] == null)
                {
                    if (!(subdomainsAndElements2[i1] == null))
                    {
                        subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                        for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                    }
                }
                if (subdomainsAndElements2[i1] == null)
                {
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                        for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                    }
                }
                if (!(subdomainsAndElements1[i1] == null))
                {
                    if (!(subdomainsAndElements2[i1] == null))
                    {
                        subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                    }
                }

            }
            return subdomainsAndElements;
        }

        public static int[][] DetermineCoheiveELementsSubdomains(Model model, int totalSubdomains)
        {
            (Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
               Dictionary<int, List<int>> hexaConnectsShells,
               Dictionary<int, List<int>> AssignedSubdomains) =
               DdmCalculationsPartb.FindEmbeddedElementsSubdomainsCorrected(model, totalSubdomains);

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

                for (int i2 = 0; i2 < connectedShellElementsLists.ElementAt(i1).Value.Count; i2++)
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

        public static void PrintSubdomainDataForPostPro(int[][] subdHexaIds, int[][] subdCohElementIds, int[][] subdShellElementIds,string generalPath)
        {
            string hexaPath = generalPath+ @"\subdomainHexas.txt";
            string cohPath = generalPath + @"\subdomainCohesiveElements.txt";
            string shellPath = generalPath + @"\subdomainShellElements.txt";

            #region hexa elements
            int hexaPrintLength = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if(subdHexaIds[i1]==null)
                {
                    hexaPrintLength += 1;
                }
                else
                {
                    hexaPrintLength += 1 + subdHexaIds[i1].Length;
                }
            }

            var hexaPrint = new int[hexaPrintLength];
            int hexaPrintCounter = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if (subdHexaIds[i1] == null)
                {
                    hexaPrint[hexaPrintCounter] = 0;
                    hexaPrintCounter += 1;
                }
                else
                {
                    hexaPrint[hexaPrintCounter] = subdHexaIds[i1].Length;
                    hexaPrintCounter += 1;
                    for (int i2 = 0; i2 < subdHexaIds[i1].Length; i2++)
                    {
                        hexaPrint[hexaPrintCounter] = subdHexaIds[i1][i2];
                        hexaPrintCounter += 1;
                    }                    
                }
            }
            WriteToFileVector(hexaPrint, hexaPath);
            #endregion

            #region cohesive elements
            int cohePrintLength = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrintLength += 1;
                }
                else
                {
                    cohePrintLength += 1 + subdCohElementIds[i1].Length;
                }
            }

            var cohePrint = new int[cohePrintLength];
            int cohePrintCounter = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrint[cohePrintCounter] = 0;
                    cohePrintCounter += 1;
                }
                else
                {
                    cohePrint[cohePrintCounter] = subdCohElementIds[i1].Length;
                    cohePrintCounter += 1;
                    for (int i2 = 0; i2 < subdCohElementIds[i1].Length; i2++)
                    {
                        cohePrint[cohePrintCounter] = subdCohElementIds[i1][i2];
                        cohePrintCounter += 1;
                    }
                }
            }
            WriteToFileVector(cohePrint, cohPath);
            #endregion

            #region shell elements
            int shellPrintLength = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrintLength += 1;
                }
                else
                {
                    shellPrintLength += 1 + subdShellElementIds[i1].Length;
                }
            }

            var shellPrint = new int[shellPrintLength];
            int shellPrintCounter = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrint[shellPrintCounter] = 0;
                    shellPrintCounter += 1;
                }
                else
                {
                    shellPrint[shellPrintCounter] = subdShellElementIds[i1].Length;
                    shellPrintCounter += 1;
                    for (int i2 = 0; i2 < subdShellElementIds[i1].Length; i2++)
                    {
                        shellPrint[shellPrintCounter] = subdShellElementIds[i1][i2];
                        shellPrintCounter += 1;
                    }
                }
            }
            WriteToFileVector(shellPrint, shellPath);
            #endregion


        }

        public static void WriteToFileVector(int[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

        }
    }
}
