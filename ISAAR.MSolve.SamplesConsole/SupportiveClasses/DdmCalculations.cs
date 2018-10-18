using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor.Elements;

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

            string print_path_gen = @"C:\Users\Dimitris Giovanis\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}globalDOFs.txt";
            string file_no = subdomainID.ToString();
            string print_path = string.Format(print_path_gen, file_no);
            Vector globalDOFS = new Vector(globalDOFs);
            globalDOFS.WriteToFile(print_path);
        }
    }
}
