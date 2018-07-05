using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;

namespace ISAAR.MSolve.FEM
{
    public static class SubdomainCalculations
    {
        public static double[][] CalculateKfreeprescribedDqMultiplications(Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            //ElementStructuralStiffnessProvider s = new ElementStructuralStiffnessProvider();
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            var keys1 = boundaryNodes.Keys;
            int[] keys2=new int [2];
            keys1.CopyTo(keys2, 0);


            double[] macroVariable = scaleTransitions.MicroToMacroTransition(boundaryNodes[keys2[0]],new double [3]);

            double[][] KfpDqVectors = new double[macroVariable.GetLength(0)][];
            for (int j1 = 0; j1 < macroVariable.GetLength(0); j1++)
            {
                KfpDqVectors[j1] = new double[subdomain.TotalDOFs]; // h allliws subdomain.Forces.GetLength(0)
            }

            // TODO: should encapsulate DOF logic into a separate entity that will manage things if embedded or not (should return element matrix and globaldofs correspondence list
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            //SkylineMatrix2D K = new SkylineMatrix2D(GlobalMatrixAssemblerSkyline.CalculateRowIndex(subdomain, nodalDOFsDictionary)); // commented out den kanoume assembly afto
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);

            foreach (Element element in subdomain.ElementsDictionary.Values) // TODOGerasimos edw mporei na xrhsimopoihthei to dictionary twn eleement pou exoun fp nodes
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    Node nodeRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                        {                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                Node nodeColumn = matrixAssemblyNodes[j];
                                //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                //{
                                //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                //    if (dofColumn != -1)
                                //    {
                                //        int height = dofRow - dofColumn;
                                //        if (height >= 0)
                                //            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
                                //    }
                                //    iElementMatrixColumn++;
                                //}
                                int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
                                if (boundaryNodes.ContainsKey(nodeColumn.ID)) 
                                {
                                    double[] element_Kfp_triplette = new double[nodalDofsNumber];
                                    for (int j1 = 0; j1 < nodalDofsNumber; j1++)
                                    {
                                       element_Kfp_triplette[j1]= ElementK[iElementMatrixRow, iElementMatrixColumn+j1];
                                    }

                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kfp_triplette);
                                    for (int j2 = 0; j2 < contribution.GetLength(0); j2++)
                                    {
                                        KfpDqVectors[j2][iElementMatrixRow] += contribution[j2];
                                    }

                                }
                                iElementMatrixColumn += nodalDofsNumber;

                            }
                        }
                        iElementMatrixRow++;
                    }
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return KfpDqVectors;
        }

        public static double[][] CalculateKpfKffinverseKfpDq(double[][] f2_vectors, Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;
            double[][] f3_vectors = new double[f2_vectors.GetLength(0)][];
            for (int i1 = 0; i1 < f2_vectors.GetLength(0); i1++)
            {
                f3_vectors[i1] = new double[scaleTransitions.PrescribedDofsPerNode() * boundaryNodes.Count];
            }
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    Node nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                Node nodeColumn = matrixAssemblyNodes[j];
                                foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                {
                                    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (dofColumn != -1)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                                        
                                        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                        {
                                            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow+i1, iElementMatrixColumn]*f2_vectors[i2][dofColumn];
                                        }

                                    }
                                    iElementMatrixColumn++;
                                }
                            }

                        }
                    }
                    iElementMatrixRow+= elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return f3_vectors;

        }


        public static Dictionary<int,int> GetNodesOrderInDictionary(Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, int> boundaryNodesOrder = new Dictionary<int, int>;
            int order = 1;
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                boundaryNodesOrder.Add(boundaryNode.ID, order);
                order += 1;
            }
            return boundaryNodesOrder;

        }

    }
}
