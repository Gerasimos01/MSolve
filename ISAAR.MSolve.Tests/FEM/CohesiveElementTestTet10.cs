using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.SamplesConsole
{
    public static class CohesiveElementTestTet10
    {
        public static void BuildTet10Interface(Model model)
        {
            BenzeggaghKenaneCohesiveMaterial material2 = new BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = 57,// N / mm2
                D_o_3 = 0.000057, // mm
                D_f_3 = 0.0098245610,  // mm

                T_o_1 = 57,// N / mm2
                D_o_1 = 0.000057, // mm
                D_f_1 = 0.0098245610,  // mm

                n_curve = 1.4
            };


            int[] nodeIds = new int[] { 5, 6, 7, 8, 17, 18, 19, 20, 38, 39, 40, 41 };

            double[,] nodeData = new double[,]
                {{-0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,1.0000000000000000},

                {0.0000000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,1.0000000000000000},

                {0.0000000000000000,0.0000000000000000,1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,1.0000000000000000} };

            int[,] elementData = new int[,]
            {
                    {6,5,38,17,39,40 },
                    {5,8,38,20,42,39 },
                    {7,6,38,18,40,41 },
                    {8,7,38,19,41,42 }
            };

            int[] upperSideSuplicateNodeIds = nodeIds.Select(x => x + 100).ToArray();

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeIds.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nodeIds[nNode], new Node(id: nodeIds[nNode], x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }
            for (int nNode = 0; nNode < upperSideSuplicateNodeIds.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(upperSideSuplicateNodeIds[nNode], new Node(id: upperSideSuplicateNodeIds[nNode], x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = 3,
                    ElementType = new cohesiveElement(material2, GaussLegendre2D.GetQuadratureWithOrder(3, 3), InterpolationQuad4.UniqueInstance) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
                for (int j = 0; j < 6; j++)
                {
                    e1.NodesDictionary.Add( elementData[nElement,j]+100, model.NodesDictionary[elementData[nElement, j]+100]);
                }
                for (int j = 0; j < 6; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j], model.NodesDictionary[elementData[nElement, j]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1.ID, e1);
            }

            
        }
        public static void BuildModel(Model model)
        {
            int[] nodeIds = new int[] { 5, 6, 7, 8, 17, 18, 19, 20, 38, 39, 40, 41 };

            double[,] nodeData = new double[,]
                {{-0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,-0.2500000000000000,1.0000000000000000},

                {0.0000000000000000,0.2500000000000000,1.0000000000000000},
                {0.2500000000000000,0.0000000000000000,1.0000000000000000},
                {0.0000000000000000,-0.2500000000000000,1.0000000000000000},
                {-0.2500000000000000,0.0000000000000000,1.0000000000000000},

                {0.0000000000000000,0.0000000000000000,1.0000000000000000},
                {-0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,0.1250000000000000,1.0000000000000000},
                {0.1250000000000000,-0.1250000000000000,1.0000000000000000} };

            int[,] elementData = new int[,]
            {
                    {6,5,38,17,39,40 },
                    {5,8,38,20,42,39 },
                    {7,6,38,18,40,41 },
                    {8,7,38,19,41,42 }
            };

            int[] upperSideSuplicateNodeIds = nodeIds.Select(x => x + 100).ToArray();



        }

        

    }
}
