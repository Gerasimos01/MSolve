using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    class Shell8disp //: IStructuralFiniteElement
    { 
        public static int gp_d1 { get; set; }
        public static int gp_d2 { get; set; }
        public static int gp_d3 { get; set; }
        public static int [] tk { get; set; }
        private static int nGaussPoints;

        private static double ksi;
        private static double heta;
        private static double zeta;
        private static int npoint;

        public static int endeixiGaussCoordinates = 1;
        private double[][] gausscoordinates;
        private double [] [] GetGaussCoordinates() //3 dianysmata me tis timew tvn ksi heta zeta se ola ta gauss points
        {
            if(endeixiGaussCoordinates == 1)               
            {   nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                gausscoordinates = new double[3][];
                for (int l = 0; l < 3; l++)
                { gausscoordinates[l] = new double[nGaussPoints]; }
                for (int l = 0; l < gp_d3; l++)
                { for (int k = 0; k < gp_d2; k++)
                    {for (int j = 0; j < gp_d1; j++)
                        {
                            npoint = l * (gp_d1 * gp_d2) + k * gp_d1 + j;
                            if (gp_d1 == 3) { ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483); }
                            if (gp_d1 == 2) { ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1); }
                            if (gp_d2 == 3) { heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483); }
                            if (gp_d2 == 2) { heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1); }
                            if (gp_d3 == 3) { zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483); }
                            if (gp_d3 == 2) { zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1); }
                            gausscoordinates[0][npoint] = ksi;
                            gausscoordinates[1][npoint] = heta;
                            gausscoordinates[2][npoint] = zeta;
                        }
                    }
                }
                endeixiGaussCoordinates = 2;
                return gausscoordinates;    
            }
            else
            { return gausscoordinates; }
        }

        private double[][] shapeFunctions;
        public static int endeixiShapeFunctions = 1;
        private double [] [] GetShapeFunctions() // 8 dianusmata me tis times twn N1....N8 se kathe gauss point
        {
            if (endeixiShapeFunctions == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                gausscoordinates = this.GetGaussCoordinates();
                shapeFunctions = new double[8][];
                for (int j = 0; j < 8; j++)
                { shapeFunctions[j] = new double[nGaussPoints]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    shapeFunctions[4][j] = 0.5 * (1 - Math.Pow(gausscoordinates[0][j], 2)) * (1 + gausscoordinates[1][j]);
                    shapeFunctions[5][j] = 0.5 * (1 - Math.Pow(gausscoordinates[1][j], 2)) * (1 - gausscoordinates[0][j]);
                    shapeFunctions[6][j] = 0.5 * (1 - Math.Pow(gausscoordinates[0][j], 2)) * (1 - gausscoordinates[1][j]);
                    shapeFunctions[7][j] = 0.5 * (1 - Math.Pow(gausscoordinates[1][j], 2)) * (1 + gausscoordinates[0][j]);
                    shapeFunctions[0][j] = 0.25 * (1 + gausscoordinates[0][j]) * (1 + gausscoordinates[1][j]) - 0.5 * shapeFunctions[4][j] - 0.5 * shapeFunctions[7][j];
                    shapeFunctions[1][j] = 0.25 * (1 - gausscoordinates[0][j]) * (1 + gausscoordinates[1][j]) - 0.5 * shapeFunctions[4][j] - 0.5 * shapeFunctions[5][j];
                    shapeFunctions[2][j] = 0.25 * (1 - gausscoordinates[0][j]) * (1 - gausscoordinates[1][j]) - 0.5 * shapeFunctions[5][j] - 0.5 * shapeFunctions[6][j];
                    shapeFunctions[3][j] = 0.25 * (1 + gausscoordinates[0][j]) * (1 - gausscoordinates[1][j]) - 0.5 * shapeFunctions[6][j] - 0.5 * shapeFunctions[7][j];
                }
                endeixiShapeFunctions = 2;
                return shapeFunctions;
            }
            else
            { return shapeFunctions; }
        }

        private double[][] shapeFunctionDerivatives;
        public static int endeixiShapeFunctionDerivatives = 1;
        private double[][] GetShapeFunctionDerivatives() // 16 dianusmata me tis times twn N1ksi....N8ksi,N1heta,....N8heta se kathe gauss point
        {
            if (endeixiShapeFunctionDerivatives == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                gausscoordinates = this.GetGaussCoordinates();
                shapeFunctionDerivatives = new double[16][];
                for (int j = 0; j < 16; j++)
                { shapeFunctionDerivatives[j] = new double[nGaussPoints]; }
                for (int j = 0; j < nGaussPoints; j++)
                {   
                    //Ni_ksi
                    shapeFunctionDerivatives[4][j] = (-gausscoordinates[0][j]) * (1 + gausscoordinates[1][j]);
                    shapeFunctionDerivatives[5][j] = -0.5 * (1 - Math.Pow(gausscoordinates[1][j], 2));
                    shapeFunctionDerivatives[6][j] = 0.5 * (-2 * gausscoordinates[0][j]) * (1 - gausscoordinates[1][j]);
                    shapeFunctionDerivatives[7][j] = 0.5 * (1 - Math.Pow(gausscoordinates[1][j], 2));
                    shapeFunctionDerivatives[0][j] = +0.25 * (1 + gausscoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[4][j] - 0.5 * shapeFunctionDerivatives[7][j];
                    shapeFunctionDerivatives[1][j] = -0.25 * (1 + gausscoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[4][j] - 0.5 * shapeFunctionDerivatives[5][j];
                    shapeFunctionDerivatives[2][j] = -0.25 * (1 - gausscoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[5][j] - 0.5 * shapeFunctionDerivatives[6][j];
                    shapeFunctionDerivatives[3][j] = +0.25 * (1 - gausscoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[6][j] - 0.5 * shapeFunctionDerivatives[7][j];
                    //Ni_heta
                    shapeFunctionDerivatives[12][j] = 0.5 * (1 - Math.Pow(gausscoordinates[0][j], 2));
                    shapeFunctionDerivatives[13][j] = 0.5 * (-2 * gausscoordinates[1][j]) * (1 - gausscoordinates[0][j]);
                    shapeFunctionDerivatives[14][j] = 0.5 * (1 - Math.Pow(gausscoordinates[0][j], 2)) * (-1);
                    shapeFunctionDerivatives[15][j] = 0.5 * (-2 * gausscoordinates[1][j]) * (1 + gausscoordinates[0][j]);
                    shapeFunctionDerivatives[8][j] = +0.25 * (1 + gausscoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[12][j] - 0.5 * shapeFunctionDerivatives[15][j];
                    shapeFunctionDerivatives[9][j] = +0.25 * (1 - gausscoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[12][j] - 0.5 * shapeFunctionDerivatives[13][j];
                    shapeFunctionDerivatives[10][j] = -0.25 * (1 - gausscoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[13][j] - 0.5 * shapeFunctionDerivatives[14][j];
                    shapeFunctionDerivatives[11][j] = -0.25 * (1 + gausscoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[14][j] - 0.5 * shapeFunctionDerivatives[15][j];
                }
                endeixiShapeFunctionDerivatives = 2;
                return shapeFunctionDerivatives;
            }
            else
            { return shapeFunctionDerivatives; }
        }

        private double[][,] ll1;
        public static int endeixill1 = 1;
        private double[][,] Getll1() //einai teliko kai oxi prok
        {
            if (endeixill1 == 1)
            {
                gausscoordinates = this.GetGaussCoordinates();
                shapeFunctionDerivatives = this.GetShapeFunctionDerivatives();
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;

                ll1 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { ll1[j] = new double[3, 24]; }
                for (int j = 0; j < nGaussPoints; j++) //dhmiourgia olklhrou tou ll1 gia kathe gauss point
                {
                    for (int k = 0; k < 8; k++)
                    {
                        ll1[j][0, 3 * k] = shapeFunctionDerivatives[k][j];
                        ll1[j][0, 3 * k + 1] = 0.5 * gausscoordinates[2][j] * tk[k] * shapeFunctionDerivatives[k][j];
                        ll1[j][0, 3 * k + 2] = -ll1[j][0, 3 * k + 1];
                        ll1[j][1, 3 * k] = shapeFunctionDerivatives[k + 8][j];
                        ll1[j][1, 3 * k + 1] = 0.5 * gausscoordinates[2][j] * tk[k] * shapeFunctionDerivatives[k + 8][j];
                        ll1[j][1, 3 * k + 2] = -ll1[j][1, 3 * k + 1];
                        ll1[j][2, 3 * k] = 0;
                        ll1[j][2, 3 * k + 1] = 0.5 * tk[k] * shapeFunctions[k][j];
                        ll1[j][2, 3 * k + 2] = -ll1[j][2, 3 * k + 1];
                    }

                }
                endeixill1 = 2;
                return ll1;
            }
            else
            { return ll1; }
          }


        private double[][,] J_0a;
        public static int endeixiJ_0a = 1;
        private double[][,] GetJ_0a() //einai teliko kai oxi prok
        {
            if (endeixiJ_0a == 1)
            {
                ll1 = this.Getll1();
                for (int j = 0; j < nGaussPoints; j++)
                { J_0a[j] = new double[3, 16]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 8; k++)
                    {
                        J_0a[j][0, 2 * k] = ll1[j][0, 3 * k];
                        J_0a[j][0, 2 * k + 1] = ll1[j][0, 3 * k + 1];
                        J_0a[j][1, 2 * k] = ll1[j][1, 3 * k];
                        J_0a[j][1, 2 * k + 1] = ll1[j][1, 3 * k + 1];
                        J_0a[j][2, 2 * k] = ll1[j][2, 3 * k];
                        J_0a[j][2, 2 * k + 1] = ll1[j][2, 3 * k + 1];
                    }
                }
                endeixiJ_0a = 2;
                return J_0a;
            }
            else
            { return J_0a; }
        }





    }
}

