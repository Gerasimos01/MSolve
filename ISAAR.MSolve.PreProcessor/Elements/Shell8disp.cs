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
        public static int endeixiGaussCoordinates = 1;
        public static int gp_d1 { get; set; }
        public static int gp_d2 { get; set; }
        public static int gp_d3 { get; set; }
        private static int nGaussPoints;

        private double[][] gausscoordinates;

        private static double ksi;
        private static double heta;
        private static double zeta;
        private static int npoint;
        private double [] [] GetGaussCoordinates()
        {
            if(endeixiGaussCoordinates == 1)               
            {   nGaussPoints = gp_d1 * gp_d2 * gp_d3;
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
    }
}
