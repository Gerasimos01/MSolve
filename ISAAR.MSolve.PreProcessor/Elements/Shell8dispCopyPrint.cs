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
    public class Shell8dispCopyPrint : IStructuralFiniteElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IFiniteElementMaterial3D[] materialsAtGaussPoints;
        protected IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw 

        public double[][] oVn_i { get; set; } 
        public double[][] oV1_i { get; set; }
        //public double[][] oV2_i { get; set; }
        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths
        public  int gp_d1  { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public  int gp_d2  { get; set; }
        public  int gp_d3  { get; set; }
        public double [] tk { get; set; } //public static int[] tk { get; set; }
        private  int nGaussPoints;

        private  double ksi;
        private  double heta;
        private  double zeta;
        private  int npoint;

        private double[] a_123g;
        private double a_1g;
        private double a_2g;
        private double a_3g;
        
        protected Shell8dispCopyPrint()//consztructor apo to hexa8
        {
        }

        public Shell8dispCopyPrint(IFiniteElementMaterial3D material,int gp_d1c,int gp_d2c,int gp_d3c)
        {
            this.gp_d1 = gp_d1c;
            this.gp_d2 = gp_d2c;
            this.gp_d3 = gp_d3c;
            this.nGaussPoints = this.gp_d1 * this.gp_d2 * this.gp_d3;
            materialsAtGaussPoints = new IFiniteElementMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IFiniteElementMaterial3D)material.Clone();

        }

        //public Shell8disp(IFiniteElementMaterial3D material, IFiniteElementDOFEnumerator dofEnumerator)//pithanotata den xreiazetai
        //    : this(material, gp_d1, gp_d2, gp_d3)
        //{
        //    this.dofEnumerator = dofEnumerator;
        //}
        // ews edw

        public int endeixiGaussCoordinates = 1;
        private double[][] gausscoordinates;
        private double[][] GetGaussCoordinates() //3 dianysmata me tis timew tvn ksi heta zeta se ola ta gauss points
        {
            if (endeixiGaussCoordinates == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                a_123g = new double[nGaussPoints];
                gausscoordinates = new double[3][];
                for (int l = 0; l < 3; l++)
                { gausscoordinates[l] = new double[nGaussPoints]; }
                for (int l = 0; l < gp_d3; l++)
                {
                    for (int k = 0; k < gp_d2; k++)
                    {
                        for (int j = 0; j < gp_d1; j++)
                        {
                            npoint = l * (gp_d1 * gp_d2) + k * gp_d1 + j;
                            if (gp_d1 == 3)
                            {
                                ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                                a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                            }
                            if (gp_d1 == 2)
                            {
                                ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                                a_1g = 1;
                            }                                            
                            if (gp_d2 == 3)
                            {
                                heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                                a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                            }
                            if (gp_d2 == 2)
                            {
                                heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                                a_2g = 1;
                            }
                            if (gp_d3 == 3)
                            {
                                zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483);
                                a_3g = 0.5 * (l - 1) * (l - 2) * (0.555555555555555) + (-1) * (l) * (l - 2) * (0.888888888888888) + 0.5 * (l) * (l - 1) * (0.555555555555555);
                            }
                            if (gp_d3 == 2)
                            {
                                zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1);
                                a_3g = 1;
                            }
                            gausscoordinates[0][npoint] = ksi;
                            gausscoordinates[1][npoint] = heta;
                            gausscoordinates[2][npoint] = zeta;

                            a_123g[npoint] = a_1g * a_2g * a_3g;
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
        public  int endeixiShapeFunctions = 1;
        private double[][] GetShapeFunctions() // 8 dianusmata me tis times twn N1....N8 se kathe gauss point
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
        public  int endeixiShapeFunctionDerivatives = 1;
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
        public  int endeixill1 = 1;
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
        public  int endeixiJ_0a = 1;
        private double[][,] GetJ_0a() //einai teliko kai oxi prok
        {
            if (endeixiJ_0a == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                J_0a = new double[nGaussPoints][,];
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

        private void GetInitialGeometricData(Element element) //TODO mhpws me endeixiInitialGeometricD...
        {
            ox_i = new double[8][];
            tx_i = new double[8][];
            tU = new double[8][];
            tUvec = new double[8][];
            oV1_i = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                tx_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                tU[j] = new double[6];
                tUvec[j] = new double[6];
                oV1_i[j] = new double[3];
                for (int k = 0; k < 3; k++) { tU[j][3 + k] = oVn_i[j][k]; }

                tUvec[j][0] = tU[j][5];
                tUvec[j][1] = 0;
                tUvec[j][2] = -tU[j][3];

                tV1norm = Math.Sqrt(tUvec[j][0] * tUvec[j][0] + tUvec[j][1] * tUvec[j][1] + tUvec[j][2] * tUvec[j][2]);

                tUvec[j][0] = tUvec[j][0] / tV1norm;
                tUvec[j][1] = tUvec[j][1] / tV1norm;
                tUvec[j][2] = tUvec[j][2] / tV1norm;

                oV1_i[j][0] = tUvec[j][0];
                oV1_i[j][1] = tUvec[j][1];
                oV1_i[j][2] = tUvec[j][2];

                tUvec[j][3] = tU[j][3 + 1] * tUvec[j][2] - tU[j][3 + 2] * tUvec[j][1];
                tUvec[j][4] = tU[j][3 + 2] * tUvec[j][0] - tU[j][3 + 0] * tUvec[j][2];
                tUvec[j][5] = tU[j][3 + 0] * tUvec[j][1] - tU[j][3 + 1] * tUvec[j][0];
            }

        }

        private double[,] J_0b;    //einai idio gia ola ta gauss points
        public int endeixiJ_0b = 1;
        private double[,] GetJ_0b(Element element)
        {
            if (endeixiJ_0b == 1)
            {
                J_0b = new double[16, 3];
                this.GetInitialGeometricData(element);
                for (int j = 0; j < 8; j++)
                {
                    J_0b[2 * j, 0] = ox_i[j][0];
                    J_0b[2 * j + 1, 0] = this.oVn_i[j][0];
                    J_0b[2 * j, 1] = ox_i[j][1];
                    J_0b[2 * j + 1, 1] = this.oVn_i[j][1];
                    J_0b[2 * j, 2] = ox_i[j][2];
                    J_0b[2 * j + 1, 2] = this.oVn_i[j][2];
                }
                endeixiJ_0b = 2;
                return J_0b;
            }
            else
            { return J_0b; }
        }

        private double[][,] J_0;       //den einai to idio gia ola ta gausspoint
        public  int endeixiJ_0 = 1;
        private double[][,] GetJ_0(Element element)   // einai teliko kai oxi prok
        {
            if (endeixiJ_0 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                J_0 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { J_0[j] = new double[3, 3]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            J_0[j][k, l] = 0;
                            for (int m = 0; m < 16; m++)
                            {
                                J_0[j][k, l] += GetJ_0a()[j][k, m] * GetJ_0b(element)[m, l];
                            }

                        }

                    }
                }
                endeixiJ_0 = 2;
                return J_0;
            }
            else
            { return J_0; }

        }

        private double[] detJ_0; //[] osa kai ta gauss points
        public  int endeixiDetJ_0 = 1;
        private double[] GetDetJ_0(Element element)
        {
            if (endeixiDetJ_0 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                detJ_0 = new double[nGaussPoints];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    double det1 = GetJ_0(element)[j][0, 0] *
                         ((GetJ_0(element)[j][1, 1] * GetJ_0(element)[j][2, 2]) - (GetJ_0(element)[j][2, 1] * GetJ_0(element)[j][1, 2]));
                    double det2 = GetJ_0(element)[j][0, 1] *
                                  ((GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][2, 2]) - (GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][1, 2]));
                    double det3 = GetJ_0(element)[j][0, 2] *
                                  ((GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][2, 1]) - (GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][1, 1]));

                    double jacobianDeterminant = det1 - det2 + det3;

                    if (jacobianDeterminant < 0)
                    {
                        throw new InvalidOperationException("The Jacobian Determinant is negative.");
                    }

                    detJ_0[j] = jacobianDeterminant;
                }
                endeixiDetJ_0 = 2;
                return detJ_0;
            }
            else
            { return detJ_0; }
        }

        private double[][,] J_0inv;
        public  int endeixiJ_0inv = 1;
        private double[][,] GetJ_0inv(Element element)
        {
            if (endeixiDetJ_0 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                J_0inv = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { J_0inv[j] = new double[3, 3]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    J_0inv[j][0, 0] = ((GetJ_0(element)[j][1, 1] * GetJ_0(element)[j][2, 2]) - (GetJ_0(element)[j][2, 1] * GetJ_0(element)[j][1, 2])) *
                                    (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][0, 1] = ((GetJ_0(element)[j][2, 1] * GetJ_0(element)[j][0, 2]) - (GetJ_0(element)[j][0, 1] * GetJ_0(element)[j][2, 2])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][0, 2] = ((GetJ_0(element)[j][0, 1] * GetJ_0(element)[j][1, 2]) - (GetJ_0(element)[j][1, 1] * GetJ_0(element)[j][0, 2])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][1, 0] = ((GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][1, 2]) - (GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][2, 2])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][1, 1] = ((GetJ_0(element)[j][0, 0] * GetJ_0(element)[j][2, 2]) - (GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][0, 2])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][1, 2] = ((GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][0, 2]) - (GetJ_0(element)[j][0, 0] * GetJ_0(element)[j][1, 2])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][2, 0] = ((GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][2, 1]) - (GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][1, 1])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][2, 1] = ((GetJ_0(element)[j][2, 0] * GetJ_0(element)[j][0, 1]) - (GetJ_0(element)[j][2, 1] * GetJ_0(element)[j][0, 0])) *
                                            (1 / GetDetJ_0(element)[j]);
                    J_0inv[j][2, 2] = ((GetJ_0(element)[j][0, 0] * GetJ_0(element)[j][1, 1]) - (GetJ_0(element)[j][1, 0] * GetJ_0(element)[j][0, 1])) *
                                            (1 / GetDetJ_0(element)[j]);
                }
                endeixiJ_0inv = 2;
                return J_0inv;
            }
            else
            { return J_0inv; }
        }

        private double[][,] BL11a;
        public  int endeixiBL11a = 1;
        private double[][,] GetBL11a(Element element)
        {
            if (endeixiBL11a == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BL11a = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BL11a[j] = new double[6, 9]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 6; k++)
                    { for (int l = 0; l < 9; l++)
                        { BL11a[j][k, l] = 0; } }

                    for (int k = 0; k < 3; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL11a[j][k, 3 * k + l] = GetJ_0inv(element)[j][k, l]; } }

                    //gemisma [4,4] ews [4,6] kai [5,7] ews [5,9]
                    for (int k = 0; k < 2; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL11a[j][3 + k, 3 + 3 * k + l] = GetJ_0inv(element)[j][k, l]; } }

                    //gemisma [4,1] ews [4,3] kai [5,4] ews [5,6]
                    for (int k = 0; k < 2; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL11a[j][3 + k, 3 * k + l] = GetJ_0inv(element)[j][1 + k, l]; } }

                    for (int l = 0; l < 3; l++)
                    { BL11a[j][5, l] = GetJ_0inv(element)[j][2, l]; }

                    for (int l = 0; l < 3; l++)
                    { BL11a[j][5, 6 + l] = GetJ_0inv(element)[j][0, l]; }
                }
                endeixiBL11a = 2;
                return BL11a;
            }
            else
            { return BL11a; }
        }

        private double[][,] BL12;
        public int endeixiBL12 = 1;
        private double[][,] GetBL12(Element element)
        {
            if (endeixiBL12 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BL12 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BL12[j] = new double[9, 9]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    { for (int l = 0; l < 9; l++)
                        { BL12[j][k, l] = 0; } }

                    for (int k = 0; k < 3; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL12[j][k, 3 * k + l] = GetJ_0inv(element)[j][0, l]; } }

                    for (int k = 0; k < 3; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL12[j][3 + k, 3 * k + l] = GetJ_0inv(element)[j][1, l]; } }

                    for (int k = 0; k < 3; k++)
                    { for (int l = 0; l < 3; l++)
                        { BL12[j][6 + k, 3 * k + l] = GetJ_0inv(element)[j][2, l]; } }
                }
                endeixiBL12 = 2;
                return BL12;
            }
            else
            { return BL12; }
        }

        private double[][,] BNL1;
        public int endeixiBNL1 = 1;
        private double[][,] GetBNL1(Element element)
        {
            if (endeixiBNL1 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BNL1 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BNL1[j] = new double[9, 9]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    { for (int l = 0; l < 9; l++)
                        { BNL1[j][k, l] = 0; } }

                    for (int m = 0; m < 3; m++)
                    { for (int k = 0; k < 3; k++)
                        { for (int l = 0; l < 3; l++)
                            { BNL1[j][3 * m + k, 3 * m + l] = GetJ_0inv(element)[j][k, l]; } } }
                }
                endeixiBNL1 = 2;
                return BNL1;
            }
            else
            { return BNL1; }
        }

        // theseis metavlhtwn pou ANANEWNONTAI
        // kai kapoies apo aftes tha xreiastoun kai upol/smous GET INITIAL....
        private double[][] tx_i; //8 arrays twn 3 stoixeiwn //den einai apo afta pou orizei o xrhsths
        private double[][] tU;   //8 arrays twn 6 stoixeiwn 
        private double[][] tUvec;//8 arrays twn 6 stoixeiwn

        // methodoi dhmiourgias pinakwn pou periexoun stoixeia pou ananewnontai

        private double[,] ll2;
        public int endeixill2 = 1;
        private void Calculatell2() //meta apo enhmerwsh h initialize
        {
            if (endeixill2 == 1)
            {
                ll2 = new double[24, 3];
                for (int j = 0; j < 8; j++)
                { for (int k = 0; k < 3; k++)
                    { ll2[3 * j + 0, k] = tU[j][k];
                        ll2[3 * j + 1, k] = tU[j][3 + k];
                        ll2[3 * j + 2, k] = oVn_i[j][k]; } }

                endeixill2 = 2;
                //return ll2;
            }
            else
            {
                for (int j = 0; j < 8; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        ll2[3 * j + 0, k] = tU[j][k];
                        ll2[3 * j + 1, k] = tU[j][3 + k];
                    }
                }
                //return ll2;
            }
        }

        private double[][,] l_circumflex;
        public int endeixil_circumflex = 1;
        private void Calculatel_circumflex() //afou periexei getll2: meta apo enhmerwsh h initialize
        {
            if (endeixil_circumflex == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                l_circumflex = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { l_circumflex[j] = new double[3, 3]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            l_circumflex[j][k, l] = 0;
                            for (int m = 0; m < 24; m++)
                            {
                                l_circumflex[j][k, l] += Getll1()[j][k, m] * ll2[m, l];
                            }

                        }

                    }

                }
                endeixil_circumflex = 2;
                //return l_circumflex;
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            l_circumflex[j][k, l] = 0;
                            for (int m = 0; m < 24; m++)
                            {
                                l_circumflex[j][k, l] += Getll1()[j][k, m] * ll2[m, l];
                            }

                        }

                    }

                }
                //return l_circumflex;
            }
        }

        private double[][,] BL11b;
        public int endeixiBL11b = 1;
        private void CalculateBL11b() //afou periexei Getl_circumflex:getll2: meta apo enhmerwsh h initialize 
        {
            if (endeixiBL11b == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BL11b = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BL11b[j] = new double[9, 9]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        { BL11b[j][k, l] = 0; }
                    }


                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; }
                        }
                    }
                }
                endeixiBL11b = 2;
                //return BL11b; 
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        { BL11b[j][k, l] = 0; }
                    }


                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            { BL11b[j][3 * k + l, 3 * k + m] = l_circumflex[j][l, m]; }
                        }
                    }
                }
                //return BL11b;
            }
        }


        private double[][,] BL11;
        public int endeixiBL11 = 1;
        private void CalculateBL11(Element element) //afou periexei BL11:Getl_circumflex:getll2: meta apo enhmerwsh h initialize 
        {
            if (endeixiBL11 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BL11 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BL11[j] = new double[6, 9]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL11[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BL11[j][k, l] += GetBL11a(element)[j][k, m] * BL11b[j][m, l];
                            }
                        }
                    }
                }
                endeixiBL11 = 2;
                //return BL11; 
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL11[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BL11[j][k, l] += GetBL11a(element)[j][k, m] * BL11b[j][m, l];
                            }
                        }
                    }
                }
                //return BL11;
            }
        }


        private double[][,] BL13;
        public int endeixilBL13 = 1; //analogws endeixi
        private void CalculateBL13() //afou periexei tVn_i kai Getll2: Xrhsimopoieitai meta apo 2)ENHMERWSH h 1)INITIALIZE 
        {
            if (endeixilBL13 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BL13 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { BL13[j] = new double[9, 40]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            BL13[j][k, l] = 0;
                        }
                    }

                    //sthles 1:3
                    for (int m = 0; m < 8; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 2; l++)
                            {
                                BL13[j][3 * k + l, 5 * m + k] = GetShapeFunctionDerivatives()[m + 8 * l][j];
                            }
                        }
                    }

                    //sthles 4:5
                    for (int m = 0; m < 8; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                BL13[j][3 * k + l, 5 * m + 3] = -GetJ_0a()[j][l, m * 2 + 1] * tUvec[m][3 + k];
                                BL13[j][3 * k + l, 5 * m + 4] = +GetJ_0a()[j][l, m * 2 + 1] * tUvec[m][k];
                            }
                        }
                    }
                }
                endeixilBL13 = 2;
                //return BL13;
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    //sthles 4:5
                    for (int m = 0; m < 8; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                BL13[j][3 * k + l, 5 * m + 3] = -GetJ_0a()[j][l, m * 2 + 1] * tUvec[m][3 + k];
                                BL13[j][3 * k + l, 5 * m + 4] = +GetJ_0a()[j][l, m * 2 + 1] * tUvec[m][k];
                            }
                        }
                    }
                }
                //return BL13;
            }
        }

        private double[,] J_1b;    //einai idio gia ola ta gauss points
        public int endeixiJ_1b = 1;
        private void CalculateJ_1b(Element element) // meta apo enhmerwsi i initialize twn tx_i,tVn_i
        {
            if (endeixiJ_1b == 1)
            {
                J_1b = new double[16, 3];
                for (int j = 0; j < 8; j++)
                {
                    J_1b[2 * j, 0] = tx_i[j][0];
                    J_1b[2 * j + 1, 0] = tU[j][3];
                    J_1b[2 * j, 1] = tx_i[j][1];
                    J_1b[2 * j + 1, 1] = tU[j][4];
                    J_1b[2 * j, 2] = tx_i[j][2];
                    J_1b[2 * j + 1, 2] = tU[j][5];
                }
                endeixiJ_1b = 2;

            }
            else
            {
                for (int j = 0; j < 8; j++)
                {
                    J_1b[2 * j, 0] = tx_i[j][0];
                    J_1b[2 * j + 1, 0] = tU[j][3];
                    J_1b[2 * j, 1] = tx_i[j][1];
                    J_1b[2 * j + 1, 1] = tU[j][4];
                    J_1b[2 * j, 2] = tx_i[j][2];
                    J_1b[2 * j + 1, 2] = tU[j][5];
                }
            }
        }

        private double[][,] J_1;       //den einai to idio gia ola ta gausspoint
        public int endeixiJ_1 = 1;
        private void CalculateJ_1(Element element)   // meta apo enhmerwsi i initialize twn tx_i,tVn_i
        {                                              //Meta Apo CALCULATE J_1b 
            if (endeixiJ_1 == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                J_1 = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { J_1[j] = new double[3, 3]; }
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            J_1[j][k, l] = 0;
                            for (int m = 0; m < 16; m++)
                            {
                                J_1[j][k, l] += GetJ_0a()[j][k, m] * J_1b[m, l];
                            }

                        }

                    }
                }
                endeixiJ_1 = 2;

            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            J_1[j][k, l] = 0;
                            for (int m = 0; m < 16; m++)
                            {
                                J_1[j][k, l] += GetJ_0a()[j][k, m] * J_1b[m, l];
                            }

                        }

                    }
                }
            }

        }

        private double[][,] DefGradTr;       //den einai to idio gia ola ta gausspoint
        public int endeixiDefGradTr = 1;
        private void CalculateDefGradTr(Element element) // Meta apo CalculateJ_1 profanws
        {
            if (endeixiDefGradTr == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                DefGradTr = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                { DefGradTr[j] = new double[3, 3]; }
                //gemisma pol/smos
                for (int j = 0; j < nGaussPoints; j++)
                {

                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            DefGradTr[j][k, l] = 0;
                            for (int m = 0; m < 3; m++)
                            {
                                DefGradTr[j][k, l] += GetJ_0inv(element)[j][k, m] * J_1[j][m, l];
                            }

                        }

                    }

                }

            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {

                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            DefGradTr[j][k, l] = 0;
                            for (int m = 0; m < 3; m++)
                            {
                                DefGradTr[j][k, l] += GetJ_0inv(element)[j][k, m] * J_1[j][m, l];
                            }

                        }

                    }
                }

            }
        }

        private double[][,] GL;       //den einai to idio gia ola ta gausspoint
        public int endeixiGL = 1;
        private void CalculateGL() // Meta apo CalculateDefGradTr profanws
        {
            if (endeixiGL == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                GL = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    GL[j] = new double[3, 3];
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            GL[j][k, l] = 0;
                            for (int m = 0; m < 3; m++)
                            {
                                GL[j][k, l] += DefGradTr[j][k, m] * DefGradTr[j][l, m];
                            }

                        }

                    }
                    for (int k = 0; k < 3; k++)
                    {
                        GL[j][k, k] = GL[j][k, k] - 1;
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        { GL[j][k, l] = 0.5 * GL[j][k, l]; }
                    }

                }
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            GL[j][k, l] = 0;
                            for (int m = 0; m < 3; m++)
                            {
                                GL[j][k, l] += DefGradTr[j][k, m] * DefGradTr[j][l, m];
                            }

                        }

                    }
                    for (int k = 0; k < 3; k++)
                    {
                        GL[j][k, k] = GL[j][k, k] - 1;
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        { GL[j][k, l] = 0.5 * GL[j][k, l]; }
                    }
                }
            }
        }

        private double[][] GLvec;
        public int endeixiGLvec = 1;
        private void CalculateGLvec()//meta apo calculate Gl
        {
            if (endeixiGLvec == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                GLvec = new double[nGaussPoints][];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    GLvec[j] = new double[6];
                    for (int k = 0; k < 3; k++)
                    { GLvec[j][k] = GL[j][k, k]; }
                    GLvec[j][3] = 2 * GL[j][0, 1];
                    GLvec[j][4] = 2 * GL[j][1, 2];
                    GLvec[j][5] = 2 * GL[j][2, 0];
                }
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 3; k++)
                    { GLvec[j][k] = GL[j][k, k]; }
                    GLvec[j][3] = 2 * GL[j][0, 1];
                    GLvec[j][4] = 2 * GL[j][1, 2];
                    GLvec[j][5] = 2 * GL[j][2, 0];
                }
            }
        }


        private double[] E;
        private double[] ni;
        private double[,] Cons;
        private double[] V3;
        private double V3_norm;
        private double[] V1;
        private double V1_norm;
        private double[] V2;
        private double[,] T_e;
        private double l1;
        private double m1;
        private double n1;
        private double l2;
        private double m2;
        private double n2;
        private double l3;
        private double m3;
        private double n3;
        private double[,] Cons_T_e;
        private double[][,] ConsCartes;
        private void CalculateCons()
        {
            V3 = new double[3];
            V1 = new double[3];
            V2 = new double[3];
            T_e = new double[6, 6];
            nGaussPoints = gp_d1 * gp_d2 * gp_d3;
            ConsCartes = new double[nGaussPoints][,];
            E = new double[nGaussPoints];
            ni = new double[nGaussPoints];
            Cons = new double[6, 6];
            Cons_T_e = new double[6, 6];
            for (int j = 0; j < nGaussPoints; j++)
            {
                E[j] = materialsAtGaussPoints[j].YoungModulus;
                ni[j] = materialsAtGaussPoints[j].PoissonRatio;
                ConsCartes[j] = new double[6, 6];
                for (int k = 0; k < 2; k++)
                { Cons[k, k] = E[j] / (1 - Math.Pow(ni[j], 2)); }
                Cons[0, 1] = ni[j] * E[j] / (1 - Math.Pow(ni[j], 2));
                Cons[1, 0] = ni[j] * E[j] / (1 - Math.Pow(ni[j], 2));
                Cons[3, 3] = (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2));
                Cons[4, 4] = (1 - ni[j]) * (0.41666666667) * E[j] / (1 - Math.Pow(ni[j], 2));
                Cons[5, 5] = (1 - ni[j]) * (0.41666666667) * E[j] / (1 - Math.Pow(ni[j], 2));
                //for (int k = 0; k < 2; k++)
                //{ Cons[4 + k, 4 + k] = (5 / 6) * (1 - ni[j]) * (0.5) * E[j] / (1 - Math.Pow(ni[j], 2)); }

                for (int k = 0; k < 3; k++)
                { V3[k] = 0; V1[k] = 0; V2[k] = 0; }

                for (int k = 0; k < 8; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        V3[l] += GetShapeFunctions()[k][j] * oVn_i[k][l];
                        V1[l] += GetShapeFunctions()[k][j] * oV1_i[k][l];
                    }
                }
                V3_norm = Math.Sqrt(V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2]);
                V1_norm = Math.Sqrt(V1[0] * V1[0] + V1[1] * V1[1] + V1[2] * V1[2]);
                for (int l = 0; l < 3; l++)
                {
                    V3[l] = V3[l] / V3_norm;
                    V1[l] = V1[l] / V1_norm;
                }

                V2[0] = V3[1] * V1[2] - V3[2] * V1[1];
                V2[1] = V3[2] * V1[0] - V3[0] * V1[2];
                V2[2] = V3[0] * V1[1] - V3[1] * V1[0];

                l1 = V1[0];
                m1 = V1[1];
                n1 = V1[2];

                l2 = V2[0];
                m2 = V2[1];
                n2 = V2[2];

                l3 = V3[0];
                m3 = V3[1];
                n3 = V3[2];

                for (int i = 0; i < 3; i++)
                {
                    T_e[0, i] = (V1[i] * V1[i]);
                    T_e[1, i] = (V2[i] * V2[i]);
                    T_e[2, i] = (V3[i] * V3[i]);

                    T_e[3, i] = (2 * V1[i] * V2[i]);
                    T_e[4, i] = (2 * V2[i] * V3[i]);
                    T_e[5, i] = (2 * V3[i] * V1[i]);

                    T_e[0, 3 + i] = (V1[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                    T_e[1, 3 + i] = (V2[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                    T_e[2, 3 + i] = (V3[i] * V3[1 + i - 3 * i * (i - 1) / 2]);

                    T_e[3, 3 + i] = (V1[i] * V2[1 + i - 3 * i * (i - 1) / 2] + V2[i] * V1[1 + i - 3 * i * (i - 1) / 2]);
                    T_e[4, 3 + i] = (V2[i] * V3[1 + i - 3 * i * (i - 1) / 2] + V3[i] * V2[1 + i - 3 * i * (i - 1) / 2]);
                    T_e[5, 3 + i] = (V3[i] * V1[1 + i - 3 * i * (i - 1) / 2] + V1[i] * V3[1 + i - 3 * i * (i - 1) / 2]);
                }

                // multiplication [Te']*[cons]*[Te];

                for (int i = 0; i < 6; i++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        Cons_T_e[i, k] = 0;
                        for (int l = 0; l < 6; l++)
                        { Cons_T_e[i, k] += Cons[i, l] * T_e[l, k]; }
                    }
                }

                for (int i = 0; i < 6; i++)
                {
                    for (int k = 0; k < 6; k++)
                    {
                        ConsCartes[j][i, k] = 0;
                        for (int l = 0; l < 6; l++)
                        { ConsCartes[j][i, k] += T_e[l, i] * Cons_T_e[l, k]; }
                    }
                }
            }

        }
       
    

        private double[][] SPKvec;
        private double[][,] SPK_circumflex;
        private int endeixiSPK = 1;
        private void CalculateSPK()
        {
            if (endeixiSPK == 1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                SPKvec = new double[nGaussPoints][];
                SPK_circumflex = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    SPKvec[j] = new double[6];
                    SPK_circumflex[j] = new double[9, 9];
                    for (int l = 0; l < 6; l++)
                        {
                            SPKvec[j][l] = 0;
                            for (int m = 0; m < 6; m++)
                            {
                            SPKvec[j][l] += ConsCartes[j][l, m] * GLvec[j][m];
                            }

                        }
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        { SPK_circumflex[j][k, l] = 0; }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            SPK_circumflex[j][3 * k + l, 3 * k + l] = SPKvec[j][l];
                        }
                        SPK_circumflex[j][3 * k , 3 * k +1] = SPKvec[j][3];
                        SPK_circumflex[j][3 * k , 3 * k+2] = SPKvec[j][5];
                        SPK_circumflex[j][3 * k + 1, 3 * k+2] = SPKvec[j][4];

                        SPK_circumflex[j][3 * k+1, 3 * k ] = SPKvec[j][3];
                        SPK_circumflex[j][3 * k+2, 3 * k] = SPKvec[j][5];
                        SPK_circumflex[j][3 * k + 2, 3 * k + 1] = SPKvec[j][4];                        
                    }

                }
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {                   
                    for (int l = 0; l < 6; l++)
                    {
                        SPKvec[j][l] = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            SPKvec[j][l] += ConsCartes[j][l, m] * GLvec[j][m];
                        }

                    }                    
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            SPK_circumflex[j][3 * k + l, 3 * k + l] = SPKvec[j][l];
                        }
                        SPK_circumflex[j][3 * k, 3 * k + 1] = SPKvec[j][3];
                        SPK_circumflex[j][3 * k, 3 * k + 2] = SPKvec[j][5];
                        SPK_circumflex[j][3 * k + 1, 3 * k + 2] = SPKvec[j][4];

                        SPK_circumflex[j][3 * k + 1, 3 * k] = SPKvec[j][3];
                        SPK_circumflex[j][3 * k + 2, 3 * k] = SPKvec[j][5];
                        SPK_circumflex[j][3 * k + 2, 3 * k + 1] = SPKvec[j][4];
                    }

                }
            }
        }

        private double[][,] ck;// 1 ana komvo kai ana gauss Point dld [GP][8komvoi,diastash9]
        private int endeixiCk = 1;
        private void CalculateCk()
        {
            if (endeixiCk == 1)
            {
                //initialize
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                ck = new double[nGaussPoints][,];
                for (int j = 0; j < nGaussPoints; j++)
                {
                    ck[j] = new double[8,9];
                }
                //tupoi
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int m = 0; m < 8; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                ck[j][m, 3 * k + l] = GetJ_0a()[j][l, 2 * m + 1] * tU[m][3 + k];
                            }
                        }
                    }
                }
            }
            else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int m = 0; m < 8; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                ck[j][m, 3 * k + l] = GetJ_0a()[j][l, 2 * m + 1] * tU[m][3 + k];
                            }
                        }
                    }
                }
            }
        }


        private double[][] kck;// 1 ana komvo kai (ana Gauss point+1 gia to athroiskma) [GP][8 vathmoi komvoi]
                               // to initialize tou einai comment out parapanw

        private double[][,] BNL;
        private double[][,] KNL;

        private double[][,] KL;
        private double[][,] BL1_2;
        //private double[][,] BL1;
        //private double[][,] BL0;
        private double[][,] BL;

        private double[][,] ConsBL;
        private double[][,] S_BNL;

        private double[][,] BL01plus1_2;
        private double[][] BL01plus1_2tSPKvec;

        private double[,] Kt = new double[40, 40];

        private double[][] Fxk;

        private int endeixiKmatrices = 1;
        private void CalculateKmatrices(Element element)
        {
        if (endeixiKmatrices==1)
            {
                nGaussPoints = gp_d1 * gp_d2 * gp_d3;
                BNL = new double[nGaussPoints][,];
                BL1_2 = new double[nGaussPoints][,];
                //BL1 = new double[nGaussPoints][,];
                //BL0 = new double[nGaussPoints][,];
                BL = new double[nGaussPoints][,];

                kck = new double[nGaussPoints+1][];
                KL = new double[nGaussPoints + 1][,];
                KNL = new double[nGaussPoints + 1][,];

                ConsBL = new double[nGaussPoints][,];
                S_BNL = new double[nGaussPoints][,];

                kck = new double[nGaussPoints + 1][];

                BL01plus1_2 = new double[nGaussPoints][,];
                BL01plus1_2tSPKvec = new double[nGaussPoints][];

                //Kt = new double[40, 40];

                Fxk = new double[nGaussPoints + 1][];

                for (int j = 0; j < nGaussPoints; j++)
                {
                    BNL[j] = new double[9, 40];
                    BL1_2[j] = new double[6, 9];
                    //BL1[j] = new double[6, 40];
                    //BL0[j] =new double[6, 40];
                    BL[j]= new double[6, 40];

                    ConsBL[j] = new double[6,40];
                    S_BNL[j] = new double[9, 40];

                    kck[j] = new double[8];

                    BL01plus1_2[j] = new double[6, 9];
                    BL01plus1_2tSPKvec[j] = new double[9];
                }

                for (int j = 0; j < nGaussPoints+1; j++)
                {
                    kck[j] = new double[8];
                    KL[j] = new double[40, 40];
                    KNL[j] = new double[40, 40];

                    Fxk[j]= new double[40];
                }

                //prepei na ginei gemisma twn parapanw mhtrwwn pollaplasiasmoi
                //athroisma olwn twn gausspoints kai prosthesi K Knl kai kck ana orous

                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            BNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BNL[j][k, l] += GetBNL1(element)[j][k, m] * BL13[j][m, l]; 
                            }

                        }

                    }

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL1_2[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BL1_2[j][k, l] += BL11[j][k, m] * GetBL12(element)[j][m, l]; //TODO BL11 keno kai BL12 null thelei getbl12 edw kai parakatw
                            }                                                   //vriskomaste sto calculate Kmatrices eprepe na trexei to calculate BL11 prwta

                        }

                    }

                    //for (int k = 0; k < 6; k++)
                    //{
                    //    for (int l = 0; l < 40; l++)
                    //    {
                    //        BL1[j][k, l] = 0;
                    //        BL0[j][k, l] = 0;
                    //        for (int m = 0; m < 9; m++) //panw apo to for BLx=BL1_2+BL11 kai mesa sto for BL=BLx*BL13
                    //        {
                    //            BL1[j][k, l] += BL1_2[j][k, m] * BL13[j][m, l];
                    //            BL0[j][k, l] += BL11[j][k, m]* BL13[j][m, l];
                    //        }
                    //        BL[j][k, l] = BL0[j][k, l] + BL1[j][k, l];
                    //    }
                    //}

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL01plus1_2[j][k, l] = BL1_2[j][k, l] + GetBL11a(element)[j][k, l];
                        }
                    }

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            BL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++) 
                            {
                                BL[j][k, l] += BL01plus1_2[j][k,m] * BL13[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 9; k++)
                    {                    
                            BL01plus1_2tSPKvec[j][k] = 0;
                            for (int m = 0; m < 6; m++)
                            {
                                BL01plus1_2tSPKvec[j][k] += BL01plus1_2[j][m, k] * SPKvec[j][m];
                            }
                    }

                    for (int k = 0; k < 8; k++)
                    {
                        kck[j][k] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            kck[j][k] += ck[j][k,m] * BL01plus1_2tSPKvec[j][m];
                        }
                    }

                    // porsthetoume kai to kck ws extra(den prokuptei apo ta comment out

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            ConsBL[j][k, l] = 0;
                            for (int m = 0; m < 6; m++) 
                            {
                                ConsBL[j][k, l] += ConsCartes[j][k, m] * BL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            S_BNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++) 
                            {
                                S_BNL[j][k, l] += SPK_circumflex[j][k, m] * BNL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            KNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++) 
                            {
                                KNL[j][k, l] += BNL[j][m, k] * S_BNL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            KL[j][k, l] = 0;
                            for (int m = 0; m < 6; m++) 
                            {
                                KL[j][k, l] += BL[j][m, k] * ConsBL[j][m, l];
                            }
                        }
                    }
                }

                //PrintUtilities.WriteToFile(KL[0], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\KL_gp1_1.txt");
                //PrintUtilities.WriteToFile(KNL[0], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\KNL_gp1_1.txt");
                //PrintUtilities.WriteToFile(KL[1], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\KL_gp2.txt");
                //PrintUtilities.WriteToFile(KNL[1], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\KNL_gp2.txt");
                //PrintUtilities.WriteToFile(BL[0], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\BL_gp1_1.txt");
                //PrintUtilities.WriteToFile(BNL[0], @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\BNL_gp1_1.txt");



                // morfwsi telikou mhtrwou
                for (int k = 0; k < 40; k++)
                {
                    for (int l = 0; l < 40; l++)
                    { Kt[k, l] = 0; }
                }
                for (int j = 0; j < nGaussPoints; j++)
                {

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        { Kt[k, l] += a_123g[j] * GetDetJ_0(element)[j]*(KL[j][k,l]+KNL[j][k,l]); }
                    }
                                        
                    for (int l = 0; l < 8; l++)
                    {
                        Kt[5 * l + 3, 5 * l + 3] += a_123g[j] * GetDetJ_0(element)[j]*kck[j][l];
                        Kt[5 * l + 4, 5 * l + 4] += a_123g[j] * GetDetJ_0(element)[j]*kck[j][l];
                    }                                                        
                }

                //mprfwsi drasewn   
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 40; k++)
                    {
                        Fxk[j][k] = 0;
                        for (int m = 0; m < 6; m++)
                        { Fxk[j][k] += BL[j][m, k] * SPKvec[j][m]; }
                    }
                }
                for (int k = 0; k < 40; k++)
                {
                    Fxk[nGaussPoints][k] = 0;
                    for (int j = 0; j < nGaussPoints; j++)
                    { Fxk[nGaussPoints][k] += a_123g[j] * GetDetJ_0(element)[j]*Fxk[j][k]; }
                }

                endeixiKmatrices = 2;
            }

        else
            {
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            BNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BNL[j][k, l] += GetBNL1(element)[j][k, m] * BL13[j][m, l];
                            }

                        }

                    }

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL1_2[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BL1_2[j][k, l] += BL11[j][k, m] * GetBL12(element)[j][m, l];
                            }

                        }

                    }

                    //for (int k = 0; k < 6; k++)
                    //{
                    //    for (int l = 0; l < 40; l++)
                    //    {
                    //        BL1[j][k, l] = 0;
                    //        BL0[j][k, l] = 0;
                    //        for (int m = 0; m < 9; m++) //panw apo to for BLx=BL1_2+BL11 kai mesa sto for BL=BLx*BL13
                    //        {
                    //            BL1[j][k, l] += BL1_2[j][k, m] * BL13[j][m, l];
                    //            BL0[j][k, l] += BL11[j][k, m]* BL13[j][m, l];
                    //        }
                    //        BL[j][k, l] = BL0[j][k, l] + BL1[j][k, l];
                    //    }
                    //}

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 9; l++)
                        {
                            BL01plus1_2[j][k, l] = BL1_2[j][k, l] + GetBL11a(element)[j][k, l];
                        }
                    }

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            BL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                BL[j][k, l] += BL01plus1_2[j][k, m] * BL13[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 9; k++)
                    {
                        BL01plus1_2tSPKvec[j][k] = 0;
                        for (int m = 0; m < 6; m++)
                        {
                            BL01plus1_2tSPKvec[j][k] += BL01plus1_2[j][m, k] * SPKvec[j][m];
                        }
                    }

                    for (int k = 0; k < 8; k++)
                    {
                        kck[j][k] = 0;
                        for (int m = 0; m < 9; m++)
                        {
                            kck[j][k] += ck[j][k, m] * BL01plus1_2tSPKvec[j][m];
                        }
                    }

                    // porsthetoume kai to kck ws extra(den prokuptei apo ta comment out

                    for (int k = 0; k < 6; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            ConsBL[j][k, l] = 0;
                            for (int m = 0; m < 6; m++)
                            {
                                ConsBL[j][k, l] += ConsCartes[j][k, m] * BL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 9; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            S_BNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                S_BNL[j][k, l] += SPK_circumflex[j][k, m] * BNL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            KNL[j][k, l] = 0;
                            for (int m = 0; m < 9; m++)
                            {
                                KNL[j][k, l] += BNL[j][m, k] * S_BNL[j][m, l];
                            }
                        }
                    }

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        {
                            KL[j][k, l] = 0;
                            for (int m = 0; m < 6; m++)
                            {
                                KL[j][k, l] += BL[j][m, k] * ConsBL[j][m, l];
                            }
                        }
                    }
                }

                // morfwsi telikou mhtrwou
                for (int k = 0; k < 40; k++)
                {
                    for (int l = 0; l < 40; l++)
                    { Kt[k, l] = 0; }
                }
                for (int j = 0; j < nGaussPoints; j++)
                {

                    for (int k = 0; k < 40; k++)
                    {
                        for (int l = 0; l < 40; l++)
                        { Kt[k, l] += a_123g[j]*GetDetJ_0(element)[j]*(KL[j][k, l] + KNL[j][k, l]); }
                    }

                    for (int l = 0; l < 8; l++)
                    {
                        Kt[5 * l + 3, 5 * l + 3] += a_123g[j] * GetDetJ_0(element)[j] * kck[j][l];
                        Kt[5 * l + 4, 5 * l + 4] += a_123g[j] * GetDetJ_0(element)[j] * kck[j][l];
                    }
                }

                //mprfwsi drasewn   
                for (int j = 0; j < nGaussPoints; j++)
                {
                    for (int k = 0; k < 40; k++)
                    {
                        Fxk[j][k] = 0;
                        for (int m = 0; m < 6; m++)
                        { Fxk[j][k] += BL[j][m, k] * SPKvec[j][m]; }
                    }
                }
                for (int k = 0; k < 40; k++)
                {
                    Fxk[nGaussPoints][k] = 0;
                    for (int j = 0; j < nGaussPoints; j++)
                    { Fxk[nGaussPoints][k] += a_123g[j] * GetDetJ_0(element)[j]*Fxk[j][k]; }
                }
            }
        }

        // ANANEWSH thw thesis tou stoixeiou-----------------------------------------
        // voithitikes metavlhtes gia upologismo strofhs-----------------------------
        private double[] ak_total = new double[8];
        private double[] bk_total = new double[8];

        // metavlhtes gia anafora stis strofes kai voithitikoi pinakes
        private double ak;
        private double bk;
        private double gk1;
        private double[,] Q = new double[3,3];
        private double[,] Q2 = new double[3, 3];
        //private double[] tdtVn = new double[3];
        private double tV1norm;

        private void UpdateCoordinateData(double[] localdisplacements)
        {
            for (int k = 0; k < 8; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    tx_i[k][l] = ox_i[k][l] + localdisplacements[5 * k + l];
                }
                //update twn tU kai tUvec 
                tU[k][0] = localdisplacements[5 * k + 0];
                tU[k][1] = localdisplacements[5 * k + 1];
                tU[k][2] = localdisplacements[5 * k + 2];
                ak = localdisplacements[5 * k + 3] - ak_total[k];
                ak_total[k] = localdisplacements[5 * k + 3];
                bk = localdisplacements[5 * k + 4] - bk_total[k];
                bk_total[k] = localdisplacements[5 * k + 4];
                this.RotateNodalDirectionVectors(ak, bk, k);
                // update twn tU kai tUvec ews edw         
            }
        }
        //prepei sto telos tou upologismou drasewn na enhmerwnontai oi ak_total kai bk_total

        // voithitikes metavlhtes gia thn peristrofh
        private double[] tdtVn = new double[3];
        private double[] tdtV1 = new double[3];
        private double[] tdtV2 = new double[3];
        private double theta;
        private double[] theta_vec = new double[3];
        private double[,] s_k = new double[3, 3];
        private void RotateNodalDirectionVectors(double ak, double bk, int n_vector)
        {
            for (int j = 0; j < 3; j++)
            {
                theta_vec[j] = ak * tUvec[n_vector][j] + bk * tUvec[n_vector][3 + j];
            }
            theta = Math.Sqrt((theta_vec[0] * theta_vec[0]) + (theta_vec[1] * theta_vec[1]) + (theta_vec[2] * theta_vec[2]));
            if (theta > 0)
            {
                s_k[0, 1] = -theta_vec[2];
                s_k[0, 2] = theta_vec[1];
                s_k[1, 0] = theta_vec[2];
                s_k[1, 2] = -theta_vec[0];
                s_k[2, 0] = -theta_vec[1];
                s_k[2, 1] = theta_vec[0];

                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] = (Math.Sin(theta) / theta) * s_k[j, m];
                    }
                }

                for (int m = 0; m < 3; m++)
                {
                    Q[m, m] += 1;
                }
                gk1 = 0.5 * ((Math.Sin(0.5 * theta) / (0.5 * theta)) * (Math.Sin(0.5 * theta) / (0.5 * theta)));
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q2[j, m] = 0;
                        for (int n = 0; n < 3; n++)
                        { Q2[j, m] += gk1 * s_k[j, n] * s_k[n, m]; }
                    }
                }
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] += Q2[j, m];
                    }
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtVn[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtVn[j] += Q[j, m] * tU[n_vector][3 + m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tU[n_vector][3 + j] = tdtVn[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV1[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV1[j] += Q[j, m] * tUvec[n_vector][m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tUvec[n_vector][j] = tdtV1[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV2[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV2[j] += Q[j, m] * tUvec[n_vector][3 + m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tUvec[n_vector][3 + j] = tdtV2[j];
                }
            }
        }


        // aparaithta tou IStructuralFiniteElement

        public int ID
        {
            get { return 12; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofTypes;
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        //aparaithta tou IstructuralElement gia to material
        public void ClearMaterialState()
        {
            foreach (IFiniteElementMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (IFiniteElementMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IFiniteElementMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (IFiniteElementMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (IFiniteElementMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return new Tuple<double[], double[]>(new double[123], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
        }

        //aparaithta tou Istructural gia th dunamiki analusi
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[123];
        }

        public virtual IMatrix2D<double> MassMatrix(Element element)
        {
            return new Matrix2D<double>(1,1);
        }

        public virtual IMatrix2D<double> DampingMatrix(Element element)
        {
            
            return new Matrix2D<double>(1, 1);            
        }

        // forces tha exei ena if me initial geometric data
        // kai to else me strofi kai meta olous tous upologismous.

        // mporei kai me ena if ean einai to trito dianusma mhdeniko. na mhn ektelountai kapoioi upologismoi.

        //implementation of basic methods

        private int endeixiArxikDianusmatwn = 1;
        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            if (endeixiArxikDianusmatwn==1)
            {                
                //this.GetInitialGeometricData(element);
                //CalculateCons();

                this.UpdateCoordinateData(localTotalDisplacements);
                Calculatell2();
                Calculatel_circumflex();
                CalculateBL11b();
                CalculateBL11(element);
                CalculateBL13();
                CalculateJ_1b(element);
                CalculateJ_1(element);
                CalculateDefGradTr(element);
                CalculateGL();
                CalculateGLvec();
                CalculateSPK();
                CalculateCk();
                CalculateKmatrices(element);

                endeixiArxikDianusmatwn = 2;
                return Fxk[nGaussPoints];
            }
            else
            {
                this.UpdateCoordinateData(localTotalDisplacements);// mporei dokimastika na ginei sxoleio meta
                Calculatell2();
                Calculatel_circumflex();
                CalculateBL11b();
                CalculateBL11(element);
                CalculateBL13();
                CalculateJ_1b(element);
                CalculateJ_1(element);
                CalculateDefGradTr(element);
                CalculateGL();
                CalculateGLvec();
                CalculateSPK();
                CalculateCk();
                CalculateKmatrices(element);

                return Fxk[nGaussPoints];
            }
        }

        private int endeixiStiffness = 1;
        public virtual IMatrix2D<double> StiffnessMatrix(Element element)
        {
            if (endeixiStiffness == 1)
            {
                this.GetInitialGeometricData(element);
                this.CalculateCons();
                
                //this.UpdateCoordinateData(localTotalDisplacements);
                Calculatell2();
                Calculatel_circumflex();
                CalculateBL11b();
                CalculateBL11(element);
                CalculateBL13();
                CalculateJ_1b(element);
                CalculateJ_1(element);
                CalculateDefGradTr(element);
                CalculateGL();
                CalculateGLvec();
                CalculateSPK();
                CalculateCk();
                CalculateKmatrices(element);
                endeixiStiffness = 2;
                PrintUtilities.WriteToFile(Kt, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\Copy_of_apo_ta_shell_New_Load_Case\Kt_1.txt");
                IMatrix2D<double> iGlobalStiffnessMatrix = new Matrix2D<double>(Kt);
                return dofEnumerator.GetTransformedMatrix(iGlobalStiffnessMatrix);
            }
            else
            {
                IMatrix2D<double> iGlobalStiffnessMatrix = new Matrix2D<double>(Kt);
                return dofEnumerator.GetTransformedMatrix(iGlobalStiffnessMatrix);
            }
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }
    }
}

