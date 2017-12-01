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
    class Hexa8NL //: IStructuralFiniteElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IFiniteElementMaterial3D[] materialsAtGaussPoints;
        protected IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw

        public int gp_d1_disp { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public int gp_d2_disp { get; set; }
        public int gp_d3_disp { get; set; }
        private int nGaussPoints;

        protected Hexa8NL()//consztructor apo to hexa8
        {
        }

        public Hexa8NL(IFiniteElementMaterial3D material, int gp_d1c, int gp_d2c, int gp_d3c)
        {
            this.gp_d1_disp = gp_d1c;
            this.gp_d2_disp = gp_d2c;
            this.gp_d3_disp = gp_d3c;
            this.nGaussPoints = this.gp_d1_disp * this.gp_d2_disp*this.gp_d3_disp;
            materialsAtGaussPoints = new IFiniteElementMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IFiniteElementMaterial3D)material.Clone();

        }

        public int endeixiShapeFunctionAndGaussPointData = 1;
        private double[] a_123g;
        private double a_1g;
        private double a_2g;
        private double a_3g;
        private double ksi;
        private double heta;
        private double zeta;
        private int npoint;

        private double[,] Ni;
        private double[,] Ni_ksi;
        private double[,] Ni_heta;
        private double[,] Ni_zeta;

        private double [] [,] ll1_hexa;
        private double[][,] BL13_hexa;

        private void CalculateShapeFunctionAndGaussPointData()
        {
            nGaussPoints = gp_d1_disp * gp_d2_disp * gp_d3_disp;
            a_123g = new double[nGaussPoints];

            Ni = new double[8,nGaussPoints]; // den sxetizetai me ta coh elements alla
            Ni_ksi = new double[8, nGaussPoints]; // me to prokat_disp (einai shapefunctionData)
            Ni_heta = new double[8, nGaussPoints]; // 
            Ni_heta = new double[8, nGaussPoints];

            ll1_hexa= new double[nGaussPoints][,];
            BL13_hexa = new double[nGaussPoints][,];

            for (int l = 0; l < gp_d3_disp; l++)
            {
                for (int k = 0; k < gp_d2_disp; k++)
                {
                    for (int j = 0; j < gp_d1_disp; j++)
                    {
                        npoint = l * (gp_d1_disp * gp_d2_disp) + k * gp_d1_disp + j;
                        if (gp_d1_disp == 3)
                        {
                            ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                            a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                        }
                        if (gp_d1_disp == 2)
                        {
                            ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                            a_1g = 1;
                        }
                        if (gp_d2_disp == 3)
                        {
                            heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                            a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                        }
                        if (gp_d2_disp == 2)
                        {
                            heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                            a_2g = 1;
                        }
                        if (gp_d3_disp == 3)
                        {
                            zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483);
                            a_3g = 0.5 * (l - 1) * (l - 2) * (0.555555555555555) + (-1) * (l) * (l - 2) * (0.888888888888888) + 0.5 * (l) * (l - 1) * (0.555555555555555);
                        }
                        if (gp_d3_disp == 2)
                        {
                            zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1);
                            a_3g = 1;
                        }                        
                        a_123g[npoint] = a_1g * a_2g * a_3g;

                        Ni[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 + zeta); //N1(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 + zeta); //N2(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 + zeta);
                        Ni[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 + zeta);
                        Ni[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 - zeta); 
                        Ni[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 - zeta);
                        Ni[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 - zeta);
                        Ni[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 - zeta);

                        Ni_ksi[0, npoint] = +0.125 * (1 + heta) * (1 + zeta); //N1_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[1, npoint] = -0.125 * (1 + heta) * (1 + zeta); //N2_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[2, npoint] = -0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[3, npoint] = +0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[4, npoint] = +0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[5, npoint] = -0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[6, npoint] = -0.125 * (1 - heta) * (1 - zeta);
                        Ni_ksi[7, npoint] = +0.125 * (1 - heta) * (1 - zeta);

                        Ni_heta[0, npoint] = 0.125 * (1 + ksi) * (+1) * (1 + zeta); //N1_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[1, npoint] = 0.125 * (1 - ksi) * (+1) * (1 + zeta); //N2_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[2, npoint] = 0.125 * (1 - ksi) * (-1) * (1 + zeta);
                        Ni_heta[3, npoint] = 0.125 * (1 + ksi) * (-1) * (1 + zeta);
                        Ni_heta[4, npoint] = 0.125 * (1 + ksi) * (+1) * (1 - zeta);
                        Ni_heta[5, npoint] = 0.125 * (1 - ksi) * (+1) * (1 - zeta);
                        Ni_heta[6, npoint] = 0.125 * (1 - ksi) * (-1) * (1 - zeta);
                        Ni_heta[7, npoint] = 0.125 * (1 + ksi) * (-1) * (1 - zeta);

                        Ni_zeta[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (+1); //N1_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (+1); //N2_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (+1);
                        Ni_zeta[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (+1);
                        Ni_zeta[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (-1);
                        Ni_zeta[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (-1);
                        Ni_zeta[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (-1);
                        Ni_zeta[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (-1);

                        ll1_hexa[npoint] = new double[3, 8];
                        for (int m = 0; m < 8; m++)
                        {
                            ll1_hexa[npoint][0, m] = Ni_ksi[m, npoint];
                            ll1_hexa[npoint][1, m] = Ni_heta[m, npoint];
                            ll1_hexa[npoint][2, m] = Ni_zeta[m, npoint];
                        }

                        BL13_hexa[npoint] = new double[9, 24];
                        for (int m = 0; m < 9; m++)
                        {
                            for (int n = 0; n < 24; n++)
                            { BL13_hexa[npoint][m, n] = 0; }
                        }

                        for (int m = 0; m < 8; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                BL13_hexa[npoint][n, 3 * m + 0] = ll1_hexa[npoint][n, m];
                                BL13_hexa[npoint][ n+3,3*m+1] = ll1_hexa[npoint][n, m];
                                BL13_hexa[npoint][ n+6,3*m+2] = ll1_hexa[npoint][n, m];
                            }
                        }


                    }
                }
            }

            endeixiShapeFunctionAndGaussPointData = 2;
        }

        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths 8 arrays twn 3 stoixeiwn
        private double[][,] J_0b_hexa; // exoume tosa [,] osa einai kai ta gpoints
        private double[][,] J_0_hexa;
        private double[][,] J_0inv_hexa;
        private double[] detJ_0; //osa kai ta gpoints
        private double[] sunt_oloklhrwmatos; // omoiws
        private double[][,] BL11a_hexa; // exoume tosa [,] osa einai kai ta gpoints
        private double[][,] BL12_hexa;
        private double[][,] BL01_hexa;
        private double[][,] BNL1_hexa;
        private double[][,] BNL_hexa;


        private void GetInitialGeometricDataAndInitializeMatrices(Element element)
        {
            ox_i = new double[8][];
            J_0b_hexa = new double[nGaussPoints][,];
            J_0_hexa = new double[nGaussPoints][,];
            J_0inv_hexa = new double[nGaussPoints][,];
            detJ_0 = new double[nGaussPoints];
            sunt_oloklhrwmatos = new double[nGaussPoints];
            BL11a_hexa = new double[nGaussPoints][,];
            BL12_hexa = new double[nGaussPoints][,];
            BL01_hexa = new double[nGaussPoints][,];
            BNL1_hexa = new double[nGaussPoints][,];
            BNL_hexa = new double[nGaussPoints][,];

            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                // initialize diastaseis twn mhtrwwn kai meta gemisma keliwn (olwn h mono oswn mporoume sthn arxh)
                J_0b_hexa[npoint] = new double[8,3];
                J_0_hexa[npoint] = new double[3,3];
                J_0inv_hexa[npoint] = new double[3,3];
                BL11a_hexa[npoint] = new double[6,9];
                BL12_hexa[npoint] = new double[9,9];
                BL01_hexa[npoint] = new double[6,9];
                BNL1_hexa[npoint] = new double[9,9];
                BNL_hexa[npoint] = new double[9,24];

                //
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0b_hexa[npoint][m, n] = ox_i[m][n];                     
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_0_hexa[npoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_0_hexa[npoint][m, n] += ll1_hexa[npoint][m, p] * J_0b_hexa[npoint][p, n];
                        }
                    }
                }

                //
                double det1 = J_0_hexa[npoint][0, 0] *
                         ((J_0_hexa[npoint][1, 1] * J_0_hexa[npoint][2, 2]) - (J_0_hexa[npoint][2, 1] * J_0_hexa[npoint][1, 2]));
                double det2 = J_0_hexa[npoint][0, 1] *
                              ((J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][2, 2]) - (J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][1, 2]));
                double det3 = J_0_hexa[npoint][0, 2] *
                              ((J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][2, 1]) - (J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][1, 1]));
                double jacobianDeterminant = det1 - det2 + det3;
                if (jacobianDeterminant < 0)
                {
                    throw new InvalidOperationException("The Jacobian Determinant is negative.");
                }
                detJ_0[npoint] = jacobianDeterminant;

                //
                J_0inv_hexa[npoint][0, 0] = ((J_0_hexa[npoint][1, 1] * J_0_hexa[npoint][2, 2]) - (J_0_hexa[npoint][2, 1] * J_0_hexa[npoint][1, 2])) *
                                    (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][0, 1] = ((J_0_hexa[npoint][2, 1] * J_0_hexa[npoint][0, 2]) - (J_0_hexa[npoint][0, 1] * J_0_hexa[npoint][2, 2])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][0, 2] = ((J_0_hexa[npoint][0, 1] * J_0_hexa[npoint][1, 2]) - (J_0_hexa[npoint][1, 1] * J_0_hexa[npoint][0, 2])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][1, 0] = ((J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][1, 2]) - (J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][2, 2])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][1, 1] = ((J_0_hexa[npoint][0, 0] * J_0_hexa[npoint][2, 2]) - (J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][0, 2])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][1, 2] = ((J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][0, 2]) - (J_0_hexa[npoint][0, 0] * J_0_hexa[npoint][1, 2])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][2, 0] = ((J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][2, 1]) - (J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][1, 1])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][2, 1] = ((J_0_hexa[npoint][2, 0] * J_0_hexa[npoint][0, 1]) - (J_0_hexa[npoint][2, 1] * J_0_hexa[npoint][0, 0])) *
                                        (1 / detJ_0[npoint]);
                J_0inv_hexa[npoint][2, 2] = ((J_0_hexa[npoint][0, 0] * J_0_hexa[npoint][1, 1]) - (J_0_hexa[npoint][1, 0] * J_0_hexa[npoint][0, 1])) *
                                        (1 / detJ_0[npoint]);

                //
                sunt_oloklhrwmatos[npoint] = detJ_0[npoint] * a_123g[npoint];

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL11a_hexa[npoint][m, n] = 0;
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL11a_hexa[npoint][m,3*m+ n] = J_0inv_hexa[npoint][m,n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL11a_hexa[npoint][3, n] = J_0inv_hexa[npoint][1, n]; // upologismos 4hs gramms
                    BL11a_hexa[npoint][3, 3+n] = J_0inv_hexa[npoint][0, n];
                    BL11a_hexa[npoint][4, 3+ n] = J_0inv_hexa[npoint][2, n]; // upologismos 5hs gramms
                    BL11a_hexa[npoint][4, 6 + n] = J_0inv_hexa[npoint][1, n];
                    BL11a_hexa[npoint][5, 0 + n] = J_0inv_hexa[npoint][2, n]; // upologismos 6hs gramms
                    BL11a_hexa[npoint][5, 6 + n] = J_0inv_hexa[npoint][0, n];
                }

                //
                for (int m = 0; m < 9; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL12_hexa[npoint][m, n] = 0;
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[npoint][m, 3 * m + n] = J_0inv_hexa[npoint][0, n];
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos grammwn 4 ews 6
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[npoint][3+m, 3 * m + n] = J_0inv_hexa[npoint][1, n];
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos grammwn 7 ews 8
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL12_hexa[npoint][6+m, 3 * m + n] = J_0inv_hexa[npoint][2, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL01_hexa[npoint][m, n] = 0;
                    }
                }
                for (int m = 0; m < 3; m++) // upologismos triwn prwtwn grammwn
                {
                    for (int n = 0; n < 3; n++)
                    {
                        BL01_hexa[npoint][m, 3 * m + n] = J_0inv_hexa[npoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    BL01_hexa[npoint][3, n] = J_0inv_hexa[npoint][1, n]; // upologismos 4hs gramms
                    BL01_hexa[npoint][3, 3 + n] = J_0inv_hexa[npoint][0, n];
                    BL01_hexa[npoint][4, 3 + n] = J_0inv_hexa[npoint][2, n]; // upologismos 5hs gramms
                    BL01_hexa[npoint][4, 6 + n] = J_0inv_hexa[npoint][1, n];
                    BL01_hexa[npoint][5, 0 + n] = J_0inv_hexa[npoint][2, n]; // upologismos 6hs gramms
                    BL01_hexa[npoint][5, 6 + n] = J_0inv_hexa[npoint][0, n];
                }

                //
                for (int m = 0; m < 9; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BNL1_hexa[npoint][m, n] = 0;
                    }
                }
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BNL1_hexa[npoint][3*m+n,3*m+p] = J_0inv_hexa[npoint][n, p];
                        }
                    }
                }

                //
                for (int m = 0; m < 9; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        BNL_hexa[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BNL_hexa[npoint][m, n] += BNL1_hexa[npoint][m, p]*BL13_hexa[npoint][p, n];
                        }
                    }
                }
            }
        }

        // implement update coordinate data
    }

    
}
