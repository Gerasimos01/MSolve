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
        private double[][] tx_i;
        private double[][] tu_i;
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
            tx_i = new double[8][];
            tu_i = new double[8][]; // apla initialized edw kai tpt allo
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
                tx_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                //tu_i[j] = new double[] { 0, 0, 0 }; den ananewnontai se afth th methodo ta mhtrwa pou periexoun tu_i
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

            this.InitializeMatrices();
        }


        private double[,] ll2; //einai anexarthto twn GP // initialize gia to update coordinate
        private double[,] J_1b;
        private double[][,] J_1;// exoume tosa [,] osa einai kai ta gpoints        
        private double[][,] BL11b;
        private double[][,] DG;
        private double[][,] GL;
        private double[][] GLvec;
        private double[][]Spkvec;

        private double[,] l_perisp; // voithitikos upologismwn // initialize gia to update forces
        private double[][] sunt_ol_Spkvec;
        private double[][,] BL11;
        private double[][,] BL1112sun01_hexa;
        private double[][,] BL;
        private double[][] fxk1; // fxk1 exei diastaseis omoiws me to shell8disp to Fxk dld mia nGausspoints+1 logw tou oti kratoume sthn
                                 // teleftaia extra thesh to athroisma

        private double[][,] sunt_ol_Spk; // initialize gia to updateKmatrices
        private double[,] sunt_ol_SPK_epi_BNL_hexa;
        private double[][,] kl_;
        private double[][,] knl_;
        private double[,] k_stoixeiou ;
        private double[][,] Cons_disp; // upologismos apo to updatematerial
        private double[,] sunt_ol_cons_disp; // voithitiko upologismwn
        private double[,] sunt_ol_cons_disp_epi_BL;


        private void InitializeMatrices()
        {
            // initialize gia to update coordinate 
            ll2 = new double[8, 3];
            J_1b = new double[8, 3];
            J_1 = new double[nGaussPoints][,];
            //BL11b= new double[nGaussPoints][,];
            DG = new double[nGaussPoints][,];
            GL = new double[nGaussPoints][,];
            GLvec = new double[nGaussPoints][];
            Spkvec = new double[nGaussPoints][];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                J_1[npoint] = new double[3, 3];
                //BL11b[npoint] = new double[9,9];                
                DG[npoint] = new double[3, 3];
                GL[npoint] = new double[3, 3];
                GLvec[npoint] = new double[6];
                Spkvec[npoint] = new double[6];
                // epiprosthtws mhdenizoueme gia efkolia edw osa apo ta parapanw the tha gemisoun plhrws 
                //for (int j = 0; j < 9; j++)
                //{
                //    for (int k = 0; k < 9; k++)
                //    {
                //        BL11b[npoint][j, k] = 0;
                //    }
                //}
            }

            // initialize gia to update forces
            l_perisp = new double[3, 3];
            sunt_ol_Spkvec = new double[nGaussPoints][];
            BL11 = new double[nGaussPoints][,];
            BL1112sun01_hexa = new double[nGaussPoints][,];
            BL = new double[nGaussPoints][,];
            fxk1 = new double[nGaussPoints+1][];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                sunt_ol_Spkvec[npoint] = new double[6];
                BL11[npoint] = new double[6, 9];
                BL1112sun01_hexa[npoint] = new double[6, 9];
                BL[npoint] = new double[6, 24];
                
            }
            for (int npoint = 0; npoint < nGaussPoints+1; npoint++)
            {
                fxk1[npoint] = new double[24];
            }

            // initialize gia to updateKmatrices
            sunt_ol_Spk = new double[nGaussPoints][,];
            sunt_ol_SPK_epi_BNL_hexa = new double[9,24];
            kl_ = new double[nGaussPoints + 1][,]; 
            knl_ = new double[nGaussPoints + 1][,]; 
            k_stoixeiou = new double[24,24]; ;
            Cons_disp = new double[nGaussPoints ][,] ; // upologismos apo to updatematerial
            sunt_ol_cons_disp = new double [6,6];
            sunt_ol_cons_disp_epi_BL = new double[6, 24];
            for (int npoint = 0; npoint < nGaussPoints ; npoint++)
            {
                sunt_ol_Spk[npoint] = new double[3,3];
                Cons_disp[npoint] = new double[6, 6];
            }
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                kl_[npoint] = new double[24,24];
                knl_[npoint] = new double[24, 24];
            }


        }



        private void UpdateCoordinateData(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    tu_i[j][k]= localdisplacements[3 * j + k];
                    tx_i[j][k] = ox_i[j][k] + tu_i[j][k];
                }
            }

            //
            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = tu_i[m][n];
                    J_1b[m,n] = tx_i[m][n];
                }
            }

            // //
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        J_1[npoint][m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            J_1[npoint][m, n] += ll1_hexa[npoint][m, p] * J_1b[p, n];
                        }
                    }
                }

                
                //
                //for (int m = 0; m < 3; m++)
                //{
                //    for (int n = 0; n < 3; n++)
                //    {
                //        for (int p = 0; p < 3; p++)
                //        {
                //            BL11b[npoint][3 * m + n, 3 * m + p] = l_perisp[n, p]; 
                //        }
                //    }
                //}

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        DG[npoint][m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            DG[npoint][m, n] += J_0inv_hexa[npoint][m, p] * J_1[npoint][n, p];
                        }
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        GL[npoint][m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            GL[npoint][m, n] += DG[npoint][p,m] * DG[npoint][p,n];
                        }
                    }
                }
                for (int m = 0; m < 3; m++)
                {
                    GL[npoint][m, m] += -1;
                }
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        GL[npoint][m, n] = 0.5 * GL[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    GLvec[npoint][m] = GL[npoint][m, m];
                }
                GLvec[npoint][3] = 2 * GL[npoint][0, 1];
                GLvec[npoint][4] = 2 * GL[npoint][1, 2];
                GLvec[npoint][5] = 2 * GL[npoint][2, 0];
            }

            }


        // apo uliko tha einai gnwsto to SpkVec[npoint][1:6] gia ola ta npoints
        // me vash afto programmatizontai oi forces


        private void UpdateForces()
        {


            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 6; m++)
                {
                    sunt_ol_Spkvec[npoint][m] = sunt_oloklhrwmatos[npoint] * Spkvec[npoint][m];
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        l_perisp[m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            l_perisp[m, n] += ll1_hexa[npoint][m, p] * ll2[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL11[npoint][m,n] = 0;
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_perisp[p, n];
                            BL11[npoint][m,3+n]+= BL11a_hexa[npoint][m, 3+p] * l_perisp[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_perisp[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL1112sun01_hexa[npoint][m, n] += BL11[npoint][m, p] * BL12_hexa[npoint][p, n];
                        }
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] += BL01_hexa[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        BL[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL[npoint][m, n] += BL1112sun01_hexa[npoint][m, p] * BL13_hexa[npoint][p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    fxk1[npoint][m] = 0;
                    for (int n = 0; n < 6; n++)
                    {
                        fxk1[npoint][m] += BL[npoint][n, m] * sunt_ol_Spkvec[npoint][n];
                    }
                }
            }

            //
            for (int m = 0; m < 24; m++)
            {
                fxk1[nGaussPoints][m] = 0;
            }
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 24; m++)
                {
                    fxk1[nGaussPoints][m] += fxk1[npoint][m];
                }
            }

        }

        private void InitializeBland_sunt_ol_Spkvec() //.first_calc_for_Kmatrices
        {
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {

                //
                for (int m = 0; m < 6; m++)
                {
                    sunt_ol_Spkvec[npoint][m] = sunt_oloklhrwmatos[npoint] * Spkvec[npoint][m];
                }

                //
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        l_perisp[m, n] = 0;
                        for (int p = 0; p < 8; p++)
                        {
                            l_perisp[m, n] += ll1_hexa[npoint][m, p] * ll2[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL11[npoint][m, n] = 0;
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            BL11[npoint][m, n] += BL11a_hexa[npoint][m, p] * l_perisp[p, n];
                            BL11[npoint][m, 3 + n] += BL11a_hexa[npoint][m, 3 + p] * l_perisp[p, n];
                            BL11[npoint][m, 6 + n] += BL11a_hexa[npoint][m, 6 + p] * l_perisp[p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL1112sun01_hexa[npoint][m, n] += BL11[npoint][m, p] * BL12_hexa[npoint][p, n];
                        }
                    }
                }
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 9; n++)
                    {
                        BL1112sun01_hexa[npoint][m, n] += BL01_hexa[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        BL[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            BL[npoint][m, n] += BL1112sun01_hexa[npoint][m, p] * BL13_hexa[npoint][p, n];
                        }
                    }
                }
            }

        }


        private void UpdateKmatrices()
        {
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //
                sunt_ol_Spk[npoint][0, 0] = sunt_ol_Spkvec[npoint][0];
                sunt_ol_Spk[npoint][0, 1] = sunt_ol_Spkvec[npoint][3];
                sunt_ol_Spk[npoint][0, 2] = sunt_ol_Spkvec[npoint][5];
                sunt_ol_Spk[npoint][1, 0] = sunt_ol_Spkvec[npoint][3];
                sunt_ol_Spk[npoint][1, 1] = sunt_ol_Spkvec[npoint][1];
                sunt_ol_Spk[npoint][1, 2] = sunt_ol_Spkvec[npoint][4];
                sunt_ol_Spk[npoint][2, 0] = sunt_ol_Spkvec[npoint][5];
                sunt_ol_Spk[npoint][2, 1] = sunt_ol_Spkvec[npoint][4];
                sunt_ol_Spk[npoint][2, 2] = sunt_ol_Spkvec[npoint][2];

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 6; n++)
                    {
                        sunt_ol_cons_disp[m, n] = sunt_oloklhrwmatos[npoint] * Cons_disp[npoint][m, n];
                    }
                }

                //
                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        sunt_ol_cons_disp_epi_BL[m, n] = 0;
                        for (int p = 0; p < 6; p++)
                        {
                            sunt_ol_cons_disp_epi_BL[m, n] += sunt_ol_cons_disp[m, p] * BL[npoint][p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        kl_[npoint][m, n] = 0;
                        for (int p = 0; p < 6; p++)
                        {
                            kl_[npoint][m, n] += BL[npoint][p, m] * sunt_ol_cons_disp_epi_BL[p, n];
                        }
                    }
                }
                //tha athroisoume meta ola ta kl- sthn teleftaia thesi

                //
                for (int m = 0; m < 3; m++) //prwtes 3x24 grammes
                {
                    for (int n = 0; n < 24; n++)
                    {
                        sunt_ol_SPK_epi_BNL_hexa[m, n] = 0;
                        sunt_ol_SPK_epi_BNL_hexa[3+m, n] = 0;
                        sunt_ol_SPK_epi_BNL_hexa[6 + m, n] = 0;
                        for (int p = 0; p < 3; p++)
                        {
                            sunt_ol_SPK_epi_BNL_hexa[m, n] += sunt_ol_Spk[npoint][m,p] * BNL_hexa[npoint][p, n];
                            sunt_ol_SPK_epi_BNL_hexa[3+m, n] += sunt_ol_Spk[npoint][m, p] * BNL_hexa[npoint][3+p, n];
                            sunt_ol_SPK_epi_BNL_hexa[6 + m, n] += sunt_ol_Spk[npoint][m, p] * BNL_hexa[npoint][6 + p, n];
                        }
                    }
                }

                //
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        knl_[npoint][m, n] = 0;
                        for (int p = 0; p < 9; p++)
                        {
                            knl_[npoint][m, n] += BNL_hexa[npoint][p, m] * sunt_ol_SPK_epi_BNL_hexa[p, n];
                        }
                    }
                }
                //tha athroisoume meta ola ta knl_ sthn teleftaia thesi i kateftheian sto k_stoixeiou

            }

            // athroisma olwn twn gpoints se k_stoixeiou kai prwta mhdenismos aftou
            for (int m = 0; m < 24; m++)
            {
                for (int n = 0; n < 24; n++)
                {
                    kl_[nGaussPoints][m, n] = 0;
                    knl_[nGaussPoints][m, n] = 0;
                }
            }
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 24; m++)
                {
                    for (int n = 0; n < 24; n++)
                    {
                        kl_[nGaussPoints][m, n] += kl_[npoint][m, n];
                        knl_[nGaussPoints][m, n] += knl_[npoint][m, n];
                    }
                }
            }
            for (int m = 0; m < 24; m++)
            {
                for (int n = 0; n < 24; n++)
                {
                    k_stoixeiou[m, n] = kl_[nGaussPoints][m, n] + knl_[nGaussPoints][m, n];
                }
            }
      }

        // telikes entoles kai mhtrwo mazas apo to hexa8

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements);
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                materialsAtGaussPoints[npoint].UpdateMaterial(GLvec[npoint]);
            }
            return new Tuple<double[], double[]>(GLvec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO mono to teleftaio dianusma tha epistrefei?
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                for (int j = 0; j < 6; j++)
                { Spkvec[npoint][j] = materialsAtGaussPoints[npoint].Stresses[j]; }
            }
            this.UpdateForces();
            return fxk1[nGaussPoints];
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public virtual IMatrix2D<double> StiffnessMatrix(Element element)
        {
            if (Cons_disp == null)
            {
                this.CalculateShapeFunctionAndGaussPointData();
                this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.UpdateCoordinateData(new double[24]);
                for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)// loop gia getfirstStressesFromMaterial
                {
                    for (int j = 0; j < 6; j++)
                    { Spkvec[npoint][j] = materialsAtGaussPoints[npoint].Stresses[j]; }
                }
                this.InitializeMatrices(); // meta to get twn stresses apo to material dioiti periexei ton pol/smo suntol epi Spkvec
                

            }
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                for (int j = 0; j < 6; j++)
                {
                    for (int k = 0; k < 6; k++)
                    { Cons_disp[npoint][j, k] = materialsAtGaussPoints[npoint].ConstitutiveMatrix[j, k]; }
                }
            }
            this.UpdateKmatrices();
            IMatrix2D<double> element_stiffnessMatrix = new Matrix2D<double>(k_stoixeiou); // TODO giati de ginetai return dof.Enumerator.GetTransformedMatrix, xrhsh symmetric
            return element_stiffnessMatrix;
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

        // omoiws me hexa 8 shell8disp implemented
        public int ID
        {
            get { return 13; }
        }
        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            return dofTypes;
        }

        //NOT IMPLEMENTED
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[24];
        }

        public virtual IMatrix2D<double> MassMatrix(Element element)
        {
            return new Matrix2D<double>(24, 24);
        }

        public virtual IMatrix2D<double> DampingMatrix(Element element)
        {

            return new Matrix2D<double>(24, 24);
        }
        //NOT IMPLEMENTED ews edw


    }

    
}
