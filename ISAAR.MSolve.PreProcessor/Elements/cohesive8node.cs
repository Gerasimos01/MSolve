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
    public class cohesive8node : IStructuralFiniteElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly IFiniteElementMaterial3D[] materialsAtGaussPoints;
        protected IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw

        public int gp_d1_coh { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public int gp_d2_coh { get; set; }
        private int nGaussPoints;

        protected cohesive8node()//consztructor apo to hexa8
        {
        }

        public cohesive8node(IFiniteElementMaterial3D material, int gp_d1c, int gp_d2c)
        {
            this.gp_d1_coh = gp_d1c;
            this.gp_d2_coh = gp_d2c;
            this.nGaussPoints = this.gp_d1_coh * this.gp_d2_coh;
            materialsAtGaussPoints = new IFiniteElementMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IFiniteElementMaterial3D)material.Clone();

        }

        public int endeixiShapeFunctionAndGaussPointData = 1;
        private double[] a_12g;
        private double a_1g;
        private double a_2g;
        private double[][] N1;
        private double[][,] N3;
        private double[][] N1_ksi;
        private double[][] N1_heta;
        private double[] N_i;
        private double[] N_i_ksi;
        private double[] N_i_heta;

        private double ksi;
        private double heta;
        private double zeta;
        private int npoint;

        private void CalculateShapeFunctionAndGaussPointData()
        {
            nGaussPoints = gp_d1_coh * gp_d2_coh ;
            a_12g = new double[nGaussPoints];
            N_i = new double[4]; 
            N_i_ksi = new double[4]; 
            N_i_heta = new double[4]; 
            N1 = new double[nGaussPoints][];
            N3 = new double[nGaussPoints][,];
            N1_ksi = new double[nGaussPoints][];
            N1_heta = new double[nGaussPoints][];
            for (int l = 0; l < nGaussPoints; l++)
            {
                N1[l] = new double[4];
                N3[l] = new double[3,12];
                N1_ksi[l] = new double[4];
                N1_heta[l] = new double[4];
            }
            for (int j = 0; j < gp_d1_coh; j++) 
            {
                    for (int k = 0; k < gp_d2_coh; k++)
                    {
                    npoint = j * gp_d1_coh + k;
                    if (gp_d1_coh == 3)
                    {
                        ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                        a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                    }
                    if (gp_d1_coh == 2)
                    {
                        ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                        a_1g = 1;
                    }
                    if (gp_d2_coh == 3)
                    {
                        heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                        a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                    }
                    if (gp_d2_coh == 2)
                    {
                        heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                        a_2g = 1;
                    }

                    a_12g[npoint] = a_1g * a_2g;

                    N_i[0] = 0.25 * (1 + ksi) * (1 + heta);
                    N_i[1] = 0.25 * (1 - ksi) * (1 + heta);
                    N_i[2] = 0.25 * (1 - ksi) * (1 - heta);
                    N_i[3] = 0.25 * (1 + ksi) * (1 - heta);

                    N_i_ksi[0] = +0.25 * (1 + heta);
                    N_i_ksi[1] = -0.25 * (1 + heta);
                    N_i_ksi[2] = -0.25 * (1 - heta);
                    N_i_ksi[3] = +0.25 * (1 - heta);

                    N_i_heta[0] = +0.25 * (1 + ksi);
                    N_i_heta[1] = +0.25 * (1 - ksi);
                    N_i_heta[2] = -0.25 * (1 - ksi);
                    N_i_heta[3] = -0.25 * (1 + ksi);

                    for (int l = 0; l < 4; l++)  // to 8 ginetai 4 gia to cohesive8node
                    { N1[npoint][l] = N_i[l]; }

                    for (int l = 0; l < 3; l++)  // arxika mhdenismos twn stoixweiwn tou pinaka
                    { for (int m = 0; m < 12; m++)
                        { N3[npoint][l, m] = 0; }
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 4; m++)
                        { N3[npoint][l, l+3*m] = N_i[m]; }
                    }

                    for (int l = 0; l < 4; l++)  
                    { N1_ksi[npoint][l] = N_i_ksi[l]; }

                    for (int l = 0; l < 4; l++)
                    { N1_heta[npoint][l] = N_i_heta[l]; }

                }
            }
            endeixiShapeFunctionAndGaussPointData = 2;
        }

        private double[] Get_a_12g()
        {
            if (endeixiShapeFunctionAndGaussPointData == 1)
            {
                CalculateShapeFunctionAndGaussPointData();
                return a_12g;
            }
            else
            { return a_12g; }
        }

        private double[][] GetN1()
        {
            if (endeixiShapeFunctionAndGaussPointData == 1)
            {
                CalculateShapeFunctionAndGaussPointData();
                return N1;
            }
            else
            { return N1; }
        }

        private double[][,] GetN3()
        {
            if (endeixiShapeFunctionAndGaussPointData == 1)
            {
                CalculateShapeFunctionAndGaussPointData();
                return N3;
            }
            else
            { return N3; }
        }

        private double[][] GetN1_ksi()
        {
            if (endeixiShapeFunctionAndGaussPointData == 1)
            {
                CalculateShapeFunctionAndGaussPointData();
                return N1_ksi;
            }
            else
            { return N1_ksi; }
        }

        private double[][] GetN1_heta()
        {
            if (endeixiShapeFunctionAndGaussPointData == 1)
            {
                CalculateShapeFunctionAndGaussPointData();
                return N1_heta;
            }
            else
            { return N1_heta; }
        }

        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths
        private double[][] tx_i; //8 arrays twn 3 stoixeiwn
        private double[] x_local; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
        private double[,] u_prok;
        private double[,] x_pavla;
        private double[,] k_stoixeiou_coh;
        private double[] fxk1_coh;
        private double[] d_trial;
        private double[] e_ksi;
        private double e_ksi_norm;
        private double[] e_heta;
        private double[] e_1;
        private double[] e_2;
        private double[] e_3;
        private double e_3_norm;
        private double[][,] R;
        private double[] u; // 3 epi 1
        private double[][] Delta; // [nGausspoints][3]
        private double[][,] D_tan; // [nGausspoints][3,3]
        private double[][] T_int;  //[nGausspoints][3]
        private double[][] c_1; // [nGausspoints][3]
        private double[] coh_det_J_t;
        private double[] sunt_olokl;
        private double[,] M; // 12 epi 12
        private double[] r_int; // 24 epi 1
        private double[] r_int_1; // to panw miso tou dianusmatos
        // gia tous pollaplasiasmous
        private double[][,] RN3;
        private double[,] D_tan_sunt_ol;
        private double[,] D_RN3_sunt_ol;
        private double[] T_int_sunt_ol;
        //temporary
        private int n_incr = 0;
        private int n_iter = 0;
        private void GetInitialGeometricDataAndInitializeMatrices(Element element)
        {
            ox_i = new double[8][];
            tx_i = new double[8][]; //mporei na mhn xreiasthei teilka!!!!
            for (int j = 0; j < 8; j++)
            {
                ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
                tx_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
            }

            x_local = new double[24];
            u_prok = new double[3, 4];
            x_pavla = new double[3, 4];
            k_stoixeiou_coh = new double[24, 24]; //allazei sto cohesive 8 node
            fxk1_coh = new double[24];
            d_trial = new double[nGaussPoints];
            e_ksi = new double[3];
            e_heta = new double[3];
            e_1 = new double[3];
            e_2 = new double[3];
            e_3 = new double[3];
            u = new double[3];
            nGaussPoints = gp_d1_coh * gp_d2_coh;
            Delta = new double[nGaussPoints][];
            R = new double[nGaussPoints][,];
            c_1 = new double[nGaussPoints][];
            coh_det_J_t = new double[nGaussPoints];
            sunt_olokl = new double[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                Delta[j] = new double[3];
                R[j] = new double[3, 3];
                c_1[j] = new double[3];
            }
            D_tan = new double[nGaussPoints][,];
            T_int = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                D_tan[j] = new double[3, 3];
                T_int[j] = new double[3];
            }
            M = new double[12, 12];
            r_int = new double[24];
            r_int_1 = new double[12];
            RN3 = new double[nGaussPoints][,];
            D_tan_sunt_ol = new double[3, 3];
            D_RN3_sunt_ol = new double[3, 12];
            T_int_sunt_ol = new double[3];
            for (int j = 0; j < nGaussPoints; j++)
            {
                RN3[j] = new double[3, 12];  //tha ginei [3,12] sto cohesive 8 node
            }

        }

        private void UpdateCoordinateData(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    x_local[3 * j + k] = ox_i[j][k] + localdisplacements[3 * j + k];
                }
            }

            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    u_prok[k, j] = x_local[k + 3 * j] - x_local[12 + k + 3 * j];
                    x_pavla[k, j] = x_local[k + 3 * j] + x_local[12 + k + 3 * j];
                }
            }
            // sunexeia ews upologismou tou Delta gia ola ta gp
            
                for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        e_ksi[l] = 0;
                        e_heta[l] = 0;
                        for (int m = 0; m < 4; m++) // tha ginei 4 sto cohesive 8 node
                        {
                            e_ksi[l] += GetN1_ksi()[npoint1][m] * x_pavla[l, m];
                            e_heta[l] += GetN1_heta()[npoint1][m] * x_pavla[l, m];
                        }
                        e_ksi[l] = 0.5 * e_ksi[l];
                        e_heta[l] = 0.5 * e_heta[l];
                    }
                    this.cross(e_ksi, e_heta, e_3);
                    e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
                    e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
                    for (int l = 0; l < 3; l++)
                    {
                        e_3[l] = e_3[l] / e_3_norm;
                        e_1[l] = e_ksi[l] / e_ksi_norm;
                    }
                    this.cross(e_1, e_3, e_2);
                    for (int l = 0; l < 3; l++)
                    {
                        R[npoint1][l, 0] = e_1[l];
                        R[npoint1][l, 1] = e_2[l];
                        R[npoint1][l, 2] = e_3[l];

                    }
                    for (int l = 0; l < 3; l++)
                    { u[l] = 0; }
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 4; m++)  // pithanws gia to cohesive 8 node na gineiews 4 to m
                        {
                            u[l] += u_prok[l, m] * GetN1()[npoint1][m];
                        }
                    }
                    for (int l = 0; l < 3; l++)
                    { Delta[npoint1][l] = 0; }
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            Delta[npoint1][l] += R[npoint1][m, l] * u[m];
                        }
                    }

                    this.cross(e_ksi, e_heta, c_1[npoint1]);
                    coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
                    sunt_olokl[npoint1] = coh_det_J_t[npoint1] * Get_a_12g()[npoint1];
                }
            


        }

        private void cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }

        private void InitializeRN3()
        {

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            RN3[npoint1][l, m] = 0;
                        }
                    }
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            { RN3[npoint1][l, m] += R[npoint1][l, n] * GetN3()[npoint1][n, m]; }
                        }
                    }
            }
        }

        private void UpdateForces()
        {
            for (int j = 0; j < 24; j++) // allagh sto cohesive 8 node
            {
                fxk1_coh[j] = 0;
            }

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                    
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            RN3[npoint1][l, m] = 0;
                        }
                    }
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            { RN3[npoint1][l, m] += R[npoint1][l, n] * GetN3()[npoint1][n, m]; }
                        }
                    }
                    for (int l = 0; l < 3; l++)
                    {
                        T_int_sunt_ol[l] = T_int[npoint1][l] * sunt_olokl[npoint1];
                    }
                    for (int l = 0; l < 12; l++)
                    { r_int_1[l] = 0; }
                    for (int l = 0; l < 12; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        { r_int_1[l] += RN3[npoint1][m, l] * T_int_sunt_ol[m]; }
                    }
                    for (int l = 0; l < 12; l++)
                    {
                        fxk1_coh[l] += r_int_1[l];
                        fxk1_coh[12 + l] += (-r_int_1[l]);
                    }
                }
        }

        private void UpdateKmatrices()
        {
            for (int k = 0; k < 24; k++) // allagh sto cohesive 8 node
            {
                for (int j = 0; j < 24; j++)
                {
                    k_stoixeiou_coh[k, j] = 0;
                }
            }
            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            D_tan_sunt_ol[l, m] = D_tan[npoint1][l, m] * sunt_olokl[npoint1];
                        }
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)  // omoiws to 24 gia ta cohesive 8 node
                        {
                            D_RN3_sunt_ol[l, m] = 0;
                        }
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            { D_RN3_sunt_ol[l, m] += D_tan_sunt_ol[l, n] * RN3[npoint1][n, m]; }
                        }
                    }
                    for (int l = 0; l < 12; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            M[l, m] = 0;
                        }
                    }
                    for (int l = 0; l < 12; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                M[l, m] += RN3[npoint1][n, l] * D_RN3_sunt_ol[n, m];
                            }
                        }
                    }
                    for (int l = 0; l < 12; l++)
                    {
                        for (int m = 0; m < 12; m++)
                        {
                            k_stoixeiou_coh[l, m] += M[l, m];
                            k_stoixeiou_coh[l, 12 + m] += -M[l, m];
                            k_stoixeiou_coh[12 + l, m] += -M[l, m];
                            k_stoixeiou_coh[12 + l, 12 + m] += M[l, m];
                        }
                    }


                }
            
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements);
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
            }
            return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO mono to teleftaio dianusma tha epistrefei?
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                for (int j = 0; j < 3; j++)
                { T_int[i][j] = materialsAtGaussPoints[i].Stresses[j]; }
            }
            this.UpdateForces();
            //temporary
            n_iter += 1;
            if (n_iter == 1)
            { n_iter += 0; }

            return fxk1_coh;           
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public virtual IMatrix2D<double> StiffnessMatrix(Element element)
        {
            if (D_tan == null)
            {
                this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.UpdateCoordinateData(new double[24]);
                this.InitializeRN3();
            }

            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    { D_tan[i][j, k] = materialsAtGaussPoints[i].ConstitutiveMatrix[j, k]; }
                }
            }
            this.UpdateKmatrices();
            IMatrix2D<double> element_stiffnessMatrix = new Matrix2D<double>(k_stoixeiou_coh); // TODO giati de ginetai return dof.Enumerator.GetTransformedMatrix, xrhsh symmetric
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
            //temporary
            n_incr += 1;
            if (n_incr == 17)
            { n_incr += 0; }
        }

        public void ClearMaterialStresses()
        {
            foreach (IFiniteElementMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        // omoiws me hexa 8 shell8disp implemented
        public int ID
        {
            get { return 14; }
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

        // aplopoihtika implemented mhdenikes masses gia cohesive - not implemented
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[24];
        }

        public virtual IMatrix2D<double> MassMatrix(Element element)
        {
            return new Matrix2D<double>(24,24);
        }

        public virtual IMatrix2D<double> DampingMatrix(Element element)
        {

            return new Matrix2D<double>(24, 24);
        }

    }

}
