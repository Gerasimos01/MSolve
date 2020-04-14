﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using System.Linq;

namespace ISAAR.MSolve.FEM.Elements
{
    public class cohesiveElement : IStructuralFiniteElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        protected readonly IDofType[][] dofTypes;
        private readonly IIsoparametricInterpolation2D interpolation;
        protected readonly ICohesiveZoneMaterial3D[] materialsAtGaussPoints;
        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator(); 
        // ews edw


        private int nGaussPoints; //TODO meta afto mporei na ginei readonly.

        public bool MatrixIsNotInitialized = true;

        protected cohesiveElement()//consztructor apo to hexa8
        {
        }

        public cohesiveElement(ICohesiveZoneMaterial3D material, IQuadrature2D quadratureForStiffness,IIsoparametricInterpolation2D interpolation)
        {
            this.QuadratureForStiffness = quadratureForStiffness;
            this.interpolation = interpolation;
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;

           materialsAtGaussPoints = new ICohesiveZoneMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = material.Clone();

            dofTypes = new IDofType[interpolation.NumFunctions*2][];
            for (int i = 0; i < interpolation.NumFunctions*2; i++)
            {
                dofTypes[i] = new IDofType[]
                {
                    StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
                };
            }
        }

        public IReadOnlyList<IFiniteElementMaterial> Materials => materialsAtGaussPoints;

        public CellType CellType { get; } = CellType.Unknown;

        public IQuadrature2D QuadratureForStiffness { get; }


        public int endeixiShapeFunctionAndGaussPointData = 1;
        

       

        

        

        //private double[][,] GetN3()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return N3;
        //    }
        //    else
        //    { return N3; }
        //}

        

        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths
        private double[] x_local; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
        //private double[][,] R;
        //private double[][,] D_tan; // [nGausspoints][3,3]
        //private double[][] T_int;  //[nGausspoints][3]
        //private double[][] c_1; // [nGausspoints][3]
        //private double[] coh_det_J_t;
        //private double[] sunt_olokl;
        //private double[,] M; // 12 epi 12
        //private double[] r_int; // 24 epi 1
        //private double[] r_int_1; // to panw miso tou dianusmatos
        // gia tous pollaplasiasmous
        //private double[][,] RN3;
        //private double[,] D_tan_sunt_ol;
        //private double[,] D_RN3_sunt_ol;
        //private double[] T_int_sunt_ol;
        //temporary
        //private int n_incr = 0;
        //private int n_iter = 0;
        private void GetInitialGeometricDataAndInitializeMatrices(IElement element)
        {
            ox_i = new double[2*interpolation.NumFunctions][];
            for (int j = 0; j < 2*interpolation.NumFunctions; j++)
            {
                ox_i[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
            }

            x_local = new double[2*3*interpolation.NumFunctions];
            

        }

        private double[][] UpdateCoordinateData(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);

            double[,] u_prok = new double[3, 2*interpolation.NumFunctions];
            double[,] x_bar = new double[3, 2*interpolation.NumFunctions];

            double[] e_ksi = new double[3];
            double e_ksi_norm;
            double[] e_heta = new double[3];
            double[] e_1 = new double[3];
            double[] e_2 = new double[3];
            double[] e_3 = new double[3];
            double e_3_norm;
            double[] u = new double[3];

            //double[] coh_det_J_t = new double[nGaussPoints];

            double[][] Delta = new double[nGaussPoints][];
            double[][] c_1 = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                Delta[j] = new double[3];
                c_1[j] = new double[3];
            }
            double[][,] R = new double[nGaussPoints][,]; //TODO: maybe cache R
            for (int j = 0; j < nGaussPoints; j++)
            {
                R[j] = new double[3, 3];
            }

            for (int j = 0; j < 2*interpolation.NumFunctions; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    x_local[3 * j + k] = ox_i[j][k] + localdisplacements[3 * j + k];
                }
            }

            for (int j = 0; j < interpolation.NumFunctions; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    u_prok[k, j] = x_local[k + 3 * j] - x_local[3*interpolation.NumFunctions + k + 3 * j];
                    x_bar[k, j] = x_local[k + 3 * j] + x_local[3*interpolation.NumFunctions + k + 3 * j];
                }
            }
            // sunexeia ews upologismou tou Delta gia ola ta gp

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                for (int l = 0; l < 3; l++)
                {
                    e_ksi[l] = 0;
                    e_heta[l] = 0;
                    for (int m = 0; m < interpolation.NumFunctions; m++) // tha ginei 4 sto cohesive 8 node
                    {
                        e_ksi[l] += shapeFunctionDerivatives[npoint1][m,0] * x_bar[l, m];
                        e_heta[l] += shapeFunctionDerivatives[npoint1][m,1] * x_bar[l, m];
                    }
                    e_ksi[l] = 0.5 * e_ksi[l];
                    e_heta[l] = 0.5 * e_heta[l];
                }
                this.Cross(e_ksi, e_heta, e_3);
                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
                e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
                for (int l = 0; l < 3; l++)
                {
                    e_3[l] = e_3[l] / e_3_norm;
                    e_1[l] = e_ksi[l] / e_ksi_norm;
                }
                this.Cross(e_1, e_3, e_2);
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
                    for (int m = 0; m < interpolation.NumFunctions; m++)  // pithanws gia to cohesive 8 node na gineiews 4 to m
                    {
                        u[l] += u_prok[l, m] * N1[npoint1][m];
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

                this.Cross(e_ksi, e_heta, c_1[npoint1]);
                //coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
                //sunt_olokl[npoint1] = coh_det_J_t[npoint1] * Get_a_12g()[npoint1];
            }

            return Delta;

        }

        private double[,] ReShapeShapeFunctionDataMatrix(double[] N1)
        {
            var N3 = new double[3, 3 * interpolation.NumFunctions];
            for (int l = 0; l < 3; l++)
            {
                for (int m = 0; m < interpolation.NumFunctions; m++)
                { N3[l, l + 3 * m] = N1[m]; }
            }
            return N3;
        }
        private Tuple<Matrix[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations()
        {
            IReadOnlyList<double[]> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix> N3 = N1.Select(x => Matrix.CreateFromArray(ReShapeShapeFunctionDataMatrix(x))).ToList();
            IReadOnlyList<Matrix> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            double[] integrationsCoeffs = new double[nGaussPoints];
            Matrix[] RtN3 = new Matrix[nGaussPoints];
            double[,] x_bar = new double[3, 8];

            double[] e_1 = new double[3];
            double[] e_2 = new double[3];
            double[] e_3 = new double[3];
            double e_3_norm;

            double[] coh_det_J_t = new double[nGaussPoints];

            double[][] c_1 = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                c_1[j] = new double[3];
            }

            Matrix[] R = new Matrix[nGaussPoints]; //TODO: perhaps cache matrices in InitializeMatrices() where RtN3 is calculated
            for (int j = 0; j < nGaussPoints; j++)
            { R[j] = Matrix.CreateZero(3, 3); }

            for (int j = 0; j < interpolation.NumFunctions; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    x_bar[k, j] = x_local[k + 3 * j] + x_local[3*interpolation.NumFunctions + k + 3 * j];
                }
            }

            // Calculate Delta for all GPs
            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                double[] e_ksi = new double[3];
                double[] e_heta = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < interpolation.NumFunctions; m++) // must be 4 in cohesive 8-nodes
                    {
                        e_ksi[l] += shapeFunctionDerivatives[npoint1][m,0] * x_bar[l, m];
                        e_heta[l] += shapeFunctionDerivatives[npoint1][m,1] * x_bar[l, m];
                    }
                    e_ksi[l] = 0.5 * e_ksi[l];
                    e_heta[l] = 0.5 * e_heta[l];
                }
                this.Cross(e_ksi, e_heta, e_3);
                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
                double e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
                for (int l = 0; l < 3; l++)
                {
                    e_3[l] = e_3[l] / e_3_norm;
                    e_1[l] = e_ksi[l] / e_ksi_norm;
                }
                this.Cross(e_1, e_3, e_2);
                for (int l = 0; l < 3; l++)
                {
                    R[npoint1][l, 0] = e_1[l];
                    R[npoint1][l, 1] = e_2[l];
                    R[npoint1][l, 2] = e_3[l];

                }

                this.Cross(e_ksi, e_heta, c_1[npoint1]);
                coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
                integrationsCoeffs[npoint1] = coh_det_J_t[npoint1] * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;

                // Calculate RtN3 here instead of in InitializeRN3() and then in UpdateForces()
                RtN3[npoint1] = R[npoint1].Transpose() * N3[npoint1];
            }
            return new Tuple<Matrix[], double[]>(RtN3, integrationsCoeffs);
        }
        private void Cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }

        

        private double[] UpdateForces(Element element, Matrix[] RtN3, double[] integrationCoeffs)
        {
            double[] fxk1_coh = new double[3*2*interpolation.NumFunctions]; // TODO: must be 24 in cohesive 8 node

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                double[] T_int_integration_coeffs = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    T_int_integration_coeffs[l] = materialsAtGaussPoints[npoint1].Tractions[l] * integrationCoeffs[npoint1];
                }

                double[] r_int_1 = new double[3*interpolation.NumFunctions];
                for (int l = 0; l < 3 * interpolation.NumFunctions; l++)
                {
                    for (int m = 0; m < 3; m++)
                    { r_int_1[l] += RtN3[npoint1][m, l] * T_int_integration_coeffs[m]; }
                }
                for (int l = 0; l < 3 * interpolation.NumFunctions; l++)
                {
                    fxk1_coh[l] += r_int_1[l];
                    fxk1_coh[3 * interpolation.NumFunctions + l] += (-r_int_1[l]);
                }
            }

            return fxk1_coh;
        }

        private double[,] UpdateKmatrices(IElement element, Matrix[] RtN3, double[] integrationCoeffs)
        {
            double[,] k_stoixeiou_coh = new double[3 * 2 * interpolation.NumFunctions, 3 * 2 * interpolation.NumFunctions];


            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                Matrix D_tan_sunt_ol = Matrix.CreateZero(3, 3);
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        D_tan_sunt_ol[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * integrationCoeffs[npoint1];// D_tan[npoint1][l, m] * integrationCoeffs[npoint1];
                    }
                }

                Matrix D_RtN3_sunt_ol = D_tan_sunt_ol * RtN3[npoint1];
                Matrix M = RtN3[npoint1].Transpose() * D_RtN3_sunt_ol;

                for (int l = 0; l < 3 *  interpolation.NumFunctions; l++)
                {
                    for (int m = 0; m < 3 * interpolation.NumFunctions; m++)
                    {
                        k_stoixeiou_coh[l, m] += M[l, m];
                        k_stoixeiou_coh[l, 3 * interpolation.NumFunctions + m] += -M[l, m];
                        k_stoixeiou_coh[3 * interpolation.NumFunctions + l, m] += -M[l, m];
                        k_stoixeiou_coh[3 * interpolation.NumFunctions + l, 3 * interpolation.NumFunctions + m] += M[l, m];
                    }
                }
            }

            return k_stoixeiou_coh;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            double[][] Delta =UpdateCoordinateData(localTotalDisplacements);
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
            }
            return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Tractions);
            //TODO mono to teleftaio dianusma tha epistrefei?
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
            Matrix[] RtN3;
            RtN3 = RtN3AndIntegrationCoeffs.Item1;
            double[] integrationCoeffs;
            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

            double[] fxk2_coh = UpdateForces(element, RtN3, integrationCoeffs);
            return fxk2_coh;           
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public virtual IMatrix StiffnessMatrix(IElement element)
        {
            double[,] k_stoixeiou_coh2;
            if (MatrixIsNotInitialized)
            {
                this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.UpdateCoordinateData(new double[2*3*interpolation.NumFunctions]); //returns Delta that can't be used for the initial material state
                MatrixIsNotInitialized = false;
            }

            Tuple<Matrix[], double[]> RtN3AndIntegrationCoeffs;
            RtN3AndIntegrationCoeffs = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
            Matrix[] RtN3;
            RtN3 = RtN3AndIntegrationCoeffs.Item1;
            double[] integrationCoeffs;
            integrationCoeffs = RtN3AndIntegrationCoeffs.Item2;

            k_stoixeiou_coh2 = this.UpdateKmatrices(element, RtN3, integrationCoeffs);
            IMatrix element_stiffnessMatrix = Matrix.CreateFromArray(k_stoixeiou_coh2);
            return element_stiffnessMatrix; // embedding
        }


        public bool MaterialModified
        {
            get
            {
                foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.SaveState();
            ////temporary
            //n_incr += 1;
            //if (n_incr == 17)
            //{ n_incr += 0; }
        }

        public void ClearMaterialStresses()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
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

        public IElementDofEnumerator DOFEnumerator
        {
            get { return DofEnumerator; }
            set { DofEnumerator = value; }
        }

        public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

        // aplopoihtika implemented mhdenikes masses gia cohesive - not implemented
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[2*3*interpolation.NumFunctions];
        }

        public virtual IMatrix MassMatrix(IElement element)
        {
            return Matrix.CreateZero(3*2*interpolation.NumFunctions,3*2*interpolation.NumFunctions);
        }

        public virtual IMatrix DampingMatrix(IElement element)
        {

            return Matrix.CreateZero(3 * 2 * interpolation.NumFunctions, 3 * 2 * interpolation.NumFunctions);
        }

    }

}
