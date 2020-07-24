using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// Isotropic elastic shell material that accounts for no out of plane shear deformation
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class ShellElasticMaterial2DtransformationbDefGrad : IShellMaterial
    {
		public double[] NormalVectorV3 { get; set; }
		public double[] TangentVectorV1 { get; set; }
		public double[] TangentVectorV2 { get; set; }
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
        Matrix transformationMatrix; // gia to shell

        private bool modified; 
		private double[,] CartesianConstitutiveMatrix;
		private double[] CartesianStresses = new double[6];

		object ICloneable.Clone() => Clone();

		public IShellMaterial Clone()
		{
			return new ShellElasticMaterial2DtransformationbDefGrad()
			{
				YoungModulus = this.YoungModulus,
				PoissonRatio = this.PoissonRatio,
			};
		}

		public void UpdateMaterial(double[] cartesianStrains) //TODO: rename cartesian strains to strains 
		{
			if (CartesianConstitutiveMatrix == null)
			{
				this.CalculateConstitutiveMatrix(Vector.CreateFromArray(TangentVectorV1), Vector.CreateFromArray(TangentVectorV2));
			}

			for (int l = 0; l < 3; l++)
			{
				CartesianStresses[l] = 0;
				for (int m = 0; m < 3; m++)
				{
					CartesianStresses[l] += CartesianConstitutiveMatrix[l, m] * cartesianStrains[m];
				}
			}
		}

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tgi">Each Column is a Vector and not normalised.</param>
        /// <param name="tG_i">Each Column is a Vector and not normalised.</param>
        /// <param name="F"></param>
        /// <returns></returns>
        public (double[,], double[]) CalculateTransformations(double[,]tgi, double[,] tG_i, double[,] F_3D)
        {
            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;

            // normalise vectorsi.
            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = tgi[0, i1] * tgi[0, i1] + tgi[1, i1] * tgi[1, i1] + tgi[2, i1] * tgi[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    tgi[i2, i1] = tgi[i2, i1] / norm;
                }

            }
            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = tG_i[0, i1] * tG_i[0, i1] + tG_i[1, i1] * tG_i[1, i1] + tG_i[2, i1] * tG_i[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    tG_i[i2, i1] = tG_i[i2, i1] / norm;
                }

            }

            var Qij = CalculateRotationMatrix(tG_i, eye);
            var Qij1 = CalculateRotationMatrix(tgi, eye);

            double[,] F_rve = Transform_F3D_to_Frve(F_3D, Qij, Qij1);

            

            double[] GLvec = Transform2DDefGradToGL(F_rve);
            (double[] SPKvec, double[,] ConsCartes )= CalculateSPK(GLvec);

            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            double[,] FPKrve = TransformSPKvecToFPK(F_rve, SPKMat );

            double[,] Cinpk = ExpandCijrs(ConsCartes);

            double[,] Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);


            FPKrve = new double[3, 3] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } };
                  
            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPKrve, Qij, Qij1);

            var Qpi = CalculateRotationMatrix(eye, tgi);
            var Qqj = CalculateRotationMatrix(eye, tG_i);
            var Qrk = Qpi;
            var Qsl = Qqj;

            var Aijkl_3D = Transform_Aijkl_rve_to_Aijkl_3D(Aijkl_rve, Qpi, Qqj, Qrk, Qsl);

            double[] FPK_3D_vec = new double[9] { FPK_3D[0, 0], FPK_3D[1, 1], FPK_3D[2, 2], FPK_3D[0, 1], FPK_3D[1, 2], FPK_3D[2, 0], FPK_3D[0, 2], FPK_3D[1, 0], FPK_3D[2, 1]  };

            return (Aijkl_3D, FPK_3D_vec);

        }

        

        private double[,] TransformCinpk(double[,] Cinpk, double[,] F_rve, double[,] SPKMat)
        {

            double[,] Aijkl = new double[4, 4];

            int[,] thesi = new int[2, 2] { { 0, 2 }, { 3, 1 } };
            int row;//anaferetai sto Aijkl mhtrwo
            int col;// anaferetai sto Aijkl mhtrwo
            int rowC;// anaferetai sto Cinpk mhtrwo
            int colC;// anaferetai sto Cinpk mhtrwo

            int[,] dlj = new int[2,2] { { 1,0 }, { 0,1 } };

            // exwterika loop einai oi theseis pou gemizoun sto Aijkl
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    row = thesi[i, j];
                    for (int k = 0; k < 2; k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            col = thesi[k, l];

                            //eswterika loop einai to athroisma logw n kai p
                            for (int n = 0; n < 2; n++)
                            {
                                rowC = thesi[i, n];
                                for (int p = 0; p < 2; p++)
                                {
                                    colC = thesi[p, k];

                                    Aijkl[row, col] += Cinpk[rowC, colC] * F_rve[j, n] * F_rve[l, p] + SPKMat[i, k] * dlj[l, j];
                                }
                            }

                        }
                    }
                }
            }

            return Aijkl;

        }

        public double[,] CalculateRotationMatrix(double[,] e_new, double[,] e_old)
        {
            double[,] Qij = new double[3, 3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        Qij[i1, i2] += + e_new[k1, i2] * e_old[k1, i1];
                    }
                }
            }

            return Qij;

        }

        public double[] Transform2DDefGradToGL(double[,] F )
        {
            //double[,] GL = new double[2, 2];

            double[,] FtrF = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    for (int k1 = 0; k1 < 2; k1++)
                    {
                        FtrF[i1, i2] += F[k1, i1] * F[k1, i2];
                    }
                }
            }

            FtrF[0, 0] += -1;
            FtrF[1, 1] += -1;

            var GL = new double[2, 2] { { 0.5 * FtrF[0, 0], 0.5 * FtrF[0, 1] }, { 0.5 * FtrF[1, 0], 0.5 * FtrF[1, 1] } };

            double[] GLvec = new double[3] { GL[0, 0], GL[1, 1], 2 * GL[1, 0], };

            return GLvec;
                       
        }

        public (double[], double[,]) CalculateSPK(double[] GLvec)
        {
            var OriginalConstitutiveMatrix = new double[3, 3];
            //if (StressState == StressState2D.PlaneStress)
            {
                double aux = YoungModulus / (1 - PoissonRatio * PoissonRatio);
                OriginalConstitutiveMatrix[0, 0] = aux;
                OriginalConstitutiveMatrix[1, 1] = aux;
                OriginalConstitutiveMatrix[0, 1] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[1, 0] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[2, 2] = (1 - PoissonRatio) / 2 * aux;
            }

            var SPKvec = new double[3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    SPKvec[i1] += OriginalConstitutiveMatrix[i1, i2] * GLvec[i2];
                }
            }

            return (SPKvec, OriginalConstitutiveMatrix);
        }

        private double[,] TransformSPKvecToFPK(double[,] F_rve, double[,] SPKMat)
        {
            var FPKmat = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    for (int k1 = 0; k1 < 2; k1++)
                    {
                        FPKmat[i1, i2] += F_rve[i1,k1] * SPKMat[k1, i2];
                    }
                }
            }
            return FPKmat;
        }

        private double[,] ExpandCijrs(double[,] consCartes)
        {
            return new double[4, 4] {{consCartes[0,0], consCartes[0, 1], 0, 0 },
                                    {consCartes[1,0], consCartes[1, 1], 0, 0 },
                                    {0,0,consCartes[2,2], consCartes[2,2] },
                                    {0,0,consCartes[2,2], consCartes[2,2] }};
        }

        private double[,] Transform_F3D_to_Frve(double[,] F_3D, double[,] Qij, double[,] Qij1)
        {
            var F_3D_times_Qij = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        F_3D_times_Qij[i1, i2] += F_3D[i1, k1] * Qij[k1, i2];
                    }
                }
            }

            var F_rve = new double[2, 2];
            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        F_rve[i1, i2] += Qij1[k1, i1] * F_3D_times_Qij[k1, i2];
                    }
                }
            }

            return F_rve;
        }

        private double[,] Transform_FPK_rve_To_FPK_3D(double[,] FPKrve, double[,] Qij, double[,] Qij1)
        {
            var FPK_rve_times_QijT = new double[3, 3]; // ta loops tou transformation tou FPK apo rve se 3D mporoun na veltistopoiithoun. to FPK_rve einai 2x2
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        FPK_rve_times_QijT[i1, i2] += FPKrve[i1, k1] * Qij[i2, k1];
                    }
                }
            }

            var FPK_3D = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        FPK_3D[i1, i2] += Qij1[i1, k1] * FPK_rve_times_QijT[k1, i2];
                    }
                }
            }

            return FPK_3D;
        }

        private double[,] Transform_Aijkl_rve_to_Aijkl_3D(double[,] Apqrs_rve, double[,] Qpi, double[,] Qqj, double[,] Qrk, double[,] Qsl)
        {
            double[,] Aijkl_3D = new double[9, 9];

            int[,] thesi_pqrs = new int[2, 2] { { 0, 2 }, { 3, 1 } };  // 11 22 12 21 
            int row_pq;//anaferetai sto Apqrs_rve mhtrwo
            int col_rs;// anaferetai sto Apqrs_rve mhtrwo

            int[,] thesi_ijkl = new int[3, 3] { { 0, 3, 6 }, { 7, 1, 4 }, { 5, 8, 2 } }; // 11 22 33 12 23 31 13 21 32
            int row_ij;//anaferetai sto Aijkl_3D mhtrwo
            int col_kl;// anaferetai sto Aijkl_3D mhtrwo

            // exwterika loop einai oi theseis pou gemizoun sto Aijkl
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    row_ij = thesi_ijkl[i, j];
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            col_kl = thesi_ijkl[k, l];

                            //eswterika loop einai to athroisma logw p,q kai r,s
                            for (int p = 0; p < 2; p++)
                            {

                                for (int q = 0; q < 2; q++)
                                {
                                    row_pq = thesi_pqrs[p, q];

                                    for (int r = 0; r < 2; r++)
                                    {
                                        for (int s = 0; s < 2; s++)
                                        {
                                            col_rs = thesi_pqrs[r, s];

                                            Aijkl_3D[row_ij, col_kl] += Qpi[p, i] * Qqj[q, j] * Qrk[r, k] * Qsl[s, l] * Apqrs_rve[row_pq, col_rs];
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }

            return Aijkl_3D;
        }










        private void CalculateConstitutiveMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
            this.CalculateTransformationMatrix(Vector.CreateFromArray(TangentVectorV1), Vector.CreateFromArray(TangentVectorV2));

            var OriginalConstitutiveMatrix = new double[3, 3];
            //if (StressState == StressState2D.PlaneStress)
            {
                double aux = YoungModulus / (1 - PoissonRatio * PoissonRatio);
                OriginalConstitutiveMatrix[0, 0] = aux;
                OriginalConstitutiveMatrix[1, 1] = aux;
                OriginalConstitutiveMatrix[0, 1] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[1, 0] = PoissonRatio * aux;
                OriginalConstitutiveMatrix[2, 2] = (1 - PoissonRatio) / 2 * aux;
            }



            //         var auxMatrix1 = new Matrix2D(2, 2);
            //auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            //auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            //auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            //auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            //(Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant();

            //var constitutiveMatrix = new Matrix2D(new double[3, 3]
            //{
            //	{
            //		inverse[0,0]*inverse[0,0],
            //		this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
            //		inverse[0,0]*inverse[1,0]
            //	},
            //	{
            //		this.PoissonRatio*inverse[0,0]*inverse[1,1]+(1-this.PoissonRatio)*inverse[1,0]*inverse[1,0],
            //		inverse[1,1]*inverse[1,1],
            //		inverse[1,1]*inverse[1,0]
            //	},
            //	{
            //		inverse[0,0]*inverse[1,0],
            //		inverse[1,1]*inverse[1,0],
            //		0.5*(1-this.PoissonRatio)*inverse[0,0]*inverse[1,1]+(1+this.PoissonRatio)*inverse[1,0]*inverse[1,0]
            //	},
            //});

            //         // Integrate over thickness takes into account multiplication *t but not (E/(1-(ni^2)) and it will be added here
            //         ConstitutiveMatrix.Scale(YoungModulus / (1 - Math.Pow(PoissonRatio, 2)));
            var constitutiveMatrix = (transformationMatrix.Transpose() * (Matrix.CreateFromArray(OriginalConstitutiveMatrix)) * transformationMatrix);
            CartesianConstitutiveMatrix = constitutiveMatrix.CopyToArray2D();
		}

        private void CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = Matrix.CreateZero(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            Matrix inverse = auxMatrix1.Invert(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)
                                                  //TODO: auxMatrix1.Invert2x2AndDeterminant(1e-20) for bad geometry

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            Ei[0] =surfaceBasisVector1.CopyToArray();
            double G1_norm = surfaceBasisVector1.Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = (Vector.CreateFromArray(Ei[1])).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = Vector.CreateFromArray(Ei[i1]).DotProduct(Vector.CreateFromArray(G_i[i2]));
                }
            }

            transformationMatrix = Matrix.CreateFromArray(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });
        }

        private bool CheckIfConstitutiveMatrixChanged()
		{
			return false;
		}

		public double[] Stresses 
		{
			get { return CartesianStresses; }
		}

		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (CartesianConstitutiveMatrix == null) UpdateMaterial(new double[6]);
				return Matrix.CreateFromArray(CartesianConstitutiveMatrix);
			}
		}

		public void SaveState()
		{
		}

		public bool Modified => modified;

		public void ResetModified()
		{
			modified = false;
		}

		public int ID
		{
			get { throw new NotImplementedException(); }
		}

		public void ClearState()
		{
		}
		public void ClearStresses()
		{

		}

		public double[] Coordinates
		{

			get { throw new NotImplementedException(); }
			set { throw new InvalidOperationException(); }
		}

	}
}
