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
        public (double[,], double[]) CalculateTransformations(double[,] tgi, double[,] tG_i, double[,] F_3D)
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

            OrthogonaliseBasisMembranePart(tgi);
            OrthogonaliseBasisMembranePart(tG_i);

            var Qij = CalculateRotationMatrix(tG_i, eye);
            var Qij1 = CalculateRotationMatrix(tgi, eye);

            double[,] F_rve = Transform_F3D_to_Frve(F_3D, Qij, Qij1);


            bool runExample = false;
            if (runExample)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };

                var exampleG1 = Vector.CreateFromArray(new double[] { 14.285714286775336, 3.2374954496774966E-16, 0 });
                var exampleG2 = Vector.CreateFromArray(new double[] { 4.4431960281171579E-17, 14.285714286775336, 0 });
                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], 0 }, { exampleG_1[1], exampleG_2[1],0 },
            { exampleG_1[2],exampleG_2[2],1} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], 0 }, { exampleG1[1], exampleG2[1],0 },
            { exampleG1[2],exampleG2[2],1} };



                double[,] exampleGlTensr = new double[3, 3] { {3.2357191827259157E-05, 0.5 * (8.3678113701926713E-06), 0  },
                { 0.5 * (8.3678113701926713E-06), 1.4505218359772698E-08, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);

                double[,] exampleSpkTensor = new double[,] { {0.0077689617554168137, 0.0010045557546931828,0    },
                { 0.0010045557546931828, 3.4827029271466839E-06, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            


            double[] GLvec = Transform2DDefGradToGL(F_rve);
            (double[] SPKvec, double[,] ConsCartes )= CalculateSPK(GLvec);

            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            double[,] FPKrve = TransformSPKvecToFPK(F_rve, SPKMat );

            double[,] Cinpk = ExpandCijrs(ConsCartes);

            double[,] Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);


            FPKrve = new double[3, 3] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } };
                  
            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPKrve, Qij, Qij1);

            bool runExample2 = true; //z=-...
            if (runExample2)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };

                var exampleG1 = Vector.CreateFromArray(new double[] { -0.010109380507998852, 1.2063529570299063, 7.1865868866394246E-10 });
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.00819119424556895, -6.8721022557950172E-05, 1.0857567150015821 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 0.99993643119452769, 0.0083795855951000336, 0.0075442769836465062 }); // a3_init dld einai hdh normalised
                
                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], exampleG3[0] }, { exampleG_1[1], exampleG_2[1],exampleG3[1] },
            { exampleG_1[2],exampleG_2[2],exampleG3[2]} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], exampleG3[0] }, { exampleG1[1], exampleG2[1],exampleG3[1] },
            { exampleG1[2],exampleG2[2],exampleG3[2]} };



                double[,] exampleGlTensr = new double[3, 3] { {-7.0829702148917644E-05, 0.5 * (1.9389935310977689E-05), 0  },
                { 0.5 * (1.9389935310977689E-05), 0.00015551997161231441, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);


                double[,] exampleSpkTensor = new double[,] { { -468.55859198619288, 296.64505487001753,0    },
                {296.64505487001753, 7463.2325991276839, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            bool runExample3 = true; //z=0
            if (runExample3)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };

                

                var exampleG1 = Vector.CreateFromArray(new double[] { -0.010125086915473506, 1.2082272041845672, 5.403943148540348E-17 });
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.0082039172057740226, -6.8828411272038153E-05, 1.0874431657141146 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 0.99993643119452769, 0.0083795855951000336, 0.0075442769836465062 });
                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], exampleG3[0] }, { exampleG_1[1], exampleG_2[1],exampleG3[1] },
            { exampleG_1[2],exampleG_2[2],exampleG3[2]} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], exampleG3[0] }, { exampleG1[1], exampleG2[1],exampleG3[1] },
            { exampleG1[2],exampleG2[2],exampleG3[2]} };


                double[,] exampleGlTensr = new double[3, 3] { {2.0737144890037307E-05, 0.5 * (1.1308041227672513E-05), 0  },
                { 0.5 * (1.1308041227672513E-05), 0.00012549972453512748, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);


                 double[,] exampleSpkTensor = new double[,] { { 2365.2521089180445, 171.93028053115697,0    },
                {171.93028053115697, 7000.4556652076062, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            bool runExample4 = true; //z=0
            if (runExample4)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };

               

                var exampleG1 = Vector.CreateFromArray(new double[] { -0.01014079332294816, 1.2101014513392281, -7.1865858058507949E-10});
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.0082166401659790958, - 6.8935799986126134E-05, 1.0891296164266471 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 0.99993643119452769, 0.0083795855951000336, 0.0075442769836465062 }); // a3_init dld einai hdh normalised

                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], exampleG3[0] }, { exampleG_1[1], exampleG_2[1],exampleG3[1] },
            { exampleG_1[2],exampleG_2[2],exampleG3[2]} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], exampleG3[0] }, { exampleG1[1], exampleG2[1],exampleG3[1] },
            { exampleG1[2],exampleG2[2],exampleG3[2]} };




                double[,] exampleGlTensr = new double[3, 3] { {0.00011230399192899226, 0.5 * (3.22614714436734E-06), 0  },
                { 0.5 * (3.22614714436734E-06), 9.5479477457940546E-05, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);

               
                double[,] exampleSpkTensor = new double[,] { { 5164.0442904259, 48.748518031151249,0    },
                {48.748518031151249, 6543.1821336104158, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            var Qpi = CalculateRotationMatrix(eye, tgi);
            var Qqj = CalculateRotationMatrix(eye, tG_i);
            var Qrk = Qpi;
            var Qsl = Qqj;

            var Aijkl_3D = Transform_Aijkl_rve_to_Aijkl_3D(Aijkl_rve, Qpi, Qqj, Qrk, Qsl);

            double[] FPK_3D_vec = new double[9] { FPK_3D[0, 0], FPK_3D[1, 1], FPK_3D[2, 2], FPK_3D[0, 1], FPK_3D[1, 2], FPK_3D[2, 0], FPK_3D[0, 2], FPK_3D[1, 0], FPK_3D[2, 1]  };

            return (Aijkl_3D, FPK_3D_vec);

        }

        private void OrthogonaliseBasisMembranePart(double[,] tgi)
        {
            double[] E1 = new double[] { tgi[0, 0], tgi[1, 0], tgi[2, 0] };
            double[] E2 = new double[] { tgi[0, 1], tgi[1, 1], tgi[2, 1] };

            double E1_dot_E2 = 0;
            for (int i1 = 0; i1 < 3; i1++) { E1_dot_E2 += E1[i1] * E2[i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = E1_dot_E2 * E1[i1]; }

            double[] E2ortho = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { E2ortho[i1] = E2[i1] - projection[i1]; }
            double norm1 = (Vector.CreateFromArray(E2ortho)).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { E2[i1] = E2ortho[i1] / norm1; }

            tgi[0, 1] = E2[0];
            tgi[1, 1] = E2[1];
            tgi[2, 1] = E2[2];

            //
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
