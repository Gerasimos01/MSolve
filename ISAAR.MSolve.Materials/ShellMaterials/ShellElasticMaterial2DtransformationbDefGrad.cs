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
    public class ShellElasticMaterial2DtransformationbDefGrad : IContinuumMaterial3DDefGrad
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

        public IContinuumMaterial3DDefGrad Clone()
        {
            return new ShellElasticMaterial2DtransformationbDefGrad()
            {
                YoungModulus = this.YoungModulus,
                PoissonRatio = this.PoissonRatio,
            };
        }

        //public void UpdateMaterial(double[] cartesianStrains) //TODO: rename cartesian strains to strains 
        //{
        //    if (CartesianConstitutiveMatrix == null)
        //    {
        //        this.CalculateConstitutiveMatrix(Vector.CreateFromArray(TangentVectorV1), Vector.CreateFromArray(TangentVectorV2));
        //    }

        //    for (int l = 0; l < 3; l++)
        //    {
        //        CartesianStresses[l] = 0;
        //        for (int m = 0; m < 3; m++)
        //        {
        //            CartesianStresses[l] += CartesianConstitutiveMatrix[l, m] * cartesianStrains[m];
        //        }
        //    }
        //}

        public void UpdateState(double[,] F_rve) //TODO: rename cartesian strains to strains 
        {
            double[] GLvec = Transform2DDefGradToGL(F_rve);
            (double[] SPKvec, double[,] ConsCartes) = CalculateSPK(GLvec);
            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            FPKrve = TransformSPKvecToFPK(F_rve, SPKMat);

            double[,] Cinpk = ExpandCijrs(ConsCartes);
            Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);
                                    
        }

        public void UpdateMaterial(double[] F_rve_vec) //TODO: rename cartesian strains to strains 
        {
            var F_rve = new double[2, 2] { { F_rve_vec[0], F_rve_vec[2] }, { F_rve_vec[3], F_rve_vec[1] } };
            double[] GLvec = Transform2DDefGradToGL(F_rve);
            (double[] SPKvec, double[,] ConsCartes) = CalculateSPK(GLvec);
            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            FPKrve = TransformSPKvecToFPK(F_rve, SPKMat);

            double[,] Cinpk = ExpandCijrs(ConsCartes);
            Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);

        }





        public double[,] FPKrve { get; set; } = new double[2, 2] { { 0, 0, }, { 0, 0 } };
        public double[,] Aijkl_rve { get; set; }

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

            double[,] F_rve = Transform_F3D_to_Frve(F_3D, Qij, Qij); // 1);


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



                //double[,] exampleGlTensr = new double[3, 3] { {3.2357191827259157E-05, 0.5 * (8.3678113701926713E-06), 0  },
                //{ 0.5 * (8.3678113701926713E-06), 1.4505218359772698E-08, 0 },  { 0,0,0} };

                //efelkusmos
                double[,] exampleGlTensr = new double[3, 3] { {3.9474348824342087E-05, 0.5 * (-8.15207639374829E-17), 0  },
                { 0.5 * (-8.15207639374829E-17), 0, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);

                double[,] exampleSpkTensor = new double[,] { {0.0094777911499087427, -9.7868119761061161E-15,0    },
                { -9.7868119761061161E-15, 5.044589476002491E-31, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }




            double[] GLvec = Transform2DDefGradToGL(F_rve);
            (double[] SPKvec, double[,] ConsCartes) = CalculateSPK(GLvec);

            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            double[,] FPKrve = TransformSPKvecToFPK(F_rve, SPKMat);

            double[,] Cinpk = ExpandCijrs(ConsCartes);

            double[,] Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);


            FPKrve = new double[3, 3] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } };

            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPKrve, Qij, Qij);// 1);

            bool runExample2 = false; //z=-...
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
                { 0.5 * (1.9389935310977689E-05), 0.00015551997161231441, 0 },  { 0,0,0} }; // afta ta values exoun allaxei logw allaghs proshmou sto class shell ele.

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);


                double[,] exampleSpkTensor = new double[,] { { -468.55859198619288, 296.64505487001753,0    },
                {296.64505487001753, 7463.2325991276839, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            bool runExample3 = false; //z=0
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

            bool runExample4 = false; //z=0
            if (runExample4)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };



                var exampleG1 = Vector.CreateFromArray(new double[] { -0.01014079332294816, 1.2101014513392281, -7.1865858058507949E-10 });
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.0082166401659790958, -6.8935799986126134E-05, 1.0891296164266471 });
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

            bool runExample5 = false; //z= -; // plate warped
            if (runExample5)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };



                var exampleG1 = Vector.CreateFromArray(new double[] { 3.6524987999804055, -0.0030524745407180458, 0 });
                var exampleG2 = Vector.CreateFromArray(new double[] { 0.31047049677853678, 13.638616788095289, 0 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 0, 00, 1 }); // a3_init dld einai hdh normalised

                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], exampleG3[0] }, { exampleG_1[1], exampleG_2[1],exampleG3[1] },
            { exampleG_1[2],exampleG_2[2],exampleG3[2]} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], exampleG3[0] }, { exampleG1[1], exampleG2[1],exampleG3[1] },
            { exampleG1[2],exampleG2[2],exampleG3[2]} };




                double[,] exampleGlTensr = new double[3, 3] { {3.513815970634937E-06, 0.5 * (1.4221485205023754E-06), 0  },
                { 0.5 * (1.4221485205023754E-06), 2.6178867074122536E-08, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);


                double[,] exampleSpkTensor = new double[,] { { 0.19715273909078354, 0.001707543307020885,0    },
                {0.001707543307020885, - 1.9278719468371267E-05, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            bool runExample6 = false; //z= -; // hemispherical
            if (runExample6)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };



                var exampleG1 = Vector.CreateFromArray(new double[] { -0.010109380507998852, 1.2063529570299063, 7.1865868866394246E-10 });
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.00819119424556895, -6.8721022557950172E-05, 1.0857567150015821 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 1.7157003609709955, 0.014377772008092392, 0.012944541613155218 }); // a3_init dld einai hdh normalised

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
                { 0.5 * (3.22614714436734E-06),9.5479477457940546E-05, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);






                double[,] exampleSpkTensor = new double[,] { { 5228.5269163275243, 49.357152115028043,0    },
                {49.357152115028043, 6624.8696580176875, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            bool runExample7 = false; //z=0; // hemispherical
            if (runExample7)
            {
                //
                //double[,] exampleBasis = new double[3, 3] {{ 14.285714922100491, -1.7441909127627026E-08, 0 }, { 8.3652124878453923E-07, 14.285714288340733,0},
                //{ 0,0,1} };



                var exampleG1 = Vector.CreateFromArray(new double[] { -0.010125086915473506, 1.2082272041845672, 5.403943148540348E-17 });
                var exampleG2 = Vector.CreateFromArray(new double[] { -0.0082039172057740226, -6.8828411272038153E-05, 1.0874431657141146 });
                var exampleG3 = Vector.CreateFromArray(new double[] { 1.7263862358694682, 0.014467320903979015, 0.013025163938266793 }); // a3_init dld einai hdh normalised

                double exampleG1_norm_sqred = exampleG1.DotProduct(exampleG1);
                double exampleG2_norm_sqred = exampleG2.DotProduct(exampleG2);
                //double G3_norm_sqred = a3.DotProduct(a3);
                double[] exampleG_1 = new double[3] { exampleG1[0] / exampleG1_norm_sqred, exampleG1[1] / exampleG1_norm_sqred, exampleG1[2] / exampleG1_norm_sqred };
                double[] exampleG_2 = new double[3] { exampleG2[0] / exampleG2_norm_sqred, exampleG2[1] / exampleG2_norm_sqred, exampleG2[2] / exampleG2_norm_sqred };
                double[,] exampleBasis = new double[3, 3] {{ exampleG_1[0], exampleG_2[0], exampleG3[0] }, { exampleG_1[1], exampleG_2[1],exampleG3[1] },
            { exampleG_1[2],exampleG_2[2],exampleG3[2]} };
                double[,] exampleCovarBasis = new double[3, 3] {{ exampleG1[0], exampleG2[0], exampleG3[0] }, { exampleG1[1], exampleG2[1],exampleG3[1] },
            { exampleG1[2],exampleG2[2],exampleG3[2]} };

                double[,] exampleGlTensr = new double[3, 3] { { 2.0737144890037307E-05, 0.5 * (1.1308041227672513E-05), 0  },
                { 0.5 * (1.1308041227672513E-05),0.00012549972453512748, 0 },  { 0,0,0} };

                var exampleQij = CalculateRotationMatrix(tG_i, exampleBasis);
                var exampleQij1 = CalculateRotationMatrix(tgi, exampleBasis);
                double[,] exampleTensorTransformed = Transform_F3D_to_Frve(exampleGlTensr, exampleQij, exampleQij1);


                double[,] exampleSpkTensor = new double[,] { {  2365.2521089180445, 171.93028053115697,0    },
                {171.93028053115697, 7000.4556652076062, 0 },  { 0,0,0} };
                var exampleQijcov = CalculateRotationMatrix(tG_i, exampleCovarBasis);
                var exampleQij1cov = CalculateRotationMatrix(tgi, exampleCovarBasis);
                double[,] exampleSpkTensorTransformed = Transform_F3D_to_Frve(exampleSpkTensor, exampleQijcov, exampleQij1cov);
            }

            var Qpi = CalculateRotationMatrix(eye, tG_i);
            var Qqj = CalculateRotationMatrix(eye, tG_i);
            var Qrk = Qpi;
            var Qsl = Qqj;

            var Aijkl_3D = Transform_Aijkl_rve_to_Aijkl_3D(Aijkl_rve, Qpi, Qqj, Qrk, Qsl);

            double[] FPK_3D_vec = new double[9] { FPK_3D[0, 0], FPK_3D[1, 1], FPK_3D[2, 2], FPK_3D[0, 1], FPK_3D[1, 2], FPK_3D[2, 0], FPK_3D[0, 2], FPK_3D[1, 0], FPK_3D[2, 1] };

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

            int[,] dlj = new int[2, 2] { { 1, 0 }, { 0, 1 } };

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

                            Aijkl[row, col] += SPKMat[i, k] * dlj[l, j];

                            //eswterika loop einai to athroisma logw n kai p
                            for (int n = 0; n < 2; n++)
                            {
                                rowC = thesi[i, n];
                                for (int p = 0; p < 2; p++)
                                {
                                    colC = thesi[p, k];

                                    Aijkl[row, col] += Cinpk[rowC, colC] * F_rve[j, n] * F_rve[l, p];
                                    //Aijkl[row, col] += Cinpk[rowC, colC] * F_rve[j, n] * F_rve[l, p] + SPKMat[i, k] * dlj[l, j];
                                }
                            }

                        }
                    }
                }
            }

            Aijkl = new double[4, 4]
            {
                {Aijkl[0,0], Aijkl[0,1], Aijkl[0,3], Aijkl[0,2] },
                {Aijkl[1,0], Aijkl[1,1], Aijkl[1,3], Aijkl[1,2]},
                {Aijkl[2,0], Aijkl[2,1], Aijkl[2,3], Aijkl[2,2]},
                {Aijkl[3,0], Aijkl[3,1], Aijkl[3,3], Aijkl[3,2]}
            };

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
                        Qij[i1, i2] += +e_new[k1, i2] * e_old[k1, i1];
                    }
                }
            }

            return Qij;

        }

        public double[] Transform2DDefGradToGL(double[,] F)
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

        public double[,] Transform3DDefGradToGL(double[,] F)
        {
            //double[,] GL = new double[2, 2];

            double[,] FtrF = new double[3, 3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        FtrF[i1, i2] += F[k1, i1] * F[k1, i2];
                    }
                }
            }

            FtrF[0, 0] += -1;
            FtrF[1, 1] += -1;

            FtrF[2, 2] += -1;

            var GL = new double[3, 3] { { 0.5 * FtrF[0, 0], 0.5 * FtrF[0, 1],  0.5 * FtrF[0, 2] }, { 0.5 * FtrF[1, 0], 0.5 * FtrF[1, 1],  0.5 * FtrF[1, 2] },
             { 0.5 * FtrF[2,0], 0.5 * FtrF[2, 1],  0.5 * FtrF[2, 2] }};

            //double[] GLvec = new double[3] { GL[0, 0], GL[1, 1], 2 * GL[1, 0], };

            return GL;

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
                        //FPKmat[i1, i2] += F_rve[i1,k1] * SPKMat[k1, i2];
                        FPKmat[i1, i2] += SPKMat[i1, k1] * F_rve[i2, k1];
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

            Ei[0] = surfaceBasisVector1.CopyToArray();
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
            get { return new double[] { FPKrve[0, 0], FPKrve[1, 1], FPKrve[0, 1], FPKrve[1, 0] }; }
        }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (Aijkl_rve == null) UpdateMaterial(new double[4] { 1, 1, 0, 0 });
                return Matrix.CreateFromArray(Aijkl_rve);
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

        public (double[,], double[], double[,], double[,], double[,], double[,], double[,], double[,], double[,], double[,], double[,], double[,], double[], double[],
            double[], tensorOrder2, tensorOrder2, tensorOrder2, tensorOrder2) CalculateTransformationsV2(Vector g1, Vector g2, Vector g3, Vector G1, Vector G2, Vector G3, double[] G_1, double[] G_2, double[] G_3)
        {
            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;
            double[,] tgi = new double[3, 3] { { g1[0], g2[0], g3[0] }, { g1[1], g2[1], g3[1] }, { g1[2], g2[2], g3[2] } };
            double[,] Gi = new double[3, 3] { { G1[0], G2[0], G3[0] }, { G1[1], G2[1], G3[1] }, { G1[2], G2[2], G3[2] } };
            double[,] G_i = new double[3, 3] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } };


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
                double norm = Gi[0, i1] * Gi[0, i1] + Gi[1, i1] * Gi[1, i1] + Gi[2, i1] * Gi[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    Gi[i2, i1] = Gi[i2, i1] / norm;
                }

            }

            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = G_i[0, i1] * G_i[0, i1] + G_i[1, i1] * G_i[1, i1] + G_i[2, i1] * G_i[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i2, i1] = G_i[i2, i1] / norm;
                }

            }

            OrthogonaliseBasisMembranePart(tgi);
            OrthogonaliseBasisMembranePart(Gi);
            OrthogonaliseBasisMembranePart(G_i);

            var ei = tgi;
            var Ei = Gi;
            var E_i = G_i;


            //[..,0] 
            double[] g1__ei = new double[] { g1[0] * ei[0, 0] + g1[1] * ei[1, 0] + g1[2] * ei[2, 0], g1[0] * ei[0, 1] + g1[1] * ei[1, 1] + g1[2] * ei[2, 1], g1[0] * ei[0, 2] + g1[1] * ei[1, 2] + g1[2] * ei[2, 2] };
            double[] g2__ei = new double[] { g2[0] * ei[0, 0] + g2[1] * ei[1, 0] + g2[2] * ei[2, 0], g2[0] * ei[0, 1] + g2[1] * ei[1, 1] + g2[2] * ei[2, 1], g2[0] * ei[0, 2] + g2[1] * ei[1, 2] + g2[2] * ei[2, 2] };
            double[] g3__ei = new double[] { g3[0] * ei[0, 0] + g3[1] * ei[1, 0] + g3[2] * ei[2, 0], g3[0] * ei[0, 1] + g3[1] * ei[1, 1] + g3[2] * ei[2, 1], g3[0] * ei[0, 2] + g3[1] * ei[1, 2] + g3[2] * ei[2, 2] };

            //double[] g1__Ei = new double[] { g1[0] * Ei[0, 0] + g1[1] * Ei[1, 0] + g1[2] * Ei[2, 0], g1[0] * Ei[0, 1] + g1[1] * Ei[1, 1] + g1[2] * Ei[2, 1], g1[0] * Ei[0, 2] + g1[1] * Ei[1, 2] + g1[2] * Ei[2, 2] };
            //double[] g2__Ei = new double[] { g2[0] * Ei[0, 0] + g2[1] * Ei[1, 0] + g2[2] * Ei[2, 0], g2[0] * Ei[0, 1] + g2[1] * Ei[1, 1] + g2[2] * Ei[2, 1], g2[0] * Ei[0, 2] + g2[1] * Ei[1, 2] + g2[2] * Ei[2, 2] };
            //double[] g3__Ei = new double[] { g3[0] * Ei[0, 0] + g3[1] * Ei[1, 0] + g3[2] * Ei[2, 0], g3[0] * Ei[0, 1] + g3[1] * Ei[1, 1] + g3[2] * Ei[2, 1], g3[0] * Ei[0, 2] + g3[1] * Ei[1, 2] + g3[2] * Ei[2, 2] };

            double[] G_1__Ei = new double[] { G_1[0] * Ei[0, 0] + G_1[1] * Ei[1, 0] + G_1[2] * Ei[2, 0], G_1[0] * Ei[0, 1] + G_1[1] * Ei[1, 1] + G_1[2] * Ei[2, 1], G_1[0] * Ei[0, 2] + G_1[1] * Ei[1, 2] + G_1[2] * Ei[2, 2] };
            double[] G_2__Ei = new double[] { G_2[0] * Ei[0, 0] + G_2[1] * Ei[1, 0] + G_2[2] * Ei[2, 0], G_2[0] * Ei[0, 1] + G_2[1] * Ei[1, 1] + G_2[2] * Ei[2, 1], G_2[0] * Ei[0, 2] + G_2[1] * Ei[1, 2] + G_2[2] * Ei[2, 2] };
            double[] G_3__Ei = new double[] { G_3[0] * Ei[0, 0] + G_3[1] * Ei[1, 0] + G_3[2] * Ei[2, 0], G_3[0] * Ei[0, 1] + G_3[1] * Ei[1, 1] + G_3[2] * Ei[2, 1], G_3[0] * Ei[0, 2] + G_3[1] * Ei[1, 2] + G_3[2] * Ei[2, 2] };

            double[] G1__Ei = new double[] { G1[0] * Ei[0, 0] + G1[1] * Ei[1, 0] + G1[2] * Ei[2, 0], G1[0] * Ei[0, 1] + G1[1] * Ei[1, 1] + G1[2] * Ei[2, 1], G1[0] * Ei[0, 2] + G1[1] * Ei[1, 2] + G1[2] * Ei[2, 2] };
            double[] G2__Ei = new double[] { G2[0] * Ei[0, 0] + G2[1] * Ei[1, 0] + G2[2] * Ei[2, 0], G2[0] * Ei[0, 1] + G2[1] * Ei[1, 1] + G2[2] * Ei[2, 1], G2[0] * Ei[0, 2] + G2[1] * Ei[1, 2] + G2[2] * Ei[2, 2] };
            double[] G3__Ei = new double[] { G3[0] * Ei[0, 0] + G3[1] * Ei[1, 0] + G3[2] * Ei[2, 0], G3[0] * Ei[0, 1] + G3[1] * Ei[1, 1] + G3[2] * Ei[2, 1], G3[0] * Ei[0, 2] + G3[1] * Ei[1, 2] + G3[2] * Ei[2, 2] };

            (double[] G_1__Ei_alte, double[] G_2__Ei_alte, double[] G_3__Ei_alte) = CalculateContravariants(Vector.CreateFromArray(G1__Ei), Vector.CreateFromArray(G2__Ei), Vector.CreateFromArray(G3__Ei));


            double[] G_1__E_i = new double[] { G_1[0] * E_i[0, 0] + G_1[1] * E_i[1, 0] + G_1[2] * E_i[2, 0], G_1[0] * E_i[0, 1] + G_1[1] * E_i[1, 1] + G_1[2] * E_i[2, 1], G_1[0] * E_i[0, 2] + G_1[1] * E_i[1, 2] + G_1[2] * E_i[2, 2] };
            double[] G_2__E_i = new double[] { G_2[0] * E_i[0, 0] + G_2[1] * E_i[1, 0] + G_2[2] * E_i[2, 0], G_2[0] * E_i[0, 1] + G_2[1] * E_i[1, 1] + G_2[2] * E_i[2, 1], G_2[0] * E_i[0, 2] + G_2[1] * E_i[1, 2] + G_2[2] * E_i[2, 2] };
            double[] G_3__E_i = new double[] { G_3[0] * E_i[0, 0] + G_3[1] * E_i[1, 0] + G_3[2] * E_i[2, 0], G_3[0] * E_i[0, 1] + G_3[1] * E_i[1, 1] + G_3[2] * E_i[2, 1], G_3[0] * E_i[0, 2] + G_3[1] * E_i[1, 2] + G_3[2] * E_i[2, 2] };

            double[,] F = CaclculateDefGrad3D(g1__ei, g2__ei, g3__ei, G_1__Ei, G_2__Ei, G_3__Ei);

            double[,] F_case_a = CaclculateDefGrad3D(g1__ei, g2__ei, g3__ei, G_1__E_i, G_2__E_i, G_3__E_i);

            //double[,] F_case_b = CaclculateDefGrad3D(g1__Ei, g2__Ei, g3__Ei, G_1__E_i, G_2__E_i, G_3__E_i);

            double[,] F_rve = new double[,]
            {
                {F[0,0], F[0,1] },
                {F[1,0], F[1,1] }
            };

            double[,] F_rve_case_a = new double[,]
            {
                {F_case_a[0,0], F_case_a[0,1] },
                {F_case_a[1,0], F_case_a[1,1] }
            };

            double[] GLvec = Transform2DDefGradToGL(F_rve);

            var GL_coeffs = new double[3, 3]
                            {
                                {GLvec[0],0.5*GLvec[2],0 },
                                {0.5*GLvec[2],GLvec[1],0 },
                                {0,0,0 },
                            };

            tensorOrder2 GLtensor = new tensorOrder2(GL_coeffs, Ei, Ei);
            var GLtensorProjected = GLtensor.ProjectIn3DCartesianBasis();

            (double[] SPKvec, double[,] ConsCartes) = CalculateSPK(GLvec);

            var SPK_coeffs = new double[3, 3]
                            {
                                {SPKvec[0],SPKvec[2],0 },
                                {SPKvec[2],SPKvec[1],0 },
                                {0,0,0 },
                            };

            tensorOrder2 SPKtensor = new tensorOrder2(SPK_coeffs, Ei, Ei);
            var SPKtensorProjected = SPKtensor.ProjectIn3DCartesianBasis();





            double[,] SPKMat = new double[2, 2] { { SPKvec[0], SPKvec[2] }, { SPKvec[2], SPKvec[1] } };
            double[,] FPKrve = TransformSPKvecToFPK(F_rve, SPKMat);

            double[,] Cinpk = ExpandCijrs(ConsCartes);

            double[,] Aijkl_rve = TransformCinpk(Cinpk, F_rve, SPKMat/*,...*/);

            #region more case a transformations
            double[] GLvec_a = Transform2DDefGradToGL(F_rve_case_a);
            (double[] SPKvec_a, double[,] ConsCartes_a) = CalculateSPK(GLvec_a);

            double[,] SPKMat_a = new double[2, 2] { { SPKvec_a[0], SPKvec_a[2] }, { SPKvec_a[2], SPKvec_a[1] } };
            double[,] FPKrve_a = TransformSPKvecToFPK(F_rve_case_a, SPKMat_a);

            double[,] Cinpk_a = ExpandCijrs(ConsCartes_a);

            double[,] Aijkl_rve_a = TransformCinpk(Cinpk_a, F_rve_case_a, SPKMat_a/*,...*/);
            #endregion

            var cartes_to_Gi = CalculateRotationMatrix(Gi, eye);
            var cartes_to_tgi = CalculateRotationMatrix(tgi, eye);

            FPKrve = new double[3, 3] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } };
            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPKrve, cartes_to_Gi, cartes_to_tgi);// 1);

            double[,] FPK_3D_tr_basis = Transform_FPK_rve_To_FPK_3D(FPKrve, cartes_to_tgi, cartes_to_Gi);// 1); 

            var FPKrve_coeffs = tensorOrder2.CopyBasis(FPKrve);

            tensorOrder2 FPKtensor = new tensorOrder2(FPKrve_coeffs, Ei, ei);
            var FPKtensorProjected = FPKtensor.ProjectIn3DCartesianBasis();



            #region more case a calculations
            var cartes_to_G_i = CalculateRotationMatrix(G_i, eye);
            FPKrve_a = new double[3, 3] { { FPKrve_a[0, 0], FPKrve_a[0, 1], 0 }, { FPKrve_a[1, 0], FPKrve_a[1, 1], 0 }, { 0, 0, 0 } };
            double[,] FPK_3D_a = Transform_FPK_rve_To_FPK_3D(FPKrve_a, cartes_to_Gi, cartes_to_tgi);// 1);

            double[,] FPK_3D_tr_basis_a = Transform_FPK_rve_To_FPK_3D(FPKrve_a, cartes_to_tgi, cartes_to_Gi);// 1); 

            #endregion


            var ch01_FPK_3D = Calculate3DtensorFrom2D(new double[2, 2] { { FPKrve[0, 0], FPKrve[0, 1] }, { FPKrve[1, 0], FPKrve[1, 1] } },
                Vector.CreateFromArray(new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }), Vector.CreateFromArray(new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }),
            Ei);//revisit this maybe we should pass dg de1_dr instead of dg1_dr.
                //var ch02_FPK_3D = Calculate3DtensorFrom2D(new double[2, 2] { { FPKrve[0, 0], FPKrve[0, 1] }, { FPKrve[1, 0], FPKrve[1, 1] } },
                //    Vector.CreateFromArray(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }), Vector.CreateFromArray(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }),
                //ei);

            double[,] ch02_FPK_3D = Calculate3DtensorFrom2Dcorrected_normaliseBothBasesCase(new double[,] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } }, Vector.CreateFromArray(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }), Vector.CreateFromArray(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }),
                Vector.CreateFromArray(new double[] { Ei[0, 2], Ei[1, 2], Ei[2, 2] }), new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }, new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }, new double[] { ei[0, 2], ei[1, 2], ei[2, 2] });
            double[] ch02_FPK_3D_vec = new double[9] { ch02_FPK_3D[0, 0], ch02_FPK_3D[1, 1], ch02_FPK_3D[2, 2], ch02_FPK_3D[0, 1], ch02_FPK_3D[1, 2], ch02_FPK_3D[2, 0], ch02_FPK_3D[0, 2], ch02_FPK_3D[1, 0], ch02_FPK_3D[2, 1] };

            tensorOrder2 ch03_FPK_3D = new tensorOrder2() { basis1 = Ei, basis2 = ei, coefficients = new double[,] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } } };

            //FPK_3D = ch01_FPK_3D;

            var ch01_GL_3D = Calculate3DtensorFrom2D(new double[2, 2] { { GLvec[0], 0.5 * GLvec[2] }, { 0.5 * GLvec[2], GLvec[1] } }, Vector.CreateFromArray(new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }), Vector.CreateFromArray(new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }),
            Ei);
            var ch01_SPKMat_3D = Calculate3DtensorFrom2D(SPKMat, Vector.CreateFromArray(new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }), Vector.CreateFromArray(new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }),
            Ei);

            var GL_exte = new double[3, 3] { { GLvec[0], 0.5 * GLvec[2], 0 }, { 0.5 * GLvec[2], GLvec[1], 0 }, { 0, 0, 0 } };

            var SPK_exte = new double[3, 3] { { SPKvec[0], SPKvec[2], 0 }, { SPKvec[2], SPKvec[1], 0 }, { 0, 0, 0 } };

            var GL3D = Transform_FPK_rve_To_FPK_3D(GL_exte, cartes_to_Gi, cartes_to_tgi);

            var SPKMat3D = Transform_FPK_rve_To_FPK_3D(SPK_exte, cartes_to_Gi, cartes_to_tgi);



            var Qpi = CalculateRotationMatrix(eye, tgi);
            var Qqj = CalculateRotationMatrix(eye, Gi);
            var Qrk = Qpi;
            var Qsl = Qqj;

            var Aijkl_3D = Transform_Aijkl_rve_to_Aijkl_3D(Aijkl_rve, Qpi, Qqj, Qrk, Qsl);
            //var Aijkl_3D = Transform_Aijkl_rve_to_Aijkl_3D(Aijkl_rve, Qqj, Qpi, Qsl, Qrk);
            double[] FPK_3D_a_vec = new double[9] { FPK_3D_a[0, 0], FPK_3D_a[1, 1], FPK_3D_a[2, 2], FPK_3D_a[0, 1], FPK_3D_a[1, 2], FPK_3D_a[2, 0], FPK_3D_a[0, 2], FPK_3D_a[1, 0], FPK_3D_a[2, 1] };

            double[] FPK_3D_vec = new double[9] { FPK_3D[0, 0], FPK_3D[1, 1], FPK_3D[2, 2], FPK_3D[0, 1], FPK_3D[1, 2], FPK_3D[2, 0], FPK_3D[0, 2], FPK_3D[1, 0], FPK_3D[2, 1] };

            double[] FPK_3D_tr_basis_a_vec = new double[9] { FPK_3D_tr_basis_a[0, 0], FPK_3D_tr_basis_a[1, 1], FPK_3D_tr_basis_a[2, 2], FPK_3D_tr_basis_a[0, 1], FPK_3D_tr_basis_a[1, 2], FPK_3D_tr_basis_a[2, 0], FPK_3D_tr_basis_a[0, 2], FPK_3D_tr_basis_a[1, 0], FPK_3D_tr_basis_a[2, 1] };

            return (Aijkl_3D, FPK_3D_vec, FPKrve, Ei, Aijkl_rve, ei, F_rve, GL3D, SPKMat3D, ch01_GL_3D, ch01_SPKMat_3D, FPK_3D_tr_basis, FPK_3D_a_vec, FPK_3D_tr_basis_a_vec,
                ch02_FPK_3D_vec, ch03_FPK_3D, GLtensorProjected, SPKtensorProjected, FPKtensorProjected);
        }

        private double[,] CaclculateDefGrad3D(double[] g1__ei, double[] g2__ei, double[] g3__ei, double[] G_1__Ei, double[] G_2__Ei, double[] G_3__Ei)
        {
            double[,] F = new double[3, 3] {
                { g1__ei[0]*G_1__Ei[0]+g2__ei[0]*G_2__Ei[0]+g3__ei[0]*G_3__Ei[0], g1__ei[0]*G_1__Ei[1]+g2__ei[0]*G_2__Ei[1]+g3__ei[0]*G_3__Ei[1], g1__ei[0]*G_1__Ei[2]+g2__ei[0]*G_2__Ei[2]+g3__ei[0]*G_3__Ei[2] },
                { g1__ei[1]*G_1__Ei[0]+g2__ei[1]*G_2__Ei[0]+g3__ei[1]*G_3__Ei[0], g1__ei[1]*G_1__Ei[1]+g2__ei[1]*G_2__Ei[1]+g3__ei[1]*G_3__Ei[1], g1__ei[1]*G_1__Ei[2]+g2__ei[1]*G_2__Ei[2]+g3__ei[1]*G_3__Ei[2] },
                { g1__ei[2]*G_1__Ei[0]+g2__ei[2]*G_2__Ei[0]+g3__ei[2]*G_3__Ei[0], g1__ei[2]*G_1__Ei[1]+g2__ei[2]*G_2__Ei[1]+g3__ei[2]*G_3__Ei[1], g1__ei[2]*G_1__Ei[2]+g2__ei[2]*G_2__Ei[2]+g3__ei[2]*G_3__Ei[2] },
                        };
            return F;
        }

        private double[,] Calculate3DtensorFrom2D(double[,] FPK_2D, Vector dg1_dr, Vector dg2_dr, double[,] Ei)
        {
            double[,] tensor3D = new double[3, 3];
            double[] leftVec = new double[3];
            double[] rigthVec = new double[3];

            for (int i1 = 0; i1 < 2; i1++)
            {
                if (i1 == 0) { leftVec = dg1_dr.CopyToArray(); }
                else if (i1 == 1) { leftVec = dg2_dr.CopyToArray(); }
                for (int i2 = 0; i2 < 2; i2++)
                {
                    rigthVec = new double[] { Ei[0, i2], Ei[1, i2], Ei[2, i2] };

                    for (int i3 = 0; i3 < 3; i3++)
                    {
                        for (int i4 = 0; i4 < 3; i4++)
                        {
                            tensor3D[i3, i4] = FPK_2D[i1, i2] * leftVec[i3] * rigthVec[i4];
                        }
                    }

                }
            }

            return tensor3D;
        }

        private (double[] G_1, double[] G_2, double[] G_3) CalculateContravariants(Vector g1, Vector g2, Vector a3)
        {
            var auxMatrix1 = Matrix.CreateZero(3, 3);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = g1.DotProduct(g1);
            auxMatrix1[0, 1] = g1.DotProduct(g2);
            auxMatrix1[0, 2] = g1.DotProduct(a3);
            auxMatrix1[1, 0] = g2.DotProduct(g1);
            auxMatrix1[1, 1] = g2.DotProduct(g2);
            auxMatrix1[1, 2] = g2.DotProduct(a3);
            auxMatrix1[2, 0] = a3.DotProduct(g1);
            auxMatrix1[2, 1] = a3.DotProduct(g2);
            auxMatrix1[2, 2] = a3.DotProduct(a3);
            Matrix inverse = auxMatrix1.Invert(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)
                                                  //TODO: auxMatrix1.Invert2x2AndDeterminant(1e-20) for bad geometry

            //Contravariant base vectors
            double[][] G_i = new double[3][];
            for (int i1 = 0; i1 < 3; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * g1[i2] + inverse[i1, 1] * g2[i2] + inverse[i1, 2] * a3[i2];
                }
            }

            return (G_i[0], G_i[1], G_i[2]);
        }

        private double[,] Calculate3DtensorFrom2Dcorrected_normaliseBothBasesCase(double[,] eye3, Vector dg1_dr, Vector dg2_dr, Vector da3_dr, double[] G_1, double[] G_2, double[] G_3)
        {
            //double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;

            #region create and normalise ei
            double[,] ei = new double[3, 3];
            double norm_e1 = dg1_dr[0] * dg1_dr[0] + dg1_dr[1] * dg1_dr[1] + dg1_dr[2] * dg1_dr[2];
            norm_e1 = Math.Sqrt(norm_e1);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 0] = dg1_dr[i2] / norm_e1;
            }

            double norm_e2 = dg2_dr[0] * dg2_dr[0] + dg2_dr[1] * dg2_dr[1] + dg2_dr[2] * dg2_dr[2];
            norm_e2 = Math.Sqrt(norm_e2);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 1] = dg2_dr[i2] / norm_e2;
            }

            double norm_e3 = da3_dr[0] * da3_dr[0] + da3_dr[1] * da3_dr[1] + da3_dr[2] * da3_dr[2];
            norm_e3 = Math.Sqrt(norm_e3);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 2] = da3_dr[i2] / norm_e3;
            }
            #endregion

            #region adapt FPK_2D for normalisation of basis vectors

            double coef1 = 0;
            double coef2 = 0;
            double[,] FPK_2D_in_normalised = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                if (i1 == 0) { coef1 = norm_e1; }
                else if (i1 == 1) { coef1 = norm_e2; }
                else if (i1 == 2) { coef1 = norm_e3; }
                for (int i2 = 0; i2 < 3; i2++)
                {
                    //if (i2 == 0) { coef2 = norm_e1; }
                    //else if (i2 == 1) { coef2 = norm_e2; }

                    //FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1 * coef2;
                    FPK_2D_in_normalised[i1, i2] = eye3[i1, i2] * coef1;// * coef2;
                }
            }
            #endregion

            #region create and normalise Ei
            double[,] Ei = new double[3, 3];
            double norm_E1 = G_1[0] * G_1[0] + G_1[1] * G_1[1] + G_1[2] * G_1[2];
            norm_E1 = Math.Sqrt(norm_E1);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 0] = G_1[i2] / norm_E1;
            }

            double norm_E2 = G_2[0] * G_2[0] + G_2[1] * G_2[1] + G_2[2] * G_2[2];
            norm_E2 = Math.Sqrt(norm_E2);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 1] = G_2[i2] / norm_E2;
            }

            double norm_E3 = G_3[0] * G_3[0] + G_3[1] * G_3[1] + G_3[2] * G_3[2];
            norm_E3 = Math.Sqrt(norm_E3);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 2] = G_3[i2] / norm_E3;
            }
            #endregion

            #region adapt FPK_2D for normalisation of basis vectors

            //double coef1 = 0;
            //double coef2 = 0;
            //double[,] FPK_2D_in_normalised = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {

                for (int i2 = 0; i2 < 3; i2++)
                {
                    if (i2 == 0) { coef2 = norm_E1; }
                    else if (i2 == 1) { coef2 = norm_E2; }
                    else if (i2 == 2) { coef2 = norm_E3; }

                    //FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1 * coef2;
                    FPK_2D_in_normalised[i1, i2] = FPK_2D_in_normalised[i1, i2] * coef2;
                }
            }
            #endregion

            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;
            var cartes_to_Gi = CalculateRotationMatrix(Ei, eye);
            var cartes_to_tgi = CalculateRotationMatrix(ei, eye);

            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPK_2D_in_normalised, cartes_to_Gi, cartes_to_tgi);// 1);



            return FPK_3D;
        }
    }
}
