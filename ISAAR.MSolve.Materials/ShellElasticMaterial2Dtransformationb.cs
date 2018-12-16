using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Materials
{
	 public class ShellElasticMaterial2Dtransformationb : IShellMaterial
	{
		public double[] NormalVectorV3 { get; set; }
		public double[] TangentVectorV1 { get; set; }
		public double[] TangentVectorV2 { get; set; }
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
        Matrix2D transformationMatrix; // gia to shell

        private bool modified; 
		private double[,] CartesianConstitutiveMatrix;
		private double[] CartesianStresses = new double[6];

		object ICloneable.Clone() => Clone();

		public IShellMaterial Clone()
		{
			return new ShellElasticMaterial2Dtransformationb()
			{
				YoungModulus = this.YoungModulus,
				PoissonRatio = this.PoissonRatio,
			};
		}

		public void UpdateMaterial(double[] cartesianStrains) //TODO: rename cartesian strains to strains 
		{
			if (CartesianConstitutiveMatrix == null)
			{
				this.CalculateConstitutiveMatrix(new Vector(TangentVectorV1), new Vector(TangentVectorV2));
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

		private void CalculateConstitutiveMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
		{
            this.CalculateTransformationMatrix(new Vector(TangentVectorV1), new Vector(TangentVectorV2));

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
            var constitutiveMatrix = (transformationMatrix.Transpose() * (new Matrix2D(OriginalConstitutiveMatrix)) * transformationMatrix);
            CartesianConstitutiveMatrix = constitutiveMatrix.Data;
		}

        private void CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = new Matrix2D(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)

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

            Ei[0] = new double[3]; Array.Copy(surfaceBasisVector1.Data, Ei[0], 3);
            double G1_norm = surfaceBasisVector1.Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = new Vector(Ei[1]).Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = new Vector(Ei[i1]).DotProduct(new Vector(G_i[i2]));
                }
            }

            transformationMatrix = new Matrix2D(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
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

		public IMatrix2D ConstitutiveMatrix
		{
			get
			{
				if (CartesianConstitutiveMatrix == null) UpdateMaterial(new double[6]);
				return new Matrix2D(CartesianConstitutiveMatrix);
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
