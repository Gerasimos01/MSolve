using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    public class tensorOrder2
    {
        public double[,] coefficients { get; set; }

        public double[,] basis1 { get; set; }

        public double[,] basis2 { get; set; }

        public bool isCartesian3Dbasis { get; set; } = false;

        public tensorOrder2() { }

        public tensorOrder2(double[,] coefficients, double[,] basis1, double[,] basis2)
        {
            this.coefficients =  CopyBasis(coefficients);
            this.basis1 =  CopyBasis(basis1);
            this.basis2 = CopyBasis(basis2);
        }

        public tensorOrder2(double[,] coefficients)
        {
            this.coefficients = coefficients;
            basis1 = new double[3, 3]; basis1[0, 0] = 1; basis1[1, 1] = 1; basis1[2, 2] = 1;
            basis2 = CopyBasis(basis1);
            isCartesian3Dbasis = true;
        }

        public void ReplaceBasisWithVector(Vector e1, Vector e2, Vector e3, bool isTheFirstBasisToBeReplaced )
        {
            if(isTheFirstBasisToBeReplaced)
            {
                basis1 = new double[3, 3]
                {
                    {e1[0],e2[0],e3[0] },
                    {e1[1],e2[1],e3[1] },
                    {e1[2],e2[2],e3[2] }
                };
            }
            else
            {
                basis2 = new double[3, 3]
                {
                    {e1[0],e2[0],e3[0] },
                    {e1[1],e2[1],e3[1] },
                    {e1[2],e2[2],e3[2] }
                };
            }
        }

        public void ReplaceBasisWithVector(double[] e1, double[] e2, double[] e3, bool isTheFirstBasisToBeReplaced)
        {
            if (isTheFirstBasisToBeReplaced)
            {
                basis1 = new double[3, 3]
                {
                    {e1[0],e2[0],e3[0] },
                    {e1[1],e2[1],e3[1] },
                    {e1[2],e2[2],e3[2] }
                };
            }
            else
            {
                basis2 = new double[3, 3]
                {
                    {e1[0],e2[0],e3[0] },
                    {e1[1],e2[1],e3[1] },
                    {e1[2],e2[2],e3[2] }
                };
            }
        }

        public double doubleContractCopy1(tensorOrder2 otherTensor)
        {
            double result = 0;

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    var dot1 = basis1[0, i1] * otherTensor.basis1[0, i1] + basis1[1, i1] * otherTensor.basis1[1, i1] + basis1[2, i1] * otherTensor.basis1[2, i1];
                    var dot2 = basis2[0, i2] * otherTensor.basis2[0, i2] + basis2[1, i2] * otherTensor.basis2[1, i2] + basis2[2, i2] * otherTensor.basis2[2, i2];
                    result += dot1 * dot2 * coefficients[i1, i2] * otherTensor.coefficients[i1, i2];
                }
            }


            return result;
        }

        public double doubleContract(tensorOrder2 otherTensor)
        {
            //Dot products

            var EiEk = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    EiEk[i1, i2] = basis1[0, i1] * otherTensor.basis1[0, i2] + basis1[1, i1] * otherTensor.basis1[1, i2] + basis1[2, i1] * otherTensor.basis1[2, i2];
                }
            }

            var EjEl = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    EjEl[i1, i2] = basis2[0, i1] * otherTensor.basis2[0, i2] + basis2[1, i1] * otherTensor.basis2[1, i2] + basis2[2, i1] * otherTensor.basis2[2, i2];
                }
            }


            double result = 0;

            //for (int i1 = 0; i1 < 3; i1++)
            //{
            //    for (int i2 = 0; i2 < 3; i2++)
            //    {
            //        var dot1 = basis1[0, i1] * otherTensor.basis1[0, i1] + basis1[1, i1] * otherTensor.basis1[1, i1] + basis1[2, i1] * otherTensor.basis1[2, i1];
            //        var dot2 = basis2[0, i2] * otherTensor.basis2[0, i2] + basis2[1, i2] * otherTensor.basis2[1, i2] + basis2[2, i2] * otherTensor.basis2[2, i2];
            //        result += dot1 * dot2 * coefficients[i1, i2] * otherTensor.coefficients[i1, i2];
            //    }
            //}

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            result+= EiEk[i,k] * EjEl[j,l]* coefficients[i,j] * otherTensor.coefficients[k,l];
                        }
                    }
                }
            }


            return result;
        }

        public tensorOrder2 TransformTensor(double[,] basis)
        {
            throw new NotImplementedException();
        }

        public tensorOrder2 SingleContract( tensorOrder2 otherTensor)
        {
            

            var EjFk = new double[3, 3]; // contracted basis dot products 

            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    // eswteriko ginomeno
                    for (int i1= 0; i1 < 3; i1++)
                    {
                        EjFk[j, k] += basis2[i1, j] * otherTensor.basis1[i1, k];
                    }
                }
            }


            var newCoeffs = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int l = 0; l < 3; l++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            newCoeffs[i, l] += coefficients[i, j] * otherTensor.coefficients[k, l] * EjFk[j, k];
                        }
                    }
                }
            }

            return new tensorOrder2()
            {
                basis1 = new double[,] { { basis1[0, 0], basis1[0, 1], basis1[0, 2] }, { basis1[1, 0], basis1[1, 1], basis1[1, 2] }, { basis1[2, 0], basis1[2, 1], basis1[2, 2] } },
                basis2 = new double[,] { { otherTensor.basis2[0, 0], otherTensor.basis2[0, 1], otherTensor.basis2[0, 2] }, { otherTensor.basis2[1, 0], otherTensor.basis2[1, 1], otherTensor.basis2[1, 2] }, { otherTensor.basis2[2, 0], otherTensor.basis2[2, 1], otherTensor.basis2[2, 2] } },
                coefficients = newCoeffs,

            };
        }

        public tensorOrder2 Transpose()
        {
            var newBasis1 = new double[,] { { basis2[0, 0], basis2[0,1], basis2[0, 2] }, { basis2[1, 0], basis2[1, 1], basis2[1, 2] }, { basis2[2, 0], basis2[2, 1], basis2[2, 2] } };
            var newBasis2 = new double[,] { { basis1[0,0], basis1[0, 1], basis1[0, 2] }, { basis1[1, 0], basis1[1, 1], basis1[1, 2] }, { basis1[2, 0], basis1[2, 1], basis1[2, 2] } };

            var newCoeffs = new double[,] { { coefficients[0, 0], coefficients[1,0], coefficients[2,0] }, { coefficients[0,1], coefficients[1, 1], coefficients[2,1] }, { coefficients[0,2], coefficients[1,2], coefficients[2, 2] } };

            return new tensorOrder2()
            {
                basis1 = newBasis1,
                basis2 = newBasis2,
                coefficients = newCoeffs
            };
        }

        public tensorOrder2 ProjectIn3DCartesianBasis()
        {
            var eye2 = new double[3, 3]; eye2[0, 0] = 1; eye2[1, 1] = 1; eye2[2, 2] = 1;
            var  eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;

            double[,] B1kEi = new double[3, 3];
            for (int k = 0; k < 3; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    //eswteriko ginomeno
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        B1kEi[k, i] += basis1[i1, k] * eye[i1, i];
                    }
                }
            }

            double[,] B2lEj = new double[3, 3];
            for (int l = 0; l < 3; l++)
            {
                for (int j = 0; j < 3; j++)
                {
                    //eswteriko ginomeno
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        B2lEj[l, j] += basis2[i1, l] * eye[i1, j];
                    }
                }
            }

            double[,] newCoeffs = new double[3,3];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    // provolh olwn twn kl sunistwsws tou parontos tanusth akl_
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            newCoeffs[i, j] += coefficients[k, l] * B1kEi[k, i] * B2lEj[l, j];
                        }
                    }
                }
            }

            return new tensorOrder2()
            {
                basis1 = eye,
                basis2 = eye2,
                coefficients = newCoeffs,
                isCartesian3Dbasis = true
            };

        }

        public tensorOrder2 Scale(double scalar)
        {
            double[,] newCoeffs1 = new double[3, 3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    newCoeffs1[i1, i2] = scalar * coefficients[i1, i2];
                }
            }

            var newTensor = new tensorOrder2()
            {
                basis1 = CopyBasis(basis1),
                basis2 = CopyBasis(basis2),
                coefficients = newCoeffs1
            };

            if(isCartesian3Dbasis)
            { newTensor.isCartesian3Dbasis = true; }

            return newTensor;
        }

        public static double[,] CopyBasis(double[,] basis)
        {
            var copiedBasis = new double[3, 3];
            for (int i1= 0; i1 < 3; i1++)
            {
                for (int i2= 0; i2 < 3; i2++)
                {
                    copiedBasis[i1, i2] = basis[i1, i2];
                }
            }

            return copiedBasis;

        }

        public tensorOrder2 AddTensor(tensorOrder2 otherTensor)
        {
            if(isCartesian3Dbasis&&otherTensor.isCartesian3Dbasis)
            {
                double[,] newCoeffs = new double[3, 3];

                for (int i1 = 0; i1 < 3; i1++)
                {
                    for (int i2 = 0; i2 < 3; i2++)
                    {
                        newCoeffs[i1, i2] = coefficients[i1, i2] + otherTensor.coefficients[i1, i2];
                    }
                }

                return new tensorOrder2()
                {
                    basis1 = CopyBasis(basis1),
                    basis2 = CopyBasis(basis2),
                    coefficients = newCoeffs
                };
            }
            else
            {
                throw new NotImplementedException();
            }
                
        }
    }
}

