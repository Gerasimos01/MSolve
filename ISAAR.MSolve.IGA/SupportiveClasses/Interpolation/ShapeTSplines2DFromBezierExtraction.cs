using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interpolation
{
    /// <summary>
	/// Two-dimensional T-spline shape functions from Bezier extraction.
	/// </summary>
	public class ShapeTSplines2DFromBezierExtraction:IShapeFunction2D
	{
        private readonly int degreeKsi;
        private readonly int degreeHeta;
        private readonly Matrix extractionOperator;
        private readonly ControlPoint[] controlPoints;

        /// <summary>
		/// Two-dimensional T-spline shape functions from Bezier extraction.
		/// </summary>
		public ShapeTSplines2DFromBezierExtraction(int degreeKsi, int degreeHeta, Matrix extractionOperator, ControlPoint[] controlPoints)
		{
            this.degreeKsi = degreeKsi;
            this.degreeHeta = degreeHeta;
            this.extractionOperator = extractionOperator;
            this.controlPoints = controlPoints;

            var gauss = new GaussQuadrature();
			var gaussPoints = gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta,
				new List<Knot>
				{
					new Knot(){ID=0,Ksi=-1,Heta = -1,Zeta = 0},
					new Knot(){ID=1,Ksi=-1,Heta = 1,Zeta = 0},
					new Knot(){ID=2,Ksi=1,Heta = -1,Zeta = 0},
					new Knot(){ID=3,Ksi=1,Heta = 1,Zeta = 0}
				});

			var parametricGaussPointKsi = new double[(degreeKsi + 1) * 2];
            for (var i = 0; i < degreeKsi + 1; i++)
			{
				parametricGaussPointKsi[i] = gaussPoints[i * (degreeHeta + 1)].Ksi;
			}

			var parametricGaussPointHeta = new double[(degreeHeta + 1) * 2];
            for (var i = 0; i < degreeHeta + 1; i++)
			{
				parametricGaussPointHeta[i] = gaussPoints[i].Heta;
			}

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) =
                CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            Values = values;
            DerivativeValuesKsi = derivativeValuesKsi;
            DerivativeValuesHeta = derivativeValuesHeta;
            SecondDerivativeValuesKsi = secondDerivativeValuesKsi;
            SecondDerivativeValuesHeta = secondDerivativeValuesHeta;
            SecondDerivativeValuesKsiHeta = secondDerivativeValuesKsiHeta;
        }

        private (double[,] values, double[,] derivativeValuesKsi, double[,] derivativeValuesHeta, double[,]
            secondDerivativeValuesKsi, double[,] secondDerivativeValuesHeta, double[,] secondDerivativeValuesKsiHeta)
            CalculateTSplines2D(double[] parametricGaussPointKsi, double[] parametricGaussPointHeta)
        {
            var knotValueVectorKsi = new double[(degreeKsi + 1) * 2];
            var knotValueVectorHeta = new double[(degreeHeta + 1) * 2];
            for (var i = 0; i < degreeKsi + 1; i++)
            {
                knotValueVectorKsi[i] = -1;
                knotValueVectorKsi[degreeKsi + 1 + i] = 1;
            }

            for (var i = 0; i < degreeHeta + 1; i++)
            {
                knotValueVectorHeta[i] = -1;
                knotValueVectorHeta[degreeHeta + 1 + i] = 1;
            }

            var bernsteinKsi = new BSplines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
            var bernsteinHeta = new BSplines1D(degreeHeta, knotValueVectorHeta, parametricGaussPointHeta);

            var supportKsi = degreeKsi + 1;
            var supportHeta = degreeHeta + 1;

            var bKsi = MatrixPart(supportKsi, bernsteinKsi.Values);
            var bdKsi = MatrixPart(supportKsi, bernsteinKsi.DerivativeValues);
            var bddKsi = MatrixPart(supportKsi, bernsteinKsi.SecondDerivativeValues);

            var bheta = MatrixPart(supportHeta, bernsteinHeta.Values);
            var bdheta = MatrixPart(supportHeta, bernsteinHeta.DerivativeValues);
            var bddheta = MatrixPart(supportHeta, bernsteinHeta.SecondDerivativeValues);

            var bernsteinShapeFunctions = KroneckerProduct(bKsi, bheta);
            var bernsteinShapeFunctionDerivativesKsi = KroneckerProduct(bheta, bdKsi);
            var bernsteinShapeFunctionDerivativesHeta = KroneckerProduct(bdheta, bKsi);
            var bernsteinShapeFunctionSecondDerivativesKsi = KroneckerProduct(bheta, bddKsi);
            var bernsteinShapeFunctionSecondDerivativesHeta = KroneckerProduct(bddheta, bKsi);
            var bernsteinShapeFunctionSecondDerivativesKsiHeta = KroneckerProduct(bdheta, bdKsi);

            var rationalTSplines = extractionOperator * bernsteinShapeFunctions;
            var rationalTSplineDerivativesKsi = extractionOperator * bernsteinShapeFunctionDerivativesKsi;
            var rationalTSplineDerivativesHeta = extractionOperator * bernsteinShapeFunctionDerivativesHeta;
            var rationalTSplineSecondDerivativesKsi = extractionOperator * bernsteinShapeFunctionSecondDerivativesKsi;
            var rationalTSplineSecondDerivativesHeta = extractionOperator * bernsteinShapeFunctionSecondDerivativesHeta;
            var rationalTSplineSecondDerivativesKsiHeta = extractionOperator * bernsteinShapeFunctionSecondDerivativesKsiHeta;

            var values = new double[controlPoints.Length, supportKsi * supportHeta];
            var derivativeValuesKsi = new double[controlPoints.Length, supportKsi * supportHeta];
            var derivativeValuesHeta = new double[controlPoints.Length, supportKsi * supportHeta];
            var secondDerivativeValuesKsi = new double[controlPoints.Length, supportKsi * supportHeta];
            var secondDerivativeValuesHeta = new double[controlPoints.Length, supportKsi * supportHeta];
            var secondDerivativeValuesKsiHeta = new double[controlPoints.Length, supportKsi * supportHeta];

            for (var i = 0; i < supportKsi; i++)
            {
                for (var j = 0; j < supportHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;
                    double sumdKsidKsi = 0;
                    double sumdHetadHeta = 0;
                    double sumdKsidHeta = 0;

                    var index = i * supportHeta + j;

                    for (var k = 0; k < controlPoints.Length; k++)
                    {
                        sumKsiHeta += rationalTSplines[k, index] * controlPoints[k].WeightFactor;
                        sumdKsiHeta += rationalTSplineDerivativesKsi[k, index] * controlPoints[k].WeightFactor;
                        sumKsidHeta += rationalTSplineDerivativesHeta[k, index] * controlPoints[k].WeightFactor;
                        sumdKsidKsi += rationalTSplineSecondDerivativesKsi[k, index] * controlPoints[k].WeightFactor;
                        sumdHetadHeta += rationalTSplineSecondDerivativesHeta[k, index] * controlPoints[k].WeightFactor;
                        sumdKsidHeta += rationalTSplineSecondDerivativesKsiHeta[k, index] * controlPoints[k].WeightFactor;
                    }

                    for (var k = 0; k < controlPoints.Length; k++)
                    {
                        values[k, index] = rationalTSplines[k, index] * controlPoints[k].WeightFactor / sumKsiHeta;
                        derivativeValuesKsi[k, index] = (rationalTSplineDerivativesKsi[k, index] * sumKsiHeta -
                                                         rationalTSplines[k, index] * sumdKsiHeta) /
                            Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
                        derivativeValuesHeta[k, index] = (rationalTSplineDerivativesHeta[k, index] * sumKsiHeta -
                                                          rationalTSplines[k, index] * sumKsidHeta) /
                            Math.Pow(sumKsiHeta, 2) * controlPoints[k].WeightFactor;
                        secondDerivativeValuesKsi[k, index] = (rationalTSplineSecondDerivativesKsi[k, index] / sumKsiHeta -
                                                               2 * rationalTSplineDerivativesKsi[k, index] * sumdKsiHeta /
                                                               Math.Pow(sumKsiHeta, 2) -
                                                               rationalTSplines[k, index] * sumdKsidKsi /
                                                               Math.Pow(sumKsiHeta, 2) +
                                                               2 * rationalTSplines[k, index] * Math.Pow(sumdKsiHeta, 2) /
                                                               Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
                        secondDerivativeValuesHeta[k, index] = (rationalTSplineSecondDerivativesHeta[k, index] / sumKsiHeta -
                                                                2 * rationalTSplineDerivativesHeta[k, index] * sumKsidHeta /
                                                                Math.Pow(sumKsiHeta, 2) -
                                                                rationalTSplines[k, index] * sumdHetadHeta /
                                                                Math.Pow(sumKsiHeta, 2) +
                                                                2 * rationalTSplines[k, index] * Math.Pow(sumKsidHeta, 2) /
                                                                Math.Pow(sumKsiHeta, 3)) * controlPoints[k].WeightFactor;
                        secondDerivativeValuesKsiHeta[k, index] =
                            (rationalTSplineSecondDerivativesKsiHeta[k, index] / sumKsiHeta -
                             rationalTSplineDerivativesKsi[k, index] * sumKsidHeta /
                             Math.Pow(sumKsiHeta, 2) -
                             rationalTSplineDerivativesHeta[k, index] * sumdKsiHeta /
                             Math.Pow(sumKsiHeta, 2) -
                             rationalTSplines[k, index] * sumdKsidHeta /
                             Math.Pow(sumKsiHeta, 2) +
                             2 * rationalTSplines[k, index] * sumdKsiHeta * sumKsidHeta /
                             Math.Pow(sumKsiHeta, 3)) *
                            controlPoints[k].WeightFactor;
                    }
                }
            }

            return (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta);
        }

        /// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function derivatives per axis Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function second derivatives per axis Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape function mixed second derivatives per axis Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiHeta { get; private set; }

        public double[,] EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Eta).Distinct().ToArray();

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return values;
        }

        public (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Eta).Distinct().ToArray();

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return (derivativeValuesKsi, derivativeValuesHeta);
        }

        public (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta) 
            EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Eta).Distinct().ToArray();

            var (_, _, _, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return (secondDerivativeValuesKsi, secondDerivativeValuesHeta,secondDerivativeValuesKsiHeta);
        }

        /// <summary>
		/// <see cref="Matrix"/> containing T-Spline shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }

		private static Matrix KroneckerProduct(Matrix A, Matrix B)
		{
			var C = Matrix.CreateZero(A.NumRows * B.NumRows, A.NumColumns * B.NumColumns);
			for (var rowAIndex = 0; rowAIndex < A.NumRows; rowAIndex++)
			{
				for (var rowBIndex = 0; rowBIndex < B.NumRows; rowBIndex++)
				{
					for (var columnAIndex = 0; columnAIndex < A.NumColumns; columnAIndex++)
					{
						for (var columnBIndex = 0; columnBIndex < B.NumColumns; columnBIndex++)
						{
							C[rowAIndex * B.NumRows + rowBIndex, columnAIndex * B.NumColumns + columnBIndex] =
								A[rowAIndex, columnAIndex] * B[rowBIndex, columnBIndex];
						}
					}
				}
			}

			return C;
		}

		private static Matrix MatrixPart(int support, double[,] matrix)
		{
			var A = Matrix.CreateZero(support, support);
			for (var i = 0; i < support; i++)
			{
				for (var j = 0; j < support; j++)
				{
					A[i, j] = matrix[i, j];
				}
			}

			return A;
		}

		private static Matrix MatrixPart(int support1, int support2, double[,] matrix)
		{
			var A = Matrix.CreateZero(support1, support2);
			for (var i = 0; i < support1; i++)
			{
				for (var j = 0; j < support2; j++)
				{
					A[i, j] = matrix[i, j];
				}
			}

			return A;
		}

        public double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint)
        {
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return values;
        }

        public (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint)
        {
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return (derivativeValuesKsi,derivativeValuesHeta);
        }

        public (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta)
            EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint)
        {
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};

            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta) = CalculateTSplines2D(parametricGaussPointKsi, parametricGaussPointHeta);

            return (secondDerivativeValuesKsi, secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta);
        }
    }
}
