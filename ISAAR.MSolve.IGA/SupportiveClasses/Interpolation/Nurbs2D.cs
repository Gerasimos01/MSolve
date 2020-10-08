using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interpolation
{
    /// <summary>
	/// Two-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs2D:IShapeFunction2D
	{
        private readonly ControlPoint[] controlPoints;
        private readonly int degreeHeta;
        private readonly int degreeKsi;
        private readonly double[] knotValueVectorHeta;
        private readonly double[] knotValueVectorKsi;
        private int numberOfControlPointsHeta;

        /// <summary>
		/// Defines a 2D NURBS shape function for an element given the per axis gauss point coordinates.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		/// <param name="parametricGaussPointKsi">An <see cref="IVector"/> containing Gauss points of axis Ksi.</param>
		/// <param name="parametricGaussPointHeta">An <see cref="IVector"/> containing Gauss points of axis Heta.</param>
		public Nurbs2D(int degreeKsi, double[] knotValueVectorKsi,
                       int degreeHeta, double[] knotValueVectorHeta,
                       ControlPoint[] controlPoints,GaussLegendrePoint3D[] gaussPoints)
		{
            this.degreeKsi = degreeKsi;
            this.knotValueVectorKsi = knotValueVectorKsi;
            this.degreeHeta = degreeHeta;
            this.knotValueVectorHeta = knotValueVectorHeta;
            this.controlPoints = controlPoints;
            var parametricPointsCount = gaussPoints.Length;
            numberOfControlPointsHeta = knotValueVectorHeta.Length - degreeHeta - 1;
            var parametricGaussPointKsi = gaussPoints.Select(x => x.Ksi).Distinct().ToArray();
            var parametricGaussPointHeta = gaussPoints.Select(x => x.Heta).Distinct().ToArray();
            (Values, DerivativeValuesKsi, DerivativeValuesHeta, SecondDerivativeValuesKsi, SecondDerivativeValuesHeta,
                SecondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
        }

        private (double[,] values, double[,] derivativeValuesKsi, double[,] derivativeValuesHeta,
            double[,] secondDerivativeValuesKsi, double[,] secondDerivativeValuesHeta, double[,] secondDerivativeValuesKsiHeta)
            CalculateNurbs(double[] parametricGaussPointKsi, double[] parametricGaussPointHeta,
            int parametricPointsCount)
        {
            var bSplinesKsi = new BSplines1D(degreeKsi, knotValueVectorKsi, parametricGaussPointKsi);
            var bSplinesHeta = new BSplines1D(degreeHeta, knotValueVectorHeta,
                parametricGaussPointHeta);

            int supportKsi = parametricGaussPointKsi.Length;
            int supportHeta = parametricGaussPointHeta.Length;
            int numberOfElementControlPoints = (degreeKsi + 1) * (degreeHeta + 1);

            var values = new double[numberOfElementControlPoints, parametricPointsCount];
            var derivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
            var derivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
            var secondDerivativeValuesKsi = new double[numberOfElementControlPoints, parametricPointsCount];
            var secondDerivativeValuesHeta = new double[numberOfElementControlPoints, parametricPointsCount];
            var secondDerivativeValuesKsiHeta = new double[numberOfElementControlPoints, parametricPointsCount];

            for (int i = 0; i < supportKsi; i++)
            {
                for (int j = 0; j < supportHeta; j++)
                {
                    double sumKsiHeta = 0;
                    double sumdKsiHeta = 0;
                    double sumKsidHeta = 0;
                    double sumdKsidKsi = 0;
                    double sumdHetadHeta = 0;
                    double sumdKsidHeta = 0;

                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
                        int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;
                        sumKsiHeta += bSplinesKsi.Values[indexKsi, i] *
                                      bSplinesHeta.Values[indexHeta, j] *
                                      controlPoints[k].WeightFactor;
                        sumdKsiHeta += bSplinesKsi.DerivativeValues[indexKsi, i] *
                                       bSplinesHeta.Values[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumKsidHeta += bSplinesKsi.Values[indexKsi, i] *
                                       bSplinesHeta.DerivativeValues[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumdKsidKsi += bSplinesKsi.SecondDerivativeValues[indexKsi, i] *
                                       bSplinesHeta.Values[indexHeta, j] *
                                       controlPoints[k].WeightFactor;
                        sumdHetadHeta += bSplinesKsi.Values[indexKsi, i] *
                                         bSplinesHeta.SecondDerivativeValues[indexHeta, j] *
                                         controlPoints[k].WeightFactor;
                        sumdKsidHeta += bSplinesKsi.DerivativeValues[indexKsi, i] *
                                        bSplinesHeta.DerivativeValues[indexHeta, j] *
                                        controlPoints[k].WeightFactor;
                    }

                    for (int k = 0; k < numberOfElementControlPoints; k++)
                    {
                        int indexKsi = controlPoints[k].ID / numberOfControlPointsHeta;
                        int indexHeta = controlPoints[k].ID % numberOfControlPointsHeta;

                        values[k, i * supportHeta + j] =
                            bSplinesKsi.Values[indexKsi, i] *
                            bSplinesHeta.Values[indexHeta, j] *
                            controlPoints[k].WeightFactor / sumKsiHeta;

                        derivativeValuesKsi[k, i * supportHeta + j] =
                            bSplinesHeta.Values[indexHeta, j] * controlPoints[k].WeightFactor *
                            (bSplinesKsi.DerivativeValues[indexKsi, i] * sumKsiHeta -
                             bSplinesKsi.Values[indexKsi, i] * sumdKsiHeta) / Math.Pow(sumKsiHeta, 2);

                        derivativeValuesHeta[k, i * supportHeta + j] =
                            bSplinesKsi.Values[indexKsi, i] * controlPoints[k].WeightFactor *
                            (bSplinesHeta.DerivativeValues[indexHeta, j] * sumKsiHeta -
                             bSplinesHeta.Values[indexHeta, j] * sumKsidHeta) / Math.Pow(sumKsiHeta, 2);

                        secondDerivativeValuesKsi[k, i * supportHeta + j] =
                            bSplinesHeta.Values[indexHeta, j] * controlPoints[k].WeightFactor *
                            (bSplinesKsi.SecondDerivativeValues[indexKsi, i] / sumKsiHeta -
                             2 * bSplinesKsi.DerivativeValues[indexKsi, i] * sumdKsiHeta /
                             Math.Pow(sumKsiHeta, 2) -
                             bSplinesKsi.Values[indexKsi, i] * sumdKsidKsi / Math.Pow(sumKsiHeta, 2) +
                             2 * bSplinesKsi.Values[indexKsi, i] * Math.Pow(sumdKsiHeta, 2) /
                             Math.Pow(sumKsiHeta, 3));

                        secondDerivativeValuesHeta[k, i * supportHeta + j] =
                            bSplinesKsi.Values[indexKsi, i] * controlPoints[k].WeightFactor *
                            (bSplinesHeta.SecondDerivativeValues[indexHeta, j] / sumKsiHeta -
                             2 * bSplinesHeta.DerivativeValues[indexHeta, j] * sumKsidHeta /
                             Math.Pow(sumKsiHeta, 2) -
                             bSplinesHeta.Values[indexHeta, j] * sumdHetadHeta / Math.Pow(sumKsiHeta, 2) +
                             2 * bSplinesHeta.Values[indexHeta, j] * Math.Pow(sumKsidHeta, 2) /
                             Math.Pow(sumKsiHeta, 3));

                        secondDerivativeValuesKsiHeta[k, i * supportHeta + j] =
                            controlPoints[k].WeightFactor *
                            (bSplinesKsi.DerivativeValues[indexKsi, i] *
                             bSplinesHeta.DerivativeValues[indexHeta, j] / sumKsiHeta -
                             bSplinesKsi.DerivativeValues[indexKsi, i] *
                             bSplinesHeta.Values[indexHeta, j] *
                             sumKsidHeta / Math.Pow(sumKsiHeta, 2) -
                             bSplinesKsi.Values[indexKsi, i] *
                             bSplinesHeta.DerivativeValues[indexHeta, j] *
                             sumdKsiHeta / Math.Pow(sumKsiHeta, 2) -
                             bSplinesKsi.Values[indexKsi, i] * bSplinesHeta.Values[indexHeta, j] *
                             sumdKsidHeta / Math.Pow(sumKsiHeta, 2) +
                             2 * bSplinesKsi.Values[indexKsi, i] * bSplinesHeta.Values[indexHeta, j] *
                             sumdKsiHeta * sumKsidHeta / Math.Pow(sumKsiHeta, 3));
                    }
                }
            }

            return (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi,
                secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta);
        }

        /// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function mixed second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesHeta { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsi { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function second derivatives per Ksi and Heta.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] SecondDerivativeValuesKsiHeta { get; private set; }
		
        public double[,] EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricPointsCount = quadrature.IntegrationPoints.Count;
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return values;
        }

        public (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricPointsCount = quadrature.IntegrationPoints.Count;
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return (derivativeValuesKsi, derivativeValuesHeta);
        }

        public (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta) 
            EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature2D quadrature)
        {
            var parametricPointsCount = quadrature.IntegrationPoints.Count;
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var parametricGaussPointHeta = quadrature.IntegrationPoints.Select(x => x.Xi).Distinct().ToArray();
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return (secondDerivativeValuesKsi, secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta);
        }

        public double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint)
        {
            var parametricPointsCount = 1;
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return values;
        }

        public (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint)
        {
            var parametricPointsCount = 1;
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return (derivativeValuesKsi, derivativeValuesHeta);
        }

        public (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta) 
            EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint)
        {
            var parametricPointsCount = 1;
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var parametricGaussPointHeta = new double[]{naturalPoint.Eta};
            var (values, derivativeValuesKsi, derivativeValuesHeta, secondDerivativeValuesKsi, secondDerivativeValuesHeta,
                secondDerivativeValuesKsiHeta) = CalculateNurbs(parametricGaussPointKsi, parametricGaussPointHeta,
                parametricPointsCount);
            return (secondDerivativeValuesKsi, secondDerivativeValuesHeta, secondDerivativeValuesKsiHeta);
        }

        /// <summary>
        /// <see cref="Matrix"/> containing NURBS shape functions.
        /// Row represent Control Points, while columns Gauss Points.
        /// </summary>
        public double[,] Values { get; private set; }
	}
}
