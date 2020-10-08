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
	/// One-dimensional NURBS shape functions.
	/// </summary>
	public class Nurbs1D:IShapeFunction1D
	{
        private readonly int degree;
        private readonly double[] knotValueVector;
        private readonly ControlPoint[] controlPoints;
        private readonly GaussLegendrePoint3D[] gaussPoints;

        /// <summary>
		/// Defines an 1D NURBS shape function for an element.
		/// </summary>
		/// <param name="element">An <see cref="Element"/> of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="controlPoints">A <see cref="List{T}"/> containing the control points of the element.</param>
		public Nurbs1D(int degree, double[] KnotValueVector, ControlPoint[] controlPoints, GaussLegendrePoint3D[] gaussPoints)
        {
            this.degree = degree;
            this.knotValueVector = KnotValueVector;
            this.controlPoints = controlPoints;

            this.gaussPoints = gaussPoints;
            var parametricGaussPointKsi = gaussPoints.Select(gp => gp.Ksi).ToArray();
            var numberOfGaussPoints = parametricGaussPointKsi.Length;
			var (values, derivativeValues) = CalculateNurbs1D(degree, KnotValueVector, controlPoints, gaussPoints, parametricGaussPointKsi, numberOfGaussPoints);
            Values = values;
            DerivativeValues = derivativeValues;
        }

        private (double[,] values, double[,] derivativeValues) CalculateNurbs1D(int degree, double[] KnotValueVector,
            ControlPoint[] controlPoints, GaussLegendrePoint3D[] gaussPoints, double[] parametricGaussPointKsi,
            int numberOfGaussPoints)
        {
            var bsplinesKsi = new BSplines1D(degree, KnotValueVector, parametricGaussPointKsi);
            var values = new double[controlPoints.Length, gaussPoints.Length];
            var derivativeValues = new double[controlPoints.Length, gaussPoints.Length];
            for (int i = 0; i < numberOfGaussPoints; i++)
            {
                double sumKsi = 0;
                double sumdKsi = 0;

                for (int j = 0; j < controlPoints.Length; j++)
                {
                    int indexKsi = controlPoints[j].ID;
                    sumKsi += bsplinesKsi.Values[indexKsi, i] * controlPoints[j].WeightFactor;
                    sumdKsi += bsplinesKsi.DerivativeValues[indexKsi, i] * controlPoints[j].WeightFactor;
                }

                for (int j = 0; j < controlPoints.Length; j++)
                {
                    int indexKsi = controlPoints[j].ID;
                    values[j, i] = bsplinesKsi.Values[indexKsi, i] * controlPoints[j].WeightFactor / sumKsi;
                    derivativeValues[j, i] = controlPoints[j].WeightFactor *
                        (bsplinesKsi.DerivativeValues[indexKsi, i] * sumKsi -
                         bsplinesKsi.Values[indexKsi, i] * sumdKsi) / Math.Pow(sumKsi, 2);
                }
            }

            return (values, derivativeValues);
        }

        /// <summary>
		/// <see cref="Matrix"/> containing NURBS shape function derivatives.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] DerivativeValues { get; private set; }

		/// <summary>
		/// <see cref="Matrix"/> containing NURBS shape functions.
		/// Row represent Control Points, while columns Gauss Points.
		/// </summary>
		public double[,] Values { get; private set; }

		public double[,] SecondDerivativeValues => throw new NotImplementedException();

        public double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint)
        {
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var (values, _) = CalculateNurbs1D(degree, knotValueVector, controlPoints, gaussPoints, parametricGaussPointKsi, parametricGaussPointKsi.Length);

            return values;
        }

        public double[,] EvaluateFunctionsAtGaussPoints(IQuadrature1D quadrature)
        {
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).ToArray();
            var (values, _) = CalculateNurbs1D(degree, knotValueVector, controlPoints, gaussPoints, parametricGaussPointKsi, parametricGaussPointKsi.Length);

            return values;
        }

        public double[,] EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint)
        {
            var parametricGaussPointKsi = new double[]{naturalPoint.Xi};
            var (_, derivativeValues) = CalculateNurbs1D(degree, knotValueVector, controlPoints, gaussPoints, parametricGaussPointKsi, parametricGaussPointKsi.Length);

            return derivativeValues;
        }

        public double[,] EvaluateNaturalDerivativesAtGaussPoints(IQuadrature1D quadrature)
        {
            var parametricGaussPointKsi = quadrature.IntegrationPoints.Select(x => x.Xi).ToArray();
            var (_, derivativeValues) = CalculateNurbs1D(degree, knotValueVector, controlPoints, gaussPoints, parametricGaussPointKsi, parametricGaussPointKsi.Length);

            return derivativeValues;
        }

        public double[,] EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint)
        {
            throw new NotImplementedException();
        }

        public double[,] EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature1D quadrature)
        {
            throw new NotImplementedException();
        }
    }
}
