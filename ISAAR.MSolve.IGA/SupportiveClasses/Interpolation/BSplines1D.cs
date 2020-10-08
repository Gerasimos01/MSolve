using System;
using System.Linq;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interpolation
{
    public class BSplines1D: IShapeFunction1D
    {
        private int numberOfControlPoints;
        public int Degree { get; private set; }

        public double[] KnotValueVector { get; private set; }

        public double[] ParametricCoordinates { get; private set; }

        public double[,] Values { get; private set; }

        public double[,] DerivativeValues { get; private set; }

		public double [,] SecondDerivativeValues { get; private set; }

        public BSplines1D(int degree, double[] knotValueVector, double[] parametricCoordinates)
        {
            if (degree <= 0)
            {
                throw new ArgumentException();
            }else if(knotValueVector == null){ throw new ArgumentNullException(); }
            this.Degree = degree;
            this.KnotValueVector = knotValueVector;
            int order = Degree + 1;
            numberOfControlPoints = KnotValueVector.Length - order;
            int numberOfGaussPoints = parametricCoordinates.Length;
            (Values, DerivativeValues, SecondDerivativeValues)=CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
        }


        private (double[,] values, double[,] derivativeValues, double[,] secondDerivativeValues)
            CalculateBSplines(int numberOfControlPoints, int numberOfGaussPoints, double[] integrationCoordinates)
        {
            var values = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
            var derivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];
            var secondDerivativeValues = new double[numberOfControlPoints + Degree, numberOfGaussPoints];

            for (int i = 0; i < numberOfGaussPoints; i++)
            for (int j = 0; j < numberOfControlPoints + Degree; j++)
                if (KnotValueVector[j] <= integrationCoordinates[i] && integrationCoordinates[i] <= KnotValueVector[j + 1])
                    values[j, i] = 1;
                else
                    values[j, i] = 0;

            for (int i = 1; i <= Degree; i++)
            {
                for (int j = 0; j < numberOfControlPoints + Degree - i; j++)
                {
                    for (int k = 0; k < numberOfGaussPoints; k++)
                    {
                        double additive1 = 0;
                        double additive2 = 0;
                        double additiveDerivative1 = 0;
                        double additiveDerivative2 = 0;
                        double additiveSecondDerivative1 = 0;
                        double additiveSecondDerivative2 = 0;
                        double denominator1 = KnotValueVector[j + i] - KnotValueVector[j];
                        double denominator2 = KnotValueVector[j + i + 1] - KnotValueVector[j + 1];
                        if (denominator1 != 0)
                        {
                            additive1 = (integrationCoordinates[k] - KnotValueVector[j])
                                / denominator1 * values[j, k];
                            additiveDerivative1 = ((integrationCoordinates[k] - KnotValueVector[j])
                                * derivativeValues[j, k] + values[j, k]) / denominator1;
                            additiveSecondDerivative1 = Degree / denominator1 * derivativeValues[j, k];
                        }

                        if (denominator2 != 0)
                        {
                            additive2 = (KnotValueVector[j + i + 1] - integrationCoordinates[k])
                                / denominator2 * values[j + 1, k];
                            additiveDerivative2 = ((KnotValueVector[j + i + 1] - integrationCoordinates[k])
                                * derivativeValues[j + 1, k] - values[j + 1, k]) / denominator2;
                            additiveSecondDerivative2 = -Degree / denominator2 * derivativeValues[j + 1, k];
                        }

                        values[j, k] = additive1 + additive2;
                        derivativeValues[j, k] = additiveDerivative1 + additiveDerivative2;
                        secondDerivativeValues[j, k] = additiveSecondDerivative1 + additiveSecondDerivative2;
                    }
                }
            }

            return (values, derivativeValues, secondDerivativeValues);
        }

        public double[,] EvaluateFunctionsAtGaussPoints(IQuadrature1D quadrature)
        {
            var numberOfGaussPoints = quadrature.IntegrationPoints.Count;
            var parametricCoordinates = quadrature.IntegrationPoints.Select(x => x.Xi).ToArray();
            var (values, derivativeValues, secondDerivativeValues)=CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return values;
        }

        public double[,] EvaluateNaturalDerivativesAtGaussPoints(IQuadrature1D quadrature)
        {
            var numberOfGaussPoints = quadrature.IntegrationPoints.Count;
            var parametricCoordinates = quadrature.IntegrationPoints.Select(x => x.Xi).ToArray();
            var (values, derivativeValues, secondDerivativeValues)=CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return derivativeValues;
        }

        public double[,] EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature1D quadrature)
        {
            var numberOfGaussPoints = quadrature.IntegrationPoints.Count;
            var parametricCoordinates = quadrature.IntegrationPoints.Select(x => x.Xi).ToArray();
            var (_, _, secondDerivativeValues)=CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return secondDerivativeValues;
        }

        public double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint)
        {
            var numberOfGaussPoints = 1;
            var parametricCoordinates = new double[] {naturalPoint.Xi};
            var (values, _, _) =CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return values;
        }

        public double[,] EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint)
        {
            var numberOfGaussPoints = 1;
            var parametricCoordinates = new double[] {naturalPoint.Xi};
            var (_, derivativeValues, _) =CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return derivativeValues;
        }

        public double[,] EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint)
        {
            var numberOfGaussPoints = 1;
            var parametricCoordinates = new double[] {naturalPoint.Xi};
            var (_, _, secondDerivativeValues) =CalculateBSplines(numberOfControlPoints, numberOfGaussPoints, parametricCoordinates);
            return secondDerivativeValues;
        }
    }
}
