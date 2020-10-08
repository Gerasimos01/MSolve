using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction3D
    {
        double[,] Values { get; }

        double[,] DerivativeValuesKsi { get; }
        double[,] DerivativeValuesHeta { get; }
        double[,] DerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesHeta { get; }
        double[,] SecondDerivativeValuesZeta { get; }

        double[,] SecondDerivativeValuesKsiHeta { get; }
        double[,] SecondDerivativeValuesHetaZeta { get; }
        double[,] SecondDerivativeValuesKsiZeta { get; }

        double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint);

        (double[,] derivativesKsi, double[,] derivativesHeta, double[,] derivativesZeta) EvaluateNaturalDerivativesAt(
            NaturalPoint naturalPoint);

        (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesZeta,
            double[,] secondDerivativesKsiHeta, double[,] secondDerivativesHetaZeta, double[,] secondDerivativesKsiZeta)
            EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint);

        double[,] EvaluateFunctionsAtGaussPoints(IQuadrature3D quadrature);
        (double[,] derivativesKsi, double[,] derivativesHeta, double[,] derivativesZeta) EvaluateNaturalDerivativesAtGaussPoints(IQuadrature3D quadrature);

        (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesZeta,
            double[,] secondDerivativesKsiHeta, double[,] secondDerivativesHetaZeta, double[,] secondDerivativesKsiZeta)
            EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature3D quadrature);

    }
}