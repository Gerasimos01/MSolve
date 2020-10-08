using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Interfaces
{
    public interface IShapeFunction2D
    {
        double[,] Values { get; }
        double[,] DerivativeValuesHeta { get;  }
        double[,] DerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesHeta { get; }
        double[,] SecondDerivativeValuesKsi { get; }
        double[,] SecondDerivativeValuesKsiHeta { get; }
        
        double[,] EvaluateFunctionsAt(NaturalPoint naturalPoint);
        (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAt(NaturalPoint naturalPoint);

        (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta)
            EvaluateNaturalSecondDerivativesAt(NaturalPoint naturalPoint);

        double[,] EvaluateFunctionsAtGaussPoints(IQuadrature2D quadrature);
        (double[,] derivativesKsi, double[,] derivativesHeta) EvaluateNaturalDerivativesAtGaussPoints(IQuadrature2D quadrature);

        (double[,] secondDerivativesKsi, double[,] secondDerivativesHeta, double[,] secondDerivativesKsiHeta)
            EvaluateNaturalSecondDerivativesAtGaussPoints(IQuadrature2D quadrature);
    }
}