using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Quadrature
{
    public class FullQuadrature1D : IQuadrature1D
    {
        public FullQuadrature1D(int degree, double minKnot, double maxKnot)
        {
            IList<GaussPoint> gaussPointsofElement = new List<GaussPoint>();
            int numberOfGp = degree + 1;
            double[] coordinatesKsi = new double[numberOfGp];
            double[] weightKsi = new double[numberOfGp];

            for (int i = 0; i < numberOfGp; i++)
            {
                coordinatesKsi[i] = 0.5 * (minKnot + maxKnot + (maxKnot - minKnot) * QuadratureRules.Coordinates[degree][i]);
                weightKsi[i] = 0.5 * (maxKnot - minKnot) * QuadratureRules.Weights[degree][i];
            }

            for (int i = 0; i < numberOfGp; i++)
            {
                gaussPointsofElement.Add(new GaussPoint(coordinatesKsi[i], 0,0,weightKsi[i]));
            }

            IntegrationPoints = gaussPointsofElement.ToArray();
        }

        public IReadOnlyList<GaussPoint> IntegrationPoints { get; private set; }

        public IList<GaussLegendrePoint3D> GaussPoints => IntegrationPoints
            .Select(x => new GaussLegendrePoint3D(x.Xi, x.Eta, x.Zeta, x.Weight)).ToList();

    }
}
