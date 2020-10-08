using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.IGA.Entities;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Quadrature
{
    public class FullQuadrature2D:IQuadrature2D
    {
        public FullQuadrature2D(int degreeKsi, int degreeHeta, IList<Knot> knotsOfElement)
        {
            IList<GaussLegendrePoint3D> gaussPointsofElement = new List<GaussLegendrePoint3D>();
            int numberOfGPKsi = degreeKsi + 1;
            int numberOfGPHeta = degreeHeta + 1;
            double[] coordinatesKsi = new double[numberOfGPKsi];
            double[] weightKsi = new double[numberOfGPKsi];
            double[] coordinatesHeta = new double[numberOfGPHeta];
            double[] weightHeta = new double[numberOfGPHeta];

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[2].Ksi
                                                                 + (knotsOfElement[2].Ksi - knotsOfElement[0].Ksi) * QuadratureRules.Coordinates[degreeKsi][i]);
                weightKsi[i] = 0.5 * ((knotsOfElement[2].Ksi - knotsOfElement[0].Ksi) * QuadratureRules.Weights[degreeKsi][i]);
            }

            for (int i = 0; i < numberOfGPHeta; i++)
            {
                coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[1].Heta
                                                                   + (knotsOfElement[1].Heta - knotsOfElement[0].Heta) * QuadratureRules.Coordinates[degreeHeta][i]);
                weightHeta[i] = 0.5 * ((knotsOfElement[1].Heta - knotsOfElement[0].Heta) * QuadratureRules.Weights[degreeHeta][i]);
            }

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    gaussPointsofElement.Add(new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], 0, weightKsi[i] * weightHeta[j]));
                }
            }

            GaussPoints = gaussPointsofElement;
            IntegrationPoints = GaussPoints.Select(x => new GaussPoint(x.Ksi, x.Heta, x.Zeta, x.WeightFactor)).ToArray();
        }

        public IReadOnlyList<GaussPoint> IntegrationPoints { get; private set; }

        public IList<GaussLegendrePoint3D> GaussPoints { get; private set; }
    }
}
