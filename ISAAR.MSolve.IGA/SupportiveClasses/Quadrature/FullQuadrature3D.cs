using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.IGA.Entities;

namespace ISAAR.MSolve.IGA.SupportiveClasses.Quadrature
{
    public class FullQuadrature3D:IQuadrature3D
    {
        public FullQuadrature3D(int degreeKsi, int degreeHeta, int degreeZeta, IList<Knot> knotsOfElement)
        {
            int numberOfGPKsi = degreeKsi + 1;
            int numberOfGPHeta = degreeHeta + 1;
            int numberOfGPZeta = degreeZeta + 1;
            var gaussPointsofElement = new GaussLegendrePoint3D[numberOfGPKsi * numberOfGPHeta * numberOfGPZeta];
            double[] coordinatesKsi = new double[numberOfGPKsi];
            double[] weightKsi = new double[numberOfGPKsi];
            double[] coordinatesHeta = new double[numberOfGPHeta];
            double[] weightHeta = new double[numberOfGPHeta];
            double[] coordinatesZeta = new double[numberOfGPZeta];
            double[] weightZeta = new double[numberOfGPZeta];

            for (int i = 0; i < numberOfGPKsi; i++)
            {
                coordinatesKsi[i] = 0.5 * (knotsOfElement[0].Ksi + knotsOfElement[4].Ksi
                                                                 + (knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * QuadratureRules.Coordinates[degreeKsi][i]);
                weightKsi[i] = 0.5 * ((knotsOfElement[4].Ksi - knotsOfElement[0].Ksi) * QuadratureRules.Weights[degreeKsi][i]);
            }

            for (int i = 0; i < numberOfGPHeta; i++)
            {
                coordinatesHeta[i] = 0.5 * (knotsOfElement[0].Heta + knotsOfElement[2].Heta
                                                                   + (knotsOfElement[2].Heta - knotsOfElement[0].Heta) * QuadratureRules.Coordinates[degreeHeta][i]);
                weightHeta[i] = 0.5 * ((knotsOfElement[2].Heta - knotsOfElement[0].Heta) * QuadratureRules.Weights[degreeHeta][i]);
            }

            for (int i = 0; i < numberOfGPZeta; i++)
            {
                coordinatesZeta[i] = 0.5 * (knotsOfElement[0].Zeta + knotsOfElement[1].Zeta
                                                                   + (knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * QuadratureRules.Coordinates[degreeZeta][i]);
                weightZeta[i] = 0.5 * ((knotsOfElement[1].Zeta - knotsOfElement[0].Zeta) * QuadratureRules.Weights[degreeZeta][i]);
            }

            var index = 0;
            for (int i = 0; i < numberOfGPKsi; i++)
            {
                for (int j = 0; j < numberOfGPHeta; j++)
                {
                    for (int k = 0; k < numberOfGPZeta; k++)
                    {
                        gaussPointsofElement[index++] = new GaussLegendrePoint3D(coordinatesKsi[i], coordinatesHeta[j], coordinatesZeta[k], weightKsi[i] * weightHeta[j] * weightZeta[k]);
                    }
                }
            }

            GaussPoints = gaussPointsofElement;
            IntegrationPoints = GaussPoints.Select(x => new GaussPoint(x.Ksi, x.Heta, x.Zeta, x.WeightFactor)).ToArray();
        }

        public IReadOnlyList<GaussPoint> IntegrationPoints { get; }

        public IList<GaussLegendrePoint3D> GaussPoints { get; }
    }
}
