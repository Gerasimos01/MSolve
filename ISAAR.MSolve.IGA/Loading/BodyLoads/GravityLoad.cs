using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.IGA.Loading.BodyLoads
{
    public class GravityLoad:IBodyLoad
    {
        private readonly double _density;
        protected double _acceleration;
        private readonly IDofType _dofType;

        public GravityLoad(double density, double acceleration, IDofType dofType)
        {
            _density = density;
            _acceleration = acceleration;
            _dofType = dofType;
        }

        public Table<INode, IDofType, double> CalculateBodyLoad(IShapeFunction3D interpolation, IQuadrature3D integration, IReadOnlyList<INode> nodes)
        {
            var loadTable = new Table<INode, IDofType, double>();
            var (derivativesKsi,derivativesHeta,derivativesZeta) =
                interpolation.EvaluateNaturalDerivativesAtGaussPoints(integration);
            var shapeFunctionNatural =
                interpolation.EvaluateFunctionsAtGaussPoints(integration);

            for (int gp = 0; gp < integration.IntegrationPoints.Count; gp++)
            {
                var jacobianMatrix = Matrix.CreateZero(3, 3);
                for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
                {
                    jacobianMatrix[0, 0] += derivativesKsi[indexNode,gp] * nodes[indexNode].X;
                    jacobianMatrix[0, 1] += derivativesKsi[indexNode,gp] * nodes[indexNode].Y;
                    jacobianMatrix[0, 2] += derivativesKsi[indexNode,gp] * nodes[indexNode].Z;

                    jacobianMatrix[1, 0] += derivativesHeta[indexNode,gp] * nodes[indexNode].X;
                    jacobianMatrix[1, 1] += derivativesHeta[indexNode,gp] * nodes[indexNode].Y;
                    jacobianMatrix[1, 2] += derivativesHeta[indexNode,gp] * nodes[indexNode].Z;

                    jacobianMatrix[2, 0] += derivativesZeta[indexNode,gp] * nodes[indexNode].X;
                    jacobianMatrix[2, 1] += derivativesZeta[indexNode,gp] * nodes[indexNode].Y;
                    jacobianMatrix[2, 2] += derivativesZeta[indexNode,gp] * nodes[indexNode].Z;
                }

                var jacdet=jacobianMatrix.CalcDeterminant();

                var weightFactor = integration.IntegrationPoints[gp].Weight;
                for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
                {
                    var node = nodes[indexNode];
                    var valueX = _acceleration * shapeFunctionNatural[indexNode,gp]* jacdet *
                                 (weightFactor*_density);
                    if (loadTable.Contains(node, _dofType))
                    {
                        loadTable[node, _dofType] += valueX;
                    }
                    else
                    {
                        loadTable.TryAdd(node, _dofType, valueX);
                    }
                }
            }

            return loadTable;
        }
    }
}
