using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Loading.LineLoads
{
    public class DistributedLineLoad : ILineLoad
    {
        protected double _loadX;
        protected double _loadY;
        protected double _loadZ;

        public DistributedLineLoad(double loadX, double loadY, double loadZ)
        {
            _loadX = loadX;
            _loadY = loadY;
            _loadZ = loadZ;
        }

        public Table<INode, IDofType, double> CalculateLineLoad(IShapeFunction1D interpolation, IQuadrature1D integration, IReadOnlyList<INode> nodes)
        {
            var loadTable = new Table<INode, IDofType, double>();
			var shapeGradientsNatural =
				interpolation.EvaluateNaturalDerivativesAtGaussPoints(integration);
			var shapeFunctionNatural =
				interpolation.EvaluateFunctionsAtGaussPoints(integration);


			for (int gp = 0; gp < integration.IntegrationPoints.Count; gp++)
			{
				var jacobianMatrix = Vector.CreateZero(3);
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					jacobianMatrix[0] += shapeGradientsNatural[indexNode,gp] * nodes[indexNode].X;
					jacobianMatrix[1] += shapeGradientsNatural[indexNode,gp] * nodes[indexNode].Y;
					jacobianMatrix[2] += shapeGradientsNatural[indexNode,gp]* nodes[indexNode].Z;
				}

				var jacdet = jacobianMatrix.Norm2();

				var weightFactor = integration.IntegrationPoints[gp].Weight;
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					var node = nodes[indexNode];
					var valueX = _loadX * shapeFunctionNatural[indexNode,gp] * jacdet * weightFactor;
					var valueY = _loadY * shapeFunctionNatural[indexNode,gp] * jacdet * weightFactor;
					var valueZ = _loadZ * shapeFunctionNatural[indexNode,gp] * jacdet * weightFactor;
					if (loadTable.Contains(node, StructuralDof.TranslationX))
					{
						loadTable[node, StructuralDof.TranslationX] += valueX;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationX, valueX);
					}

					if (loadTable.Contains(node, StructuralDof.TranslationY))
					{
						loadTable[node, StructuralDof.TranslationY] += valueY;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationY, valueY);
					}

					if (loadTable.Contains(node, StructuralDof.TranslationZ))
					{
						loadTable[node, StructuralDof.TranslationZ] += valueZ;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationZ, valueZ);
					}
				}
			}

            return loadTable;
        }
    }
}
