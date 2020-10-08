using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using Element = ISAAR.MSolve.IGA.Entities.Element;


namespace ISAAR.MSolve.IGA.Elements.Continuum
{
	public class ContinuumElement2D : Element, IStructuralIsogeometricElement
	{
        public ContinuumElement2D(IContinuumMaterial2D material, 
            IShapeFunction2D shapeFunctions, GaussLegendrePoint3D[] gaussPoints,
            double thickness)
        {
            _material = material;
            _shapeFunctions = shapeFunctions;
            _gaussPoints = gaussPoints;
            _thickness = thickness;
        }


		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private IDofType[][] _dofTypes;
        private readonly IContinuumMaterial2D _material;
        public readonly IShapeFunction2D _shapeFunctions;
        private readonly GaussLegendrePoint3D[] _gaussPoints;
        private readonly double _thickness;

        /// <summary>
        /// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
        /// </summary>
        public CellType CellType { get; } = CellType.Unknown;

		#region IStructuralIsogeometricElement

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Retrieves the number of Dimensions of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
            throw new NotImplementedException();
			//Contract.Requires(element != null, "The element cannot be null");
   //         var nurbsElement = (NURBSElement2D)element;
			//var elementControlPoints = nurbsElement.ControlPoints.ToArray();
			//var elemenetKnots = nurbsElement.Knots.ToArray();
			//var knotParametricCoordinatesKsi = Vector.CreateFromArray(new double[] { elemenetKnots[0].Ksi, elemenetKnots[2].Ksi });
			//var knotParametricCoordinatesHeta = Vector.CreateFromArray(new double[] { elemenetKnots[0].Heta, elemenetKnots[1].Heta });
			//Nurbs2D nurbs = new Nurbs2D(nurbsElement, elementControlPoints, knotParametricCoordinatesKsi, knotParametricCoordinatesHeta);
			//var knotDisplacements = new double[4, 2];
			//var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
			//for (int j = 0; j < elemenetKnots.Length; j++)
			//{
			//	for (int i = 0; i < elementControlPoints.Length; i++)
			//	{
			//		knotDisplacements[paraviewKnotRenumbering[j], 0] += nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
			//		knotDisplacements[paraviewKnotRenumbering[j], 1] += nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
			//	}
			//}

			//return knotDisplacements;
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom</returns>
		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="NURBSElement2D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="ControlPoint"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");
			_dofTypes = new IDofType[element.Nodes.Count][];
			for (var i = 0; i < element.Nodes.Count; i++)
			{
				_dofTypes[i] = ControlPointDofTypes;
			}

			return _dofTypes;
		}

		/// <summary>
		/// Calculates the mass matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NURBSElement2D"/>.</returns>
		public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NURBSElement2D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element)
		{
			Contract.Requires(element != null, "The element cannot be null");
            var numberOfCP = element.Nodes.Count;
			var stiffnessMatrixElement = Matrix.CreateZero(
                numberOfCP * 2, numberOfCP * 2);
            var elementControlPoints = element.Nodes.ToArray();

			for (var j = 0; j < _gaussPoints.Length; j++)
			{
				var jacobianMatrix = Matrix.CreateZero(2, 2);

				for (int k = 0; k < numberOfCP; k++)
				{
					jacobianMatrix[0, 0] += _shapeFunctions.DerivativeValuesKsi[k, j] * elementControlPoints[k].X;
					jacobianMatrix[0, 1] += _shapeFunctions.DerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
					jacobianMatrix[1, 0] += _shapeFunctions.DerivativeValuesHeta[k, j] * elementControlPoints[k].X;
					jacobianMatrix[1, 1] += _shapeFunctions.DerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				}

				double jacdet = (jacobianMatrix[0, 0] * jacobianMatrix[1, 1])
								- (jacobianMatrix[1, 0] * jacobianMatrix[0, 1]);

				Matrix B1 = Matrix.CreateZero(3, 4);

				B1[0, 0] += jacobianMatrix[1, 1] / jacdet;
				B1[0, 1] += -jacobianMatrix[0, 1] / jacdet;
				B1[1, 2] += -jacobianMatrix[1, 0] / jacdet;
				B1[1, 3] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 0] += -jacobianMatrix[1, 0] / jacdet;
				B1[2, 1] += jacobianMatrix[0, 0] / jacdet;
				B1[2, 2] += jacobianMatrix[1, 1] / jacdet;
				B1[2, 3] += -jacobianMatrix[0, 1] / jacdet;

				Matrix B2 = Matrix.CreateZero(4, 2 * numberOfCP);
				for (int column = 0; column < 2 * numberOfCP; column += 2)
				{
					B2[0, column] += _shapeFunctions.DerivativeValuesKsi[column / 2, j];
					B2[1, column] += _shapeFunctions.DerivativeValuesHeta[column / 2, j];
					B2[2, column + 1] += _shapeFunctions.DerivativeValuesKsi[column / 2, j];
					B2[3, column + 1] += _shapeFunctions.DerivativeValuesHeta[column / 2, j];
				}

				Matrix B = B1 * B2;
				IMatrixView elasticityMatrix = _material.ConstitutiveMatrix;
				Matrix stiffnessMatrixGaussPoint = B.ThisTransposeTimesOtherTimesThis(elasticityMatrix);
				stiffnessMatrixGaussPoint *= jacdet * _gaussPoints[j].WeightFactor * _thickness;

				for (int m = 0; m < numberOfCP * 2; m++)
				{
					for (int n = 0; n < numberOfCP * 2; n++)
					{
						stiffnessMatrixElement[m, n] += stiffnessMatrixGaussPoint[m, n];
					}
				}
			}

			return stiffnessMatrixElement;
		}

		private static Matrix JacobianMatrixForLoadCalculation(Element element, IShapeFunction2D nurbs, int j, out double xGaussPoint,
			out double yGaussPoint, out double zGaussPoint)
		{
			var elementControlPoints = element.ControlPoints.ToArray();
			Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
			xGaussPoint = 0;
			yGaussPoint = 0;
			zGaussPoint = 0;
			for (int k = 0; k < elementControlPoints.Length; k++)
			{
				xGaussPoint += nurbs.Values[k, j] * elementControlPoints[k].X;
				yGaussPoint += nurbs.Values[k, j] * elementControlPoints[k].Y;
				zGaussPoint += nurbs.Values[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[0, 0] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].X;
				jacobianMatrix[0, 1] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[0, 2] += nurbs.DerivativeValuesKsi[k, j] * elementControlPoints[k].Z;
				jacobianMatrix[1, 0] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].X;
				jacobianMatrix[1, 1] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].Y;
				jacobianMatrix[1, 2] += nurbs.DerivativeValuesHeta[k, j] * elementControlPoints[k].Z;
			}

			return jacobianMatrix;
		}

		#endregion IStructuralIsogeometricElement

		/// <summary>
		/// Calculates the forces applies to an <see cref="NURBSElement2D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NURBSElement2D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		//public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => throw new NotImplementedException();

		/// <summary>
		/// Clear any saved material states of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() => throw new NotImplementedException();

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() => throw new NotImplementedException();
        
		private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element, int degreeKsi, int degreeHeta)
		{
			GaussQuadrature gauss = new GaussQuadrature();
			return gauss.CalculateElementGaussPoints(degreeKsi, degreeHeta, element.Knots.ToArray());
		}

        public double[,] CalculatePointsForPostProcessing(Element element)
        {
            throw new NotImplementedException();
        }
    }
}
