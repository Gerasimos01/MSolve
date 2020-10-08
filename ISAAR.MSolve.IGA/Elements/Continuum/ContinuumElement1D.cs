using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using Element = ISAAR.MSolve.IGA.Entities.Element;

namespace ISAAR.MSolve.IGA.Elements.Continuum
{
	public class ContinuumElement1D : Element, IStructuralIsogeometricElement
	{

        public ContinuumElement1D(IShapeFunction1D shapeFunctions, GaussLegendrePoint3D[] gaussPoints,
            int numberOfCPHeta, int numberOfCpZeta)
        {
            _numberOfCPHeta = numberOfCPHeta;
            _numberOfCPZeta = numberOfCpZeta;
            _shapeFunctions = shapeFunctions;
            _gaussPoints = gaussPoints;
        }

		protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        private readonly IShapeFunction1D _shapeFunctions;
		private readonly GaussLegendrePoint3D[] _gaussPoints;
		private IDofType[][] _dofTypes;
        private int _numberOfCPHeta;
        private int _numberOfCPZeta;

        /// <summary>
        /// Retrieves the type of Finite Element used. Since the element is Isogeometric its type is defined as unknown.
        /// </summary>
        public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Property that Polynomial degree of the NURBS shape functions.
		/// </summary>
		public int Degree { get; set; }

		#region IStructuralIsogeometricElement

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		/// <summary>
		/// Retrieves the number of Dimensions of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.OneD;

		/// <summary>
		/// Boolean property that determines whether the material used for this elements has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Calculates the forces applies to an <see cref="NurbsElement1D"/> due to <see cref="FEM.Entities.MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <param name="loads">A list of <see cref="FEM.Entities.MassAccelerationLoad"/>. For more info see <seealso cref="FEM.Entities.MassAccelerationLoad"/></param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		//public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) => Array.Empty<double>();

		/// <summary>
		/// Calculates displacements of knots for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements) =>
			Array.Empty<double>();

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(Element element, double[] localDisplacements) =>
			Array.Empty<double>();

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements) => throw new NotImplementedException();

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Clear any saved material states of the element.
		/// </summary>
		public void ClearMaterialStresses() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/></param>
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
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() => throw new NotImplementedException();

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() => throw new NotImplementedException();

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="NurbsElement1D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="NurbsElement1D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element) => throw new NotImplementedException();
		

        public double[,] CalculatePointsForPostProcessing(Element element)
        {
            throw new NotImplementedException();
        }

        #endregion IStructuralIsogeometricElement
    }
}
