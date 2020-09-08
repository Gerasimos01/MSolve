using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.IGA.Interfaces
{
    /// <summary>
	/// Isogeometric element interface. Implements <see cref="IElementType"/>.
	/// </summary>
	public interface IIsogeometricElement : IElementType
	{
		/// <summary>
		/// Retries the element dimensions.
		/// </summary>
		ElementDimensions ElementDimensions { get; }

		/// <summary>
		/// Retrieves the element ID.
		/// </summary>
		int ID { get; }

		/// <summary>
		/// Calculates the knot displacements for post-processing with Paraview.
		/// </summary>
		/// <param name="element">An isogeometric <see cref="Element"/>.</param>
		/// <param name="localDisplacements">A <see cref="Matrix"/> containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array calculating the displacement of the element Knots'.
		/// The rows of the matrix denote the knot numbering while the columns the displacements for each degree of freedom.</returns>
		double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements);

        double[,] CalculatePointsForPostProcessing(Element element);
    }
}
