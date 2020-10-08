using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.Interfaces
{
    public interface IBodyLoad
    {
        Table<INode, IDofType, double> CalculateBodyLoad(IShapeFunction3D interpolation,
            IQuadrature3D integration, IReadOnlyList<INode> nodes);
    }
}