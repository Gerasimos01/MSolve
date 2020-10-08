using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.LoadElements
{
    public class SurfaceLoadElement:ILoad
    {
        private readonly ISurfaceLoad surfaceLoad;
        private readonly IShapeFunction2D isoparametricInterpolation2D;
        private readonly IQuadrature2D quadrature2D;
        private readonly IReadOnlyList<INode> nodes;

        public SurfaceLoadElement(ISurfaceLoad surfaceLoad, IShapeFunction2D isoparametricInterpolation2D,
            IQuadrature2D quadrature2D, IReadOnlyList<INode> nodes)
        {
            this.surfaceLoad = surfaceLoad;
            this.isoparametricInterpolation2D = isoparametricInterpolation2D;
            this.quadrature2D = quadrature2D;
            this.nodes = nodes;
        }

        public Table<INode, IDofType, double> CalculateLoad() =>
            surfaceLoad.CalculateSurfaceLoad(isoparametricInterpolation2D, quadrature2D, nodes);

    }
}
