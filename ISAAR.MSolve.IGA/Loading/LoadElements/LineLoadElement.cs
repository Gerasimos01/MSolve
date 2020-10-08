using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.LoadElements
{
    public class LineLoadElement:ILoad
    {
        private readonly ILineLoad _lineLoad;
        private readonly IShapeFunction1D _interpolation1D;
        private readonly IQuadrature1D _quadrature1D;
        private readonly IReadOnlyList<INode> _nodes;

        public LineLoadElement(ILineLoad lineLoad, IShapeFunction1D interpolation1D,
            IQuadrature1D quadrature1D, IReadOnlyList<INode> nodes)
        {
            _lineLoad = lineLoad;
            _interpolation1D = interpolation1D;
            _quadrature1D = quadrature1D;
            _nodes = nodes;
        }

        public Table<INode, IDofType, double> CalculateLoad()
            => _lineLoad.CalculateLineLoad(_interpolation1D, _quadrature1D, _nodes);
    }
}
