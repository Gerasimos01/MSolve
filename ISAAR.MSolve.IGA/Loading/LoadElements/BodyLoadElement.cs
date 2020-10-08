using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.LoadElements
{
    public class BodyLoadElement:ILoad
    {
        private readonly IBodyLoad _bodyLoad;
        private readonly IShapeFunction3D _interpolation;
        private readonly IQuadrature3D _quadrature;
        private readonly IReadOnlyList<INode> _nodes;

        public BodyLoadElement(IBodyLoad bodyLoad, IShapeFunction3D interpolation3D, 
            IQuadrature3D quadrature3D, IReadOnlyList<INode> nodes)
        {
            _bodyLoad = bodyLoad;
            _interpolation = interpolation3D;
            _quadrature = quadrature3D;
            _nodes = nodes;
        }

        public Table<INode, IDofType, double> CalculateLoad()=>
            _bodyLoad.CalculateBodyLoad(_interpolation, _quadrature, _nodes);
        
    }
}
