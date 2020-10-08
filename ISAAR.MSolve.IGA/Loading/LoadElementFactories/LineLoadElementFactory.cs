using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.Loading.LoadElements;

namespace ISAAR.MSolve.IGA.Loading.LoadElementFactories
{

    public enum NurbsSurfaceEdges
    {
        Left = 0,
        Right = 1,
        Bottom = 2,
        Top = 3
    }

    public enum NurbsSolidEdges
    {
        Edge1,
        Edge2,
        Edge3,
        Edge4,
        Edge5,
        Edge6,
        Edge7,
        Edge8,
        Edge9,
        Edge10,
        Edge11,
        Edge12,
    }


    public class LineLoadElementFactory
    {
        private readonly ILineLoad load;
        private readonly Model model;

        public LineLoadElementFactory(ILineLoad load, Model model)
        {
            this.load = load;
            this.model = model;
        }

        public List<LineLoadElement> CreateElements(NurbsSurfaceGeometryDto nurbsSurfaceGeometry,
            NurbsSurfaceEdges edge)
        {
            throw new NotImplementedException();
        }

        public List<LineLoadElement> CreatElements(NurbsSolidGeometryDto nurbsSolidGeometry,
            NurbsSolidEdges edge)
        {
            throw new NotImplementedException();
        }
    }
}
