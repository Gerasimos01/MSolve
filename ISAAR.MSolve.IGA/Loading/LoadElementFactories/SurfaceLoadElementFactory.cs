using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.LoadElementFactories
{
    public enum NurbsSolidFaces
    {
        Left,
        Right,
        Bottom,
        Top,
        Front,
        Back
    }

    public class SurfaceLoadElementFactory
    {
        private readonly ISurfaceLoad load;
        private readonly Model model;

        public SurfaceLoadElementFactory(ISurfaceLoad load, Model model)
        {
            this.load = load;
            this.model = model;
        }

        public List<SurfaceLoadElementFactory> CreateElements(NurbsSolidGeometryDto solidGeometry, NurbsSolidFaces face)
        {
            throw new NotImplementedException();
        }

        public List<SurfaceLoadElementFactory> CreateElements(NurbsSurfaceGeometryDto surfaceGeometry)
        {
            throw new NotImplementedException();
        }
    }
}
