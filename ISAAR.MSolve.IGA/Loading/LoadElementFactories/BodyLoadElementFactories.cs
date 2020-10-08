using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.Loading.LoadElements;

namespace ISAAR.MSolve.IGA.Loading.LoadElementFactories
{
    public class BodyLoadElementFactories
    {
        private readonly IBodyLoad load;
        private readonly Model model;

        public BodyLoadElementFactories(IBodyLoad load, Model model)
        {
            this.load = load;
            this.model = model;
        }

        public List<BodyLoadElement> CreateElements(NurbsSolidGeometryDto nurbsGeometry)
        {
            throw  new NotImplementedException();
        }
    }
}
