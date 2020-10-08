using System.Collections.Generic;

namespace ISAAR.MSolve.IGA.Geometry.NurbsMesh
{
    public class ModelDto
    {
        public List<NurbsSurfaceGeometryDto> NurbsSurfacePatches { get; set; }
        public List<NurbsSolidGeometryDto> NurbsSolidPatches { get; set; }
        public List<ControlPointDto> ControlPoints { get; set; }
    }
}
