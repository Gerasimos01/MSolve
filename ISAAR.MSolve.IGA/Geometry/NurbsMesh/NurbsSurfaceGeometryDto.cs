namespace ISAAR.MSolve.IGA.Geometry.NurbsMesh
{
    public enum NurbsGeometryType
    {
        Plane=0,
        ShellLinear=1,
        ShellNonLinear=2,
    }

    public class NurbsSurfaceGeometryDto
    {
        public int ID { get; set; }

        public double Thickness { get; set; }

        public int NumberOfCpKsi { get; set; }

        public int NumberOfCpHeta { get; set; }

        public int DegreeKsi { get; set; }

        public int DegreeHeta { get; set; }

        public double[] KnotValueVectorKsi { get; set; }

        public double[] KnotValueVectorHeta { get; set; }

        public int[] ControlPointIDs { get; set; }

        public NurbsGeometryType GeometryType { get; set; }
    }
}
