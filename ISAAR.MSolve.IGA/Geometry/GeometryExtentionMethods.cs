using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.IGA.Loading.LoadElementFactories;
using ISAAR.MSolve.IGA.Loading.LoadElements;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.IGA.SupportiveClasses.Quadrature;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Geometry
{
    public static class GeometryExtentionMethods
    {
        public static List<ILoad> CreateLoadForEdge(this NurbsSurfaceGeometryDto surfaceGeometry,Model model, NurbsSurfaceEdges edge,
            ILineLoad lineLoad)
        {
            var (KnotValueVector, degree, controlPoints) = FindEdgeData2D(surfaceGeometry, model, edge);

            #region Knots

            var singleKnotValuesKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[0];
            var knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = 0.0, Zeta = 0.0 });
                id++;
            }

            #endregion Knots

            #region Elements

            Vector multiplicityKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            if (numberOfElementsKsi == 0)
            {
                throw new ArgumentNullException("Number of Elements should be defined before Element Connectivity");
            }

            var edgeLoads= new List<ILoad>();
            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                IList<Knot> knotsOfElement = new List<Knot>
                {
                    knots[i],
                    knots[i + 1]
                };

                int multiplicityElementKsi = 0;
                if (multiplicityKsi[i + 1] - degree > 0)
                {
                    multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degree;
                }

                int nurbsSupportKsi = degree + 1;

                IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                for (int k = 0; k < nurbsSupportKsi; k++)
                {
                    int controlPointID = i + multiplicityElementKsi + k;
                    elementControlPoints.Add(controlPoints[controlPointID]);
                }

                var quadrature= new FullQuadrature1D(degree, knots[i].Ksi, knots[i+1].Ksi);

                var nurbs = new Nurbs1D(degree, KnotValueVector.CopyToArray(), elementControlPoints.Select(controlPoint => new ControlPoint()
                {
                    ID = (edge == NurbsSurfaceEdges.Left || edge == NurbsSurfaceEdges.Right)
                        ? controlPoint.ID % surfaceGeometry.NumberOfCpHeta
                        : controlPoint.ID / surfaceGeometry.NumberOfCpHeta,
                    Ksi = controlPoint.Ksi,
                    Heta = controlPoint.Heta,
                    Zeta = controlPoint.Zeta,
                    X = controlPoint.X,
                    Y = controlPoint.Y,
                    Z = controlPoint.Z,
                    WeightFactor = controlPoint.WeightFactor,
                }).ToArray(), quadrature.GaussPoints.ToArray());
                var lineLoadElement = new LineLoadElement(lineLoad, nurbs, quadrature, elementControlPoints.ToArray());
                edgeLoads.Add(lineLoadElement);
            }

            return edgeLoads;

            #endregion Elements
        }

        private static (Vector KnotValueVector, int degree, List<ControlPoint> controlPoints) FindEdgeData2D(NurbsSurfaceGeometryDto surfaceGeometry,
            Model model, NurbsSurfaceEdges edge)
        {
            Vector KnotValueVector = null;
            int degree = 0;
            int numberOfControlPoints;
            var controlPoints = new List<ControlPoint>();
            switch (edge)
            {
                case NurbsSurfaceEdges.Left:
                    KnotValueVector = LinearAlgebra.Vectors.Vector.CreateFromArray(surfaceGeometry.KnotValueVectorHeta);
                    degree = surfaceGeometry.DegreeHeta;
                    numberOfControlPoints = surfaceGeometry.NumberOfCpHeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                        controlPoints.Add(model.ControlPointsDictionary[surfaceGeometry.ControlPointIDs[i]]);
                    break;
                case NurbsSurfaceEdges.Right:
                    KnotValueVector = LinearAlgebra.Vectors.Vector.CreateFromArray(surfaceGeometry.KnotValueVectorHeta);
                    degree = surfaceGeometry.DegreeHeta;
                    numberOfControlPoints = surfaceGeometry.NumberOfCpHeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i + surfaceGeometry.NumberOfCpHeta * (surfaceGeometry.NumberOfCpKsi - 1);
                        controlPoints.Add(model.ControlPointsDictionary[surfaceGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSurfaceEdges.Bottom:
                    KnotValueVector = LinearAlgebra.Vectors.Vector.CreateFromArray(surfaceGeometry.KnotValueVectorKsi);
                    degree = surfaceGeometry.DegreeKsi;
                    numberOfControlPoints = surfaceGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * surfaceGeometry.NumberOfCpHeta;
                        controlPoints.Add(model.ControlPointsDictionary[surfaceGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSurfaceEdges.Top:
                    KnotValueVector = LinearAlgebra.Vectors.Vector.CreateFromArray(surfaceGeometry.KnotValueVectorKsi);
                    degree = surfaceGeometry.DegreeKsi;
                    numberOfControlPoints = surfaceGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * surfaceGeometry.NumberOfCpHeta + surfaceGeometry.NumberOfCpHeta - 1;
                        controlPoints.Add(model.ControlPointsDictionary[surfaceGeometry.ControlPointIDs[index]]);
                    }

                    break;
            }

            return (KnotValueVector, degree, controlPoints);
        }

        public static void ConstraintDofsOfEdge(this NurbsSurfaceGeometryDto surfaceGeometry,Model model, NurbsSurfaceEdges edge,
            List<IDofType> constrainedDofs)
        {
            var (_, _, controlPoints) = FindEdgeData2D(surfaceGeometry, model, edge);

            foreach (var controlPoint in controlPoints)
            {
                controlPoint.Constraints.AddRange(constrainedDofs.Select(x=>new Constraint()
                {
                    Amount = 0,
                    DOF = x
                }));
            }
        }

        public static List<ILoad> CreateSurfaceLoad(this NurbsSurfaceGeometryDto surfaceGeometry,Model model, ISurfaceLoad surfaceLoad)
        {
            int degree1 = surfaceGeometry.DegreeKsi;
            int degree2 = surfaceGeometry.DegreeHeta;

            var knotValueVector1 = Vector.CreateFromArray(surfaceGeometry.KnotValueVectorKsi);
            var knotValueVector2 = Vector.CreateFromArray(surfaceGeometry.KnotValueVectorHeta);
            
            #region Knots

			Vector singleKnotValuesKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[0];

			List<Knot> knots = new List<Knot>();

			int id = 0;
			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
			{
				for (int j = 0; j < singleKnotValuesHeta.Length; j++)
				{
					knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
					id++;
				}
			}

			#endregion Knots

			#region Elements

			Vector multiplicityKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[1];
			Vector multiplicityHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
			if (numberOfElementsKsi * numberOfElementsHeta == 0)
			{
				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
			}

            var surfaceLoads = new List<ILoad>();
			for (int i = 0; i < numberOfElementsKsi; i++)
			{
				for (int j = 0; j < numberOfElementsHeta; j++)
				{
					IList<Knot> knotsOfElement = new List<Knot>
					{
						knots[i * singleKnotValuesHeta.Length + j],
						knots[i * singleKnotValuesHeta.Length + j + 1],
						knots[(i + 1) * singleKnotValuesHeta.Length + j],
						knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]
					};

					int multiplicityElementKsi = 0;
					if (multiplicityKsi[i + 1] - degree1 > 0)
					{
						multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degree1;
					}

					int multiplicityElementHeta = 0;
					if (multiplicityHeta[j + 1] - degree2 > 0)
					{
						multiplicityElementHeta = (int)multiplicityHeta[j + 1] - degree2;
					}

					int nurbsSupportKsi = degree1 + 1;
					int nurbsSupportHeta = degree2 + 1;

					int NumberOfControlPointsHeta = knotValueVector2.Length - degree2 - 1;

					IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

					for (int k = 0; k < nurbsSupportKsi; k++)
					{
						for (int l = 0; l < nurbsSupportHeta; l++)
						{
							int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta +
								(j + multiplicityElementHeta) + k * NumberOfControlPointsHeta + l;

							elementControlPoints.Add(model.ControlPointsDictionary[controlPointID]);
						}
					}
					int elementID = i * numberOfElementsHeta + j;

                    var quadrature = new FullQuadrature2D(degree1, degree2, knotsOfElement);
                    var nurbs = new Nurbs2D(degree1, knotValueVector1.CopyToArray(), degree2,
                        knotValueVector2.CopyToArray(), elementControlPoints.ToArray(), quadrature.GaussPoints.ToArray());
                    var surfaceLoadElement = new SurfaceLoadElement(surfaceLoad, nurbs, quadrature, model.ControlPointsDictionary.Values.ToList());
                    surfaceLoads.Add(surfaceLoadElement);
				}
			}

            return surfaceLoads;

            #endregion Elements
        }

        public static List<ILoad> CreateLoadForEdge(this NurbsSolidGeometryDto solidGeometry,Model model, NurbsSolidEdges edge,
            ILineLoad lineLoad)
        {
            var (KnotValueVector, degree, controlPoints) = FindEdgeData3D(solidGeometry, model, edge);
            
            #region Knots

            var singleKnotValuesKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[0];
            var knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = 0.0, Zeta = 0.0 });
                id++;
            }

            #endregion Knots

            #region Elements

            Vector multiplicityKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            if (numberOfElementsKsi == 0)
            {
                throw new ArgumentNullException("Number of Elements should be defined before Element Connectivity");
            }

            var edgeLoads= new List<ILoad>();
            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                IList<Knot> knotsOfElement = new List<Knot>
                {
                    knots[i],
                    knots[i + 1]
                };

                int multiplicityElementKsi = 0;
                if (multiplicityKsi[i + 1] - degree > 0)
                {
                    multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degree;
                }

                int nurbsSupportKsi = degree + 1;

                IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                for (int k = 0; k < nurbsSupportKsi; k++)
                {
                    int controlPointID = i + multiplicityElementKsi + k;
                    elementControlPoints.Add(controlPoints[controlPointID]);
                }

                var quadrature= new FullQuadrature1D(degree, knots[i].Ksi, knots[i+1].Ksi);
                var nurbs = new Nurbs1D(degree, KnotValueVector.CopyToArray(), controlPoints.ToArray(), quadrature.GaussPoints.ToArray());
                var lineLoadElement = new LineLoadElement(lineLoad, nurbs, quadrature, controlPoints);
                edgeLoads.Add(lineLoadElement);
            }

            return edgeLoads;

            #endregion Elements
        }

        private static (Vector KnotValueVector, int degree, List<ControlPoint> controlPoints) FindEdgeData3D(
            NurbsSolidGeometryDto solidGeometry, Model model, NurbsSolidEdges edge)
        {
            Vector KnotValueVector = null;
            int degree = 0;
            int numberOfControlPoints;
            var controlPoints = new List<ControlPoint>();
            switch (edge)
            {
                case NurbsSolidEdges.Edge1:
                    degree = solidGeometry.DegreeZeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpZeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[i]]);
                    break;
                case NurbsSolidEdges.Edge2:
                    degree = solidGeometry.DegreeZeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpZeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i + solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpZeta - 1);
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge3:
                    degree = solidGeometry.DegreeHeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpHeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge4:
                    degree = solidGeometry.DegreeHeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpHeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta + solidGeometry.NumberOfCpZeta - 1;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge5:
                    degree = solidGeometry.DegreeZeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    int offset = solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta *
                                 (solidGeometry.NumberOfCpKsi - 1);
                    numberOfControlPoints = solidGeometry.NumberOfCpZeta;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i + offset;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge6:
                    degree = solidGeometry.DegreeZeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpZeta;
                    offset = solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta *
                             (solidGeometry.NumberOfCpKsi - 1);
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i + solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpHeta - 1) + offset;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge7:
                    degree = solidGeometry.DegreeHeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpHeta;
                    offset = solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta *
                             (solidGeometry.NumberOfCpKsi - 1);
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta + offset;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge8:
                    degree = solidGeometry.DegreeHeta;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    numberOfControlPoints = solidGeometry.NumberOfCpHeta;
                    offset = solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta *
                             (solidGeometry.NumberOfCpKsi - 1);
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta + solidGeometry.NumberOfCpZeta - 1 + offset;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge9:
                    degree = solidGeometry.DegreeKsi;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    numberOfControlPoints = solidGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge10:
                    degree = solidGeometry.DegreeKsi;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    numberOfControlPoints = solidGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta
                            + solidGeometry.NumberOfCpZeta - 1;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge11:
                    degree = solidGeometry.DegreeKsi;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    numberOfControlPoints = solidGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta +
                                    solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpHeta - 1);
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
                case NurbsSolidEdges.Edge12:
                    degree = solidGeometry.DegreeKsi;
                    KnotValueVector = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    numberOfControlPoints = solidGeometry.NumberOfCpKsi;
                    for (var i = 0; i < numberOfControlPoints; i++)
                    {
                        var index = i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta +
                            solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpHeta - 1) +
                            solidGeometry.NumberOfCpZeta - 1;
                        controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                    }

                    break;
            }

            return (KnotValueVector, degree, controlPoints);
        }

        public static void ConstraintDofsOfEdge(this NurbsSolidGeometryDto solidGeometry,Model model, NurbsSolidEdges edge,
            List<IDofType> constrainedDofs)
        {
            var (_, degree, controlPoints) = FindEdgeData3D(solidGeometry, model, edge);

            foreach (var controlPoint in controlPoints)
            {
                controlPoint.Constraints.AddRange(constrainedDofs.Select(x=>new Constraint()
                {
                    Amount = 0,
                    DOF = x
                }));
            }
        }

        public static void ConstraintDofsOfFace(this NurbsSolidGeometryDto solidGeometry,Model model, NurbsSolidFaces face,
            List<IDofType> constrainedDofs)
        {
            var (_, _, _, _, controlPoints) = FindDataOfFace(solidGeometry, model, face);

            foreach (var controlPoint in controlPoints)
            {
                controlPoint.Constraints.AddRange(constrainedDofs.Select(x=>new Constraint()
                {
                    Amount = 0,
                    DOF = x
                }));
            }
        }

        public static List<ILoad> CreateLoadForFace(this NurbsSolidGeometryDto solidGeometry,Model model, NurbsSolidFaces face, ISurfaceLoad surfaceLoad)
        {
            var (degree1, degree2, knotValueVector1, knotValueVector2, controlPoints) = FindDataOfFace(solidGeometry, model, face);
            
            #region Knots

			Vector singleKnotValuesKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[0];

			List<Knot> knots = new List<Knot>();

			int id = 0;
			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
			{
				for (int j = 0; j < singleKnotValuesHeta.Length; j++)
				{
					knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
					id++;
				}
			}

			#endregion Knots

			#region Elements

			Vector multiplicityKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[1];
			Vector multiplicityHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
			if (numberOfElementsKsi * numberOfElementsHeta == 0)
			{
				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
			}

            var surfaceLoads = new List<ILoad>();
			for (int i = 0; i < numberOfElementsKsi; i++)
			{
				for (int j = 0; j < numberOfElementsHeta; j++)
				{
					IList<Knot> knotsOfElement = new List<Knot>
					{
						knots[i * singleKnotValuesHeta.Length + j],
						knots[i * singleKnotValuesHeta.Length + j + 1],
						knots[(i + 1) * singleKnotValuesHeta.Length + j],
						knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]
					};

					int multiplicityElementKsi = 0;
					if (multiplicityKsi[i + 1] - degree1 > 0)
					{
						multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degree1;
					}

					int multiplicityElementHeta = 0;
					if (multiplicityHeta[j + 1] - degree2 > 0)
					{
						multiplicityElementHeta = (int)multiplicityHeta[j + 1] - degree2;
					}

					int nurbsSupportKsi = degree1 + 1;
					int nurbsSupportHeta = degree2 + 1;

					int NumberOfControlPointsHeta = knotValueVector2.Length - degree2 - 1;

					IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

					for (int k = 0; k < nurbsSupportKsi; k++)
					{
						for (int l = 0; l < nurbsSupportHeta; l++)
						{
							int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta +
								(j + multiplicityElementHeta) + k * NumberOfControlPointsHeta + l;

							elementControlPoints.Add(controlPoints[controlPointID]);
						}
					}
					int elementID = i * numberOfElementsHeta + j;

                    var quadrature = new FullQuadrature2D(degree1, degree2, knotsOfElement);
                    var nurbs = new Nurbs2D(degree1, knotValueVector1.CopyToArray(), degree2,
                        knotValueVector2.CopyToArray(), elementControlPoints.ToArray(), quadrature.GaussPoints.ToArray());
                    var surfaceLoadElement = new SurfaceLoadElement(surfaceLoad, nurbs, quadrature, controlPoints);
                    surfaceLoads.Add(surfaceLoadElement);
				}
			}

            return surfaceLoads;

            #endregion Elements
        }

        private static (int degree1, int degree2, Vector knotValueVector1, Vector knotValueVector2, List<ControlPoint> controlPoints)
            FindDataOfFace(NurbsSolidGeometryDto solidGeometry, Model model, NurbsSolidFaces face)
        {
            int degree1 = 0;
            int degree2 = 0;
            Vector knotValueVector1 = null;
            Vector knotValueVector2 = null;
            var controlPoints = new List<ControlPoint>();

            switch (face)
            {
                case NurbsSolidFaces.Left:
                    degree1 = solidGeometry.DegreeHeta;
                    degree2 = solidGeometry.DegreeZeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpHeta; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpZeta; j++)
                        {
                            controlPoints.Add(
                                model.ControlPointsDictionary[
                                    solidGeometry.ControlPointIDs[j + solidGeometry.NumberOfCpZeta * i]]);
                        }
                    }

                    break;
                case NurbsSolidFaces.Right:
                    degree1 = solidGeometry.DegreeHeta;
                    degree2 = solidGeometry.DegreeZeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpHeta; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpZeta; j++)
                        {
                            var index = j + solidGeometry.NumberOfCpZeta * i + solidGeometry.NumberOfCpHeta *
                                solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpKsi - 1);
                            controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                        }
                    }

                    break;
                case NurbsSolidFaces.Bottom:
                    degree1 = solidGeometry.DegreeKsi;
                    degree2 = solidGeometry.DegreeHeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpKsi; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpHeta; j++)
                        {
                            var index = j * solidGeometry.NumberOfCpZeta +
                                        i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta;
                            controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                        }
                    }

                    break;
                case NurbsSolidFaces.Top:
                    degree1 = solidGeometry.DegreeKsi;
                    degree2 = solidGeometry.DegreeHeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpKsi; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpHeta; j++)
                        {
                            var index = solidGeometry.NumberOfCpZeta - 1 + j * solidGeometry.NumberOfCpZeta +
                                        i * solidGeometry.NumberOfCpZeta * solidGeometry.NumberOfCpHeta;
                            controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                        }
                    }

                    break;
                case NurbsSolidFaces.Front:
                    degree1 = solidGeometry.DegreeKsi;
                    degree2 = solidGeometry.DegreeZeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpKsi; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpZeta; j++)
                        {
                            var index = j + i * solidGeometry.NumberOfCpHeta * solidGeometry.NumberOfCpZeta;
                            controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                        }
                    }

                    break;
                case NurbsSolidFaces.Back:
                    degree1 = solidGeometry.DegreeKsi;
                    degree2 = solidGeometry.DegreeZeta;
                    knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
                    knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);
                    for (int i = 0; i < solidGeometry.NumberOfCpKsi; i++)
                    {
                        for (int j = 0; j < solidGeometry.NumberOfCpZeta; j++)
                        {
                            var index = j + i * solidGeometry.NumberOfCpHeta * solidGeometry.NumberOfCpZeta +
                                        solidGeometry.NumberOfCpZeta * (solidGeometry.NumberOfCpHeta - 1);
                            controlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[index]]);
                        }
                    }

                    break;
            }

            return (degree1, degree2, knotValueVector1, knotValueVector2, controlPoints);
        }

        public static List<ILoad> CreateBodyLoad(this NurbsSolidGeometryDto solidGeometry,Model model, IBodyLoad bodyLoad)
        {
            var degree1 = solidGeometry.DegreeKsi;
            var degree2 = solidGeometry.DegreeHeta;
            var degree3 = solidGeometry.DegreeZeta;

            var knotValueVector1 = Vector.CreateFromArray(solidGeometry.KnotValueVectorKsi);
            var knotValueVector2 = Vector.CreateFromArray(solidGeometry.KnotValueVectorHeta);
            var knotValueVector3 = Vector.CreateFromArray(solidGeometry.KnotValueVectorZeta);

            
            #region Knots

            var singleKnotValuesKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[0];
            var singleKnotValuesHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[0];
            var singleKnotValuesZeta = knotValueVector3.RemoveDuplicatesFindMultiplicity()[0];

            var knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
                    {
                        knots.Add(new Knot()
                        {
                            ID = id,
                            Ksi = singleKnotValuesKsi[i],
                            Heta = singleKnotValuesHeta[j],
                            Zeta = singleKnotValuesZeta[k]
                        });
                        id++;
                    }
                }
            }

            #endregion Knots

            #region Elements
            var singlesKnotValuesKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityKsi = knotValueVector1.RemoveDuplicatesFindMultiplicity()[1];
            var singlesKnotValuesHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityHeta = knotValueVector2.RemoveDuplicatesFindMultiplicity()[1];
            var singlesKnotValuesZeta =knotValueVector3.RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityZeta = knotValueVector3.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
            {
                throw new ArgumentException("Number of Elements should be defined before Element Connectivity");
            }

            var bodyLoads= new List<ILoad>();
            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    for (int k = 0; k < numberOfElementsZeta; k++)
                    {
                        IList<Knot> knotsOfElement = new List<Knot>
                        {
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length +
                            k],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length +
                            k + 1],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k],
                            knots[
                            i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k + 1],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            j * singlesKnotValuesZeta.Length + k],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            j * singlesKnotValuesZeta.Length + k + 1],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k],
                            knots[
                            (i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length +
                            (j + 1) * singlesKnotValuesZeta.Length + k + 1]
                        };

                        int multiplicityElementKsi = 0;
                        if (multiplicityKsi[i + 1] - degree1 > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degree1;
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - degree2 > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - degree2;
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - degree3 > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - degree3;
                        }

                        int nurbsSupportKsi = degree1 + 1;
                        int nurbsSupportHeta = degree2 + 1;
                        int nurbsSupportZeta = degree3 + 1;

                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    int controlPointID = (i + multiplicityElementKsi) * solidGeometry.NumberOfCpHeta*
                                                         solidGeometry.NumberOfCpZeta+
                                                         (j + multiplicityElementHeta) * solidGeometry.NumberOfCpZeta+
                                                         (k + multiplicityElementZeta) +
                                                         l * solidGeometry.NumberOfCpHeta* solidGeometry.NumberOfCpZeta+
                                                         m * solidGeometry.NumberOfCpZeta+ n;

                                    elementControlPoints.Add(model.ControlPointsDictionary[solidGeometry.ControlPointIDs[controlPointID]]);
                                }
                            }
                        }

                       
                        var quadrature = new FullQuadrature3D(degree1, degree2, degree3, knotsOfElement);
                        var nurbs = new Nurbs3D(solidGeometry.NumberOfCpKsi, solidGeometry.NumberOfCpHeta,
                            solidGeometry.NumberOfCpZeta,
                            degree1, degree2, degree3, solidGeometry.KnotValueVectorKsi,
                            solidGeometry.KnotValueVectorHeta,
                            solidGeometry.KnotValueVectorZeta, elementControlPoints.ToArray(), quadrature.GaussPoints.ToArray());
                        var bodyLoadElement = new BodyLoadElement(bodyLoad, nurbs, quadrature, elementControlPoints.ToList());
                        bodyLoads.Add(bodyLoadElement);
                    }
                }
            }

            return bodyLoads;

            #endregion
        }
    }
}
