using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.IGA.Elements.Continuum;
using ISAAR.MSolve.IGA.Elements.Structural;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Geometry.NurbsMesh;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interpolation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;
using Newtonsoft.Json;

namespace ISAAR.MSolve.IGA.Readers
{
    public class JsonModelReader
    {
        private readonly string filePath;
        private readonly IContinuumMaterial2D planeMaterial;
        private readonly IContinuumMaterial3D solidMaterial;
        private readonly IShellMaterial shellMaterial;
        private readonly IShellSectionMaterial shellSectionMaterial;

        public JsonModelReader(string filePath, IContinuumMaterial2D planeMaterial=null, IContinuumMaterial3D solidMaterial=null,
            IShellMaterial shellMaterial=null,IShellSectionMaterial shellSectionMaterial=null)
        {
            this.filePath = filePath;
            this.planeMaterial = planeMaterial;
            this.solidMaterial = solidMaterial;
            this.shellMaterial = shellMaterial;
            this.shellSectionMaterial = shellSectionMaterial;
        }

        public (ModelDto geometry, Model model) ReadGeometryAndCreateModel()
        {
            var jsonFile=File.ReadAllText(filePath);
            var geometry = JsonConvert.DeserializeObject<ModelDto>(jsonFile);

            var model=CreateModelData(geometry);

            return (geometry, model);
        }

        private Model CreateModelData(ModelDto geometry)
        {
            var model= new Model();

            foreach (var controlPoint in geometry.ControlPoints)
            {
                model.ControlPointsDictionary.Add(controlPoint.ID, new ControlPoint()
                {
                    ID = controlPoint.ID,
                    X = controlPoint.X,
                    Y = controlPoint.Y,
                    Z = controlPoint.Z,
                    WeightFactor = controlPoint.Weight
                });
            }

            if (geometry.NurbsSolidPatches!=null)
                foreach (var solidPatch in geometry.NurbsSolidPatches)
                    CreateSolidPatchData(solidPatch, model);

            if (geometry.NurbsSurfacePatches!=null)
                foreach (var surfacePatch in geometry.NurbsSurfacePatches)
                    CreateSurfacePatchData(surfacePatch,model);

            return model;
        }

        private void CreateSurfacePatchData(NurbsSurfaceGeometryDto surfacePatch, Model model)
        {
            model.PatchesDictionary.Add(surfacePatch.ID, new Patch());
            model.PatchesDictionary[0].NumberOfDimensions = 3;

            foreach (var controlPointID in surfacePatch.ControlPointIDs)
            {
                model.PatchesDictionary[surfacePatch.ID].ControlPoints.Add(model.ControlPointsDictionary[controlPointID]);
            }

            #region Knots
            var singleKnotValuesKsi = Vector.CreateFromArray(surfacePatch.KnotValueVectorKsi).RemoveDuplicatesFindMultiplicity()[0];
            var singleKnotValuesHeta = Vector.CreateFromArray(surfacePatch.KnotValueVectorHeta).RemoveDuplicatesFindMultiplicity()[0];

            var knots = new List<Knot>();

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

            var multiplicityKsi = Vector.CreateFromArray(surfacePatch.KnotValueVectorKsi).RemoveDuplicatesFindMultiplicity()[1];
            var multiplicityHeta = Vector.CreateFromArray(surfacePatch.KnotValueVectorHeta).RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
            if (numberOfElementsKsi * numberOfElementsHeta == 0)
            {
                throw new ArgumentException("Number of Elements should be defined before Element Connectivity");
            }

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
                    if (multiplicityKsi[i + 1] - surfacePatch.DegreeKsi> 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - surfacePatch.DegreeKsi;
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - surfacePatch.DegreeHeta> 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - surfacePatch.DegreeHeta;
                    }

                    int nurbsSupportKsi = surfacePatch.DegreeKsi + 1;
                    int nurbsSupportHeta = surfacePatch.DegreeHeta+ 1;

                    var elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * surfacePatch.NumberOfCpHeta+
                                                 (j + multiplicityElementHeta) + k * surfacePatch.NumberOfCpHeta+ l;
                            elementControlPoints.Add(model.ControlPointsDictionary[surfacePatch.ControlPointIDs[controlPointID]]);
                        }
                    }

                    int elementID = i * numberOfElementsHeta + j;


                    var gauss = new GaussQuadrature();
                    var gaussPoints = gauss.CalculateElementGaussPoints(surfacePatch.DegreeKsi, 
                        surfacePatch.DegreeHeta, knotsOfElement).ToArray();
                    var nurbs = new Nurbs2D(surfacePatch.DegreeKsi, surfacePatch.KnotValueVectorKsi,
                        surfacePatch.DegreeHeta, surfacePatch.KnotValueVectorHeta,
                        elementControlPoints.ToArray(), gaussPoints);
                    var element = new Element
                    {
                        ID = elementID,
                        Patch = model.PatchesDictionary[surfacePatch.ID],
                    };

                    switch (surfacePatch.GeometryType)
                    {
                        case NurbsGeometryType.Plane:
                            element.ElementType = new ContinuumElement2D(planeMaterial, nurbs, gaussPoints,
                                surfacePatch.Thickness);
                            break;
                        case NurbsGeometryType.ShellLinear:
                            element.ElementType = new KirchhoffLoveShell(shellSectionMaterial, nurbs, gaussPoints,
                                surfacePatch.Thickness);
                            break;
                        case NurbsGeometryType.ShellNonLinear:
                            element.ElementType = new KirchhoffLoveShellNL(shellMaterial,
                                knotsOfElement, nurbs, elementControlPoints, model.PatchesDictionary[0],
                                surfacePatch.Thickness, surfacePatch.DegreeKsi, surfacePatch.DegreeHeta);
                            break;
                        default:
                            throw new NotImplementedException("The nurbs geometry type given in the input file is not implemented");
                    }

                    element.AddKnots(knotsOfElement);
                    element.AddControlPoints(elementControlPoints.ToList());
                    model.PatchesDictionary[surfacePatch.ID].Elements.Add(element);
                    model.ElementsDictionary.Add(elementID, element);
                }
            }

            #endregion Elements
        }
        
        private void CreateSolidPatchData(NurbsSolidGeometryDto solidPatch, Model model)
        {
            model.PatchesDictionary.Add(solidPatch.ID, new Patch());
            model.PatchesDictionary[0].NumberOfDimensions = 3;

            foreach (var controlPointID in solidPatch.ControlPointIDs)
            {
                model.PatchesDictionary[solidPatch.ID].ControlPoints.Add(model.ControlPointsDictionary[controlPointID]);
            }

            #region Knots

            var singleKnotValuesKsi = Vector.CreateFromArray(solidPatch.KnotValueVectorKsi).RemoveDuplicatesFindMultiplicity()[0];
            var singleKnotValuesHeta = Vector.CreateFromArray(solidPatch.KnotValueVectorHeta).RemoveDuplicatesFindMultiplicity()[0];
            var singleKnotValuesZeta = Vector.CreateFromArray(solidPatch.KnotValueVectorZeta).RemoveDuplicatesFindMultiplicity()[0];

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
            var singlesKnotValuesKsi = Vector.CreateFromArray(solidPatch.KnotValueVectorKsi).RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityKsi = Vector.CreateFromArray(solidPatch.KnotValueVectorKsi).RemoveDuplicatesFindMultiplicity()[1];
            var singlesKnotValuesHeta = Vector.CreateFromArray(solidPatch.KnotValueVectorHeta).RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityHeta = Vector.CreateFromArray(solidPatch.KnotValueVectorHeta).RemoveDuplicatesFindMultiplicity()[1];
            var singlesKnotValuesZeta = Vector.CreateFromArray(solidPatch.KnotValueVectorZeta).RemoveDuplicatesFindMultiplicity()[0];
            var multiplicityZeta = Vector.CreateFromArray(solidPatch.KnotValueVectorZeta).RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
            {
                throw new ArgumentException("Number of Elements should be defined before Element Connectivity");
            }

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
                        if (multiplicityKsi[i + 1] - solidPatch.DegreeKsi > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - solidPatch.DegreeKsi;
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - solidPatch.DegreeHeta > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - solidPatch.DegreeHeta;
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - solidPatch.DegreeZeta > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - solidPatch.DegreeZeta;
                        }

                        int nurbsSupportKsi = solidPatch.DegreeKsi + 1;
                        int nurbsSupportHeta = solidPatch.DegreeHeta + 1;
                        int nurbsSupportZeta = solidPatch.DegreeZeta + 1;

                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    int controlPointID = (i + multiplicityElementKsi) * solidPatch.NumberOfCpHeta*
                                                         solidPatch.NumberOfCpZeta+
                                                         (j + multiplicityElementHeta) * solidPatch.NumberOfCpZeta+
                                                         (k + multiplicityElementZeta) +
                                                         l * solidPatch.NumberOfCpHeta* solidPatch.NumberOfCpZeta+
                                                         m * solidPatch.NumberOfCpZeta+ n;

                                    elementControlPoints.Add(model.ControlPointsDictionary[solidPatch.ControlPointIDs[controlPointID]]);
                                }
                            }
                        }

                        int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
                        var gauss = new GaussQuadrature();
                        var gaussPoints = gauss.CalculateElementGaussPoints(solidPatch.DegreeKsi,
                            solidPatch.DegreeHeta, solidPatch.DegreeZeta, knotsOfElement);
                        var nurbs= new Nurbs3D(solidPatch.NumberOfCpKsi, 
                            solidPatch.NumberOfCpHeta,solidPatch.NumberOfCpZeta,
                            solidPatch.DegreeKsi, solidPatch.DegreeHeta,solidPatch.DegreeZeta,
                            solidPatch.KnotValueVectorKsi, solidPatch.KnotValueVectorHeta,
                            solidPatch.KnotValueVectorZeta, elementControlPoints.ToArray(), gaussPoints);

                        var element = new Element
                        {
                            ID = elementID,
                            Patch = model.PatchesDictionary[solidPatch.ID],
                            ElementType = new ContinuumElement3D(solidMaterial, nurbs, gaussPoints)
                        };
                        element.AddKnots(knotsOfElement);
                        element.AddControlPoints(elementControlPoints.ToList());
                        model.PatchesDictionary[solidPatch.ID].Elements.Add(element);
                        model.ElementsDictionary.Add(elementID, element);
                    }
                }
            }
            

            #endregion
        }

    }
}
