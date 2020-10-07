using System;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.IGA.SupportiveClasses.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Elements.Structural
{
    public struct a3r
    {
        public double a3r00;
        public double a3r01;
        public double a3r02;

        public double a3r10;
        public double a3r11;
        public double a3r12;

        public double a3r20;
        public double a3r21;
        public double a3r22;
    }

    public struct a3rs
    {
        public double a3rs00_0;
        public double a3rs00_1;
        public double a3rs00_2;

        public double a3rs01_0;
        public double a3rs01_1;
        public double a3rs01_2;

        public double a3rs02_0;
        public double a3rs02_1;
        public double a3rs02_2;

        public double a3rs10_0;
        public double a3rs10_1;
        public double a3rs10_2;

        public double a3rs11_0;
        public double a3rs11_1;
        public double a3rs11_2;

        public double a3rs12_0;
        public double a3rs12_1;
        public double a3rs12_2;

        public double a3rs20_0;
        public double a3rs20_1;
        public double a3rs20_2;

        public double a3rs21_0;
        public double a3rs21_1;
        public double a3rs21_2;

        public double a3rs22_0;
        public double a3rs22_1;
        public double a3rs22_2;
    }

    public struct Bab_rs
    {
        public double Bab_rs00_0;
        public double Bab_rs00_1;
        public double Bab_rs00_2;

        public double Bab_rs01_0;
        public double Bab_rs01_1;
        public double Bab_rs01_2;

        public double Bab_rs02_0;
        public double Bab_rs02_1;
        public double Bab_rs02_2;

        public double Bab_rs10_0;
        public double Bab_rs10_1;
        public double Bab_rs10_2;

        public double Bab_rs11_0;
        public double Bab_rs11_1;
        public double Bab_rs11_2;

        public double Bab_rs12_0;
        public double Bab_rs12_1;
        public double Bab_rs12_2;

        public double Bab_rs20_0;
        public double Bab_rs20_1;
        public double Bab_rs20_2;

        public double Bab_rs21_0;
        public double Bab_rs21_1;
        public double Bab_rs21_2;

        public double Bab_rs22_0;
        public double Bab_rs22_1;
        public double Bab_rs22_2;
    }

    public struct Forces
    {
        public double v0;
        public double v1;
        public double v2;
    }

    public class KirchhoffLoveShellNL : Element, IStructuralIsogeometricElement
    {
        internal IShapeFunction2D _shapeFunctions;

        public double[] _solution;

        internal Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
            materialsAtThicknessGP =
                new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();

        protected static readonly IDofType[] ControlPointDofTypes = { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };

        private const int ThicknessIntegrationDegree = 2;

        private readonly ControlPoint[] _controlPoints;

        private readonly int _degreeHeta;

        private readonly int _degreeKsi;

        private readonly double[,] jacobianMatrix = new double[2, 3];

        private a3rs a3rs = new a3rs();

        private Bab_rs Bab_rs = new Bab_rs();

        private IDofType[][] dofTypes;

        private double[] InitialJ1;

        private double[][] initialSurfaceBasisVectorDerivative1;

        private double[][] initialSurfaceBasisVectorDerivative12;

        private double[][] initialSurfaceBasisVectorDerivative2;

        private double[][] initialSurfaceBasisVectors1;

        private double[][] initialSurfaceBasisVectors2;

        private double[][] initialUnitSurfaceBasisVectors3;

        private bool isInitialized;

        private Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints =
                    new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();

        public KirchhoffLoveShellNL(IShellMaterial shellMaterial, IList<Knot> elementKnots,
            IShapeFunction2D shapeFunctions, IList<ControlPoint> elementControlPoints, Patch patch, double thickness,
            int degreeKsi, int degreeHeta)
        {
            this.Patch = patch;
            this.Thickness = thickness;
            _degreeKsi = degreeKsi;
            _degreeHeta = degreeHeta;
            foreach (var knot in elementKnots)
            {
                if (!KnotsDictionary.ContainsKey(knot.ID))
                    this.KnotsDictionary.Add(knot.ID, knot);
            }

            _shapeFunctions = shapeFunctions;
            _solution = new double[3 * elementControlPoints.Count];

            CreateElementGaussPoints(this);
            foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
            {
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
                foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
                {
                    materialsAtThicknessGP[medianSurfaceGP].Add(point, shellMaterial.Clone());
                }
            }

            _controlPoints = elementControlPoints.ToArray();
        }

        public KirchhoffLoveShellNL(List<IShellMaterial> shellMaterials, IList<Knot> elementKnots,
            IShapeFunction2D shapeFunctions, IList<ControlPoint> elementControlPoints, Patch patch, double thickness,
            int degreeKsi, int degreeHeta)
        {
            this.Patch = patch;
            this.Thickness = thickness;
            _degreeKsi = degreeKsi;
            _degreeHeta = degreeHeta;
            foreach (var knot in elementKnots)
            {
                if (!KnotsDictionary.ContainsKey(knot.ID))
                    this.KnotsDictionary.Add(knot.ID, knot);
            }

            _shapeFunctions = shapeFunctions;
            _solution = new double[3 * elementControlPoints.Count];

            CreateElementGaussPoints(this);
            foreach (var medianSurfaceGP in thicknessIntegrationPoints.Keys)
            {
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IShellMaterial>());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][0], shellMaterials[0].Clone());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][1], shellMaterials[1].Clone());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][2], shellMaterials[2].Clone());
            }

            _controlPoints = elementControlPoints.ToArray();
        }

        public CellType CellType { get; } = CellType.Unknown;

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public bool MaterialModified => false;

        public double Thickness { get; }

        //public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) =>
        //    throw new NotImplementedException();

        public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
        {
            var knots = element.Knots.ToList();
            var localCoordinates = new double[4, 2]
            {
                {knots[0].Ksi, knots[0].Heta},
                {knots[1].Ksi, knots[1].Heta},
                {knots[2].Ksi, knots[2].Heta},
                {knots[3].Ksi, knots[3].Heta}
            };
            
            var shapeFunctionsK0 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[0, 0], localCoordinates[0, 1]));
            var shapeFunctionsK1 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[1, 0], localCoordinates[1, 1]));
            var shapeFunctionsK2 =_shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[2, 0], localCoordinates[2, 1]));
            var shapeFunctionsK3 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[3, 0], localCoordinates[3, 1]));

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };

            for (int i = 0; i < element.ControlPoints.Count(); i++)
            {
                knotDisplacements[paraviewKnotRenumbering[0], 0] += shapeFunctionsK0[i,0] * localDisplacements[i, 0];
                knotDisplacements[paraviewKnotRenumbering[0], 1] += shapeFunctionsK0[i,0] * localDisplacements[i, 1];
                knotDisplacements[paraviewKnotRenumbering[0], 2] += shapeFunctionsK0[i,0] * localDisplacements[i, 2];
                                                                                      
                knotDisplacements[paraviewKnotRenumbering[1], 0] += shapeFunctionsK1[i,0] * localDisplacements[i, 0];
                knotDisplacements[paraviewKnotRenumbering[1], 1] += shapeFunctionsK1[i,0] * localDisplacements[i, 1];
                knotDisplacements[paraviewKnotRenumbering[1], 2] += shapeFunctionsK1[i,0] * localDisplacements[i, 2];
                                                                                      
                knotDisplacements[paraviewKnotRenumbering[2], 0] += shapeFunctionsK2[i,0] * localDisplacements[i, 0];
                knotDisplacements[paraviewKnotRenumbering[2], 1] += shapeFunctionsK2[i,0] * localDisplacements[i, 1];
                knotDisplacements[paraviewKnotRenumbering[2], 2] += shapeFunctionsK2[i,0] * localDisplacements[i, 2];

                knotDisplacements[paraviewKnotRenumbering[3], 0] += shapeFunctionsK3[i,0] * localDisplacements[i, 0];
                knotDisplacements[paraviewKnotRenumbering[3], 1] += shapeFunctionsK3[i,0] * localDisplacements[i, 1];
                knotDisplacements[paraviewKnotRenumbering[3], 2] += shapeFunctionsK3[i,0] * localDisplacements[i, 2];
            }

            return knotDisplacements;
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            var shellElement = (KirchhoffLoveShellNL)element;
            var elementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(_controlPoints);
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            var Bmembrane = new double[3, _controlPoints.Length * 3];
            var Bbending = new double[3, _controlPoints.Length * 3];
            var numberOfControlPoints = _controlPoints.Length;
            var MembraneForces = new Forces();
            var BendingMoments = new Forces();

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                CalculateJacobian(newControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, _shapeFunctions, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;


                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector(hessianMatrix, 2);

                CalculateMembraneDeformationMatrix(numberOfControlPoints, _shapeFunctions, j, surfaceBasisVector1,
                    surfaceBasisVector2, Bmembrane);
                CalculateBendingDeformationMatrix(numberOfControlPoints, surfaceBasisVector3, _shapeFunctions, j,
                    surfaceBasisVector2,
                    surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, Bbending);


                // Bbending = Matrix.CreateFromArray(Bbending).Scale(-1).CopytoArray2D();
                IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

                var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                for (int i = 0; i < Bmembrane.GetLength(1); i++)
                {
                    elementNodalForces[i] +=
                        (Bmembrane[0, i] * MembraneForces.v0 * wfactor + Bbending[0, i] * BendingMoments.v0 * wfactor) +
                        (Bmembrane[1, i] * MembraneForces.v1 * wfactor + Bbending[1, i] * BendingMoments.v1 * wfactor) +
                        (Bmembrane[2, i] * MembraneForces.v2 * wfactor + Bbending[2, i] * BendingMoments.v2 * wfactor);
                }
            }
            return elementNodalForces;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[,] CalculatePointsForPostProcessing(Element element)
        {
            var knots = element.Knots.ToList();
            var localCoordinates = new double[4, 2]
            {
                {knots[0].Ksi, knots[0].Heta},
                {knots[1].Ksi, knots[1].Heta},
                {knots[2].Ksi, knots[2].Heta},
                {knots[3].Ksi, knots[3].Heta}
            };

            
            var shapeFunctionsK0 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[0, 0], localCoordinates[0, 1]));
            var shapeFunctionsK1 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[1, 0], localCoordinates[1, 1]));
            var shapeFunctionsK2 =_shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[2, 0], localCoordinates[2, 1]));
            var shapeFunctionsK3 = _shapeFunctions.EvaluateFunctionsAt(new NaturalPoint(localCoordinates[3, 0], localCoordinates[3, 1]));

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] { 0, 3, 1, 2 };
            var elementControlPoints = element.ControlPoints.ToArray();
            for (int i = 0; i < elementControlPoints.Length; i++)
            {
                knotDisplacements[paraviewKnotRenumbering[0], 0] += shapeFunctionsK0[i,0] * elementControlPoints[i].X;
                knotDisplacements[paraviewKnotRenumbering[0], 1] += shapeFunctionsK0[i,0] * elementControlPoints[i].Y;
                knotDisplacements[paraviewKnotRenumbering[0], 2] += shapeFunctionsK0[i,0] * elementControlPoints[i].Z;
                                                                                      
                knotDisplacements[paraviewKnotRenumbering[1], 0] += shapeFunctionsK1[i,0] * elementControlPoints[i].X;
                knotDisplacements[paraviewKnotRenumbering[1], 1] += shapeFunctionsK1[i,0] * elementControlPoints[i].Y;
                knotDisplacements[paraviewKnotRenumbering[1], 2] += shapeFunctionsK1[i,0] * elementControlPoints[i].Z;
                                                                                      
                knotDisplacements[paraviewKnotRenumbering[2], 0] += shapeFunctionsK2[i,0] * elementControlPoints[i].X;
                knotDisplacements[paraviewKnotRenumbering[2], 1] += shapeFunctionsK2[i,0] * elementControlPoints[i].Y;
                knotDisplacements[paraviewKnotRenumbering[2], 2] += shapeFunctionsK2[i,0] * elementControlPoints[i].Z;
                                                                                      
                knotDisplacements[paraviewKnotRenumbering[3], 0] += shapeFunctionsK3[i,0] * elementControlPoints[i].X;
                knotDisplacements[paraviewKnotRenumbering[3], 1] += shapeFunctionsK3[i,0] * elementControlPoints[i].Y;
                knotDisplacements[paraviewKnotRenumbering[3], 2] += shapeFunctionsK3[i,0] * elementControlPoints[i].Z;
            }

            return knotDisplacements;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            var elementControlPoints = ((KirchhoffLoveShellNL)element).ControlPoints.ToArray();


            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(elementControlPoints);
            var midsurfaceGP = materialsAtThicknessGP.Keys.ToArray();

            for (var j = 0; j < midsurfaceGP.Length; j++)
            {
                CalculateJacobian(newControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, _shapeFunctions, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0]
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;

                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector(hessianMatrix, 2);

                var A11 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors1[j][0] +
                          initialSurfaceBasisVectors1[j][1] * initialSurfaceBasisVectors1[j][1] +
                          initialSurfaceBasisVectors1[j][2] * initialSurfaceBasisVectors1[j][2];

                var A22 = initialSurfaceBasisVectors2[j][0] * initialSurfaceBasisVectors2[j][0] +
                          initialSurfaceBasisVectors2[j][1] * initialSurfaceBasisVectors2[j][1] +
                          initialSurfaceBasisVectors2[j][2] * initialSurfaceBasisVectors2[j][2];

                var A12 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors2[j][0] +
                          initialSurfaceBasisVectors1[j][1] * initialSurfaceBasisVectors2[j][1] +
                          initialSurfaceBasisVectors1[j][2] * initialSurfaceBasisVectors2[j][2];

                var a11 = surfaceBasisVector1[0] * surfaceBasisVector1[0] +
                          surfaceBasisVector1[1] * surfaceBasisVector1[1] +
                          surfaceBasisVector1[2] * surfaceBasisVector1[2];

                var a22 = surfaceBasisVector2[0] * surfaceBasisVector2[0] +
                          surfaceBasisVector2[1] * surfaceBasisVector2[1] +
                          surfaceBasisVector2[2] * surfaceBasisVector2[2];

                var a12 = surfaceBasisVector1[0] * surfaceBasisVector2[0] +
                          surfaceBasisVector1[1] * surfaceBasisVector2[1] +
                          surfaceBasisVector1[2] * surfaceBasisVector2[2];


                var membraneStrain = new double[] { 0.5 * (a11 - A11), 0.5 * (a22 - A22), a12 - A12 };


                var B11 = initialSurfaceBasisVectorDerivative1[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
                          initialSurfaceBasisVectorDerivative1[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
                          initialSurfaceBasisVectorDerivative1[j][2] * initialUnitSurfaceBasisVectors3[j][2];

                var B22 = initialSurfaceBasisVectorDerivative2[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
                          initialSurfaceBasisVectorDerivative2[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
                          initialSurfaceBasisVectorDerivative2[j][2] * initialUnitSurfaceBasisVectors3[j][2];

                var B12 = initialSurfaceBasisVectorDerivative12[j][0] * initialUnitSurfaceBasisVectors3[j][0] +
                          initialSurfaceBasisVectorDerivative12[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
                          initialSurfaceBasisVectorDerivative12[j][2] * initialUnitSurfaceBasisVectors3[j][2];

                var b11 = surfaceBasisVectorDerivative1[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative1[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative1[2] * surfaceBasisVector3[2];

                var b22 = surfaceBasisVectorDerivative2[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative2[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative2[2] * surfaceBasisVector3[2];

                var b12 = surfaceBasisVectorDerivative12[0] * surfaceBasisVector3[0] +
                          surfaceBasisVectorDerivative12[1] * surfaceBasisVector3[1] +
                          surfaceBasisVectorDerivative12[2] * surfaceBasisVector3[2];

                var bendingStrain = new double[] { b11 - B11, b22 - B22, 2 * b12 - 2 * B12 };

                foreach (var keyValuePair in materialsAtThicknessGP[midsurfaceGP[j]])
                {
                    var thicknessPoint = keyValuePair.Key;
                    var material = keyValuePair.Value;
                    var gpStrain = new double[bendingStrain.Length];
                    var z = thicknessPoint.Zeta;
                    for (var i = 0; i < bendingStrain.Length; i++)
                    {
                        gpStrain[i] += membraneStrain[i] + bendingStrain[i] * z;
                    }

                    material.UpdateMaterial(gpStrain);
                }
            }

            return new Tuple<double[], double[]>(new double[0], new double[0]);
        }
        
        public void ClearMaterialState()
        {
        }

        public void ClearMaterialStresses() => throw new NotImplementedException();

        public IMatrix DampingMatrix(IElement element) => throw new NotImplementedException();

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
        {
            dofTypes = new IDofType[element.Nodes.Count][];
            for (var i = 0; i < element.Nodes.Count; i++)
            {
                dofTypes[i] = ControlPointDofTypes;
            }

            return dofTypes;
        }

        public IMatrix MassMatrix(IElement element) => throw new NotImplementedException();

        public void ResetMaterialModified() => throw new NotImplementedException();

        public void SaveMaterialState()
        {
            foreach (var gp in materialsAtThicknessGP.Keys)
            {
                foreach (var material in materialsAtThicknessGP[gp].Values)
                {
                    material.SaveState();
                }
            }
        }
        public IMatrix StiffnessMatrix(IElement element)
        {
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();
            if (!isInitialized)
            {
                CalculateInitialConfigurationData(_controlPoints, _shapeFunctions, gaussPoints);
                isInitialized = true;
            }

            var elementControlPoints = CurrentControlPoint(_controlPoints);

            var bRows = 3;
            var bCols = elementControlPoints.Length * 3;
            var stiffnessMatrix = new double[bCols, bCols];
            var KmembraneL = new double[bCols, bCols];
            var KbendingL = new double[bCols, bCols];
            var KmembraneNL = new double[bCols, bCols];
            var KbendingNL = new double[bCols, bCols];
            var Bmembrane = new double[bRows, bCols];
            var Bbending = new double[bRows, bCols];

            var BmTranspose = new double[bCols, bRows];
            var BbTranspose = new double[bCols, bRows];

            var BmTransposeMultStiffness = new double[bCols, bRows];
            var BbTransposeMultStiffness = new double[bCols, bRows];
            var BmbTransposeMultStiffness = new double[bCols, bRows];
            var BbmTransposeMultStiffness = new double[bCols, bRows];
            var MembraneForces = new Forces();
            var BendingMoments = new Forces();

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                CalculateJacobian(elementControlPoints, _shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(elementControlPoints, _shapeFunctions, j);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0] * surfaceBasisVector3[0] +
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1] +
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;

                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector(hessianMatrix, 2);

                double wFactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                CalculateLinearStiffness(elementControlPoints, _shapeFunctions, j, surfaceBasisVector1, surfaceBasisVector2,
                    Bmembrane, surfaceBasisVector3, surfaceBasisVectorDerivative1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, Bbending, gaussPoints, BmTranspose, bRows, bCols, BbTranspose,
                    wFactor, BmTransposeMultStiffness, BbTransposeMultStiffness, BmbTransposeMultStiffness,
                    BbmTransposeMultStiffness, stiffnessMatrix, KmembraneL, KbendingL);
                CalculateNonLinearStiffness(gaussPoints, j, KmembraneNL, bCols, KbendingNL, elementControlPoints, _shapeFunctions,
                    surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, surfaceBasisVectorDerivative1,
                    surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, J1, stiffnessMatrix, wFactor,
                    ref MembraneForces, ref BendingMoments);

            }


            return Matrix.CreateFromArray(stiffnessMatrix);
        }

        internal static (a3rs, double[,][], double[][], double[][], double[], double[], double[,], double[], double[,][]) Calculate_a3rs(Vector surfaceBasisVector1, Vector surfaceBasisVector2,
            Vector surfaceBasisVector3, double J1,
            double dKsi_r, double dKsi_s, double dHeta_r, double dHeta_s)
        {
            var da3_drds = new double[3, 3][];
            var da3tilde_drds = new double[3, 3][];
            var da3tilde_dr = new double[3][];
            var da3tilde_ds = new double[3][];
            double[] a3_tilde;
            double[] dnorma3_dr = new double[3];
            double[] dnorma3_ds = new double[3];
            double[,] dnorma3_drds = new double[3, 3];

            //5.30
            Calculate_da3tilde_drds(dKsi_r, dKsi_s, dHeta_r, dHeta_s, da3tilde_drds);

            //5.24
            Calculate_da3tilde_dr(surfaceBasisVector1, surfaceBasisVector2, dKsi_r, dHeta_r, da3tilde_dr);

            //5.24
            Calculate_da3tilde_dr(surfaceBasisVector1, surfaceBasisVector2, dKsi_s, dHeta_s, da3tilde_ds);


            //5.25
            a3_tilde = CalculateTerm525(surfaceBasisVector3, J1, dnorma3_dr, da3tilde_dr, dnorma3_ds, da3tilde_ds);

            //5.31
            CalculateTerm531(J1, da3tilde_drds, a3_tilde, da3tilde_dr, da3tilde_ds, dnorma3_drds);

            //5.32
            CalculateTerm532(J1, da3tilde_drds, dnorma3_ds, da3tilde_dr, dnorma3_dr, da3tilde_ds, dnorma3_drds, a3_tilde, da3_drds);

            a3rs a3rsAlternative = new a3rs();
            a3rsAlternative.a3rs00_0 = da3_drds[0, 0][0]; a3rsAlternative.a3rs00_1 = da3_drds[0, 0][1]; a3rsAlternative.a3rs00_2 = da3_drds[0, 0][2];
            a3rsAlternative.a3rs01_0 = da3_drds[0, 1][0]; a3rsAlternative.a3rs01_1 = da3_drds[0, 1][1]; a3rsAlternative.a3rs01_2 = da3_drds[0, 1][2];
            a3rsAlternative.a3rs02_0 = da3_drds[0, 2][0]; a3rsAlternative.a3rs02_1 = da3_drds[0, 2][1]; a3rsAlternative.a3rs02_2 = da3_drds[0, 2][2];

            a3rsAlternative.a3rs10_0 = da3_drds[1, 0][0]; a3rsAlternative.a3rs10_1 = da3_drds[1, 0][1]; a3rsAlternative.a3rs10_2 = da3_drds[1, 0][2];
            a3rsAlternative.a3rs11_0 = da3_drds[1, 1][0]; a3rsAlternative.a3rs11_1 = da3_drds[1, 1][1]; a3rsAlternative.a3rs11_2 = da3_drds[1, 1][2];
            a3rsAlternative.a3rs12_0 = da3_drds[1, 2][0]; a3rsAlternative.a3rs12_1 = da3_drds[1, 2][1]; a3rsAlternative.a3rs12_2 = da3_drds[1, 2][2];

            a3rsAlternative.a3rs20_0 = da3_drds[2, 0][0]; a3rsAlternative.a3rs20_1 = da3_drds[2, 0][1]; a3rsAlternative.a3rs20_2 = da3_drds[2, 0][2];
            a3rsAlternative.a3rs21_0 = da3_drds[2, 1][0]; a3rsAlternative.a3rs21_1 = da3_drds[2, 1][1]; a3rsAlternative.a3rs21_2 = da3_drds[2, 1][2];
            a3rsAlternative.a3rs22_0 = da3_drds[2, 2][0]; a3rsAlternative.a3rs22_1 = da3_drds[2, 2][1]; a3rsAlternative.a3rs22_2 = da3_drds[2, 2][2];

            return (a3rsAlternative, da3tilde_drds, da3tilde_dr, da3tilde_ds, dnorma3_dr, dnorma3_ds, dnorma3_drds, a3_tilde, da3_drds);

        }

        internal static void CalculateA3r(double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double dksi_r, double dheta_r, double J1, ref a3r da3_unit_dr_out)
        {
            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var da3tilde_dr = new double[3][];
            //5.24
            var a1r = new double[3];
            var a2r = new double[3];
            var sum1 = new double[3];
            var sum2 = new double[3];

            // da3_tilde_dr0     ...  da3tilde_dr[0][...]; r1=0
            a1r[0] = dksi_r;
            a2r[0] = dheta_r;
            CalculateCrossProduct(a1r, surfaceBasisVector2, sum1);
            CalculateCrossProduct(surfaceBasisVector1, a2r, sum2);

            sum2[0] += sum1[0];
            sum2[1] += sum1[1];
            sum2[2] += sum1[2];

            var da3_tilde_dr00 = sum2[0];
            var da3_tilde_dr10 = sum2[1];
            var da3_tilde_dr20 = sum2[2];

            // da3_tilde_dr1     ...  da3tilde_dr[1][...]; r1=1
            a1r[0] = 0.0;
            a2r[0] = 0.0;
            a1r[1] = dksi_r;
            a2r[1] = dheta_r;
            sum1[0] = 0.0; sum1[1] = 0.0; sum1[2] = 0.0;
            sum2[0] = 0.0; sum2[1] = 0.0; sum2[2] = 0.0;
            CalculateCrossProduct(a1r, surfaceBasisVector2, sum1);
            CalculateCrossProduct(surfaceBasisVector1, a2r, sum2);

            sum2[0] += sum1[0];
            sum2[1] += sum1[1];
            sum2[2] += sum1[2];

            var da3_tilde_dr01 = sum2[0];
            var da3_tilde_dr11 = sum2[1];
            var da3_tilde_dr21 = sum2[2];


            // da3_tilde_dr2     ...  da3tilde_dr[2][...]; r1=2
            a1r[1] = 0.0;
            a2r[1] = 0.0;
            a1r[2] = dksi_r;
            a2r[2] = dheta_r;
            sum1[0] = 0.0; sum1[1] = 0.0; sum1[2] = 0.0;
            sum2[0] = 0.0; sum2[1] = 0.0; sum2[2] = 0.0;
            CalculateCrossProduct(a1r, surfaceBasisVector2, sum1);
            CalculateCrossProduct(surfaceBasisVector1, a2r, sum2);

            sum2[0] += sum1[0];
            sum2[1] += sum1[1];
            sum2[2] += sum1[2];

            var da3_tilde_dr02 = sum2[0];
            var da3_tilde_dr12 = sum2[1];
            var da3_tilde_dr22 = sum2[2];


            var dnorma3_dr0 = s30 * da3_tilde_dr00 +
                              s31 * da3_tilde_dr10 +
                              s32 * da3_tilde_dr20;

            var dnorma3_dr1 = s30 * da3_tilde_dr01 +
                              s31 * da3_tilde_dr11 +
                              s32 * da3_tilde_dr21;

            var dnorma3_dr2 = s30 * da3_tilde_dr02 +
                              s31 * da3_tilde_dr12 +
                              s32 * da3_tilde_dr22;

            da3_unit_dr_out.a3r00 = (da3_tilde_dr00 - s30 * dnorma3_dr0) / J1;
            da3_unit_dr_out.a3r10 = (da3_tilde_dr10 - s31 * dnorma3_dr0) / J1;
            da3_unit_dr_out.a3r20 = (da3_tilde_dr20 - s32 * dnorma3_dr0) / J1;

            da3_unit_dr_out.a3r01 = (da3_tilde_dr01 - s30 * dnorma3_dr1) / J1;
            da3_unit_dr_out.a3r11 = (da3_tilde_dr11 - s31 * dnorma3_dr1) / J1;
            da3_unit_dr_out.a3r21 = (da3_tilde_dr21 - s32 * dnorma3_dr1) / J1;

            da3_unit_dr_out.a3r02 = (da3_tilde_dr02 - s30 * dnorma3_dr2) / J1;
            da3_unit_dr_out.a3r12 = (da3_tilde_dr12 - s31 * dnorma3_dr2) / J1;
            da3_unit_dr_out.a3r22 = (da3_tilde_dr22 - s32 * dnorma3_dr2) / J1;
        }

        internal static Bab_rs CalculateBab_rs(double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12,
            double ddksi_r, double ddksi_s, double ddheta_r, double ddheta_s, double dksidheta_r,
            double dksidheta_s, a3rs a3rsAlternative, a3r a3r, a3r a3s, double[,][] da3_drds)
        {
            var s11 = surfaceBasisVectorDerivative1;
            var s22 = surfaceBasisVectorDerivative2;
            var s12 = surfaceBasisVectorDerivative12;

            // r= 0 sunistwses _0, _1, kai _2
            var a3rVecs_r0_0 = a3r.a3r00;
            var a3rVecs_r0_1 = a3r.a3r10;
            var a3rVecs_r0_2 = a3r.a3r20;

            // r= 1 sunistwses _0, _1, kai _2
            var a3rVecs_r1_0 = a3r.a3r01;
            var a3rVecs_r1_1 = a3r.a3r11;
            var a3rVecs_r1_2 = a3r.a3r21;

            // r= 2 sunistwses _0, _1, kai _2
            var a3rVecs_r2_0 = a3r.a3r02;
            var a3rVecs_r2_1 = a3r.a3r12;
            var a3rVecs_r2_2 = a3r.a3r22;

            // s= 0 sunistwses _0, _1, kai _2
            var a3sVecs_s0_0 = a3s.a3r00;
            var a3sVecs_s0_1 = a3s.a3r10;
            var a3sVecs_s0_2 = a3s.a3r20;

            var a3sVecs_s1_0 = a3s.a3r01;
            var a3sVecs_s1_1 = a3s.a3r11;
            var a3sVecs_s1_2 = a3s.a3r21;

            var a3sVecs_s2_0 = a3s.a3r02;
            var a3sVecs_s2_1 = a3s.a3r12;
            var a3sVecs_s2_2 = a3s.a3r22;

            var Bab_rsAlternative = new Bab_rs();

            //...............................[r, s].............[r, s]
            var da3 = new double[3] { da3_drds[0, 0][0], da3_drds[0, 0][1], da3_drds[0, 0][2] };
            //[A]: 0--> b11rs , a11r ,a11s surfaceBasisVectorDerivative1,   
            //     1--> b22rs , a22r ,a22s surfaceBasisVectorDerivative2, 
            //     2--> b12rs , a12r ,a12s surfaceBasisVectorDerivative12, times x2
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s..........................................
            Bab_rsAlternative.Bab_rs00_0 = ddksi_r * a3sVecs_s0_0 + ddksi_s * a3rVecs_r0_0 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs00_1 = ddheta_r * a3sVecs_s0_0 + ddheta_s * a3rVecs_r0_0 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs00_2 = dksidheta_r * a3sVecs_s0_0 + dksidheta_s * a3rVecs_r0_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s0_0 + dksidheta_s * a3rVecs_r0_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //...............[r, s]
            da3[0] = da3_drds[0, 1][0]; da3[1] = da3_drds[0, 1][1]; da3[2] = da3_drds[0, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s.......................................................................[r,s]
            Bab_rsAlternative.Bab_rs01_0 = ddksi_r * a3sVecs_s1_0 + ddksi_s * a3rVecs_r0_1 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs01_1 = ddheta_r * a3sVecs_s1_0 + ddheta_s * a3rVecs_r0_1 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs01_2 = dksidheta_r * a3sVecs_s1_0 + dksidheta_s * a3rVecs_r0_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s1_0 + dksidheta_s * a3rVecs_r0_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //...............[r, s]
            da3[0] = da3_drds[0, 2][0]; da3[1] = da3_drds[0, 2][1]; da3[2] = da3_drds[0, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s]..........r_s............................................................................[r,s]
            Bab_rsAlternative.Bab_rs02_0 = ddksi_r * a3sVecs_s2_0 + ddksi_s * a3rVecs_r0_2 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs02_1 = ddheta_r * a3sVecs_s2_0 + ddheta_s * a3rVecs_r0_2 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs02_2 = dksidheta_r * a3sVecs_s2_0 + dksidheta_s * a3rVecs_r0_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s2_0 + dksidheta_s * a3rVecs_r0_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 0][0]; da3[1] = da3_drds[1, 0][1]; da3[2] = da3_drds[1, 0][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s..........................................................................[r,s]
            Bab_rsAlternative.Bab_rs10_0 = ddksi_r * a3sVecs_s0_1 + ddksi_s * a3rVecs_r1_0 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs10_1 = ddheta_r * a3sVecs_s0_1 + ddheta_s * a3rVecs_r1_0 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs10_2 = dksidheta_r * a3sVecs_s0_1 + dksidheta_s * a3rVecs_r1_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s0_1 + dksidheta_s * a3rVecs_r1_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 1][0]; da3[1] = da3_drds[1, 1][1]; da3[2] = da3_drds[1, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]
            Bab_rsAlternative.Bab_rs11_0 = ddksi_r * a3sVecs_s1_1 + ddksi_s * a3rVecs_r1_1 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs11_1 = ddheta_r * a3sVecs_s1_1 + ddheta_s * a3rVecs_r1_1 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs11_2 = dksidheta_r * a3sVecs_s1_1 + dksidheta_s * a3rVecs_r1_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s1_1 + dksidheta_s * a3rVecs_r1_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 2][0]; da3[1] = da3_drds[1, 2][1]; da3[2] = da3_drds[1, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s].....
            Bab_rsAlternative.Bab_rs12_0 = ddksi_r * a3sVecs_s2_1 + ddksi_s * a3rVecs_r1_2 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs12_1 = ddheta_r * a3sVecs_s2_1 + ddheta_s * a3rVecs_r1_2 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs12_2 = dksidheta_r * a3sVecs_s2_1 + dksidheta_s * a3rVecs_r1_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s2_1 + dksidheta_s * a3rVecs_r1_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 0][0]; da3[1] = da3_drds[2, 0][1]; da3[2] = da3_drds[2, 0][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]..
            Bab_rsAlternative.Bab_rs20_0 = ddksi_r * a3sVecs_s0_2 + ddksi_s * a3rVecs_r2_0 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs20_1 = ddheta_r * a3sVecs_s0_2 + ddheta_s * a3rVecs_r2_0 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs20_2 = dksidheta_r * a3sVecs_s0_2 + dksidheta_s * a3rVecs_r2_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s0_2 + dksidheta_s * a3rVecs_r2_0 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 1][0]; da3[1] = da3_drds[2, 1][1]; da3[2] = da3_drds[2, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s].
            Bab_rsAlternative.Bab_rs21_0 = ddksi_r * a3sVecs_s1_2 + ddksi_s * a3rVecs_r2_1 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs21_1 = ddheta_r * a3sVecs_s1_2 + ddheta_s * a3rVecs_r2_1 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs21_2 = dksidheta_r * a3sVecs_s1_2 + dksidheta_s * a3rVecs_r2_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s1_2 + dksidheta_s * a3rVecs_r2_1 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 2][0]; da3[1] = da3_drds[2, 2][1]; da3[2] = da3_drds[2, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]..
            Bab_rsAlternative.Bab_rs22_0 = ddksi_r * a3sVecs_s2_2 + ddksi_s * a3rVecs_r2_2 + s11[0] * da3[0] + s11[1] * da3[1] + s11[2] * da3[2];
            Bab_rsAlternative.Bab_rs22_1 = ddheta_r * a3sVecs_s2_2 + ddheta_s * a3rVecs_r2_2 + s22[0] * da3[0] + s22[1] * da3[1] + s22[2] * da3[2];
            Bab_rsAlternative.Bab_rs22_2 = dksidheta_r * a3sVecs_s2_2 + dksidheta_s * a3rVecs_r2_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2]
                                         + dksidheta_r * a3sVecs_s2_2 + dksidheta_s * a3rVecs_r2_2 + s12[0] * da3[0] + s12[1] * da3[1] + s12[2] * da3[2];

            return Bab_rsAlternative;
        }

        internal static void CalculateCrossProduct(double[] vector1, double[] vector2, double[] result)
        {
            result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
            result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
            result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
        }

        internal void CalculateBendingDeformationMatrix(int controlPointsCount, double[] surfaceBasisVector3,
            IShapeFunction2D shapeFunctions, int j, double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVector1, double J1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double[,] BbendingOut)
        {
            var s10 = surfaceBasisVector1[0];
            var s11 = surfaceBasisVector1[1];
            var s12 = surfaceBasisVector1[2];

            var s20 = surfaceBasisVector2[0];
            var s21 = surfaceBasisVector2[1];
            var s22 = surfaceBasisVector2[2];

            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var s11_0 = surfaceBasisVectorDerivative1[0];
            var s11_1 = surfaceBasisVectorDerivative1[1];
            var s11_2 = surfaceBasisVectorDerivative1[2];

            var s22_0 = surfaceBasisVectorDerivative2[0];
            var s22_1 = surfaceBasisVectorDerivative2[1];
            var s22_2 = surfaceBasisVectorDerivative2[2];

            var s12_0 = surfaceBasisVectorDerivative12[0];
            var s12_1 = surfaceBasisVectorDerivative12[1];
            var s12_2 = surfaceBasisVectorDerivative12[2];

            for (int column = 0; column < controlPointsCount * 3; column += 3)
            {
                var dksi = shapeFunctions.DerivativeValuesKsi[column / 3, j];
                var dheta = shapeFunctions.DerivativeValuesHeta[column / 3, j];
                var d2Ksi = shapeFunctions.SecondDerivativeValuesKsi[column / 3, j];
                var d2Heta = shapeFunctions.SecondDerivativeValuesHeta[column / 3, j];
                var d2KsiHeta = shapeFunctions.SecondDerivativeValuesKsiHeta[column / 3, j];

                BbendingOut[0, column] = -d2Ksi * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s11 * s11_2 - s12 * s11_1) + dksi * (s21 * s11_2 - s22 * s11_1)) / J1;
                BbendingOut[0, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_2 - s12 * s11_0) + dksi * (s20 * s11_2 - s22 * s11_0)) / J1 - d2Ksi * s31;
                BbendingOut[0, column + 2] = -d2Ksi * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_1 - s11 * s11_0) + dksi * (s20 * s11_1 - s21 * s11_0)) / J1;

                BbendingOut[1, column] = -d2Heta * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s11 * s22_2 - s12 * s22_1) + dksi * (s21 * s22_2 - s22 * s22_1)) / J1;
                BbendingOut[1, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_2 - s12 * s22_0) + dksi * (s20 * s22_2 - s22 * s22_0)) / J1 - d2Heta * s31;
                BbendingOut[1, column + 2] = -d2Heta * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_1 - s11 * s22_0) + dksi * (s20 * s22_1 - s21 * s22_0)) / J1;

                BbendingOut[2, column] = -2 * d2KsiHeta * s30 - (2 * ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s11 * s12_2 - s12 * s12_1) + dksi * (s21 * s12_2 - s22 * s12_1))) / J1;
                BbendingOut[2, column + 1] = (2 * ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_2 - s12 * s12_0) + dksi * (s20 * s12_2 - s22 * s12_0))) / J1 - 2 * d2KsiHeta * s31;
                BbendingOut[2, column + 2] = -2 * d2KsiHeta * s32 - (2 * ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_1 - s11 * s12_0) + dksi * (s20 * s12_1 - s21 * s12_0))) / J1;

            }
        }

        internal double[,] CalculateHessian(ControlPoint[] controlPoints, IShapeFunction2D shapeFunctions, int j)
        {
            var hessianMatrix = new double[3, 3];
            for (var k = 0; k < controlPoints.Length; k++)
            {
                hessianMatrix[0, 0] += shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].X;
                hessianMatrix[0, 1] += shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Y;
                hessianMatrix[0, 2] += shapeFunctions.SecondDerivativeValuesKsi[k, j] * controlPoints[k].Z;
                hessianMatrix[1, 0] += shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].X;
                hessianMatrix[1, 1] += shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[1, 2] += shapeFunctions.SecondDerivativeValuesHeta[k, j] * controlPoints[k].Z;
                hessianMatrix[2, 0] += shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].X;
                hessianMatrix[2, 1] += shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[2, 2] += shapeFunctions.SecondDerivativeValuesKsiHeta[k, j] * controlPoints[k].Z;
            }

            return hessianMatrix;
        }

        internal void CalculateInitialConfigurationData(ControlPoint[] controlPoints,
            IShapeFunction2D shapeFunctions, IList<GaussLegendrePoint3D> gaussPoints)
        {
            var numberOfGP = gaussPoints.Count;
            InitialJ1 = new double[numberOfGP];
            initialSurfaceBasisVectors1 = new double[numberOfGP][];
            initialSurfaceBasisVectors2 = new double[numberOfGP][];
            initialUnitSurfaceBasisVectors3 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative1 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative2 = new double[numberOfGP][];
            initialSurfaceBasisVectorDerivative12 = new double[numberOfGP][];

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(controlPoints, shapeFunctions, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(controlPoints, shapeFunctions, j);
                initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector(jacobianMatrix, 0);
                initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector(jacobianMatrix, 1);
                var s3 = new double[3];
                CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j], s3);
                var norm = s3.Sum(t => t * t);
                InitialJ1[j] = Math.Sqrt(norm);
                var vector3 = new double[3];
                CalculateCrossProduct(initialSurfaceBasisVectors1[j], initialSurfaceBasisVectors2[j], vector3);
                initialUnitSurfaceBasisVectors3[j] = new double[]
                {
                    vector3[0] / InitialJ1[j],
                    vector3[1] / InitialJ1[j],
                    vector3[2] / InitialJ1[j],
                };

                initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector(hessianMatrix, 0);
                initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector(hessianMatrix, 1);
                initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector(hessianMatrix, 2);

                foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
                {
                    integrationPointMaterial.TangentVectorV1 = initialSurfaceBasisVectors1[j];
                    integrationPointMaterial.TangentVectorV2 = initialSurfaceBasisVectors2[j];
                    integrationPointMaterial.NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
                }
            }
        }

        internal void CalculateJacobian(ControlPoint[] controlPoints, IShapeFunction2D shapeFunctions, int j, double[,] jacobianOut)
        {
            jacobianOut[0, 0] = jacobianOut[0, 1] = jacobianOut[0, 2] =
                jacobianOut[1, 0] = jacobianOut[1, 1] = jacobianOut[1, 2] = 0.0;
            for (var k = 0; k < controlPoints.Length; k++)
            {
                jacobianOut[0, 0] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].X;
                jacobianOut[0, 1] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].Y;
                jacobianOut[0, 2] += shapeFunctions.DerivativeValuesKsi[k, j] * controlPoints[k].Z;
                jacobianOut[1, 0] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].X;
                jacobianOut[1, 1] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].Y;
                jacobianOut[1, 2] += shapeFunctions.DerivativeValuesHeta[k, j] * controlPoints[k].Z;
            }
        }

        internal void CalculateKbendingNL(ControlPoint[] controlPoints,
            ref Forces bendingMoments, IShapeFunction2D shapeFunctions, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double J1, int j, double[,] KbendingNLOut)
        {
            var a3rArray = new a3r[controlPoints.Length];
            for (var i = 0; i < controlPoints.Length; i++)
            {
                var a3r = new a3r();
                var dksi_r = shapeFunctions.DerivativeValuesKsi[i, j];
                var dheta_r = shapeFunctions.DerivativeValuesHeta[i, j];
                CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                    dheta_r, J1, ref a3r);
                a3rArray[i] = a3r;
            }

            for (int i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = shapeFunctions.DerivativeValuesKsi[i, j];
                var dheta_r = shapeFunctions.DerivativeValuesHeta[i, j];
                var d2Ksi_dr2 = shapeFunctions.SecondDerivativeValuesKsi[i, j];
                var d2Heta_dr2 = shapeFunctions.SecondDerivativeValuesHeta[i, j];
                var d2KsiHeta_dr2 = shapeFunctions.SecondDerivativeValuesKsiHeta[i, j];

                var a3r = a3rArray[i];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var d2Ksi_ds2 = shapeFunctions.SecondDerivativeValuesKsi[k, j];
                    var d2Heta_ds2 = shapeFunctions.SecondDerivativeValuesHeta[k, j];
                    var d2KsiHeta_ds2 = shapeFunctions.SecondDerivativeValuesKsiHeta[k, j];

                    var dksi_s = shapeFunctions.DerivativeValuesKsi[k, j];
                    var dheta_s = shapeFunctions.DerivativeValuesHeta[k, j];

                    var a3s = a3rArray[k];
                    a3rs = new a3rs();//Clear struct values
                    (a3rs a3rsAlternative, var da3tilde_drds, var da3tilde_dr, var da3tilde_ds,
                        double[] dnorma3_dr, double[] dnorma3_ds, double[,] dnorma3_drds, var a3_tilde, var da3_drds) =
                        Calculate_a3rs(Vector.CreateFromArray(surfaceBasisVector1), Vector.CreateFromArray(surfaceBasisVector2),
                            Vector.CreateFromArray(surfaceBasisVector3), J1, dksi_r, dksi_s, dheta_r, dheta_s);

                    Bab_rs Bab_rsAlternative = CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                        surfaceBasisVectorDerivative12, d2Ksi_dr2, d2Ksi_ds2, d2Heta_dr2, d2Heta_ds2, d2KsiHeta_dr2, d2KsiHeta_ds2,
                        a3rsAlternative, a3r, a3s, da3_drds);
                    Bab_rs = Bab_rsAlternative;

                    KbendingNLOut[i * 3 + 0, k * 3 + 0] -= (Bab_rs.Bab_rs00_0 * bendingMoments.v0 + Bab_rs.Bab_rs00_1 * bendingMoments.v1 + Bab_rs.Bab_rs00_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 0, k * 3 + 1] -= (Bab_rs.Bab_rs01_0 * bendingMoments.v0 + Bab_rs.Bab_rs01_1 * bendingMoments.v1 + Bab_rs.Bab_rs01_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 0, k * 3 + 2] -= (Bab_rs.Bab_rs02_0 * bendingMoments.v0 + Bab_rs.Bab_rs02_1 * bendingMoments.v1 + Bab_rs.Bab_rs02_2 * bendingMoments.v2);

                    KbendingNLOut[i * 3 + 1, k * 3 + 0] -= (Bab_rs.Bab_rs10_0 * bendingMoments.v0 + Bab_rs.Bab_rs10_1 * bendingMoments.v1 + Bab_rs.Bab_rs10_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 1, k * 3 + 1] -= (Bab_rs.Bab_rs11_0 * bendingMoments.v0 + Bab_rs.Bab_rs11_1 * bendingMoments.v1 + Bab_rs.Bab_rs11_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 1, k * 3 + 2] -= (Bab_rs.Bab_rs12_0 * bendingMoments.v0 + Bab_rs.Bab_rs12_1 * bendingMoments.v1 + Bab_rs.Bab_rs12_2 * bendingMoments.v2);

                    KbendingNLOut[i * 3 + 2, k * 3 + 0] -= (Bab_rs.Bab_rs20_0 * bendingMoments.v0 + Bab_rs.Bab_rs20_1 * bendingMoments.v1 + Bab_rs.Bab_rs20_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 2, k * 3 + 1] -= (Bab_rs.Bab_rs21_0 * bendingMoments.v0 + Bab_rs.Bab_rs21_1 * bendingMoments.v1 + Bab_rs.Bab_rs21_2 * bendingMoments.v2);
                    KbendingNLOut[i * 3 + 2, k * 3 + 2] -= (Bab_rs.Bab_rs22_0 * bendingMoments.v0 + Bab_rs.Bab_rs22_1 * bendingMoments.v1 + Bab_rs.Bab_rs22_2 * bendingMoments.v2);

                }
            }
        }

        internal void CalculateKmembraneNL(ControlPoint[] controlPoints, ref Forces membraneForces, IShapeFunction2D shapeFunctions,
            int j, double[,] KmembraneNLOut)
        {
            for (var i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = shapeFunctions.DerivativeValuesKsi[i, j];
                var dheta_r = shapeFunctions.DerivativeValuesHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var dksi_s = shapeFunctions.DerivativeValuesKsi[k, j];
                    var dheta_s = shapeFunctions.DerivativeValuesHeta[k, j];


                    var aux = membraneForces.v0 * dksi_r * dksi_s +
                              membraneForces.v1 * dheta_r * dheta_s +
                              membraneForces.v2 * (dksi_r * dheta_s + dksi_s * dheta_r);

                    KmembraneNLOut[i * 3, k * 3] += aux;
                    KmembraneNLOut[i * 3 + 1, k * 3 + 1] += aux;
                    KmembraneNLOut[i * 3 + 2, k * 3 + 2] += aux;
                }
            }
        }

        internal void CalculateMembraneDeformationMatrix(int controlPointsCount, IShapeFunction2D shapeFunctions, int j,
            double[] surfaceBasisVector1, double[] surfaceBasisVector2, double[,] BmembraneOut)
        {
            var s1_0 = surfaceBasisVector1[0];
            var s1_1 = surfaceBasisVector1[1];
            var s1_2 = surfaceBasisVector1[2];

            var s2_0 = surfaceBasisVector2[0];
            var s2_1 = surfaceBasisVector2[1];
            var s2_2 = surfaceBasisVector2[2];

            for (int column = 0; column < controlPointsCount * 3; column += 3)
            {
                var dKsi = shapeFunctions.DerivativeValuesKsi[column / 3, j];
                var dHeta = shapeFunctions.DerivativeValuesHeta[column / 3, j];

                BmembraneOut[0, column] = dKsi * s1_0;
                BmembraneOut[0, column + 1] = dKsi * s1_1;
                BmembraneOut[0, column + 2] = dKsi * s1_2;

                BmembraneOut[1, column] = dHeta * s2_0;
                BmembraneOut[1, column + 1] = dHeta * s2_1;
                BmembraneOut[1, column + 2] = dHeta * s2_2;

                BmembraneOut[2, column] = dHeta * s1_0 + dKsi * s2_0;
                BmembraneOut[2, column + 1] = dHeta * s1_1 + dKsi * s2_1;
                BmembraneOut[2, column + 2] = dHeta * s1_2 + dKsi * s2_2;
            }
        }

        internal double[] CalculateSurfaceBasisVector(double[,] Matrix, int row)
        {
            var surfaceBasisVector1 = new double[3];
            surfaceBasisVector1[0] = Matrix[row, 0];
            surfaceBasisVector1[1] = Matrix[row, 1];
            surfaceBasisVector1[2] = Matrix[row, 2];
            return surfaceBasisVector1;
        }

        internal ControlPoint[] CurrentControlPoint(ControlPoint[] controlPoints)
        {
            var cp = new ControlPoint[controlPoints.Length];

            for (int i = 0; i < controlPoints.Length; i++)
            {
                cp[i] = new ControlPoint()
                {
                    ID = controlPoints[i].ID,
                    X = controlPoints[i].X + _solution[i * 3],
                    Y = controlPoints[i].Y + _solution[i * 3 + 1],
                    Z = controlPoints[i].Z + _solution[i * 3 + 2],
                    Ksi = controlPoints[i].Ksi,
                    Heta = controlPoints[i].Heta,
                    Zeta = controlPoints[i].Zeta,
                    WeightFactor = controlPoints[i].WeightFactor
                };
            }

            return cp;
        }

        internal (double[,] MembraneConstitutiveMatrix, double[,] BendingConstitutiveMatrix, double[,]
            CouplingConstitutiveMatrix) IntegratedConstitutiveOverThickness(GaussLegendrePoint3D midSurfaceGaussPoint)
        {
            var MembraneConstitutiveMatrix = new double[3, 3];
            var BendingConstitutiveMatrix = new double[3, 3];
            var CouplingConstitutiveMatrix = new double[3, 3];

            foreach (var keyValuePair in materialsAtThicknessGP[midSurfaceGaussPoint])
            {
                var thicknessPoint = keyValuePair.Key;
                var material = keyValuePair.Value;
                var constitutiveMatrixM = material.ConstitutiveMatrix;
                double tempc = 0;
                double w = thicknessPoint.WeightFactor;
                double z = thicknessPoint.Zeta;
                for (int i = 0; i < 3; i++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        tempc = constitutiveMatrixM[i, k];
                        MembraneConstitutiveMatrix[i, k] += tempc * w;
                        CouplingConstitutiveMatrix[i, k] += tempc * w * z;
                        BendingConstitutiveMatrix[i, k] += tempc * w * z * z;
                    }
                }
            }

            return (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix);
        }

        internal void IntegratedStressesOverThickness(
            GaussLegendrePoint3D midSurfaceGaussPoint, ref Forces MembraneForces, ref Forces BendingMoments)
        {
            MembraneForces = new Forces();
            BendingMoments = new Forces();
            var thicknessPoints = thicknessIntegrationPoints[midSurfaceGaussPoint];

            for (int i = 0; i < thicknessPoints.Count; i++)
            {
                var thicknessPoint = thicknessPoints[i];
                var material = materialsAtThicknessGP[midSurfaceGaussPoint][thicknessPoints[i]];
                var w = thicknessPoint.WeightFactor;
                var z = thicknessPoint.Zeta;
                MembraneForces.v0 += material.Stresses[0] * w;
                MembraneForces.v1 += material.Stresses[1] * w;
                MembraneForces.v2 += material.Stresses[2] * w;

                BendingMoments.v0 -= material.Stresses[0] * w * z;
                BendingMoments.v1 -= material.Stresses[1] * w * z;
                BendingMoments.v2 -= material.Stresses[2] * w * z;
            }
        }

        private static void Calculate_da3tilde_dr(Vector surfaceBasisVector1, Vector surfaceBasisVector2, double dksi_r,
             double dHeta_r, double[][] da3tilde_dr)
        {
            //da3tilde_dr[r1] = a1r.GetColumn(r1).CrossProduct(surfaceBasisVector2) + surfaceBasisVector1.CrossProduct(a2r.GetColumn(r1));

            da3tilde_dr[0] = new double[]
            {
                 0,
                 -dksi_r*surfaceBasisVector2[2]+surfaceBasisVector1[2]*dHeta_r,
                 dksi_r*surfaceBasisVector2[1]-surfaceBasisVector1[1]*dHeta_r
            };

            da3tilde_dr[1] = new double[]
            {
                 dksi_r*surfaceBasisVector2[2]-surfaceBasisVector1[2]*dHeta_r,
                 0,
                 -dksi_r*surfaceBasisVector2[0]+surfaceBasisVector1[0]*dHeta_r
            };

            da3tilde_dr[2] = new double[]
            {
                 -dksi_r*surfaceBasisVector2[1]+dHeta_r*surfaceBasisVector1[1],
                 dksi_r*surfaceBasisVector2[0]-dHeta_r*surfaceBasisVector1[0],
                 0
            };
        }

        private static void Calculate_da3tilde_drds(double dKsi_r, double dKsi_s, double dHeta_r, double dHeta_s,
             double[,][] da3tilde_drds)
        {
            //da3tilde_drds[r1, s1] = a1r.GetColumn(r1).CrossProduct(a2s.GetColumn(s1)) +
            //                        a1s.GetColumn(s1).CrossProduct(a2r.GetColumn(r1));

            var dksiRxdHetaS = dKsi_r * dHeta_s;
            var dHetaRxdKsiS = dHeta_r * dKsi_s;
            da3tilde_drds[0, 0] = new double[3];
            da3tilde_drds[0, 1] = new double[] { 0, 0, dksiRxdHetaS - dHetaRxdKsiS };
            da3tilde_drds[0, 2] = new double[] { 0, dHetaRxdKsiS - dksiRxdHetaS, 0 };

            da3tilde_drds[1, 0] = new double[] { 0, 0, dHetaRxdKsiS - dksiRxdHetaS };
            da3tilde_drds[1, 1] = new double[3];
            da3tilde_drds[1, 2] = new double[] { dksiRxdHetaS - dHetaRxdKsiS, 0, 0 };

            da3tilde_drds[2, 0] = new double[] { 0, dksiRxdHetaS - dHetaRxdKsiS, 0 };
            da3tilde_drds[2, 1] = new double[] { dHetaRxdKsiS - dksiRxdHetaS, 0, 0 };
            da3tilde_drds[2, 2] = new double[3];
        }

        private static double[] CalculateTerm525(Vector surfaceBasisVector3, double J1, double[] dnorma3_dr,
             double[][] da3tilde_dr, double[] dnorma3_ds, double[][] da3tilde_ds)
        {
            double[] a3_tilde;
            a3_tilde = new double[]
            {
                 surfaceBasisVector3[0] * J1,
                 surfaceBasisVector3[1] * J1,
                 surfaceBasisVector3[2] * J1,
            };
            for (int r1 = 0; r1 < 3; r1++)
            {
                //dnorma3_dr[r1] = (a3_tilde.DotProduct(da3tilde_dr[r1])) / J1;
                dnorma3_dr[r1] = (a3_tilde[0] * da3tilde_dr[r1][0] + a3_tilde[1] * da3tilde_dr[r1][1] +
                                  a3_tilde[2] * da3tilde_dr[r1][2]) / J1;
            }

            for (int s1 = 0; s1 < 3; s1++)
            {
                //dnorma3_ds[s1] = (a3_tilde.DotProduct(da3tilde_ds[s1])) / J1;
                dnorma3_ds[s1] = (a3_tilde[0] * da3tilde_ds[s1][0] + a3_tilde[1] * da3tilde_ds[s1][1] +
                                  a3_tilde[2] * da3tilde_ds[s1][2]) / J1;
            }

            return a3_tilde;
        }

        private static void CalculateTerm531(double J1, double[,][] da3tilde_drds, double[] a3_tilde, double[][] da3tilde_dr,
             double[][] da3tilde_ds, double[,] dnorma3_drds)
        {
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    //double firstNumerator = da3tilde_drds[r1, s1].DotProduct(a3_tilde) + da3tilde_dr[r1].DotProduct(da3tilde_ds[s1]);
                    double firstNumerator = da3tilde_drds[r1, s1][0] * a3_tilde[0] + da3tilde_drds[r1, s1][1] * a3_tilde[1] +
                                            da3tilde_drds[r1, s1][2] * a3_tilde[2] +
                                            da3tilde_dr[r1][0] * da3tilde_ds[s1][0] + da3tilde_dr[r1][1] * da3tilde_ds[s1][1] +
                                            da3tilde_dr[r1][2] * da3tilde_ds[s1][2];
                    double firstDenominator = J1;
                    //double secondNumerator = (da3tilde_dr[r1].DotProduct(a3_tilde)) * (da3tilde_ds[s1].DotProduct(a3_tilde));
                    double secondNumerator = (da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] +
                                              da3tilde_dr[r1][2] * a3_tilde[2]) *
                                             (da3tilde_ds[s1][0] * a3_tilde[0] + da3tilde_ds[s1][1] * a3_tilde[1] +
                                              da3tilde_ds[s1][2] * a3_tilde[2]);
                    double secondDenominator = Math.Pow(J1, 3);

                    dnorma3_drds[r1, s1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);
                }
            }
        }

        private static void CalculateTerm532(double J1, double[,][] da3tilde_drds, double[] dnorma3_ds, double[][] da3tilde_dr,
            double[] dnorma3_dr, double[][] da3tilde_ds, double[,] dnorma3_drds, double[] a3_tilde, double[,][] da3_drds)
        {
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    var firstVec_0 = da3tilde_drds[r1, s1][0] / J1;
                    var firstVec_1 = da3tilde_drds[r1, s1][1] / J1;
                    var firstVec_2 = da3tilde_drds[r1, s1][2] / J1;

                    double scale2 = -((double)1 / (Math.Pow(J1, 2))); //denominator of vectors 2 3 and 4 and a minus.

                    var scale3 = dnorma3_ds[s1] * scale2;
                    var secondVec_0 = da3tilde_dr[r1][0] * scale3;
                    var secondVec_1 = da3tilde_dr[r1][1] * scale3;
                    var secondVec_2 = da3tilde_dr[r1][2] * scale3;

                    var scale4 = dnorma3_dr[r1] * scale2;
                    var thirdVec_0 = da3tilde_ds[s1][0] * scale4;
                    var thirdVec_1 = da3tilde_ds[s1][1] * scale4;
                    var thirdVec_2 = da3tilde_ds[s1][2] * scale4;

                    var scale6 = dnorma3_drds[r1, s1] * scale2;
                    var fourthVec_0 = a3_tilde[0] * scale6;
                    var fourthVec_1 = a3_tilde[1] * scale6;
                    var fourthVec_2 = a3_tilde[2] * scale6;

                    double scale5 = ((double)1) / Math.Pow(J1, 3);

                    var scale7 = 2 * dnorma3_dr[r1] * dnorma3_ds[s1] * scale5;
                    var fifthvector_0 = a3_tilde[0] * scale7;
                    var fifthvector_1 = a3_tilde[1] * scale7;
                    var fifthvector_2 = a3_tilde[2] * scale7;

                    da3_drds[r1, s1] = new double[]
                    {
                        firstVec_0 + secondVec_0 + thirdVec_0 + fourthVec_0 + fifthvector_0,
                        firstVec_1 + secondVec_1 + thirdVec_1 + fourthVec_1 + fifthvector_1,
                        firstVec_2 + secondVec_2 + thirdVec_2 + fourthVec_2 + fifthvector_2,
                    };
                }
            }
        }

        private void CalculateLinearStiffness(ControlPoint[] elementControlPoints, IShapeFunction2D shapeFunctions, int j,
            double[] surfaceBasisVector1, double[] surfaceBasisVector2, double[,] Bmembrane, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double J1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double[,] Bbending, GaussLegendrePoint3D[] gaussPoints,
            double[,] BmTranspose, int bRows, int bCols, double[,] BbTranspose, double wFactor,
            double[,] BmTransposeMultStiffness, double[,] BbTransposeMultStiffness, double[,] BmbTransposeMultStiffness,
            double[,] BbmTransposeMultStiffness, double[,] stiffnessMatrix, double[,] KmembraneL, double[,] KbendingL)
        {
            CalculateMembraneDeformationMatrix(elementControlPoints.Length, shapeFunctions, j, surfaceBasisVector1,
                surfaceBasisVector2, Bmembrane);
            CalculateBendingDeformationMatrix(elementControlPoints.Length, surfaceBasisVector3, shapeFunctions, j,
                surfaceBasisVector2, surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, Bbending);

            var (MembraneConstitutiveMatrix, BendingConstitutiveMatrix, CouplingConstitutiveMatrix) =
                IntegratedConstitutiveOverThickness(gaussPoints[j]);


            double tempb = 0;
            double tempm = 0;
            Array.Clear(BmTranspose, 0, bRows * bCols);
            Array.Clear(BbTranspose, 0, bRows * bCols);
            for (int i = 0; i < bRows; i++)
            {
                for (int k = 0; k < bCols; k++)
                {
                    BmTranspose[k, i] = Bmembrane[i, k] * wFactor;
                    BbTranspose[k, i] = Bbending[i, k] * wFactor;
                }
            }

            double tempcm = 0;
            double tempcb = 0;
            double tempcc = 0;
            Array.Clear(BmTransposeMultStiffness, 0, bRows * bCols);
            Array.Clear(BbTransposeMultStiffness, 0, bRows * bCols);
            Array.Clear(BmbTransposeMultStiffness, 0, bRows * bCols);
            Array.Clear(BbmTransposeMultStiffness, 0, bRows * bCols);
            for (int i = 0; i < bCols; i++)
            {
                for (int k = 0; k < bRows; k++)
                {
                    tempm = BmTranspose[i, k];
                    tempb = BbTranspose[i, k];
                    for (int m = 0; m < bRows; m++)
                    {
                        tempcm = MembraneConstitutiveMatrix[k, m];
                        tempcb = BendingConstitutiveMatrix[k, m];
                        tempcc = CouplingConstitutiveMatrix[k, m];

                        BmTransposeMultStiffness[i, m] += tempm * tempcm;
                        BbTransposeMultStiffness[i, m] += tempb * tempcb;
                        BmbTransposeMultStiffness[i, m] += tempm * tempcc;
                        BbmTransposeMultStiffness[i, m] += tempb * tempcc;
                    }
                }
            }

            double tempmb = 0;
            double tempbm = 0;
            double mem = 0;
            double ben = 0;
            for (int i = 0; i < bCols; i++)
            {
                for (int k = 0; k < bRows; k++)
                {
                    tempm = BmTransposeMultStiffness[i, k];
                    tempb = BbTransposeMultStiffness[i, k];
                    tempmb = BmbTransposeMultStiffness[i, k];
                    tempbm = BbmTransposeMultStiffness[i, k];

                    for (int m = 0; m < bCols; m++)
                    {
                        mem = Bmembrane[k, m];
                        ben = Bbending[k, m];
                        stiffnessMatrix[i, m] += tempm * mem + tempb * ben + tempmb * ben + tempbm * mem;
                    }
                }
            }
        }

        private void CalculateNonLinearStiffness(GaussLegendrePoint3D[] gaussPoints, int j, double[,] KmembraneNL, int bCols,
            double[,] KbendingNL, ControlPoint[] elementControlPoints, IShapeFunction2D shapeFunctions, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, double J1,
            double[,] stiffnessMatrix, double wFactor, ref Forces MembraneForces, ref Forces BendingMoments)
        {
            IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

            Array.Clear(KmembraneNL, 0, bCols * bCols);
            Array.Clear(KbendingNL, 0, bCols * bCols);

            CalculateKmembraneNL(elementControlPoints, ref MembraneForces, shapeFunctions, j, KmembraneNL);
            CalculateKbendingNL(elementControlPoints, ref BendingMoments, shapeFunctions,
                surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3,
                surfaceBasisVectorDerivative1,
                surfaceBasisVectorDerivative2,
                surfaceBasisVectorDerivative12, J1, j, KbendingNL);

            for (var i = 0; i < stiffnessMatrix.GetLength(0); i++)
            {
                for (var k = 0; k < stiffnessMatrix.GetLength(1); k++)
                {
                    stiffnessMatrix[i, k] += (KmembraneNL[i, k] + KbendingNL[i, k]) * wFactor;
                }
            }
        }
        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(KirchhoffLoveShellNL shellElement)
        {
            var gauss = new GaussQuadrature();
            var medianSurfaceGP = gauss.CalculateElementGaussPoints(_degreeKsi,
                _degreeHeta, shellElement.Knots.ToList());
            foreach (var point in medianSurfaceGP)
            {
                var gp = gauss.CalculateElementGaussPoints(ThicknessIntegrationDegree,
                    new List<Knot>
                    {
                        new Knot() {ID = 0, Ksi = -shellElement.Thickness / 2, Heta = point.Heta},
                        new Knot() {ID = 1, Ksi = shellElement.Thickness / 2, Heta = point.Heta},
                    }).ToList();

                thicknessIntegrationPoints.Add(point,
                    gp.Select(g => new GaussLegendrePoint3D(point.Ksi, point.Heta, g.Ksi, g.WeightFactor))
                        .ToList());
            }

            return medianSurfaceGP;
        }
    }
}