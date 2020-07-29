using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using Element = ISAAR.MSolve.IGA.Entities.Element;

[assembly: InternalsVisibleTo("ISAAR.MSolve.IGA.Tests")]

namespace ISAAR.MSolve.IGA.Elements
{
    public class NurbsKirchhoffLoveShellElementNLDevelop : Element, IStructuralIsogeometricElement, ISurfaceLoadedElement
    {
        protected static readonly IDofType[] ControlPointDofTypes =
            {StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ};

        private IDofType[][] dofTypes;

        private Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>> thicknessIntegrationPoints =
            new Dictionary<GaussLegendrePoint3D, List<GaussLegendrePoint3D>>();

        internal Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>
            materialsAtThicknessGP =
                new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IShellMaterial>>();

        private bool isInitialized;
        internal double[] _solution;

        public NurbsKirchhoffLoveShellElementNLDevelop(IShellMaterial shellMaterial, IList<Knot> elementKnots,
            IList<ControlPoint> elementControlPoints, Patch patch, double thickness)
        {
            Contract.Requires(shellMaterial != null);
            this.Patch = patch;
            this.Thickness = thickness;
            foreach (var knot in elementKnots)
            {
                if (!KnotsDictionary.ContainsKey(knot.ID))
                    this.KnotsDictionary.Add(knot.ID, knot);
            }

            _solution = new double[3 * elementControlPoints.Count];

            foreach (var controlPoint in elementControlPoints)
            {
                if (!ControlPointsDictionary.ContainsKey(controlPoint.ID))
                    ControlPointsDictionary.Add(controlPoint.ID, controlPoint);
            }

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

        public CellType CellType { get; } = CellType.Unknown;

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public bool MaterialModified => false;

        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads) =>
            throw new NotImplementedException();

        public double[,] CalculateDisplacementsForPostProcessing(Element element, Matrix localDisplacements)
        {
            var nurbsElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var knotParametricCoordinatesKsi = Vector.CreateFromArray(Knots.Select(k => k.Ksi).ToArray());
            var knotParametricCoordinatesHeta = Vector.CreateFromArray(Knots.Select(k => k.Heta).ToArray());

            var nurbs = new Nurbs2D(nurbsElement, nurbsElement.ControlPoints.ToArray(), knotParametricCoordinatesKsi,
                knotParametricCoordinatesHeta);

            var knotDisplacements = new double[4, 3];
            var paraviewKnotRenumbering = new int[] {0, 3, 1, 2};
            for (var j = 0; j < knotDisplacements.GetLength(0); j++)
            {
                for (int i = 0; i < element.ControlPoints.Count(); i++)
                {
                    knotDisplacements[paraviewKnotRenumbering[j], 0] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 0];
                    knotDisplacements[paraviewKnotRenumbering[j], 1] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 1];
                    knotDisplacements[paraviewKnotRenumbering[j], 2] +=
                        nurbs.NurbsValues[i, j] * localDisplacements[i, 2];
                }
            }

            return knotDisplacements;
        }

        public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var elementNodalForces = new double[shellElement.ControlPointsDictionary.Count * 3];
            var elementNodalMembraneForces = new double[shellElement.ControlPointsDictionary.Count * 3];
            var elementNodalBendingForces = new double[shellElement.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(_controlPoints);
            var nurbs = CalculateShapeFunctions(shellElement, _controlPoints);
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

            var Bmembrane = new double[3, _controlPoints.Length * 3];
            var Bbending = new double[3, _controlPoints.Length * 3];
            var numberOfControlPoints = _controlPoints.Length;
            var MembraneForces = new Forces();
            var BendingMoments = new Forces();

            var forcesDevelop = new double[shellElement.ControlPointsDictionary.Count * 3];

            bool runNewForces = true;

            for (int j = 0; j < gaussPoints.Length; j++)
            {

                CalculateJacobian(newControlPoints, nurbs, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, nurbs, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

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


                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                if (!runNewForces)
                {
                    CalculateMembraneDeformationMatrix(numberOfControlPoints, nurbs, j, surfaceBasisVector1,
                        surfaceBasisVector2, Bmembrane);
                    CalculateBendingDeformationMatrix(numberOfControlPoints, surfaceBasisVector3, nurbs, j,
                        surfaceBasisVector2,
                        surfaceBasisVectorDerivative1, surfaceBasisVector1, J1, surfaceBasisVectorDerivative2,
                        surfaceBasisVectorDerivative12, Bbending);


                    // Bbending = Matrix.CreateFromArray(Bbending).Scale(-1).CopytoArray2D();
                    IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

                    if (j == ElementStiffnesses.gpNumberToCheck)
                    {
                        var Bmem_matrix = Matrix.CreateFromArray(Bmembrane);
                        var Bben_matrix = Matrix.CreateFromArray(Bbending);
                        for (int i1 = 0; i1 < Bmembrane.GetLength(1); i1++)
                        {
                            ElementStiffnesses.ProccessVariable(4, Bmem_matrix.GetColumn(i1).CopyToArray(), true, i1);
                            ElementStiffnesses.ProccessVariable(5, Bben_matrix.GetColumn(i1).CopyToArray(), true, i1);

                        }
                    }


                    if (j == ElementStiffnesses.gpNumberToCheck)
                    {
                        for (int i = 0; i < _controlPoints.Length; i++)
                        {
                            var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                            var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);

                            for (int i1 = 0; i1 < 3; i1++)
                            {
                                ElementStiffnesses.ProccessVariable(6, a1r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                                ElementStiffnesses.ProccessVariable(7, a2r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                            }
                        }
                    }




                    for (int i = 0; i < Bmembrane.GetLength(1); i++)
                    {
                        elementNodalForces[i] +=
                            (Bmembrane[0, i] * MembraneForces.v0 * wfactor + Bbending[0, i] * BendingMoments.v0 * wfactor) +
                            (Bmembrane[1, i] * MembraneForces.v1 * wfactor + Bbending[1, i] * BendingMoments.v1 * wfactor) +
                            (Bmembrane[2, i] * MembraneForces.v2 * wfactor + Bbending[2, i] * BendingMoments.v2 * wfactor);

                        if ((ElementStiffnesses.saveForcesState1 | ElementStiffnesses.saveForcesState2s) | ElementStiffnesses.saveForcesState0)
                        {
                            elementNodalMembraneForces[i] +=
                             (Bmembrane[0, i] * MembraneForces.v0 * wfactor) +
                             (Bmembrane[1, i] * MembraneForces.v1 * wfactor) +
                             (Bmembrane[2, i] * MembraneForces.v2 * wfactor);

                            elementNodalBendingForces[i] +=
                             (Bbending[0, i] * BendingMoments.v0 * wfactor) +
                             (Bbending[1, i] * BendingMoments.v1 * wfactor) +
                             (Bbending[2, i] * BendingMoments.v2 * wfactor);
                        }
                    }

                    if (j == ElementStiffnesses.gpNumberToCheck)
                    {
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }

                        ElementStiffnesses.ProccessVariable(8, surfaceBasisVector3, false);
                        //ElementStiffnesses.ProccessVariable(9, surfaceBasisVector2, false);
                        ElementStiffnesses.ProccessVariable(10, new double[] { surfaceBasisVector3[0] * J1, surfaceBasisVector3[1] * J1, surfaceBasisVector3[2] * J1 }, false);

                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                    }

                }


                #region develop formulation
                var forcesDevelopGp = new double[shellElement.ControlPointsDictionary.Count * 3];

                var thicknessGPoints = thicknessIntegrationPoints[gaussPoints[j]];
                var materialpoint = materialsAtThicknessGP[gaussPoints[j]][thicknessGPoints[0]];
                var transformations = new ShellElasticMaterial2DtransformationbDefGrad() { YoungModulus = materialpoint.YoungModulus, PoissonRatio = materialpoint.PoissonRatio };


                var elementControlPoints = CurrentControlPoint(_controlPoints);

                var a11 = Vector.CreateFromArray(surfaceBasisVectorDerivative1);
                var a22 = Vector.CreateFromArray(surfaceBasisVectorDerivative2);
                var a12 = Vector.CreateFromArray(surfaceBasisVectorDerivative12);
                var a1 = Vector.CreateFromArray(surfaceBasisVector1);
                var a2 = Vector.CreateFromArray(surfaceBasisVector2);
                var a3 = Vector.CreateFromArray(surfaceBasisVector3); // einai to mono pou einai normalised
                Vector a3_tilde = a3.Scale(J1);

                (Vector da3tilde_dksi, Vector da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, Vector da3_dksi, Vector da3_dheta) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1, a2, a11, a22, a12, a3, J1);

                if (j == ElementStiffnesses.gpNumberToCheck)
                {
                    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                    ElementStiffnesses.ProccessVariable(11, new double[] { J1 }, false);
                    ElementStiffnesses.ProccessVariable(12, da3tilde_dksi.CopyToArray(), false);
                    ElementStiffnesses.ProccessVariable(13, da3tilde_dheta.CopyToArray(), false);
                    ElementStiffnesses.ProccessVariable(14, new double[] { da3norm_dksi }, false);
                    ElementStiffnesses.ProccessVariable(15, new double[] { da3norm_dheta }, false);
                    ElementStiffnesses.ProccessVariable(16, da3_dksi.CopyToArray(), false);
                    ElementStiffnesses.ProccessVariable(17, da3_dheta.CopyToArray(), false);
                    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                }

                #region original Config
                var originalControlPoints = shellElement.ControlPoints.ToArray();
                var originalHessianMatrix = CalculateHessian(originalControlPoints, nurbs, j);
                double[,] jacobian_init = new double[3, 3];
                CalculateJacobian(originalControlPoints, nurbs, j, jacobian_init);
                var a1_init = Vector.CreateFromArray(CalculateSurfaceBasisVector1(jacobian_init, 0));
                var a2_init = Vector.CreateFromArray(CalculateSurfaceBasisVector1(jacobian_init, 1));
                var a3_init = Vector.CreateFromArray(new double[3]
                {
                    a1_init[1] * a2_init[2] - a1_init[2] * a2_init[1],
                    a1_init[2] * a2_init[0] - a1_init[0] * a2_init[2],
                    a1_init[0] * a2_init[1] - a1_init[1] * a2_init[0]
                });

                var J1_init = Math.Sqrt(a3_init[0] * a3_init[0] +
                                        a3_init[1] * a3_init[1] +
                                        a3_init[2] * a3_init[2]);

                a3_init[0] /= J1_init;
                a3_init[1] /= J1_init;
                a3_init[2] /= J1_init;

                var a11_init = Vector.CreateFromArray(CalculateSurfaceBasisVector1(originalHessianMatrix, 0));
                var a22_init = Vector.CreateFromArray(CalculateSurfaceBasisVector1(originalHessianMatrix, 1));
                var a12_init = Vector.CreateFromArray(CalculateSurfaceBasisVector1(originalHessianMatrix, 2));
                #endregion

                (Vector da3tilde_dksi_init, Vector da3tilde_dheta_init, double da3norm_dksi_init, double da3norm_dheta_init, Vector da3_dksi_init, Vector da3_dheta_init) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1_init, a2_init, a11_init, a22_init, a12_init, a3_init, J1_init);


                for (int i = 0; i < elementControlPoints.Length; i++)
                {
                    var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                    var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                    var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                    var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                    var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);

                    //var a3r = new a3r();
                    var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                    var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];
                    //CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                    //    dheta_r, J1, ref a3r); //gia to a3r


                    //5.24
                    var da3tilde_dr = new double[3][];
                    Calculate_da3tilde_dr(a1, a2, dksi_r, dheta_r, da3tilde_dr);

                    //5.25
                    double[] dnorma3_dr = new double[3];
                    a3_tilde = Vector.CreateFromArray(CalculateTerm525(a3, J1, dnorma3_dr, da3tilde_dr));

                    //5.30 b
                    (Vector[] da3tilde_dksidr, Vector[] da3tilde_dhetadr) = Calculate_da3tilde_dksidr(a1r, a2r, a11r, a22r, a12r, a1, a2, a11, a22, a12);

                    //5.31 b
                    (double[] da3norm_dksidr, double[] da3norm_dhetadr) = Calculate_da3norm_dksidr(da3tilde_dksidr, da3tilde_dhetadr,
                        a3_tilde, da3tilde_dksi, da3tilde_dheta, da3tilde_dr, J1);

                    //5.32 b
                    (Vector[] da3_dksidr, Vector[] da3_dhetadr) = Calculate_da3_dksidr(da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksi, da3tilde_dheta,
                        dnorma3_dr, a3_tilde, da3norm_dksidr, da3norm_dhetadr, da3norm_dksi, da3norm_dheta, J1, da3tilde_dr);

                    if (j == ElementStiffnesses.gpNumberToCheck)
                    {
                        // conrol point einai to iota
                        //var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                        //var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);



                        for (int i1 = 0; i1 < 3; i1++)
                        {
                            ElementStiffnesses.ProccessVariable(11, new double[1] { dnorma3_dr[i1] }, true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(12, da3tilde_dksidr[i1].CopyToArray(), true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(13, da3tilde_dhetadr[i1].CopyToArray(), true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(14, new double[1] { da3norm_dksidr[i1] }, true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(15, new double[1] { da3norm_dhetadr[i1] }, true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(16, da3_dksidr[i1].CopyToArray(), true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(17, da3_dhetadr[i1].CopyToArray(), true, 3 * i + i1);
                        }

                    }




                    var thicknessPoints = thicknessIntegrationPoints[gaussPoints[j]];
                    double[,] forceIntegration = new double[thicknessGPoints.Count(), 3]; //[thickness,dof]
                    double[,] thicknesCoeffs = new double[thicknessGPoints.Count(), 3]; //[thickness,dof]
                    for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                    {
                        var thicknessPoint = thicknessPoints[i1];
                        var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                        var w = thicknessPoint.WeightFactor;
                        var z = thicknessPoint.Zeta;
                        //MembraneForces.v0 += material.Stresses[0] * w;
                        //MembraneForces.v1 += material.Stresses[1] * w;
                        //MembraneForces.v2 += material.Stresses[2] * w;

                        //BendingMoments.v0 -= material.Stresses[0] * w * z;
                        //BendingMoments.v1 -= material.Stresses[1] * w * z;
                        //BendingMoments.v2 -= material.Stresses[2] * w * z;

                        //for (int r1 = 0; r1 < 3; r1++)
                        //{
                        //    Vector dg1_dr = a1r.GetColumn(r1) - da3_dksidr[r1].Scale(z);
                        //    Vector dg2_dr = a2r.GetColumn(r1) - da3_dhetadr[r1].Scale(z);
                        //}

                        // (14) 
                        Vector G1 = a1_init + da3_dksi_init.Scale(z);
                        Vector G2 = a2_init + da3_dheta_init.Scale(z);

                        double G1_norm_sqred = G1.DotProduct(G1);
                        double G2_norm_sqred = G2.DotProduct(G2);
                        double G3_norm_sqred = a3.DotProduct(a3);

                        double[] G_1 = new double[3] { G1[0] / G1_norm_sqred, G1[1] / G1_norm_sqred, G1[2] / G1_norm_sqred };
                        double[] G_2 = new double[3] { G2[0] / G2_norm_sqred, G2[1] / G2_norm_sqred, G2[2] / G2_norm_sqred };
                        double[] G_3 = new double[3] { a3[0] / G3_norm_sqred, a3[1] / G3_norm_sqred, a3[2] / G3_norm_sqred };

                        Vector g1 = a1 + da3_dksi.Scale(z);
                        Vector g2 = a2 + da3_dheta.Scale(z);

                        double[,] F_3D = new double[3, 3] { { g1[0]*G_1[0]+g2[0]*G_2[0], g1[0]*G_1[1]+g2[0]*G_2[1], g1[0]*G_1[2]+g2[0]*G_2[2] },
                                                            { g1[1]*G_1[0]+g2[1]*G_2[0], g1[1]*G_1[1]+g2[1]*G_2[1], g1[1]*G_1[2]+g2[1]*G_2[2] },
                                                            { g1[2]*G_1[0]+g2[2]*G_2[0], g1[2]*G_1[1]+g2[2]*G_2[1], g1[2]*G_1[2]+g2[2]*G_2[2] },
                        };

                        if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0))
                        {
                            if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                            double[] F_3D_vec = { F_3D[0, 0], F_3D[1, 1], F_3D[2, 2], F_3D[0, 1], F_3D[1, 2], F_3D[2, 0], F_3D[0, 2], F_3D[1, 0], F_3D[2, 1] };//  .
                            ElementStiffnesses.ProccessVariable(18, F_3D_vec, false);
                            if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                        }

                        double[,] tgi = new double[3, 3] { { g1[0], g2[0], a3[0] }, { g1[1], g2[1], a3[1] }, { g1[2], g2[2], a3[2] } };
                        double[,] G_i = new double[3, 3] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } };

                        (double[,] Aijkl_3D, double[] FPK_3D_vec) = transformations.CalculateTransformations(tgi, G_i, F_3D);

                        for (int r1 = 0; r1 < 3; r1++)
                        {
                            //(31)
                            Vector dg1_dr = a1r.GetColumn(r1) + da3_dksidr[r1] * z;
                            Vector dg2_dr = a2r.GetColumn(r1) + da3_dhetadr[r1] * z;

                            //Vector dg3_dr = a3r. ....

                            //(39)
                            double[,] dF_3D_dr = new double[3, 3] { { dg1_dr[0]*G_1[0]+dg2_dr[0]*G_2[0], dg1_dr[0]*G_1[1]+dg2_dr[0]*G_2[1], dg1_dr[0]*G_1[2]+dg2_dr[0]*G_2[2] },
                                                                 { dg1_dr[1]*G_1[0]+dg2_dr[1]*G_2[0], dg1_dr[1]*G_1[1]+dg2_dr[1]*G_2[1], dg1_dr[1]*G_1[2]+dg2_dr[1]*G_2[2] },
                                                                 { dg1_dr[2]*G_1[0]+dg2_dr[2]*G_2[0], dg1_dr[2]*G_1[1]+dg2_dr[2]*G_2[1], dg1_dr[2]*G_1[2]+dg2_dr[2]*G_2[2] }, };

                            double[] dF_3D_dr_vec = { dF_3D_dr[0, 0], dF_3D_dr[1, 1], dF_3D_dr[2, 2], dF_3D_dr[0, 1], dF_3D_dr[1, 2], dF_3D_dr[2, 0], dF_3D_dr[0, 2], dF_3D_dr[1, 0], dF_3D_dr[2, 1] };

                            if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0))
                            { ElementStiffnesses.ProccessVariable(18, dF_3D_dr_vec, true, 3 * i + r1); }

                            for (int i2 = 0; i2 < 9; i2++)
                            {
                                forceIntegration[i1, r1] += FPK_3D_vec[i2] * dF_3D_dr_vec[i2] * wfactor;
                                thicknesCoeffs[i1, r1] = w;

                                forcesDevelop[3 * i + r1] += FPK_3D_vec[i2] * dF_3D_dr_vec[i2] * wfactor * w;
                                forcesDevelopGp[3 * i + r1] += FPK_3D_vec[i2] * dF_3D_dr_vec[i2] * wfactor * w;
                            }


                        }
                        //  (31) 



                    }
                }

                #endregion

            }

            if ((ElementStiffnesses.saveForcesState1 | ElementStiffnesses.saveForcesState2s) | ElementStiffnesses.saveForcesState0)
            {
                ElementStiffnesses.SaveNodalForces(elementNodalMembraneForces, elementNodalBendingForces, element);
            }

            if (runNewForces)
            { return forcesDevelop; }
            else
            {
                return elementNodalForces;
            }



        }

        private (Vector[] da3_dksidr, Vector[] da3_dhetadr) Calculate_da3_dksidr(Vector[] da3tilde_dksidr, Vector[] da3tilde_dhetadr, Vector da3tilde_dksi,
            Vector da3tilde_dheta, double[] dnorma3_dr, Vector a3_tilde, double[] da3norm_dksidr, double[] da3norm_dhetadr, double da3norm_dksi,
            double da3norm_dheta, double J1, double[][] da3tilde_dr)
        {
            Vector[] da3_dksidr = new Vector[3];
            Vector[] da3_dhetadr = new Vector[3];

            for (int r1 = 0; r1 < 3; r1++)
            {
                var firstVec_0 = da3tilde_dksidr[r1][0] / J1;
                var firstVec_1 = da3tilde_dksidr[r1][1] / J1;
                var firstVec_2 = da3tilde_dksidr[r1][2] / J1;

                double scale2 = -((double)1 / (Math.Pow(J1, 2))); //denominator of vectors 2 3 and 4 and a minus.

                var scale3 = dnorma3_dr[r1] * scale2;
                var secondVec_0 = da3tilde_dksi[0] * scale3;
                var secondVec_1 = da3tilde_dksi[1] * scale3;
                var secondVec_2 = da3tilde_dksi[2] * scale3;

                var scale4 = da3norm_dksi * scale2;
                var thirdVec_0 = da3tilde_dr[r1][0] * scale4;
                var thirdVec_1 = da3tilde_dr[r1][1] * scale4;
                var thirdVec_2 = da3tilde_dr[r1][2] * scale4;

                var scale6 = da3norm_dksidr[r1] * scale2;
                var fourthVec_0 = a3_tilde[0] * scale6;
                var fourthVec_1 = a3_tilde[1] * scale6;
                var fourthVec_2 = a3_tilde[2] * scale6;

                double scale5 = ((double)1) / Math.Pow(J1, 3);
                var scale7 = 2 * da3norm_dksi * dnorma3_dr[r1] * scale5;
                var fifthvector_0 = a3_tilde[0] * scale7;
                var fifthvector_1 = a3_tilde[1] * scale7;
                var fifthvector_2 = a3_tilde[2] * scale7;

                da3_dksidr[r1] = Vector.CreateFromArray(new double[]
                {
                    firstVec_0 + secondVec_0 + thirdVec_0 + fourthVec_0 + fifthvector_0,
                    firstVec_1 + secondVec_1 + thirdVec_1 + fourthVec_1 + fifthvector_1,
                    firstVec_2 + secondVec_2 + thirdVec_2 + fourthVec_2 + fifthvector_2,
                });

            }

            for (int r1 = 0; r1 < 3; r1++)
            {
                var firstVec_0 = da3tilde_dhetadr[r1][0] / J1;
                var firstVec_1 = da3tilde_dhetadr[r1][1] / J1;
                var firstVec_2 = da3tilde_dhetadr[r1][2] / J1;

                double scale2 = -((double)1 / (Math.Pow(J1, 2))); //denominator of vectors 2 3 and 4 and a minus.

                var scale3 = dnorma3_dr[r1] * scale2;
                var secondVec_0 = da3tilde_dheta[0] * scale3;
                var secondVec_1 = da3tilde_dheta[1] * scale3;
                var secondVec_2 = da3tilde_dheta[2] * scale3;

                var scale4 = da3norm_dheta * scale2;
                var thirdVec_0 = da3tilde_dr[r1][0] * scale4;
                var thirdVec_1 = da3tilde_dr[r1][1] * scale4;
                var thirdVec_2 = da3tilde_dr[r1][2] * scale4;

                var scale6 = da3norm_dhetadr[r1] * scale2;
                var fourthVec_0 = a3_tilde[0] * scale6;
                var fourthVec_1 = a3_tilde[1] * scale6;
                var fourthVec_2 = a3_tilde[2] * scale6;

                double scale5 = ((double)1) / Math.Pow(J1, 3);
                var scale7 = 2 * da3norm_dheta * dnorma3_dr[r1] * scale5;
                var fifthvector_0 = a3_tilde[0] * scale7;
                var fifthvector_1 = a3_tilde[1] * scale7;
                var fifthvector_2 = a3_tilde[2] * scale7;

                da3_dhetadr[r1] = Vector.CreateFromArray(new double[]
                {
                    firstVec_0 + secondVec_0 + thirdVec_0 + fourthVec_0 + fifthvector_0,
                    firstVec_1 + secondVec_1 + thirdVec_1 + fourthVec_1 + fifthvector_1,
                    firstVec_2 + secondVec_2 + thirdVec_2 + fourthVec_2 + fifthvector_2,
                });


            }

            return (da3_dksidr, da3_dhetadr);
        }

        private (double[] da3norm_dksidr, double[] da3norm_dhetadr) Calculate_da3norm_dksidr(Vector[] da3tilde_dksidr, Vector[] da3tilde_dhetadr, Vector a3_tilde, Vector da3tilde_dksi, Vector da3tilde_dheta, double[][] da3tilde_dr, double J1)
        {
            var da3norm_dksidr = new double[3];
            var da3norm_dhetadr = new double[3];

            for (int r1 = 0; r1 < 3; r1++)
            {
                double firstNumerator = da3tilde_dksidr[r1].DotProduct(a3_tilde) + da3tilde_dksi.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double firstDenominator = J1;

                double secondNumerator = da3tilde_dksi.DotProduct(a3_tilde) * da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray());
                double secondDenominator = Math.Pow(J1, 3);

                da3norm_dksidr[r1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);

            }

            for (int r1 = 0; r1 < 3; r1++)
            {
                double firstNumerator = da3tilde_dhetadr[r1].DotProduct(a3_tilde) + da3tilde_dheta.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double firstDenominator = J1;

                double secondNumerator = da3tilde_dheta.DotProduct(a3_tilde) * da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray());
                double secondDenominator = Math.Pow(J1, 3);

                da3norm_dhetadr[r1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);

            }

            return (da3norm_dksidr, da3norm_dhetadr);
        }

        private (Vector[] da3tilde_dksidr, Vector[] da3tilde_dhetadr) Calculate_da3tilde_dksidr(Matrix3by3 a1r, Matrix3by3 a2r, Matrix3by3 a11r, Matrix3by3 a22r, Matrix3by3 a12r,
            Vector a1, Vector a2, Vector a11, Vector a22, Vector a12)
        {
            Vector[] da3tilde_dksidr = new Vector[3];
            Vector[] da3tilde_dhetadr = new Vector[3];

            for (int r1 = 0; r1 < 3; r1++)
            {
                da3tilde_dksidr[r1] = a11r.GetColumn(r1).CrossProduct(a2) + a11.CrossProduct(a2r.GetColumn(r1)) + a1r.GetColumn(r1).CrossProduct(a12) + a1.CrossProduct(a12r.GetColumn(r1));
                da3tilde_dhetadr[r1] = a12r.GetColumn(r1).CrossProduct(a2) + a12.CrossProduct(a2r.GetColumn(r1)) + a1r.GetColumn(r1).CrossProduct(a22) + a1.CrossProduct(a22r.GetColumn(r1));
            }

            return (da3tilde_dksidr, da3tilde_dhetadr);
        }

        private (Vector da3tilde_dksi, Vector da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, Vector da3_dksi, Vector da3_dheta) 
            Calculate_da3tilde_dksi_524_525_526_b(Vector a1, Vector a2, Vector a11, Vector a22, Vector a12, Vector a3, double normA3)
        {
            var da3tilde_dksi = a11.CrossProduct(a2) + a1.CrossProduct(a12);
            var da3tilde_dheta = a12.CrossProduct(a2) + a1.CrossProduct(a22);

            var da3norm_dksi = a3.DotProduct(da3tilde_dksi);
            var da3norm_dheta = a3.DotProduct(da3tilde_dheta);

            double scaleFactor = (double)1 / normA3;
            var da3_dksi = da3tilde_dksi.Scale(scaleFactor) - a3.Scale(da3norm_dksi).Scale(scaleFactor);
            var da3_dheta = da3tilde_dheta.Scale(scaleFactor) - a3.Scale(da3norm_dheta).Scale(scaleFactor);


            return (da3tilde_dksi, da3tilde_dheta, da3norm_dksi, da3norm_dheta, da3_dksi, da3_dheta);

        }


        #region develop formulation

        #endregion




        public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            NeumannBoundaryCondition neumann) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            NeumannBoundaryCondition neumann) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge,
            PressureBoundaryCondition pressure) => throw new NotImplementedException();

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face,
            PressureBoundaryCondition pressure) => throw new NotImplementedException();

        double[,] jacobianMatrix = new double[2, 3];
        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var nurbs = CalculateShapeFunctions(shellElement, elementControlPoints);

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(elementControlPoints);
            var midsurfaceGP = materialsAtThicknessGP.Keys.ToArray();

            for (var j = 0; j < midsurfaceGP.Length; j++)
            {
                CalculateJacobian(newControlPoints, nurbs, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, nurbs, j);

                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

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
                
                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                if (j == ElementStiffnesses.gpNumberToCheck)
                {
                    ElementStiffnesses.ProccessVariable(1, surfaceBasisVectorDerivative1, false);
                    ElementStiffnesses.ProccessVariable(2, surfaceBasisVectorDerivative2, false);
                    ElementStiffnesses.ProccessVariable(3, surfaceBasisVectorDerivative12, false);

                    ElementStiffnesses.ProccessVariable(6, surfaceBasisVector1, false);
                    ElementStiffnesses.ProccessVariable(7, surfaceBasisVector2, false);
                }

                var A11 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors1[j][0] + //E
                          initialSurfaceBasisVectors1[j][1] * initialSurfaceBasisVectors1[j][1] +
                          initialSurfaceBasisVectors1[j][2] * initialSurfaceBasisVectors1[j][2];

                var A22 = initialSurfaceBasisVectors2[j][0] * initialSurfaceBasisVectors2[j][0] + //G
                          initialSurfaceBasisVectors2[j][1] * initialSurfaceBasisVectors2[j][1] +
                          initialSurfaceBasisVectors2[j][2] * initialSurfaceBasisVectors2[j][2];

                var A12 = initialSurfaceBasisVectors1[j][0] * initialSurfaceBasisVectors2[j][0] + //F
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


                var membraneStrain = new double[] {0.5 * (a11 - A11), 0.5 * (a22 - A22), a12 - A12};

                if (j == ElementStiffnesses.gpNumberToCheck)
                {
                    ElementStiffnesses.ProccessVariable(4, membraneStrain, false);
                }

                var B11 = initialSurfaceBasisVectorDerivative1[j][0] * initialUnitSurfaceBasisVectors3[j][0] + //L
                          initialSurfaceBasisVectorDerivative1[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
                          initialSurfaceBasisVectorDerivative1[j][2] * initialUnitSurfaceBasisVectors3[j][2];

                var B22 = initialSurfaceBasisVectorDerivative2[j][0] * initialUnitSurfaceBasisVectors3[j][0] + //N
                          initialSurfaceBasisVectorDerivative2[j][1] * initialUnitSurfaceBasisVectors3[j][1] +
                          initialSurfaceBasisVectorDerivative2[j][2] * initialUnitSurfaceBasisVectors3[j][2];

                var B12 = initialSurfaceBasisVectorDerivative12[j][0] * initialUnitSurfaceBasisVectors3[j][0] + //M
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

                //var bendingStrain = new double[] {b11 - B11, b22 - B22, 2 * b12 - 2 * B12};
                var bendingStrain = new double[] { -(b11 - B11), -(b22 - B22), -(2 * b12 - 2 * B12) };

                //double du = 0.1;
                //double dv = 0.1;
                //double kn = (B11 * du * du + 2 * B12 * du * dv + B22 * dv * dv) / (A11 * du * du + 2 * A12 * du * dv + A22 * dv * dv);
                if (j == ElementStiffnesses.gpNumberToCheck)
                {
                    ElementStiffnesses.ProccessVariable(5, bendingStrain, false);
                }

                

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

        public Dictionary<int, double> CalculateSurfaceDistributedLoad(Element element, IDofType loadedDof,
            double loadMagnitude)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var distributedLoad = new Dictionary<int, double>();
            var nurbs = new Nurbs2D(shellElement, elementControlPoints);

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(elementControlPoints, nurbs, j, jacobianMatrix);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                var J1 = surfaceBasisVector3.Norm2();
                surfaceBasisVector3.ScaleIntoThis(1 / J1);

                for (int i = 0; i < elementControlPoints.Length; i++)
                {
                    var loadedDofIndex = ControlPointDofTypes.FindFirstIndex(loadedDof);
                    if (!element.Model.GlobalDofOrdering.GlobalFreeDofs.Contains(elementControlPoints[i], loadedDof))
                        continue;
                    var dofId = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[i], loadedDof];

                    if (distributedLoad.ContainsKey(dofId))
                    {
                        distributedLoad[dofId] += loadMagnitude * J1 *
                                                  nurbs.NurbsValues[i, j] * gaussPoints[j].WeightFactor;
                    }
                    else
                    {
                        distributedLoad.Add(dofId,
                            loadMagnitude * nurbs.NurbsValues[i, j] * J1 * gaussPoints[j].WeightFactor);
                    }
                }
            }

            return distributedLoad;
        }

        public Dictionary<int, double> CalculateSurfacePressure(Element element, double pressureMagnitude)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var elementControlPoints = shellElement.ControlPoints.ToArray();
            var gaussPoints = CreateElementGaussPoints(shellElement);
            var pressureLoad = new Dictionary<int, double>();
            var nurbs = new Nurbs2D(shellElement, elementControlPoints);

            for (var j = 0; j < gaussPoints.Count; j++)
            {
                CalculateJacobian(elementControlPoints, nurbs, j, jacobianMatrix);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
                var surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                var J1 = surfaceBasisVector3.Norm2();
                surfaceBasisVector3.ScaleIntoThis(1 / J1);

                for (int i = 0; i < elementControlPoints.Length; i++)
                {
                    for (int k = 0; k < ControlPointDofTypes.Length; k++)
                    {
                        int dofId = element.Model.GlobalDofOrdering.GlobalFreeDofs[elementControlPoints[i],
                            ControlPointDofTypes[k]];

                        if (pressureLoad.ContainsKey(dofId))
                        {
                            pressureLoad[dofId] += pressureMagnitude * surfaceBasisVector3[k] *
                                                   nurbs.NurbsValues[i, j] * gaussPoints[j].WeightFactor;
                        }
                        else
                        {
                            pressureLoad.Add(dofId,
                                pressureMagnitude * surfaceBasisVector3[k] * nurbs.NurbsValues[i, j] *
                                gaussPoints[j].WeightFactor);
                        }
                    }
                }
            }

            return pressureLoad;
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
            MembraneForces= new Forces();
            BendingMoments= new Forces();
            var thicknessPoints = thicknessIntegrationPoints[midSurfaceGaussPoint];

            for (int i = 0; i < thicknessPoints.Count; i++)
            {
                var thicknessPoint = thicknessPoints[i];
                var material = materialsAtThicknessGP[midSurfaceGaussPoint][thicknessPoints[i]];
                var w = thicknessPoint.WeightFactor;
                var z = thicknessPoint.Zeta;
                MembraneForces.v0+= material.Stresses[0] * w;
                MembraneForces.v1 += material.Stresses[1] * w;
                MembraneForces.v2 += material.Stresses[2] * w;

                BendingMoments.v0 += material.Stresses[0] * w * z;
                BendingMoments.v1 += material.Stresses[1] * w * z;
                BendingMoments.v2 += material.Stresses[2] * w * z;
            }
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

        private Nurbs2D _nurbs;

        public IMatrix StiffnessMatrix(IElement element)
        {
            var shellElement = (NurbsKirchhoffLoveShellElementNLDevelop) element;
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();
            var nurbs = CalculateShapeFunctions(shellElement, _controlPoints);

            if (!isInitialized)
            {
                CalculateInitialConfigurationData(_controlPoints, nurbs, gaussPoints);
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

            var KmembraneNL_total = new double[bCols, bCols];
            var KbendingNL_total = new double[bCols, bCols];

            //var KL_total = new double[bCols, bCols];

            for (int j = 0; j < gaussPoints.Length; j++)
            {
                ElementStiffnesses.gpNumber = j;

                CalculateJacobian(elementControlPoints, nurbs, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(elementControlPoints, nurbs, j);
                var surfaceBasisVector1 = CalculateSurfaceBasisVector1(jacobianMatrix, 0);

                var surfaceBasisVector2 = CalculateSurfaceBasisVector1(jacobianMatrix, 1);

                var surfaceBasisVector3 = new[]
                {
                    surfaceBasisVector1[1] * surfaceBasisVector2[2] - surfaceBasisVector1[2] * surfaceBasisVector2[1],
                    surfaceBasisVector1[2] * surfaceBasisVector2[0] - surfaceBasisVector1[0] * surfaceBasisVector2[2],
                    surfaceBasisVector1[0] * surfaceBasisVector2[1] - surfaceBasisVector1[1] * surfaceBasisVector2[0],
                };

                var J1 = Math.Sqrt(surfaceBasisVector3[0]* surfaceBasisVector3[0]+
                                   surfaceBasisVector3[1] * surfaceBasisVector3[1]+
                                   surfaceBasisVector3[2] * surfaceBasisVector3[2]);

                surfaceBasisVector3[0] /= J1;
                surfaceBasisVector3[1] /= J1;
                surfaceBasisVector3[2] /= J1;

                var surfaceBasisVectorDerivative1 = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                var surfaceBasisVectorDerivative2 = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                var surfaceBasisVectorDerivative12 = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                double wFactor = InitialJ1[j] * gaussPoints[j].WeightFactor;

                CalculateLinearStiffness(elementControlPoints, nurbs, j, surfaceBasisVector1, surfaceBasisVector2,
                    Bmembrane, surfaceBasisVector3, surfaceBasisVectorDerivative1, J1, surfaceBasisVectorDerivative2,
                    surfaceBasisVectorDerivative12, Bbending, gaussPoints, BmTranspose, bRows, bCols, BbTranspose,
                    wFactor, BmTransposeMultStiffness, BbTransposeMultStiffness, BmbTransposeMultStiffness,
                    BbmTransposeMultStiffness, stiffnessMatrix, KmembraneL, KbendingL);
                CalculateNonLinearStiffness(gaussPoints, j, KmembraneNL, bCols, KbendingNL, elementControlPoints, nurbs,
                    surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, surfaceBasisVectorDerivative1,
                    surfaceBasisVectorDerivative2, surfaceBasisVectorDerivative12, J1, stiffnessMatrix, wFactor,
                    ref MembraneForces, ref BendingMoments);
                if (ElementStiffnesses.saveStiffnessMatrixState)
                {
                    for (var i = 0; i < stiffnessMatrix.GetLength(0); i++)
                    {
                        for (var k = 0; k < stiffnessMatrix.GetLength(1); k++)
                        {


                            KmembraneNL_total[i, k] += KmembraneNL[i, k] * wFactor;
                            KbendingNL_total[i, k] += KbendingNL[i, k] * wFactor;
                        }


                    }
                }

            }

            if (ElementStiffnesses.saveStiffnessMatrixState)
            {
                ElementStiffnesses.SaveStiffnessMatrixes((Matrix.CreateFromArray(stiffnessMatrix)- Matrix.CreateFromArray(KmembraneNL_total) - Matrix.CreateFromArray(KbendingNL_total)).CopytoArray2D(),
                    KmembraneNL_total, KbendingNL_total, Matrix.CreateFromArray(KbendingNL_total) + Matrix.CreateFromArray(KmembraneNL_total),
                     KmembraneL, KbendingL,element);
            }
            
            return Matrix.CreateFromArray(stiffnessMatrix);
        }

        private void CalculateNonLinearStiffness(GaussLegendrePoint3D[] gaussPoints, int j, double[,] KmembraneNL, int bCols,
            double[,] KbendingNL, ControlPoint[] elementControlPoints, Nurbs2D nurbs, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3, double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12, double J1,
            double[,] stiffnessMatrix, double wFactor, ref Forces MembraneForces, ref Forces BendingMoments)
        {
            IntegratedStressesOverThickness(gaussPoints[j], ref MembraneForces, ref BendingMoments);

            Array.Clear(KmembraneNL, 0, bCols * bCols);
            Array.Clear(KbendingNL, 0, bCols * bCols);

            CalculateKmembraneNL(elementControlPoints, ref MembraneForces, nurbs, j, KmembraneNL);
            CalculateKbendingNL(elementControlPoints, ref BendingMoments, nurbs,
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

        private void CalculateLinearStiffness(ControlPoint[] elementControlPoints, Nurbs2D nurbs, int j,
            double[] surfaceBasisVector1, double[] surfaceBasisVector2, double[,] Bmembrane, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double J1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double[,] Bbending, GaussLegendrePoint3D[] gaussPoints,
            double[,] BmTranspose, int bRows, int bCols, double[,] BbTranspose, double wFactor,
            double[,] BmTransposeMultStiffness, double[,] BbTransposeMultStiffness, double[,] BmbTransposeMultStiffness,
            double[,] BbmTransposeMultStiffness, double[,] stiffnessMatrix, double[,] KmembraneL, double[,] KbendingL)
        {
            CalculateMembraneDeformationMatrix(elementControlPoints.Length, nurbs, j, surfaceBasisVector1,
                surfaceBasisVector2, Bmembrane);
            CalculateBendingDeformationMatrix(elementControlPoints.Length, surfaceBasisVector3, nurbs, j,
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
                for (int k = 0; k < bRows; k++) // to k sumvolizei to athroisma 
                {
                    tempm = BmTranspose[i, k]; // kai tp , m sumvolizei th thesi pou gemizetai
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
                        if(ElementStiffnesses.saveStiffnessMatrixState)
                        {
                            KbendingL[i, m] += tempb * ben + tempbm * mem;
                            KmembraneL[i, m] += tempm * mem + tempmb * ben;
                        }
                    }
                }
            }
        }

        private Nurbs2D CalculateShapeFunctions(NurbsKirchhoffLoveShellElementNLDevelop shellElement,
            ControlPoint[] controlPoints)
        {
            return _nurbs ?? (_nurbs = new Nurbs2D(shellElement, controlPoints));
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

        internal double[,] CalculateHessian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j)
        {
            var hessianMatrix = new double[3, 3];
            for (var k = 0; k < controlPoints.Length; k++)
            {
                hessianMatrix[0, 0] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].X;
                hessianMatrix[0, 1] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Y;
                hessianMatrix[0, 2] +=
                    nurbs.NurbsSecondDerivativeValueKsi[k, j] * controlPoints[k].Z;
                hessianMatrix[1, 0] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].X;
                hessianMatrix[1, 1] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[1, 2] +=
                    nurbs.NurbsSecondDerivativeValueHeta[k, j] * controlPoints[k].Z;
                hessianMatrix[2, 0] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].X;
                hessianMatrix[2, 1] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Y;
                hessianMatrix[2, 2] +=
                    nurbs.NurbsSecondDerivativeValueKsiHeta[k, j] * controlPoints[k].Z;
            }

            return hessianMatrix;
        }

        internal void CalculateJacobian(ControlPoint[] controlPoints, Nurbs2D nurbs, int j, double[,] jacobianOut)
        {
            jacobianOut[0, 0] = jacobianOut[0, 1] = jacobianOut[0, 2] =
                jacobianOut[1, 0] = jacobianOut[1, 1] = jacobianOut[1, 2] = 0.0;
            for (var k = 0; k < controlPoints.Length; k++)
            {
                jacobianOut[0, 0] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].X;
                jacobianOut[0, 1] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].Y;
                jacobianOut[0, 2] += nurbs.NurbsDerivativeValuesKsi[k, j] * controlPoints[k].Z;
                jacobianOut[1, 0] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].X;
                jacobianOut[1, 1] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].Y;
                jacobianOut[1, 2] += nurbs.NurbsDerivativeValuesHeta[k, j] * controlPoints[k].Z;
            }
        }

        internal double[] CalculateSurfaceBasisVector1(double[,] Matrix, int row)
        {
            var surfaceBasisVector1 = new double[3];
            surfaceBasisVector1[0] = Matrix[row, 0];
            surfaceBasisVector1[1] = Matrix[row, 1];
            surfaceBasisVector1[2] = Matrix[row, 2];
            return surfaceBasisVector1;
        }

        internal void CalculateBendingDeformationMatrix(int controlPointsCount, double[] surfaceBasisVector3,
            Nurbs2D nurbs, int j, double[] surfaceBasisVector2, double[] surfaceBasisVectorDerivative1,
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
                var dksi = nurbs.NurbsDerivativeValuesKsi[column / 3, j];
                var dheta = nurbs.NurbsDerivativeValuesHeta[column / 3, j];
                var d2Ksi = nurbs.NurbsSecondDerivativeValueKsi[column / 3, j];
                var d2Heta = nurbs.NurbsSecondDerivativeValueHeta[column / 3, j];
                var d2KsiHeta = nurbs.NurbsSecondDerivativeValueKsiHeta[column / 3, j];

                BbendingOut[0, column] = -d2Ksi * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s11 * s11_2 - s12 * s11_1) + dksi * (s21 * s11_2 - s22 * s11_1)) / J1;
                BbendingOut[0, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_2 - s12 * s11_0) + dksi * (s20 * s11_2 - s22 * s11_0)) / J1 - d2Ksi * s31;
                BbendingOut[0, column + 2] = -d2Ksi * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s11_0 * s30 + s11_1 * s31 + s11_2 * s32) - dheta * (s10 * s11_1 - s11 * s11_0) + dksi * (s20 * s11_1 - s21 * s11_0)) / J1;

                BbendingOut[1, column] = -d2Heta * s30 - ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s11 * s22_2 - s12 * s22_1) + dksi * (s21 * s22_2 - s22 * s22_1)) / J1;
                BbendingOut[1, column + 1] = ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_2 - s12 * s22_0) + dksi * (s20 * s22_2 - s22 * s22_0)) / J1 - d2Heta * s31;
                BbendingOut[1, column + 2] = -d2Heta * s32 - ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s22_0 * s30 + s22_1 * s31 + s22_2 * s32) - dheta * (s10 * s22_1 - s11 * s22_0) + dksi * (s20 * s22_1 - s21 * s22_0)) / J1;

                BbendingOut[2, column] = -2 * d2KsiHeta * s30 - (2 * ((dheta * (s11 * s32 - s12 * s31) - dksi * (s21 * s32 - s22 * s31)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s11 * s12_2 - s12 * s12_1) + dksi * (s21 * s12_2 - s22 * s12_1))) / J1;
                BbendingOut[2, column + 1] = (2 * ((dheta * (s10 * s32 - s12 * s30) - dksi * (s20 * s32 - s22 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_2 - s12 * s12_0) + dksi * (s20 * s12_2 - s22 * s12_0))) / J1 - 2 * d2KsiHeta * s31;
                BbendingOut[2, column + 2] = -2 * d2KsiHeta * s32 - (2 * ((dheta * (s10 * s31 - s11 * s30) - dksi * (s20 * s31 - s21 * s30)) * (s12_0 * s30 + s12_1 * s31 + s12_2 * s32) - dheta * (s10 * s12_1 - s11 * s12_0) + dksi * (s20 * s12_1 - s21 * s12_0))) / J1;

                //BbendingOut[0, column] = -BbendingOut[0, column];
                //BbendingOut[0, column + 1] = -BbendingOut[0, column + 1];
                //BbendingOut[0, column + 2] = -BbendingOut[0, column + 2];

                //BbendingOut[1, column] = -BbendingOut[1, column];
                //BbendingOut[1, column + 1] = -BbendingOut[1, column + 1];
                //BbendingOut[1, column + 2] = -BbendingOut[1, column + 2];

                //BbendingOut[2, column] = -BbendingOut[2, column];
                //BbendingOut[2, column + 1] = -BbendingOut[2, column + 1];
                //BbendingOut[2, column + 2] = -BbendingOut[2, column + 2];
            }
        }

        internal static void CalculateCrossProduct(double[] vector1, double[] vector2, double[] result )
        {
            result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
            result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
            result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
        }

        private double[][] initialSurfaceBasisVectors1;
        private double[][] initialSurfaceBasisVectors2;
        private double[][] initialUnitSurfaceBasisVectors3;

        private double[][] initialSurfaceBasisVectorDerivative1;
        private double[][] initialSurfaceBasisVectorDerivative2;
        private double[][] initialSurfaceBasisVectorDerivative12;

        private double[] InitialJ1;

        internal void CalculateInitialConfigurationData(ControlPoint[] controlPoints,
            Nurbs2D nurbs, IList<GaussLegendrePoint3D> gaussPoints)
        {
            bool useDifferenetTicknessCovariants = true; 
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
                CalculateJacobian(controlPoints, nurbs, j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(controlPoints, nurbs, j);
                initialSurfaceBasisVectors1[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 0);
                initialSurfaceBasisVectors2[j] = CalculateSurfaceBasisVector1(jacobianMatrix, 1);
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

                initialSurfaceBasisVectorDerivative1[j] = CalculateSurfaceBasisVector1(hessianMatrix, 0);
                initialSurfaceBasisVectorDerivative2[j] = CalculateSurfaceBasisVector1(hessianMatrix, 1);
                initialSurfaceBasisVectorDerivative12[j] = CalculateSurfaceBasisVector1(hessianMatrix, 2);

                var a1_init = Vector.CreateFromArray(initialSurfaceBasisVectors1[j]);
                var a2_init = Vector.CreateFromArray(initialSurfaceBasisVectors2[j]);
                var a3_init = Vector.CreateFromArray(initialUnitSurfaceBasisVectors3[j]);
                var a11_init = Vector.CreateFromArray(initialSurfaceBasisVectorDerivative1[j]);
                var a22_init = Vector.CreateFromArray(initialSurfaceBasisVectorDerivative2[j]);
                var a12_init = Vector.CreateFromArray(initialSurfaceBasisVectorDerivative12[j]);

                (Vector da3tilde_dksi_init, Vector da3tilde_dheta_init, double da3norm_dksi_init, double da3norm_dheta_init, Vector da3_dksi_init, Vector da3_dheta_init) =
            Calculate_da3tilde_dksi_524_525_526_b(a1_init, a2_init, a11_init, a22_init, a12_init, a3_init, InitialJ1[j]);

                var thicknessPoints = thicknessIntegrationPoints[gaussPoints[j]];
                for (int i1 = 0; i1 < thicknessPoints.Count ; i1++)
                {
                    var integrationPointMaterial = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    //}
                    //foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
                    //{

                    var thicknessPoint = thicknessPoints[i1];
                    if (!useDifferenetTicknessCovariants)
                    {
                        integrationPointMaterial.TangentVectorV1 = initialSurfaceBasisVectors1[j];
                        integrationPointMaterial.TangentVectorV2 = initialSurfaceBasisVectors2[j];
                        integrationPointMaterial.NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
                    }
                    else
                    {
                        var z = thicknessPoint.Zeta;
                        Vector G1 = a1_init + da3_dksi_init.Scale(z);
                        Vector G2 = a2_init + da3_dheta_init.Scale(z);

                        integrationPointMaterial.TangentVectorV1 = G1.CopyToArray();
                        integrationPointMaterial.TangentVectorV2 = G2.CopyToArray();
                        integrationPointMaterial.NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
                    }
                }
            }
        }

        private a3rs a3rs=new a3rs();
        private Bab_rs Bab_rs = new Bab_rs();
        private a3r a3r= new a3r();
        private a3r a3s= new a3r();
        private ControlPoint[] _controlPoints;

        internal void CalculateKbendingNL(ControlPoint[] controlPoints,
            ref Forces bendingMoments, Nurbs2D nurbs, double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double[] surfaceBasisVectorDerivative1, double[] surfaceBasisVectorDerivative2,
            double[] surfaceBasisVectorDerivative12, double J1, int j, double[,] KbendingNLOut)
        {
            var a3rArray = new a3r[controlPoints.Length];
            for (var i = 0; i < controlPoints.Length; i++)
            {
                var a3r= new a3r();
                var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];
                CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                    dheta_r, J1, ref a3r);
                a3rArray[i] = a3r;
            }
            
            for (int i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];
                var d2Ksi_dr2 = nurbs.NurbsSecondDerivativeValueKsi[i, j];
                var d2Heta_dr2 = nurbs.NurbsSecondDerivativeValueHeta[i, j];
                var d2KsiHeta_dr2 = nurbs.NurbsSecondDerivativeValueKsiHeta[i, j];

                var a3r = a3rArray[i];

                if (ElementStiffnesses.performCalculations)
                {
                    var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                    var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                    var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                    var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                    var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);
                    if (ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck)
                    {
                        for (int i1 = 0; i1 < 3; i1++)
                        {
                            ElementStiffnesses.ProccessVariable(1, a11r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(2, a22r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                            ElementStiffnesses.ProccessVariable(3, a12r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);

                        }
                    }
                }

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var d2Ksi_ds2 = nurbs.NurbsSecondDerivativeValueKsi[k, j];
                    var d2Heta_ds2 = nurbs.NurbsSecondDerivativeValueHeta[k, j];
                    var d2KsiHeta_ds2 = nurbs.NurbsSecondDerivativeValueKsiHeta[k, j];

                    var dksi_s = nurbs.NurbsDerivativeValuesKsi[k, j];
                    var dheta_s = nurbs.NurbsDerivativeValuesHeta[k, j];
                    
                    var a3s = a3rArray[k];
                    a3rs=new a3rs();//Clear struct values
                    //var a1s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[k, j]);
                    //var a2s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[k, j]);
                    (a3rs a3rsAlternative, var da3tilde_drds, var da3tilde_dr, var da3tilde_ds,
                        double[] dnorma3_dr, double[] dnorma3_ds, double[,] dnorma3_drds, var a3_tilde, var da3_drds) =
                        Calculate_a3rs(Vector.CreateFromArray(surfaceBasisVector1), Vector.CreateFromArray(surfaceBasisVector2), 
                            Vector.CreateFromArray(surfaceBasisVector3),J1, dksi_r,dksi_s, dheta_r, dheta_s);
                    
                    Bab_rs Bab_rsAlternative = CalculateBab_rs(surfaceBasisVectorDerivative1, surfaceBasisVectorDerivative2,
                        surfaceBasisVectorDerivative12,d2Ksi_dr2,d2Ksi_ds2, d2Heta_dr2,d2Heta_ds2, d2KsiHeta_dr2,d2KsiHeta_ds2,
                        a3rsAlternative, a3r, a3s, da3_drds);
                    Bab_rs = Bab_rsAlternative;

                    if (ElementStiffnesses.performCalculations)
                    {
                        if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveStiffnessMatrixState)
                        {
                            ElementStiffnesses.saveOriginalState = true;
                            ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 }, true, 3 * i+ 0);
                            ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r01, a3r.a3r11, a3r.a3r21 }, true, 3 * i + 1);
                            ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r02, a3r.a3r12, a3r.a3r22 }, true, 3 * i + 2); ;

                            if (i == 0) // ennoume thn paragwgo ws pros r gia k=0 kai gai to x dof
                            {
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs00_0, a3rs.a3rs00_1, a3rs.a3rs00_2 }, true, 3 * k + 0);
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs01_0, a3rs.a3rs01_1, a3rs.a3rs01_2 }, true, 3 * k + 1);
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs02_0, a3rs.a3rs02_1, a3rs.a3rs02_2 }, true, 3 * k + 2);
                                ElementStiffnesses.ProccessVariable(9, new double[] { a3rsAlternative.a3rs00_0, a3rsAlternative.a3rs00_1, a3rsAlternative.a3rs00_2 }, true, 3 * k + 0);
                                ElementStiffnesses.ProccessVariable(9, new double[] { a3rsAlternative.a3rs01_0, a3rsAlternative.a3rs01_1, a3rsAlternative.a3rs01_2 }, true, 3 * k + 1);
                                ElementStiffnesses.ProccessVariable(9, new double[] { a3rsAlternative.a3rs02_0, a3rsAlternative.a3rs02_1, a3rsAlternative.a3rs02_2 }, true, 3 * k + 2);
                                ElementStiffnesses.ProccessVariable(9, new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 }, false);

                            }


                            ElementStiffnesses.saveOriginalState = false;

                        }

                        if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveVariationStates)
                        {
                            //ElementStiffnesses.saveOriginalState = true;
                            //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r00, a3r.a3r01, a3r.a3r02 }, true, 3 * i + 0);
                            //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r10, a3r.a3r11, a3r.a3r12 }, true, 3 * i + 1);
                            //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r20, a3r.a3r21, a3r.a3r22 }, true, 3 * i + 2);

                            //ElementStiffnesses.ProccessVariable(8, surfaceBasisVector3, false);

                            if (i == 0) // ennoume thn paragwgo ws pros r gia k=0 kai gai to x dof
                            {
                                ElementStiffnesses.ProccessVariable(9,new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 }, false);
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs00_0, a3rs.a3rs00_1, a3rs.a3rs00_2 }, true, 3 * k + 0);
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs01_0, a3rs.a3rs01_1, a3rs.a3rs01_2 }, true, 3 * k + 1);
                                //ElementStiffnesses.ProccessVariable(9, new double[] { a3rs.a3rs02_0, a3rs.a3rs02_1, a3rs.a3rs02_2 }, true, 3 * k + 2);
                            }

                            ElementStiffnesses.saveOriginalState = false;

                        }
                    }


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

        /// <summary>
        /// 
        /// </summary>
        /// <param name="surfaceBasisVectorDerivative1"></param>
        /// <param name="surfaceBasisVectorDerivative2"></param>
        /// <param name="surfaceBasisVectorDerivative12"></param>
        /// <param name="ddksi_r">Second derivative ksi[i,j]</param>
        /// <param name="ddksi_s">Second derivative ksi[k,j]</param>
        
        /// <param name="ddheta_r">Second derivative heta[i,j]</param>
        /// <param name="ddheta_s">Second derivative heta[i,j]</param>
        /// <param name="dksidheta_r">Second derivative ksi and heta[i,j]</param>
        /// <param name="dksidheta_s">Second derivative ksi and heta[k,j]</param>
        /// <param name="a3rsAlternative"></param>
        /// <param name="a3r"></param>
        /// <param name="a3s"></param>
        /// <param name="da3_drds"></param>
        /// <returns></returns>
        internal static Bab_rs CalculateBab_rs(double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12,
            double ddksi_r, double ddksi_s, double ddheta_r, double ddheta_s,  double dksidheta_r,
            double dksidheta_s,a3rs a3rsAlternative, a3r a3r, a3r a3s, double[,][] da3_drds)
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
            var da3 = new double[3] {da3_drds[0, 0][0], da3_drds[0, 0][1], da3_drds[0, 0][2]};
            //[A]: 0--> b11rs , a11r ,a11s surfaceBasisVectorDerivative1,   
            //     1--> b22rs , a22r ,a22s surfaceBasisVectorDerivative2, 
            //     2--> b12rs , a12r ,a12s surfaceBasisVectorDerivative12, times x2
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s..........................................
            Bab_rsAlternative.Bab_rs00_0 = ddksi_r*a3sVecs_s0_0 + ddksi_s*a3rVecs_r0_0 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs00_1 = ddheta_r*a3sVecs_s0_0 + ddheta_s*a3rVecs_r0_0 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2]; 
            Bab_rsAlternative.Bab_rs00_2 = dksidheta_r*a3sVecs_s0_0 + dksidheta_s*a3rVecs_r0_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2]
                                         + dksidheta_r*a3sVecs_s0_0 + dksidheta_s*a3rVecs_r0_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //...............[r, s]
            da3[0] = da3_drds[0, 1][0]; da3[1] = da3_drds[0, 1][1]; da3[2] = da3_drds[0, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s.......................................................................[r,s]
            Bab_rsAlternative.Bab_rs01_0 = ddksi_r*a3sVecs_s1_0 + ddksi_s*a3rVecs_r0_1 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs01_1 = ddheta_r*a3sVecs_s1_0 + ddheta_s*a3rVecs_r0_1 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2]; 
            Bab_rsAlternative.Bab_rs01_2 = dksidheta_r*a3sVecs_s1_0 + dksidheta_s*a3rVecs_r0_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2]
                                         + dksidheta_r*a3sVecs_s1_0 + dksidheta_s*a3rVecs_r0_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //...............[r, s]
            da3[0] = da3_drds[0, 2][0]; da3[1] = da3_drds[0, 2][1]; da3[2] = da3_drds[0, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s]..........r_s............................................................................[r,s]
            Bab_rsAlternative.Bab_rs02_0 = ddksi_r*a3sVecs_s2_0 + ddksi_s*a3rVecs_r0_2 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs02_1 = ddheta_r*a3sVecs_s2_0 + ddheta_s*a3rVecs_r0_2 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs02_2 = dksidheta_r*a3sVecs_s2_0 + dksidheta_s*a3rVecs_r0_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s2_0 + dksidheta_s*a3rVecs_r0_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 0][0]; da3[1] = da3_drds[1, 0][1]; da3[2] = da3_drds[1, 0][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s..........................................................................[r,s]
            Bab_rsAlternative.Bab_rs10_0 = ddksi_r*a3sVecs_s0_1 + ddksi_s*a3rVecs_r1_0 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs10_1 = ddheta_r*a3sVecs_s0_1 + ddheta_s*a3rVecs_r1_0 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2]; 
            Bab_rsAlternative.Bab_rs10_2 = dksidheta_r*a3sVecs_s0_1 + dksidheta_s*a3rVecs_r1_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s0_1 + dksidheta_s*a3rVecs_r1_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 1][0]; da3[1] = da3_drds[1, 1][1]; da3[2] = da3_drds[1, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]
            Bab_rsAlternative.Bab_rs11_0 = ddksi_r*a3sVecs_s1_1 + ddksi_s*a3rVecs_r1_1 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2]; 
            Bab_rsAlternative.Bab_rs11_1 = ddheta_r*a3sVecs_s1_1 + ddheta_s*a3rVecs_r1_1 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs11_2 = dksidheta_r*a3sVecs_s1_1 + dksidheta_s*a3rVecs_r1_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s1_1 + dksidheta_s*a3rVecs_r1_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[1, 2][0]; da3[1] = da3_drds[1, 2][1]; da3[2] = da3_drds[1, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s].....
            Bab_rsAlternative.Bab_rs12_0 = ddksi_r*a3sVecs_s2_1 + ddksi_s*a3rVecs_r1_2 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs12_1 = ddheta_r*a3sVecs_s2_1 + ddheta_s*a3rVecs_r1_2 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs12_2 = dksidheta_r*a3sVecs_s2_1 + dksidheta_s*a3rVecs_r1_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s2_1 + dksidheta_s*a3rVecs_r1_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 0][0]; da3[1] = da3_drds[2, 0][1]; da3[2] = da3_drds[2, 0][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]..
            Bab_rsAlternative.Bab_rs20_0 = ddksi_r*a3sVecs_s0_2 + ddksi_s*a3rVecs_r2_0 +s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs20_1 = ddheta_r*a3sVecs_s0_2 + ddheta_s*a3rVecs_r2_0 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs20_2 = dksidheta_r*a3sVecs_s0_2 + dksidheta_s*a3rVecs_r2_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s0_2 + dksidheta_s*a3rVecs_r2_0 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 1][0]; da3[1] = da3_drds[2, 1][1]; da3[2] = da3_drds[2, 1][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s].
            Bab_rsAlternative.Bab_rs21_0 = ddksi_r*a3sVecs_s1_2 + ddksi_s*a3rVecs_r2_1 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs21_1 = ddheta_r*a3sVecs_s1_2 + ddheta_s*a3rVecs_r2_1 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs21_2 = dksidheta_r*a3sVecs_s1_2 + dksidheta_s*a3rVecs_r2_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s1_2 + dksidheta_s*a3rVecs_r2_1 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];
            //............[r, s]
            da3[0] = da3_drds[2, 2][0]; da3[1] = da3_drds[2, 2][1]; da3[2] = da3_drds[2, 2][2];
            //......................rs[A]..[A}..[r].........s_r.+..[A].[s].........r_s...........................................................................................[r,s]..
            Bab_rsAlternative.Bab_rs22_0 = ddksi_r*a3sVecs_s2_2 + ddksi_s*a3rVecs_r2_2 + s11[0]*da3[0]+s11[1]*da3[1]+s11[2]*da3[2];
            Bab_rsAlternative.Bab_rs22_1 = ddheta_r*a3sVecs_s2_2 + ddheta_s*a3rVecs_r2_2 + s22[0]*da3[0]+s22[1]*da3[1]+s22[2]*da3[2];
            Bab_rsAlternative.Bab_rs22_2 = dksidheta_r*a3sVecs_s2_2 + dksidheta_s*a3rVecs_r2_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2] 
                                         + dksidheta_r*a3sVecs_s2_2 + dksidheta_s*a3rVecs_r2_2 + s12[0]*da3[0]+s12[1]*da3[1]+s12[2]*da3[2];

            return Bab_rsAlternative;
        }


        private Bab_rs CalculateBab_rs_OLD(double[] surfaceBasisVectorDerivative1,
            double[] surfaceBasisVectorDerivative2, double[] surfaceBasisVectorDerivative12,
            Nurbs2D nurbs, int i, int k, int j, a3rs a3rsAlternative, a3r a3r, a3r a3s, Vector[,] da3_drds)
        {
            var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
            var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
            var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);

            var a11s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[k, j]);
            var a22s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[k, j]);
            var a12s = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[k, j]);

            Vector[] a3rVecs = new Vector[3];
            a3rVecs[0] = Vector.CreateFromArray(new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 });
            a3rVecs[1] = Vector.CreateFromArray(new double[] { a3r.a3r01, a3r.a3r11, a3r.a3r21 });
            a3rVecs[2] = Vector.CreateFromArray(new double[] { a3r.a3r02, a3r.a3r12, a3r.a3r22 });

            Vector[] a3sVecs = new Vector[3];
            a3sVecs[0] = Vector.CreateFromArray(new double[] { a3s.a3r00, a3s.a3r10, a3s.a3r20 });
            a3sVecs[1] = Vector.CreateFromArray(new double[] { a3s.a3r01, a3s.a3r11, a3s.a3r21 });
            a3sVecs[2] = Vector.CreateFromArray(new double[] { a3s.a3r02, a3s.a3r12, a3s.a3r22 });

            var Bab_rsAlternative = new Bab_rs();

            double[,] b11_rs = new double[3, 3];
            double[,] b22_rs = new double[3, 3];
            double[,] b12_rs = new double[3, 3];
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    b11_rs[r1, s1] = a11r.GetColumn(r1).DotProduct(a3sVecs[s1]) + a11s.GetColumn(s1).DotProduct(a3rVecs[r1]) + Vector.CreateFromArray(surfaceBasisVectorDerivative1).DotProduct(da3_drds[r1, s1]);
                    b22_rs[r1, s1] = a22r.GetColumn(r1).DotProduct(a3sVecs[s1]) + a22s.GetColumn(s1).DotProduct(a3rVecs[r1]) + Vector.CreateFromArray(surfaceBasisVectorDerivative2).DotProduct(da3_drds[r1, s1]);
                    b12_rs[r1, s1] = a12r.GetColumn(r1).DotProduct(a3sVecs[s1]) + a12s.GetColumn(s1).DotProduct(a3rVecs[r1]) + Vector.CreateFromArray(surfaceBasisVectorDerivative12).DotProduct(da3_drds[r1, s1])
                        + a12r.GetColumn(r1).DotProduct(a3sVecs[s1]) + a12s.GetColumn(s1).DotProduct(a3rVecs[r1]) + Vector.CreateFromArray(surfaceBasisVectorDerivative12).DotProduct(da3_drds[r1, s1]);
                }
            }

            Bab_rsAlternative.Bab_rs00_0 = b11_rs[0, 0]; Bab_rsAlternative.Bab_rs00_1 = b22_rs[0, 0]; Bab_rsAlternative.Bab_rs00_2 = b12_rs[0, 0];
            Bab_rsAlternative.Bab_rs01_0 = b11_rs[0, 1]; Bab_rsAlternative.Bab_rs01_1 = b22_rs[0, 1]; Bab_rsAlternative.Bab_rs01_2 = b12_rs[0, 1];
            Bab_rsAlternative.Bab_rs02_0 = b11_rs[0, 2]; Bab_rsAlternative.Bab_rs02_1 = b22_rs[0, 2]; Bab_rsAlternative.Bab_rs02_2 = b12_rs[0, 2];

            Bab_rsAlternative.Bab_rs10_0 = b11_rs[1, 0]; Bab_rsAlternative.Bab_rs10_1 = b22_rs[1, 0]; Bab_rsAlternative.Bab_rs10_2 = b12_rs[1, 0];
            Bab_rsAlternative.Bab_rs11_0 = b11_rs[1, 1]; Bab_rsAlternative.Bab_rs11_1 = b22_rs[1, 1]; Bab_rsAlternative.Bab_rs11_2 = b12_rs[1, 1];
            Bab_rsAlternative.Bab_rs12_0 = b11_rs[1, 2]; Bab_rsAlternative.Bab_rs12_1 = b22_rs[1, 2]; Bab_rsAlternative.Bab_rs12_2 = b12_rs[1, 2];

            Bab_rsAlternative.Bab_rs20_0 = b11_rs[2, 0]; Bab_rsAlternative.Bab_rs20_1 = b22_rs[2, 0]; Bab_rsAlternative.Bab_rs20_2 = b12_rs[2, 0];
            Bab_rsAlternative.Bab_rs21_0 = b11_rs[2, 1]; Bab_rsAlternative.Bab_rs21_1 = b22_rs[2, 1]; Bab_rsAlternative.Bab_rs21_2 = b12_rs[2, 1];
            Bab_rsAlternative.Bab_rs22_0 = b11_rs[2, 2]; Bab_rsAlternative.Bab_rs22_1 = b22_rs[2, 2]; Bab_rsAlternative.Bab_rs22_2 = b12_rs[2, 2];

            return Bab_rsAlternative;
        }

        private static void Calculate_a3rs_OLD_Dimitris(double[] surfaceBasisVector1, double[] surfaceBasisVector2,
            double[] surfaceBasisVector3, double J1, double dksi_r, double dheta_r, double dksi_s, double dheta_s, ref a3rs a3rsOut)
        {
            #region Initializations
            var s10 = surfaceBasisVector1[0];
            var s11 = surfaceBasisVector1[1];
            var s12 = surfaceBasisVector1[2];

            var s20 = surfaceBasisVector2[0];
            var s21 = surfaceBasisVector2[1];
            var s22 = surfaceBasisVector2[2];

            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var aux1Term1 = (dheta_s * dksi_r - dheta_r * dksi_s) * J1;
            var aux2Term1 = dheta_r * dksi_s - dheta_s * dksi_r;
            var aux3Term1 = aux2Term1 * J1;

            var aux1Term2 = dheta_s * s11 - dksi_s * s21;
            var aux2Term2 = dheta_s * s12 - dksi_s * s22;
            var aux3Term2 = dheta_r * s12 - dksi_r * s22;
            var aux4Term2 = dheta_r * s11 - dksi_r * s21;
            var aux5Term2 = dheta_s * s10 - dksi_s * s20;
            var aux6Term2 = dheta_r * s10 - dksi_r * s20;
            var aux7Term2 = s32 * aux1Term2 - s31 * aux2Term2;
            var aux8Term2 = s32 * aux5Term2 - s30 * aux2Term2;
            var aux9Term2 = s31 * aux5Term2 - s30 * aux1Term2;
            var J1squared = J1 * J1;

            var aux1Term3 = s32 * aux4Term2 - s31 * aux3Term2;
            var aux2Term3 = s32 * aux6Term2 - s30 * aux3Term2;
            var aux3Term3 = s31 * aux6Term2 - s30 * aux4Term2;

            var aux1 = aux4Term2 * aux1Term2;
            var aux2 = aux3Term2 * aux2Term2;
            var aux3 = aux1Term3 * aux7Term2;
            var aux4 = aux7Term2 * aux3Term2;
            var aux5 = aux1Term3 * aux2Term2;
            var aux6 = aux7Term2 * aux4Term2;
            var aux7 = aux1Term3 * aux1Term2;
            var aux8 = aux1 + aux2 - aux3;
            var aux9 = aux1Term3 * aux8Term2;
            var aux10 = aux4Term2 * aux5Term2;
            var aux11 = J1 * s32 * J1 * aux2Term1;
            var aux12 = aux8Term2 * aux4Term2;
            var aux13 = aux8Term2 * aux3Term2;
            var aux14 = aux1Term3 * aux5Term2;
            var aux15 = (aux10 + aux11 - aux9);
            var aux16 = aux9Term2 * aux1Term3;
            var aux17 = aux3Term2 * aux5Term2;
            var aux18 = J1 * s31 * J1 * aux2Term1;
            var aux19 = aux9Term2 * aux3Term2;
            var aux20 = aux9Term2 * aux4Term2;
            var aux21 = (aux17 - aux18 + aux16);
            var aux22 = aux2Term3 * aux7Term2;
            var aux23 = aux6Term2 * aux1Term2;
            var aux24 = aux2Term3 * aux2Term2;
            var aux25 = aux7Term2 * aux6Term2;
            var aux26 = aux2Term3 * aux1Term2;
            var aux27 = aux23 - aux11 - aux22;
            var aux28 = aux2Term3 * aux8Term2;
            var aux29 = aux6Term2 * aux5Term2;
            var aux30 = aux8Term2 * aux6Term2;
            var aux31 = aux2Term3 * aux5Term2;
            var aux32 = (aux29 + aux2 - aux28);
            var aux33 = aux2Term3 * aux9Term2;
            var aux34 = J1 * s30 * J1 * aux2Term1;
            var aux35 = aux3Term2 * aux1Term2;
            var aux36 = aux9Term2 * aux6Term2;
            var aux37 = (aux35 + aux34 - aux33);
            var aux38 = aux3Term3 * aux7Term2;
            var aux39 = aux6Term2 * aux2Term2;
            var aux40 = aux3Term3 * aux2Term2;
            var aux41 = aux3Term3 * aux1Term2;
            var aux42 = (aux39 + aux18 + aux38);
            var aux43 = aux3Term3 * aux8Term2;
            var aux44 = aux4Term2 * aux2Term2;
            var aux45 = aux3Term3 * aux5Term2;
            var aux46 = (aux44 - aux34 - aux43);
            var aux47 = aux3Term3 * aux9Term2;
            var aux48 = (aux29 + aux1 - aux47);
            #endregion

            #region Term1



            a3rsOut.a3rs00_0 = (2 * s30 * aux3 - s30 * aux8) / J1squared;
            a3rsOut.a3rs00_1 = (aux4 + aux5 + 2 * s31 * aux3 - s31 * aux8) / J1squared;
            a3rsOut.a3rs00_2 = -(aux6 + aux7 - 2 * s32 * aux3 + s32 * aux8) / J1squared;
            
            a3rsOut.a3rs01_0 = -(aux5 + 2 * s30 * aux9- s30 * aux15) / J1squared;
            a3rsOut.a3rs01_1 = -(aux13 + 2 * s31 * aux9 - s31 * aux15) / J1squared;
            a3rsOut.a3rs01_2 = aux1Term1 + (aux12 + aux14 - 2 * s32 * aux9 + s32 * aux15) / J1squared;
            
            a3rsOut.a3rs02_0 = (aux7 + 2 * s30 * aux16 + s30 * aux21) / J1squared;
            a3rsOut.a3rs02_1 = aux3Term1 + (aux19 - aux14 + 2 * s31 * aux16 + s31 * aux21) / J1squared;
            a3rsOut.a3rs02_2 = -(aux20 - 2 * s32 * aux16 - s32 * aux21) / J1squared;
            
            a3rsOut.a3rs10_0 = -(aux4 + 2 * s30 * aux22 - s30 * aux27) / J1squared;
            a3rsOut.a3rs10_1 = -(aux24 + 2 * s31 * aux22 - s31 * aux27) / J1squared;
            a3rsOut.a3rs10_2 = aux3Term1 + (aux25 + aux26 - 2 * s32 * aux22 + s32 * aux27) / J1squared;
            
            a3rsOut.a3rs11_0 = (aux13 + aux24 + 2 * s30 * aux28 - s30 * aux32) / J1squared;
            a3rsOut.a3rs11_1 = (2 * s31 * aux28 - s31 * aux32) / J1squared;
            a3rsOut.a3rs11_2 = -(aux30 + aux31 - 2 * s32 * aux28 + s32 * aux32) /J1squared;
            
            a3rsOut.a3rs12_0 = aux1Term1 - (aux19 + aux26 + 2 * s30 * aux33 - s30 * aux37) / J1squared;
            a3rsOut.a3rs12_1 = (aux31 - 2 * s31 * aux33 + s31 * aux37) / J1squared;
            a3rsOut.a3rs12_2 = (aux36 - 2 * s32 * aux33 + s32 * aux37) / J1squared;
            
            a3rsOut.a3rs20_0 = (aux6 + 2 * s30 * aux38 + s30 * aux42) / J1squared;
            a3rsOut.a3rs20_1 = aux1Term1 - (aux25 - aux40 - 2 * s31 * aux38 - s31 * aux42) / J1squared;
            a3rsOut.a3rs20_2 = -(aux41 - 2 * s32 * aux38 - s32 * aux42) / J1squared;
            
            a3rsOut.a3rs21_0 = aux3Term1 - (aux12 + aux40 + 2 * s30 * aux43 - s30 * aux46) / J1squared;
            a3rsOut.a3rs21_1 = (aux30 - 2 * s31 * aux43 + s31 * aux46) / J1squared;
            a3rsOut.a3rs21_2 = (aux45 - 2 * s32 * aux43 + s32 * aux46) / J1squared;
            
            a3rsOut.a3rs22_0 = (aux20 + aux41 + 2 * s30 * aux47 - s30 * aux48) / J1squared;
            a3rsOut.a3rs22_1 = -(aux36 + aux45 - 2 * s31 * aux47 + s31 * aux48) / J1squared;
            a3rsOut.a3rs22_2 = (2 * s32 * aux47 - s32 * aux48) / J1squared;

            #endregion
        }

        internal static (a3rs ,double[,][], double[][], double[][], double[], double[], double[,], double[], double[,][] ) Calculate_a3rs(Vector surfaceBasisVector1, Vector surfaceBasisVector2,
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

                    double scale2 = -((double) 1 / (Math.Pow(J1, 2))); //denominator of vectors 2 3 and 4 and a minus.

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

                    double scale5 = ((double) 1) / Math.Pow(J1, 3);

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

        private static double[] CalculateTerm525(Vector surfaceBasisVector3, double J1, double[] dnorma3_dr,
            double[][] da3tilde_dr)
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

            return a3_tilde;
        }

        private static void Calculate_da3tilde_drds(double dKsi_r, double dKsi_s, double dHeta_r, double dHeta_s,
            double[,][] da3tilde_drds)
        {
            //da3tilde_drds[r1, s1] = a1r.GetColumn(r1).CrossProduct(a2s.GetColumn(s1)) +
            //                        a1s.GetColumn(s1).CrossProduct(a2r.GetColumn(r1));
            
            var dksiRxdHetaS = dKsi_r*dHeta_s;
            var dHetaRxdKsiS = dHeta_r*dKsi_s;
            da3tilde_drds[0, 0]= new double[3];
            da3tilde_drds[0, 1]= new double[]{0,0,dksiRxdHetaS-dHetaRxdKsiS};
            da3tilde_drds[0, 2]= new double[]{0,dHetaRxdKsiS-dksiRxdHetaS,0};

            da3tilde_drds[1, 0]= new double[]{0,0,dHetaRxdKsiS-dksiRxdHetaS};
            da3tilde_drds[1, 1]= new double[3];
            da3tilde_drds[1, 2]= new double[]{dksiRxdHetaS-dHetaRxdKsiS,0,0};

            da3tilde_drds[2, 0]= new double[]{0,dksiRxdHetaS-dHetaRxdKsiS,0};
            da3tilde_drds[2, 1]= new double[]{dHetaRxdKsiS-dksiRxdHetaS,0,0};
            da3tilde_drds[2, 2]= new double[3];
        }

        private static void Calculate_da3tilde_dr(Vector surfaceBasisVector1, Vector surfaceBasisVector2, double dksi_r,
            double dHeta_r, double[][] da3tilde_dr)
        {
            //da3tilde_dr[r1] = a1r.GetColumn(r1).CrossProduct(surfaceBasisVector2) + surfaceBasisVector1.CrossProduct(a2r.GetColumn(r1));

            da3tilde_dr[0]=new double[]
            {
                0,
                -dksi_r*surfaceBasisVector2[2]+surfaceBasisVector1[2]*dHeta_r,
                dksi_r*surfaceBasisVector2[1]-surfaceBasisVector1[1]*dHeta_r
            };

            da3tilde_dr[1]=new double[]
            {
               dksi_r*surfaceBasisVector2[2]-surfaceBasisVector1[2]*dHeta_r,
               0,
               -dksi_r*surfaceBasisVector2[0]+surfaceBasisVector1[0]*dHeta_r
            };

            da3tilde_dr[2]=new double[]
            {
               -dksi_r*surfaceBasisVector2[1]+dHeta_r*surfaceBasisVector1[1],
               dksi_r*surfaceBasisVector2[0]-dHeta_r*surfaceBasisVector1[0],
               0
            };
        }

        internal static (a3rs ,Vector[,], Vector[], Vector[], double[], double[], double[,], Vector, Vector[,] ) Calculate_a3rs_OLD(Vector surfaceBasisVector1, Vector surfaceBasisVector2,
            Vector surfaceBasisVector3, double J1, Matrix3by3 a1r, Matrix3by3 a2s, Matrix3by3 a1s, Matrix3by3 a2r)
        {
            var da3_drds = new Vector[3, 3];
            Vector[,] da3tilde_drds = new Vector[3, 3];
            Vector[] da3tilde_dr = new Vector[3];
            Vector[] da3tilde_ds = new Vector[3];
            Vector a3_tilde;
            double[] dnorma3_dr = new double[3];
            double[] dnorma3_ds = new double[3];
            double[,] dnorma3_drds = new double[3, 3];


            //5.30
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    da3tilde_drds[r1, s1] = a1r.GetColumn(r1).CrossProduct(a2s.GetColumn(s1)) + a1s.GetColumn(s1).CrossProduct(a2r.GetColumn(r1));
                }
            }

            //5.24
            for (int r1 = 0; r1 < 3; r1++)
            {
                da3tilde_dr[r1] = a1r.GetColumn(r1).CrossProduct(surfaceBasisVector2) + surfaceBasisVector1.CrossProduct(a2r.GetColumn(r1));
            }

            //5.24
            for (int s1 = 0; s1 < 3; s1++)
            {
                da3tilde_ds[s1] = a1s.GetColumn(s1).CrossProduct(surfaceBasisVector2) + surfaceBasisVector1.CrossProduct(a2s.GetColumn(s1));
            }

            //5.25
            a3_tilde = surfaceBasisVector3.Scale(J1);
            for (int r1 = 0; r1 < 3; r1++)
            {
                dnorma3_dr[r1] = (a3_tilde.DotProduct(da3tilde_dr[r1])) / J1;
            }
            for (int s1 = 0; s1 < 3; s1++)
            {
                dnorma3_ds[s1] = (a3_tilde.DotProduct(da3tilde_ds[s1])) / J1;
            }

            //5.31
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    double firstNumerator = da3tilde_drds[r1, s1].DotProduct(a3_tilde) + da3tilde_dr[r1].DotProduct(da3tilde_ds[s1]);
                    double firstDenominator = J1;
                    double secondNumerator = (da3tilde_dr[r1].DotProduct(a3_tilde)) * (da3tilde_ds[s1].DotProduct(a3_tilde));
                    double secondDenominator = Math.Pow(J1, 3);

                    dnorma3_drds[r1, s1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);
                }
            }

            //5.32
            for (int r1 = 0; r1 < 3; r1++)
            {
                for (int s1 = 0; s1 < 3; s1++)
                {
                    Vector firstVec = da3tilde_drds[r1, s1].Scale(((double)1 / J1));

                    double scale2 =-( (double)1 / (Math.Pow(J1, 2))); //denominator of vectors 2 3 and 4 and a minus.

                    Vector secondVec = da3tilde_dr[r1].Scale(dnorma3_ds[s1]).Scale(scale2);

                    Vector thirdVec = da3tilde_ds[s1].Scale(dnorma3_dr[r1]).Scale(scale2);

                    Vector fourthVec = a3_tilde.Scale(dnorma3_drds[r1, s1]).Scale(scale2);

                    double scale5 = ((double)1) / Math.Pow(J1, 3);

                    Vector fifthvector = a3_tilde.Scale(2 * dnorma3_dr[r1] * dnorma3_ds[s1]).Scale(scale5);

                    da3_drds[r1, s1] = firstVec + secondVec + thirdVec + fourthVec + fifthvector;

                }
            }

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
        
        private void CalculateA3r_OLD(double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double dksi_r, double dheta_r, double J1, ref a3r da3_unit_dr_out)
        {
            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];

            var a1r = Matrix3by3.CreateIdentity().Scale(dksi_r);
            var a2r = Matrix3by3.CreateIdentity().Scale(dheta_r);
            Vector[] da3tilde_dr = new Vector[3];
            //5.24
            for (int r1 = 0; r1 < 3; r1++)
            {
                da3tilde_dr[r1] = a1r.GetColumn(r1).CrossProduct(Vector.CreateFromArray(surfaceBasisVector2)) +
                    Vector.CreateFromArray(surfaceBasisVector1).CrossProduct(a2r.GetColumn(r1));
            }

            var da3_tilde_dr00 = da3tilde_dr[0][0];
            var da3_tilde_dr10 = da3tilde_dr[0][1];
            var da3_tilde_dr20 = da3tilde_dr[0][2];

            var da3_tilde_dr01 = da3tilde_dr[1][0];
            var da3_tilde_dr11 = da3tilde_dr[1][1];
            var da3_tilde_dr21 = da3tilde_dr[1][2];

            var da3_tilde_dr02 = da3tilde_dr[2][0];
            var da3_tilde_dr12 = da3tilde_dr[2][1];
            var da3_tilde_dr22 = da3tilde_dr[2][2];

            //var da3_tilde_dr10 = dheta_r * surfaceBasisVector1[2] - dksi_r * surfaceBasisVector2[2];
            //var da3_tilde_dr20 = dksi_r * surfaceBasisVector2[1] - dheta_r * surfaceBasisVector1[1];

            //var da3_tilde_dr01 = dksi_r * surfaceBasisVector2[2] - dheta_r * surfaceBasisVector1[2];
            //var da3_tilde_dr21 = dheta_r * surfaceBasisVector1[0] - dksi_r * surfaceBasisVector2[0];

            //var da3_tilde_dr02 = dheta_r * surfaceBasisVector1[1] - dksi_r * surfaceBasisVector2[1];
            //var da3_tilde_dr12 = dksi_r * surfaceBasisVector2[0] - dheta_r * surfaceBasisVector1[0];


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

        internal static void CalculateA3r(double[] surfaceBasisVector1,
            double[] surfaceBasisVector2, double[] surfaceBasisVector3,
            double dksi_r, double dheta_r, double J1, ref a3r da3_unit_dr_out)
        {
            var s30 = surfaceBasisVector3[0];
            var s31 = surfaceBasisVector3[1];
            var s32 = surfaceBasisVector3[2];
            
            var da3tilde_dr = new double[3][];
            //5.24
            var a1r= new double[3];
            var a2r= new double[3];
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
            sum1[0] = 0.0;sum1[1] = 0.0;sum1[2] = 0.0;
            sum2[0] = 0.0;sum2[1] = 0.0;sum2[2] = 0.0;
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
            sum1[0] = 0.0;sum1[1] = 0.0;sum1[2] = 0.0;
            sum2[0] = 0.0;sum2[1] = 0.0;sum2[2] = 0.0;
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

        internal void CalculateKmembraneNL(ControlPoint[] controlPoints, ref Forces membraneForces, Nurbs2D nurbs,
            int j, double[,] KmembraneNLOut)
        {
            for (var i = 0; i < controlPoints.Length; i++)
            {
                var dksi_r = nurbs.NurbsDerivativeValuesKsi[i, j];
                var dheta_r = nurbs.NurbsDerivativeValuesHeta[i, j];

                for (int k = 0; k < controlPoints.Length; k++)
                {
                    var dksi_s = nurbs.NurbsDerivativeValuesKsi[k, j];
                    var dheta_s = nurbs.NurbsDerivativeValuesHeta[k, j];

                    
                    var aux = membraneForces.v0 * dksi_r * dksi_s +
                              membraneForces.v1 * dheta_r * dheta_s +
                              membraneForces.v2 * (dksi_r * dheta_s + dksi_s * dheta_r);

                    KmembraneNLOut[i * 3, k * 3] += aux;
                    KmembraneNLOut[i * 3 + 1, k * 3 + 1] += aux;
                    KmembraneNLOut[i * 3 + 2, k * 3 + 2] += aux;
                }
            }
        }
        
        internal void CalculateMembraneDeformationMatrix(int controlPointsCount, Nurbs2D nurbs, int j,
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
                var dKsi = nurbs.NurbsDerivativeValuesKsi[column / 3, j];
                var dHeta = nurbs.NurbsDerivativeValuesHeta[column / 3, j];
                
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
        
        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(NurbsKirchhoffLoveShellElementNLDevelop shellElement)
        {
            var gauss = new GaussQuadrature();
            var medianSurfaceGP = gauss.CalculateElementGaussPoints(shellElement.Patch.DegreeKsi,
                shellElement.Patch.DegreeHeta, shellElement.Knots.ToList());
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

        private const int ThicknessIntegrationDegree = 2;

        public double Thickness { get; set; }
    }

    
}