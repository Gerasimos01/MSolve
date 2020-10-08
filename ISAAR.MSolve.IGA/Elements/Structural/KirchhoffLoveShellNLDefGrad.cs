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
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.IGA.Elements.Structural
{
    
    public class KirchhoffLoveShellNLDefGrad : Element, IStructuralIsogeometricElement
    {
        internal IShapeFunction2D _shapeFunctions;

        public double[] _solution;

        internal Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IContinuumMaterial3DDefGrad>>
            materialsAtThicknessGP =
                new Dictionary<GaussLegendrePoint3D, Dictionary<GaussLegendrePoint3D, IContinuumMaterial3DDefGrad>>();

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

        public KirchhoffLoveShellNLDefGrad(IContinuumMaterial3DDefGrad shellMaterial, IList<Knot> elementKnots,
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
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IContinuumMaterial3DDefGrad>());
                foreach (var point in thicknessIntegrationPoints[medianSurfaceGP])
                {
                    materialsAtThicknessGP[medianSurfaceGP].Add(point, (IContinuumMaterial3DDefGrad)shellMaterial.Clone());
                }
            }

            _controlPoints = elementControlPoints.ToArray();
        }

        public KirchhoffLoveShellNLDefGrad(List<IContinuumMaterial3DDefGrad> shellMaterials, IList<Knot> elementKnots,
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
                materialsAtThicknessGP.Add(medianSurfaceGP, new Dictionary<GaussLegendrePoint3D, IContinuumMaterial3DDefGrad>());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][0], (IContinuumMaterial3DDefGrad)shellMaterials[0].Clone());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][1], (IContinuumMaterial3DDefGrad)shellMaterials[1].Clone());
                materialsAtThicknessGP[medianSurfaceGP].Add(thicknessIntegrationPoints[medianSurfaceGP][2], (IContinuumMaterial3DDefGrad)shellMaterials[2].Clone());
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

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)//.
        {
            //var shellElement = (KirchhoffLoveShellNLDefGrad)element;
            var elementNodalForces = new double[element.ControlPointsDictionary.Count * 3];
            var elementNodalMembraneForces = new double[element.ControlPointsDictionary.Count * 3];
            var elementNodalBendingForces = new double[element.ControlPointsDictionary.Count * 3];

            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(_controlPoints);
            
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();


            var numberOfControlPoints = _controlPoints.Length;



            var forcesDevelop_v3 = new double[element.ControlPointsDictionary.Count * 3];
            double[,] dF3Dtensor_dr_Tr_check = new double[3, 3];
            for (int j = 0; j < gaussPoints.Length; j++)
            {

                CalculateJacobian(newControlPoints, _shapeFunctions , j, jacobianMatrix);

                var hessianMatrix = CalculateHessian(newControlPoints, _shapeFunctions , j);

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

                var wfactor = InitialJ1[j] * gaussPoints[j].WeightFactor;




                #region develop formulation

                var thicknessGPoints = thicknessIntegrationPoints[gaussPoints[j]];
                var materialpoint = materialsAtThicknessGP[gaussPoints[j]][thicknessGPoints[0]];


                var elementControlPoints = CurrentControlPoint(_controlPoints);

                var a11 = surfaceBasisVectorDerivative1;
                var a22 = surfaceBasisVectorDerivative2;
                var a12 = surfaceBasisVectorDerivative12;
                var a1 = surfaceBasisVector1;
                var a2 = surfaceBasisVector2;
                var a3 = surfaceBasisVector3; // einai to mono pou einai normalised
                double[] a3_tilde = new double[] { a3[0] * (J1), a3[1] * (J1), a3[2] * (J1) };

                (double[] da3tilde_dksi, double[] da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, double[] da3_dksi, double[] da3_dheta) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1, a2, a11, a22, a12, a3, J1);

                if (j == ElementStiffnesses.gpNumberToCheck && ElementStiffnesses.performCalculations)
                {
                    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                    //ElementStiffnesses.ProccessVariable(11, new double[] { J1 }, false);
                    //ElementStiffnesses.ProccessVariable(12, da3tilde_dksi.CopyToArray(), false);
                    //ElementStiffnesses.ProccessVariable(13, da3tilde_dheta.CopyToArray(), false);
                    //ElementStiffnesses.ProccessVariable(14, new double[] { da3norm_dksi }, false);
                    //ElementStiffnesses.ProccessVariable(15, new double[] { da3norm_dheta }, false);
                    ElementStiffnesses.ProccessVariable(16, da3_dksi, false);
                    ElementStiffnesses.ProccessVariable(17, da3_dheta, false);
                    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                }

                #region original Config
                var originalControlPoints = element.ControlPoints.ToArray();
                var originalHessianMatrix = CalculateHessian(originalControlPoints, _shapeFunctions , j);
                double[,] jacobian_init = new double[3, 3];
                CalculateJacobian(originalControlPoints, _shapeFunctions , j, jacobian_init);
                var a1_init = CalculateSurfaceBasisVector(jacobian_init, 0);
                var a2_init = CalculateSurfaceBasisVector(jacobian_init, 1);
                var a3_init = new double[3]
                {
                    a1_init[1] * a2_init[2] - a1_init[2] * a2_init[1],
                    a1_init[2] * a2_init[0] - a1_init[0] * a2_init[2],
                    a1_init[0] * a2_init[1] - a1_init[1] * a2_init[0]
                };

                var J1_init = Math.Sqrt(a3_init[0] * a3_init[0] +
                                        a3_init[1] * a3_init[1] +
                                        a3_init[2] * a3_init[2]);

                a3_init[0] /= J1_init;
                a3_init[1] /= J1_init;
                a3_init[2] /= J1_init;

                var a11_init = CalculateSurfaceBasisVector(originalHessianMatrix, 0);
                var a22_init = CalculateSurfaceBasisVector(originalHessianMatrix, 1);
                var a12_init = CalculateSurfaceBasisVector(originalHessianMatrix, 2);
                #endregion

                (double[] da3tilde_dksi_init, double[] da3tilde_dheta_init, double da3norm_dksi_init, double da3norm_dheta_init, double[] da3_dksi_init, double[] da3_dheta_init) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1_init, a2_init, a11_init, a22_init, a12_init, a3_init, J1_init);
                var thicknessPoints = thicknessIntegrationPoints[gaussPoints[j]];
                double[,] forceIntegration = new double[thicknessGPoints.Count(), 3]; //[thickness,dof]
                double[,] thicknesCoeffs = new double[thicknessGPoints.Count(), 3]; //[thickness,dof]
                for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                {
                    var materialDevelop = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    var thicknessPoint = thicknessPoints[i1];
                    var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    var w = thicknessPoint.WeightFactor;
                    var z = thicknessPoint.Zeta;


                    //Vector G1 = a1_init + da3_dksi_init.Scale(z);
                    //Vector G2 = a2_init + da3_dheta_init.Scale(z);
                    double[] G1 = new double[] { a1_init[0] + da3_dksi_init[0] * z, a1_init[1] + da3_dksi_init[1] * z, a1_init[2] + da3_dksi_init[2] * z };
                    double[] G2 = new double[] { a2_init[0] + da3_dheta_init[0] * z, a2_init[1] + da3_dheta_init[1] * z, a2_init[2] + da3_dheta_init[2] * z };

                    (double[] G_1, double[] G_2, double[] G_3) = CalculateContravariants(G1, G2, a3_init);

                    //Vector g1 = a1 + da3_dksi.Scale(z);
                    //Vector g2 = a2 + da3_dheta.Scale(z);
                    double[] g1 = new double[] { a1[0] + da3_dksi[0] * z, a1[1] + da3_dksi[1] * z, a1[2] + da3_dksi[2] * z };
                    double[] g2 = new double[] { a2[0] + da3_dheta[0] * z, a2[1] + da3_dheta[1] * z, a2[2] + da3_dheta[2] * z };




                    if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0) && ElementStiffnesses.performCalculations)
                    {
                        double[,] F_3D = new double[3, 3] { { g1[0]*G_1[0]+g2[0]*G_2[0], g1[0]*G_1[1]+g2[0]*G_2[1], g1[0]*G_1[2]+g2[0]*G_2[2] },
                                                            { g1[1]*G_1[0]+g2[1]*G_2[0], g1[1]*G_1[1]+g2[1]*G_2[1], g1[1]*G_1[2]+g2[1]*G_2[2] },
                                                            { g1[2]*G_1[0]+g2[2]*G_2[0], g1[2]*G_1[1]+g2[2]*G_2[1], g1[2]*G_1[2]+g2[2]*G_2[2] },
                        };
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                        double[] F_3D_vec = { F_3D[0, 0], F_3D[1, 1], F_3D[2, 2], F_3D[0, 1], F_3D[1, 2], F_3D[2, 0], F_3D[0, 2], F_3D[1, 0], F_3D[2, 1] };//  .
                        ElementStiffnesses.ProccessVariable(18, F_3D_vec, false);
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                    }



                    //(double[,] Aijkl_3D, double[] FPK_3D_vec, double[,] FPK_2D, double[,] Ei, _, double[,] ei, double[,] F_rve, double[,] GL3D, double[,] SPKMat3D,
                    //        double[,] ch01_GL_3D, double[,] ch01_SPKMat_3D, double[,] FPK_3D_tr_basis, double[] FPK_3D_a_vec, double[] FPK_3D_tr_basis_a_vec,
                    //        double[] ch02_FPK_3D_vec, tensorOrder2 ch03_FPK_3D, tensorOrder2 GLtensorProjected2, tensorOrder2 SPKtensorProjected2, tensorOrder2 FPKtensorProjected2
                    //        ) = transformations.CalculateTransformationsV2(g1, g2, a3, G1, G2, a3_init, G_1, G_2, G_3);

                    double[,] ei = new double[3, 3] { { g1[0], g2[0], a3[0] }, { g1[1], g2[1], a3[1] }, { g1[2], g2[2], a3[2] } };
                    double[,] Ei = new double[3, 3] { { G1[0], G2[0], a3_init[0] }, { G1[1], G2[1], a3_init[1] }, { G1[2], G2[2], a3_init[2] } };

                    NormaliseAndOrthogonaliseBasis(ei);
                    NormaliseAndOrthogonaliseBasis(Ei);



                    double[,] FPKrve = new double[2, 2] { { materialDevelop.Stresses[0], materialDevelop.Stresses[2] }, { materialDevelop.Stresses[3], materialDevelop.Stresses[1] } };
                    //double[,] FPKrve_coeffs = new double[3, 3] { { FPKrve[0, 0], FPKrve[0, 1], 0 }, { FPKrve[1, 0], FPKrve[1, 1], 0 }, { 0, 0, 0 } }; //TODO speed up this operation since only 2x2  are nonzero
                    //tensorOrder2 FPKtensor = new tensorOrder2(FPKrve_coeffs, Ei, ei);
                    //var FPKtensorProjected2 = FPKtensor.ProjectIn3DCartesianBasis(); //TODO: Speed up this operation.

                    var FPKtensorProjected2_check = new double[3, 3];
                    for (int j1 = 0; j1 < 3; j1++)
                    {
                        for (int j2 = 0; j2 < 3; j2++)
                        {
                            FPKtensorProjected2_check[j1, j2] += FPKrve[0, 0] * Ei[j1, 0] * ei[j2, 0] + FPKrve[0, 1] * Ei[j1, 0] * ei[j2, 1] +
                                                                 FPKrve[1, 0] * Ei[j1, 1] * ei[j2, 0] + FPKrve[1, 1] * Ei[j1, 1] * ei[j2, 1];
                        }
                    }


                    if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0) && ElementStiffnesses.performCalculations)
                    {
                        tensorOrder2 GLtensorProjected;
                        tensorOrder2 SPKtensorProjected = new tensorOrder2();
                        tensorOrder2 defGradTensor = new tensorOrder2();
                        var defGradTensorTr = new tensorOrder2();
                        tensorOrder2 FPKtensorProjected = new tensorOrder2();
                        tensorOrder2 SPKtensor = new tensorOrder2();

                        double[,] ch_F_3D = Calculate3DtensorFrom2Dcorrected_normaliseBothBasesCase(new double[,] { { 1, 0, 0 }, { 0, 1, 0, }, { 0, 0, 1 } }, Vector.CreateFromArray(g1), Vector.CreateFromArray(g2), Vector.CreateFromArray(a3), G_1, G_2, G_3);


                        bool moreChecks = false;
                        if (moreChecks)
                        {
                            //var GL_coeffs = new double[3, 3]
                            //{
                            //        {g1.DotProduct(g1)-G1.DotProduct(G1),g1.DotProduct(g2)-G1.DotProduct(G2),g1.DotProduct(a3)-G1.DotProduct(a3_init) },
                            //        {g2.DotProduct(g1)-G2.DotProduct(G1),g2.DotProduct(g2)-G2.DotProduct(G2),g2.DotProduct(a3)-G2.DotProduct(a3_init) },
                            //        {a3.DotProduct(g1)-a3_init.DotProduct(G1),a3.DotProduct(g2)-a3_init.DotProduct(G2),a3.DotProduct(a3)-a3_init.DotProduct(a3_init) },
                            //};

                            var GL_coeffs = new double[3, 3]
                            {
                                {g1.DotProduct(g1)-G1.DotProduct(G1),g1.DotProduct(g2)-G1.DotProduct(G2),g1.DotProduct(a3)-G1.DotProduct(a3_init) },
                                {g2.DotProduct(g1)-G2.DotProduct(G1),g2.DotProduct(g2)-G2.DotProduct(G2),g2.DotProduct(a3)-G2.DotProduct(a3_init) },
                                {a3.DotProduct(g1)-a3_init.DotProduct(G1),a3.DotProduct(g2)-a3_init.DotProduct(G2),a3.DotProduct(a3)-a3_init.DotProduct(a3_init) },
                            };

                            var corrections = new double[2, 2]
                            {
                            {-da3_dksi.Scale(z).DotProduct(da3_dksi.Scale(z))+da3_dksi_init.Scale(z).DotProduct(da3_dksi_init.Scale(z)),-da3_dksi.Scale(z).DotProduct(da3_dheta.Scale(z))+da3_dksi_init.Scale(z).DotProduct(da3_dheta_init.Scale(z)) },
                            {-da3_dheta.Scale(z).DotProduct(da3_dksi.Scale(z))+da3_dheta_init.Scale(z).DotProduct(da3_dksi_init.Scale(z)),-da3_dheta.Scale(z).DotProduct(da3_dheta.Scale(z))+da3_dheta_init.Scale(z).DotProduct(da3_dheta_init.Scale(z)) }
                            };

                            //for (int j1 = 0; j1 < 2; j1++)
                            //{
                            //    for (int j2 = 0; j2 < 2; j2++)
                            //    {
                            //        GL_coeffs[j1, j2] += corrections[j1, j2];
                            //    }
                            //}
                            tensorOrder2 GLtensor = new tensorOrder2(GL_coeffs);
                            GLtensor = GLtensor.Scale(0.5);

                            GLtensor.ReplaceBasisWithVector(G_1, G_2, G_3, true);
                            GLtensor.ReplaceBasisWithVector(G_1, G_2, G_3, false);

                            GLtensorProjected = GLtensor.ProjectIn3DCartesianBasis();

                            var materialaux_ = new ShellElasticMaterial2Dtransformationb()
                            {
                                PoissonRatio = material.PoissonRatio,
                                YoungModulus = material.YoungModulus,
                            }; 

                            materialaux_.TangentVectorV1 = G1;
                            materialaux_.TangentVectorV2 = G2;
                            materialaux_.NormalVectorV3 = a3_init;

                            materialaux_.UpdateMaterial(new double[] { GLtensor.coefficients[0, 0], GLtensor.coefficients[1, 1], 2 * GLtensor.coefficients[0, 1] });

                            var stresses = materialaux_.Stresses;

                            var SPK_coeffs = new double[3, 3]
                            {
                                {stresses[0], stresses[2], 0 },
                                {stresses[2], stresses[1], 0 },
                                {0,0,0 },
                            };

                            SPKtensor = new tensorOrder2(SPK_coeffs);

                            SPKtensor.ReplaceBasisWithVector(G1, G2, a3_init, true);
                            SPKtensor.ReplaceBasisWithVector(G1, G2, a3_init, false);

                            SPKtensorProjected = SPKtensor.ProjectIn3DCartesianBasis();

                            defGradTensor = new tensorOrder2(ch_F_3D);

                            defGradTensorTr = defGradTensor.Transpose();

                            var FPKtensor_ = SPKtensorProjected.SingleContract(defGradTensorTr);

                            FPKtensorProjected = FPKtensor_.ProjectIn3DCartesianBasis();



                        }

                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }

                        //ElementStiffnesses.ProccessVariable(28, /*dF2D_coefs_dr_vec  */ new double[] { F_rve[0, 0], F_rve[1, 1], F_rve[0, 1], F_rve[1, 0] }, false); //TODO: restore this maybe in CalculateStressesDevelop//.
                        ElementStiffnesses.ProccessVariable(29, /*dFPK2D_coefs_dr_vec*/ new double[] { FPKrve[0, 0], FPKrve[1, 1], FPKrve[0, 1], FPKrve[1, 0] }, false); //TODO: restore this maybe in CalculateStressesDevelop No this will be kept//.
                        //ElementStiffnesses.ProccessVariable(30, /*dFPK_3D_dr_vec*/ FPK_3D_vec, false); // TODO: probably we don't need this
                        //ElementStiffnesses.ProccessVariable(33, /*dFPK_3D_dr_vec*/ FPK_3D_vec, false); 

                        ElementStiffnesses.ProccessVariable(31, new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }, false);

                        ElementStiffnesses.ProccessVariable(32, new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }, false);

                        double[] dFPK_tensor_vec = { FPKtensorProjected.coefficients[0, 0], FPKtensorProjected.coefficients[1, 1], FPKtensorProjected.coefficients[2, 2], FPKtensorProjected.coefficients[0, 1], FPKtensorProjected.coefficients[1, 2], FPKtensorProjected.coefficients[2, 0], FPKtensorProjected.coefficients[0, 2], FPKtensorProjected.coefficients[1, 0], FPKtensorProjected.coefficients[2, 1] };// 
                        ElementStiffnesses.ProccessVariable(34, dFPK_tensor_vec, false);



                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                    }

                    for (int i = 0; i < elementControlPoints.Length; i++)
                    {
                        //var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                        //var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                        //var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                        //var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                        //var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);



                        //var a3r = new a3r();
                        var dksi_r = _shapeFunctions.DerivativeValuesKsi[i, j];
                        var dheta_r = _shapeFunctions.DerivativeValuesHeta[i, j];
                        //CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                        //    dheta_r, J1, ref a3r); //gia to a3r
                        var d2Ksi_dr = _shapeFunctions.SecondDerivativeValuesKsi[i, j];
                        var d2Heta_dr = _shapeFunctions.SecondDerivativeValuesHeta[i, j];
                        var d2KsiHeta_dr = _shapeFunctions.SecondDerivativeValuesKsiHeta[i, j];


                        ////5.24
                        //var da3tilde_dr = new double[3][];
                        //Calculate_da3tilde_dr(a1, a2, dksi_r, dheta_r, da3tilde_dr);

                        ////5.25
                        //double[] dnorma3_dr = new double[3];
                        //a3_tilde = Vector.CreateFromArray(CalculateTerm525(a3, J1, dnorma3_dr, da3tilde_dr));

                        ////5.30 b
                        //(Vector[] da3tilde_dksidr, Vector[] da3tilde_dhetadr) = Calculate_da3tilde_dksidr(a1r, a2r, a11r, a22r, a12r, a1, a2, a11, a22, a12);

                        ////5.31 b
                        //(double[] da3norm_dksidr, double[] da3norm_dhetadr) = Calculate_da3norm_dksidr(da3tilde_dksidr, da3tilde_dhetadr,
                        //    a3_tilde, da3tilde_dksi, da3tilde_dheta, da3tilde_dr, J1);

                        ////5.32 b
                        //(Vector[] da3_dksidr, Vector[] da3_dhetadr) = Calculate_da3_dksidr(da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksi, da3tilde_dheta,
                        //    dnorma3_dr, a3_tilde, da3norm_dksidr, da3norm_dhetadr, da3norm_dksi, da3norm_dheta, J1, da3tilde_dr);



                        //5.24
                        var da3tilde_dr = new double[3][];
                        Calculate_da3tilde_dr_Array(a1, a2, dksi_r, dheta_r, da3tilde_dr);

                        //5.25
                        var dnorma3_dr = new double[3];
                        a3_tilde = CalculateTerm525_Array(a3, J1, dnorma3_dr, da3tilde_dr);

                        //5.30 b
                        //(da3tilde_dksidrArray[i], da3tilde_dhetadrArray[i]) = Calculate_da3tilde_dksidr(a1rArray[i], a2rArray[i], a11rArray[i], a22rArray[i], a12rArray[i], a1, a2, a11, a22, a12);
                        (double[][] da3tilde_dksidr, double[][] da3tilde_dhetadr) = Calculate_da3tilde_dksidrDvelop(a1, a2, a11, a22, a12,
                            dksi_r, dheta_r, d2Ksi_dr, d2Heta_dr, d2KsiHeta_dr);


                        //5.31 b
                        (double[] da3norm_dksidr, double[] da3norm_dhetadr) = Calculate_da3norm_dksidr_Develop(da3tilde_dksidr, da3tilde_dhetadr,
                            a3_tilde, da3tilde_dksi, da3tilde_dheta, da3tilde_dr, J1);

                        //5.32 b
                        (double[][] da3_dksidr, double[][] da3_dhetadr) = Calculate_da3_dksidrDevelop(da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksi, da3tilde_dheta,
                            dnorma3_dr, a3_tilde, da3norm_dksidr, da3norm_dhetadr, da3norm_dksi, da3norm_dheta, J1, da3tilde_dr);





                        for (int r1 = 0; r1 < 3; r1++)
                        {
                            //(31)
                            //Vector dg1_dr = a1rArray[i].GetColumn(r1) + da3_dksidrArray[i][r1] * z;
                            double[] dg1_dr = new double[3] { da3_dksidr[r1][0] * z, da3_dksidr[r1][1] * z, da3_dksidr[r1][2] * z };
                            dg1_dr[r1] += dksi_r;

                            //Vector dg2_dr = a2rArray[i].GetColumn(r1) + da3_dhetadrArray[i][r1] * z;
                            double[] dg2_dr = new double[3] { da3_dhetadr[r1][0] * z, da3_dhetadr[r1][1] * z, da3_dhetadr[r1][2] * z };
                            dg2_dr[r1] += dheta_r;



                            //Vector dg3_dr = a3r. ....

                            //(39)


                            double[] da3_dr = CalculateDerivativeOfVectorNormalisedArray(a3_tilde, da3tilde_dr[r1]);



                            //tensorOrder2 dF3Dtensor_dr = new tensorOrder2()
                            //{
                            //    basis1 = new double[,] { { dg1_dr[0], dg2_dr[0], da3_dr[0] }, { dg1_dr[1], dg2_dr[1], da3_dr[1] }, { dg1_dr[2], dg2_dr[2], da3_dr[2] } },
                            //    basis2 = new double[,] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } },
                            //    coefficients = new double[,] { { 1, 0, 0 }, { 0, 1, 0, }, { 0, 0, 1 } }
                            //};

                            //tensorOrder2 dF3Dtensor_dr_Tr = dF3Dtensor_dr.Transpose();


                            Array.Clear(dF3Dtensor_dr_Tr_check, 0, 9);
                            for (int j1 = 0; j1 < 3; j1++)
                            {
                                for (int j2 = 0; j2 < 3; j2++)
                                {
                                    dF3Dtensor_dr_Tr_check[j2, j1] += dg1_dr[j1] * G_1[j2] + dg2_dr[j1] * G_2[j2] + da3_dr[j1] * G_3[j2];
                                }
                            }


                            //forcesDevelop_v3[3 * i + r1] += FPKtensorProjected2.doubleContract(dF3Dtensor_dr_Tr) * wfactor * w;// this is going to be replaced by simplified tensor expressions speed
                            for (int j1 = 0; j1 < 3; j1++)
                            {
                                for (int j2 = 0; j2 < 3; j2++)
                                {
                                    forcesDevelop_v3[3 * i + r1] += FPKtensorProjected2_check[j1, j2] * dF3Dtensor_dr_Tr_check[j1, j2] * wfactor * w;
                                }
                            }




                            if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0) && ElementStiffnesses.performCalculations)
                            {
                                double[,] dF_3D_dr = new double[3, 3] { { dg1_dr[0]*G_1[0]+dg2_dr[0]*G_2[0], dg1_dr[0]*G_1[1]+dg2_dr[0]*G_2[1], dg1_dr[0]*G_1[2]+dg2_dr[0]*G_2[2] },
                                                                 { dg1_dr[1]*G_1[0]+dg2_dr[1]*G_2[0], dg1_dr[1]*G_1[1]+dg2_dr[1]*G_2[1], dg1_dr[1]*G_1[2]+dg2_dr[1]*G_2[2] },
                                                                 { dg1_dr[2]*G_1[0]+dg2_dr[2]*G_2[0], dg1_dr[2]*G_1[1]+dg2_dr[2]*G_2[1], dg1_dr[2]*G_1[2]+dg2_dr[2]*G_2[2] }, };
                                double[] dF_3D_dr_vec = { dF_3D_dr[0, 0], dF_3D_dr[1, 1], dF_3D_dr[2, 2], dF_3D_dr[0, 1], dF_3D_dr[1, 2], dF_3D_dr[2, 0], dF_3D_dr[0, 2], dF_3D_dr[1, 0], dF_3D_dr[2, 1] };

                                ElementStiffnesses.ProccessVariable(18, dF_3D_dr_vec, true, 3 * i + r1);
                            }


                            if ((j == ElementStiffnesses.gpNumberToCheck) && (i == 0) && (r1 == 1) && (i1 == 0) && ElementStiffnesses.performCalculations)
                            {
                                if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                                double[,] dF_3D_dr = new double[3, 3] { { dg1_dr[0]*G_1[0]+dg2_dr[0]*G_2[0], dg1_dr[0]*G_1[1]+dg2_dr[0]*G_2[1], dg1_dr[0]*G_1[2]+dg2_dr[0]*G_2[2] },
                                                                 { dg1_dr[1]*G_1[0]+dg2_dr[1]*G_2[0], dg1_dr[1]*G_1[1]+dg2_dr[1]*G_2[1], dg1_dr[1]*G_1[2]+dg2_dr[1]*G_2[2] },
                                                                 { dg1_dr[2]*G_1[0]+dg2_dr[2]*G_2[0], dg1_dr[2]*G_1[1]+dg2_dr[2]*G_2[1], dg1_dr[2]*G_1[2]+dg2_dr[2]*G_2[2] }, };
                                double[] dF_3D_dr_vec = { dF_3D_dr[0, 0], dF_3D_dr[1, 1], dF_3D_dr[2, 2], dF_3D_dr[0, 1], dF_3D_dr[1, 2], dF_3D_dr[2, 0], dF_3D_dr[0, 2], dF_3D_dr[1, 0], dF_3D_dr[2, 1] };
                                ElementStiffnesses.ProccessVariable(27, dF_3D_dr_vec, false);

                                if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                            }

                            if (j == ElementStiffnesses.gpNumberToCheck)
                            {

                                //ElementStiffnesses.ProccessVariable(11, new double[1] { dnorma3_dr[r1] }, true, 3 * i + r1);
                                //ElementStiffnesses.ProccessVariable(12, da3tilde_dksidr[r1].CopyToArray(), true, 3 * i + r1);
                                //ElementStiffnesses.ProccessVariable(13, da3tilde_dhetadr[r1].CopyToArray(), true, 3 * i + r1);
                                //ElementStiffnesses.ProccessVariable(14, new double[1] { da3norm_dksidr[r1] }, true, 3 * i + r1);
                                //ElementStiffnesses.ProccessVariable(15, new double[1] { da3norm_dhetadr[r1] }, true, 3 * i + r1);
                                ElementStiffnesses.ProccessVariable(16, da3_dksidr[r1], true, 3 * i + r1);
                                ElementStiffnesses.ProccessVariable(17, da3_dhetadr[r1], true, 3 * i + r1);

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

            
             return forcesDevelop_v3;
            


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

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {


            _solution = localDisplacements;

            var newControlPoints = CurrentControlPoint(_controlPoints);
            var gaussPoints = materialsAtThicknessGP.Keys.ToArray();

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





                #region develop formulation
                var thicknessGPoints = thicknessIntegrationPoints[gaussPoints[j]];
                var materialpoint = materialsAtThicknessGP[gaussPoints[j]][thicknessGPoints[0]];
                //var transformations = new ShellElasticMaterial2DtransformationbDefGrad() { YoungModulus = materialpoint.YoungModulus, PoissonRatio = materialpoint.PoissonRatio };



                var a11 = surfaceBasisVectorDerivative1;
                var a22 = surfaceBasisVectorDerivative2;
                var a12 = surfaceBasisVectorDerivative12;
                var a1 =  surfaceBasisVector1;
                var a2 =  surfaceBasisVector2;
                var a3 =  surfaceBasisVector3; // einai to mono pou einai normalised
                double[] a3_tilde = new double[] { a3[0] * (J1), a3[1] * (J1), a3[2] * (J1) };

                (double[] da3tilde_dksi, double[] da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, double[] da3_dksi, double[] da3_dheta) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1, a2, a11, a22, a12, a3, J1);

                //if (j == ElementStiffnesses.gpNumberToCheck)
                //{
                //    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                //    ElementStiffnesses.ProccessVariable(11, new double[] { J1 }, false);
                //    ElementStiffnesses.ProccessVariable(12, da3tilde_dksi.CopyToArray(), false);
                //    ElementStiffnesses.ProccessVariable(13, da3tilde_dheta.CopyToArray(), false);
                //    ElementStiffnesses.ProccessVariable(14, new double[] { da3norm_dksi }, false);
                //    ElementStiffnesses.ProccessVariable(15, new double[] { da3norm_dheta }, false);
                //    ElementStiffnesses.ProccessVariable(16, da3_dksi.CopyToArray(), false);
                //    ElementStiffnesses.ProccessVariable(17, da3_dheta.CopyToArray(), false);
                //    if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                //}

                #region original Config
                var originalControlPoints = _controlPoints.ToArray();
                var originalHessianMatrix = CalculateHessian(originalControlPoints, _shapeFunctions, j);
                double[,] jacobian_init = new double[3, 3];
                CalculateJacobian(originalControlPoints, _shapeFunctions, j, jacobian_init);
                var a1_init = CalculateSurfaceBasisVector(jacobian_init, 0);
                var a2_init = CalculateSurfaceBasisVector(jacobian_init, 1);
                var a3_init = new double[3]
                {
                    a1_init[1] * a2_init[2] - a1_init[2] * a2_init[1],
                    a1_init[2] * a2_init[0] - a1_init[0] * a2_init[2],
                    a1_init[0] * a2_init[1] - a1_init[1] * a2_init[0]
                };

                var J1_init = Math.Sqrt(a3_init[0] * a3_init[0] +
                                        a3_init[1] * a3_init[1] +
                                        a3_init[2] * a3_init[2]);

                a3_init[0] /= J1_init;
                a3_init[1] /= J1_init;
                a3_init[2] /= J1_init;

                var a11_init = CalculateSurfaceBasisVector(originalHessianMatrix, 0);
                var a22_init = CalculateSurfaceBasisVector(originalHessianMatrix, 1);
                var a12_init = CalculateSurfaceBasisVector(originalHessianMatrix, 2);
                #endregion

                (double[] da3tilde_dksi_init, double[] da3tilde_dheta_init, double da3norm_dksi_init, double da3norm_dheta_init, double[] da3_dksi_init, double[] da3_dheta_init) =
                    Calculate_da3tilde_dksi_524_525_526_b(a1_init, a2_init, a11_init, a22_init, a12_init, a3_init, J1_init);
                var thicknessPoints = thicknessIntegrationPoints[gaussPoints[j]];
                for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                {
                    var materialDevelop = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    var thicknessPoint = thicknessPoints[i1];
                    var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]]; // TODO: comment out mallon afto..//.
                    var w = thicknessPoint.WeightFactor;
                    var z = thicknessPoint.Zeta;

                    //Vector G1 = a1_init + da3_dksi_init.Scale(z);
                    //Vector G2 = a2_init + da3_dheta_init.Scale(z);
                    double[] G1 = new double[] { a1_init[0] + da3_dksi_init[0] * z, a1_init[1] + da3_dksi_init[1] * z, a1_init[2] + da3_dksi_init[2] * z };
                    double[] G2 = new double[] { a2_init[0] + da3_dheta_init[0] * z, a2_init[1] + da3_dheta_init[1] * z, a2_init[2] + da3_dheta_init[2] * z };

                    (double[] G_1, double[] G_2, double[] G_3) = CalculateContravariants(G1, G2, a3_init);

                    //Vector g1 = a1 + da3_dksi.Scale(z);
                    //Vector g2 = a2 + da3_dheta.Scale(z);
                    double[] g1 = new double[] { a1[0] + da3_dksi[0] * z, a1[1] + da3_dksi[1] * z, a1[2] + da3_dksi[2] * z };
                    double[] g2 = new double[] { a2[0] + da3_dheta[0] * z, a2[1] + da3_dheta[1] * z, a2[2] + da3_dheta[2] * z };



                    if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0) && ElementStiffnesses.performCalculations)
                    {
                        double[,] F_3D = new double[3, 3] { { g1[0]*G_1[0]+g2[0]*G_2[0], g1[0]*G_1[1]+g2[0]*G_2[1], g1[0]*G_1[2]+g2[0]*G_2[2] },
                                                            { g1[1]*G_1[0]+g2[1]*G_2[0], g1[1]*G_1[1]+g2[1]*G_2[1], g1[1]*G_1[2]+g2[1]*G_2[2] },
                                                            { g1[2]*G_1[0]+g2[2]*G_2[0], g1[2]*G_1[1]+g2[2]*G_2[1], g1[2]*G_1[2]+g2[2]*G_2[2] },
                        };
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; }
                        double[] F_3D_vec = { F_3D[0, 0], F_3D[1, 1], F_3D[2, 2], F_3D[0, 1], F_3D[1, 2], F_3D[2, 0], F_3D[0, 2], F_3D[1, 0], F_3D[2, 1] };//  .
                        ElementStiffnesses.ProccessVariable(18, F_3D_vec, false);
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                    }


                    (double[,] ei, double[,] F_rve) = CalculateTransformationsV2_FrveOnly(g1, g2, a3, G1, G2, a3_init, G_1, G_2, G_3);


                    materialDevelop.UpdateMaterial(new double[] { F_rve[0, 0], F_rve[1, 1], F_rve[0, 1], F_rve[1, 0] });


                    if ((j == ElementStiffnesses.gpNumberToCheck) && (i1 == 0) && ElementStiffnesses.performCalculations)
                    {
                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = true; } // afou tha eimaste entos twn stresses afto den xreiazetai poia
                                                                                                                    // dioti to Stresses kaleitai entos saveVariationStates=true;

                        ElementStiffnesses.ProccessVariable(28, /*dF2D_coefs_dr_vec  */ new double[] { F_rve[0, 0], F_rve[1, 1], F_rve[0, 1], F_rve[1, 0] }, false);

                        //ElementStiffnesses.ProccessVariable(29, /*dFPK2D_coefs_dr_vec*/ new double[] { FPK_2D[0, 0], FPK_2D[1, 1], FPK_2D[0, 1], FPK_2D[1, 0] }, false);
                        //ElementStiffnesses.ProccessVariable(30, /*dFPK_3D_dr_vec*/ FPK_3D_vec, false);
                        //ElementStiffnesses.ProccessVariable(33, /*dFPK_3D_dr_vec*/ FPK_3D_vec, false);

                        ElementStiffnesses.ProccessVariable(31, new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }, false); // TODO: mporei na mpei entos CalculateTransformationsV2only

                        ElementStiffnesses.ProccessVariable(32, new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }, false);

                        //ElementStiffnesses.ProccessVariable(34, dFPK_tensor_vec, false);



                        if (ElementStiffnesses.saveForcesState1) { ElementStiffnesses.saveVariationStates = false; }
                    }


                }

                #endregion

            }
            return new Tuple<double[], double[]>(new double[0], new double[0]);
        }

        #region Calculate stresses (and Forces) supportive
        private (double[] da3tilde_dksi, double[] da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, double[] da3_dksi, double[] da3_dheta)
            Calculate_da3tilde_dksi_524_525_526_b(double[] a1, double[] a2, double[] a11, double[] a22, double[] a12, double[] a3, double normA3)
        {
            //var da3tilde_dksi = a11.CrossProduct(a2) + a1.CrossProduct(a12);
            double[] da3tilde_dksi = new double[]
            {
                a11[1]*a2[2] - a11[2]*a2[1] - a12[1]*a1[2] + a12[2]*a1[1],
                a11[2]*a2[0] - a11[0]*a2[2] + a12[0]*a1[2] - a12[2]*a1[0],
                a11[0]*a2[1] - a11[1]*a2[0] - a12[0]*a1[1] + a12[1]*a1[0]
            };

            //var da3tilde_dheta = a12.CrossProduct(a2) + a1.CrossProduct(a22);//1
            double[] da3tilde_dheta = new double[]
            {
                a12[1]*a2[2] - a12[2]*a2[1] - a22[1]*a1[2] + a22[2]*a1[1],
                 a12[2]*a2[0] - a12[0]*a2[2] + a22[0]*a1[2] - a22[2]*a1[0],
                 a12[0]*a2[1] - a12[1]*a2[0] - a22[0]*a1[1] + a22[1]*a1[0],
            };

            var da3norm_dksi = a3.DotProduct(da3tilde_dksi);
            var da3norm_dheta = a3.DotProduct(da3tilde_dheta);

            double scaleFactor = (double)1 / normA3;
            //var da3_dksi = da3tilde_dksi.Scale(scaleFactor) - a3.Scale(da3norm_dksi).Scale(scaleFactor);
            var da3_dksi = new double[]{ da3tilde_dksi[0]*scaleFactor - a3[0]*da3norm_dksi*scaleFactor,
                da3tilde_dksi[1]*scaleFactor - a3[1]*da3norm_dksi*scaleFactor,
                da3tilde_dksi[2]*scaleFactor - a3[2]*da3norm_dksi*scaleFactor };

            //var da3_dheta = da3tilde_dheta.Scale(scaleFactor) - a3.Scale(da3norm_dheta).Scale(scaleFactor);
            var da3_dheta = new double[] { da3tilde_dheta[0]*scaleFactor - a3[0]*da3norm_dheta*scaleFactor,
                da3tilde_dheta[1]*scaleFactor - a3[1]*da3norm_dheta*scaleFactor,
                da3tilde_dheta[2]*scaleFactor - a3[2]*da3norm_dheta*scaleFactor };




            return (da3tilde_dksi, da3tilde_dheta, da3norm_dksi, da3norm_dheta, da3_dksi, da3_dheta);

        }

        private (double[] G_1, double[] G_2, double[] G_3) CalculateContravariants(double[] g1, double[] g2, double[] a3)
        {
            var auxMatrix1 = Matrix.CreateZero(3, 3);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = g1.DotProduct(g1);
            auxMatrix1[0, 1] = g1.DotProduct(g2);
            auxMatrix1[0, 2] = g1.DotProduct(a3);
            auxMatrix1[1, 0] = g2.DotProduct(g1);
            auxMatrix1[1, 1] = g2.DotProduct(g2);
            auxMatrix1[1, 2] = g2.DotProduct(a3);
            auxMatrix1[2, 0] = a3.DotProduct(g1);
            auxMatrix1[2, 1] = a3.DotProduct(g2);
            auxMatrix1[2, 2] = a3.DotProduct(a3);
            Matrix inverse = auxMatrix1.Invert(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)
                                                  //TODO: auxMatrix1.Invert2x2AndDeterminant(1e-20) for bad geometry

            //Contravariant base vectors
            double[][] G_i = new double[3][];
            for (int i1 = 0; i1 < 3; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * g1[i2] + inverse[i1, 1] * g2[i2] + inverse[i1, 2] * a3[i2];
                }
            }

            return (G_i[0], G_i[1], G_i[2]);
        }

        public (double[,], double[,]) CalculateTransformationsV2_FrveOnly(double[] g1, double[] g2, double[] g3, double[] G1, double[] G2, double[] G3, double[] G_1, double[] G_2, double[] G_3)
        {
            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;
            double[,] tgi = new double[3, 3] { { g1[0], g2[0], g3[0] }, { g1[1], g2[1], g3[1] }, { g1[2], g2[2], g3[2] } };
            double[,] Gi = new double[3, 3] { { G1[0], G2[0], G3[0] }, { G1[1], G2[1], G3[1] }, { G1[2], G2[2], G3[2] } };
            double[,] G_i = new double[3, 3] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } };


            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = tgi[0, i1] * tgi[0, i1] + tgi[1, i1] * tgi[1, i1] + tgi[2, i1] * tgi[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    tgi[i2, i1] = tgi[i2, i1] / norm;
                }

            }
            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = Gi[0, i1] * Gi[0, i1] + Gi[1, i1] * Gi[1, i1] + Gi[2, i1] * Gi[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    Gi[i2, i1] = Gi[i2, i1] / norm;
                }

            }

            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = G_i[0, i1] * G_i[0, i1] + G_i[1, i1] * G_i[1, i1] + G_i[2, i1] * G_i[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i2, i1] = G_i[i2, i1] / norm;
                }

            }

            OrthogonaliseBasisMembranePart(tgi);
            OrthogonaliseBasisMembranePart(Gi);
            OrthogonaliseBasisMembranePart(G_i);

            var ei = tgi;
            var Ei = Gi;


            double[] g1__ei = new double[] { g1[0] * ei[0, 0] + g1[1] * ei[1, 0] + g1[2] * ei[2, 0], g1[0] * ei[0, 1] + g1[1] * ei[1, 1] + g1[2] * ei[2, 1], g1[0] * ei[0, 2] + g1[1] * ei[1, 2] + g1[2] * ei[2, 2] };
            double[] g2__ei = new double[] { g2[0] * ei[0, 0] + g2[1] * ei[1, 0] + g2[2] * ei[2, 0], g2[0] * ei[0, 1] + g2[1] * ei[1, 1] + g2[2] * ei[2, 1], g2[0] * ei[0, 2] + g2[1] * ei[1, 2] + g2[2] * ei[2, 2] };
            double[] g3__ei = new double[] { g3[0] * ei[0, 0] + g3[1] * ei[1, 0] + g3[2] * ei[2, 0], g3[0] * ei[0, 1] + g3[1] * ei[1, 1] + g3[2] * ei[2, 1], g3[0] * ei[0, 2] + g3[1] * ei[1, 2] + g3[2] * ei[2, 2] };


            double[] G_1__Ei = new double[] { G_1[0] * Ei[0, 0] + G_1[1] * Ei[1, 0] + G_1[2] * Ei[2, 0], G_1[0] * Ei[0, 1] + G_1[1] * Ei[1, 1] + G_1[2] * Ei[2, 1], G_1[0] * Ei[0, 2] + G_1[1] * Ei[1, 2] + G_1[2] * Ei[2, 2] };
            double[] G_2__Ei = new double[] { G_2[0] * Ei[0, 0] + G_2[1] * Ei[1, 0] + G_2[2] * Ei[2, 0], G_2[0] * Ei[0, 1] + G_2[1] * Ei[1, 1] + G_2[2] * Ei[2, 1], G_2[0] * Ei[0, 2] + G_2[1] * Ei[1, 2] + G_2[2] * Ei[2, 2] };
            double[] G_3__Ei = new double[] { G_3[0] * Ei[0, 0] + G_3[1] * Ei[1, 0] + G_3[2] * Ei[2, 0], G_3[0] * Ei[0, 1] + G_3[1] * Ei[1, 1] + G_3[2] * Ei[2, 1], G_3[0] * Ei[0, 2] + G_3[1] * Ei[1, 2] + G_3[2] * Ei[2, 2] };


            double[,] F = CaclculateDefGrad3D(g1__ei, g2__ei, g3__ei, G_1__Ei, G_2__Ei, G_3__Ei);





            double[,] F_rve = new double[,]
            {
                {F[0,0], F[0,1] },
                {F[1,0], F[1,1] }
            };




            return (ei, F_rve);
        }

        private void OrthogonaliseBasisMembranePart(double[,] tgi)
        {
            double[] E1 = new double[] { tgi[0, 0], tgi[1, 0], tgi[2, 0] };
            double[] E2 = new double[] { tgi[0, 1], tgi[1, 1], tgi[2, 1] };

            double E1_dot_E2 = 0;
            for (int i1 = 0; i1 < 3; i1++) { E1_dot_E2 += E1[i1] * E2[i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = E1_dot_E2 * E1[i1]; }

            double[] E2ortho = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { E2ortho[i1] = E2[i1] - projection[i1]; }
            double norm1 = (Vector.CreateFromArray(E2ortho)).Norm2();
            for (int i1 = 0; i1 < 3; i1++) { E2[i1] = E2ortho[i1] / norm1; }

            tgi[0, 1] = E2[0];
            tgi[1, 1] = E2[1];
            tgi[2, 1] = E2[2];

            //
        }

        private double[,] CaclculateDefGrad3D(double[] g1__ei, double[] g2__ei, double[] g3__ei, double[] G_1__Ei, double[] G_2__Ei, double[] G_3__Ei)
        {

            //TODO: only the 2x2 is needed.
            double[,] F = new double[3, 3] {
                { g1__ei[0]*G_1__Ei[0]+g2__ei[0]*G_2__Ei[0]+g3__ei[0]*G_3__Ei[0], g1__ei[0]*G_1__Ei[1]+g2__ei[0]*G_2__Ei[1]+g3__ei[0]*G_3__Ei[1], g1__ei[0]*G_1__Ei[2]+g2__ei[0]*G_2__Ei[2]+g3__ei[0]*G_3__Ei[2] },
                { g1__ei[1]*G_1__Ei[0]+g2__ei[1]*G_2__Ei[0]+g3__ei[1]*G_3__Ei[0], g1__ei[1]*G_1__Ei[1]+g2__ei[1]*G_2__Ei[1]+g3__ei[1]*G_3__Ei[1], g1__ei[1]*G_1__Ei[2]+g2__ei[1]*G_2__Ei[2]+g3__ei[1]*G_3__Ei[2] },
                { g1__ei[2]*G_1__Ei[0]+g2__ei[2]*G_2__Ei[0]+g3__ei[2]*G_3__Ei[0], g1__ei[2]*G_1__Ei[1]+g2__ei[2]*G_2__Ei[1]+g3__ei[2]*G_3__Ei[1], g1__ei[2]*G_1__Ei[2]+g2__ei[2]*G_2__Ei[2]+g3__ei[2]*G_3__Ei[2] },
                        };
            return F;
        }

        private void NormaliseAndOrthogonaliseBasis(double[,] tgi)
        {
            for (int i1 = 0; i1 < 3; i1++)
            {
                double norm = tgi[0, i1] * tgi[0, i1] + tgi[1, i1] * tgi[1, i1] + tgi[2, i1] * tgi[2, i1];
                norm = Math.Sqrt(norm);
                for (int i2 = 0; i2 < 3; i2++)
                {
                    tgi[i2, i1] = tgi[i2, i1] / norm;
                }

            }

            OrthogonaliseBasisMembranePart(tgi);

        }

        private double[,] Calculate3DtensorFrom2Dcorrected_normaliseBothBasesCase(double[,] eye3, Vector dg1_dr, Vector dg2_dr, Vector da3_dr, double[] G_1, double[] G_2, double[] G_3)
        {
            //double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;

            #region create and normalise ei
            double[,] ei = new double[3, 3];
            double norm_e1 = dg1_dr[0] * dg1_dr[0] + dg1_dr[1] * dg1_dr[1] + dg1_dr[2] * dg1_dr[2];
            norm_e1 = Math.Sqrt(norm_e1);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 0] = dg1_dr[i2] / norm_e1;
            }

            double norm_e2 = dg2_dr[0] * dg2_dr[0] + dg2_dr[1] * dg2_dr[1] + dg2_dr[2] * dg2_dr[2];
            norm_e2 = Math.Sqrt(norm_e2);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 1] = dg2_dr[i2] / norm_e2;
            }

            double norm_e3 = da3_dr[0] * da3_dr[0] + da3_dr[1] * da3_dr[1] + da3_dr[2] * da3_dr[2];
            norm_e3 = Math.Sqrt(norm_e3);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 2] = da3_dr[i2] / norm_e3;
            }
            #endregion

            #region adapt FPK_2D for normalisation of basis vectors

            double coef1 = 0;
            double coef2 = 0;
            double[,] FPK_2D_in_normalised = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                if (i1 == 0) { coef1 = norm_e1; }
                else if (i1 == 1) { coef1 = norm_e2; }
                else if (i1 == 2) { coef1 = norm_e3; }
                for (int i2 = 0; i2 < 3; i2++)
                {
                    //if (i2 == 0) { coef2 = norm_e1; }
                    //else if (i2 == 1) { coef2 = norm_e2; }

                    //FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1 * coef2;
                    FPK_2D_in_normalised[i1, i2] = eye3[i1, i2] * coef1;// * coef2;
                }
            }
            #endregion

            #region create and normalise Ei
            double[,] Ei = new double[3, 3];
            double norm_E1 = G_1[0] * G_1[0] + G_1[1] * G_1[1] + G_1[2] * G_1[2];
            norm_E1 = Math.Sqrt(norm_E1);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 0] = G_1[i2] / norm_E1;
            }

            double norm_E2 = G_2[0] * G_2[0] + G_2[1] * G_2[1] + G_2[2] * G_2[2];
            norm_E2 = Math.Sqrt(norm_E2);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 1] = G_2[i2] / norm_E2;
            }

            double norm_E3 = G_3[0] * G_3[0] + G_3[1] * G_3[1] + G_3[2] * G_3[2];
            norm_E3 = Math.Sqrt(norm_E3);
            for (int i2 = 0; i2 < 3; i2++)
            {
                Ei[i2, 2] = G_3[i2] / norm_E3;
            }
            #endregion

            #region adapt FPK_2D for normalisation of basis vectors

            //double coef1 = 0;
            //double coef2 = 0;
            //double[,] FPK_2D_in_normalised = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {

                for (int i2 = 0; i2 < 3; i2++)
                {
                    if (i2 == 0) { coef2 = norm_E1; }
                    else if (i2 == 1) { coef2 = norm_E2; }
                    else if (i2 == 2) { coef2 = norm_E3; }

                    //FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1 * coef2;
                    FPK_2D_in_normalised[i1, i2] = FPK_2D_in_normalised[i1, i2] * coef2;
                }
            }
            #endregion

            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;
            var cartes_to_Gi = CalculateRotationMatrix(Ei, eye);
            var cartes_to_tgi = CalculateRotationMatrix(ei, eye);

            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPK_2D_in_normalised, cartes_to_Gi, cartes_to_tgi);// 1);



            return FPK_3D;
        }

        private static void Calculate_da3tilde_dr_Array(double[] surfaceBasisVector1, double[] surfaceBasisVector2, double dksi_r,
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

        private (double[][] da3tilde_dksidr, double[][] da3tilde_dhetadr) Calculate_da3tilde_dksidrDvelop(
            double[] a1, double[] a2, double[] a11, double[] a22, double[] a12,
            double dksi_r, double dheta_r, double d2Ksi_dr, double d2Heta_dr, double d2KsiHeta_dr)
        {
            double[][] da3tilde_dksidr = new double[3][];
            double[][] da3tilde_dhetadr = new double[3][];

            //for (int r1 = 0; r1 < 3; r1++)
            //{
            //    da3tilde_dksidr[r1] = a11r.GetColumn(r1).CrossProduct(a2) + a11.CrossProduct(a2r.GetColumn(r1)) + a1r.GetColumn(r1).CrossProduct(a12) + a1.CrossProduct(a12r.GetColumn(r1));
            //    da3tilde_dhetadr[r1] = a12r.GetColumn(r1).CrossProduct(a2) + a12.CrossProduct(a2r.GetColumn(r1)) + a1r.GetColumn(r1).CrossProduct(a22) + a1.CrossProduct(a22r.GetColumn(r1));
            //}

            da3tilde_dksidr[0] = new double[3] { 0, a1[2] * d2KsiHeta_dr - a2[2] * d2Ksi_dr + a11[2] * dheta_r - a12[2] * dksi_r, a2[1] * d2Ksi_dr - a1[1] * d2KsiHeta_dr - a11[1] * dheta_r + a12[1] * dksi_r };
            da3tilde_dksidr[1] = new double[3] { a2[2] * d2Ksi_dr - a1[2] * d2KsiHeta_dr - a11[2] * dheta_r + a12[2] * dksi_r, 0, a1[0] * d2KsiHeta_dr - a2[0] * d2Ksi_dr + a11[0] * dheta_r - a12[0] * dksi_r };
            da3tilde_dksidr[2] = new double[3] { a1[1] * d2KsiHeta_dr - a2[1] * d2Ksi_dr + a11[1] * dheta_r - a12[1] * dksi_r, a2[0] * d2Ksi_dr - a1[0] * d2KsiHeta_dr - a11[0] * dheta_r + a12[0] * dksi_r, 0 };

            da3tilde_dhetadr[0] = new double[3] { 0, a1[2] * d2Heta_dr - a2[2] * d2KsiHeta_dr + a12[2] * dheta_r - a22[2] * dksi_r, a2[1] * d2KsiHeta_dr - a1[1] * d2Heta_dr - a12[1] * dheta_r + a22[1] * dksi_r };
            da3tilde_dhetadr[1] = new double[3] { a2[2] * d2KsiHeta_dr - a1[2] * d2Heta_dr - a12[2] * dheta_r + a22[2] * dksi_r, 0, a1[0] * d2Heta_dr - a2[0] * d2KsiHeta_dr + a12[0] * dheta_r - a22[0] * dksi_r };
            da3tilde_dhetadr[2] = new double[3] { a1[1] * d2Heta_dr - a2[1] * d2KsiHeta_dr + a12[1] * dheta_r - a22[1] * dksi_r, a2[0] * d2KsiHeta_dr - a1[0] * d2Heta_dr - a12[0] * dheta_r + a22[0] * dksi_r, 0 };



            return (da3tilde_dksidr, da3tilde_dhetadr);
        }

        private (double[] da3norm_dksidr, double[] da3norm_dhetadr) Calculate_da3norm_dksidr_Develop(double[][] da3tilde_dksidr, double[][] da3tilde_dhetadr, double[] a3_tilde, double[] da3tilde_dksi, double[] da3tilde_dheta, double[][] da3tilde_dr, double J1)
        {
            var da3norm_dksidr = new double[3];
            var da3norm_dhetadr = new double[3];

            for (int r1 = 0; r1 < 3; r1++)
            {
                //double firstNumerator = da3tilde_dksidr[r1].DotProduct(a3_tilde) + da3tilde_dksi.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double firstNumerator = da3tilde_dksidr[r1][0] * a3_tilde[0] + da3tilde_dksidr[r1][1] * a3_tilde[1] + da3tilde_dksidr[r1][2] * a3_tilde[2] +
                                        da3tilde_dksi[0] * da3tilde_dr[r1][0] + da3tilde_dksi[1] * da3tilde_dr[r1][1] + da3tilde_dksi[2] * da3tilde_dr[r1][2];

                double firstDenominator = J1;

                //double secondNumerator = da3tilde_dksi.DotProduct(a3_tilde) * da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray());
                double secondNumerator = (da3tilde_dksi[0] * a3_tilde[0] + da3tilde_dksi[1] * a3_tilde[1] + da3tilde_dksi[2] * a3_tilde[2]) *
                    (da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] + da3tilde_dr[r1][2] * a3_tilde[2]);

                double secondDenominator = Math.Pow(J1, 3);

                da3norm_dksidr[r1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);

            }

            for (int r1 = 0; r1 < 3; r1++)
            {
                //double firstNumerator = da3tilde_dhetadr[r1].DotProduct(a3_tilde) + da3tilde_dheta.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double firstNumerator = da3tilde_dhetadr[r1][0] * a3_tilde[0] + da3tilde_dhetadr[r1][1] * a3_tilde[1] + da3tilde_dhetadr[r1][2] * a3_tilde[2] +
                    +da3tilde_dheta[0] * da3tilde_dr[r1][0] + da3tilde_dheta[1] * da3tilde_dr[r1][1] + da3tilde_dheta[2] * da3tilde_dr[r1][2];

                double firstDenominator = J1;

                double secondNumerator = da3tilde_dheta.DotProduct(a3_tilde) * da3tilde_dr[r1].DotProduct(a3_tilde);
                double secondDenominator = Math.Pow(J1, 3);

                da3norm_dhetadr[r1] = (firstNumerator / firstDenominator) - (secondNumerator / secondDenominator);

            }

            return (da3norm_dksidr, da3norm_dhetadr);
        }

        private static double[] CalculateTerm525_Array(double[] surfaceBasisVector3, double J1, double[] dnorma3_dr,
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

        private (double[][] da3_dksidr, double[][] da3_dhetadr) Calculate_da3_dksidrDevelop(double[][] da3tilde_dksidr, double[][] da3tilde_dhetadr, double[] da3tilde_dksi,
            double[] da3tilde_dheta, double[] dnorma3_dr, double[] a3_tilde, double[] da3norm_dksidr, double[] da3norm_dhetadr, double da3norm_dksi,
            double da3norm_dheta, double J1, double[][] da3tilde_dr)
        {
            double[][] da3_dksidr = new double[3][];
            double[][] da3_dhetadr = new double[3][];

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

                da3_dksidr[r1] = new double[]
                {
                    firstVec_0 + secondVec_0 + thirdVec_0 + fourthVec_0 + fifthvector_0,
                    firstVec_1 + secondVec_1 + thirdVec_1 + fourthVec_1 + fifthvector_1,
                    firstVec_2 + secondVec_2 + thirdVec_2 + fourthVec_2 + fifthvector_2,
                };

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

                da3_dhetadr[r1] = new double[]
                {
                    firstVec_0 + secondVec_0 + thirdVec_0 + fourthVec_0 + fifthvector_0,
                    firstVec_1 + secondVec_1 + thirdVec_1 + fourthVec_1 + fifthvector_1,
                    firstVec_2 + secondVec_2 + thirdVec_2 + fourthVec_2 + fifthvector_2,
                };


            }

            return (da3_dksidr, da3_dhetadr);
        }
        public double[,] CalculateRotationMatrix(double[,] e_new, double[,] e_old)
        {
            double[,] Qij = new double[3, 3];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        Qij[i1, i2] += +e_new[k1, i2] * e_old[k1, i1];
                    }
                }
            }

            return Qij;

        }

        private double[,] Transform_FPK_rve_To_FPK_3D(double[,] FPKrve, double[,] Qij, double[,] Qij1)
        {
            var FPK_rve_times_QijT = new double[3, 3]; // ta loops tou transformation tou FPK apo rve se 3D mporoun na veltistopoiithoun. to FPK_rve einai 2x2
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        FPK_rve_times_QijT[i1, i2] += FPKrve[i1, k1] * Qij[i2, k1];
                    }
                }
            }

            var FPK_3D = new double[3, 3];
            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int i2 = 0; i2 < 3; i2++)
                {
                    for (int k1 = 0; k1 < 3; k1++)
                    {
                        FPK_3D[i1, i2] += Qij1[i1, k1] * FPK_rve_times_QijT[k1, i2];
                    }
                }
            }

            return FPK_3D;
        }

        private double[] CalculateDerivativeOfVectorNormalisedArray(double[] g1, double[] dg1_dr)
        {
            double norm = Math.Sqrt(g1[0] * g1[0] + g1[1] * g1[1] + g1[2] * g1[2]);

            double dnorm_dr = (g1[0] * dg1_dr[0] + g1[1] * dg1_dr[1] + g1[2] * dg1_dr[2]) * ((double)1 / norm);

            double[] firstTerm1 = new double[] { dg1_dr[0] * ((double)1 / norm), dg1_dr[1] * ((double)1 / norm), dg1_dr[2] * ((double)1 / norm) }; // 5.26


            double coeff = -dnorm_dr / (Math.Pow(norm, 2));
            double[] secondTerm1 = new double[] { g1[0] * coeff, g1[1] * coeff, g1[2] * coeff };

            return new double[] { firstTerm1[0] + secondTerm1[0], firstTerm1[1] + secondTerm1[1], firstTerm1[2] + secondTerm1[2] };
        }

        #endregion


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





            double[,] StiffnessDevelop_v2 =  new double[_controlPoints.Length * 3, _controlPoints.Length * 3];

                double[] dFPK2D_coefs_dr_vec_prealloc = new double[4];

                var da3_drds_prealloc = new double[3, 3][];
                var da3tilde_drds_prealloc = new double[3, 3][];
                var dF2D_coefs_dr_vec_prealloc = new double[4];
                for (int j = 0; j < gaussPoints.Length; j++)
                {
                    var thicknessGPoints = thicknessIntegrationPoints[gaussPoints[j]];
                    var materialpoint = materialsAtThicknessGP[gaussPoints[j]][thicknessGPoints[0]];
                    //var transformations = new ShellElasticMaterial2DtransformationbDefGrad() { YoungModulus = materialpoint.YoungModulus, PoissonRatio = materialpoint.PoissonRatio };

                    ElementStiffnesses.gpNumber = j;


                    #region current configuration
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


                    var a11 = surfaceBasisVectorDerivative1;
                    var a22 = surfaceBasisVectorDerivative2;
                    var a12 = surfaceBasisVectorDerivative12;
                    var a1 = surfaceBasisVector1;
                    var a2 = surfaceBasisVector2;
                    var a3 = surfaceBasisVector3; // einai to mono pou einai normalised
                    double[] a3_tilde = new double[] { a3[0] * (J1), a3[1] * (J1), a3[2] * (J1) };

                    #endregion

                    (double[] da3tilde_dksi, double[] da3tilde_dheta, double da3norm_dksi, double da3norm_dheta, double[] da3_dksi, double[] da3_dheta) =
                        Calculate_da3tilde_dksi_524_525_526_b(a1, a2, a11, a22, a12, a3, J1);

                    #region original Config
                    var originalControlPoints = _controlPoints.ToArray();
                    var originalHessianMatrix = CalculateHessian(originalControlPoints, _shapeFunctions, j);
                    double[,] jacobian_init = new double[3, 3];
                    CalculateJacobian(originalControlPoints, _shapeFunctions, j, jacobian_init);
                    var a1_init = CalculateSurfaceBasisVector(jacobian_init, 0);
                    var a2_init = CalculateSurfaceBasisVector(jacobian_init, 1);
                    var a3_init = new double[3]
                    {
                    a1_init[1] * a2_init[2] - a1_init[2] * a2_init[1],
                    a1_init[2] * a2_init[0] - a1_init[0] * a2_init[2],
                    a1_init[0] * a2_init[1] - a1_init[1] * a2_init[0]
                    };

                    var J1_init = Math.Sqrt(a3_init[0] * a3_init[0] +
                                            a3_init[1] * a3_init[1] +
                                            a3_init[2] * a3_init[2]);

                    a3_init[0] /= J1_init;
                    a3_init[1] /= J1_init;
                    a3_init[2] /= J1_init;

                    var a11_init = CalculateSurfaceBasisVector(originalHessianMatrix, 0);
                    var a22_init = CalculateSurfaceBasisVector(originalHessianMatrix, 1);
                    var a12_init = CalculateSurfaceBasisVector(originalHessianMatrix, 2);
                    #endregion

                    (double[] da3tilde_dksi_init, double[] da3tilde_dheta_init, double da3norm_dksi_init, double da3norm_dheta_init, double[] da3_dksi_init, double[] da3_dheta_init) =
                        Calculate_da3tilde_dksi_524_525_526_b(a1_init, a2_init, a11_init, a22_init, a12_init, a3_init, J1_init);

                    #region materials at thickness points. (this will be unnecessar if the material is used correctly)
                    var thicknessPoints = thicknessIntegrationPoints[gaussPoints[j]];
                    //double[][,] Aijkl_3D_ofGPs = new double[thicknessGPoints.Count()][,];
                    //double[][,] FPK_2D_ofGPs = new double[thicknessGPoints.Count()][,];
                    //double[][] FPK_3D_vec_ofGPs = new double[thicknessGPoints.Count()][];
                    double[][] G_1_ofGPs = new double[thicknessGPoints.Count()][];
                    double[][] G_2_ofGPs = new double[thicknessGPoints.Count()][];
                    double[][] G_3_ofGPs = new double[thicknessGPoints.Count()][];
                    //tensorOrder2[] FPK_3D_tensor_ofGPs = new tensorOrder2[thicknessGPoints.Count()];
                    double[][,] FPK_3D_tensor_ofGPs_check = new double[thicknessGPoints.Count()][,];
                    //tensorOrder2[] FPKtensorProjected_ofGPs = new tensorOrder2[thicknessGPoints.Count()];
                    double[][,] Ei_of_Gps = new double[thicknessGPoints.Count()][,];
                    //double[][,] Aijkl_2D_ofGPs = new double[thicknessGPoints.Count()][,];
                    double[][,] ei_of_Gps = new double[thicknessGPoints.Count()][,];

                for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                {
                    var materialDevelop = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    var thicknessPoint = thicknessPoints[i1];
                    //var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                    //var w = thicknessPoint.WeightFactor;
                    var z = thicknessPoint.Zeta;

                    //Vector G1 = a1_init + da3_dksi_init.Scale(z);
                    double[] G1 = new double[] { a1_init[0] + da3_dksi_init[0] * z, a1_init[1] + da3_dksi_init[1] * z, a1_init[2] + da3_dksi_init[2] * z };
                    //Vector G2 = a2_init + da3_dheta_init.Scale(z);
                    double[] G2 = new double[] { a2_init[0] + da3_dheta_init[0] * z, a2_init[1] + da3_dheta_init[1] * z, a2_init[2] + da3_dheta_init[2] * z };

                    //double G1_norm_sqred = G1.DotProduct(G1);
                    //double G2_norm_sqred = G2.DotProduct(G2);
                    //double G3_norm_sqred = a3.DotProduct(a3);

                    //double[] G_1 = new double[3] { G1[0] / G1_norm_sqred, G1[1] / G1_norm_sqred, G1[2] / G1_norm_sqred };
                    //double[] G_2 = new double[3] { G2[0] / G2_norm_sqred, G2[1] / G2_norm_sqred, G2[2] / G2_norm_sqred };
                    //(double[] G_1, double[] G_2) = CalculateContravariants(G1, G2);
                    //double[] G_3 = new double[3] { a3[0] / G3_norm_sqred, a3[1] / G3_norm_sqred, a3[2] / G3_norm_sqred };
                    (double[] G_1, double[] G_2, double[] G_3) = CalculateContravariants(G1, G2, a3_init);

                    G_1_ofGPs[i1] = G_1;
                    G_2_ofGPs[i1] = G_2;
                    G_3_ofGPs[i1] = G_3;



                    //Vector g1 = a1 + da3_dksi.Scale(z);
                    double[] g1 = new double[] { a1[0] + da3_dksi[0] * z, a1[1] + da3_dksi[1] * z, a1[2] + da3_dksi[2] * z };
                    //Vector g2 = a2 + da3_dheta.Scale(z);
                    double[] g2 = new double[] { a2[0] + da3_dheta[0] * z, a2[1] + da3_dheta[1] * z, a2[2] + da3_dheta[2] * z };

                    //double[,] F_3D = new double[3, 3] { { g1[0]*G_1[0]+g2[0]*G_2[0], g1[0]*G_1[1]+g2[0]*G_2[1], g1[0]*G_1[2]+g2[0]*G_2[2] },
                    //                                    { g1[1]*G_1[0]+g2[1]*G_2[0], g1[1]*G_1[1]+g2[1]*G_2[1], g1[1]*G_1[2]+g2[1]*G_2[2] },
                    //                                    { g1[2]*G_1[0]+g2[2]*G_2[0], g1[2]*G_1[1]+g2[2]*G_2[1], g1[2]*G_1[2]+g2[2]*G_2[2] },
                    //};
                    //double[,] tgi = new double[3, 3] { { g1[0], g2[0], a3[0] }, { g1[1], g2[1], a3[1] }, { g1[2], g2[2], a3[2] } };
                    //double[,] G_i = new double[3, 3] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } };
                    //double[,] Gi = new double[3, 3] { { G1[0], G2[0], a3_init[0] }, { G1[1], G2[1], a3_init[1] }, { G1[2], G2[2], a3_init[2] } };

                    double[,] ei = new double[3, 3] { { g1[0], g2[0], a3[0] }, { g1[1], g2[1], a3[1] }, { g1[2], g2[2], a3[2] } };
                    double[,] Ei = new double[3, 3] { { G1[0], G2[0], a3_init[0] }, { G1[1], G2[1], a3_init[1] }, { G1[2], G2[2], a3_init[2] } };

                    NormaliseAndOrthogonaliseBasis(ei);
                    NormaliseAndOrthogonaliseBasis(Ei);
                    Ei_of_Gps[i1] = Ei;
                    ei_of_Gps[i1] = ei;

                    //(Aijkl_3D_ofGPs[i1], FPK_3D_vec_ofGPs[i1], FPK_2D_ofGPs[i1], Ei_of_Gps[i1], Aijkl_2D_ofGPs[i1], ei_of_Gps[i1],
                    //    _, /*GL3D, SPKMat3D,*/_, _, _, _ /* ch01_GL_3D, ch01_SPKMat_3D*/, _, _, _, _, _, _, _, _) = transformations.CalculateTransformationsV2(g1, g2, a3, G1, G2, a3_init, G_1, G_2, G_3);


                    #region for higher order terms
                    //var FPK_2D = FPK_2D_ofGPs[i1];
                    var FPK_2D = new double[2, 2] { { materialDevelop.Stresses[0], materialDevelop.Stresses[2] }, { materialDevelop.Stresses[3], materialDevelop.Stresses[1] } };
                    //var FPK3Dcoeffs = new double[3, 3]
                    //        {
                    //                        {FPK_2D[0,0] ,FPK_2D[0,1],0 },
                    //                        {FPK_2D[1,0] ,FPK_2D[1,1],0  },
                    //                        {0,0,0 },
                    //        };
                    //tensorOrder2 FPK_3D_tensor = new tensorOrder2(FPK3Dcoeffs, Ei_of_Gps[i1], ei_of_Gps[i1]);
                    //FPK_3D_tensor_ofGPs[i1] = FPK_3D_tensor.ProjectIn3DCartesianBasis();
                    FPK_3D_tensor_ofGPs_check[i1] = new double[3, 3];
                    for (int j1 = 0; j1 < 3; j1++)
                    {
                        for (int j2 = 0; j2 < 3; j2++)
                        {
                            FPK_3D_tensor_ofGPs_check[i1][j1, j2] += FPK_2D[0, 0] * Ei_of_Gps[i1][j1, 0] * ei_of_Gps[i1][j2, 0] + FPK_2D[0, 1] * Ei_of_Gps[i1][j1, 0] * ei_of_Gps[i1][j2, 1] + FPK_2D[1, 0] * Ei_of_Gps[i1][j1, 1] * ei_of_Gps[i1][j2, 0] + FPK_2D[1, 1] * Ei_of_Gps[i1][j1, 1] * ei_of_Gps[i1][j2, 1];
                        }
                    }


                    #endregion

                    #region higher order for check
                    /*var GL_coeffs = new double[3, 3]
                    {
                            {g1.DotProduct(g1)-G1.DotProduct(G1),g1.DotProduct(g2)-G1.DotProduct(G2),g1.DotProduct(a3)-G1.DotProduct(a3_init) },
                            {g2.DotProduct(g1)-G2.DotProduct(G1),g2.DotProduct(g2)-G2.DotProduct(G2),g2.DotProduct(a3)-G2.DotProduct(a3_init) },
                            {a3.DotProduct(g1)-a3_init.DotProduct(G1),a3.DotProduct(g2)-a3_init.DotProduct(G2),a3.DotProduct(a3)-a3_init.DotProduct(a3_init) },
                    };

                    var corrections = new double[2, 2]
                    {
                        {-da3_dksi.Scale(z).DotProduct(da3_dksi.Scale(z))+da3_dksi_init.Scale(z).DotProduct(da3_dksi_init.Scale(z)),-da3_dksi.Scale(z).DotProduct(da3_dheta.Scale(z))+da3_dksi_init.Scale(z).DotProduct(da3_dheta_init.Scale(z)) },
                        {-da3_dheta.Scale(z).DotProduct(da3_dksi.Scale(z))+da3_dheta_init.Scale(z).DotProduct(da3_dksi_init.Scale(z)),-da3_dheta.Scale(z).DotProduct(da3_dheta.Scale(z))+da3_dheta_init.Scale(z).DotProduct(da3_dheta_init.Scale(z)) }
                    };

                    //for (int j1 = 0; j1 < 2; j1++)
                    //{
                    //    for (int j2 = 0; j2 < 2; j2++)
                    //    {
                    //        GL_coeffs[j1, j2] += corrections[j1, j2];
                    //    }
                    //}
                    tensorOrder2 GLtensor = new tensorOrder2(GL_coeffs);
                    GLtensor = GLtensor.Scale(0.5);

                    GLtensor.ReplaceBasisWithVector(G_1, G_2, G_3, true);
                    GLtensor.ReplaceBasisWithVector(G_1, G_2, G_3, false);

                    var GLtensorProjected = GLtensor.ProjectIn3DCartesianBasis();

                    var materialaux_ = material.Clone();
                    materialaux_.TangentVectorV1 = material.TangentVectorV1;
                    materialaux_.TangentVectorV2 = material.TangentVectorV2;
                    materialaux_.NormalVectorV3 = material.NormalVectorV3;

                    materialaux_.UpdateMaterial(new double[] { GLtensor.coefficients[0, 0], GLtensor.coefficients[1, 1], 2 * GLtensor.coefficients[0, 1] });

                    var stresses = materialaux_.Stresses;

                    var SPK_coeffs = new double[3, 3]
                    {
                            {stresses[0], stresses[2], 0 },
                            {stresses[2], stresses[1], 0 },
                            {0,0,0 },
                    };

                    var SPKtensor = new tensorOrder2(SPK_coeffs);

                    SPKtensor.ReplaceBasisWithVector(G1, G2, a3_init, true);
                    SPKtensor.ReplaceBasisWithVector(G1, G2, a3_init, false);

                    var SPKtensorProjected = SPKtensor.ProjectIn3DCartesianBasis();
                    double[,] ch_F_3D = Calculate3DtensorFrom2Dcorrected_normaliseBothBasesCase(new double[,] { { 1, 0, 0 }, { 0, 1, 0, }, { 0, 0, 1 } }, g1, g2, a3, G_1, G_2, G_3);
                    var defGradTensor = new tensorOrder2(ch_F_3D);

                    var defGradTensorTr = defGradTensor.Transpose();

                    var FPKtensor = SPKtensorProjected.SingleContract(defGradTensorTr);

                    FPKtensorProjected_ofGPs[i1] = FPKtensor.ProjectIn3DCartesianBasis();*/
                    #endregion


                }
                    #endregion

                    #region small loop precalculated values
                    //315 line and 
                    //1442
                    var a3rArray = new a3r[elementControlPoints.Length];
                    //var a1rArray = new Matrix3by3[elementControlPoints.Length];
                    //var a2rArray = new Matrix3by3[elementControlPoints.Length];
                    //var a11rArray = new Matrix3by3[elementControlPoints.Length];
                    //var a22rArray = new Matrix3by3[elementControlPoints.Length];
                    //var a12rArray = new Matrix3by3[elementControlPoints.Length];

                    var da3tilde_drArray = new double[elementControlPoints.Length][][]; //now an index will be necessary
                    var dnorma3_drArray = new double[elementControlPoints.Length][];
                    var da3tilde_dksidrArray = new double[elementControlPoints.Length][][];
                    var da3tilde_dhetadrArray = new double[elementControlPoints.Length][][];
                    var da3norm_dksidrArray = new double[elementControlPoints.Length][];
                    var da3norm_dhetadrArray = new double[elementControlPoints.Length][];
                    var da3_dksidrArray = new double[elementControlPoints.Length][][];
                    var da3_dhetadrArray = new double[elementControlPoints.Length][][];
                    double[,][] dF_3D_dr_vecArray = new double[3 * elementControlPoints.Length, thicknessPoints.Count()][];
                    double[,][] dF_3D_dr_tensor_tr_vecArray = new double[3 * elementControlPoints.Length, thicknessPoints.Count()][];

                    double[,] dF3Dtensor_dr_Tr_check = new double[3, 3];
                    for (int i = 0; i < elementControlPoints.Length; i++)
                    {
                        var a3r = new a3r();
                        var dksi_r = _shapeFunctions.DerivativeValuesKsi[i, j];
                        var dheta_r = _shapeFunctions.DerivativeValuesHeta[i, j];
                        var d2Ksi_dr = _shapeFunctions.SecondDerivativeValuesKsi[i, j];
                        var d2Heta_dr = _shapeFunctions.SecondDerivativeValuesHeta[i, j];
                        var d2KsiHeta_dr = _shapeFunctions.SecondDerivativeValuesKsiHeta[i, j];

                        CalculateA3r(surfaceBasisVector1, surfaceBasisVector2, surfaceBasisVector3, dksi_r,
                            dheta_r, J1, ref a3r);
                        a3rArray[i] = a3r;

                        //a1rArray[i] = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                        //a2rArray[i] = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                        //a11rArray[i] = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                        //a22rArray[i] = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                        //a12rArray[i] = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);

                        //5.24
                        da3tilde_drArray[i] = new double[3][];
                        Calculate_da3tilde_dr_Array(a1, a2, dksi_r, dheta_r, da3tilde_drArray[i]);

                        //5.25
                        dnorma3_drArray[i] = new double[3];
                        a3_tilde = CalculateTerm525_Array(a3, J1, dnorma3_drArray[i], da3tilde_drArray[i]);

                        //5.30 b
                        //(da3tilde_dksidrArray[i], da3tilde_dhetadrArray[i]) = Calculate_da3tilde_dksidr(a1rArray[i], a2rArray[i], a11rArray[i], a22rArray[i], a12rArray[i], a1, a2, a11, a22, a12);
                        (da3tilde_dksidrArray[i], da3tilde_dhetadrArray[i]) = Calculate_da3tilde_dksidrDvelop(a1, a2, a11, a22, a12,
                            dksi_r, dheta_r, d2Ksi_dr, d2Heta_dr, d2KsiHeta_dr);


                        //5.31 b
                        (da3norm_dksidrArray[i], da3norm_dhetadrArray[i]) = Calculate_da3norm_dksidr_Develop(da3tilde_dksidrArray[i], da3tilde_dhetadrArray[i],
                            a3_tilde, da3tilde_dksi, da3tilde_dheta, da3tilde_drArray[i], J1);

                        //5.32 b
                        (da3_dksidrArray[i], da3_dhetadrArray[i]) = Calculate_da3_dksidrDevelop(da3tilde_dksidrArray[i], da3tilde_dhetadrArray[i], da3tilde_dksi, da3tilde_dheta,
                            dnorma3_drArray[i], a3_tilde, da3norm_dksidrArray[i], da3norm_dhetadrArray[i], da3norm_dksi, da3norm_dheta, J1, da3tilde_drArray[i]);

                        for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                        {
                            var thicknessPoint = thicknessPoints[i1];
                            //var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                            //var w = thicknessPoint.WeightFactor;
                            var z = thicknessPoint.Zeta;

                            for (int r1 = 0; r1 < 3; r1++)
                            {
                                //(31)
                                //Vector dg1_dr = a1rArray[i].GetColumn(r1) + da3_dksidrArray[i][r1] * z;
                                double[] dg1_dr = new double[3] { da3_dksidrArray[i][r1][0] * z, da3_dksidrArray[i][r1][1] * z, da3_dksidrArray[i][r1][2] * z };
                                dg1_dr[r1] += dksi_r;

                                //Vector dg2_dr = a2rArray[i].GetColumn(r1) + da3_dhetadrArray[i][r1] * z;
                                double[] dg2_dr = new double[3] { da3_dhetadrArray[i][r1][0] * z, da3_dhetadrArray[i][r1][1] * z, da3_dhetadrArray[i][r1][2] * z };
                                dg2_dr[r1] += dheta_r;




                                var G_1 = G_1_ofGPs[i1];
                                var G_2 = G_2_ofGPs[i1];
                                var G_3 = G_3_ofGPs[i1];

                                //(39)
                                double[,] dF_3D_dr = new double[3, 3] { { dg1_dr[0]*G_1[0]+dg2_dr[0]*G_2[0], dg1_dr[0]*G_1[1]+dg2_dr[0]*G_2[1], dg1_dr[0]*G_1[2]+dg2_dr[0]*G_2[2] },
                                                                 { dg1_dr[1]*G_1[0]+dg2_dr[1]*G_2[0], dg1_dr[1]*G_1[1]+dg2_dr[1]*G_2[1], dg1_dr[1]*G_1[2]+dg2_dr[1]*G_2[2] },
                                                                 { dg1_dr[2]*G_1[0]+dg2_dr[2]*G_2[0], dg1_dr[2]*G_1[1]+dg2_dr[2]*G_2[1], dg1_dr[2]*G_1[2]+dg2_dr[2]*G_2[2] }, };

                                //double[] dF_3D_dr_vec =
                                dF_3D_dr_vecArray[3 * i + r1, i1] = new double[] { dF_3D_dr[0, 0], dF_3D_dr[1, 1], dF_3D_dr[2, 2], dF_3D_dr[0, 1], dF_3D_dr[1, 2], dF_3D_dr[2, 0], dF_3D_dr[0, 2], dF_3D_dr[1, 0], dF_3D_dr[2, 1] };

                                var da3_dr = new double[3];
                                if (r1 == 0) { da3_dr[0] = a3r.a3r00; da3_dr[1] = a3r.a3r10; da3_dr[2] = a3r.a3r20; }
                                else if (r1 == 1) { da3_dr[0] = a3r.a3r01; da3_dr[1] = a3r.a3r11; da3_dr[2] = a3r.a3r21; }
                                else if (r1 == 2) { da3_dr[0] = a3r.a3r02; da3_dr[1] = a3r.a3r12; da3_dr[2] = a3r.a3r22; }//.

                                //tensorOrder2 dF3Dtensor_dr = new tensorOrder2()
                                //{
                                //    basis1 = new double[,] { { dg1_dr[0], dg2_dr[0], da3_dr[0] }, { dg1_dr[1], dg2_dr[1], da3_dr[1] }, { dg1_dr[2], dg2_dr[2], da3_dr[2] } },
                                //    basis2 = new double[,] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } },
                                //    coefficients = new double[,] { { 1, 0, 0 }, { 0, 1, 0, }, { 0, 0, 1 } }
                                //};

                                //tensorOrder2 dF3Dtensor_dr_Tr = dF3Dtensor_dr.Transpose().ProjectIn3DCartesianBasis();

                                Array.Clear(dF3Dtensor_dr_Tr_check, 0, 9);
                                for (int j1 = 0; j1 < 3; j1++)
                                {
                                    for (int j2 = 0; j2 < 3; j2++)
                                    {
                                        dF3Dtensor_dr_Tr_check[j2, j1] += dg1_dr[j1] * G_1[j2] + dg2_dr[j1] * G_2[j2] + da3_dr[j1] * G_3[j2];
                                    }
                                }


                                dF_3D_dr_tensor_tr_vecArray[3 * i + r1, i1] = new double[] { dF3Dtensor_dr_Tr_check[0, 0], dF3Dtensor_dr_Tr_check[1, 1], dF3Dtensor_dr_Tr_check[2, 2], dF3Dtensor_dr_Tr_check[0, 1], dF3Dtensor_dr_Tr_check[1, 2], dF3Dtensor_dr_Tr_check[2, 0], dF3Dtensor_dr_Tr_check[0, 2], dF3Dtensor_dr_Tr_check[1, 0], dF3Dtensor_dr_Tr_check[2, 1] };

                            }
                            //  (31) 



                        }
                    }
                    #endregion

                    #region actual loop, total calculations
                    var dFPK_tensor_dr_check = new double[3, 3];
                    var dF3Dtensor_drds_Tr_check = new double[3, 3];
                    for (int i = 0; i < elementControlPoints.Length; i++)
                    {
                        var dksi_r = _shapeFunctions.DerivativeValuesKsi[i, j];
                        var dheta_r = _shapeFunctions.DerivativeValuesHeta[i, j];
                        var d2Ksi_dr = _shapeFunctions.SecondDerivativeValuesKsi[i, j];
                        var d2Heta_dr = _shapeFunctions.SecondDerivativeValuesHeta[i, j];
                        var d2KsiHeta_dr = _shapeFunctions.SecondDerivativeValuesKsiHeta[i, j];
                        var a3r = a3rArray[i];

                        //var a1r = a1rArray[i];
                        //var a2r = a2rArray[i];
                        //var a11r = a11rArray[i];
                        //var a22r = a22rArray[i];
                        //var a12r = a12rArray[i];

                        //var da3tilde_dr = da3tilde_drArray[i];
                        //var dnorma3_dr = dnorma3_drArray[i];
                        var da3tilde_dksidr = da3tilde_dksidrArray[i];
                        var da3tilde_dhetadr = da3tilde_dhetadrArray[i];
                        var da3norm_dksidr = da3norm_dksidrArray[i];
                        var da3norm_dhetadr = da3norm_dhetadrArray[i];
                        var da3_dksidr = da3_dksidrArray[i];
                        var da3_dhetadr = da3_dhetadrArray[i];

                        //if (ElementStiffnesses.performCalculations)
                        //{
                        //    var a1r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesKsi[i, j]);
                        //    var a2r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsDerivativeValuesHeta[i, j]);
                        //    var a11r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsi[i, j]);
                        //    var a22r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueHeta[i, j]);
                        //    var a12r = Matrix3by3.CreateIdentity().Scale(nurbs.NurbsSecondDerivativeValueKsiHeta[i, j]);
                        //    if (ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck)
                        //    {
                        //        for (int i1 = 0; i1 < 3; i1++)
                        //        {
                        //            ElementStiffnesses.ProccessVariable(1, a11r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                        //            ElementStiffnesses.ProccessVariable(2, a22r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);
                        //            ElementStiffnesses.ProccessVariable(3, a12r.GetColumn(i1).CopyToArray(), true, 3 * i + i1);

                        //        }
                        //    }
                        //}

                        for (int k = i; k < elementControlPoints.Length; k++)
                        {
                            var d2Ksi_ds = _shapeFunctions.SecondDerivativeValuesKsi[k, j];
                            var d2Heta_ds = _shapeFunctions.SecondDerivativeValuesHeta[k, j];
                            var d2KsiHeta_ds = _shapeFunctions.SecondDerivativeValuesKsiHeta[k, j];
                            var dksi_s = _shapeFunctions.DerivativeValuesKsi[k, j];
                            var dheta_s = _shapeFunctions.DerivativeValuesHeta[k, j];

                            //var a3s = a3rArray[k];

                            //var a1s = a1rArray[k];
                            //var a2s = a2rArray[k];
                            //var a11s = a11rArray[k];
                            //var a22s = a22rArray[k];
                            //var a12s = a12rArray[k];

                            //var da3tilde_ds = da3tilde_drArray[k];
                            //var dnorma3_ds = dnorma3_drArray[k];
                            var da3tilde_dksids = da3tilde_dksidrArray[k];
                            var da3tilde_dhetads = da3tilde_dhetadrArray[k];
                            var da3norm_dksids = da3norm_dksidrArray[k];
                            var da3norm_dhetads = da3norm_dhetadrArray[k];
                            //var da3_dksids = da3_dksidrArray[k];
                            //var da3_dhetads = da3_dhetadrArray[k];

                            //// =apo prohgoiumena 
                            //(a3rs a3rsAlternative, var da3tilde_drds, _, _,
                            //_, _, double[,] dnorma3_drds, _, var da3_drds) =
                            //Calculate_a3rs(Vector.CreateFromArray(surfaceBasisVector1), Vector.CreateFromArray(surfaceBasisVector2),
                            //    Vector.CreateFromArray(surfaceBasisVector3), J1, dksi_r, dksi_s, dheta_r, dheta_s);

                            //// develop
                            //(Vector[,] da3tilde_dksidrds, Vector[,] da3tilde_dhetadrds) = Calculate_530_c(a1, a2,
                            //    a11r, a12r, a22r, a1r, a2r, a11s, a12s, a22s, a1s, a2s, a11, a12, a22);

                            //(double[,] da3norm_dksidrds, double[,] da3norm_dhetadrds) = Calculate_531_c(J1, da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksids, da3tilde_dhetads, a3_tilde,
                            //    da3tilde_dksidrds, da3tilde_dhetadrds, a3r, a3s, da3_dksi, da3_dheta, a3, da3tilde_drds, da3tilde_dr, da3tilde_ds, dnorma3_dr, dnorma3_ds, da3tilde_dksi, da3tilde_dheta);

                            //(Vector[,] da3_dksidrds, Vector[,] da3_dhetadrds) = Calculate_532_c(J1, da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksids, da3tilde_dhetads, a3_tilde,
                            //    da3tilde_dksidrds, da3tilde_dhetadrds, a3r, a3s, da3_dksi, da3_dheta, a3, da3tilde_drds, da3tilde_dr, da3tilde_ds, dnorma3_dr, dnorma3_ds, da3tilde_dksi, da3tilde_dheta,
                            //    da3norm_dksidrds, da3norm_dhetadrds, dnorma3_drds, da3norm_dksi, da3norm_dheta, da3norm_dksids, da3norm_dhetads, da3norm_dhetadr, da3norm_dksidr);


                            (double[,][] da3_dksidrds, double[,][] da3_dhetadrds) = Calculate_dvariable_drds_second_derivatives(surfaceBasisVector1, surfaceBasisVector2,
                                surfaceBasisVector3, J1, dksi_r, dksi_s, dheta_r, dheta_s, da3tilde_dksidr, da3tilde_dhetadr, da3tilde_dksids, da3tilde_dhetads, a3_tilde,
                                  da3tilde_dksi, da3tilde_dheta,
                                da3norm_dksi, da3norm_dheta, da3norm_dksids, da3norm_dhetads, da3norm_dhetadr, da3norm_dksidr, dksi_r, dheta_r, d2Ksi_dr, d2Heta_dr, d2KsiHeta_dr,
                                d2Ksi_ds, d2Heta_ds, d2KsiHeta_ds, dksi_s, dheta_s, da3_drds_prealloc, da3tilde_drds_prealloc);



                            if (ElementStiffnesses.performCalculations)
                            {
                                if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveStiffnessMatrixState)
                                {
                                    ElementStiffnesses.saveOriginalState = true;

                                    if (i == 0) // ennoume thn paragwgo ws pros r gia k=0 kai gai to x dof
                                    {


                                        // stathero to r=0 allazei to s afou edw vlepoume oti to exoume exarthsei apo to k
                                        //ElementStiffnesses.ProccessVariable(21, da3tilde_dksidrds[0, 0].CopyToArray(), true, 3 * k + 0);
                                        //ElementStiffnesses.ProccessVariable(21, da3tilde_dksidrds[0, 1].CopyToArray(), true, 3 * k + 1);
                                        //ElementStiffnesses.ProccessVariable(21, da3tilde_dksidrds[0, 2].CopyToArray(), true, 3 * k + 2);
                                        //ElementStiffnesses.ProccessVariable(21, da3tilde_dksidr[0].CopyToArray(), false);

                                        //ElementStiffnesses.ProccessVariable(22, da3tilde_dhetadrds[0, 0].CopyToArray(), true, 3 * k + 0);
                                        //ElementStiffnesses.ProccessVariable(22, da3tilde_dhetadrds[0, 1].CopyToArray(), true, 3 * k + 1);
                                        //ElementStiffnesses.ProccessVariable(22, da3tilde_dhetadrds[0, 2].CopyToArray(), true, 3 * k + 2);
                                        //ElementStiffnesses.ProccessVariable(22, da3tilde_dhetadr[0].CopyToArray(), false);

                                        //ElementStiffnesses.ProccessVariable(23, new double[] { da3norm_dksidrds[0, 0] }, true, 3 * k + 0);
                                        //ElementStiffnesses.ProccessVariable(23, new double[] { da3norm_dksidrds[0, 1] }, true, 3 * k + 1);
                                        //ElementStiffnesses.ProccessVariable(23, new double[] { da3norm_dksidrds[0, 2] }, true, 3 * k + 2);
                                        //ElementStiffnesses.ProccessVariable(23, new double[] { da3norm_dksidr[0] }, false);

                                        //ElementStiffnesses.ProccessVariable(24, new double[] { da3norm_dhetadrds[0, 1] }, true, 3 * k + 1);
                                        //ElementStiffnesses.ProccessVariable(24, new double[] { da3norm_dhetadrds[0, 2] }, true, 3 * k + 2);
                                        //ElementStiffnesses.ProccessVariable(24, new double[] { da3norm_dhetadrds[0, 0] }, true, 3 * k + 0);
                                        //ElementStiffnesses.ProccessVariable(24, new double[] { da3norm_dhetadr[0] }, false);

                                        ElementStiffnesses.ProccessVariable(25, da3_dksidrds[0, 0], true, 3 * k + 0);
                                        ElementStiffnesses.ProccessVariable(25, da3_dksidrds[0, 1], true, 3 * k + 1);
                                        ElementStiffnesses.ProccessVariable(25, da3_dksidrds[0, 2], true, 3 * k + 2);
                                        ElementStiffnesses.ProccessVariable(25, da3_dksidr[0], false);

                                        ElementStiffnesses.ProccessVariable(26, da3_dhetadrds[0, 0], true, 3 * k + 0);
                                        ElementStiffnesses.ProccessVariable(26, da3_dhetadrds[0, 1], true, 3 * k + 1);
                                        ElementStiffnesses.ProccessVariable(26, da3_dhetadrds[0, 2], true, 3 * k + 2);
                                        ElementStiffnesses.ProccessVariable(26, da3_dhetadr[0], false);

                                    }


                                    ElementStiffnesses.saveOriginalState = false;

                                }

                                if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveVariationStates)
                                {
                                    if (i == 0) // ennoume thn paragwgo ws pros r gia k=0 kai gai to x dof
                                    {
                                        //ElementStiffnesses.ProccessVariable(9, new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 }, false);

                                        ElementStiffnesses.ProccessVariable(21, da3tilde_dksidr[0], false);

                                        ElementStiffnesses.ProccessVariable(22, da3tilde_dhetadr[0], false);

                                        ElementStiffnesses.ProccessVariable(23, new double[] { da3norm_dksidr[0] }, false);

                                        ElementStiffnesses.ProccessVariable(24, new double[] { da3norm_dhetadr[0] }, false);

                                        ElementStiffnesses.ProccessVariable(25, da3_dksidr[0], false);

                                        ElementStiffnesses.ProccessVariable(26, da3_dhetadr[0], false);
                                    }

                                    ElementStiffnesses.saveOriginalState = false;

                                }
                            }

                            for (int i1 = 0; i1 < thicknessPoints.Count; i1++)
                            {
                                var materialDevelop = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                                var thicknessPoint = thicknessPoints[i1];
                                //var material = materialsAtThicknessGP[gaussPoints[j]][thicknessPoints[i1]];
                                var w = thicknessPoint.WeightFactor;
                                var z = thicknessPoint.Zeta;

                                var G_1 = G_1_ofGPs[i1];
                                var G_2 = G_2_ofGPs[i1];
                                var G_3 = G_3_ofGPs[i1];

                                //var FPK_2D = FPK_2D_ofGPs[i1];
                                var FPK_2D = new double[2, 2] { { materialDevelop.Stresses[0], materialDevelop.Stresses[2] }, { materialDevelop.Stresses[3], materialDevelop.Stresses[1] } };

                                var Ei = Ei_of_Gps[i1];
                                //var Aijkl_2D = Aijkl_2D_ofGPs[i1];
                                var ei = ei_of_Gps[i1];
                                double[] g1 = new double[] { a1[0] + da3_dksi[0] * z, a1[1] + da3_dksi[1] * z, a1[2] + da3_dksi[2] * z };
                                double[] g2 = new double[] { a2[0] + da3_dheta[0] * z, a2[1] + da3_dheta[1] * z, a2[2] + da3_dheta[2] * z };

                                double norm_g2 = Math.Sqrt(g2[0] * g2[0] + g2[1] * g2[1] + g2[2] * g2[2]);
                                double[] e2_init = new double[] { g2[0] * ((double)1 / norm_g2), g2[1] * ((double)1 / norm_g2), g2[2] * ((double)1 / norm_g2) }; // einai to g2 unit 

                                #region higher order terms FPK:F3D_drds

                                //var FPK3Dcoeffs = new double[3, 3]
                                //{
                                //                {FPK_2D[0,0] ,FPK_2D[0,1],0 },
                                //                {FPK_2D[1,0] ,FPK_2D[1,1],0  },
                                //                {0,0,0 },
                                //};
                                //tensorOrder2 FPK_3D_tensor = new tensorOrder2(FPK3Dcoeffs, Ei, ei);
                                //FPK_3D_tensor = FPK_3D_tensor.ProjectIn3DCartesianBasis();

                                //tensorOrder2 FPK_3D_tensor = FPK_3D_tensor_ofGPs[i1];
                                var FPK_3D_tensor_check = FPK_3D_tensor_ofGPs_check[i1];
                                //tensorOrder2 FPKtensorProjected = FPKtensorProjected_ofGPs[i1];
                                #endregion


                                dFrvedri dFrve_dr = new dFrvedri();

                                for (int r1 = 0; r1 < 3; r1++)
                                {
                                    //var dF_3D_dr1 = dF_3D_dr_vecArray[3 * i + r1, i1];
                                    //(31)
                                    //Vector dg1_dr = a1rArray[i].GetColumn(r1) + da3_dksidrArray[i][r1] * z;
                                    double[] dg1_dr = new double[3] { da3_dksidrArray[i][r1][0] * z, da3_dksidrArray[i][r1][1] * z, da3_dksidrArray[i][r1][2] * z };
                                    dg1_dr[r1] += dksi_r;

                                    //Vector dg2_dr = a2rArray[i].GetColumn(r1) + da3_dhetadrArray[i][r1] * z;
                                    double[] dg2_dr = new double[3] { da3_dhetadrArray[i][r1][0] * z, da3_dhetadrArray[i][r1][1] * z, da3_dhetadrArray[i][r1][2] * z };
                                    dg2_dr[r1] += dheta_r;

                                    var da3_dr = new double[3];
                                    if (r1 == 0) { da3_dr[0] = a3r.a3r00; da3_dr[1] = a3r.a3r10; da3_dr[2] = a3r.a3r20; }
                                    else if (r1 == 1) { da3_dr[0] = a3r.a3r01; da3_dr[1] = a3r.a3r11; da3_dr[2] = a3r.a3r21; }
                                    else if (r1 == 2) { da3_dr[0] = a3r.a3r02; da3_dr[1] = a3r.a3r12; da3_dr[2] = a3r.a3r22; }//.

                                    double[] de1_dr = CalculateDerivativeOfVectorNormalisedArray(g1, dg1_dr);
                                    double[] de2_init_dr = CalculateDerivativeOfVectorNormalisedArray(g2, dg2_dr); // einai to dg2unit_dr
                                    double[] de2_dr = CalculateDerivativeOfOrthogonalisedNormalisedCovariantBasisArray(ei, de1_dr, e2_init, de2_init_dr);




                                    for (int s1 = 0; s1 < 3; s1++)
                                    {
                                        // linear stiffness
                                        var dF_3D_ds1 = dF_3D_dr_vecArray[3 * k + s1, i1];
                                        var dF_3D_tr_ds1 = dF_3D_dr_tensor_tr_vecArray[3 * k + s1, i1];



                                        var dg1_drds = new double[] { da3_dksidrds[r1, s1][0] * z, da3_dksidrds[r1, s1][1] * z, da3_dksidrds[r1, s1][2] * z }; //Scale(z);
                                        var dg2_drds = new double[] { da3_dhetadrds[r1, s1][0] * z, da3_dhetadrds[r1, s1][1] * z, da3_dhetadrds[r1, s1][2] * z };//.





                                        //(_,_/*double[,] Pcontributions_term_1, double[] dF2D_coefs_dr_vec*/, double[] dFPK2D_coefs_dr_vec, double[,] dFPK2D_coefs_dr) = Calculate_FPK_variation_term1(
                                        //materialDevelop.ConstitutiveMatrix.CopytoArray2D(), dg1_dr, dg2_dr, de1_dr, de2_dr, G_1, G_2, a3r, Ei, ei, g1, g2, a3);// to do anti gia copy to array i methodos na mporei na doulevei me IMatrixView
                                        Calculate_FPK_variation_term1_Short(
                                            materialDevelop.ConstitutiveMatrix, dg1_dr, dg2_dr, de1_dr, de2_dr, G_1, G_2, a3r, Ei, ei, g1, g2, dF2D_coefs_dr_vec_prealloc, dFPK2D_coefs_dr_vec_prealloc, ref dFrve_dr);//. to do anti gia copy to array i methodos na mporei na doulevei me IMatrixView
                                        //var p_contrib_coeffs_ = new double[3, 3]
                                        //{
                                        //        {FPK_2D[0,0] ,FPK_2D[0,1],0 },
                                        //        {FPK_2D[1,0] ,FPK_2D[1,1],0  },
                                        //        {0,0,0 },
                                        //};
                                        //tensorOrder2 P_contrib_tensor = new tensorOrder2(p_contrib_coeffs_, Ei, ei);
                                        //P_contrib_tensor.ReplaceBasisWithVector(de1_dr, de2_dr, da3_dr, false);
                                        //P_contrib_tensor = P_contrib_tensor.ProjectIn3DCartesianBasis();


                                        //var p_contrib_term_1_coeffs_ = new double[3, 3]
                                        //{
                                        //        {dFPK2D_coefs_dr[0,0] ,dFPK2D_coefs_dr[0,1],0 },
                                        //        {dFPK2D_coefs_dr[1,0] ,dFPK2D_coefs_dr[1,1],0 },
                                        //        {0,0,0 },
                                        //};
                                        //tensorOrder2 P_contrib_term1_tensor = new tensorOrder2(p_contrib_term_1_coeffs_, Ei, ei);
                                        //P_contrib_term1_tensor = P_contrib_term1_tensor.ProjectIn3DCartesianBasis();

                                        //tensorOrder2 dFPK_tensor_dr = P_contrib_tensor.AddTensor(P_contrib_term1_tensor);
                                        //dFPK_tensor_dr.CorrectLimitsToZero();
                                        //double[] dFPK_tensor_dr_vec = { dFPK_tensor_dr.coefficients[0, 0], dFPK_tensor_dr.coefficients[1, 1], dFPK_tensor_dr.coefficients[2, 2], dFPK_tensor_dr.coefficients[0, 1], dFPK_tensor_dr.coefficients[1, 2], dFPK_tensor_dr.coefficients[2, 0], dFPK_tensor_dr.coefficients[0, 2], dFPK_tensor_dr.coefficients[1, 0], dFPK_tensor_dr.coefficients[2, 1] };// 

                                        Array.Clear(dFPK_tensor_dr_check, 0, 9);
                                        for (int j1 = 0; j1 < 3; j1++)
                                        {
                                            for (int j2 = 0; j2 < 3; j2++)
                                            {
                                                dFPK_tensor_dr_check[j1, j2] += FPK_2D[0, 0] * Ei[j1, 0] * de1_dr[j2] + FPK_2D[0, 1] * Ei[j1, 0] * de2_dr[j2] + FPK_2D[1, 0] * Ei[j1, 1] * de1_dr[j2] + FPK_2D[1, 1] * Ei[j1, 1] * de2_dr[j2] +
                                                    dFrve_dr.r00 * Ei[j1, 0] * ei[j2, 0] + dFrve_dr.r01 * Ei[j1, 0] * ei[j2, 1] + dFrve_dr.r10 * Ei[j1, 1] * ei[j2, 0] + dFrve_dr.r11 * Ei[j1, 1] * ei[j2, 1];
                                            }
                                        }

                                        double[] dFPK_tensor_dr_vec = new double[] { dFPK_tensor_dr_check[0, 0], dFPK_tensor_dr_check[1, 1], dFPK_tensor_dr_check[2, 2], dFPK_tensor_dr_check[0, 1], dFPK_tensor_dr_check[1, 2], dFPK_tensor_dr_check[2, 0], dFPK_tensor_dr_check[0, 2], dFPK_tensor_dr_check[1, 0], dFPK_tensor_dr_check[2, 1] };// 


                                        //var test01 = 0.01;
                                        //for (int j1 = 0; j1 < 3; j1++)
                                        //{
                                        //    for (int j2 = 0; j2 < 3; j2++)
                                        //    {
                                        //       test01 += Ei[j1, 0] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 0] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 0] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 0] * Ei[j1, 1] * Ei[j1, 1] * Ei[j1, 0] * Ei[j1, 1] * Ei[j1, 1];
                                        //            //dFPK2D_coefs_dr_prealloc[0, 0] * Ei[j1, 0] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[0, 1] * Ei[j1, 0] * ei[j2, 1] + dFPK2D_coefs_dr_prealloc[1, 0] * Ei[j1, 1] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[1, 1] * Ei[j1, 1] * ei[j2, 1];
                                        //    }
                                        //}

                                        //for (int j1 = 0; j1 < 3; j1++)
                                        //{
                                        //    for (int j2 = 0; j2 < 3; j2++)
                                        //    {
                                        //        test01 += de1_dr[j2]* de2_dr[j2]* de1_dr[j2] * de2_dr[j2]* de1_dr[j2] * de2_dr[j2] * de1_dr[j2] * de2_dr[j2]* de1_dr[j2] * de2_dr[j2] * de1_dr[j2] * de2_dr[j2] * de2_dr[j2] * de1_dr[j2] * de2_dr[j2] * de1_dr[j2];
                                        //        //dFPK2D_coefs_dr_prealloc[0, 0] * Ei[j1, 0] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[0, 1] * Ei[j1, 0] * ei[j2, 1] + dFPK2D_coefs_dr_prealloc[1, 0] * Ei[j1, 1] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[1, 1] * Ei[j1, 1] * ei[j2, 1];
                                        //    }
                                        //}

                                        //for (int j1 = 0; j1 < 3; j1++)
                                        //{
                                        //    for (int j2 = 0; j2 < 3; j2++)
                                        //    {
                                        //        test01 += dFrve_dr[j1, 0] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 0] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 0] *
                                        //            dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 0] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1] * dFrve_dr[j1, 0] * dFrve_dr[j1, 1] * dFrve_dr[j1, 1];
                                        //        //dFPK2D_coefs_dr_prealloc[0, 0] * Ei[j1, 0] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[0, 1] * Ei[j1, 0] * ei[j2, 1] + dFPK2D_coefs_dr_prealloc[1, 0] * Ei[j1, 1] * ei[j2, 0] + dFPK2D_coefs_dr_prealloc[1, 1] * Ei[j1, 1] * ei[j2, 1];
                                        //    }
                                        //}


                                        //#region Higher order terms FPK:Ft_drds
                                        //tensorOrder2 dF3Dtensor_drds = new tensorOrder2()
                                        //{
                                        //    basis1 = new double[,] { { dg1_drds[0], dg2_drds[0], da3_drds[r1, s1][0] }, { dg1_drds[1], dg2_drds[1], da3_drds[r1, s1][1] }, { dg1_drds[2], dg2_drds[2], da3_drds[r1, s1][2] } },
                                        //    basis2 = new double[,] { { G_1[0], G_2[0], G_3[0] }, { G_1[1], G_2[1], G_3[1] }, { G_1[2], G_2[2], G_3[2] } },
                                        //    coefficients = new double[,] { { 1, 0, 0 }, { 0, 1, 0, }, { 0, 0, 1 } }
                                        //};

                                        //tensorOrder2 dF3Dtensor_drds_Tr = dF3Dtensor_drds.Transpose().ProjectIn3DCartesianBasis();


                                        Array.Clear(dF3Dtensor_drds_Tr_check, 0, 9);
                                        for (int j1 = 0; j1 < 3; j1++)
                                        {
                                            for (int j2 = 0; j2 < 3; j2++)
                                            {
                                                dF3Dtensor_drds_Tr_check[j2, j1] += dg1_drds[j1] * G_1[j2] + dg2_drds[j1] * G_2[j2] + da3_drds_prealloc[r1, s1][j1] * G_3[j2];
                                            }
                                        }
                                        #endregion



                                        #region teliko assembly tou v1
                                        for (int i2 = 0; i2 < 9; i2++)
                                        {
                                            //StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += FPK_3D_vec[i2] * dF_3D_drds_vecArray[i2] * w * wFactor;
                                        }

                                    //StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += FPK_3D_tensor.doubleContract(dF3Dtensor_drds_Tr) * w * wFactor;//this is going to be replaced
                                    for (int j1 = 0; j1 < 3; j1++)
                                    {
                                        for (int j2 = 0; j2 < 3; j2++)
                                        {
                                            StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += FPK_3D_tensor_check[j1, j2] * dF3Dtensor_drds_Tr_check[j1, j2] * w * wFactor;
                                         }

                                    }

                                        //stiffnessDevelop_v2_part1[3 * i + r1, 3 * k + s1] += FPK_3D_tensor.doubleContract(dF3Dtensor_drds_Tr) * w * wFactor;
                                        //for (int p1 = 0; p1 < 9; p1++)
                                        //{
                                        //    StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += Pcontributions_dr_vec[p1] * dF_3D_ds1[p1] * w * wFactor;
                                        //    StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += Pcontributions_term_1_dr_vec[p1] * dF_3D_ds1[p1] * w * wFactor;

                                        //}

                                        //double[] dF_3D_ds1_tr = new double[9]
                                        //    {dF_3D_ds1[0],dF_3D_ds1[1],dF_3D_ds1[2],dF_3D_ds1[7],dF_3D_ds1[8],dF_3D_ds1[6],dF_3D_ds1[5],dF_3D_ds1[3],dF_3D_ds1[4]};
                                        for (int p1 = 0; p1 < 9; p1++)
                                        {
                                            //StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += dFPK_tensor_dr_vec[p1] * dF_3D_ds1[p1] * w * wFactor;
                                            StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += dFPK_tensor_dr_vec[p1] * dF_3D_tr_ds1[p1] * w * wFactor;
                                            //StiffnessDevelop_v2[3 * i + r1, 3 * k + s1] += Pcontributions_term_1_dr_vec[p1] * dF_3D_ds1[p1] * w * wFactor;


                                            //stiffnessDevelop_v2_part2[3 * i + r1, 3 * k + s1] += dFPK_tensor_dr_vec[p1] * dF_3D_tr_ds1[p1] * w * wFactor;
                                        }
                                        #endregion


                                        if (ElementStiffnesses.performCalculations)
                                        {

                                            if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveStiffnessMatrixState && (i == 0) && (r1 == 1) && (i1 == 0))
                                            {
                                                double[,] dF_3D_drds = new double[3, 3] { { dg1_drds[0]*G_1[0]+dg2_drds[0]*G_2[0], dg1_drds[0]*G_1[1]+dg2_drds[0]*G_2[1], dg1_drds[0]*G_1[2]+dg2_drds[0]*G_2[2] },
                                                                 { dg1_drds[1]*G_1[0]+dg2_drds[1]*G_2[0], dg1_drds[1]*G_1[1]+dg2_drds[1]*G_2[1], dg1_drds[1]*G_1[2]+dg2_drds[1]*G_2[2] },
                                                                 { dg1_drds[2]*G_1[0]+dg2_drds[2]*G_2[0], dg1_drds[2]*G_1[1]+dg2_drds[2]*G_2[1], dg1_drds[2]*G_1[2]+dg2_drds[2]*G_2[2] }, };

                                                //double[] dF_3D_dr_vec =
                                                var dF_3D_drds_vecArray = new double[] { dF_3D_drds[0, 0], dF_3D_drds[1, 1], dF_3D_drds[2, 2], dF_3D_drds[0, 1], dF_3D_drds[1, 2], dF_3D_drds[2, 0], dF_3D_drds[0, 2], dF_3D_drds[1, 0], dF_3D_drds[2, 1] };

                                                ElementStiffnesses.saveOriginalState = true;
                                                ElementStiffnesses.ProccessVariable(27, dF_3D_drds_vecArray, true, 3 * k + s1);


                                                ElementStiffnesses.saveOriginalState = false;

                                            }

                                            if ((ElementStiffnesses.gpNumber == ElementStiffnesses.gpNumberToCheck) && ElementStiffnesses.saveStiffnessMatrixState && (i1 == 0))
                                            {

                                                double[,] Pcontributions = Calculate3DtensorFrom2Dcorrected(FPK_2D, de1_dr, de2_dr, da3_dr, Ei);//revisit this maybe we should pass dg de1_dr instead of dg1_dr.

                                                (double[,] Pcontributions_term_1, double[] dF2D_coefs_dr_vec, double[] dFPK2D_coefs_dr_vec_1, double[,] dFPK2D_coefs_dr_v2) = Calculate_FPK_variation_term1(materialDevelop.ConstitutiveMatrix.CopytoArray2D(), Vector.CreateFromArray(dg1_dr), Vector.CreateFromArray(dg2_dr),
                                                    Vector.CreateFromArray(de1_dr), Vector.CreateFromArray(de2_dr), G_1, G_2, a3r, Ei, ei, Vector.CreateFromArray(g1), Vector.CreateFromArray(g2), Vector.CreateFromArray(a3));//revisit this maybe we should pass dg de1_dr instead of dg1_dr.



                                                ElementStiffnesses.saveOriginalState = true;
                                                //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r00, a3r.a3r10, a3r.a3r20 }, true, 3 * i + 0);
                                                //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r01, a3r.a3r11, a3r.a3r21 }, true, 3 * i + 1);
                                                //ElementStiffnesses.ProccessVariable(8, new double[] { a3r.a3r02, a3r.a3r12, a3r.a3r22 }, true, 3 * i + 2); ;
                                                ElementStiffnesses.ProccessVariable(28, dF2D_coefs_dr_vec, true, 3 * i + r1);
                                                ElementStiffnesses.ProccessVariable(29, dFPK2D_coefs_dr_vec_1, true, 3 * i + r1);

                                                var Pcontributions_dr_vec = new double[] { Pcontributions[0, 0], Pcontributions[1, 1], Pcontributions[2, 2], Pcontributions[0, 1], Pcontributions[1, 2], Pcontributions[2, 0], Pcontributions[0, 2], Pcontributions[1, 0], Pcontributions[2, 1] };
                                                var Pcontributions_term_1_dr_vec = new double[] { Pcontributions_term_1[0, 0], Pcontributions_term_1[1, 1], Pcontributions_term_1[2, 2], Pcontributions_term_1[0, 1], Pcontributions_term_1[1, 2], Pcontributions_term_1[2, 0], Pcontributions_term_1[0, 2], Pcontributions_term_1[1, 0], Pcontributions_term_1[2, 1] };


                                                double[] dFPK_3D_dr_vec = new double[9];
                                                for (int j1 = 0; j1 < 9; j1++)
                                                {
                                                    dFPK_3D_dr_vec[j1] = Pcontributions_term_1_dr_vec[j1] + Pcontributions_dr_vec[j1];
                                                }

                                                ElementStiffnesses.ProccessVariable(30, new double[9], true, 3 * i + r1);
                                                ElementStiffnesses.ProccessVariable(33, dFPK_3D_dr_vec, true, 3 * i + r1);
                                                ElementStiffnesses.ProccessVariable(31, de1_dr, true, 3 * i + r1);

                                                ElementStiffnesses.ProccessVariable(32, de2_dr, true, 3 * i + r1);
                                                ElementStiffnesses.ProccessVariable(34, dFPK_tensor_dr_vec, true, 3 * i + r1);
                                                ElementStiffnesses.saveOriginalState = false;

                                            }


                                        }








                                    }
                                }


                            }
                        }
                    }



                }

            



            SymmetrizeMatrix(StiffnessDevelop_v2);
            return Matrix.CreateFromArray(StiffnessDevelop_v2);

        }

        #region stiffness matrix supportive methods implementation
        private (double[,][] da3_dksidrds, double[,][] da3_dhetadrds) Calculate_dvariable_drds_second_derivatives(double[] surfaceBasisVector1, double[] surfaceBasisVector2,
            double[] surfaceBasisVector3, double J1, double dKsi_r, double dKsi_s, double dHeta_r, double dHeta_s,
             double[][] da3tilde_dksidr, double[][] da3tilde_dhetadr, double[][] da3tilde_dksids, double[][] da3tilde_dhetads, double[] a3_tilde, double[] da3tilde_dksi, double[] da3tilde_dheta, double da3norm_dksi,
            double da3norm_dheta, double[] da3norm_dksids, double[] da3norm_dhetads, double[] da3norm_dhetadr, double[] da3norm_dksidr,
            double dksi_r, double dheta_r, double d2Ksi_dr, double d2Heta_dr, double d2KsiHeta_dr,
            double d2Ksi_ds, double d2Heta_ds, double d2KsiHeta_ds, double dksi_s, double dheta_s, double[,][] da3_drds, double[,][] da3tilde_drds)
        {
            double J1_squared = Math.Pow(J1, 2);
            var J1_cubed = Math.Pow(J1, 3);
            double J1_fourth = Math.Pow(J1, 4);

            #region Calculate_a3rs
            //var da3_drds = new double[3, 3][];
            //var da3tilde_drds = new double[3, 3][];
            var da3tilde_dr = new double[3][];
            var da3tilde_ds = new double[3][];
            //double[] a3_tilde;
            double[] dnorma3_dr = new double[3];
            double[] dnorma3_ds = new double[3];
            double[,] dnorma3_drds = new double[3, 3];

            //5.30
            Calculate_da3tilde_drds(dKsi_r, dKsi_s, dHeta_r, dHeta_s, da3tilde_drds);

            //5.24
            Calculate_da3tilde_dr_Array(surfaceBasisVector1, surfaceBasisVector2, dKsi_r, dHeta_r, da3tilde_dr);

            //5.24
            Calculate_da3tilde_dr_Array(surfaceBasisVector1, surfaceBasisVector2, dKsi_s, dHeta_s, da3tilde_ds);


            //5.25
            CalculateTerm525(surfaceBasisVector3, J1, dnorma3_dr, da3tilde_dr, dnorma3_ds, da3tilde_ds);

            //5.31
            CalculateTerm531(J1, da3tilde_drds, a3_tilde, da3tilde_dr, da3tilde_ds, dnorma3_drds);

            //5.32
            CalculateTerm532(J1, da3tilde_drds, dnorma3_ds, da3tilde_dr, dnorma3_dr, da3tilde_ds, dnorma3_drds, a3_tilde, da3_drds);//.

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

            #endregion

            #region 530_c


            var da3tilde_dksidrds = new double[3, 3][];
            var da3tilde_dhetadrds = new double[3, 3][];
            //for (int r1 = 0; r1 < 3; r1++)
            //{
            //    for (int s1 = 0; s1 < 3; s1++)
            //    {
            //        da3tilde_dksidrds[r1, s1] = a11r.GetColumn(r1).CrossProduct(a2s.GetColumn(s1)) + a11s.GetColumn(s1).CrossProduct(a2r.GetColumn(r1)) +
            //            a1r.GetColumn(r1).CrossProduct(a12s.GetColumn(s1)) + a1s.GetColumn(s1).CrossProduct(a12r.GetColumn(r1));

            //        da3tilde_dhetadrds[r1, s1] = a12r.GetColumn(r1).CrossProduct(a2s.GetColumn(s1)) + a12s.GetColumn(s1).CrossProduct(a2r.GetColumn(r1)) +
            //            a1r.GetColumn(r1).CrossProduct(a22s.GetColumn(s1)) + a1s.GetColumn(s1).CrossProduct(a22r.GetColumn(r1));
            //    }
            //}
            double aux1 = d2Ksi_dr * dheta_s - d2Ksi_ds * dheta_r - d2KsiHeta_dr * dksi_s + d2KsiHeta_ds * dksi_r;
            double aux2 = d2Ksi_ds * dheta_r - d2Ksi_dr * dheta_s + d2KsiHeta_dr * dksi_s - d2KsiHeta_ds * dksi_r;
            da3tilde_dksidrds[0, 0] = new double[3];
            da3tilde_dksidrds[0, 1] = new double[3] { 0, 0, aux1 };
            da3tilde_dksidrds[0, 2] = new double[3] { 0, aux2, 0 };

            da3tilde_dksidrds[1, 0] = new double[3] { 0, 0, aux2 };
            da3tilde_dksidrds[1, 1] = new double[3] { 0, 0, 0 };
            da3tilde_dksidrds[1, 2] = new double[3] { aux1, 0, 0 };

            da3tilde_dksidrds[2, 0] = new double[3] { 0, aux1, 0 };
            da3tilde_dksidrds[2, 1] = new double[3] { aux2, 0, 0 };
            da3tilde_dksidrds[2, 2] = new double[3] { 0, 0, 0 };

            double aux3 = d2KsiHeta_dr * dheta_s - d2KsiHeta_ds * dheta_r - d2Heta_dr * dksi_s + d2Heta_ds * dksi_r;
            double aux4 = d2KsiHeta_ds * dheta_r - d2KsiHeta_dr * dheta_s + d2Heta_dr * dksi_s - d2Heta_ds * dksi_r;
            da3tilde_dhetadrds[0, 0] = new double[3] { 0, 0, 0 };
            da3tilde_dhetadrds[0, 1] = new double[3] { 0, 0, aux3 };
            da3tilde_dhetadrds[0, 2] = new double[3] { 0, aux4, 0 };

            da3tilde_dhetadrds[1, 0] = new double[3] { 0, 0, aux4 };
            da3tilde_dhetadrds[1, 1] = new double[3] { 0, 0, 0 };
            da3tilde_dhetadrds[1, 2] = new double[3] { aux3, 0, 0 };

            da3tilde_dhetadrds[2, 0] = new double[3] { 0, aux3, 0 };
            da3tilde_dhetadrds[2, 1] = new double[3] { aux4, 0, 0 };
            da3tilde_dhetadrds[2, 2] = new double[3] { 0, 0, 0 };

            #endregion

            #region 531_c
            var da3norm_dksidrds = new double[3, 3];
            var da3norm_dhetadrds = new double[3, 3];

            //double term05a = da3tilde_dksi.DotProduct(a3_tilde);
            double term05a = da3tilde_dksi[0] * a3_tilde[0] + da3tilde_dksi[1] * a3_tilde[1] + da3tilde_dksi[2] * a3_tilde[2];
            for (int r1 = 0; r1 < 3; r1++)
            {
                //double term01 = da3tilde_dksidr[r1].DotProduct(a3_tilde) + da3tilde_dksi.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double term01 = da3tilde_dksidr[r1][0] * a3_tilde[0] + da3tilde_dksidr[r1][1] * a3_tilde[1] + da3tilde_dksidr[r1][2] * a3_tilde[2] +
                    da3tilde_dksi[0] * da3tilde_dr[r1][0] + da3tilde_dksi[1] * da3tilde_dr[r1][1] + da3tilde_dksi[2] * da3tilde_dr[r1][2];


                //double TERM03 = (da3tilde_dksi.DotProduct(a3_tilde)) * (da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray()));
                double TERM03 = (da3tilde_dksi[0] * a3_tilde[0] + da3tilde_dksi[1] * a3_tilde[1] + da3tilde_dksi[2] * a3_tilde[2]) *
                    (da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] + da3tilde_dr[r1][2] * a3_tilde[2]);


                //double term04b = da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray());
                double term04b = da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] + da3tilde_dr[r1][2] * a3_tilde[2];

                for (int s1 = 0; s1 < 3; s1++)
                {
                    //double term02 = da3tilde_dksidrds[r1, s1].DotProduct(a3_tilde) + da3tilde_dksidr[r1].CopyToArray().DotProduct(da3tilde_ds[s1]) +
                    //da3tilde_dksids[s1].CopyToArray().DotProduct(da3tilde_dr[r1]) + da3tilde_dksi.CopyToArray().DotProduct(da3tilde_drds[r1, s1]);
                    double term02 = da3tilde_dksidrds[r1, s1][0] * a3_tilde[0] + da3tilde_dksidrds[r1, s1][1] * a3_tilde[1] + da3tilde_dksidrds[r1, s1][2] * a3_tilde[2] +
                        da3tilde_dksidr[r1][0] * da3tilde_ds[s1][0] + da3tilde_dksidr[r1][1] * da3tilde_ds[s1][1] + da3tilde_dksidr[r1][2] * da3tilde_ds[s1][2] +
                        da3tilde_dksids[s1][0] * da3tilde_dr[r1][0] + da3tilde_dksids[s1][1] * da3tilde_dr[r1][1] + da3tilde_dksids[s1][2] * da3tilde_dr[r1][2] +
                        da3tilde_dksi[0] * da3tilde_drds[r1, s1][0] + da3tilde_dksi[1] * da3tilde_drds[r1, s1][1] + da3tilde_dksi[2] * da3tilde_drds[r1, s1][2];

                    //double term04a = da3tilde_dksids[s1].DotProduct(a3_tilde) + da3tilde_dksi.CopyToArray().DotProduct(da3tilde_ds[s1]);
                    double term04a = da3tilde_dksids[s1][0] * a3_tilde[0] + da3tilde_dksids[s1][1] * a3_tilde[1] + da3tilde_dksids[s1][2] * a3_tilde[2] +
                        da3tilde_dksi[0] * da3tilde_ds[s1][0] + da3tilde_dksi[1] * da3tilde_ds[s1][1] + da3tilde_dksi[2] * da3tilde_ds[s1][2];

                    //double term05b = da3tilde_drds[r1, s1].DotProduct(a3_tilde.CopyToArray()) + da3tilde_dr[r1].DotProduct(da3tilde_ds[s1]);
                    double term05b = da3tilde_drds[r1, s1][0] * a3_tilde[0] + da3tilde_drds[r1, s1][1] * a3_tilde[1] + da3tilde_drds[r1, s1][2] * a3_tilde[2] +
                        da3tilde_dr[r1][0] * da3tilde_ds[s1][0] + da3tilde_dr[r1][1] * da3tilde_ds[s1][1] + da3tilde_dr[r1][2] * da3tilde_ds[s1][2];

                    da3norm_dksidrds[r1, s1] = (term02 / J1) - (term01 * dnorma3_ds[s1] / (J1_squared)) -
                        ((double)1 / (J1_cubed)) * ((term04a * term04b) + (term05a * term05b)) + TERM03 * 3 * dnorma3_ds[s1] * ((double)1 / J1_fourth);
                }
            }

            //term05a = da3tilde_dheta.DotProduct(a3_tilde);
            term05a = da3tilde_dheta[0] * a3_tilde[0] + da3tilde_dheta[1] * a3_tilde[1] + da3tilde_dheta[2] * a3_tilde[2];
            for (int r1 = 0; r1 < 3; r1++)
            {
                //double term01 = da3tilde_dhetadr[r1].DotProduct(a3_tilde) + da3tilde_dheta.DotProduct(Vector.CreateFromArray(da3tilde_dr[r1]));
                double term01 = da3tilde_dhetadr[r1][0] * a3_tilde[0] + da3tilde_dhetadr[r1][1] * a3_tilde[1] + da3tilde_dhetadr[r1][2] * a3_tilde[2] +
                    da3tilde_dheta[0] * da3tilde_dr[r1][0] + da3tilde_dheta[1] * da3tilde_dr[r1][1] + da3tilde_dheta[2] * da3tilde_dr[r1][2];
                //double TERM03 = (da3tilde_dheta.DotProduct(a3_tilde)) * (da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray()));
                double TERM03 = da3tilde_dheta[0] * a3_tilde[0] + da3tilde_dheta[1] * a3_tilde[1] + da3tilde_dheta[2] * a3_tilde[2] +
                    da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] + da3tilde_dr[r1][2] * a3_tilde[2];
                //double term04b = da3tilde_dr[r1].DotProduct(a3_tilde.CopyToArray());
                double term04b = da3tilde_dr[r1][0] * a3_tilde[0] + da3tilde_dr[r1][1] * a3_tilde[1] + da3tilde_dr[r1][2] * a3_tilde[2];
                for (int s1 = 0; s1 < 3; s1++)
                {
                    //double term02 = da3tilde_dhetadrds[r1, s1].DotProduct(a3_tilde) + da3tilde_dhetadr[r1].CopyToArray().DotProduct(da3tilde_ds[s1]) +
                    //    da3tilde_dhetads[s1].CopyToArray().DotProduct(da3tilde_dr[r1]) + da3tilde_dheta.CopyToArray().DotProduct(da3tilde_drds[r1, s1]);

                    double term02 = da3tilde_dhetadrds[r1, s1][0] * a3_tilde[0] + da3tilde_dhetadrds[r1, s1][1] * a3_tilde[1] + da3tilde_dhetadrds[r1, s1][2] * a3_tilde[2] +
                        da3tilde_dhetadr[r1][0] * da3tilde_ds[s1][0] + da3tilde_dhetadr[r1][1] * da3tilde_ds[s1][1] + da3tilde_dhetadr[r1][2] * da3tilde_ds[s1][2] +
                        da3tilde_dhetads[s1][0] * da3tilde_dr[r1][0] + da3tilde_dhetads[s1][1] * da3tilde_dr[r1][1] + da3tilde_dhetads[s1][2] * da3tilde_dr[r1][2] +
                        da3tilde_dheta[0] * da3tilde_drds[r1, s1][0] + da3tilde_dheta[1] * da3tilde_drds[r1, s1][1] + da3tilde_dheta[2] * da3tilde_drds[r1, s1][2];



                    //double term04a = da3tilde_dhetads[s1].DotProduct(a3_tilde) + da3tilde_dheta.CopyToArray().DotProduct(da3tilde_ds[s1]);
                    double term04a = da3tilde_dhetads[s1][0] * a3_tilde[0] + da3tilde_dhetads[s1][1] * a3_tilde[1] + da3tilde_dhetads[s1][2] * a3_tilde[2] +
                        da3tilde_dheta[0] * da3tilde_ds[s1][0] + da3tilde_dheta[1] * da3tilde_ds[s1][1] + da3tilde_dheta[2] * da3tilde_ds[s1][2];



                    //double term05b = da3tilde_drds[r1, s1].DotProduct(a3_tilde.CopyToArray()) + da3tilde_dr[r1].DotProduct(da3tilde_ds[s1]);
                    double term05b = (da3tilde_drds[r1, s1][0] * a3_tilde[0] + da3tilde_drds[r1, s1][1] * a3_tilde[1] + da3tilde_drds[r1, s1][2] * a3_tilde[2]) +
                                     (da3tilde_dr[r1][0] * da3tilde_ds[s1][0] + da3tilde_dr[r1][1] * da3tilde_ds[s1][1] + da3tilde_dr[r1][2] * da3tilde_ds[s1][2]);

                    da3norm_dhetadrds[r1, s1] = (term02 / J1) - (term01 * dnorma3_ds[s1] / (J1_squared)) -
                        ((double)1 / (J1_cubed)) * ((term04a * term04b) + (term05a * term05b)) + TERM03 * 3 * dnorma3_ds[s1] * ((double)1 / J1_fourth);
                }
            }
            #endregion

            #region 532_c
            double[,][] da3_dksidrds = new double[3, 3][];
            double[,][] da3_dhetadrds = new double[3, 3][];//.





            double t_ars_coeff = (double)1 / J1;
            for (int r1 = 0; r1 < 3; r1++)
            {
                double t_as_coeff = -dnorma3_dr[r1] / J1_squared;
                double t_s_coeff = (2 * da3norm_dksi * dnorma3_dr[r1] / J1_cubed) - (da3norm_dksidr[r1] / J1_squared);
                for (int s1 = 0; s1 < 3; s1++)
                {
                    double t_ar_coeff = -dnorma3_ds[s1] / J1_squared;
                    double t_a_coeff = (2 * dnorma3_dr[r1] * dnorma3_ds[s1] / J1_cubed) - dnorma3_drds[r1, s1] / J1_squared;
                    double t_rs_coeff = -da3norm_dksi / J1_squared;
                    double t_r_coeff = (2 * da3norm_dksi * dnorma3_ds[s1] / J1_cubed) - (da3norm_dksids[s1] / J1_squared);

                    double t_coef = -(6 * da3norm_dksi * dnorma3_dr[r1] * dnorma3_ds[s1] / J1_fourth) +
                                    (2 * da3norm_dksids[s1] * dnorma3_dr[r1] / J1_cubed) +
                                    (2 * da3norm_dksi * dnorma3_drds[r1, s1] / J1_cubed) +
                                    (2 * da3norm_dksidr[r1] * dnorma3_ds[s1] / J1_cubed) -
                                    (da3norm_dksidrds[r1, s1] / J1_squared);

                    da3_dksidrds[r1, s1] = new double[3]{
                                        da3tilde_dksidrds[r1, s1][0]*(t_ars_coeff) +
                                        da3tilde_dksidr[r1][0]*(t_ar_coeff) + da3tilde_dksids[s1][0]*(t_as_coeff) + da3tilde_dksi[0]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][0]*(t_rs_coeff) + da3tilde_dr[r1][0]*(t_r_coeff) + da3tilde_ds[s1][0]*(t_s_coeff) +
                                        a3_tilde[0]*(t_coef),

                                        da3tilde_dksidrds[r1, s1][1]*(t_ars_coeff) +
                                        da3tilde_dksidr[r1][1]*(t_ar_coeff) + da3tilde_dksids[s1][1]*(t_as_coeff) + da3tilde_dksi[1]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][1]*(t_rs_coeff) + da3tilde_dr[r1][1]*(t_r_coeff) + da3tilde_ds[s1][1]*(t_s_coeff) +
                                        a3_tilde[1]*(t_coef),
                                        da3tilde_dksidrds[r1, s1][2]*(t_ars_coeff) +
                                        da3tilde_dksidr[r1][2]*(t_ar_coeff) + da3tilde_dksids[s1][2]*(t_as_coeff) + da3tilde_dksi[2]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][2]*(t_rs_coeff) + da3tilde_dr[r1][2]*(t_r_coeff) + da3tilde_ds[s1][2]*(t_s_coeff) +
                                        a3_tilde[2]*(t_coef)
                    };


                }
            }

            t_ars_coeff = (double)1 / J1;
            for (int r1 = 0; r1 < 3; r1++)
            {
                double t_as_coeff = -dnorma3_dr[r1] / J1_squared;
                double t_s_coeff = (2 * da3norm_dheta * dnorma3_dr[r1] / J1_cubed) - (da3norm_dhetadr[r1] / J1_squared);
                for (int s1 = 0; s1 < 3; s1++)
                {
                    double t_ar_coeff = -dnorma3_ds[s1] / J1_squared;
                    double t_a_coeff = (2 * dnorma3_dr[r1] * dnorma3_ds[s1] / J1_cubed) - dnorma3_drds[r1, s1] / J1_squared;
                    double t_rs_coeff = -da3norm_dheta / J1_squared;
                    double t_r_coeff = (2 * da3norm_dheta * dnorma3_ds[s1] / J1_cubed) - (da3norm_dhetads[s1] / J1_squared);

                    double t_coef = -(6 * da3norm_dheta * dnorma3_dr[r1] * dnorma3_ds[s1] / J1_fourth) +
                                    (2 * da3norm_dhetads[s1] * dnorma3_dr[r1] / J1_cubed) +
                                    (2 * da3norm_dheta * dnorma3_drds[r1, s1] / J1_cubed) +
                                    (2 * da3norm_dhetadr[r1] * dnorma3_ds[s1] / J1_cubed) -
                                    (da3norm_dhetadrds[r1, s1] / J1_squared);

                    da3_dhetadrds[r1, s1] = new double[3]{ da3tilde_dhetadrds[r1, s1][0]*(t_ars_coeff) +
                                        da3tilde_dhetadr[r1][0]*(t_ar_coeff) + da3tilde_dhetads[s1][0]*(t_as_coeff) + da3tilde_dheta[0]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][0]*(t_rs_coeff) + da3tilde_dr[r1][0]*(t_r_coeff) + da3tilde_ds[s1][0]*(t_s_coeff) +
                                        a3_tilde[0]*(t_coef),
                                        da3tilde_dhetadrds[r1, s1][1]*(t_ars_coeff) +
                                        da3tilde_dhetadr[r1][1]*(t_ar_coeff) + da3tilde_dhetads[s1][1]*(t_as_coeff) + da3tilde_dheta[1]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][1]*(t_rs_coeff) + da3tilde_dr[r1][1]*(t_r_coeff) + da3tilde_ds[s1][1]*(t_s_coeff) +
                                        a3_tilde[1]*(t_coef),
                                        da3tilde_dhetadrds[r1, s1][2]*(t_ars_coeff) +
                                        da3tilde_dhetadr[r1][2]*(t_ar_coeff) + da3tilde_dhetads[s1][2]*(t_as_coeff) + da3tilde_dheta[2]*(t_a_coeff) +
                                        da3tilde_drds[r1, s1][2]*(t_rs_coeff) + da3tilde_dr[r1][2]*(t_r_coeff) + da3tilde_ds[s1][2]*(t_s_coeff) +
                                        a3_tilde[2]*(t_coef) };
                }
            }
            #endregion

            return (da3_dksidrds, da3_dhetadrds);

        }

        private double[] CalculateDerivativeOfOrthogonalisedNormalisedCovariantBasisArray(double[,] ei, double[] de1_dr, double[] e2_init, double[] de2_init_dr)
        {
            double e1_dot_e2_init = ei[0, 0] * e2_init[0] + ei[1, 0] * e2_init[1] + ei[2, 0] * e2_init[2];
            double e1_dot_de2_init_dr = ei[0, 0] * de2_init_dr[0] + ei[1, 0] * de2_init_dr[1] + ei[2, 0] * de2_init_dr[2];
            double de1_dr_dot_e2_init = de1_dr[0] * e2_init[0] + de1_dr[1] * e2_init[1] + de1_dr[2] * e2_init[2];

            //Vector projection = e1.Scale(e1_dot_e2_init);
            double[] projection = new double[] { ei[0, 0] * e1_dot_e2_init, ei[1, 0] * e1_dot_e2_init, ei[2, 0] * e1_dot_e2_init };

            //Vector dprojection_dr = e1.Scale(e1_dot_de2_init_dr + de1_dr_dot_e2_init) + de1_dr.Scale(e1_dot_e2_init);
            double e1coeff = e1_dot_de2_init_dr + de1_dr_dot_e2_init;
            double[] dprojection_dr = new double[] { ei[0, 0] * e1coeff + de1_dr[0] * e1_dot_e2_init, ei[1, 0] * e1coeff + de1_dr[1] * e1_dot_e2_init, ei[2, 0] * e1coeff + de1_dr[2] * e1_dot_e2_init };


            double[] e2_tilde = new double[] { e2_init[0] - projection[0], e2_init[1] - projection[1], e2_init[2] - projection[2] };

            double[] de2tilde_dr = new double[] { de2_init_dr[0] - dprojection_dr[0], de2_init_dr[1] - dprojection_dr[1], de2_init_dr[2] - dprojection_dr[2] };

            double[] de2_dr = CalculateDerivativeOfVectorNormalisedArray(e2_tilde, de2tilde_dr);

            return de2_dr;
        }

        public struct dFrvedri
        {
            public double r00;
            public double r01;
            public double r02;

            public double r10;
            public double r11;
            public double r12;

            public double r20;
            public double r21;
            public double r22;
        }

        private void Calculate_FPK_variation_term1_Short(IMatrixView Aijkl_2D, double[] dg1_dr, double[] dg2_dr, double[] de1_dr, double[] de2_dr, double[] G_1, double[] G_2, a3r a3r, double[,] Ei, double[,] ei, double[] g1, double[] g2, double[] dF2D_coefs_dr_vec,
            double[] dFPK2D_coefs_dr_vec,
            ref dFrvedri dFPK2D_coefs_dr)//.
        {
            //double[,] dF2D_coefs_dr = new double[2, 2];

            double dg1_dr__ei0 = (dg1_dr[0] * ei[0, 0] + dg1_dr[1] * ei[1, 0] + dg1_dr[2] * ei[2, 0] + g1[0] * de1_dr[0] + g1[1] * de1_dr[1] + g1[2] * de1_dr[2]);
            double G_1__Ei0 = (G_1[0] * Ei[0, 0] + G_1[1] * Ei[1, 0] + G_1[2] * Ei[2, 0]);
            double dg2_dr__ei0 = (dg2_dr[0] * ei[0, 0] + dg2_dr[1] * ei[1, 0] + dg2_dr[2] * ei[2, 0] + g2[0] * de1_dr[0] + g2[1] * de1_dr[1] + g2[2] * de1_dr[2]);
            double G_2__Ei0 = (G_2[0] * Ei[0, 0] + G_2[1] * Ei[1, 0] + G_2[2] * Ei[2, 0]);
            double dg1_dr_ei1 = (dg1_dr[0] * ei[0, 1] + dg1_dr[1] * ei[1, 1] + dg1_dr[2] * ei[2, 1] + g1[0] * de2_dr[0] + g1[1] * de2_dr[1] + g1[2] * de2_dr[2]);
            double dg2_dr_ei1 = (dg2_dr[0] * ei[0, 1] + dg2_dr[1] * ei[1, 1] + dg2_dr[2] * ei[2, 1] + g2[0] * de2_dr[0] + g2[1] * de2_dr[1] + g2[2] * de2_dr[2]);
            double G_2__Ei1 = (G_2[0] * Ei[0, 1] + G_2[1] * Ei[1, 1] + G_2[2] * Ei[2, 1]);
            double G_1__Ei1 = (G_1[0] * Ei[0, 1] + G_1[1] * Ei[1, 1] + G_1[2] * Ei[2, 1]);

            //dF2D_coefs_dr[0, 0] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g1.DotProduct(de1_dr)) * (G_1.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] })) +
            //                        (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g2.DotProduct(de1_dr)) * (G_2.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }));

            dF2D_coefs_dr_vec[0] = dg1_dr__ei0 * G_1__Ei0 +
                dg2_dr__ei0 * G_2__Ei0;

            //dF2D_coefs_dr[0, 1] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g1.DotProduct(de1_dr)) * (G_1.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] })) +
            //                        (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g2.DotProduct(de1_dr)) * (G_2.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }));

            dF2D_coefs_dr_vec[2] = dg1_dr__ei0 * G_1__Ei1 +
                dg2_dr__ei0 * G_2__Ei1;

            //dF2D_coefs_dr[1, 0] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g1.DotProduct(de2_dr)) * (G_1.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] })) +
            //                        (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g2.DotProduct(de2_dr)) * (G_2.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }));

            dF2D_coefs_dr_vec[3] = dg1_dr_ei1 * G_1__Ei0 +
                                  dg2_dr_ei1 * G_2__Ei0;

            //dF2D_coefs_dr[1, 1] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g1.DotProduct(de2_dr)) * (G_1.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] })) +
            //                        (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g2.DotProduct(de2_dr)) * (G_2.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }));

            dF2D_coefs_dr_vec[1] = dg1_dr_ei1 * G_1__Ei1 +
                                  dg2_dr_ei1 * G_2__Ei1;

            //double[] dF2D_coefs_dr_vec = new double[] { dF2D_coefs_dr[0, 0], dF2D_coefs_dr[1, 1], dF2D_coefs_dr[0, 1], dF2D_coefs_dr[1, 0] }; // to [1,0] vgainei iso  me to [0,1]


            Array.Clear(dFPK2D_coefs_dr_vec, 0, 4);
            for (int i1 = 0; i1 < 4; i1++)
            {

                for (int i2 = 0; i2 < 4; i2++)
                {
                    dFPK2D_coefs_dr_vec[i1] += Aijkl_2D[i1, i2] * dF2D_coefs_dr_vec[i2];
                }
            }

            dFPK2D_coefs_dr.r00 = dFPK2D_coefs_dr_vec[0]; dFPK2D_coefs_dr.r01 = dFPK2D_coefs_dr_vec[2]; dFPK2D_coefs_dr.r10 = dFPK2D_coefs_dr_vec[3]; dFPK2D_coefs_dr.r11 = dFPK2D_coefs_dr_vec[1];


            //.return dFPK2D_coefs_dr;
        }

        private double[,] Calculate3DtensorFrom2Dcorrected(double[,] FPK_2D, double[] dg1_dr, double[] dg2_dr, double[] da3_dr, double[,] Ei)
        {
            double[,] eye = new double[3, 3]; eye[0, 0] = 1; eye[1, 1] = 1; eye[2, 2] = 1;

            #region create and normalise ei
            double[,] ei = new double[3, 3];
            double norm_e1 = dg1_dr[0] * dg1_dr[0] + dg1_dr[1] * dg1_dr[1] + dg1_dr[2] * dg1_dr[2];
            norm_e1 = Math.Sqrt(norm_e1);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 0] = dg1_dr[i2] / norm_e1;
            }

            double norm_e2 = dg2_dr[0] * dg2_dr[0] + dg2_dr[1] * dg2_dr[1] + dg2_dr[2] * dg2_dr[2];
            norm_e2 = Math.Sqrt(norm_e2);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 1] = dg2_dr[i2] / norm_e2;
            }

            double norm_e3 = da3_dr[0] * da3_dr[0] + da3_dr[1] * da3_dr[1] + da3_dr[2] * da3_dr[2];
            norm_e3 = Math.Sqrt(norm_e3);
            for (int i2 = 0; i2 < 3; i2++)
            {
                ei[i2, 2] = da3_dr[i2] / norm_e3;
            }
            #endregion

            #region adapt FPK_2D for normalisation of basis vectors

            double coef1 = 0;
            double coef2 = 0;
            double[,] FPK_2D_in_normalised = new double[3, 3];
            for (int i1 = 0; i1 < 2; i1++)
            {
                if (i1 == 0) { coef1 = norm_e1; }
                else if (i1 == 1) { coef1 = norm_e2; }
                for (int i2 = 0; i2 < 2; i2++)
                {
                    if (i2 == 0) { coef2 = norm_e1; }
                    else if (i2 == 1) { coef2 = norm_e2; }

                    //FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1 * coef2;
                    FPK_2D_in_normalised[i1, i2] = FPK_2D[i1, i2] * coef1;// * coef2;
                }
            }
            #endregion

            var cartes_to_Gi = CalculateRotationMatrix(Ei, eye);
            var cartes_to_tgi = CalculateRotationMatrix(ei, eye);

            double[,] FPK_3D = Transform_FPK_rve_To_FPK_3D(FPK_2D_in_normalised, cartes_to_Gi, cartes_to_tgi);// 1);



            return FPK_3D;
        }

        private (double[,], double[], double[], double[,]) Calculate_FPK_variation_term1(double[,] Aijkl_2D, Vector dg1_dr, Vector dg2_dr, Vector de1_dr, Vector de2_dr, double[] G_1, double[] G_2, a3r a3r, double[,] Ei, double[,] ei, Vector g1, Vector g2, Vector a3)
        {
            double[,] dF2D_coefs_dr = new double[2, 2];

            dF2D_coefs_dr[0, 0] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g1.DotProduct(de1_dr)) * (G_1.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] })) +
                                    (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g2.DotProduct(de1_dr)) * (G_2.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }));

            dF2D_coefs_dr[0, 1] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g1.DotProduct(de1_dr)) * (G_1.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] })) +
                                    (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 0], ei[1, 0], ei[2, 0] }) + g2.DotProduct(de1_dr)) * (G_2.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }));

            dF2D_coefs_dr[1, 0] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g1.DotProduct(de2_dr)) * (G_1.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] })) +
                                    (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g2.DotProduct(de2_dr)) * (G_2.DotProduct(new double[] { Ei[0, 0], Ei[1, 0], Ei[2, 0] }));

            dF2D_coefs_dr[1, 1] = (dg1_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g1.DotProduct(de2_dr)) * (G_1.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] })) +
                                    (dg2_dr.CopyToArray().DotProduct(new double[3] { ei[0, 1], ei[1, 1], ei[2, 1] }) + g2.DotProduct(de2_dr)) * (G_2.DotProduct(new double[] { Ei[0, 1], Ei[1, 1], Ei[2, 1] }));

            double[] dF2D_coefs_dr_vec = new double[] { dF2D_coefs_dr[0, 0], dF2D_coefs_dr[1, 1], dF2D_coefs_dr[0, 1], dF2D_coefs_dr[1, 0] }; // to [1,0] vgainei iso  me to [0,1]

            double[] dFPK2D_coefs_dr_vec = new double[4];


            for (int i1 = 0; i1 < 4; i1++)
            {

                for (int i2 = 0; i2 < 4; i2++)
                {
                    dFPK2D_coefs_dr_vec[i1] += Aijkl_2D[i1, i2] * dF2D_coefs_dr_vec[i2];
                }
            }

            var dFPK2D_coefs_dr = new double[,] { { dFPK2D_coefs_dr_vec[0], dFPK2D_coefs_dr_vec[2] }, { dFPK2D_coefs_dr_vec[3], dFPK2D_coefs_dr_vec[1] } };

            //double[,] Pcontributions = Calculate3DtensorFrom2D(dFPK2D_coefs_dr, Vector.CreateFromArray(new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }), Vector.CreateFromArray(new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }),
            //   Ei);
            double[,] Pcontributions = Calculate3DtensorFrom2Dcorrected(dFPK2D_coefs_dr, new double[] { ei[0, 0], ei[1, 0], ei[2, 0] }, new double[] { ei[0, 1], ei[1, 1], ei[2, 1] }, a3.CopyToArray(), Ei);

            return (Pcontributions, dF2D_coefs_dr_vec, dFPK2D_coefs_dr_vec, dFPK2D_coefs_dr);
        }

        private void SymmetrizeMatrix(double[,] matrix)
        {
            for (int i1 = 1; i1 < matrix.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < i1; i2++)
                {
                    matrix[i1, i2] = matrix[i2, i1];

                    if (double.IsNaN(matrix[i1, i2]))
                    {
                        var breakpoint1 = "here";
                    }
                }
            }
        }

        private static double[] CalculateTerm525(double[] surfaceBasisVector3, double J1, double[] dnorma3_dr,
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
        #endregion


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

                //foreach (var integrationPointMaterial in materialsAtThicknessGP[gaussPoints[j]].Values)
                //{
                //    integrationPointMaterial.TangentVectorV1 = initialSurfaceBasisVectors1[j];
                //    integrationPointMaterial.TangentVectorV2 = initialSurfaceBasisVectors2[j];
                //    integrationPointMaterial.NormalVectorV3 = initialUnitSurfaceBasisVectors3[j];
                //}
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
        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(KirchhoffLoveShellNLDefGrad shellElement)
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