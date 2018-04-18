﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Tensors;

//TODO: Replace double[] with Vector
namespace ISAAR.MSolve.XFEM.CrackGeometry
{
    //TODO: decide what data structures (arrays, tuples, matrix & vector classes I will use as arguments, return types 
    // and for operations. Implement convenience methods for those operations on these data types.
    // Perhaps vector-vector operations could be abstracted. 
    // Actually wouldn't all methods be clearer if I operated directly with cosa, sina, instead of rotation matrices?
    class TipCoordinateSystem
    {
        private readonly double[] localCoordinatesOfGlobalOrigin;

        public double RotationAngle { get; }
        public Matrix RotationMatrixGlobalToLocal { get; }
        public Matrix TransposeRotationMatrixGlobalToLocal { get; } // cache this for efficiency

        /// <summary>
        /// det(J_globToLoc) = det(Q) = (cosa)^2 + (sina)^2 = 1
        /// </summary>
        public double DeterminantOfJacobianGlobalToLocalCartesian { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tipCoordinates">Coordinates of the crack tip in the global cartesian system.</param>
        /// <param name="tipRotationAngle">Counter-clockwise angle from the O-x axis of the global cartesian system to  
        ///     the T-x1 axis of the local corrdinate system of the tip (T being the tip point)</param>
        public TipCoordinateSystem(ICartesianPoint2D tipCoordinates, double tipRotationAngle)
        {
            this.RotationAngle = tipRotationAngle;

            double cosa = Math.Cos(tipRotationAngle);
            double sina = Math.Sin(tipRotationAngle);
            RotationMatrixGlobalToLocal = Matrix.CreateFromArray(new double[,] { { cosa, sina }, { -sina, cosa } });
            TransposeRotationMatrixGlobalToLocal = RotationMatrixGlobalToLocal.Transpose();
            localCoordinatesOfGlobalOrigin = RotationMatrixGlobalToLocal.
                MultiplyRight(VectorMKL.CreateFromArray(new double[] { -tipCoordinates.X, -tipCoordinates.Y }), false).
                CopyToArray(); //TODO: Replace double[] with Vector

            DeterminantOfJacobianGlobalToLocalCartesian = 1.0; // det = (cosa)^2 +(sina)^2 = 1
        }

        public ICartesianPoint2D TransformPointGlobalCartesianToLocalCartesian(ICartesianPoint2D cartesianGlobalPoint)
        {
            double[] rotated = (RotationMatrixGlobalToLocal * VectorMKL.CreateFromArray(cartesianGlobalPoint.Coordinates, false)).
                CopyToArray();
            double x1 = rotated[0] + localCoordinatesOfGlobalOrigin[0];
            double x2 = rotated[1] + localCoordinatesOfGlobalOrigin[1];
            return new CartesianPoint2D(x1, x2);
        }

        public PolarPoint2D TransformPointLocalCartesianToLocalPolar(ICartesianPoint2D cartesianLocalPoint)
        {
            double x1 = cartesianLocalPoint.X;
            double x2 = cartesianLocalPoint.Y;
            double r = Math.Sqrt(x1 * x1 + x2 * x2);
            double theta = Math.Atan2(x2, x1);
            return new PolarPoint2D(r, theta);
        }

        public PolarPoint2D TransformPointGlobalCartesianToLocalPolar(ICartesianPoint2D cartesianGlobalPoint)
        {
            double[] rotated = (RotationMatrixGlobalToLocal * VectorMKL.CreateFromArray(cartesianGlobalPoint.Coordinates, false)).
                CopyToArray();
            double x1 = rotated[0] + localCoordinatesOfGlobalOrigin[0];
            double x2 = rotated[1] + localCoordinatesOfGlobalOrigin[1];
            double r = Math.Sqrt(x1 * x1 + x2 * x2);
            double theta = Math.Atan2(x2, x1);
            return new PolarPoint2D(r, theta);
        }

        public TipJacobians CalculateJacobiansAt(PolarPoint2D polarCoordinates)
        {
            return new TipJacobians(this, polarCoordinates);
        }

        public double[] TransformScalarFieldDerivativesGlobalCartesianToLocalCartesian(double[] gradient)

        {
            return (VectorMKL.CreateFromArray(gradient, false) * TransposeRotationMatrixGlobalToLocal).CopyToArray();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="gradient">A 2x2 matrix, for which: Row i is the gradient of the ith component of the vector 
        ///     field, thus:    gradient = [Fx,x Fx,y; Fy,x Fy,y],
        ///     where Fi,j is the derivative of component i w.r.t. coordinate j</param>
        /// <returns></returns>
        public Matrix TransformVectorFieldDerivativesGlobalCartesianToLocalCartesian(Matrix gradient)
        {
            return RotationMatrixGlobalToLocal * (gradient * TransposeRotationMatrixGlobalToLocal);
        }

        public Tensor2D TransformTensorGlobalCartesianToLocalCartesian(Tensor2D tensor)
        {
            return tensor.Rotate(RotationAngle);
        }
    }
}
