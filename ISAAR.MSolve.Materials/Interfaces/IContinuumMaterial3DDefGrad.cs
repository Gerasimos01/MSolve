using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial3DDefGrad : ICloneable//,IFiniteElementMaterial
    {
        Vector Stresses { get; }
        Matrix2D ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);
        void ClearState();
        void SaveState();
        void ClearStresses();

        //copied from IfiniteElementMaterial anti comment out youngmodulus kai poissonRation
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        double[] Coordinates { get; set; }
    }
}
