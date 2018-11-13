using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;


namespace ISAAR.MSolve.Materials
{
    public class BondSlipCohMatCopy //: IFiniteElementMaterial3D // TODOGerasimos
    {
        private bool modified; // opws sto MohrCoulomb gia to modified

        public double k_elastic { get; set; } // opws sto elastic 3d 
        public double k_elastic2 { get; set; }
        public double t_max { get; set; }
        public double[] s_0 { get; set; }
        public double[] a_0 { get; set; }
        public double tol { get; set; }

        public object Clone()
        {
            return new BondSlipCohMatCopy()
            { k_elastic = this.k_elastic,
                k_elastic2 = this.k_elastic2,
                t_max = this.t_max,
                s_0 = new double[2] { this.s_0[0], this.s_0[1] },
                a_0 = new double[2] { this.a_0[0], this.a_0[1] },
                tol = this.tol
            };
        }

        private double c1;
        private double[] sigma;

        private bool matrices_not_initialized = true;
        public void InitializeMatrices()
        {
            sigma = new double[3];
            //tan_stif = new double[3, 3]; // TODOGerasimos commented out einai temporarily
            //tan_stif_prev = new double[3, 3];
            //alpha = new double[3];
            c1 = (k_elastic * k_elastic2) / (k_elastic - k_elastic2);
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
        }


    }
}
