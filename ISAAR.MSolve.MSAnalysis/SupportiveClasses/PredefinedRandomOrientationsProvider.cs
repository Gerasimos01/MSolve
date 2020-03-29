using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.MSAnalysis.SupportiveClasses
{
    public static class PredefinedRandomOrientationsProvider
    {
        public static int ExampleNo;

        private static Dictionary<int, (double[], double[], double[][])> rot_phi_1_rot_phi_2_ekk_xyzData = new Dictionary<int, (double[], double[], double[][])>()
        {
            {42,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[2][]{new double [] {0,0,0 }, new double[]{0,0,0} }) },
            {43,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[2][]{new double [] {0,0,0 }, new double[]{0,0,0} }) },

        };
        public static Tuple<double[], double[], double[][]> GetPredefinedOrienntationData()
        {
            return new Tuple<double[], double[], double[][]>(rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item1,rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item2,
                rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item3 );
        }
    }
}
