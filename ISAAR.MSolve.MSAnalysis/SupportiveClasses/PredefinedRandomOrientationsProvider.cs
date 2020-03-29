using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.MSAnalysis.SupportiveClasses
{
    public static class PredefinedRandomOrientationsProvider
    {
        public static int ExampleNo;

        public static double L01 = 120;

        public static double elem_len = L01 / 18;

        public static double a= 38;
        public static double b=L01/3;

        public static double term_1 = -0.5*a;

        public static double s_term = 1.5 * b / 18;

        public static double term = -0.5 * b;

        public static double x_mikro = term_1 + s_term + term;

        public static double x_megalo = term;

        public static double y_m = -b;

        private static Dictionary<int, (double[], double[], double[][])> rot_phi_1_rot_phi_2_ekk_xyzData = new Dictionary<int, (double[], double[], double[][])>()
        {
            {42,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_mikro,y_m,4.5*elem_len }, 
                    new double [] {x_mikro,y_m,5.5*elem_len }, 
                    new double [] {x_mikro,y_m,6.5*elem_len },
                    new double [] {x_mikro,y_m,7.5*elem_len },
                    new double [] {x_mikro,y_m,8.5*elem_len },})},

            {43,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_megalo,y_m,4.5*elem_len },
                    new double [] {x_megalo,y_m,5.5*elem_len },
                    new double [] {x_megalo,y_m,6.5*elem_len },
                    new double [] {x_megalo,y_m,7.5*elem_len },
                    new double [] { x_megalo, y_m,8.5*elem_len },})},

            {44,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_mikro,y_m,4.5*elem_len } }) },

            {45,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_megalo,y_m,4.5*elem_len } }) },

            {46,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {y_m,y_m,4.5*elem_len },
                    new double [] {y_m,y_m,5.5*elem_len },
                    new double [] {y_m,y_m,6.5*elem_len },
                    new double [] {y_m,y_m,7.5*elem_len },
                    new double [] {y_m,y_m,8.5*elem_len },})},

            {47,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_mikro, x_mikro, 4.5 *elem_len },
                    new double [] {x_mikro,x_mikro,5.5*elem_len },
                    new double [] {x_mikro,x_mikro,6.5*elem_len },
                    new double [] {x_mikro,x_mikro,7.5*elem_len },
                    new double [] {x_mikro, x_mikro, 8.5 *elem_len },})},

            {48,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_megalo,x_megalo,4.5*elem_len },
                    new double [] {x_megalo,x_megalo,5.5*elem_len },
                    new double [] {x_megalo,x_megalo,6.5*elem_len },
                    new double [] {x_megalo, x_megalo, 7.5 *elem_len },
                    new double [] { x_megalo, x_megalo, 8.5 *elem_len },})},

            {49,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_mikro, x_mikro, 4.5 *elem_len } }) },

            {50,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_megalo, x_megalo, 4.5 *elem_len } }) },


            //
            //
            //Z_ekkentrothtes

            {51,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_mikro,y_m,4.5*elem_len-b },
                    new double [] {x_mikro,y_m,5.5*elem_len-b },
                    new double [] {x_mikro,y_m,6.5*elem_len-b },
                    new double [] {x_mikro,y_m,7.5*elem_len-b },
                    new double [] {x_mikro,y_m,8.5*elem_len-b },})},

            {52,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_megalo,y_m,4.5*elem_len-b },
                    new double [] {x_megalo,y_m,5.5*elem_len-b },
                    new double [] {x_megalo,y_m,6.5*elem_len-b },
                    new double [] {x_megalo,y_m,7.5*elem_len-b },
                    new double [] { x_megalo, y_m,8.5*elem_len-b },})},

            {53,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_mikro,y_m,4.5*elem_len-b } }) },

            {54,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_megalo,y_m,4.5*elem_len-b } }) },

            {55,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {y_m,y_m,4.5*elem_len-b },
                    new double [] {y_m,y_m,5.5*elem_len-b },
                    new double [] {y_m,y_m,6.5*elem_len-b },
                    new double [] {y_m,y_m,7.5*elem_len-b },
                    new double [] {y_m,y_m,8.5*elem_len-b },})},

            {56,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_mikro, x_mikro, 4.5 *elem_len-b },
                    new double [] {x_mikro,x_mikro,5.5*elem_len-b },
                    new double [] {x_mikro,x_mikro,6.5*elem_len-b },
                    new double [] {x_mikro,x_mikro,7.5*elem_len-b },
                    new double [] {x_mikro, x_mikro, 8.5 *elem_len-b },})},

            {57,(first: new double[]{0,0,0,0,0 }, second: new double[]{0,0,0,0,0 },
                third: new double[5][]{new double [] {x_megalo,x_megalo,4.5*elem_len-b },
                    new double [] {x_megalo,x_megalo,5.5*elem_len-b },
                    new double [] {x_megalo,x_megalo,6.5*elem_len-b },
                    new double [] {x_megalo, x_megalo, 7.5 *elem_len-b },
                    new double [] { x_megalo, x_megalo, 8.5 *elem_len-b },})},

            {58,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_mikro, x_mikro, 4.5 *elem_len-b } }) },

            {59,(first: new double[]{0,0 }, second: new double[]{0,0 }, third: new double[1][]{new double [] {x_megalo, x_megalo, 4.5 *elem_len-b } }) },
        };
        public static Tuple<double[], double[], double[][]> GetPredefinedOrienntationData()
        {
            return new Tuple<double[], double[], double[][]>(rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item1,rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item2,
                rot_phi_1_rot_phi_2_ekk_xyzData[ExampleNo].Item3 );
        }
    }
}
