using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.IGA.Tests
{
    public static class SupportiveMethods
    {
        public static double[] modify_ox_sunol_forRotationAndTranslation(double[] o_xsunol, double rot_phi_1, double rot_phi_2, double[] ekk_xyz)
        {

            // rotation data
            double e1_new_z = Math.Sin(rot_phi_2);
            double e1_new_y = Math.Sin(rot_phi_1) * Math.Cos(rot_phi_2); // e1_new_xy = cos(rot_phi_2);
            double e1_new_x = Math.Cos(rot_phi_1) * Math.Cos(rot_phi_2);

            double e2_new_y = Math.Cos(rot_phi_1);
            double e2_new_x = -Math.Sin(rot_phi_1);
            double e2_new_z = 0;

            double[,] e_new = new double[3, 3] { { e1_new_x, e2_new_x, 0 }, { e1_new_y, e2_new_y, 0 }, { e1_new_z, e2_new_z, 0 } };
            double[] e_new_cross = new double[3];
            cross(new double[3] { e1_new_x, e1_new_y, e1_new_z }, new double[3] { e2_new_x, e2_new_y, e2_new_z }, e_new_cross);
            for (int i1 = 0; i1 < 3; i1++)
            { e_new[i1, 2] = e_new_cross[i1]; }

            double[,] e_old = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

            double[,] Qij = new double[3, 3];
            for (int q1 = 0; q1 < 3; q1++)
            {
                for (int q2 = 0; q2 < 3; q2++)
                {
                    Qij[q1, q2] = dot_product(new double[3] { e_old[0, q1], e_old[1, q1], e_old[2, q1] }, new double[3] { e_new[0, q2], e_new[1, q2], e_new[2, q2] });
                }
            }

            // dhmiourgia neou dianusmatos me ROTATION tou initial
            double[] o_xsunol_new = new double[o_xsunol.GetLength(0)];

            for (int komvos = 0; komvos < o_xsunol.GetLength(0) / 6; komvos++)
            {
                double[] product;
                product = MatVecMult(Qij, new double[3] { o_xsunol[6 * (komvos) + 0], o_xsunol[6 * (komvos) + 1], o_xsunol[6 * (komvos) + 2] });
                for (int q2 = 0; q2 < 3; q2++) { o_xsunol_new[6 * (komvos) + q2] = product[q2]; }
                product = MatVecMult(Qij, new double[3] { o_xsunol[6 * (komvos) + 3], o_xsunol[6 * (komvos) + 4], o_xsunol[6 * (komvos) + 5] });
                for (int q2 = 0; q2 < 3; q2++) { o_xsunol_new[6 * (komvos) + 3 + q2] = product[q2]; }
            }

            // TRANSLATION dlh h ekkentrothta
            for (int komvos = 0; komvos < o_xsunol.GetLength(0) / 6; komvos++)
            {
                for (int q2 = 0; q2 < 3; q2++) { o_xsunol_new[6 * (komvos) + q2] = o_xsunol_new[6 * (komvos) + q2] + ekk_xyz[q2]; }
            }

            return o_xsunol_new;
        }

        public static double dot_product(double[] vec1, double[] vec2)
        {
            double dot_product = 0;
            for (int i1 = 0; i1 < vec1.GetLength(0); i1++)
            {
                dot_product += vec1[i1] * vec2[i1];
            }
            return dot_product;
        }

        public static double[] MatVecMult(double[,] matrix, double[] vec)
        {
            double[] product = new double[matrix.GetLength(0)];
            for (int q1 = 0; q1 < matrix.GetLength(0); q1++)
            {
                for (int q2 = 0; q2 < matrix.GetLength(1); q2++)
                {
                    product[q1] += matrix[q1, q2] * vec[q2];
                }
            }
            return product;
        }

        public static void cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }
    }
}
