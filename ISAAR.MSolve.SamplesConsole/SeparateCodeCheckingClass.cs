using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;



namespace ISAAR.MSolve.SamplesConsole
{
    class SeparateCodeCheckingClass
    {
        public static void Check01()
        {
            double[,] Aijkl = new double[9, 9];
            for (int k = 0; k < 9; k++)
            {
                for (int l = 0; l < 9; l++)
                {
                    Aijkl[k, l] = 100 * (k+1) + (l+1);
                }
            }

            double[,] SPK = new double[3, 3] { { 0.5, -0.1, -0.0 },{ -0.1, 0.2, 0.1 },{ 0, 0.1, 0.15 } } ;
            double[,] F = new double[3, 3] { { 1.2, -0.1, -0.0 }, { -0.1,1.1, 0.1 }, { 0, 0.1, 1.15 } };

            double[,] Cinpk = Transform_d2Wdfdf_to_Cijrs(Aijkl, SPK, F);
        }

        public static void Check02()
        {
            double[,] F__F__ = new double[9, 9];
            

            for (int i1 = 0; i1 < 9; i1++)
            {
                F__F__[i1, i1] = 1;
            }

            Vector solution = new Vector(new double[9]);

            SkylineMatrix2D F__F__Mat = new SkylineMatrix2D(F__F__);
            int linearsystemID = 1;
            SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[9]);
            var solver = new SolverSkyline(linearSystem);
            // BuildMatrices();
            linearSystem.Matrix = F__F__Mat;
            //solver.Initialize();
            solver.Initialize(); // dld factorize

            double[] rhs = new double[9];
            for (int l = 0; l < 9; l++)
            {
                rhs[l] = 10;
            }
            Vector RHS = new Vector(rhs);
            SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
            k.Solve(RHS, solution);

        }

        public static void Check03()
        {
            Dictionary<int, double> dictionary1 = new Dictionary<int, double>();
            dictionary1.Add(4, 0.111112);
            dictionary1.Add(8, 0.222223);
            var dok1 = dictionary1.Keys;
            var dok2 = dictionary1.GetEnumerator();
        }

        private  static double[,] Transform_d2Wdfdf_to_Cijrs(double[,] Aijkl, double[,] SPK, double[,] F)
        {
            int[,] i_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };
            int[,] k_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };

            double[,] Cinpk = new double[9, 9];

            double[,] F__F__ = new double[9, 9];
            //[F(:, 1) * F(1,:), F(:, 2) * F(1,:), F(:, 3) * F(1,:);
            //F(:, 1) * F(2,:),F(:, 2) * F(2,:),F(:, 3) * F(2,:);
            //F(:, 1) * F(3,:),F(:, 2) * F(3,:),F(:, 3) * F(3,:)];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            F__F__[3 * i1 + k, 3 * j1 + l] = F[k, j1] * F[i1, l];
                        }
                    }

                }
            }

            //TODO upologismos F_inv (F__F__^(-1))

            double[,] F_inv = new double[9, 9];
            double[,] eye9 = new double[9, 9];
            for (int i1 = 0; i1 < 9; i1++)
            { eye9[i1, i1] = 1; }
            Vector solution = new Vector(new double[9]);

            SkylineMatrix2D F__F__Mat = new SkylineMatrix2D(F__F__);
            F__F__Mat.WriteToFile(@"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\F__F__.txt");
            int linearsystemID = 1;
            SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[9]);
            var solver = new SolverSkyline(linearSystem);
            // BuildMatrices();
            linearSystem.Matrix = F__F__Mat;
            //solver.Initialize();
            solver.Initialize(); // dld factorize

            for (int j1 = 0; j1 < 9; j1++)
            {
                Vector RHS = new Vector(new double[9] { eye9[0, j1], eye9[1, j1], eye9[2, j1], eye9[3, j1], eye9[4, j1], eye9[5, j1], eye9[6, j1], eye9[7, j1], eye9[8, j1] });
                SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
                k.Solve(RHS, solution);
                for (int i1 = 0; i1 < 9; i1++)
                {
                    F_inv[i1, j1] = solution[i1];
                }
            }


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    double[] A_j_l = new double[9] { Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1]};

                    double[] sec_term = new double[9] { -SPK[i1 - 1, k1 - 1], 0, 0, 0, -SPK[i1 - 1, k1 - 1], 0, 0, 0, -SPK[i1 - 1, k1 - 1] };

                    Matrix2D F_invMat = new Matrix2D(F_inv);
                    Vector A_j_lVec = new Vector(A_j_l);
                    Vector sec_termVec = new Vector(sec_term);
                    Vector C_np_ = F_invMat * (new Vector(A_j_lVec + sec_termVec));

                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[0];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[1];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[2];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[3];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[4];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[5];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[6];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[7];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[8];

                }
            }

            return Cinpk;
        }



    }
}
