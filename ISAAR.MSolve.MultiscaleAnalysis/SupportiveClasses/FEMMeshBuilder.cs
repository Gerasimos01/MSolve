using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    public static class FEMMeshBuilder
    {
        public static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
        {
            int elem;
            int[,] t_shell = new int[elements, 8];
            for (int nrow = 0; nrow < elem1; nrow++)
            {
                for (int nline = 0; nline < elem2; nline++)
                {
                    elem = (nrow + 1 - 1) * elem2 + nline + 1;//nrow+ 1 nline+1 einai zero based 
                    t_shell[elem - 1, -1 + 1] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 8] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 4] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;

                    t_shell[elem - 1, -1 + 5] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 2;
                    t_shell[elem - 1, -1 + 7] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 1;

                    t_shell[elem - 1, -1 + 2] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 6] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 3] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;
                }
            }
            return t_shell;

        }

        public static int Topol_rve(int h1, int h2, int h3, int hexa1, int hexa2, int hexa3, int kuvos, int endiam_plaka, int katw_plaka)
        {
            int arith;
            if (h3 == 1)
            { arith = h1 + (h2 - 1) * (hexa1 + 1) + kuvos; }
            else
            {
                if (h3 == hexa3 + 1)
                { arith = hexa3 * (hexa1 + 1) * (hexa2 + 1) + h1 + (h2 - 1) * (hexa1 + 1); }
                else
                {
                    if (h2 == 1)
                    { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + h1; }
                    else
                    {
                        if (h2 == hexa2 + 1)
                        { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + (hexa1 + 1) + 2 * (hexa2 - 1) + h1; }
                        else
                        {
                            if (h1 == 1)
                            { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 1; }
                            else
                            {
                                if (h1 == hexa1 + 1)
                                { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 2; }
                                else
                                { arith = (h1 - 1) + (h2 - 2) * (hexa1 - 1) + (h3 - 2) * (hexa1 - 1) * (hexa2 - 1); }
                            }
                        }
                    }

                }
            }
            return arith;
        }

    }
}
