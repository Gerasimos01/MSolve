using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.SamplesConsole
{
    class ProgramRotateCheck
    {
        

        private static void Example()
        {

                  
            int paradeigma = 1;
            if (paradeigma == 1)
            {
                //tVn = new double[] { 0, 0, 1 };
                //tV1 = new double[] { 1, 0, 0 };
                //tV2 = new double[] { 0, 1, 0 };
                calculation calculation1 = new calculation
                {
                    tVn = new double[] { 0, 0, 1 },
                    tV1 = new double[] { 1, 0, 0 },
                    tV2 = new double[] { 0, 1, 0 }
                };
                for (int j = 0; j < 180; j++)
                {
                      calculation1.RotateNodalDirectionVectors(2*0.01234134149, 2*0.01234134149);                    
                }
            }


        }

        
        

        //static void Main(string[] args)
        //{
        //    ProgramRotateCheck.Example();
        //}
    }

    class calculation
    {
        public double[] tVn;
        public double[] tV1;
        public double[] tV2;
        public double ak;
        public double bk;
        public int paradeigma;

        // voithitikes metavlhtes gia thn peristrofh
        private double[] tdtVn = new double[3];
        private double[] tdtV1 = new double[3];
        private double[] tdtV2 = new double[3];
        private double theta;
        private double[] theta_vec = new double[3];
        private double[,] s_k = new double[3, 3];
        private double gk1;
        private double[,] Q = new double[3, 3];
        private double[,] Q2 = new double[3, 3];
        public void RotateNodalDirectionVectors(double ak, double bk)
        {
            for (int j = 0; j < 3; j++)
            {
                theta_vec[j] = ak * tV1[j] + bk * tV2[j];
            }
            theta = Math.Sqrt((theta_vec[0] * theta_vec[0]) + (theta_vec[1] * theta_vec[1]) + (theta_vec[2] * theta_vec[2]));
            if (theta > 0)
            {
                s_k[0, 1] = -theta_vec[2];
                s_k[0, 2] = theta_vec[1];
                s_k[1, 0] = theta_vec[2];
                s_k[1, 2] = -theta_vec[0];
                s_k[2, 0] = -theta_vec[1];
                s_k[2, 1] = theta_vec[0];

                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] = (Math.Sin(theta) / theta) * s_k[j, m];
                    }
                }

                for (int m = 0; m < 3; m++)
                {
                    Q[m, m] += 1;
                }
                gk1 = 0.5 * ((Math.Sin(0.5 * theta) / (0.5 * theta)) * (Math.Sin(0.5 * theta) / (0.5 * theta)));
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q2[j, m] = 0;
                        for (int n = 0; n < 3; n++)
                        { Q2[j, m] += gk1 * s_k[j, n] * s_k[n, m]; }
                    }
                }
                for (int j = 0; j < 3; j++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Q[j, m] += Q2[j, m];
                    }
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtVn[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtVn[j] += Q[j, m] * tVn[m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tVn[j] = tdtVn[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV1[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV1[j] += Q[j, m] * tV1[m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tV1[j] = tdtV1[j];
                }
                //
                for (int j = 0; j < 3; j++)
                {
                    tdtV2[j] = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        tdtV2[j] += Q[j, m] * tV2[m];
                    }
                }

                for (int j = 0; j < 3; j++)
                {
                    tV2[j] = tdtV2[j];
                }
            }
        }
    }
}
