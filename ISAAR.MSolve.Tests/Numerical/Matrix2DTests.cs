﻿using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Numerical.Tests.LinearAlgebra
{
    public static class Matrix2DTests
    {
        [Fact]
        private static void TestLU()
        {
            var A = new Matrix2D(new double[,]
            {
                {1.0201,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,1.0100,0.0000,0.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0100,0.0000,0.0000},
                {0.0000,1.0100,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,0.0000,1.0000,0.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000,0.0000},
                {0.0000,0.0000,1.0100,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,0.0000,0.0000,1.0000,0.0000,0.0000,0.0000},
                {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,1.0000},

            });

            double [,] MultipleRHSs = new double[9, 9];
            for (int i1=0; i1<9; i1++)
            {
                for (int j1 = 0; j1 < 9; j1++)
                {
                    MultipleRHSs[i1, j1] = j1 + 1;
                }
            }

            Matrix2D solution = A.SolveLU(new Matrix2D(MultipleRHSs), true);

            double[,] expectedSolution = new double[9, 9]
            {
                {0.980296049407,1.960592098814,2.940888148221,3.921184197628,4.901480247035,5.881776296442,6.862072345848,7.842368395255,8.822664444662},
                {0.990099009901,1.980198019802,2.970297029703,3.960396039604,4.950495049505,5.940594059406,6.930693069307,7.920792079208,8.910891089109},
                {0.990099009901,1.980198019802,2.970297029703,3.960396039604,4.950495049505,5.940594059406,6.930693069307,7.920792079208,8.910891089109},
                {0.990099009901,1.980198019802,2.970297029703,3.960396039604,4.950495049505,5.940594059406,6.930693069307,7.920792079208,8.910891089109},
                {1.000000000000,2.000000000000,3.000000000000,4.000000000000,5.000000000000,6.000000000000,7.000000000000,8.000000000000,9.000000000000},
                {1.000000000000,2.000000000000,3.000000000000,4.000000000000,5.000000000000,6.000000000000,7.000000000000,8.000000000000,9.000000000000},
                {0.990099009901,1.980198019802,2.970297029703,3.960396039604,4.950495049505,5.940594059406,6.930693069307,7.920792079208,8.910891089109},
                {1.000000000000,2.000000000000,3.000000000000,4.000000000000,5.000000000000,6.000000000000,7.000000000000,8.000000000000,9.000000000000},
                {1.000000000000,2.000000000000,3.000000000000,4.000000000000,5.000000000000,6.000000000000,7.000000000000,8.000000000000,9.000000000000},
            };

            double tolerance = 1E-8;
            Assert.True(Utilities.AreEqual(new Matrix2D(expectedSolution), solution, tolerance));
        }
    }
}