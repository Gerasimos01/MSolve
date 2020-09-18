﻿using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Tests.TestData;
using ISAAR.MSolve.LinearAlgebra.Tests.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using Xunit;

namespace ISAAR.MSolve.LinearAlgebra.Tests.Matrices
{
    /// <summary>
    /// Tests for <see cref="Matrix"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class TensorsTests1
    {
        [Fact]
        private static void TestContraction()
        {
            // gia to paradeigma i vash einai monadiaia kai ta dianusmata ths einai orthogonika metaxi tous
            var Ei = new double[3, 3] { { 0.750000000000000, -0.500000000000000, -0.433012701892219 }, { 0.433012701892219, 0.866025403784439, -0.250000000000000 }, { 0.500000000000000, 0, 0.866025403784439, } };
            var ei = new double[3, 3] { { 0.709406479916223 ,- 0.555128922733495, - 0.434251454111122 }, { 0.496731764892154,    0.830902789308865, - 0.250715193927043 }, { 0.500000000000000, - 0.0378475080223054, 0.865197992448261 } };

            var ci = new double[3, 3]; ci[0, 0] = 1; ci[1, 1] = 1; ci[2, 2] = 1;

            tensorOrder2 tensor1 = new tensorOrder2() { basis1 = Ei, basis2 = ei, coefficients = new double[,] { { 10, 20, 30 }, { 40, 50, 60 }, { 70, 80, 90 } } };
            tensorOrder2 tensor2 = new tensorOrder2() { basis1 = ei, basis2 = Ei, coefficients = new double[,] { { 10, 20, 30 }, { 40, 50, 60 }, { 70, 80, 90 } } };

            var tensor3 = tensor1.SingleContract(tensor2);
            var tensor4 = tensor3.ProjectIn3DCartesianBasis();

            var tensor5 = tensor1.ProjectIn3DCartesianBasis();
            var tensor6 = tensor2.ProjectIn3DCartesianBasis();

            var multiplicationData = (Matrix.CreateFromArray(tensor5.coefficients)* Matrix.CreateFromArray(tensor6.coefficients)).CopyToArray2D();

            var tensor7 = tensor5.SingleContract(tensor6);

            //comparer.AssertEqual()
            comparer.AssertEqual(multiplicationData, tensor4.coefficients);

            comparer.AssertEqual(multiplicationData, tensor7.coefficients);

        }

        [Fact]
        private static void TestScale()
        {
            // gia to paradeigma i vash einai monadiaia kai ta dianusmata ths einai orthogonika metaxi tous
            var ci = new double[3, 3]; ci[0, 0] = 1; ci[1, 1] = 1; ci[2, 2] = 1;

            tensorOrder2 tensor1 = new tensorOrder2(new double[,] { { 10, 20, 30 }, { 40, 50, 60 }, { 70, 80, 90 } } );
            tensorOrder2 tensor2 = new tensorOrder2(new double[,] { { 10, 20, 30 }, { 40, 50, 60 }, { 70, 80, 90 } });

            var tensor3 = tensor2.Scale(2);

            var tensor4 = tensor1.AddTensor(tensor3);


            comparer.AssertEqual(Matrix.CreateFromArray(tensor1.coefficients).Scale(3).CopyToArray2D(), tensor4.coefficients);

        }




        private static readonly MatrixComparer comparer = new MatrixComparer(1E-13);

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        internal static void TestAddition(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var expected = Matrix.CreateFromArray(
                    MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, 1.0, SymmPosDef10by10.Matrix));

                // operator+
                comparer.AssertEqual(expected, A1 + A2);
            });
        }

        [Fact]
        private static void TestClear()
        {
            var zero = Matrix.CreateZero(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols);
            var matrix = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            matrix.Clear();
            comparer.AssertEqual(zero, matrix);
        }

        [Fact]
        private static void TestEquality()
        {
            // Equals(SkylineMatrix)
            var full1 = Matrix.CreateFromArray(SparsePosDef10by10.Matrix);
            var skyline1 = SkylineMatrix.CreateFromArrays(SparsePosDef10by10.Order,
                SparsePosDef10by10.SkylineValues, SparsePosDef10by10.SkylineDiagOffsets, true, true);
            Assert.True(full1.Equals(skyline1));

            // Equals(CsrMatrix)
            var full2 = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csr2 = CsrMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CsrValues, SparseRectangular10by5.CsrColIndices, SparseRectangular10by5.CsrRowOffsets,
                true);
            Assert.True(full2.Equals(csr2));

            // Equals(CscMatrix)
            var full3 = Matrix.CreateFromArray(SparseRectangular10by5.Matrix);
            var csc3 = CscMatrix.CreateFromArrays(SparseRectangular10by5.NumRows, SparseRectangular10by5.NumCols,
                SparseRectangular10by5.CscValues, SparseRectangular10by5.CscRowIndices, SparseRectangular10by5.CscColOffsets,
                true);
            Assert.True(full3.Equals(csc3));
        }

        [Fact]
        private static void TestGetColumn()
        {
            var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
            for (int j = 0; j < RectangularFullRank10by5.NumCols; ++j)
            {
                Vector colExpected = DenseStrategies.GetColumn(matrix, j);
                Vector colComputed = matrix.GetColumn(j);
                comparer.AssertEqual(colExpected, colComputed);
            }
        }

        [Fact]
        private static void TestGetRow()
        {
            var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
            for (int i = 0; i < RectangularFullRank10by5.NumRows; ++i)
            {
                Vector rowExpected = DenseStrategies.GetRow(matrix, i);
                Vector rowComputed = matrix.GetRow(i);
                comparer.AssertEqual(rowExpected, rowComputed);
            }
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestLinearCombination(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                double scalar1 = 2.0;
                var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                double scalar2 = 3.5;
                var expected = Matrix.CreateFromArray(
                    MatrixOperations.LinearCombination(scalar1, SquareSingular10by10.Matrix, scalar2, SymmPosDef10by10.Matrix));

                // LinearCombination()
                comparer.AssertEqual(expected, A1.LinearCombination(scalar1, A2, scalar2));

                // LinearCombinationIntoThis()
                Matrix temp = A1.Copy();
                temp.LinearCombinationIntoThis(scalar1, A2, scalar2);
                comparer.AssertEqual(expected, temp);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixMatrixMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                var expectedA1TimesA2 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(SquareSingular10by10.Matrix, RectangularFullRank10by5.Matrix));
                var expectedTransposeA2TimesA1 = Matrix.CreateFromArray(
                    MatrixOperations.MatrixTimesMatrix(
                        MatrixOperations.Transpose(RectangularFullRank10by5.Matrix), SquareSingular10by10.Matrix));

                // MultiplyRight() without transposition
                comparer.AssertEqual(expectedA1TimesA2, A1.MultiplyRight(A2, false, false));

                // operator*
                comparer.AssertEqual(expectedA1TimesA2, A1 * A2);

                // MultiplyRight() with transposition
                comparer.AssertEqual(expectedTransposeA2TimesA1, A2.MultiplyRight(A1, true, false));

                // MultiplyRight() with incorrect dimensions
                Assert.Throws<NonMatchingDimensionsException>(() => A2.MultiplyRight(A1, false, false));
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplication(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // rectangular 10-by-5
                var A1 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                var x1 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs5);
                var b1Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs10);
                Vector b1Computed = A1.Multiply(x1, false);
                comparer.AssertEqual(b1Expected, b1Computed);

                // rectangular 5-by-10
                double[,] fullRank5by10 = MatrixOperations.Transpose(RectangularFullRank10by5.Matrix);
                var x2 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs10);
                var b2Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs5);
                Vector b2Computed = A1.Multiply(x2, true);
                comparer.AssertEqual(b2Expected, b2Computed);

                // square invertible 10-by-10
                var A3 = Matrix.CreateFromArray(SquareInvertible10by10.Matrix);
                var x3 = Vector.CreateFromArray(SquareInvertible10by10.Lhs);
                var b3Expected = Vector.CreateFromArray(SquareInvertible10by10.Rhs);
                Vector b3Computed = A3.Multiply(x3, false);
                comparer.AssertEqual(b3Expected, b3Computed);

                // square singular 10-by-10 (rank = 8)
                var A4 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var x4 = Vector.CreateFromArray(SquareSingular10by10.Lhs);
                var b4Expected = Vector.CreateFromArray(SquareSingular10by10.Rhs);
                Vector b4Computed = A4.Multiply(x4, false);
                comparer.AssertEqual(b4Expected, b4Computed);

                // square singular 10-by-10 (rank = 9)
                var A5 = Matrix.CreateFromArray(SquareSingularSingleDeficiency10by10.Matrix);
                var x5 = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.Lhs);
                var b5Expected = Vector.CreateFromArray(SquareSingularSingleDeficiency10by10.Rhs);
                Vector b5Computed = A5.Multiply(x5, false);
                comparer.AssertEqual(b5Expected, b5Computed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestMatrixVectorMultiplicationIntoResult(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                // The result vectors will first be set to some non zero values to make sure that the result overwrites 
                // them instead of being added to them.

                // MultiplyIntoResult() - untransposed 
                var A1 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                var x1 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs5);
                var b1Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs10);
                Vector b1Computed = Vector.CreateWithValue(A1.NumRows, 1.0);
                A1.MultiplyIntoResult(x1, b1Computed, false);
                comparer.AssertEqual(b1Expected, b1Computed);

                // MultiplyIntoResult() - transposed
                var x2 = Vector.CreateFromArray(RectangularFullRank10by5.Lhs10);
                var b2Expected = Vector.CreateFromArray(RectangularFullRank10by5.Rhs5);
                Vector b2Computed = Vector.CreateWithValue(A1.NumColumns, 1.0);
                A1.MultiplyIntoResult(x2, b2Computed, true);
                comparer.AssertEqual(b2Expected, b2Computed);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestScaling(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var matrix = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
                double scalar = 5.0;
                var expected = Matrix.CreateFromArray(MatrixOperations.Scale(scalar, RectangularFullRank10by5.Matrix));

                // Scale()
                comparer.AssertEqual(expected, matrix.Scale(scalar));

                // ScaleIntoThis()
                Matrix temp = matrix.Copy();
                temp.ScaleIntoThis(scalar);
                comparer.AssertEqual(expected, temp);

                // operator*
                comparer.AssertEqual(expected, scalar * matrix);
            });
        }

        [Theory]
        [MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
        private static void TestSubtraction(LinearAlgebraProviderChoice providers)
        {
            TestSettings.RunMultiproviderTest(providers, delegate ()
            {
                var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
                var A2 = Matrix.CreateFromArray(SymmPosDef10by10.Matrix);
                var expected = Matrix.CreateFromArray(
                    MatrixOperations.LinearCombination(1.0, SquareSingular10by10.Matrix, -1.0, SymmPosDef10by10.Matrix));

                // operator+
                comparer.AssertEqual(expected, A1 - A2);
            });
        }

        [Fact]
        private static void TestTransposition()
        {
            // square
            var A1 = Matrix.CreateFromArray(SquareSingular10by10.Matrix);
            var A1TransposeExpected = MatrixOperations.Transpose(SquareSingular10by10.Matrix);
            Matrix A1TransposeComputed = A1.Transpose();
            comparer.AssertEqual(A1TransposeExpected, A1TransposeComputed.CopyToArray2D());

            // rectangular
            var A2 = Matrix.CreateFromArray(RectangularFullRank10by5.Matrix);
            var A2TransposeExpected = MatrixOperations.Transpose(RectangularFullRank10by5.Matrix);
            Matrix A2TransposeComputed = A2.Transpose();
            comparer.AssertEqual(A2TransposeExpected, A2TransposeComputed.CopyToArray2D());
        }
    }
}
