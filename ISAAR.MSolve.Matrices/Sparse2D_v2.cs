using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISAAR.MSolve.Matrices
{
    public class Sparse2D_v2<T> : IMatrix2D<T>
    {
        private int rows, columns;
        private readonly Dictionary<int, Dictionary<int, int>> rowColDataPositions;
        private readonly List<T> data;

        public Sparse2D_v2(int rows, int columns)
        {
            this.rows = rows;
            this.columns = columns;
            rowColDataPositions = new Dictionary<int, Dictionary<int, int>>(rows);
            data = new List<T>(rows);
        }

        private Sparse2D_v2(int rows, int columns, List<T> data, Dictionary<int, Dictionary<int, int>> rowColDataPositions)
        {
            this.rows = rows;
            this.columns = columns;
            this.rowColDataPositions = rowColDataPositions;
        }

        private T GetValueFromRowCol(int row, int col)
        {
            T value = default(T);
            if (rowColDataPositions.ContainsKey(row))
                if (rowColDataPositions[row].ContainsKey(col)) value = data[rowColDataPositions[row][col]];
            return value;
        }

        private void SetValueAtRowCol(int row, int col, T value)
        {
            int pos = data.Count;
            if (rowColDataPositions.ContainsKey(row))
            {
                if (rowColDataPositions[row].ContainsKey(col))
                {
                    data[rowColDataPositions[row][col]] = value;
                    return;
                }
                else
                    rowColDataPositions[row].Add(col, pos);
            }
            else
            {
                rowColDataPositions.Add(row, new Dictionary<int, int>());
                rowColDataPositions[row].Add(col, pos);
            }
            data.Add(value);
        }

        public double[] MultiplyRight(double[] vector)
        {
            if (Rows != vector.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");
            List<double> d = data as List<double>;
            double[] result = new double[vector.Length];
            foreach (int row in rowColDataPositions.Keys)
            {
                double sum = 0.0;
                foreach (int col in rowColDataPositions[row].Keys)
                    sum += d[rowColDataPositions[row][col]] * vector[col];
                result[row] = sum;
            }
            return result;
        }

        public Matrix2D<double> MultiplyRight(IMatrix2D<double> matrix)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            if (this.Columns != matrix.Rows) throw new InvalidOperationException("Matrix sizes mismatch.");

            double[,] c = new double[this.Rows, matrix.Columns];
            var AA = new Sparse2D_v2<double>(rows, columns, this.data as List<double>, rowColDataPositions);

            for (int i = 0; i < this.Rows; i++)
                for (int k = 0; k < matrix.Columns; k++)
                    for (int j = 0; j < matrix.Rows; j++)
                        c[i, k] += AA[i, j] * matrix[j, k];

            return new Matrix2D<double>(c);
        }

        public Matrix2D<double> MultiplyLeft(IMatrix2D<double> matrix)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            if (this.Rows != matrix.Columns) throw new InvalidOperationException("Matrix sizes mismatch.");

            double[,] c = new double[matrix.Rows, this.Columns];
            var AA = new Sparse2D_v2<double>(rows, columns, this.data as List<double>, rowColDataPositions);

            for (int i = 0; i < matrix.Rows; i++)
                for (int k = 0; k < this.Columns; k++)
                    for (int j = 0; j < matrix.Columns; j++)
                        c[i, k] += matrix[i, j] * this[j, k];

            return new Matrix2D<double>(c);
        }

        public Sparse2D_v2<T> Transpose()
        {
            Sparse2D_v2<T> transpose = new Sparse2D_v2<T>(rows, columns);
            foreach (int row in rowColDataPositions.Keys)
                foreach (int col in rowColDataPositions[row].Keys)
                    transpose[col, row] = data[rowColDataPositions[row][col]];
            return transpose;
        }


        #region IMatrix2D<double> Members

        public int Rows
        {
            get { return rows; }
        }

        public int Columns
        {
            get { return columns; }
        }

        public T this[int x, int y]
        {
            get { return GetValueFromRowCol(x, y); }
            set { SetValueAtRowCol(x, y, value); }
        }

        public void Scale(double scale)
        {
            throw new NotImplementedException();
        }

        public void Multiply(IVector<double> vIn, double[] vOut)
        {
            if (!(typeof(T) == typeof(double))) throw new InvalidOperationException("Cannot multiply for types other than double");
            if (Rows != vIn.Length) throw new InvalidOperationException("Matrix and vector size mismatch.");
            Array.Clear(vOut, 0, vOut.Length);
            List<double> d = data as List<double>;

            foreach (int row in rowColDataPositions.Keys)
                foreach (int col in rowColDataPositions[row].Keys)
                    vOut[row] += d[rowColDataPositions[row][col]] * vIn[col];
        }

        public void LinearCombination(IList<T> coefficients, IList<IMatrix2D<T>> matrices)
        {
            throw new NotImplementedException();
        }

        public void Solve(IVector<double> f, double[] result)
        {
            throw new NotImplementedException();
        }

        public void WriteToFile(string name)
        {
            throw new NotImplementedException();
        }

        public void ReadFromFile(string name)
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
