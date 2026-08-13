// <copyright file="NumericMatrix.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Collections;
using System.Numerics;

namespace NumericsLab.Internals;

internal sealed class NumericMatrix<T>
    where T : INumber<T>
{
    private readonly T[,] _matrix;

    public NumericMatrix(int numberOfRows, int numberOfColumns)
    {
        _matrix = new T[numberOfRows, numberOfColumns];
        NumberOfRows = numberOfRows;
        NumberOfColumns = numberOfColumns;
    }

    public int NumberOfRows { get; }

    public int NumberOfColumns { get; }

    /// <summary>
    /// Gets or sets the value of the element at the specified one-based row and column indices in the matrix.
    /// </summary>
    /// <param name="row">The one-based row index of the element.</param>
    /// <param name="column">The one-based column index of the element.</param>
    /// <returns>The value of the element at the specified one-based row and column indices.</returns>
    public T this[int row, int column]
    {
        get => _matrix[row - 1, column - 1];
        set => _matrix[row - 1, column - 1] = value;
    }

    internal INumericVector<T> GetColumnVector(int j, int startRowIndex = 1) => new ColumnVector(this, j, startRowIndex);

    private sealed class ColumnVector(NumericMatrix<T> numericMatrix, int columnIndex, int startRowIndex) : INumericVector<T>
    {
        public int Length => numericMatrix.NumberOfRows - startRowIndex + 1;

        public T this[int index]
        {
            get => numericMatrix._matrix[index + startRowIndex - 2, columnIndex - 1];
            set => numericMatrix._matrix[index + startRowIndex - 2, columnIndex - 1] = value;
        }

        public IReadOnlyList<T> AsReadOnlyList() => new ColumnList(this);

        private readonly struct ColumnList(NumericMatrix<T>.ColumnVector columnVector) : IReadOnlyList<T>
        {
            public int Count => columnVector.Length;

            public readonly T this[int index] => columnVector[index + 1];

            public IEnumerator<T> GetEnumerator() => new ColumnEnumerator(columnVector);

            IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

            private struct ColumnEnumerator(ColumnVector columnVector) : IEnumerator<T>
            {
                private int _currentIndex = 1;

                public readonly T Current => columnVector[_currentIndex];

#pragma warning disable IDE0251 // Make member 'readonly'
                object IEnumerator.Current => Current;
#pragma warning restore IDE0251 // Make member 'readonly'

                public bool MoveNext()
                {
                    if (_currentIndex < columnVector.Length)
                    {
                        _currentIndex++;
                        return true;
                    }

                    return false;
                }

                public void Reset() => _currentIndex = 1;

#pragma warning disable IDE0251 // Make member 'readonly'
                public void Dispose()
                {
                }
#pragma warning restore IDE0251 // Make member 'readonly'
            }
        }
    }
}