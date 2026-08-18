// <copyright file="NumericMatrix.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Collections;
using System.Numerics;

namespace NumericsLab.Internals;

internal sealed class NumericMatrix<T>(T[,] matrix)
    where T : INumber<T>
{
    private readonly T[,] _matrix = matrix;

    public NumericMatrix(int numberOfRows, int numberOfColumns)
        : this(new T[numberOfRows, numberOfColumns])
    {
    }

    public int NumberOfRows { get; } = matrix.GetLength(0);

    public int NumberOfColumns { get; } = matrix.GetLength(1);

    /// <summary>
    /// Gets or sets the value of the element at the specified one-based row and column indices in the matrix.
    /// </summary>
    /// <param name="oneBasedRowIndex">The one-based row index of the element.</param>
    /// <param name="oneBasedColumnIndex">The one-based column index of the element.</param>
    /// <returns>The value of the element at the specified one-based row and column indices.</returns>
    public T this[int oneBasedRowIndex, int oneBasedColumnIndex]
    {
        get => _matrix[oneBasedRowIndex - 1, oneBasedColumnIndex - 1];
        set => _matrix[oneBasedRowIndex - 1, oneBasedColumnIndex - 1] = value;
    }

    internal INumericVector<T> GetColumnVector(int oneBasedColumnIndex, int oneBasedStartRowIndex = 1) => new ColumnVector(this, oneBasedColumnIndex, oneBasedStartRowIndex);

    private sealed class ColumnVector(NumericMatrix<T> numericMatrix, int oneBasedColumnIndex, int oneBasedStartRowIndex) : INumericVector<T>
    {
        public int Length => numericMatrix.NumberOfRows - oneBasedStartRowIndex + 1;

        public T this[int oneBasedIndex]
        {
            get => numericMatrix._matrix[oneBasedIndex + oneBasedStartRowIndex - 2, oneBasedColumnIndex - 1];
            set => numericMatrix._matrix[oneBasedIndex + oneBasedStartRowIndex - 2, oneBasedColumnIndex - 1] = value;
        }

        public IReadOnlyList<T> AsReadOnlyList() => new ReadOnlyColumnList(this);

        private readonly struct ReadOnlyColumnList(ColumnVector columnVector) : IReadOnlyList<T>
        {
            public int Count => columnVector.Length;

            public readonly T this[int zeroBasedIndex] => columnVector[zeroBasedIndex + 1];

            public IEnumerator<T> GetEnumerator() => new ColumnEnumerator(columnVector);

            IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

            private struct ColumnEnumerator(ColumnVector columnVector) : IEnumerator<T>
            {
                private int _currentZeroBasedIndex = 0;

                public readonly T Current => columnVector[_currentZeroBasedIndex];

#pragma warning disable IDE0251 // Make member 'readonly'
                object IEnumerator.Current => Current;
#pragma warning restore IDE0251 // Make member 'readonly'

                public bool MoveNext()
                {
                    if (_currentZeroBasedIndex < columnVector.Length)
                    {
                        _currentZeroBasedIndex++;
                        return true;
                    }

                    return false;
                }

                public void Reset() => _currentZeroBasedIndex = 0;

#pragma warning disable IDE0251 // Make member 'readonly'
                public void Dispose()
                {
                }
#pragma warning restore IDE0251 // Make member 'readonly'
            }
        }
    }
}