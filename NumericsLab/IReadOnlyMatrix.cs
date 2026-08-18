// <copyright file="IReadOnlyMatrix.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Numerics;

namespace NumericsLab;

/// <summary>
/// Represents a read-only matrix with numeric elements.
/// </summary>
/// <typeparam name="T">The type of the numeric elements.</typeparam>
public interface IReadOnlyMatrix<T>
    where T : INumber<T>
{
    /// <summary>
    /// Gets the number of rows in the matrix.
    /// </summary>
    int NumberOfRows { get; }

    /// <summary>
    /// Gets the number of columns in the matrix.
    /// </summary>
    int NumberOfColumns { get; }

    /// <summary>
    /// Gets the value of the element at the specified row and column indices in the matrix.
    /// </summary>
    /// <param name="rowIndex">The zero-based index of the row.</param>
    /// <param name="columnIndex">The zero-based index of the column.</param>
    /// <returns>The value of the element at the specified row and column indices.</returns>
    T this[int rowIndex, int columnIndex] { get; }
}