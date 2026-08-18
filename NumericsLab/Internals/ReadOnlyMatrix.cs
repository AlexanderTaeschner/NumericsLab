// <copyright file="ReadOnlyMatrix.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Numerics;

namespace NumericsLab.Internals;

internal sealed class ReadOnlyMatrix<T>(T[,] matrix) : IReadOnlyMatrix<T>
    where T : INumber<T>
{
    public int NumberOfRows => matrix.GetLength(0);

    public int NumberOfColumns => matrix.GetLength(1);

    public T this[int rowIndex, int columnIndex] => matrix[rowIndex, columnIndex];
}