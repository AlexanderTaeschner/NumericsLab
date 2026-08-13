// <copyright file="INumericVector.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Numerics;

namespace NumericsLab.Internals;

/// <summary>
/// Represents a one-based numeric vector with numeric elements.
/// </summary>
/// <typeparam name="T">The numeric type of the elements in the vector.</typeparam>
internal interface INumericVector<T>
    where T : INumber<T>
{
    /// <summary>
    /// Gets the number of elements in the vector.
    /// </summary>
    int Length { get; }

    /// <summary>
    /// Gets or sets the value of the element at the specified one-based index in the vector.
    /// </summary>
    /// <param name="index">The one-based index of the element to get or set.</param>
    /// <returns>The value of the element at the specified one-based index.</returns>
    T this[int index] { get; set; }

    /// <summary>
    /// Returns a read-only one-based list representation of the vector.
    /// </summary>
    /// <returns>A read-only one-based list containing the elements of the vector.</returns>
    IReadOnlyList<T> AsReadOnlyList();
}