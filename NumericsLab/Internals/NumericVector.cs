// <copyright file="NumericVector.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using System.Numerics;

namespace NumericsLab.Internals;

/// <summary>
/// Represents a one-based numeric vector with numeric elements.
/// </summary>
internal sealed class NumericVector<T> : INumericVector<T>
    where T : INumber<T>
{
    private readonly T[] _array;

    public NumericVector(T[] array) => _array = array;

    public NumericVector(int length) => _array = new T[length];

    /// <inheritdoc/>
    public int Length => _array.Length;

    /// <inheritdoc/>
    public T this[int index]
    {
        get => _array[index - 1];
        set => _array[index - 1] = value;
    }

    /// <inheritdoc/>
    public IReadOnlyList<T> AsReadOnlyList() => _array;

    /// <summary>
    /// Gets the underlying array.
    /// </summary>
    /// <returns>The internal array.</returns>
    internal T[] AsArray() => _array;

    /// <summary>
    /// Sets the value of the element at index i in the underlying array to v.
    /// </summary>
    /// <param name="i">The zero-based index of the element to set.</param>
    /// <param name="value">The value to set the element to.</param>
    internal void SetArrayElement(int i, T value) => _array[i] = value;
}