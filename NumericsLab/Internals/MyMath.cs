// <copyright file="MyMath.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

namespace NumericsLab.Internals;

internal static class MyMath
{
    /// <summary>
    /// The difference between 1 and the smallest value greater than 1 that
    /// is representable for the data type double.
    /// Value taken from float.h.
    /// </summary>
    public const double Epsilon = 2.2204460492503131e-016;
}