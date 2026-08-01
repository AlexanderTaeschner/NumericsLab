// <copyright file="Differentiate.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using NumericsLab.Internals;

namespace NumericsLab;

/// <summary>
/// Provides methods for numerical differentiation of one-dimensional functions.
/// </summary>
public sealed class Differentiate
{
    private static readonly double BaseStepSize = Math.Pow(MyMath.Epsilon, 1.0 / 3.0);

    /// <summary>
    /// Computes the first derivative of a one-dimensional function at a given point using central difference approximation.
    /// </summary>
    /// <param name="function">The function to differentiate.</param>
    /// <param name="x">The point at which to evaluate the derivative.</param>
    /// <returns>The first derivative of the function at the given point.</returns>
    public static double FirstDerivative(Func<double, double> function, double x)
    {
        double h = BaseStepSize * (Math.Abs(x) + 1.0);
        double fxph = function(x + h);
        double fxmh = function(x - h);
        return (fxph - fxmh) / (2.0 * h);
    }
}