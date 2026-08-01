// <copyright file="CubicSpline.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using NumericsLab.Internals;

namespace NumericsLab.Interpolation;

/// <summary>
/// Represents a cubic spline interpolation of a set of data points.
/// </summary>
public sealed class CubicSpline
{
    private readonly double[] _x;
    private readonly double[] _c0;
    private readonly double[] _c1;
    private readonly double[] _c2;
    private readonly double[] _c3;

    private CubicSpline(double[] x, double[] c0, double[] c1, double[] c2, double[] c3)
    {
        _x = x;
        _c0 = c0;
        _c1 = c1;
        _c2 = c2;
        _c3 = c3;
    }

    private CubicSpline() => throw new NotSupportedException("Default constructor is not supported.");

    /// <summary>
    /// Creates an Akima spline interpolation for the given data points.
    /// </summary>
    /// <param name="x">The x-coordinates of the data points.</param>
    /// <param name="y">The y-coordinates of the data points.</param>
    /// <returns>A cubic spline representing the Akima interpolation of the data points.</returns>
    /// <exception cref="ArgumentException">Thrown when the input arrays have different lengths or contain fewer than 5 points.</exception>
    public static CubicSpline CreateAkimaSpline(double[] x, double[] y)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Input arrays must have the same length.");
        }

        int numPoints = x.Length;

        if (numPoints < 5)
        {
            throw new ArgumentException("At least 5 data points are required for Akima interpolation.");
        }

        int numSegments = numPoints - 1;

        // Calculate the slopes of the segments - see https://en.wikipedia.org/wiki/Akima_spline for details
        double[] m = new double[numSegments];
        for (int i = 0; i < numSegments; i++)
        {
            m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]); // m_i is the slope of the segment between (x_i, y_i) and (x_{i+1}, y_{i+1})
        }

        double mm1 = (2 * m[0]) - m[1]; // m[-1], see 2.3 of Akima's paper
        double mm2 = (2 * mm1) - m[0]; // m[-2]
        double mN = (2 * m[^1]) - m[^2]; // m[N]
        double mNp1 = (2 * mN) - m[^1]; // m[N+1]

        double[] w = new double[m.Length - 1];
        for (int i = 0; i < w.Length; i++)
        {
            w[i] = Math.Abs(m[i + 1] - m[i]);
        }

        double Slope(int i)
        {
            if (i == -2)
            {
                return mm2;
            }

            if (i == -1)
            {
                return mm1;
            }

            if (i == numSegments)
            {
                return mN;
            }

            if (i == numSegments + 1)
            {
                return mNp1;
            }

            return m[i];
        }

        double Weight(int i)
        {
            if (i < 0 || i >= w.Length)
            {
                return Math.Abs(Slope(i + 1) - Slope(i));
            }

            return w[i];
        }

        double[] t = new double[numPoints];
        for (int i = 0; i < numPoints; i++)
        {
            double sum = Weight(i) + Weight(i - 2);
            if (sum < MyMath.Epsilon)
            {
                t[i] = (Slope(i - 1) + Slope(i)) / 2;
            }
            else
            {
                t[i] = ((Weight(i) * Slope(i - 1)) + (Weight(i - 2) * Slope(i))) / sum;
            }
        }

        // Implementation of Akima interpolation
        return CreateSplineViaHermiteInterpolation(x, y, t, m);
    }

    /// <summary>
    /// Interpolates the value at the given x-coordinate using the cubic spline.
    /// </summary>
    /// <param name="x">The x-coordinate at which to interpolate the value.</param>
    /// <returns>The interpolated value at the given x-coordinate.</returns>
    public double Interpolate(double x)
    {
        int i = IndexOfKnot(x);
        double dx = x - _x[i];
        return _c0[i] + (dx * (_c1[i] + (dx * (_c2[i] + (dx * _c3[i])))));
    }

    private static CubicSpline CreateSplineViaHermiteInterpolation(double[] x, double[] y, double[] firstDerivatives, double[] segmentSlopes)
    {
        double[] c0 = new double[x.Length - 1];
        double[] c1 = new double[x.Length - 1];
        double[] c2 = new double[x.Length - 1];
        double[] c3 = new double[x.Length - 1];

        for (int i = 0; i < c0.Length; i++)
        {
            double w = x[i + 1] - x[i];
            double m = segmentSlopes[i];
            double w2 = w * w;
            c0[i] = y[i];
            c1[i] = firstDerivatives[i];
            c2[i] = ((3 * m) - (2 * firstDerivatives[i]) - firstDerivatives[i + 1]) / w;
            c3[i] = (firstDerivatives[i] + firstDerivatives[i + 1] - (2 * m)) / w2;
        }

        return new CubicSpline(x, c0, c1, c2, c3);
    }

    private int IndexOfKnot(double x)
    {
        int i = Array.BinarySearch(_x, x);
        if (i < 0)
        {
            i = ~i - 1;
        }

        return Math.Clamp(i, 0, _x.Length - 2);
    }
}