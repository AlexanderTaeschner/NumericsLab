//--------------------------------------------------------------------------
// <copyright file="Brent.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>
//--------------------------------------------------------------------------

using NumericsLab.Internals;

namespace NumericsLab.RootFinding;

/// <summary>
/// Class used to find the root of a function using Brent's method.
/// </summary>
public static class Brent
{
    /*
     * Source of the C code: https://netlib.org/c/brent.shar
    file brent.shar
    for Brent's univariate minimizer and zero finder.
    by Oleg Keselyov <oleg@ponder.csci.unt.edu, oleg@unt.edu> May 23, 1991
    ref G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical computations.
    # Contains the source code for the program fminbr.c and
    # zeroin.c, test drivers for both, and verifivation protocols.
    */

    /*
     ************************************************************************
     *           C math library
     * function ZEROIN - obtain a function zero within the given range
     *
     * Input
     * double zeroin(ax,bx,f,tol)
     * double ax;    Root will be seeked for within
     * double bx;     a range [ax,bx]
     * double (*f)(double x);  Name of the function whose zero
     *     will be seeked for
     * double tol;   Acceptable tolerance for the root
     *     value.
     *     May be specified as 0.0 to cause
     *     the program to find the root as
     *     accurate as possible
     *
     * Output
     * Zeroin returns an estimate for the root with accuracy
     * 4*EPSILON*abs(x) + tol
     *
     * Algorithm
     * G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
     * computations. M., Mir, 1980, p.180 of the Russian edition
     *
     * The function makes use of the bissection procedure combined with
     * the linear or quadric inverse interpolation.
     * At every step program operates on three abscissae - a, b, and c.
     * b - the last and the best approximation to the root
     * a - the last but one approximation
     * c - the last but one or even earlier approximation than a that
     *  1) |f(b)| <= |f(c)|
     *  2) f(b) and f(c) have opposite signs, i.e. b and c confine
     *     the root
     * At every step Zeroin selects one of the two new approximations, the
     * former being obtained by the bissection procedure and the latter
     * resulting in the interpolation (if a,b, and c are all different
     * the quadric interpolation is utilized, otherwise the linear one).
     * If the latter (i.e. obtained by the interpolation) point is
     * reasonable (i.e. lies within the current interval [b,c] not being
     * too close to the boundaries) it is accepted. The bissection result
     * is used in the other case. Therefore, the range of uncertainty is
     * ensured to be reduced at least by the factor 1.6
     *
     ************************************************************************
     */
    /*
    double zeroin(ax,bx,f,tol)  /* An estimate to the root * /
    double ax;    /* Left border | of the range * /
    double bx;      /* Right border| the root is seeked* /
    double (*f)(double x);   /* Function under investigation * /
    double tol;    /* Acceptable tolerance  * /
    */

    /// <summary>
    /// Find the root of the specified function.
    /// </summary>
    /// <param name="function">The function.</param>
    /// <param name="lowerLimit">The lower limit.</param>
    /// <param name="upperLimit">The upper limit.</param>
    /// <param name="tolerance">The tolerance.</param>
    /// <returns>The root of the function.</returns>
    public static double FindRoot(
        Func<double, double> function,
        double lowerLimit,
        double upperLimit,
        double tolerance)
    {
        if (!TryFindRoot(function, lowerLimit, upperLimit, tolerance, int.MaxValue, out double root))
        {
            throw new InvalidOperationException("Root finding failed.");
        }

        return root;
    }

    /// <summary>
    /// Find the root of the specified function.
    /// </summary>
    /// <param name="function">The function.</param>
    /// <param name="lowerLimit">The lower limit.</param>
    /// <param name="upperLimit">The upper limit.</param>
    /// <param name="tolerance">The tolerance.</param>
    /// <param name="maxIterations">The maximum number of iterations.</param>
    /// <param name="root">When this method returns, contains the root of the function if found; otherwise, 0.</param>
    /// <returns>True if the root was found; otherwise, false.</returns>
    public static bool TryFindRoot(
        Func<double, double> function,
        double lowerLimit,
        double upperLimit,
        double tolerance,
        int maxIterations,
        out double root)
    {
        double a, b, c;    /* Abscissae, descr. see above */
        double fa;    /* f(a)    */
        double fb;    /* f(b)    */
        double fc;    /* f(c)    */

        a = lowerLimit;
        b = upperLimit;
        fa = function(a);
        fb = function(b);
        c = a;
        fc = fa;

        int iterations = 0;

        while (true)
        {
            /* Main iteration loop */
            double prev_step = b - a;  /* Distance from the last but one*/
            /* to the last approximation */
            double tol_act;   /* Actual tolerance  */
            double p;         /* Interpolation step is calcu- */
            double q;         /* lated in the form p/q; divi- */
            /* sion operations is delayed   */
            /* until the last moment */
            double new_step;        /* Step at this iteration       */

            if (Math.Abs(fc) < Math.Abs(fb))
            { /* Swap data for b to be the  */
                a = b; /* best approximation  */
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            tol_act = (2 * MyMath.Epsilon * Math.Abs(b)) + (tolerance / 2);
            new_step = (c - b) / 2;

            if (Math.Abs(new_step) <= tol_act || fb == 0)
            {
                root = b;
                return true;    /* Acceptable approx. is found */
            }

            iterations++;
            if (iterations >= maxIterations)
            {
                root = b;
                return false;   /* Maximum iterations reached */
            }

            /* Decide if the interpolation can be tried */
            /* If prev_step was large enough and was in true direction, */
            /* Interpolatiom may be tried */
            if (Math.Abs(prev_step) >= tol_act
            && Math.Abs(fa) > Math.Abs(fb))
            {
                double t1, cb, t2;
                cb = c - b;
                if (a == c)
                {
                    /* If we have only two distinct */
                    /* points linear interpolation  */
                    t1 = fb / fa;   /* can only be applied  */
                    p = cb * t1;
                    q = 1.0 - t1;
                }
                else
                {
                    /* Quadric inverse interpolation*/
                    q = fa / fc;
                    t1 = fb / fc;
                    t2 = fb / fa;
                    p = t2 * ((cb * q * (q - t1)) - ((b - a) * (t1 - 1.0)));
                    q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
                }

                if (p > 0)
                {
                    /* p was calculated with the op-*/
                    q = -q;   /* posite sign; make p positive */
                }
                else
                {
                    /* and assign possible minus to */
                    p = -p;   /* q    */
                }

                /* If b+p/q falls in [b,c] and isn't too large it is accepted */
                if (p < ((0.75 * cb * q) - (Math.Abs(tol_act * q) / 2))
                    && p < Math.Abs(prev_step * q / 2))
                {
                    new_step = p / q;
                }

                /* If p/q is too large then the */
                /* bissection procedure can  */
                /* reduce [b,c] range to more */
                /* extent   */
            }

            if (Math.Abs(new_step) < tol_act)
            {
                /* Adjust the step to be not less than tolerance*/
                new_step = new_step > 0 ? tol_act : -tol_act;
            }

            a = b;  /* Save the previous approx. */
            fa = fb;
            b += new_step;
            fb = function(b); /* Do step to a new approxim. */
            if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
            { /* Adjust c for it to have a sign*/
                c = a; /* opposite to that of b */
                fc = fa;
            }
        }
    }
}