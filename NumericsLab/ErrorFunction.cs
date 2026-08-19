//--------------------------------------------------------------------------
// <copyright file="ErrorFunction.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>
//--------------------------------------------------------------------------

/* Define this macro to suppress error propagation in exp(x^2)
   by using the expx2 function.  The tradeoff is that doing so
   generates two calls to the exponential function instead of one.  */
#define USE_EXPXSQ

namespace NumericsLab;

/// <summary>
/// Class providing special mathematical functions.
/// </summary>
public static partial class SpecialFunctions
{
    /* http://www.netlib.org/cephes/readme
    Some software in this archive may be from the book _Methods and
    Programs for Mathematical Functions_ (Prentice-Hall or Simon & Schuster
    International, 1989) or from the Cephes Mathematical Library, a
    commercial product. In either event, it is copyrighted by the author.
    What you see here may be used freely but it comes with no support or
    guarantee.

       The two known misprints in the book are repaired here in the
    source listings for the gamma function and the incomplete beta
    integral.


       Stephen L. Moshier
       moshier@na-net.ornl.gov
    */

    /*       ndtr.c (http://www.netlib.org/cephes/cprob.tgz)
     *
     * Normal distribution function
     *
     *
     *
     * SYNOPSIS:
     *
     * double x, y, ndtr();
     *
     * y = ndtr( x );
     *
     *
     *
     * DESCRIPTION:
     *
     * Returns the area under the Gaussian probability density
     * function, integrated from minus infinity to x:
     *
     *                            x
     *                             -
     *                   1        | |          2
     *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
     *               sqrt(2pi)  | |
     *                           -
     *                          -inf.
     *
     *             =  ( 1 + erf(z) ) / 2
     *             =  erfc(z) / 2
     *
     * where z = x/sqrt(2). Computation is via the functions
     * erf and erfc with care to avoid error amplification in computing exp(-x^2).
     *
     *
     * ACCURACY:
     *
     *                      Relative error:
     * arithmetic   domain     # trials      peak         rms
     *    IEEE     -13,0        30000       1.3e-15     2.2e-16
     *
     *
     * ERROR MESSAGES:
     *
     *   message         condition         value returned
     * erfc underflow    x > 37.519379347       0.0
     *
     */
    /*       erf.c
    *
    * Error function
    *
    *
    *
    * SYNOPSIS:
    *
    * double x, y, erf();
    *
    * y = erf( x );
    *
    *
    *
    * DESCRIPTION:
    *
    * The integral is
    *
    *                           x
    *                            -
    *                 2         | |          2
    *   erf(x)  =  --------     |    exp( - t  ) dt.
    *              sqrt(pi)   | |
    *                          -
    *                           0
    *
    * The magnitude of x is limited to 9.231948545 for DEC
    * arithmetic; 1 or -1 is returned outside this range.
    *
    * For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
    * erf(x) = 1 - erfc(x).
    *
    *
    *
    * ACCURACY:
    *
    *                      Relative error:
    * arithmetic   domain     # trials      peak         rms
    *    DEC       0,1         14000       4.7e-17     1.5e-17
    *    IEEE      0,1         30000       3.7e-16     1.0e-16
    *
    */
    /*       erfc.c
    *
    * Complementary error function
    *
    *
    *
    * SYNOPSIS:
    *
    * double x, y, erfc();
    *
    * y = erfc( x );
    *
    *
    *
    * DESCRIPTION:
    *
    *
    *  1 - erf(x) =
    *
    *                           inf.
    *                             -
    *                  2         | |          2
    *   erfc(x)  =  --------     |    exp( - t  ) dt
    *               sqrt(pi)   | |
    *                           -
    *                            x
    *
    *
    * For small x, erfc(x) = 1 - erf(x); otherwise rational
    * approximations are computed.
    *
    * A special function expx2.c is used to suppress error amplification
    * in computing exp(-x^2).
    *
    *
    * ACCURACY:
    *
    *                      Relative error:
    * arithmetic   domain     # trials      peak         rms
    *    IEEE      0,26.6417   30000       1.3e-15     2.2e-16
    *
    *
    * ERROR MESSAGES:
    *
    *   message         condition              value returned
    * erfc underflow    x > 9.231948545 (DEC)       0.0
    *
    *
    */

    /*
    Cephes Math Library Release 2.9:  November, 2000
    Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
    */

    /* Constants from const.c (http://www.netlib.org/cephes/cprob.tgz) */

    /// <summary>
    /// Logarithm of the maximum number.
    /// </summary>
    private const double MAXLOG = 7.09782712893383996732E2;     /* log(MAXNUM) */

    private const double MINLOG = -7.08396418532264106224E2;         /* log(2**-1022) */

    /* ndtr.c: */

    /// <summary>
    /// Constant used to calculate expx2.
    /// </summary>
    private const double M = 128.0;

    /// <summary>
    /// Constant used to calculate expx2.
    /// </summary>
    private const double MINV = 0.0078125;

    /// <summary>
    /// Constant used to calculate erfc.
    /// </summary>
    private static readonly double[] s_parameterComplementaryErrorFunction1 =
    [
        2.46196981473530512524E-10,
        5.64189564831068821977E-1,
        7.46321056442269912687E0,
        4.86371970985681366614E1,
        1.96520832956077098242E2,
        5.26445194995477358631E2,
        9.34528527171957607540E2,
        1.02755188689515710272E3,
        5.57535335369399327526E2,
    ];

    /// <summary>
    /// Constant used to calculate erfc.
    /// </summary>
    private static readonly double[] s_parameterComplementaryErrorFunction2 =
    [
        /*1.00000000000000000000E0,*/
        1.32281951154744992508E1,
        8.67072140885989742329E1,
        3.54937778887819891062E2,
        9.75708501743205489753E2,
        1.82390916687909736289E3,
        2.24633760818710981792E3,
        1.65666309194161350182E3,
        5.57535340817727675546E2,
    ];

    /// <summary>
    /// Constant used to calculate erfc.
    /// </summary>
    private static readonly double[] s_parameterComplementaryErrorFunction3 =
    [
        5.64189583547755073984E-1,
        1.27536670759978104416E0,
        5.01905042251180477414E0,
        6.16021097993053585195E0,
        7.40974269950448939160E0,
        2.97886665372100240670E0,
    ];

    /// <summary>
    /// Constant used to calculate erfc.
    /// </summary>
    private static readonly double[] s_parameterComplementaryErrorFunction4 =
    [
        /*1.00000000000000000000E0,*/
        2.26052863220117276590E0,
        9.39603524938001434673E0,
        1.20489539808096656605E1,
        1.70814450747565897222E1,
        9.60896809063285878198E0,
        3.36907645100081516050E0,
    ];

    /// <summary>
    /// Constant used to calculate erf.
    /// </summary>
    private static readonly double[] s_parameterErrorFunction1 =
    [
        9.60497373987051638749E0,
        9.00260197203842689217E1,
        2.23200534594684319226E3,
        7.00332514112805075473E3,
        5.55923013010394962768E4,
    ];

    /// <summary>
    /// Constant used to calculate erf.
    /// </summary>
    private static readonly double[] s_parameterErrorFunction2 =
    [
        /*1.00000000000000000000E0,*/
        3.35617141647503099647E1,
        5.21357949780152679795E2,
        4.59432382970980127987E3,
        2.26290000613890934246E4,
        4.92673942608635921086E4,
    ];

    /// <summary>
    /// Calculates the complementary error function.
    /// </summary>
    /// <param name="argument">Argument of the function.</param>
    /// <returns>The complementary error function.</returns>
    public static double ComplementaryErrorFunction(double argument)
    {
        double p, q, x, y, z;

        x = argument < 0.0 ? -argument : argument;

        if (x < 1.0)
        {
            return 1.0 - ErrorFunction(argument);
        }

        z = -argument * argument;

        if (z < -MAXLOG)
        {
            // under:
            // mtherr( "erfc", UNDERFLOW );
            return argument < 0 ? 2.0 : 0.0;
        }

#if USE_EXPXSQ
        /* Compute z = exp(z).  */
        z = ExpXSquared(argument, -1);
#else
            z = exp(z);
#endif
        if (x < 8.0)
        {
            p = PolynomialEvaluation(x, s_parameterComplementaryErrorFunction1, 8);
            q = PolynomialEvaluationOne(x, s_parameterComplementaryErrorFunction2, 8);
        }
        else
        {
            p = PolynomialEvaluation(x, s_parameterComplementaryErrorFunction3, 5);
            q = PolynomialEvaluationOne(x, s_parameterComplementaryErrorFunction4, 6);
        }

        y = z * p / q;

        if (argument < 0)
        {
            y = 2.0 - y;
        }

        if (y == 0.0)
        {
            // goto under;
            return argument < 0 ? 2.0 : 0.0;
        }

        return y;
    }

    /* Exponentially scaled erfc function
       exp(x^2) erfc(x)
       valid for x > 1.
       Use with ndtr and expx2.  * /
    static double erfce(x)
    double x;
    {
    double p,q;

    if( x < 8.0 )
        {
        p = polevl( x, P, 8 );
        q = p1evl( x, Q, 8 );
        }
    else
        {
        p = polevl( x, R, 5 );
        q = p1evl( x, S, 6 );
        }
    return (p/q);
    }
    */

    /*
    double erf(x)
    double x;
    */

    /// <summary>
    /// Calculates the error function.
    /// </summary>
    /// <param name="argument">Argument of the function.</param>
    /// <returns>The error function.</returns>
    public static double ErrorFunction(double argument)
    {
        double y, z;

        if (Math.Abs(argument) > 1.0)
        {
            return 1.0 - ComplementaryErrorFunction(argument);
        }

        z = argument * argument;
        double h = PolynomialEvaluationOne(z, s_parameterErrorFunction2, 5);
        y = argument * PolynomialEvaluation(z, s_parameterErrorFunction1, 4) / h;
        return y;
    }

    /*       expx2.c
     *
     * Exponential of squared argument
     *
     *
     *
     * SYNOPSIS:
     *
     * double x, y, expx2();
     * int sign;
     *
     * y = expx2( x, sign );
     *
     *
     *
     * DESCRIPTION:
     *
     * Computes y = exp(x*x) while suppressing error amplification
     * that would ordinarily arise from the inexactness of the
     * exponential argument x*x.
     *
     * If sign < 0, the result is inverted; i.e., y = exp(-x*x) .
     *
     *
     * ACCURACY:
     *
     *                      Relative error:
     * arithmetic    domain     # trials      peak         rms
     *   IEEE      -26.6, 26.6    10^7       3.9e-16     8.9e-17
     *
     */

    /*
    Cephes Math Library Release 2.9:  June, 2000
    Copyright 2000 by Stephen L. Moshier
    */

    /// <summary>
    /// Computes y = exp(x*x) while suppressing error amplification.
    /// </summary>
    /// <param name="x">The argument.</param>
    /// <param name="sign">If sign &lt; 0, the result is inverted; i.e., y = exp(-x*x).</param>
    /// <returns>The function result.</returns>
    private static double ExpXSquared(double x, int sign)
    {
        double u, u1, m, f;

        x = Math.Abs(x);
        if (sign < 0)
        {
            x = -x;
        }

        /* Represent x as an exact multiple of M plus a residual.
           M is a power of 2 chosen so that exp(m * m) does not overflow
           or underflow and so that |x - m| is small.  */
        m = MINV * Math.Floor((M * x) + 0.5);
        f = x - m;

        /* x^2 = m^2 + 2mf + f^2 */
        u = m * m;
        u1 = (2 * m * f) + (f * f);

        if (sign < 0)
        {
            u = -u;
            u1 = -u1;
        }

        if ((u + u1) > MAXLOG)
        {
            return double.PositiveInfinity;
        }

        /* u is exact, u1 is small.  */
        u = Math.Exp(u) * Math.Exp(u1);
        return u;
    }

    /*       polevl.c
     *       p1evl.c
     *
     * Evaluate polynomial
     *
     *
     *
     * SYNOPSIS:
     *
     * int N;
     * double x, y, coef[N+1], polevl[];
     *
     * y = polevl( x, coef, N );
     *
     *
     *
     * DESCRIPTION:
     *
     * Evaluates polynomial of degree N:
     *
     *                     2          N
     * y  =  C  + C x + C x  +...+ C x
     *        0    1     2          N
     *
     * Coefficients are stored in reverse order:
     *
     * coef[0] = C  , ..., coef[N] = C  .
     *            N                   0
     *
     *  The function p1evl() assumes that coef[N] = 1.0 and is
     * omitted from the array.  Its calling arguments are
     * otherwise the same as polevl().
     *
     *
     * SPEED:
     *
     * In the interest of speed, there are no checks for out
     * of bounds arithmetic.  This routine is used by most of
     * the functions in the library.  Depending on available
     * equipment features, the user may wish to rewrite the
     * program in microcode or assembly language.
     *
     */

    /*
    Cephes Math Library Release 2.1:  December, 1988
    Copyright 1984, 1987, 1988 by Stephen L. Moshier
    Direct inquiries to 30 Frost Street, Cambridge, MA 02140
    */

    /// <summary>
    /// Evaluates polynomial of degree N.
    /// </summary>
    /// <param name="x">The argument.</param>
    /// <param name="coef">The coefficients of the polynomial in reversed order.</param>
    /// <param name="n">The degree of the polynomial.</param>
    /// <returns>The function result.</returns>
    private static double PolynomialEvaluation(double x, double[] coef, int n)
    {
        int p = 0;
        double ans = coef[p];
        p++;

        do
        {
            ans = (ans * x) + coef[p];
            p++;
        }
        while (p <= n);

        return ans;
    }

    /*       p1evl() */
    /*                                          N
     * Evaluate polynomial when coefficient of x  is 1.0.
     * Otherwise same as polevl.
     */

    /// <summary>
    /// Evaluates polynomial when coefficient of x^N  is 1.0.
    /// </summary>
    /// <param name="x">The argument.</param>
    /// <param name="coef">The coefficients of the polynomial in reversed order.</param>
    /// <param name="n">The degree of the polynomial.</param>
    /// <returns>The function result.</returns>
    private static double PolynomialEvaluationOne(double x, double[] coef, int n)
    {
        double ans;
        int p;
        int i;

        p = 0;
        ans = x + coef[p];
        p++;
        i = n - 1;

        do
        {
            ans = (ans * x) + coef[p];
            p++;
            i--;
        }
        while (i != 0);

        return ans;
    }
}
