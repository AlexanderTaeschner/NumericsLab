// <copyright file="NonlinearLeastSquares.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using NumericsLab.Internals;
using System.Diagnostics.CodeAnalysis;

namespace NumericsLab.Fitting;

/// <summary>
/// Represents a model function for nonlinear least squares fitting.
/// </summary>
/// <param name="x">The independent variable.</param>
/// <param name="parameters">The parameters of the model function.</param>
/// <returns>The value of the model function at the given x with the specified parameters.</returns>
public delegate double ModelFunction(double x, IReadOnlyList<double> parameters);

/// <summary>
/// Provides methods for performing nonlinear least squares fitting using the Levenberg-Marquardt algorithm.
/// </summary>
public static class NonlinearLeastSquares
{
    private static readonly double EpsForwardDifferenceJacobian = Math.Sqrt(MyMath.Epsilon); // lmdif1 uses epsfcn=0, so we simpilify the max(epsfcn, epsmch) to just epsmch here.

    /// <summary>
    /// Represents a function that computes the residuals for nonlinear least squares fitting.
    /// F_i(parameters) = y_i - g(x_i, parameters) where g is the model function.
    /// </summary>
    /// <param name="parameters">The parameters at which to evaluate the residuals.</param>
    /// <param name="functionValues">The computed function values.
    /// In case the function values vector is null at the time of the call, a new NumericVector will be created inside this function.</param>
    private delegate void NLSFunction(NumericVector<double> parameters, [NotNull] ref NumericVector<double>? functionValues);

    /// <summary>
    /// Fits the specified model function to the given data using the Levenberg-Marquardt algorithm.
    /// </summary>
    /// <param name="model">The model function to fit.</param>
    /// <param name="xdata">The independent variable data.</param>
    /// <param name="ydata">The dependent variable data.</param>
    /// <param name="initialParameters">The initial guess for the parameters.</param>
    /// <param name="tolerance">The tolerance for the fitting algorithm.</param>
    /// <returns>The result of the fit.</returns>
    public static FitResult Fit(ModelFunction model, IReadOnlyList<double> xdata, IReadOnlyList<double> ydata, IReadOnlyList<double> initialParameters, double tolerance = MyMath.EpsilonSqrt)
    {
        void ResidualFunction(NumericVector<double> parameters, [NotNull] ref NumericVector<double>? functionValues)
        {
            functionValues ??= new NumericVector<double>(xdata.Count);

            for (int i = 0; i < xdata.Count; i++)
            {
                functionValues.SetArrayElement(i, ydata[i] - model(xdata[i], parameters.AsReadOnlyList()));
            }
        }

        NumericVector<double> x = new([.. initialParameters]); // initialize x with a copy of the initial parameters
        int maxfev = 200 * (x.Length + 1);
        return LevenbergMarquardt(ResidualFunction, x, ftol: tolerance, xtol: tolerance, gtol: 0.0, maxfev);
    }

    internal static double EuclidianNorm(int n, INumericVector<double> x)
    {
        // Based on the subroutine enorm from MINPACK (see https://netlib.org/minpack/), ported to C#

        // given an n-vector x, this function calculates the
        // euclidean norm of x.
        //
        // the euclidean norm is computed by accumulating the sum of
        // squares in three different sums.the sums of squares for the
        // small and large components are scaled so that no overflows
        // occur.non - destructive underflows are permitted.underflows
        // and overflows do not occur in the computation of the unscaled
        // sum of squares for the intermediate components.
        // the definitions of small, intermediate and large components
        // depend on two constants, rdwarf and rgiant.the main
        // restrictions on these constants are that rdwarf * *2 not
        // underflow and rgiant * *2 not overflow.the constants
        // given here are suitable for every known computer.
        //
        // argonne national laboratory.minpack project.march 1980.
        // burton s.garbow, kenneth e.hillstrom, jorge j.more
        const double rdwarf = 3.834e-20;
        const double rgiant = 1.304e19;

        double s1 = 0.0;
        double s2 = 0.0;
        double s3 = 0.0;
        double x1max = 0.0;
        double x3max = 0.0;
        double floatn = n;
        double agiant = rgiant / floatn;

        for (int i = 1; i <= n; i++)
        {
            double xabs = Math.Abs(x[i]);
            if (xabs > rdwarf && xabs < agiant)
            {
                // sum for intermediate components
                s2 += xabs * xabs;
            }
            else
            {
                if (xabs <= rdwarf)
                {
                    // sum for small components.
                    if (xabs <= x3max)
                    {
                        if (xabs != 0.0)
                        {
                            s3 += Math.Pow(xabs / x3max, 2);
                        }
                    }
                    else
                    {
                        s3 = 1.0 + (s3 * Math.Pow(x3max / xabs, 2));
                        x3max = xabs;
                    }
                }
                else
                {
                    // sum for large components.
                    if (xabs <= x1max)
                    {
                        s1 += Math.Pow(xabs / x1max, 2);
                    }
                    else
                    {
                        s1 = 1.0 + (s1 * Math.Pow(x1max / xabs, 2));
                        x1max = xabs;
                    }
                }
            }
        }

        // calculation of norm.
        double enorm;
        if (s1 != 0.0)
        {
            enorm = x1max * Math.Sqrt(s1 + (s2 / x1max / x1max));
        }
        else
        {
            if (s2 == 0.0)
            {
                enorm = x3max * Math.Sqrt(s3);
                return enorm;
            }
            else
            {
                if (s2 >= x3max)
                {
                    enorm = Math.Sqrt(s2 * (1.0 + (x3max / s2 * (x3max * s3))));
                }
                else
                {
                    enorm = Math.Sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
                }
            }
        }

        return enorm;
    }

    private static FitResult LevenbergMarquardt(
        NLSFunction function,
        NumericVector<double> x,
        double ftol,
        double xtol,
        double gtol,
        int maxfev,
        int mode = 1,
        double factor = 100.0,
        NumericVector<double> diag = null!)
    {
        // Based on the LMDIF subroutine from MINPACK (see https://netlib.org/minpack/), ported to C#
        /*
      subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *                 diag,mode,factor,nprint,info,nfev,fjac,ldfjac,
     *                 ipvt,qtf,wa1,wa2,wa3,wa4)
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      external fcn
c     **********
c
c     subroutine lmdif
c
c     the purpose of lmdif is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
c                        diag,mode,factor,nprint,info,nfev,fjac,
c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,
     *                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
     *                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
         */

        // epsmch is the machine precision.
        double epsmch = MyMath.Epsilon;

        int info = 0; // info is the return status of the algorithm
        int nfev; // number of function evaluations

        // evaluate the function at the starting point
        // and calculate its norm.
        int n = x.Length;
        NumericVector<double> fvec = null!;
        function(x, ref fvec);
        int m = fvec.Length;
        nfev = 1;
        NumericMatrix<double> fjac = new(m, n);
        double fnorm = EuclidianNorm(m, fvec);

        NumericVector<int> ipvt = new(n);
        NumericVector<double> qtf = new(n);
        NumericVector<double> wa1 = new(n);
        NumericVector<double> wa2 = new(n);
        NumericVector<double> wa3 = new(n);
        NumericVector<double> wa4 = new(m);
        diag ??= new(n);

        double xnorm = double.NaN;
        double delta = double.NaN;

        // initialize levenberg-marquardt parameter and iteration counter.
        double par = 0.0;
        int iter = 1;

        // beginning of the outer loop.
        do
        {
            // calculate the jacobian matrix.
            ForwardDifferenceJacobian(function, x, fvec, ref fjac);
            nfev += n;

            // compute the qr factorization of the jacobian.
            QRFactorization(m, n, ref fjac, true, ref ipvt, ref wa1, ref wa2, ref wa3);

            // on the first iteration and if mode is 1, scale according
            // to the norms of the columns of the initial jacobian.
            if (iter == 1)
            {
                if (mode != 2)
                {
                    for (int j = 1; j <= n; j++)
                    {
                        diag[j] = wa2[j];
                        if (wa2[j] == 0)
                        {
                            diag[j] = 1.0;
                        }
                    }
                }

                // on the first iteration, calculate the norm of the scaled x
                // and initialize the step bound delta.
                for (int j = 1; j <= n; j++)
                {
                    wa3[j] = diag[j] * x[j];
                }

                xnorm = EuclidianNorm(n, wa3);
                delta = factor * xnorm;
                if (delta == 0.0)
                {
                    delta = factor;
                }
            }

            // form(q transpose) * fvec and store the first n components in qtf.
            for (int i = 1; i <= m; i++)
            {
                wa4[i] = fvec[i];
            }

            for (int j = 1; j <= n; j++)
            {
                if (fjac[j, j] != 0.0)
                {
                    double sum = 0.0;
                    for (int i = j; i <= m; i++)
                    {
                        sum += fjac[i, j] * wa4[i];
                    }

                    double temp = -sum / fjac[j, j];
                    for (int i = j; i <= m; i++)
                    {
                        wa4[i] = wa4[i] + (fjac[i, j] * temp);
                    }
                }

                fjac[j, j] = wa1[j];
                qtf[j] = wa4[j];
            }

            // compute the norm of the scaled gradient.
            double gnorm = 0.0;
            if (fnorm != 0.0)
            {
                for (int j = 1; j <= n; j++)
                {
                    int l = ipvt[j];
                    if (wa2[l] != 0.0)
                    {
                        double sum = 0.0;
                        for (int i = 1; i <= j; i++)
                        {
                            sum += fjac[i, j] * (qtf[i] / fnorm);
                        }

                        gnorm = Math.Max(gnorm, Math.Abs(sum / wa2[l]));
                    }
                }
            }

            // test for convergence of the gradient norm.
            if (gnorm <= gtol)
            {
                info = 4;
            }

            if (info != 0)
            {
                return new FitResult(info, x, fvec, fjac, ipvt);
            }

            // rescale if necessary.
            if (mode != 2)
            {
                for (int j = 1; j <= n; j++)
                {
                    diag[j] = Math.Max(diag[j], wa2[j]);
                }
            }

            double ratio;

            // beginning of the inner loop.
            do
            {
                // determine the levenberg - marquardt parameter.
                LevenbergMarquardtParameter(n, fjac, ipvt, diag, qtf, delta, ref par, wa1, wa2, wa3, wa4);

                // store the direction p and x +p.calculate the norm of p.
                for (int j = 1; j <= n; j++)
                {
                    wa1[j] = -wa1[j];
                    wa2[j] = x[j] + wa1[j];
                    wa3[j] = diag[j] * wa1[j];
                }

                double pnorm = EuclidianNorm(n, wa3);

                // on the first iteration, adjust the initial step bound.
                if (iter == 1)
                {
                    delta = Math.Min(delta, pnorm);
                }

                // evaluate the function at x + p and calculate its norm.
                function(wa2, ref wa4);
                nfev++;

                // if (iflag.lt. 0) go to 300
                double fnorm1 = EuclidianNorm(m, wa4);

                // compute the scaled actual reduction.
                double actred = -1.0;
                if (0.1 * fnorm1 < fnorm)
                {
                    actred = 1.0 - ((fnorm1 / fnorm) * (fnorm1 / fnorm));
                }

                // compute the scaled predicted reduction and
                // the scaled directional derivative.
                for (int j = 1; j <= n; j++)
                {
                    wa3[j] = 0.0;
                    int l = ipvt[j];
                    double temp = wa1[l];
                    for (int i = 1; i <= j; i++)
                    {
                        wa3[i] = wa3[i] + (fjac[i, j] * temp);
                    }
                }

                double temp1 = EuclidianNorm(n, wa3) / fnorm;
                double temp2 = (Math.Sqrt(par) * pnorm) / fnorm;
                double prered = (temp1 * temp1) + ((temp2 * temp2) / 0.5);
                double dirder = -((temp1 * temp1) + (temp2 * temp2));

                // compute the ratio of the actual to the predicted
                // reduction.
                ratio = 0.0;
                if (prered != 0.0)
                {
                    ratio = actred / prered;
                }

                // update the step bound.
                if (ratio <= 0.25)
                {
                    double temp;
                    if (actred >= 0.0)
                    {
                        temp = 0.5;
                    }
                    else
                    {
                        temp = 0.5 * dirder / (dirder + (0.5 * actred));
                    }

                    if (0.1 * fnorm1 >= fnorm || temp < 0.1)
                    {
                        temp = 0.1;
                    }

                    delta = temp * Math.Min(delta, pnorm / 0.1);
                    par /= temp;
                }
                else
                {
                    if (par == 0.0 || ratio >= 0.75)
                    {
                        delta = pnorm / 0.5;
                        par = 0.5 * par;
                    }
                }

                // test for successful iteration.
                if (ratio >= 0.0001)
                {
                    // successful iteration. update x, fvec, and their norms.
                    for (int j = 1; j <= n; j++)
                    {
                        x[j] = wa2[j];
                        wa2[j] = diag[j] * x[j];
                    }

                    for (int i = 1; i <= m; i++)
                    {
                        fvec[i] = wa4[i];
                    }

                    xnorm = EuclidianNorm(n, wa2);
                    fnorm = fnorm1;
                    iter++;
                }

                // tests for convergence.
                if (Math.Abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1.0)
                {
                    info = 1;
                }

                if (delta <= xtol * xnorm)
                {
                    info = 2;
                }

                if (Math.Abs(actred) <= ftol && prered <= ftol
                  && 0.5 * ratio <= 1.0 && info == 2)
                {
                    info = 3;
                }

                if (info != 0)
                {
                    return new(info, x, fvec, fjac, ipvt);
                }

                // tests for termination and stringent tolerances.
                if (nfev >= maxfev)
                {
                    info = 5;
                }

                if (Math.Abs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1.0)
                {
                    info = 6;
                }

                if (delta <= epsmch * xnorm)
                {
                    info = 7;
                }

                if (gnorm <= epsmch)
                {
                    info = 8;
                }

                if (info != 0)
                {
                    return new(info, x, fvec, fjac, ipvt);
                }

                // end of the inner loop.repeat if iteration unsuccessful.
            }
            while (ratio < 0.0001);

            // end of the outer loop.
        }
        while (true);
    }

    private static void LevenbergMarquardtParameter(
        int n,
        NumericMatrix<double> r,
        NumericVector<int> ipvt,
        NumericVector<double> diag,
        NumericVector<double> qtb,
        double delta,
        ref double par,
        NumericVector<double> x,
        NumericVector<double> sdiag,
        NumericVector<double> wa1,
        NumericVector<double> wa2)
    {
        // subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,
        // *                 wa2)
        // integer n,ldr
        // integer ipvt(n)
        // double precision delta,par
        // double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),
        // *                 wa2(n)
        // c     **********
        // c
        // c     subroutine lmpar
        // c
        // c     given an m by n matrix a, an n by n nonsingular diagonal
        // c     matrix d, an m-vector b, and a positive number delta,
        // c     the problem is to determine a value for the parameter
        // c     par such that if x solves the system
        // c
        // c           a*x = b ,     sqrt(par)*d*x = 0 ,
        // c
        // c     in the least squares sense, and dxnorm is the euclidean
        // c     norm of d*x, then either par is zero and
        // c
        // c           (dxnorm-delta) .le. 0.1*delta ,
        // c
        // c     or par is positive and
        // c
        // c           abs(dxnorm-delta) .le. 0.1*delta .
        // c
        // c     this subroutine completes the solution of the problem
        // c     if it is provided with the necessary information from the
        // c     qr factorization, with column pivoting, of a. that is, if
        // c     a*p = q*r, where p is a permutation matrix, q has orthogonal
        // c     columns, and r is an upper triangular matrix with diagonal
        // c     elements of nonincreasing magnitude, then lmpar expects
        // c     the full upper triangle of r, the permutation matrix p,
        // c     and the first n components of (q transpose)*b. on output
        // c     lmpar also provides an upper triangular matrix s such that
        // c
        // c            t   t                   t
        // c           p *(a *a + par*d*d)*p = s *s .
        // c
        // c     s is employed within lmpar and may be of separate interest.
        // c
        // c     only a few iterations are generally needed for convergence
        // c     of the algorithm. if, however, the limit of 10 iterations
        // c     is reached, then the output par will contain the best
        // c     value obtained so far.
        // c
        // c     the subroutine statement is
        // c
        // c       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
        // c                        wa1,wa2)
        // c
        // c     where
        // c
        // c       n is a positive integer input variable set to the order of r.
        // c
        // c       r is an n by n array. on input the full upper triangle
        // c         must contain the full upper triangle of the matrix r.
        // c         on output the full upper triangle is unaltered, and the
        // c         strict lower triangle contains the strict upper triangle
        // c         (transposed) of the upper triangular matrix s.
        // c
        // c       ldr is a positive integer input variable not less than n
        // c         which specifies the leading dimension of the array r.
        // c
        // c       ipvt is an integer input array of length n which defines the
        // c         permutation matrix p such that a*p = q*r. column j of p
        // c         is column ipvt(j) of the identity matrix.
        // c
        // c       diag is an input array of length n which must contain the
        // c         diagonal elements of the matrix d.
        // c
        // c       qtb is an input array of length n which must contain the first
        // c         n elements of the vector (q transpose)*b.
        // c
        // c       delta is a positive input variable which specifies an upper
        // c         bound on the euclidean norm of d*x.
        // c
        // c       par is a nonnegative variable. on input par contains an
        // c         initial estimate of the levenberg-marquardt parameter.
        // c         on output par contains the final estimate.
        // c
        // c       x is an output array of length n which contains the least
        // c         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
        // c         for the output par.
        // c
        // c       sdiag is an output array of length n which contains the
        // c         diagonal elements of the upper triangular matrix s.
        // c
        // c       wa1 and wa2 are work arrays of length n.
        // c
        // c     subprograms called
        // c
        // c       minpack-supplied ... dpmpar,enorm,qrsolv
        // c
        // c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
        // c
        // c     argonne national laboratory. minpack project. march 1980.
        // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
        // c
        // c     **********
        // integer i,iter,j,jm1,jp1,k,l,nsing
        // double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,
        // *                 sum,temp,zero
        // double precision dpmpar,enorm
        // data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/

        // dwarf is the smallest positive magnitude.
        double dwarf = MyMath.SmallestMagnitude;

        // compute and store in x the gauss-newton direction. if the
        // jacobian is rank-deficient, obtain a least squares solution.
        int nsing = n;
        for (int j = 1; j <= n; j++)
        {
            wa1[j] = qtb[j];
            if (r[j, j] == 0.0 && nsing == n)
            {
                nsing = j - 1;
            }

            if (nsing < n)
            {
                wa1[j] = 0.0;
            }
        }

        if (nsing >= 1)
        {
            for (int k = 1; k <= nsing; k++)
            {
                int j = nsing - k + 1;
                wa1[j] = wa1[j] / r[j, j];
                double temp = wa1[j];
                int jm1 = j - 1;
                if (jm1 >= 1)
                {
                    for (int i = 1; i <= jm1; i++)
                    {
                        wa1[i] = wa1[i] - (r[i, j] * temp);
                    }
                }
            }
        }

        for (int j = 1; j <= n; j++)
        {
            int l = ipvt[j];
            x[l] = wa1[j];
        }

        // initialize the iteration counter.
        // evaluate the function at the origin, and test
        // for acceptance of the gauss-newton direction.
        int iter = 0;
        for (int j = 1; j <= n; j++)
        {
            wa2[j] = diag[j] * x[j];
        }

        double dxnorm = EuclidianNorm(n, wa2);
        double fp = dxnorm - delta;
        if (fp > 0.1 * delta)
        {
            // if the jacobian is not rank deficient, the newton
            // step provides a lower bound, parl, for the zero of
            // the function. otherwise set this bound to zero.
            double parl = 0.0;
            if (nsing >= n)
            {
                for (int j = 1; j <= n; j++)
                {
                    int l = ipvt[j];
                    wa1[j] = diag[l] * (wa2[l] / dxnorm);
                }

                for (int j = 1; j <= n; j++)
                {
                    double sum = 0.0;
                    int jm1 = j - 1;
                    if (jm1 >= 1)
                    {
                        for (int i = 1; i <= jm1; i++)
                        {
                            sum += r[i, j] * wa1[i];
                        }
                    }

                    wa1[j] = (wa1[j] - sum) / r[j, j];
                }

                double temp = EuclidianNorm(n, wa1);
                parl = ((fp / delta) / temp) / temp;
            }

            // calculate an upper bound, paru, for the zero of the function.
            for (int j = 1; j <= n; j++)
            {
                double sum = 0.0;
                for (int i = 1; i <= j; i++)
                {
                    sum += r[i, j] * qtb[i];
                }

                int l = ipvt[j];
                wa1[j] = sum / diag[l];
            }

            double gnorm = EuclidianNorm(n, wa1);
            double paru = gnorm / delta;
            if (paru == 0.0)
            {
                paru = dwarf / Math.Min(delta, 0.1);
            }

            // if the input par lies outside of the interval (parl,paru),
            // set par to the closer endpoint.
            par = Math.Max(par, parl);
            par = Math.Min(par, paru);
            if (par == 0.0)
            {
                par = gnorm / dxnorm;
            }

            // beginning of an iteration.
            do
            {
                iter++;

                // evaluate the function at the current value of par.
                if (par == 0.0)
                {
                    par = Math.Max(dwarf, 1.0e-3 * paru);
                }

                double temp = Math.Sqrt(par);
                for (int j = 1; j <= n; j++)
                {
                    wa1[j] = temp * diag[j];
                }

                QRSolve(n, r, ipvt, wa1, qtb, x, sdiag, wa2);
                for (int j = 1; j <= n; j++)
                {
                    wa2[j] = diag[j] * x[j];
                }

                dxnorm = EuclidianNorm(n, wa2);
                temp = fp;
                fp = dxnorm - delta;

                // if the function is small enough, accept the current value
                // of par. also test for the exceptional cases where parl
                // is zero or the number of iterations has reached 10.
                if (Math.Abs(fp) <= 0.1 * delta ||
                    (parl == 0.0 && fp <= temp && temp < 0.0) ||
                    iter == 10)
                {
                    break;
                }

                // compute the newton correction.
                for (int j = 1; j <= n; j++)
                {
                    int l = ipvt[j];
                    wa1[j] = diag[l] * (wa2[l] / dxnorm);
                }

                for (int j = 1; j <= n; j++)
                {
                    wa1[j] = wa1[j] / sdiag[j];
                    temp = wa1[j];
                    int jp1 = j + 1;
                    if (n >= jp1)
                    {
                        for (int i = jp1; i <= n; i++)
                        {
                            wa1[i] = wa1[i] - (r[i, j] * temp);
                        }
                    }
                }

                temp = EuclidianNorm(n, wa1);
                double parc = ((fp / delta) / temp) / temp;

                // depending on the sign of the function, update parl or paru.
                if (fp > 0.0)
                {
                    parl = Math.Max(parl, par);
                }

                if (fp < 0.0)
                {
                    paru = Math.Min(paru, par);
                }

                // compute an improved estimate for par.
                par = Math.Max(parl, par + parc);

                // end of an iteration.
            }
            while (true);
        }

        // termination.
        if (iter == 0)
        {
            par = 0.0;
        }
    }

    private static void QRSolve(
        int n,
        NumericMatrix<double> r,
        NumericVector<int> ipvt,
        NumericVector<double> diag,
        NumericVector<double> qtb,
        NumericVector<double> x,
        NumericVector<double> sdiag,
        NumericVector<double> wa)
    {
        // subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
        // integer n,ldr
        // integer ipvt(n)
        // double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
        // c     **********
        // c
        // c     subroutine qrsolv
        // c
        // c     given an m by n matrix a, an n by n diagonal matrix d,
        // c     and an m-vector b, the problem is to determine an x which
        // c     solves the system
        // c
        // c           a*x = b ,     d*x = 0 ,
        // c
        // c     in the least squares sense.
        // c
        // c     this subroutine completes the solution of the problem
        // c     if it is provided with the necessary information from the
        // c     qr factorization, with column pivoting, of a. that is, if
        // c     a*p = q*r, where p is a permutation matrix, q has orthogonal
        // c     columns, and r is an upper triangular matrix with diagonal
        // c     elements of nonincreasing magnitude, then qrsolv expects
        // c     the full upper triangle of r, the permutation matrix p,
        // c     and the first n components of (q transpose)*b. the system
        // c     a*x = b, d*x = 0, is then equivalent to
        // c
        // c                  t       t
        // c           r*z = q *b ,  p *d*p*z = 0 ,
        // c
        // c     where x = p*z. if this system does not have full rank,
        // c     then a least squares solution is obtained. on output qrsolv
        // c     also provides an upper triangular matrix s such that
        // c
        // c            t   t               t
        // c           p *(a *a + d*d)*p = s *s .
        // c
        // c     s is computed within qrsolv and may be of separate interest.
        // c
        // c     the subroutine statement is
        // c
        // c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
        // c
        // c     where
        // c
        // c       n is a positive integer input variable set to the order of r.
        // c
        // c       r is an n by n array. on input the full upper triangle
        // c         must contain the full upper triangle of the matrix r.
        // c         on output the full upper triangle is unaltered, and the
        // c         strict lower triangle contains the strict upper triangle
        // c         (transposed) of the upper triangular matrix s.
        // c
        // c       ldr is a positive integer input variable not less than n
        // c         which specifies the leading dimension of the array r.
        // c
        // c       ipvt is an integer input array of length n which defines the
        // c         permutation matrix p such that a*p = q*r. column j of p
        // c         is column ipvt(j) of the identity matrix.
        // c
        // c       diag is an input array of length n which must contain the
        // c         diagonal elements of the matrix d.
        // c
        // c       qtb is an input array of length n which must contain the first
        // c         n elements of the vector (q transpose)*b.
        // c
        // c       x is an output array of length n which contains the least
        // c         squares solution of the system a*x = b, d*x = 0.
        // c
        // c       sdiag is an output array of length n which contains the
        // c         diagonal elements of the upper triangular matrix s.
        // c
        // c       wa is a work array of length n.
        // c
        // c     subprograms called
        // c
        // c       fortran-supplied ... dabs,dsqrt
        // c
        // c     argonne national laboratory. minpack project. march 1980.
        // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
        // c
        // c     **********
        // integer i,j,jp1,k,kp1,l,nsing
        // double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
        // data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/

        // copy r and (q transpose)*b to preserve input and initialize s.
        // in particular, save the diagonal elements of r in x.
        for (int j = 1; j <= n; j++)
        {
            for (int i = j; i <= n; i++)
            {
                r[i, j] = r[j, i];
            }

            x[j] = r[j, j];
            wa[j] = qtb[j];
        }

        // eliminate the diagonal matrix d using a givens rotation.
        for (int j = 1; j <= n; j++)
        {
            // prepare the row of d to be eliminated, locating the
            // diagonal element using p from the qr factorization.
            int l = ipvt[j];
            if (diag[l] != 0.0)
            {
                for (int k = j; k <= n; k++)
                {
                    sdiag[k] = 0.0;
                }

                sdiag[j] = diag[l];

                // the transformations to eliminate the row of d
                // modify only a single element of (q transpose)*b
                // beyond the first n, which is initially zero.
                double qtbpj = 0.0;
                for (int k = j; k <= n; k++)
                {
                    // determine a givens rotation which eliminates the
                    // appropriate element in the current row of d.
                    if (sdiag[k] != 0.0)
                    {
                        double sin, cos;
                        if (Math.Abs(r[k, k]) < Math.Abs(sdiag[k]))
                        {
                            double cotan = r[k, k] / sdiag[k];
                            sin = 0.5 / Math.Sqrt(0.25 + (0.25 * cotan * cotan));
                            cos = sin * cotan;
                        }
                        else
                        {
                            double tan = sdiag[k] / r[k, k];
                            cos = 0.5 / Math.Sqrt(0.25 + (0.25 * tan * tan));
                            sin = cos * tan;
                        }

                        // compute the modified diagonal element of r and
                        // the modified element of ((q transpose)*b,0).
                        r[k, k] = (cos * r[k, k]) + (sin * sdiag[k]);
                        double temp = (cos * wa[k]) + (sin * qtbpj);
                        qtbpj = (-sin * wa[k]) + (cos * qtbpj);
                        wa[k] = temp;

                        // accumulate the tranformation in the row of s.
                        int kp1 = k + 1;
                        if (n >= kp1)
                        {
                            for (int i = kp1; i <= n; i++)
                            {
                                temp = (cos * r[i, k]) + (sin * sdiag[i]);
                                sdiag[i] = (-sin * r[i, k]) + (cos * sdiag[i]);
                                r[i, k] = temp;
                            }
                        }
                    }
                }
            }

            // store the diagonal element of s and restore
            // the corresponding diagonal element of r.
            sdiag[j] = r[j, j];
            r[j, j] = x[j];
        }

        // solve the triangular system for z. if the system is
        // singular, then obtain a least squares solution.
        int nsing = n;
        for (int j = 1; j <= n; j++)
        {
            if (sdiag[j] == 0.0 && nsing == n)
            {
                nsing = j - 1;
            }

            if (nsing < n)
            {
                wa[j] = 0.0;
            }
        }

        if (nsing >= 1)
        {
            for (int k = 1; k <= nsing; k++)
            {
                int j = nsing - k + 1;
                double sum = 0.0;
                int jp1 = j + 1;
                if (nsing >= jp1)
                {
                    for (int i = jp1; i <= nsing; i++)
                    {
                        sum += r[i, j] * wa[i];
                    }
                }

                wa[j] = (wa[j] - sum) / sdiag[j];
            }
        }

        // permute the components of z back to components of x.
        for (int j = 1; j <= n; j++)
        {
            int l = ipvt[j];
            x[l] = wa[j];
        }
    }

    private static void ForwardDifferenceJacobian(NLSFunction function, NumericVector<double> x, NumericVector<double> fvec, ref NumericMatrix<double> fjac)
    {
        // Based on the subroutine fdjac2 from MINPACK (see https://netlib.org/minpack/), ported to C#

        // subroutine fdjac2
        //
        // this subroutine computes a forward-difference approximation
        // to the m by n jacobian matrix associated with a specified
        // problem of m functions in n variables.

        // x is an input array of length n.

        // fvec is an input array of length m which must contain the
        // functions evaluated at x.

        // fjac is an output m by n array which contains the
        // approximation to the jacobian matrix evaluated at x.

        // argonne national laboratory. minpack project. march 1980.
        // burton s. garbow, kenneth e. hillstrom, jorge j. more
        int n = x.Length;
        int m = fvec.Length;
        NumericVector<double> wa = new(m);
        for (int j = 1; j <= n; j++)
        {
            double temp = x[j];
            double h = EpsForwardDifferenceJacobian * Math.Abs(temp);
            if (h == 0)
            {
                h = EpsForwardDifferenceJacobian;
            }

            x[j] = temp + h;
            function(x, ref wa);
            x[j] = temp;
            for (int i = 1; i <= m; i++)
            {
                fjac[i, j] = (wa[i] - fvec[i]) / h;
            }
        }
    }

    private static void QRFactorization(
        int m,
        int n,
        ref NumericMatrix<double> a,
        bool pivot,
        ref NumericVector<int> ipvt,
        ref NumericVector<double> rdiag,
        ref NumericVector<double> acnorm,
        ref NumericVector<double> wa)
    {
        // subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
        // integer m,n,lda,lipvt
        // integer ipvt(lipvt)
        // logical pivot
        // double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
        // c     **********
        // c
        // c     subroutine qrfac
        // c
        // c     this subroutine uses householder transformations with column
        // c     pivoting (optional) to compute a qr factorization of the
        // c     m by n matrix a. that is, qrfac determines an orthogonal
        // c     matrix q, a permutation matrix p, and an upper trapezoidal
        // c     matrix r with diagonal elements of nonincreasing magnitude,
        // c     such that a*p = q*r. the householder transformation for
        // c     column k, k = 1,2,...,min(m,n), is of the form
        // c
        // c                           t
        // c           i - (1/u(k))*u*u
        // c
        // c     where u has zeros in the first k-1 positions. the form of
        // c     this transformation and the method of pivoting first
        // c     appeared in the corresponding linpack subroutine.
        // c
        // c     the subroutine statement is
        // c
        // c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
        // c
        // c     where

        // m is a positive integer input variable set to the number
        // of rows of a.

        // n is a positive integer input variable set to the number
        // of columns of a.

        // a is an m by n array. on input a contains the matrix for
        // which the qr factorization is to be computed. on output
        // the strict upper trapezoidal part of a contains the strict
        // upper trapezoidal part of r, and the lower trapezoidal
        // part of a contains a factored form of q (the non-trivial
        // elements of the u vectors described above).

        // c       lda is a positive integer input variable not less than m
        // c         which specifies the leading dimension of the array a.
        // c
        // c       pivot is a logical input variable. if pivot is set true,
        // c         then column pivoting is enforced. if pivot is set false,
        // c         then no column pivoting is done.
        // c
        // c       ipvt is an integer output array of length lipvt. ipvt
        // c         defines the permutation matrix p such that a*p = q*r.
        // c         column j of p is column ipvt(j) of the identity matrix.
        // c         if pivot is false, ipvt is not referenced.
        // c
        // c       lipvt is a positive integer input variable. if pivot is false,
        // c         then lipvt may be as small as 1. if pivot is true, then
        // c         lipvt must be at least n.
        // c
        // c       rdiag is an output array of length n which contains the
        // c         diagonal elements of r.
        // c
        // c       acnorm is an output array of length n which contains the
        // c         norms of the corresponding columns of the input matrix a.
        // c         if this information is not needed, then acnorm can coincide
        // c         with rdiag.
        // c
        // c       wa is a work array of length n. if pivot is false, then wa
        // c         can coincide with rdiag.
        // c
        // c     subprograms called
        // c
        // c       minpack-supplied ... dpmpar,enorm
        // c
        // c       fortran-supplied ... dmax1,dsqrt,min0
        // c
        // c     argonne national laboratory. minpack project. march 1980.
        // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
        // c
        // c     **********
        // integer i,j,jp1,k,kmax,minmn
        // double precision ajnorm,epsmch,one,p05,sum,temp,zero
        // double precision dpmpar,enorm
        // data one,p05,zero /1.0d0,5.0d-2,0.0d0/
        // c

        // epsmch is the machine precision.
        double epsmch = MyMath.Epsilon;

        // compute the initial column norms and initialize several arrays.
        for (int j = 1; j <= n; j++)
        {
            acnorm[j] = EuclidianNorm(m, a.GetColumnVector(j));
            rdiag[j] = acnorm[j];
            wa[j] = rdiag[j];
            if (pivot)
            {
                ipvt[j] = j;
            }
        }

        // reduce a to r with householder transformations.
        int minmn = Math.Min(m, n);
        for (int j = 1; j <= minmn; j++)
        {
            if (pivot)
            {
                // bring the column of largest norm into the pivot position.
                int kmax = j;
                for (int k = j; k <= n; k++)
                {
                    if (rdiag[k] > rdiag[kmax])
                    {
                        kmax = k;
                    }
                }

                if (kmax != j)
                {
                    for (int i = 1; i <= m; i++)
                    {
                        (a[i, kmax], a[i, j]) = (a[i, j], a[i, kmax]);
                    }

                    rdiag[kmax] = rdiag[j];
                    wa[kmax] = wa[j];
                    (ipvt[kmax], ipvt[j]) = (ipvt[j], ipvt[kmax]);
                }
            }

            // compute the householder transformation to reduce the
            // j-th column of a to a multiple of the j-th unit vector.
            double ajnorm = EuclidianNorm(m - j + 1, a.GetColumnVector(j, j));
            if (ajnorm != 0.0)
            {
                if (a[j, j] < 0.0)
                {
                    ajnorm = -ajnorm;
                }

                for (int i = j; i <= m; i++)
                {
                    a[i, j] /= ajnorm;
                }

                a[j, j] = a[j, j] + 1.0;

                // apply the transformation to the remaining columns
                // and update the norms.
                int jp1 = j + 1;
                if (n >= jp1)
                {
                    for (int k = jp1; k <= n; k++)
                    {
                        double sum = 0.0;
                        for (int i = j; i <= m; i++)
                        {
                            sum += a[i, j] * a[i, k];
                        }

                        double temp = sum / a[j, j];
                        for (int i = j; i <= m; i++)
                        {
                            a[i, k] = a[i, k] - (temp * a[i, j]);
                        }

                        if (pivot && rdiag[k] != 0.0)
                        {
                            temp = a[j, k] / rdiag[k];
                            rdiag[k] = rdiag[k] * Math.Sqrt(Math.Max(0.0, 1.0 - (temp * temp)));
                            if (5.0e-2 * Math.Pow(rdiag[k] / wa[k], 2) <= epsmch)
                            {
                                rdiag[k] = EuclidianNorm(m - j, a.GetColumnVector(k, jp1));
                                wa[k] = rdiag[k];
                            }
                        }
                    }
                }
            }

            rdiag[j] = -ajnorm;
        }
    }
}