// <copyright file="FitResult.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

using NumericsLab.Internals;

namespace NumericsLab.Fitting;

/// <summary>
/// Represents the result of a nonlinear least squares fitting operation, including fitted parameters, residuals, and statistical information about the fit.
/// </summary>
public class FitResult
{
    private readonly double[,] _covarianceMatrix;
    private readonly double[] _parameterStandardDeviation;

    internal FitResult(int info, NumericVector<double> x, NumericVector<double> fvec, NumericMatrix<double> fjac, NumericVector<int> ipvt)
    {
        Info = (FitResultStatus)info;
        Parameters = x.AsReadOnlyList();
        Residuals = fvec.AsReadOnlyList();
        EuclidianNormOfResiduals = NonlinearLeastSquares.EuclidianNorm(fvec.Length, fvec);
        ChiSquared = EuclidianNormOfResiduals * EuclidianNormOfResiduals;
        NumberOfDegreesOfFreedom = fvec.Length - x.Length;
        EstimatedStandardDeviationOfFit = Math.Sqrt(ChiSquared / NumberOfDegreesOfFreedom);
        _covarianceMatrix = CalculateCovarianceMatrix(x.Length, fjac, ipvt);
        _parameterStandardDeviation = CalculateParameterStandardDeviation(_covarianceMatrix);
    }

    /// <summary>
    /// Gets the status of the fitting operation, indicating whether it was successful or if any issues were encountered.
    /// </summary>
    public FitResultStatus Info { get; }

    /// <summary>
    /// Gets the fitted parameters resulting from the nonlinear least squares fitting operation.
    /// </summary>
    public IReadOnlyList<double> Parameters { get; }

    /// <summary>
    /// Gets the residuals resulting from the nonlinear least squares fitting operation.
    /// </summary>
    public IReadOnlyList<double> Residuals { get; }

    /// <summary>
    /// Gets the Euclidian norm of the residuals.
    /// </summary>
    public double EuclidianNormOfResiduals { get; }

    /// <summary>
    /// Gets the chi-squared value of the fit.
    /// </summary>
    public double ChiSquared { get; }

    /// <summary>
    /// Gets the number of degrees of freedom of the fit.
    /// </summary>
    public int NumberOfDegreesOfFreedom { get; }

    /// <summary>
    /// Gets the estimated standard deviation of the fit.
    /// </summary>
    public double EstimatedStandardDeviationOfFit { get; }

    /// <summary>
    /// Gets the covariance matrix of the fitted parameters.
    /// </summary>
    /// <param name="errorScaled">If true, scales the covariance matrix by the estimated standard deviation of the fit.</param>
    /// <returns>The covariance matrix.</returns>
    public IReadOnlyMatrix<double> GetCovarianceMatrix(bool errorScaled)
    {
        if (errorScaled)
        {
            double scale = EstimatedStandardDeviationOfFit * EstimatedStandardDeviationOfFit;
            double[,] scaledCovarianceMatrix = new double[_covarianceMatrix.GetLength(0), _covarianceMatrix.GetLength(1)];
            for (int i = 0; i < _covarianceMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < _covarianceMatrix.GetLength(1); j++)
                {
                    scaledCovarianceMatrix[i, j] = _covarianceMatrix[i, j] * scale;
                }
            }

            return new ReadOnlyMatrix<double>(scaledCovarianceMatrix);
        }
        else
        {
            return new ReadOnlyMatrix<double>(_covarianceMatrix);
        }
    }

    /// <summary>
    /// Gets the standard deviation of the fitted parameters.
    /// </summary>
    /// <param name="errorScaled">If true, scales the standard deviation by the estimated standard deviation of the fit.</param>
    /// <returns>The standard deviation of the fitted parameters.</returns>
    public IReadOnlyList<double> GetParameterStandardDeviation(bool errorScaled)
    {
        if (errorScaled)
        {
            double scale = EstimatedStandardDeviationOfFit;
            double[] scaledParameterStandardDeviation = new double[_parameterStandardDeviation.Length];
            for (int i = 0; i < _parameterStandardDeviation.Length; i++)
            {
                scaledParameterStandardDeviation[i] = _parameterStandardDeviation[i] * scale;
            }

            return scaledParameterStandardDeviation;
        }
        else
        {
            return _parameterStandardDeviation;
        }
    }

    private static double[,] CalculateCovarianceMatrix(int n, NumericMatrix<double> fjac, NumericVector<int> ipvt)
    {
        double[,] covarianceMatrix = new double[n, n];
        Covariance(n, fjac, ipvt, MyMath.EpsilonSqrt);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                covarianceMatrix[i, j] = fjac[i + 1, j + 1]; // Adjusting for 1-based indexing in the original Fortran code
            }
        }

        return covarianceMatrix;
    }

    private static double[] CalculateParameterStandardDeviation(double[,] covarianceMatrix)
    {
        double[] parameterStandardDeviation = new double[covarianceMatrix.GetLength(0)];
        for (int i = 0; i < covarianceMatrix.GetLength(0); i++)
        {
            parameterStandardDeviation[i] = Math.Sqrt(covarianceMatrix[i, i]);
        }

        return parameterStandardDeviation;
    }

    private static void Covariance(int n, NumericMatrix<double> r, NumericVector<int> ipvt, double tol)
    {
        // SUBROUTINE COVAR(N,R,LDR,IPVT,TOL,WA)
        // INTEGER N,LDR
        // INTEGER IPVT(N)
        // DOUBLE PRECISION TOL
        // DOUBLE PRECISION R(LDR,N),WA(N)
        // C     **********
        //
        // C     SUBROUTINE COVAR
        // C     GIVEN AN M BY N MATRIX A, THE PROBLEM IS TO DETERMINE
        // C     THE COVARIANCE MATRIX CORRESPONDING TO A, DEFINED AS
        // C
        // C                     T
        // C           INVERSE (A *A) .
        // C
        // C     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
        // C     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
        // C     QR FACTORIZATION, WITH COLUMN PIVOTING, OF A. THAT IS, IF
        // C     A*P = Q*R, WHERE P IS A PERMUTATION MATRIX, Q HAS ORTHOGONAL
        // C     COLUMNS, AND R IS AN UPPER TRIANGULAR MATRIX WITH DIAGONAL
        // C     ELEMENTS OF NONINCREASING MAGNITUDE, THEN COVAR EXPECTS
        // C     THE FULL UPPER TRIANGLE OF R AND THE PERMUTATION MATRIX P.
        // C     THE COVARIANCE MATRIX IS THEN COMPUTED AS
        // C
        // C                      T     T
        // C           P*INVERSE(R *R)*P
        // C
        // C     IF A IS NEARLY RANK DEFICIENT, IT MAY BE DESIRABLE TO COMPUTE
        // C     THE COVARIANCE MATRIX CORRESPONDING TO THE LINEARLY INDEPENDENT
        // C     COLUMNS OF A. TO DEFINE THE NUMERICAL RANK OF A, COVAR USES
        // C     THE TOLERANCE TOL. IF L IS THE LARGEST INTEGER SUCH THAT
        // C
        // C     ABS(R(L,L)) .GT. TOL*ABS(R(1,1)) ,
        // C
        // C     THEN COVAR COMPUTES THE COVARIANCE MATRIX CORRESPONDING TO
        // C     THE FIRST L COLUMNS OF R. FOR K GREATER THAN L, COLUMN
        // C     AND ROW IPVT(K) OF THE COVARIANCE MATRIX ARE SET TO ZERO.
        // C
        // C     THE SUBROUTINE STATEMENT IS
        // C
        // C       SUBROUTINE COVAR(N,R,LDR,IPVT,TOL,WA)
        // C
        // C     WHERE
        // C
        // C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
        // C
        // C       R IS AN N BY N ARRAY. ON INPUT THE FULL UPPER TRIANGLE MUST
        // C         CONTAIN THE FULL UPPER TRIANGLE OF THE MATRIX R. ON OUTPUT
        // C         R CONTAINS THE SQUARE SYMMETRIC COVARIANCE MATRIX.
        // C
        // C       LDR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
        // C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY R.
        // C
        // C       IPVT IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH DEFINES THE
        // C         PERMUTATION MATRIX P SUCH THAT A*P = Q*R. COLUMN J OF P
        // C         IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
        // C
        // C       TOL IS A NONNEGATIVE INPUT VARIABLE USED TO DEFINE THE
        // C         NUMERICAL RANK OF A IN THE MANNER DESCRIBED ABOVE.
        // C
        // C       WA IS A WORK ARRAY OF LENGTH N.
        // C
        // C     SUBPROGRAMS CALLED
        // C
        // C       FORTRAN-SUPPLIED ... DABS
        // C
        // C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. AUGUST 1980.
        // C     BURTON S. GARBOW , KENNETH E. HILLSTROM, JORGE J. MORE
        // C
        // C     **********
        // INTEGER I,II,J,JJ,K,KM1,L
        // LOGICAL SING
        // DOUBLE PRECISION ONE,TEMP,TOLR,ZERO
        // DATA ONE,ZERO /1.0D0,0.0D0/

        // FORM THE INVERSE OF'R IN THE FULL UPPER TRIANGLE OF R.
        NumericVector<double> wa = new(n);
        double tolr = tol * Math.Abs(r[1, 1]);
        int l = 0;
        for (int k = 1; k <= n; k++)
        {
            if (Math.Abs(r[k, k]) <= tolr)
            {
                break;
            }

            r[k, k] = 1.0 / r[k, k];
            int km1 = k - 1;
            if (km1 >= 1)
            {
                for (int j = 1; j <= km1; j++)
                {
                    double temp = r[k, k] * r[j, k];
                    r[j, k] = 0.0;
                    for (int i = 1; i <= j; i++)
                    {
                        r[i, k] = r[i, k] - (temp * r[i, j]);
                    }
                }
            }

            l = k;
        }

        // FORM THE FULL UPPER TRIANGLE OF THE INVERSE OF (R TRANSPOSE)*R
        // IN THE FULL UPPER TRIANGLE OF R.
        if (l >= 1)
        {
            for (int k = 1; k <= l; k++)
            {
                double temp;
                int km1 = k - 1;
                if (km1 >= 1)
                {
                    for (int j = 1; j <= km1; j++)
                    {
                        temp = r[j, k];
                        for (int i = 1; i <= j; i++)
                        {
                            r[i, j] = r[i, j] + (temp * r[i, k]);
                        }
                    }
                }

                temp = r[k, k];
                for (int i = 1; i <= k; i++)
                {
                    r[i, k] = temp * r[i, k];
                }
            }
        }

        // FORM THE FULL LOWER TRIANGLE OF THE COVARIANCE MATRIX
        // IN THE STRICT LOWER TRIANGLE OF R AND IN WA.
        for (int j = 1; j <= n; j++)
        {
            int jj = ipvt[j];
            bool sing = j > l;
            for (int i = 1; i <= j; i++)
            {
                if (sing)
                {
                    r[i, j] = 0.0;
                }

                int ii = ipvt[i];
                if (ii > jj)
                {
                    r[ii, jj] = r[i, j];
                }

                if (ii < jj)
                {
                    r[jj, ii] = r[i, j];
                }
            }

            wa[jj] = r[j, j];
        }

        // SYMMETRIZE THE COVARIANCE MATRIX IN R.
        for (int j = 1; j <= n; j++)
        {
            for (int i = 1; i <= j; i++)
            {
                r[i, j] = r[j, i];
            }

            r[j, j] = wa[j];
        }

        // C
        // C     LAST CARD OF SUBROUTINE COVAR.
        // C
        // END
    }
}