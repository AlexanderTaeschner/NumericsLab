// <copyright file="FitResultStatus.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

namespace NumericsLab.Fitting;

/// <summary>
/// Specifies the termination status of a nonlinear least squares fitting operation.
/// </summary>
public enum FitResultStatus
{
    /*
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
     */

    /// <summary>
    /// Improper input parameters.
    /// </summary>
    ImproperInputParameters = 0,

    /// <summary>
    /// Both actual and predicted relative reductions in the sum of squares are at most ftol.
    /// </summary>
    RelativeReductionsInSumOfSquaresAtMostFtoll = 1,

    /// <summary>
    /// Relative error between two consecutive iterates is at most xtol.
    /// </summary>
    RelativeErrorBetweenConsecutiveIteratesAtMostXtoll = 2,

    /// <summary>
    /// Conditions for RelativeReductionsInSumOfSquaresAtMostFtoll and RelativeErrorBetweenConsecutiveIteratesAtMostXtoll both hold.
    /// </summary>
    RelativeReductionsInSumOfSquaresAndRelativeErrorBetweenConsecutiveIteratesAtMostFtollAndXtoll = 3,

    /// <summary>
    /// The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.
    /// </summary>
    CosineOfAngleBetweenFvecAndJacobianColumnsAtMostGtol = 4,

    /// <summary>
    /// Number of calls to fcn has reached or exceeded maxfev.
    /// </summary>
    NumberOfCallsToFcnReachedOrExceededMaxfev = 5,

    /// <summary>
    /// ftol is too small. No further reduction in the sum of squares is possible.
    /// </summary>
    FtollIsTooSmall = 6,

    /// <summary>
    /// xtol is too small. No further improvement in the approximate solution x is possible.
    /// </summary>
    XtollIsTooSmall = 7,

    /// <summary>
    /// gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.
    /// </summary>
    GtolIsTooSmall = 8,
}