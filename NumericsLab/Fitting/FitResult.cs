// <copyright file="FitResult.cs" company="Alexander Täschner">
// Copyright (c) Alexander Täschner. All rights reserved.
// </copyright>

namespace NumericsLab.Fitting;

public class FitResult
{
    internal FitResult(int info, Internals.NumericVector<double> x, Internals.NumericVector<double> fvec)
    {
        Info = (FitResultStatus)info;
        Parameters = x.AsReadOnlyList();
        ////Residuals = fvec.AsReadOnlyList();
        EuclidianNormOfResiduals = NonlinearLeastSquares.EuclidianNorm(fvec.Length, fvec);
    }

    public FitResultStatus Info { get; }

    public IReadOnlyList<double> Parameters { get; }

    ////public IReadOnlyList<double> Residuals { get; }

    public double EuclidianNormOfResiduals { get; set; }
}