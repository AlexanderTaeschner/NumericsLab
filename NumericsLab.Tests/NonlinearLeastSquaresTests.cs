using NumericsLab.Fitting;

namespace NumericsLab.Tests;

public class NonlinearLeastSquaresTests
{
    [Fact]
    public void Test1()
    {
        static double ExampleFunction(double x, IReadOnlyList<double> parameters)
        {
            double u = x;
            double v = 16.0 - x;
            double w = Math.Min(u, v);
            return parameters[0] + (u / ((v * parameters[1]) + (w * parameters[2])));
        }

        double[] xdata = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        double[] ydata = [0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39];
        double[] initialParams = [1.0, 1.0, 1.0];

        FitResult result = NonlinearLeastSquares.Fit(ExampleFunction, xdata, ydata, initialParams);
        Assert.Equal(FitResultStatus.RelativeReductionsInSumOfSquaresAtMostFtoll, result.Info);
        Assert.Equal(0.9063596e-01,result.EuclidianNormOfResiduals,9); // number from file06 of https://netlib.org/minpack/ex/
        Assert.Equal(0.8241057e-01, result.Parameters[0],7); // number from file06 of https://netlib.org/minpack/ex/
        Assert.Equal(0.1133037e+01, result.Parameters[1],6); // number from file06 of https://netlib.org/minpack/ex/
        Assert.Equal(0.2343695e+01, result.Parameters[2],4); // number from file06 of https://netlib.org/minpack/ex/
        Assert.Equal(0.0261643480503795, result.EstimatedStandardDeviationOfFit, 9); // number from gnuplot fit
        IReadOnlyMatrix<double> covarianceMatrix = result.GetCovarianceMatrix(true);
        Assert.Equal(3, covarianceMatrix.NumberOfRows);
        Assert.Equal(3, covarianceMatrix.NumberOfColumns);
        Assert.Equal(0.000153120891007205, covarianceMatrix[0, 0], 7); // number from gnuplot fit
        Assert.Equal(0.00287232957086124, covarianceMatrix[0, 1], 5); // number from gnuplot fit
        Assert.Equal(0.00287232957086124, covarianceMatrix[1, 0], 5); // number from gnuplot fit
        Assert.Equal(0.0949673479311945, covarianceMatrix[1, 1], 3); // number from gnuplot fit
        Assert.Equal(-0.00265990524494864, covarianceMatrix[0, 2], 5); // number from gnuplot fit
        Assert.Equal(-0.00265990524494864, covarianceMatrix[2, 0], 5); // number from gnuplot fit
        Assert.Equal(-0.0911752995397185, covarianceMatrix[2, 1], 3); // number from gnuplot fit
        Assert.Equal(-0.0911752995397185, covarianceMatrix[1, 2], 3); // number from gnuplot fit
        Assert.Equal(0.0879981606994613, covarianceMatrix[2, 2], 3); // number from gnuplot fit
        IReadOnlyList<double> parameterStdDev = result.GetParameterStandardDeviation(true);
        Assert.Equal(3, parameterStdDev.Count);
        Assert.Equal(0.0123742026412697, parameterStdDev[0], 7); // number from gnuplot fit
        Assert.Equal(0.308167726946211, parameterStdDev[1], 3); // number from gnuplot fit
        Assert.Equal(0.296644839327202, parameterStdDev[2], 2); // number from gnuplot fit
    }
}
