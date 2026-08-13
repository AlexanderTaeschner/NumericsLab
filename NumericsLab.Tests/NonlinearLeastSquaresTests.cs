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

        var result = NonlinearLeastSquares.Fit(ExampleFunction, xdata, ydata, initialParams);
    }
}
