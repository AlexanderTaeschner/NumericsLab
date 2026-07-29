using NumericsLab;

namespace NumericsLab.Tests;

public class DifferentiateTests
{
    [Fact]
    public void TestFirstDerivative()
    {
        // Test the first derivative of f(x) = x^2 at x = 2, which should be 4.
        double result = Differentiate.FirstDerivative((double x, double[] parameters) => Math.Pow(x, 2), [], 2.0);
        Assert.Equal(4.0, result, tolerance: 1e-11);
    }

    [Fact]
    public void TestFirstDerivative2()
    {
        // Test the first derivative of f(x) = sin(x) at x = 0, which should be 1.
        double result = Differentiate.FirstDerivative((double x, double[] parameters) => Math.Sin(x), [], 0.0);
        Assert.Equal(1.0, result, tolerance: 1e-11);
    }
}
