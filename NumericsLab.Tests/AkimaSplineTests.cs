using NumericsLab.Interpolation;

namespace NumericsLab.Tests;

public class AkimaSplineTests
{
    [Fact]
    public void TestAkimaSpline()
    {
        // Sample data points
        double[] x = { 0, 1, 2, 3, 4, 5 };
        double[] y = { 0, 1, 2, 3, 4, 5 };
        // Create Akima spline
        CubicSpline spline = CubicSpline.CreateAkimaSpline(x, y);
        // Test interpolation at various points
        double[] testPoints = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
        double[] expectedValues = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
        for (int i = 0; i < testPoints.Length; i++)
        {
            double interpolatedValue = spline.Interpolate(testPoints[i]);
            Assert.Equal(expectedValues[i], interpolatedValue);
        }
    }

    [Fact]
    public void TestAkimaSpline2()
    {
        // Sample data points
        double[] x = { -4, -3, -2, -1, 0, 1, 2, 3, 4 };
        double[] y = { 16, 9, 4, 1, 0, 1, 4, 9, 16 };
        // Create Akima spline
        CubicSpline spline = CubicSpline.CreateAkimaSpline(x, y);
        // Test interpolation at various points
        double[] testPoints = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };
        for (int i = 0; i < testPoints.Length; i++)
        {
            double interpolatedValue = spline.Interpolate(testPoints[i]);
            double expectedValue = testPoints[i] * testPoints[i];
            Assert.Equal(expectedValue, interpolatedValue);
        }
    }

    [Fact]
    public void TestAkimaSpline3()
    {
        // Sample data points
        double[] x = { -4, -3, -2, -1, 0, 1, 2, 3, 4 };
        double[] y = { -64, -27, -8, -1, 0, 1, 8, 27, 64 };
        // Create Akima spline
        CubicSpline spline = CubicSpline.CreateAkimaSpline(x, y);
        // Test interpolation at various points
        double[] testPoints = { -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 };
        double[] expectedValues = { -64, -43.025, -27, -15.475, -8, -3.375, -1, -0.5, 0, 0.5, 1, 3.375, 8, 15.4750, 27, 43.025, 64 };
        for (int i = 0; i < testPoints.Length; i++)
        {
            double interpolatedValue = spline.Interpolate(testPoints[i]);
            Assert.Equal(expectedValues[i], interpolatedValue, 1e-14);
        }
    }
}
