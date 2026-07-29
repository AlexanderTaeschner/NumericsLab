using NumericsLab.RootFinding;

namespace NumericsLab.Tests;

/// <summary>
/// Tests for the Brent root-finding algorithm. Taken from Source of the C code: https://netlib.org/c/brent.shar.
/// </summary>
public class RootFindingBrentTests
{
    [Fact]
    public void Test1()
    {
        int counter = 0;
        double result = Brent.FindRoot((double x, double[] parameters) => { counter++; return ((Math.Pow(x, 2) - 2.0) * x) - 5.0; }, [], 2.0, 3.0, 0.0);
        Assert.Equal(2.0945514815423265, result);
        Assert.Equal(10, counter);
    }

    [Fact]
    public void Test2()
    {
        int counter = 0;
        double result = Brent.FindRoot((double x, double[] parameters) => { counter++; return Math.Cos(x) - x; }, [], 2.0, 3.0, 0.0);
        Assert.Equal(2.0, result);
        Assert.Equal(51, counter);
    }

    [Fact]
    public void Test3()
    {
        int counter = 0;
        double result = Brent.FindRoot((double x, double[] parameters) => { counter++; return Math.Cos(x) - x; }, [], -1.0, 3.0, 0.0);
        Assert.Equal(0.73908513321516067, result);
        Assert.Equal(10, counter);
    }

    [Fact]
    public void Test4()
    {
        int counter = 0;
        double result = Brent.FindRoot((double x, double[] parameters) => { counter++; return Math.Sin(x) - x; }, [], -1.0, 3.0, 0.0);
        Assert.Equal(-1.6437373573790809E-08, result);
        Assert.Equal(57, counter);
    }

    [Fact]
    public void Test5()
    {
        int counter = 0;
        bool success = Brent.TryFindRoot((double x, double[] parameters) => { counter++; return Math.Cos(x) - x; }, [], -1.0, 3.0, 0.0, 9, out double result);
        Assert.True(success);
        Assert.Equal(0.73908513321516067, result);
        Assert.Equal(10, counter);
    }

    [Fact]
    public void Test6()
    {
        int counter = 0;
        bool success = Brent.TryFindRoot((double x, double[] parameters) => { counter++; return Math.Cos(x) - x; }, [], -1.0, 3.0, 0.0, 8, out double result);
        Assert.False(success);
        Assert.Equal(0.73908513321526159, result);
        Assert.Equal(9, counter);
    }
}
