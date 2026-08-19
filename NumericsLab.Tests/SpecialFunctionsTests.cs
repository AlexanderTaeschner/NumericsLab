namespace NumericsLab.Tests;

public class SpecialFunctionsTests
{
    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the incomplete beta function.
    [InlineData(1.0, 1.0, 0.0, 0.0)]
    [InlineData(1.0, 1.0, 1.0, 1.0)]
    [InlineData(0.1, 0.1, 1.0, 1.0)]
    [InlineData(1.0, 1.0, 0.5, 0.5)]
    [InlineData(0.1, 1.0, 0.5, 0.9330329915368074160, 1e-15)]
    [InlineData(10.0, 1.0, 0.5, 0.0009765625000000000000)]
    [InlineData(50.0, 1.0, 0.5, 8.881784197001252323e-16)]
    [InlineData(1.0, 0.1, 0.5, 0.06696700846319258402)]
    [InlineData(1.0, 10.0, 0.5, 0.99902343750000000000)]
    [InlineData(1.0, 50.0, 0.5, 0.99999999999999911180)]
    [InlineData(1.0, 1.0, 0.1, 0.10)]
    [InlineData(1.0, 2.0, 0.1, 0.19)]
    [InlineData(1.0, 2.0, 0.9, 0.99)]
    [InlineData(50.0, 60.0, 0.5, 0.8309072939016694143, 1e-15)]
    [InlineData(90.0, 90.0, 0.5, 0.5, 1e-12)]
    [InlineData(500.0, 500.0, 0.6, 0.9999999999157549630)]
    [InlineData(5000.0, 5000.0, 0.4, 4.518543727260666383e-91)]
    [InlineData(5000.0, 5000.0, 0.6, 1.0, 1e-15)]
    [InlineData(5000.0, 2000.0, 0.6, 8.445388773903332659e-89)]
    public void IncompleteBetaFunctionTest(double a, double b, double x, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.IncompleteBetaIntegral(a, b, x);
        Assert.Equal(expected, result, tolerance);
    }

    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the incomplete gamma function (gsl_sf_gamma_inc_P_e).
    [InlineData(1e-100, 0.001, 1.0, 1e-13)]
    [InlineData(0.001, 0.001, 0.9936876467088602902, 1e-15)]
    [InlineData(0.001, 1.0, 0.9997803916424144436, 1e-15)]
    [InlineData(0.001, 10.0, 0.9999999958306921828)]
    [InlineData(1.0, 0.001, 0.0009995001666250083319)]
    [InlineData(1.0, 1.01, 0.6357810204284766802, 1e-15)]
    [InlineData(1.0, 10.0, 0.9999546000702375151)]
    [InlineData(10.0, 10.01, 0.5433207586693410570, 1e-15)]
    [InlineData(10.0, 20.0, 0.9950045876916924128)]
    [InlineData(1000.0, 1000.1, 0.5054666401440661753, 1e-12)]
    [InlineData(1000.0, 2000.0, 1.0)]
    [InlineData(34.0, 32.0, 0.3849626436463866776322932129, 1e-14)]
    [InlineData(37.0, 3.499999999999999289e+01, 0.3898035054195570860969333039, 1e-15)]
    [InlineData(10, 1e-16, 2.755731922398588814734648067e-167)]
    [InlineData(1263131.0, 1261282.3637, 0.04994777516935182963821362168, 1e-10)]
    [InlineData(1263131.0, 1263131.0, 0.500118321758657770672882362502514254, 1e-9)]
    [InlineData(100, 99.0, 0.4733043303994607, 1e-13)]
    [InlineData(200, 199.0, 0.4811585880878718, 1e-14)]
    [InlineData(5670, 4574, 3.063972328743934e-55)]
    public void IncompleteGammaIntegralTest(double a, double x, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.IncompleteGammaIntegral(a, x);
        Assert.Equal(expected, result, tolerance);
    }

    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the gamma function (gsl_sf_gamma_e).
    [InlineData(1.0 + 1.0 / 4096.0, 0.9998591371459403421, 1e-15)]
    [InlineData(1.0 + 1.0 / 32.0, 0.9829010992836269148, 1e-15)]
    [InlineData(2.0 + 1.0 / 256.0, 1.0016577903733583299)]
    [InlineData(9.0, 40320.0)]
    [InlineData(10.0, 362880.0)]
    [InlineData(100.0, 9.332621544394415268e+155, 1e141)]
    [InlineData(170.0, 4.269068009004705275e+304, 1e289)]
    [InlineData(171.0, 7.257415615307998967e+306, 1e292)]
    [InlineData(-10.5, -2.640121820547716316e-07)]
    [InlineData(-11.25, 6.027393816261931672e-08)]
    [InlineData(-1.0 + 1.0 / 65536.0, -65536.42280587818970)]
    public void GammaFunctionTest(double input, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.GammaFunction(input);
        Assert.Equal(expected, result, tolerance);
    }

    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the logarithmic gamma function (gsl_sf_lngamma_e).
    [InlineData(-0.1, 2.368961332728788655, 1e-15)]
    [InlineData(-1.0 / 256.0, 5.547444766967471595)]
    [InlineData(1.0e-08, 18.420680738180208905)]
    [InlineData(0.1, 2.252712651734205, 1e-15)]
    [InlineData(1.0 + 1.0 / 256.0, -0.0022422226599611501448)]
    [InlineData(2.0 + 1.0 / 256.0, 0.0016564177556961728692)]
    [InlineData(100.0, 359.1342053695753, 1e-12)]
    [InlineData(-1.0 - 1.0 / 65536.0, 11.090348438090047844, 1e-12)]
    [InlineData(-1.0 - 1.0 / 268435456.0, 19.408121054103474300)]
    [InlineData(-100.5, -364.9009683094273518, 1e-13)]
    [InlineData(-100 - 1.0 / 65536.0, -352.6490910117097874, 1e-13)]
    public void LogarithmicGammaFunctionTest(double input, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.LogarithmicGammaFunction(input);
        Assert.Equal(expected, result, tolerance);
    }

    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the error function (gsl_sf_erf_e)
    [InlineData(-10.0, -1.0000000000000000000)]
    [InlineData(0.5, 0.5204998778130465377)]
    [InlineData(1.0, 0.8427007929497148693, 1e-15)]
    [InlineData(10.0, 1.0000000000000000000)]
    public void ErrorFunctionTest(double input, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.ErrorFunction(input);
        Assert.Equal(expected, result, tolerance);
    }

    [Theory]
    // Test cases taken from GNU Scientific Library (GSL) tests for the complementary error function (gsl_sf_erfc_e)
    [InlineData(-10.0, 2.0)]
    [InlineData(-5.0000002, 1.9999999999984625433)]
    [InlineData(-5.0, 1.9999999999984625402)]
    [InlineData(-1.0, 1.8427007929497148693)]
    [InlineData(-0.5, 1.5204998778130465377)]
    [InlineData(1.0, 0.15729920705028513066)]
    [InlineData(3.0, 0.000022090496998585441373)]
    [InlineData(7.0, 4.183825607779414399e-23)]
    [InlineData(10.0, 2.0884875837625447570e-45)]
    public void ComplementaryErrorFunctionTest(double input, double expected, double tolerance = 1.0e-16)
    {
        double result = SpecialFunctions.ComplementaryErrorFunction(input);
        Assert.Equal(expected, result, tolerance);
    }
}
