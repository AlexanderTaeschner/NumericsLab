using System;
using System.Collections.Generic;
using System.Text;

namespace NumericsLab.Fitting;

public static class NonlinearLeastSquares
{
    /// <summary>
    /// Represents a model function for nonlinear least squares fitting.
    /// </summary>
    /// <param name="x">The independent variable.</param>
    /// <param name="parameters">The parameters of the model function.</param>
    /// <returns>The value of the model function at the given x with the specified parameters.</returns>
    public delegate double ModelFunction(double x, double[] parameters);

    // F_i(parameters) = y_i - g(x_i, parameters); g is the model function
    private delegate double[] NLSFunction(double[] parameters);

    public static double[] Fit(ModelFunction model, double[] xdata, double[] ydata, double[] initialParams, int maxIterations = 1000, double tolerance = 1e-6)
    {
        NLSFunction function = (double[] parameters) =>
        {
            double[] residuals = new double[xdata.Length];
            for (int i = 0; i < xdata.Length; i++)
            {
                residuals[i] = ydata[i] - model(xdata[i], parameters);
            }

            return residuals;
        };

        return LevenbergMarquardt(function, initialParams, maxIterations, tolerance);
    }

    private static double[] LevenbergMarquardt(NLSFunction function, double[] initialParams, int maxIterations, double tolerance)
    {
// KI generated code below - TO BE REPLACED
        int n = initialParams.Length;
        double[] parameters = (double[])initialParams.Clone();
        double lambda = 0.001;
        for (int iteration = 0; iteration < maxIterations; iteration++)
        {
            double[] residuals = function(parameters);
            double[][] jacobian = ComputeJacobian(function, parameters);
            // Compute the normal equations
            double[][] jTj = Multiply(Transpose(jacobian), jacobian);
            double[] jTr = Multiply(Transpose(jacobian), residuals);
            // Add damping factor
            for (int i = 0; i < n; i++)
            {
                jTj[i][i] += lambda;
            }
            // Solve for parameter update
            double[] deltaParams = SolveLinearSystem(jTj, jTr);
            // Update parameters
            for (int i = 0; i < n; i++)
            {
                parameters[i] += deltaParams[i];
            }
            // Check for convergence
            if (Norm(deltaParams) < tolerance)
            {
                break;
            }
        }
        return parameters;
    }

    private static double[][] ComputeJacobian(NLSFunction function, double[] parameters)
    {
        int n = parameters.Length;
        double[] f0 = function(parameters);
        int m = f0.Length;
        double[][] jacobian = new double[m][];
        for (int i = 0; i < m; i++) jacobian[i] = new double[n];
        double eps = 1e-8;
        for (int j = 0; j < n; j++)
        {
            double temp = parameters[j];
            parameters[j] = temp + eps;
            double[] f1 = function(parameters);
            parameters[j] = temp;
            for (int i = 0; i < m; i++)
            {
                jacobian[i][j] = (f1[i] - f0[i]) / eps;
            }
        }
        return jacobian;
    }

    private static double Norm(double[] vector)
    {
        double sum = 0.0;
        for (int i = 0; i < vector.Length; i++)
        {
            sum += vector[i] * vector[i];
        }
        return Math.Sqrt(sum);
    }

    private static double[][] Transpose(double[][] matrix)
    {
        int rows = matrix.Length;
        int cols = matrix[0].Length;
        var result = new double[cols][];
        for (int i = 0; i < cols; i++)
        {
            result[i] = new double[rows];
            for (int j = 0; j < rows; j++)
            {
                result[i][j] = matrix[j][i];
            }
        }
        return result;
    }

    private static double[][] Multiply(double[][] a, double[][] b)
    {
        int rows = a.Length;
        int cols = b[0].Length;
        int n = b.Length;
        var result = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            result[i] = new double[cols];
            for (int j = 0; j < cols; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

    private static double[] Multiply(double[][] a, double[] b)
    {
        int rows = a.Length;
        int cols = a[0].Length;
        var result = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                result[i] += a[i][j] * b[j];
            }
        }
        return result;
    }

    private static double[] SolveLinearSystem(double[][] a, double[] b)
    {
        int n = a.Length;
        var x = new double[n];
        var lu = (double[][])a.Clone();
        var indx = new int[n];
        LUDecomposition(lu, indx);
        LUSolve(lu, indx, b, x);
        return x;
    }

    private static void LUDecomposition(double[][] a, int[] indx)
    {
        int n = a.Length;
        var vv = new double[n];
        for (int i = 0; i < n; i++)
        {
            double big = 0.0;
            for (int j = 0; j < n; j++)
            {
                double temp = Math.Abs(a[i][j]);
                if (temp > big) big = temp;
            }
            if (big == 0.0) throw new ArgumentException("Singular matrix");
            vv[i] = 1.0 / big;
        }
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < j; i++)
            {
                double sum = a[i][j];
                for (int k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
                a[i][j] = sum;
            }
            double big = 0.0;
            int imax = j;
            for (int i = j; i < n; i++)
            {
                double sum = a[i][j];
                for (int k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
                a[i][j] = sum;
                double dum = vv[i] * Math.Abs(sum);
                if (dum >= big)
                {
                    big = dum;
                    imax = i;
                }
            }
            if (j != imax)
            {
                var temp = a[imax];
                a[imax] = a[j];
                a[j] = temp;
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if (a[j][j] == 0.0) a[j][j] = 1e-20;
            if (j != n - 1)
            {
                double dum = 1.0 / a[j][j];
                for (int i = j + 1; i < n; i++) a[i][j] *= dum;
            }
        }
    }

    private static void LUSolve(double[][] a, int[] indx, double[] b, double[] x)
    {
        int n = a.Length;
        b.CopyTo(x, 0);
        int ii = -1;
        for (int i = 0; i < n; i++)
        {
            int ip = indx[i];
            double sum = x[ip];
            x[ip] = x[i];
            if (ii != -1)
            {
                for (int j = ii; j <= i - 1; j++) sum -= a[i][j] * x[j];
            }
            else if (sum != 0.0)
            {
                ii = i;
            }
            x[i] = sum;
        }
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = x[i];
            for (int j = i + 1; j < n; j++) sum -= a[i][j] * x[j];
            x[i] = sum / a[i][i];
        }
    }

    private static double Norm(double[] v)
    {
        double sum = 0.0;
        for (int i = 0; i < v.Length; i++) sum += v[i] * v[i];
        return Math.Sqrt(sum);
    }
}
