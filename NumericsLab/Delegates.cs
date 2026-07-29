namespace NumericsLab;

/// <summary>
/// Function to be used for finding its root.
/// </summary>
/// <param name="argument">Argument of the function.</param>
/// <param name="parameters">Array of constant parameters for the function.</param>
/// <returns>The functional value.</returns>
public delegate double OneDimensionalFunction(double argument, double[] parameters);