using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Helpers for evaluating the Duhamel integral via the spectral contour formula.
///
/// <para>For a non-Hermitian generator L = R · diag(λ) · R⁻¹, the integral
/// ∫_0^t e^(L(t-s)) V e^(L s) ds appears as ∂ρ/∂J in the J-perturbation. In the eigenbasis
/// it reduces to F_{j,k}(t) · X_{j,k} where X = R⁻¹ V R and F_{j,k}(t) = (e^(λ_k t) − e^(λ_j t)) / (λ_k − λ_j),
/// with the degenerate fallback F_{j,j}(t) = t · e^(λ_j t).</para>
/// </summary>
public static class SpectralContourIntegral
{
    /// <summary>F_{j,k}(t) = (exp(λ_k·t) − exp(λ_j·t)) / (λ_k − λ_j), with the degenerate
    /// fallback t·exp(λ_j·t) when |λ_k − λ_j| is below <paramref name="tolerance"/>.</summary>
    public static Complex Element(Complex lamJ, Complex lamK, Complex eJ, Complex eK,
        double t, double tolerance = 1e-10)
    {
        Complex diff = lamK - lamJ;
        if (diff.Magnitude > tolerance)
            return (eK - eJ) / diff;
        return t * eJ;
    }

    /// <summary>Build the full F matrix at a single t given eigenvalues and the precomputed
    /// e_i = exp(λ_i · t) vector.</summary>
    public static Complex[,] FullMatrix(IReadOnlyList<Complex> eigenvalues, Complex[] expLambdaT,
        double t, double tolerance = 1e-10)
    {
        int n = eigenvalues.Count;
        var f = new Complex[n, n];
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                f[j, k] = Element(eigenvalues[j], eigenvalues[k], expLambdaT[j], expLambdaT[k], t, tolerance);
            }
        }
        return f;
    }
}
