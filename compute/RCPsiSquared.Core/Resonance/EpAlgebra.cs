namespace RCPsiSquared.Core.Resonance;

/// <summary>Closed-form constants from the heuristic 2-level EP algebra (F86 Statement 1).
///
/// <para>For two adjacent rate channels at HD = 2k−1 and HD = 2k+1 with inter-channel
/// coupling J·g_eff, the effective Liouvillian has eigenvalues
/// λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²). The discriminant vanishes at the exceptional
/// point Q_EP = 2/g_eff, where the two eigenvalues coalesce at λ_± = −4γ₀·k.</para>
///
/// <para>The slowest mode k=1 sets a universal e-folding time t_peak = 1/(4γ₀),
/// independent of c, N, n, and bond position. Bit-exact verified against full block-L
/// numerics across all tested cases.</para>
/// </summary>
public static class EpAlgebra
{
    /// <summary>Universal EP timescale at the slowest channel pair: t_peak = 1/(4γ₀).</summary>
    public static double TPeak(double gammaZero)
    {
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero),
            $"γ₀ must be > 0; got {gammaZero}.");
        return 1.0 / (4.0 * gammaZero);
    }

    /// <summary>Q_EP = 2/g_eff, the dimensionless coupling at which the slowest pair
    /// coalesces. g_eff is the heuristic 2-level inter-channel coupling.</summary>
    public static double QEp(double gEff)
    {
        if (gEff <= 0) throw new ArgumentOutOfRangeException(nameof(gEff),
            $"g_eff must be > 0; got {gEff}.");
        return 2.0 / gEff;
    }

    /// <summary>Pre-EP (real) eigenvalue of the slowest pair k=1 in the 2-level model
    /// (relative to the trace midpoint −4γ₀).</summary>
    public static (double LamPlus, double LamMinus) SlowestPairEigenvalues(double gammaZero, double j, double gEff)
    {
        double discriminant = 4.0 * gammaZero * gammaZero - j * j * gEff * gEff;
        double centre = -4.0 * gammaZero;
        if (discriminant >= 0)
        {
            double sqrt = Math.Sqrt(discriminant);
            return (centre + sqrt, centre - sqrt);
        }
        // Post-EP: complex pair around centre. We return the real centre for both — the
        // imaginary component must be queried separately.
        return (centre, centre);
    }
}
