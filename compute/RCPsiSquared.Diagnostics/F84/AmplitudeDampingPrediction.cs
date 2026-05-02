using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Diagnostics.F84;

/// <summary>F84 closed-form prediction and inversion of the F81 violation under thermal
/// amplitude damping (cooling σ⁻ + heating σ⁺).
///
/// Theorem F84 (proved in PROOF_F84_AMPLITUDE_DAMPING.md) generalises F82:
/// <code>
///   ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)
/// </code>
///
/// <para>The F81 violation depends only on the net cooling rate Δγ_l = γ_↓_l − γ_↑_l. At
/// thermal equilibrium γ_↓ = γ_↑ (detailed balance), the violation vanishes even though
/// both channels are active — the F81 violation isolates the vacuum-fluctuation contribution
/// (temperature-independent).</para>
///
/// <para>Pauli-channel dissipators D[Z], D[X], D[Y] do NOT contribute to the violation
/// (Pauli-Channel Cancellation Lemma in PROOF_F84): they are Π²-symmetric. F84 violation is
/// exclusive to σ⁻/σ⁺ channels.</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F84 entry.</para>
/// </summary>
public static class AmplitudeDampingPrediction
{
    /// <summary>‖D_{AmplDamp, odd}‖_F = √(Σ (γ_↓ − γ_↑)²) · 2^(N−1).
    /// <paramref name="gammaPumpPerSite"/> defaults to zero (pure cooling, F82 special case).</summary>
    public static double PredictViolation(ChainSystem chain,
        IReadOnlyList<double> gammaT1PerSite, IReadOnlyList<double>? gammaPumpPerSite = null)
    {
        if (gammaT1PerSite.Count != chain.N)
            throw new ArgumentException($"gamma_t1 list length {gammaT1PerSite.Count} != N {chain.N}");
        var pump = gammaPumpPerSite ?? Enumerable.Repeat(0.0, chain.N).ToArray();
        if (pump.Count != chain.N)
            throw new ArgumentException($"gamma_pump list length {pump.Count} != N {chain.N}");
        double sumSq = 0;
        for (int l = 0; l < chain.N; l++)
        {
            double diff = gammaT1PerSite[l] - pump[l];
            sumSq += diff * diff;
        }
        return Math.Sqrt(sumSq) * Math.Pow(2, chain.N - 1);
    }

    /// <summary>Inverse: |Δγ|_RMS = √(Σ (γ_↓ − γ_↑)² / N) = f81_violation / (√N · 2^(N−1)).
    /// At T=0 (γ_↑=0) this reduces to <see cref="F82.T1DissipatorPrediction.EstimateRmsT1FromViolation"/>.</summary>
    public static double EstimateRmsNetCoolingFromViolation(ChainSystem chain, double f81Violation)
    {
        if (f81Violation < 0)
            throw new ArgumentOutOfRangeException(nameof(f81Violation), $"must be non-negative; got {f81Violation}");
        return f81Violation / (Math.Sqrt(chain.N) * Math.Pow(2, chain.N - 1));
    }
}
