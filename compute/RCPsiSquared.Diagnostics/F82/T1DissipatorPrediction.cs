using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Diagnostics.F82;

/// <summary>F82 closed-form prediction and inversion of the F81 violation under T1 amplitude damping.
///
/// Theorem F82 (proved in PROOF_F82_T1_DISSIPATOR_CORRECTION.md):
/// <code>
///   ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1)
/// </code>
///
/// <para>The F81 violation magnitude. Hamiltonian-independent and γ_z-independent (Master
/// Lemma propagates through F82). For uniform γ_T1: simplifies to γ_T1 · √N · 2^(N−1).</para>
///
/// <para>Use case: predict the F81 violation expected for a given hardware T1 profile,
/// compare with <see cref="F81.PiDecomposition"/> output to validate the F82 model.</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F82 entry.</para>
/// </summary>
public static class T1DissipatorPrediction
{
    /// <summary>‖D_{T1, odd}‖_F = √(Σ γ_T1²) · 2^(N−1).</summary>
    public static double PredictViolation(ChainSystem chain, IReadOnlyList<double> gammaT1PerSite)
    {
        if (gammaT1PerSite.Count != chain.N)
            throw new ArgumentException($"gamma_t1 list length {gammaT1PerSite.Count} != N {chain.N}");
        double sumSq = gammaT1PerSite.Sum(g => g * g);
        return Math.Sqrt(sumSq) * Math.Pow(2, chain.N - 1);
    }

    /// <summary>Inverse: γ_T1, RMS = f81_violation / (√N · 2^(N−1)). For uniform per-site T1
    /// returns the actual γ_T1; for non-uniform returns the root-mean-square.</summary>
    public static double EstimateRmsT1FromViolation(ChainSystem chain, double f81Violation)
    {
        if (f81Violation < 0)
            throw new ArgumentOutOfRangeException(nameof(f81Violation), $"must be non-negative; got {f81Violation}");
        return f81Violation / (Math.Sqrt(chain.N) * Math.Pow(2, chain.N - 1));
    }
}
