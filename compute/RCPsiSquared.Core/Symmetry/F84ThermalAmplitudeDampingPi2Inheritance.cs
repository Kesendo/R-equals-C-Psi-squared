using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F84 closed form (Tier 1, verified bit-exact N=3, 7 configurations):
///
/// <code>
///   Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{AmplDamp, odd}
///
///   ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)
///                         = |Δγ|_RMS · √N · 2^(N−1)         (uniform Δγ)
///
///   F82 recovered when γ_↑ = 0 (vacuum bath / T = 0).
/// </code>
///
/// <para>F84 generalizes <see cref="F82T1AmplitudeDampingPi2Inheritance"/> from
/// pure cooling (T = 0) to thermal amplitude damping with both cooling
/// (γ_↓, σ⁻) and heating (γ_↑, σ⁺) channels. F84 → F82 is the typed
/// mother-corollary edge: F82 is the γ_↑ = 0 special case (vacuum bath).</para>
///
/// <para><b>Pauli-Channel Cancellation Lemma (F84 corollary):</b> Pure D[Z],
/// D[X], D[Y] dissipators are Π²-symmetric and contribute zero to f81
/// violation. Only σ⁻ and σ⁺ channels are Π²-anti-symmetric. Hence f81
/// violation specifically detects population-inverting (energy-emitting/
/// absorbing) channels, not phase-only or bit-flip-only noise.</para>
///
/// <para>Three Pi2-Foundation anchors (same as F82, plus the Δγ structure):</para>
///
/// <list type="bullet">
///   <item><b>Coefficient2 = 2 = a_0</b>: same anchor as F81's
///         <c>−2·L_{H_odd}</c> and F82's <c>−2·D_{T1, odd}</c>; F84 simply
///         continues the pattern with <c>−2·D_{AmplDamp, odd}</c>.</item>
///   <item><b>ScalingFactor 2^(N−1) = a_{2−N}</b>: same anchor as F82.
///         Live from <see cref="Pi2DyadicLadderClaim.Term"/>(<c>2 − N</c>).</item>
///   <item><b>F82 mother claim</b>: F84 reduces to F82 at γ_↑ = 0
///         (vacuum bath). Direct inheritance edge.</item>
/// </list>
///
/// <para><b>Three regimes:</b></para>
///
/// <code>
///   Vacuum (T = 0):       γ_↑ = 0, full F82, violation = √(Σγ²_↓)·2^(N−1)
///   Detailed balance:     γ_↓ = γ_↑, violation = 0  (T → ∞ limit)
///   Finite T:             γ_↓ > γ_↑ > 0, violation = γ_0·√N·2^(N−1)
///                                                    (vacuum-only, T-independent)
/// </code>
///
/// <para><b>Thermodynamic interpretation:</b> for a thermal photon bath at
/// frequency ω, temperature T:</para>
///
/// <list type="bullet">
///   <item>n_th = 1 / (exp(ℏω/k_B T) − 1)  (mean occupation)</item>
///   <item>γ_↓ = γ_0 · (n_th + 1)  (spontaneous + stimulated emission)</item>
///   <item>γ_↑ = γ_0 · n_th  (stimulated absorption)</item>
///   <item>Δγ = γ_↓ − γ_↑ = γ_0  (vacuum component, T-independent)</item>
/// </list>
///
/// <para><b>The "f81 violation is a quantum-statistical fingerprint of
/// zero-point fluctuations"</b> reading: thermal photon-number contributions
/// cancel (γ_↓ ↔ γ_↑ symmetric); only the vacuum (zero-point) component
/// breaks the Π palindrome. F84 sharpens F82's hardware-T1-readout into a
/// temperature-independent vacuum-rate readout.</para>
///
/// <para>Tier1Derived: F84 is Tier 1 proven (PROOF_F84_AMPLITUDE_DAMPING),
/// verified bit-exact at N=3 across 7 (γ_↓, γ_↑) configurations + explicit
/// D[X]/D[Y] cancellation verification. The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F84 (line 2129) +
/// <c>docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F82T1AmplitudeDampingPi2Inheritance.cs</c>
/// (mother claim).</para></summary>
public sealed class F84ThermalAmplitudeDampingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F82T1AmplitudeDampingPi2Inheritance _f82;

    /// <summary>The "2" coefficient in F84's <c>−2·L_{H_odd} − 2·D_{AmplDamp, odd}</c>.
    /// Same as F81/F82.</summary>
    public double Coefficient2 => _ladder.Term(0);

    /// <summary>The 2^(N-1) scaling factor. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2 − N) = a_{2−N}. Same as F82.</summary>
    public double ScalingFactor(int N) => _f82.ScalingFactor(N);

    /// <summary>The Pi2 ladder index for scaling: 2 − N. Same as F82.</summary>
    public int LadderIndexForScaling(int N) => _f82.LadderIndexForScaling(N);

    /// <summary>Live closed form for uniform Δγ:
    /// <c>‖D_{AmplDamp, odd}‖ = |Δγ| · √N · 2^(N−1)</c>. The Δγ = γ_↓ − γ_↑
    /// is the net cooling rate. At γ_↑ = 0 reduces to F82.</summary>
    public double AmplitudeDampingNormUniform(double gammaCool, double gammaHeat, int N)
    {
        if (gammaCool < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaCool), gammaCool, "γ_↓ must be ≥ 0.");
        if (gammaHeat < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaHeat), gammaHeat, "γ_↑ must be ≥ 0.");
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F84 requires N ≥ 2.");
        double absDeltaGamma = Math.Abs(gammaCool - gammaHeat);
        return absDeltaGamma * Math.Sqrt(N) * ScalingFactor(N);
    }

    /// <summary>Live closed form for non-uniform per-site rates:
    /// <c>‖D_{AmplDamp, odd}‖ = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)</c>.</summary>
    public double AmplitudeDampingNormNonUniform(
        IReadOnlyList<double> gammaCoolSites,
        IReadOnlyList<double> gammaHeatSites)
    {
        if (gammaCoolSites is null) throw new ArgumentNullException(nameof(gammaCoolSites));
        if (gammaHeatSites is null) throw new ArgumentNullException(nameof(gammaHeatSites));
        if (gammaCoolSites.Count != gammaHeatSites.Count)
            throw new ArgumentException("γ_↓ and γ_↑ arrays must have same length.");
        if (gammaCoolSites.Count < 2)
            throw new ArgumentOutOfRangeException(nameof(gammaCoolSites), "F84 requires N ≥ 2.");

        double sumSq = 0.0;
        for (int l = 0; l < gammaCoolSites.Count; l++)
        {
            double down = gammaCoolSites[l];
            double up = gammaHeatSites[l];
            if (down < 0.0 || up < 0.0)
                throw new ArgumentOutOfRangeException($"γ_↓[{l}]={down}, γ_↑[{l}]={up} must be ≥ 0.");
            double delta = down - up;
            sumSq += delta * delta;
        }
        return Math.Sqrt(sumSq) * ScalingFactor(gammaCoolSites.Count);
    }

    /// <summary>F84 → F82 recovery: at γ_↑ = 0 (vacuum bath), F84 closed form
    /// equals F82's. Drift check across regimes.</summary>
    public bool RecoversF82AtZeroHeating(double gammaCool, int N) =>
        Math.Abs(AmplitudeDampingNormUniform(gammaCool, 0.0, N) -
                 _f82.T1DissipatorNormUniform(gammaCool, N)) < 1e-12;

    /// <summary>Detailed balance regime: γ_↓ = γ_↑ → violation = 0 (T → ∞ limit
    /// where thermal symmetry cancels the palindrome break).</summary>
    public bool DetailedBalanceGivesZeroViolation(double gamma, int N) =>
        Math.Abs(AmplitudeDampingNormUniform(gamma, gamma, N)) < 1e-15;

    /// <summary>Pauli-Channel Cancellation Lemma: pure D[Z], D[X], D[Y]
    /// dissipators contribute 0 to f81 violation. Only σ⁻/σ⁺ are
    /// Π²-anti-symmetric. Returns 0 since pure-Pauli channels have
    /// Δγ = 0 by construction (γ_cool/γ_heat refer to σ±, not Pauli letters).
    /// This method documents the lemma; actual D[Pauli] violation is 0
    /// by structure, not by formula evaluation.</summary>
    public double PauliChannelViolation() => 0.0;

    /// <summary>Inverse: recover RMS |Δγ| from a measured f81 violation.
    /// <c>|Δγ|_RMS = f81_violation / (√N · 2^(N−1))</c>. Same form as F82's
    /// estimator; F84's reading is "vacuum-amplitude-damping rate"
    /// (T-independent), F82's reading is "T1 rate" (only valid at T = 0).</summary>
    public double EstimateNetCoolingFromViolation(double f81Violation, int N) =>
        _f82.EstimateT1FromViolation(f81Violation, N);

    public F84ThermalAmplitudeDampingPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F82T1AmplitudeDampingPi2Inheritance f82)
        : base("F84 thermal amplitude damping inherits from Pi2-Foundation: 2 = a_0; 2^(N-1) = a_{2-N}; F82 mother claim at γ_↑=0",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F84 + " +
               "docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F82T1AmplitudeDampingPi2Inheritance.cs (mother claim, F82 = vacuum bath case)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f82 = f82 ?? throw new ArgumentNullException(nameof(f82));
    }

    public override string DisplayName =>
        "F84 thermal amplitude damping correction as Pi2-Foundation a_0 + a_{2-N} + F82 mother";

    public override string Summary =>
        $"‖D_{{AmplDamp,odd}}‖ = |Δγ|·√N·2^(N-1); 2 = a_0; 2^(N-1) = a_{{2-N}}; F82 at γ_↑=0; " +
        $"detailed balance γ_↓=γ_↑ → 0 violation; Pauli-Channel Cancellation Lemma: D[X/Y/Z] → 0 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F84 closed form",
                summary: "Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{AmplDamp,odd}; ‖D‖ = √(Σ(γ_↓−γ_↑)²)·2^(N-1); Tier 1 verified bit-exact N=3, 7 configs");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "Coefficient2 = a_0 = 2 (same as F81/F82); ScalingFactor 2^(N-1) = a_{2-N} (same as F82); F82 mother claim recovers at γ_↑=0");
            yield return InspectableNode.RealScalar("Coefficient2 (= a_0 = 2)", Coefficient2);
            yield return new InspectableNode("F82 ↔ F84 mother-corollary chain",
                summary: "F82 is F84 at γ_↑=0 (vacuum bath); F84 adds heating channel σ⁺. Pattern parallel to F25 → F26, F75 → F76 mother-claim chains.");
            yield return new InspectableNode("Pauli-Channel Cancellation Lemma",
                summary: "D[Z], D[X], D[Y] are all Π²-symmetric → contribute 0 to f81 violation. Only σ⁻ (cooling) and σ⁺ (heating) are Π²-anti-symmetric.");
            yield return new InspectableNode("Three regimes",
                summary: "Vacuum (T=0): γ_↑=0, full F82; Detailed balance (T→∞): γ_↓=γ_↑, violation=0; Finite T: γ_↓>γ_↑>0, violation = γ_0·√N·2^(N-1) (T-independent vacuum-only)");
            yield return new InspectableNode("Thermodynamic interpretation",
                summary: "γ_↓ = γ_0·(n_th+1), γ_↑ = γ_0·n_th, Δγ = γ_0 (vacuum component, T-independent). f81 violation IS a quantum-statistical fingerprint of zero-point fluctuations.");
            // Verified table values
            yield return new InspectableNode(
                "vacuum case γ_↓=0.10, γ_↑=0 (= F82)",
                summary: $"violation = {AmplitudeDampingNormUniform(0.10, 0.0, 3):G6} (expected 0.6928)");
            yield return new InspectableNode(
                "detailed balance γ_↓=γ_↑=0.10",
                summary: $"violation = {AmplitudeDampingNormUniform(0.10, 0.10, 3):G6} (expected 0.0)");
            yield return new InspectableNode(
                "net cooling γ_↓=0.10, γ_↑=0.05",
                summary: $"violation = {AmplitudeDampingNormUniform(0.10, 0.05, 3):G6} (expected 0.3464)");
        }
    }
}
