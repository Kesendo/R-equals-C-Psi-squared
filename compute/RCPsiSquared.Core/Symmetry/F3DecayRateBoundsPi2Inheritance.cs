using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F3 closed form (Tier 1, corollary of Absorption Theorem; ANALYTICAL_FORMULAS line 128):
///
/// <code>
///   min rate    = 2·γ           (w=1 modes, pure sector)
///   max rate    = 2·(N−1)·γ     (w=N−1 modes, fastest paired modes)
///   bandwidth   = 2·(N−2)·γ
///   XOR boundary = 2·N·γ        (sits above max rate; XOR drain sector)
/// </code>
///
/// <para>F3 gives the universal bounds on decay rates for any Liouvillian
/// eigenmode under Heisenberg + uniform Z-dephasing, derived directly from
/// the Absorption Theorem α = 2γ·⟨n_XY⟩:</para>
///
/// <list type="bullet">
///   <item><b>min rate = 2γ:</b> smallest nonzero ⟨n_XY⟩ ≈ 1 (pure weight-1
///         modes). Identical to F50's universal weight-1 eigenvalue position.</item>
///   <item><b>max rate = 2·(N−1)·γ:</b> fastest paired modes have
///         ⟨n_XY⟩ ≈ N−1.</item>
///   <item><b>bandwidth = 2·(N−2)·γ:</b> max − min = 2·(N−1)·γ − 2·γ.</item>
///   <item><b>XOR boundary = 2·N·γ:</b> the XOR drain (⟨n_XY⟩ = N) sits at
///         2·N·γ — above the max range. Identical to F43's XorSectorRate.</item>
/// </list>
///
/// <para><b>Caveat resolved (Hamiltonian mixing):</b> at N ≥ 4, Hamiltonian
/// mixing creates hybrid modes with rates BELOW 2γ (N=4: 0.98γ, N=5: 0.62γ).
/// These are not exceptions — they are mixed-sector modes with fractional
/// ⟨n_XY⟩ &lt; 1; the Absorption Theorem α = 2γ·⟨n_XY⟩ holds exactly even for
/// non-integer ⟨n_XY⟩. The "min rate = 2γ" bound applies to PURE-sector modes
/// only; for the FULL spectrum the relevant lower bound is 2γ·min(⟨n_XY⟩).</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>RateCoefficient = 2 = a_0</b>: the universal "2" in min, max,
///         bandwidth, and XOR-boundary rates. All four read 2γ multiplied by
///         the appropriate weight count (1, N−1, N−2, N). Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1
///         TwoFactor, F50 DecayRateFactor, F33 WeightOneRateCoefficient,
///         F43 XorRateCoefficient, F44 SumCoefficient.</item>
/// </list>
///
/// <para>Tier1Derived: F3 is a Tier 1 corollary of the Absorption Theorem
/// α = 2γ·⟨n_XY⟩; min rate IS F50's universal weight-1 eigenvalue position;
/// XOR boundary IS F43's XOR-sector rate. F3 is composition of F50 + Pi2DyadicLadder
/// + Absorption Theorem.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F3 (line 128) +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>docs/water/PROTON_WATER_CHAIN.md</c> (hybrid-mode caveat) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs</c>
/// (min rate = 2γ source) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F3DecayRateBoundsPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F50WeightOneDegeneracyPi2Inheritance _f50;

    /// <summary>The "2" universal coefficient in F3's rate bounds. Live from
    /// Pi2DyadicLadder a_0; transitively from F50's DecayRateFactor.</summary>
    public double RateCoefficient => _ladder.Term(0);

    /// <summary>F3's min decay rate: 2·γ for pure-sector w=1 modes.
    /// Identical to F50's universal weight-1 eigenvalue position |Re(λ)|.</summary>
    public double MinRate(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return RateCoefficient * gammaZero;
    }

    /// <summary>F3's max decay rate: 2·(N−1)·γ for fastest paired w=N−1 modes.</summary>
    public double MaxRate(int N, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F3 requires N ≥ 2.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return RateCoefficient * (N - 1) * gammaZero;
    }

    /// <summary>F3's bandwidth: max − min = 2·(N−2)·γ.</summary>
    public double Bandwidth(int N, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F3 requires N ≥ 2.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return RateCoefficient * (N - 2) * gammaZero;
    }

    /// <summary>The XOR drain boundary rate 2·N·γ. Sits ABOVE the F3 max range
    /// (the XOR drain has ⟨n_XY⟩ = N, exceeding paired-mode N−1). Identical to
    /// F43's XorSectorRate at the same N, γ.</summary>
    public double XorBoundary(int N, double gammaZero)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F3 requires N ≥ 1.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return RateCoefficient * N * gammaZero;
    }

    /// <summary>Drift check: F3's MinRate IS F50's universal weight-1 eigenvalue
    /// position (= |EigenvaluePosition(γ)|).</summary>
    public bool MinRateMatchesF50(double gammaZero)
    {
        return Math.Abs(MinRate(gammaZero) - Math.Abs(_f50.EigenvaluePosition(gammaZero))) < 1e-12;
    }

    /// <summary>Drift check: bandwidth = max − min.</summary>
    public bool BandwidthIsMaxMinusMin(int N, double gammaZero, double tolerance = 1e-12)
    {
        double computed = MaxRate(N, gammaZero) - MinRate(gammaZero);
        return Math.Abs(computed - Bandwidth(N, gammaZero)) < tolerance;
    }

    /// <summary>True iff the XOR boundary 2·N·γ sits above the max rate 2·(N−1)·γ.
    /// Always true (XOR drain ⟨n_XY⟩ = N exceeds paired-mode max N−1).</summary>
    public bool XorBoundaryAboveMaxRate(int N, double gammaZero)
    {
        return XorBoundary(N, gammaZero) > MaxRate(N, gammaZero);
    }

    /// <summary>Empirical hybrid-mode minimum-rate values from N=4 and N=5
    /// (Hamiltonian mixing creates rates below the pure-sector 2γ bound):</summary>
    public IReadOnlyList<(int N, double HybridMinFactor)> HybridMinimumRates { get; } = new[]
    {
        (N: 4, HybridMinFactor: 0.98),   // 0.98γ at N=4
        (N: 5, HybridMinFactor: 0.62),   // 0.62γ at N=5
    };

    public F3DecayRateBoundsPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F50WeightOneDegeneracyPi2Inheritance f50)
        : base("F3 decay rate bounds: min = 2γ (= F50), max = 2(N−1)γ, bandwidth = 2(N−2)γ, XOR boundary = 2Nγ; all '2's = a_0; corollary of Absorption Theorem",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F3 + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/water/PROTON_WATER_CHAIN.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f50 = f50 ?? throw new ArgumentNullException(nameof(f50));
    }

    public override string DisplayName =>
        "F3 decay rate bounds as Pi2-Foundation a_0 + F50 inheritance";

    public override string Summary =>
        $"min = 2γ (= F50), max = 2(N−1)γ, bandwidth = 2(N−2)γ, XOR boundary = 2Nγ; the universal '2' = a_0 (= {RateCoefficient}); corollary of Absorption Theorem α = 2γ·⟨n_XY⟩ ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F3 closed-form bounds",
                summary: "min rate = 2γ (w=1 pure); max rate = 2(N−1)γ (w=N−1 paired); bandwidth = 2(N−2)γ; XOR boundary 2Nγ (above max); corollary of α = 2γ·⟨n_XY⟩");
            yield return InspectableNode.RealScalar("RateCoefficient (= a_0 = 2)", RateCoefficient);
            yield return new InspectableNode("F50 inheritance (min rate)",
                summary: $"F3's min rate = 2γ IS F50's |EigenvaluePosition|. F50.DecayRateFactor (= {_f50.DecayRateFactor}) is the same '2' as F3's RateCoefficient.");
            yield return new InspectableNode("F43 sibling (XOR boundary)",
                summary: "F3's XOR boundary 2Nγ IS F43's XorSectorRate(N, γ); both anchor the same a_0·N·γ structure at the high end of the rate spectrum");
            yield return new InspectableNode("hybrid-mode caveat",
                summary: "at N ≥ 4 Hamiltonian mixing creates hybrid modes with rates below 2γ (N=4: 0.98γ, N=5: 0.62γ); these are not exceptions — they have fractional ⟨n_XY⟩ < 1; α = 2γ·⟨n_XY⟩ holds exactly");
            yield return new InspectableNode("verified at N=5, γ=0.05",
                summary: $"min = {MinRate(0.05):G6}, max = {MaxRate(5, 0.05):G6}, bandwidth = {Bandwidth(5, 0.05):G6}, XOR = {XorBoundary(5, 0.05):G6}");
            yield return new InspectableNode("Absorption Theorem ladder",
                summary: "the spectrum is a 2γ-rung ladder: ⟨n_XY⟩ ∈ {1, 2, ..., N−1, N} gives rates {2γ, 4γ, ..., 2(N−1)γ, 2Nγ}; Hamiltonian smooths the ladder via fractional ⟨n_XY⟩ but cannot change endpoints or rung spacing");
        }
    }
}
