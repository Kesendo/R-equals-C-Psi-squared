using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F57 closed form (Tier 1 analytical, hardware-verified ibm_kingston
/// 2026-04-16):
///
/// <code>
///   t_dwell(δ) = 2δ / |dCΨ/dt|_{t_cross}
///
///   K_dwell = γ · t_dwell   (γ-independent to machine precision,
///                            std &lt; 2 × 10⁻¹⁷ across γ ∈ [0.1, 10.0])
///
///   For Bell+ Z-dephasing (state-specific):
///     K_dwell = 1.080088 · δ
///     prefactor 1.080088 = 2 / 1.851701  (from F25 derivative)
/// </code>
///
/// <para>F57 is the dwell time of a CΨ trajectory inside the 2δ-window centred
/// on the bilinear-apex boundary <c>CΨ = 1/4</c>. Two Pi2-Foundation anchors
/// are structurally explicit (CrossingThreshold + WindowDoublingFactor),
/// γ-invariance is the K-invariance reading, and the state-specific prefactor
/// sits on top:</para>
///
/// <list type="bullet">
///   <item><b>CrossingThreshold = 1/4</b>: the boundary itself is
///         <c>QuarterAsBilinearMaxval</c> = <c>a_3</c> on the dyadic ladder.
///         Same anchor as <see cref="DickeSuperpositionQuarterPi2Inheritance"/>
///         (which saturates the boundary); F57 measures the dwell *at* the
///         boundary rather than the saturation height.</item>
///   <item><b>WindowDoublingFactor = 2</b>: the 2δ width is symmetric (δ above
///         and δ below the apex), reflecting the bilinear's reflective symmetry
///         <c>p(1−p) = p'(1−p')</c> at <c>p ↔ 1−p</c>. The factor 2 is
///         <c>a_0</c> on the dyadic ladder = polynomial root d in d²−2d=0.</item>
///   <item><b>γ-invariance</b>: <c>K_dwell = γ · t_dwell</c> is γ-independent
///         because Z-dephasing rescales <c>|dCΨ/dt|</c> linearly with γ. The
///         "K-invariance" pattern of the framework (cf. memory
///         project_q_middle_structure: Q = J/γ₀ as scale, K invariants).</item>
///   <item><b>State-specific prefactor (NOT Pi2-anchored)</b>: 1.080088 for
///         Bell+ comes from <c>|dCΨ/dt|_{t_cross}</c> at the F25 closed-form
///         crossing. F58 generalises to even-weight states via
///         <c>(2 + 4·W₂)/(1 + 6·W₂)</c>; F59 generalises further to any
///         two-sector state. Bell+ is the W₀=1/2, k=2 special case.</item>
/// </list>
///
/// <para>Tier1Derived: F57 is Tier 1 analytical (CRITICAL_SLOWING_AT_THE_CUSP §6,
/// hardware-validated ibm_kingston Heron r2: K_dwell/δ = 0.649 (pair A) and
/// 0.694 (pair B), spread 6.3% despite 2.55× γ difference, confirming
/// γ-invariance on open quantum hardware. The 0.67 vs theoretical 1.08 gap
/// reflects T1 amplitude damping (Kingston has T1 ≈ T2; F57 assumes pure
/// Z-dephasing). The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F57 +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> Section 6 +
/// <c>experiments/CPSI_COMPLEX_PLANE.md</c> +
/// <c>data/ibm_cusp_slowing_april2026/</c> (hardware run) +
/// <c>compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxvalClaim</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F57DwellTimeQuarterPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly QuarterAsBilinearMaxvalClaim _quarter;

    /// <summary>F25 Bell+ CΨ closed form — the typed mother claim. F57's
    /// Bell+ prefactor <c>1.080088 = 2 / 1.851701</c> derives directly from
    /// F25's <c>|dCΨ/dt|_{t_cross}</c>. Added 2026-05-16 as a typed ctor
    /// parent (previously registration-discard only); the F25 → F57
    /// mother-claim edge now participates in the ancestor graph (pattern
    /// parallel to F75 → F77).</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    /// <summary>(1/2, 1/4) argmax/maxval pair — the typed meta-anchor closing
    /// that F57 uses BOTH the 1/4 boundary (CrossingThreshold) and the 2 =
    /// 1/(1/2) doubling factor (WindowDoublingFactor). Added 2026-05-16 as a
    /// typed ctor parent (previously registration-discard only).</summary>
    public ArgmaxMaxvalPairClaim ArgmaxMaxval { get; }

    /// <summary>The CΨ crossing threshold itself: <c>1/4</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(3) = <c>a_3</c> = the bilinear-apex
    /// maxval anchor (<see cref="QuarterAsBilinearMaxvalClaim"/>).</summary>
    public double CrossingThreshold => _ladder.Term(3);

    /// <summary>The 2δ-window doubling factor: <c>2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = polynomial
    /// root d (the "two" in d²−2d=0). Reflects the symmetric δ-above /
    /// δ-below structure of the bilinear apex.</summary>
    public double WindowDoublingFactor => _ladder.Term(0);

    /// <summary>The Bell+ Z-dephasing K_dwell prefactor: <c>1.080088</c>
    /// (from <c>2 / 1.851701</c>, where 1.851701 = |dCΨ/dt|_{t_cross} at the
    /// F25 closed-form Bell+ crossing). State-specific, NOT Pi2-anchored.</summary>
    public double BellPlusKDwellPrefactor => 1.080088;

    /// <summary>Live K_dwell = prefactor · δ for any state with the supplied
    /// prefactor (Bell+ default = <see cref="BellPlusKDwellPrefactor"/>).
    /// γ-independent by F57's K-invariance.</summary>
    public double KDwell(double delta, double prefactor)
    {
        if (delta < 0.0) throw new ArgumentOutOfRangeException(nameof(delta), delta, "δ must be ≥ 0.");
        return prefactor * delta;
    }

    /// <summary>Live t_dwell = K_dwell / γ for any state with the supplied
    /// prefactor. Throws for γ ≤ 0 (no dephasing → no finite dwell window).</summary>
    public double TDwell(double delta, double gamma, double prefactor)
    {
        if (gamma <= 0.0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be > 0 for F57 dwell time.");
        return KDwell(delta, prefactor) / gamma;
    }

    /// <summary>The F58 generalisation of the K_dwell prefactor for even-weight
    /// states: <c>prefactor = (2 + 4·W₂) / (1 + 6·W₂)</c>, where W₂ is the
    /// light-face Pauli sector weight at the crossing moment. Bell+ has
    /// W₂ = 0.3709 → prefactor = 1.080088, matching
    /// <see cref="BellPlusKDwellPrefactor"/> exactly.</summary>
    public double EvenWeightPrefactor(double w2)
    {
        if (w2 < 0.0 || w2 > 1.0)
            throw new ArgumentOutOfRangeException(nameof(w2), w2, "W₂ must be in [0, 1].");
        return (2.0 + 4.0 * w2) / (1.0 + 6.0 * w2);
    }

    /// <summary>The F59 generalisation to any two-sector state:
    /// <c>prefactor = (4/k) · (W₀ + W_k) / (W₀ + 3·W_k)</c>. Bell+ recovers as
    /// <c>k = 2, W₀ = 1/2</c>. Verified Bell+ and W₃ per ANALYTICAL_FORMULAS.</summary>
    public double TwoSectorPrefactor(int k, double w0, double wk)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 1.");
        if (w0 < 0.0 || wk < 0.0)
            throw new ArgumentOutOfRangeException("Sector weights must be ≥ 0.");
        double denom = w0 + 3.0 * wk;
        if (Math.Abs(denom) < 1e-15)
            throw new ArgumentException("Two-sector prefactor undefined: W₀ + 3·W_k ≈ 0.");
        return (4.0 / k) * (w0 + wk) / denom;
    }

    /// <summary>Cross-check: the live <see cref="CrossingThreshold"/> from the
    /// dyadic ladder equals the QuarterAsBilinearMaxval anchor's pinned value
    /// (1/4) bit-exactly. Drift indicator.</summary>
    public bool ThresholdMatchesQuarterAnchor() =>
        Math.Abs(CrossingThreshold - 0.25) < 1e-15;

    public F57DwellTimeQuarterPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QuarterAsBilinearMaxvalClaim quarter,
        F25CPsiBellPlusPi2Inheritance f25,
        ArgmaxMaxvalPairClaim argmaxMaxval)
        : base("F57 t_dwell = 2δ/|dCΨ/dt| at CΨ=1/4 inherits from Pi2-Foundation: 1/4 = a_3 (Quarter), 2 = a_0 (root d), mother claim F25",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F57 + " +
               "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md (Section 6) + " +
               "experiments/CPSI_COMPLEX_PLANE.md + " +
               "data/ibm_cusp_slowing_april2026/ (hardware run, ibm_kingston Heron r2) + " +
               "compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxvalClaim + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs (mother claim, typed) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (ArgmaxMaxvalPair)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
        ArgmaxMaxval = argmaxMaxval ?? throw new ArgumentNullException(nameof(argmaxMaxval));
    }

    public override string DisplayName =>
        "F57 trajectory dwell time at CΨ = 1/4 as Pi2-Foundation inheritance";

    public override string Summary =>
        $"t_dwell = 2δ / |dCΨ/dt|_{{t_cross}}; K_dwell = γ·t_dwell γ-invariant; 1/4 boundary = a_3 (Quarter), " +
        $"2δ-window factor = a_0 = root d; Bell+ prefactor 1.080088 state-specific (F58 even-weight: (2+4W₂)/(1+6W₂)) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F57 closed form",
                summary: "t_dwell(δ) = 2δ/|dCΨ/dt|_{t_cross}; K_dwell = γ·t_dwell γ-independent (std < 2×10⁻¹⁷ across γ ∈ [0.1, 10.0])");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "two anchors: CrossingThreshold = 1/4 = a_3 (QuarterAsBilinearMaxval); WindowDoublingFactor = 2 = a_0 (polynomial root d)");
            yield return InspectableNode.RealScalar("CrossingThreshold (= a_3 = 1/4)", CrossingThreshold);
            yield return InspectableNode.RealScalar("WindowDoublingFactor (= a_0 = 2)", WindowDoublingFactor);
            yield return InspectableNode.RealScalar("BellPlusKDwellPrefactor (state-specific, NOT Pi2)", BellPlusKDwellPrefactor);
            yield return new InspectableNode("Hardware verification",
                summary: "ibm_kingston Heron r2 2026-04-16: K_dwell/δ = 0.649 (pair A) / 0.694 (pair B), spread 6.3% despite 2.55× γ difference; γ-invariance confirmed (per data/ibm_cusp_slowing_april2026/)");
            yield return new InspectableNode("F58/F59 sibling readings",
                summary: "F58 even-weight prefactor (2+4·W₂)/(1+6·W₂); F59 two-sector (4/k)·(W₀+W_k)/(W₀+3·W_k); Bell+ at W₂=0.3709 / k=2,W₀=1/2 reproduces 1.080088");
            // Sample state-specific prefactors via F58
            yield return new InspectableNode(
                "F58 reading at W₂=0.3709 (Bell+)",
                summary: $"prefactor = (2 + 4·0.3709)/(1 + 6·0.3709) = {EvenWeightPrefactor(0.3709):G6} (matches Bell+ K_dwell prefactor 1.080088)");
        }
    }
}
