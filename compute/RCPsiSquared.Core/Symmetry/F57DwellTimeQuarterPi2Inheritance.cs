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
///     K_dwell = 1.08008787 · δ
///     prefactor = 2 / |dCΨ/dt|_{t_cross}/γ  (solved live off F25)
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
///   <item><b>State-specific prefactor (NOT Pi2-anchored)</b>: Bell+'s comes
///         from <c>|dCΨ/dt|_{t_cross}</c> at the F25 closed-form crossing.
///         F58 generalises to even-weight states via
///         <c>(2 + 4·W₂)/(1 + 6·W₂)</c>; F59 generalises further to any
///         two-sector state. Bell+ is the W₀ = 1/2, k = 2 special case, and
///         its W₂ = f*²/2 makes F58 and F57 the same identity rather than two
///         numbers that agree to four digits.</item>
/// </list>
///
/// <para>Tier1Derived: F57 is Tier 1 analytical (CRITICAL_SLOWING_AT_THE_CUSP §6). The
/// γ-invariance is exact for the ideal Bell+/pure-Z trajectory the formula models.
/// The ibm_kingston Heron r2 comparison is approximate and two-pair: K_dwell/δ = 0.6492
/// (pair A) and 0.6937 (pair B), agreeing to 6.4% despite a 2.55× γ difference, at an
/// absolute prefactor 0.67 rather than 1.0801. Pair identity and T1 vary with γ (Kingston
/// has T1 ≈ T2; F57 assumes pure Z-dephasing), so the comparison neither establishes
/// γ-invariance beyond that model nor isolates the cause of the prefactor gap.
/// The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F57 +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> Section 6 +
/// <c>experiments/CPSI_COMPLEX_PLANE.md</c> +
/// <c>data/ibm_cusp_slowing_april2026/</c> (hardware run) +
/// <c>compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxvalClaim</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F57DwellTimeQuarterPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    public Pi2DyadicLadderClaim Ladder { get; }
    public QuarterAsBilinearMaxvalClaim Quarter { get; }
    /// <summary>F25 Bell+ CΨ closed form — the typed mother claim. F57's
    /// Bell+ prefactor is <c>2 / |dCΨ/dt|_{t_cross}/γ</c>, read straight off
    /// F25 rather than re-tabulated here. Added 2026-05-16 as a typed ctor
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
    public double CrossingThreshold => Ladder.Term(3);

    /// <summary>The 2δ-window doubling factor: <c>2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = polynomial
    /// root d (the "two" in d²−2d=0). Reflects the symmetric δ-above /
    /// δ-below structure of the bilinear apex.</summary>
    public double WindowDoublingFactor => Ladder.Term(0);

    /// <summary>The Bell+ Z-dephasing K_dwell prefactor, live off the mother
    /// claim: <c>2 / |dCΨ/dt|_{t_cross} per γ = 1.0800878671056402</c>.
    /// State-specific, NOT Pi2-anchored. Documents quote it rounded to
    /// 1.080088.</summary>
    public double BellPlusKDwellPrefactor => F25.BellPlusF57Prefactor;

    /// <summary>Bell+'s light-face weight at the crossing, <c>W₂ = f*²/2</c>,
    /// the argument F58 takes. It is not a fitted number: under Z-dephasing
    /// Bell⁺ is <c>ρ = ¼(II + ZZ + f·XX − f·YY)</c>, so with the sector weight
    /// read as <c>Σ_P c_P²/4</c> the frozen diagonal carries <c>W₀ = (1+1)/4
    /// = 1/2</c> at every time and the light face carries
    /// <c>W₂ = (f² + f²)/4 = f²/2</c>. At the crossing f = f*, giving
    /// <c>0.3708534749861196</c>; documents quote it rounded to 0.3709.</summary>
    public double BellPlusW2AtCrossing
    {
        get { double f = F25.BellPlusFCross; return f * f / 2.0; }
    }

    /// <summary>Bell+'s frozen-diagonal weight <c>W₀ = 1/2</c>, the other
    /// argument F59 takes. Z-dephasing leaves the II and ZZ coefficients at 1,
    /// so this is a constant of the trajectory, not a reading at the
    /// crossing.</summary>
    public double BellPlusW0 => 0.5;

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
    /// <see cref="BellPlusW2AtCrossing"/> = f*²/2 → exactly
    /// <see cref="BellPlusKDwellPrefactor"/>, and not approximately: substituting
    /// W₂ = f²/2 turns <c>(2+4W₂)/(1+6W₂)</c> into <c>(2+2f²)/(1+3f²)</c>, which
    /// equals F57's own <c>3/(f(1+3f²))</c> precisely when <c>f(1+f²) = 3/2</c>,
    /// the F25 crossing equation. F58 and F57 are one identity at the crossing,
    /// so the agreement gates that equation rather than a tabulated weight
    /// (see <see cref="EvenWeightPrefactorReducesToF57"/>).</summary>
    public double EvenWeightPrefactor(double w2)
    {
        if (w2 < 0.0 || w2 > 1.0)
            throw new ArgumentOutOfRangeException(nameof(w2), w2, "W₂ must be in [0, 1].");
        return (2.0 + 4.0 * w2) / (1.0 + 6.0 * w2);
    }

    /// <summary>The F59 generalisation to any two-sector state:
    /// <c>prefactor = (4/k) · (W₀ + W_k) / (W₀ + 3·W_k)</c>. Bell+ recovers as
    /// <c>k = 2</c>, <see cref="BellPlusW0"/>, <see cref="BellPlusW2AtCrossing"/>,
    /// where it collapses term for term onto F58's even-weight form. Verified
    /// Bell+ and W₃ per ANALYTICAL_FORMULAS.</summary>
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

    /// <summary>The identity above, computed: F58 evaluated at the exact
    /// W₂ = f*²/2 against <see cref="BellPlusKDwellPrefactor"/>. Both routes are
    /// float, so the residual is rounding only: measured 4.441·10⁻¹⁶, two eps,
    /// and gated at 8 eps. Feeding the rounded 0.3709 instead misses by
    /// 3.578·10⁻⁵, a factor 8.06·10¹⁰ larger, which is what the round-trip through
    /// a four-digit weight costs.</summary>
    public bool EvenWeightPrefactorReducesToF57() =>
        Math.Abs(EvenWeightPrefactor(BellPlusW2AtCrossing) - BellPlusKDwellPrefactor)
            < 8.0 * F25CPsiBellPlusPi2Inheritance.MachineEps;

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
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
        ArgmaxMaxval = argmaxMaxval ?? throw new ArgumentNullException(nameof(argmaxMaxval));
    }

    public override string DisplayName =>
        "F57 trajectory dwell time at CΨ = 1/4 as Pi2-Foundation inheritance";

    public override string Summary =>
        $"t_dwell = 2δ / |dCΨ/dt|_{{t_cross}}; K_dwell = γ·t_dwell γ-invariant; 1/4 boundary = a_3 (Quarter), " +
        $"2δ-window factor = a_0 = root d; Bell+ prefactor {BellPlusKDwellPrefactor:G8} state-specific, and F58's (2+4W₂)/(1+6W₂) at W₂ = f*²/2 IS that prefactor by the F25 crossing equation ({Tier.Label()})";

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
            yield return new InspectableNode("Hardware comparison",
                summary: "ibm_kingston Heron r2 2026-04-16: K_dwell/δ = 0.6492 (pair A) / 0.6937 (pair B), agreeing to 6.4% despite 2.55× γ difference, at prefactor 0.67 rather than 1.0801; an approximate two-pair check, not a verification of γ-invariance (per data/ibm_cusp_slowing_april2026/)");
            yield return InspectableNode.RealScalar("BellPlusW2AtCrossing (= f*²/2)", BellPlusW2AtCrossing);
            yield return new InspectableNode("F58/F59 sibling readings",
                summary: "F58 even-weight prefactor (2+4·W₂)/(1+6·W₂); F59 two-sector (4/k)·(W₀+W_k)/(W₀+3·W_k); Bell+ enters both with W₀ = 1/2 (frozen diagonal) and W₂ = f*²/2 (light face at the crossing), k = 2");
            yield return new InspectableNode(
                "F58 at Bell+'s exact W₂",
                summary: $"prefactor = (2 + 4·W₂)/(1 + 6·W₂) at W₂ = {BellPlusW2AtCrossing:G10} gives {EvenWeightPrefactor(BellPlusW2AtCrossing):G12}, against the F25 route's {BellPlusKDwellPrefactor:G12}; identity holds (drift check: {EvenWeightPrefactorReducesToF57()})");
        }
    }
}
