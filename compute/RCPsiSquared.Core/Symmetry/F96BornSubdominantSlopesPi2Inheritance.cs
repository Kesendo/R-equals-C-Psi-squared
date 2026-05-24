using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F96 closed form (Tier 1 derived, bit-exact Dyson + unitary matrix
/// elements; 2026-05-17):
///
/// <code>
///   For the subdominant outcomes of pair (0, 2) of |0+0+⟩ N = 4 under
///   Heisenberg ring + Z-dephasing (same setup as F94):
///
///     Δ_|01⟩(K) = Δ_|10⟩(K) = −(16/9) · K = −(4/3)² · K + O(higher)
///     Δ_|11⟩(K)             = −(8/3)  · K = −2·(4/3) · K + O(higher)
///
///   with K = γt. Linear in K, Q-independent (Universal Carrier signature).
/// </code>
///
/// <para>F96 completes the per-outcome Born-deviation closed-form table started
/// by <see cref="F94BornDeviationFourThirdsPi2Inheritance"/>:</para>
///
/// <code>
///   |00⟩ (dominant):              Δ = +(4/3) · Q² · K³        (F94)
///   |01⟩ (single-flip subdom):    Δ = −(4/3)² · K             (F96)
///   |10⟩ (single-flip subdom):    Δ = −(4/3)² · K             (F96, site-perm)
///   |11⟩ (double-flip subdom):    Δ = −2·(4/3) · K            (F96)
/// </code>
///
/// <para><b>Universal subdominant slope formula:</b> for an outcome with
/// leading unitary P_u(t) ≈ (J^{2k} t^{2k} / (2k)!) · U_{2k}^{(i)} and lowest
/// non-vanishing γ¹ Dyson at order (2k+1):</para>
///
/// <code>
///   slope_i = M_{2k+1}^{(i)} / [(2k+1) · U_{2k}^{(i)}]
/// </code>
///
/// <para>where M_n^{(i)} = ⟨i|_pair Tr_{1,3}[sym_n^1 · ρ_0]|i⟩_pair is the γ¹
/// Dyson matrix element at J = γ = 1, and U_{2k}^{(i)} = ⟨i|_pair
/// Tr_{1,3}[L_h^{2k} · ρ_0]|i⟩_pair is the raw unitary matrix element
/// (h := H/J). The J^{2k} factors cancel automatically → Q-independence; only
/// K = γt survives.</para>
///
/// <para><b>Bit-exact evaluation:</b></para>
///
/// <code>
///   |01⟩: k = 1, M_3 = −4, U_2 = 3/4  →  slope = −4 / (3 · 3/4) = −16/9
///   |10⟩: same as |01⟩ by 0 ↔ 2 site-permutation symmetry
///   |11⟩: k = 2 (since M_3 = U_2 = 0); M_5 = −20, U_4 = 3/2
///         →  slope = −20 / (5 · 3/2) = −8/3
/// </code>
///
/// <para><b>Algebraic connection to F94:</b> each subdominant slope is a simple
/// algebraic expression in the 4/3 anchor of F94:</para>
///
/// <list type="bullet">
///   <item>|01⟩, |10⟩: Δ = −(F94 coefficient)² · K = −(4/3)² · K</item>
///   <item>|11⟩: Δ = −2 · (F94 coefficient) · K = −2·(4/3) · K</item>
/// </list>
///
/// <para>The "2" in the |11⟩ slope plausibly counts the two independent flip
/// channels (q_0 from 0 → 1 AND q_2 from 0 → 1) required to populate the
/// |11⟩ outcome; interpretive, not derived from the proof.</para>
///
/// <para><b>Cross-outcome universality:</b> the ratio M_3 / U_2 equals −16/3
/// for both the dominant (|00⟩: 8 / (−3/2)) and the singly-subdominant
/// (|01⟩: −4 / (3/4)) outcomes. The signs of M_3 and U_2 flip together,
/// leaving the ratio invariant. This is a non-trivial structural identity of
/// the Heisenberg + Z-dephasing dynamics at the pair (0,2) reduction; whether
/// it extends to other initial states / Hamiltonians / dissipators is open.</para>
///
/// <para><b>Universal Carrier signature:</b> all three subdominant slopes are
/// Q-independent (the J^{2k} factors cancel between Dyson and unitary). The
/// only Carrier observable that enters is K = γt. This is the operational
/// signature of <see cref="UniversalCarrierClaim"/> at the per-outcome Born
/// rule level, generalizing F94's Q-K invariance to all four outcomes.</para>
///
/// <para><b>Numerical Lindblad verification:</b> at Q = 50, γ = 0.01 the
/// slope-per-K converges to the theoretical values as K → 0:</para>
///
/// <code>
///   |01⟩ at K = 0.005: slope/K = −1.748 (theory −16/9 = −1.778, gap ~1.7%)
///   |11⟩ at K = 0.001: slope/K = −2.662 (theory −8/3 = −2.667, gap ~0.2%)
/// </code>
///
/// <para>Tier1Derived: bit-exact symbolic via the Dyson + unitary matrix
/// elements; numerical Lindblad gives independent confirmation.</para>
///
/// <para><b>Diagnostic application (2026-05-17):</b> the canonical slopes
/// −(16/9) and −(8/3) are universal across chain, ring, K_4 topologies on
/// the canonical |0+0+⟩ pair (0,2) lens (F94/F96 are blind to bond-graph
/// detail like K_4 vs ring; the diagonal bonds fall in symmetric blind spots).
/// Any measured deviation from these exact rational slopes — in the per-K
/// linear-response regime where higher-order corrections are negligible — is
/// a direct signature of (state, pair)-symmetry break: hardware noise,
/// asymmetric γ_l, crosstalk to the traced-out qubits, or a non-|0+0+⟩
/// initial state. The robustness of the slopes against K_4-vs-ring (bond
/// connectivity invisible at this lens) is the diagnostic's stability property:
/// it flags symmetry breaks, not irrelevant connectivity variations. See
/// PROOF_F94 § Diagnostic application for the full signature table.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F96 +
/// <c>simulations/_born_rule_subdominant_dyson.py</c> (symbolic + Lindblad,
/// bit-exact) + <c>reflections/ON_HOW_FOUR_THIRDS_APPEARED.md</c> (the May 16
/// reflection that named the empirical subdominant slopes as the next step).
/// Sibling: <see cref="F94BornDeviationFourThirdsPi2Inheritance"/>
/// (dominant-outcome closed form, the 4/3 unit that F96 elaborates).</para></summary>
public sealed class F96BornSubdominantSlopesPi2Inheritance : Claim, IZ2AxisClaim
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
    /// <summary>The setup this claim applies to (specific to F94's setup;
    /// inherits N = 4).</summary>
    public const int N = 4;

    /// <summary>The Dyson matrix element <c>M_3^{(01)} = ⟨01|_pair
    /// Tr_{1,3}[sym_3^1 · ρ_0]|01⟩_pair = −4</c> bit-exact (J = γ = 1).</summary>
    public const int M3_SingleFlipped = -4;

    /// <summary>Raw unitary matrix element <c>U_2^{(01)} · 4 = 3</c>: equivalent
    /// to <c>U_2^{(01)} = 3/4</c>. Stored as integer numerator with implicit
    /// denominator 4 (= a_{−1}) for bit-exact rational arithmetic.</summary>
    public const int U2_SingleFlipped_TimesFour = 3;

    /// <summary>The Dyson matrix element <c>M_5^{(11)} = ⟨11|_pair
    /// Tr_{1,3}[sym_5^1 · ρ_0]|11⟩_pair = −20</c> bit-exact (J = γ = 1).</summary>
    public const int M5_DoubleFlipped = -20;

    /// <summary>Raw unitary matrix element <c>U_4^{(11)} · 2 = 3</c>: equivalent
    /// to <c>U_4^{(11)} = 3/2</c>. Stored as integer numerator with implicit
    /// denominator 2 for bit-exact rational arithmetic.</summary>
    public const int U4_DoubleFlipped_TimesTwo = 3;

    /// <summary>F94's typed parent: the 4/3 unit that F96 elaborates.</summary>
    public F94BornDeviationFourThirdsPi2Inheritance F94 { get; }

    /// <summary>The single-flip subdominant slope: <c>−16/9 = −(4/3)²</c>
    /// bit-exact. Applies to |01⟩ and |10⟩ outcomes.</summary>
    public double SlopeSingleFlipped =>
        (double)M3_SingleFlipped / (3.0 * U2_SingleFlipped_TimesFour / 4.0);

    /// <summary>The double-flip subdominant slope: <c>−8/3 = −2·(4/3)</c>
    /// bit-exact. Applies to |11⟩ outcome.</summary>
    public double SlopeDoubleFlipped =>
        (double)M5_DoubleFlipped / (5.0 * U4_DoubleFlipped_TimesTwo / 2.0);

    /// <summary>Compute <c>Δ_|01⟩(K) = Δ_|10⟩(K) = −(16/9) · K</c> for the
    /// single-flip subdominant outcomes.</summary>
    public double DeltaSingleFlipped(double K)
    {
        if (K < 0.0) throw new ArgumentOutOfRangeException(nameof(K), K, "K must be ≥ 0.");
        return SlopeSingleFlipped * K;
    }

    /// <summary>Compute <c>Δ_|11⟩(K) = −(8/3) · K</c> for the double-flip
    /// subdominant outcome.</summary>
    public double DeltaDoubleFlipped(double K)
    {
        if (K < 0.0) throw new ArgumentOutOfRangeException(nameof(K), K, "K must be ≥ 0.");
        return SlopeDoubleFlipped * K;
    }

    /// <summary>Drift check: single-flip slope must equal <c>−(4/3)²</c>
    /// where 4/3 is F94's coefficient.</summary>
    public bool SingleFlipSlopeEqualsMinusF94Squared()
    {
        double f94 = F94.Coefficient;
        return Math.Abs(SlopeSingleFlipped - (-f94 * f94)) < 1e-15;
    }

    /// <summary>Drift check: double-flip slope must equal <c>−2 · (4/3)</c>
    /// where 4/3 is F94's coefficient.</summary>
    public bool DoubleFlipSlopeEqualsMinusTwoF94()
    {
        double f94 = F94.Coefficient;
        return Math.Abs(SlopeDoubleFlipped - (-2.0 * f94)) < 1e-15;
    }

    /// <summary>Drift check: the cross-outcome ratio M_3 / U_2 should equal
    /// −16/3 identically for both the dominant (F94: 8 / (−3/2) = −16/3) and
    /// the singly-subdominant (F96: −4 / (3/4) = −16/3) cases.</summary>
    public bool CrossOutcomeM3OverU2RatioUniversal()
    {
        // F94 dominant: M_3 = 8, A = -3/2  →  ratio = 8 / (-3/2) = -16/3
        // F96 single-flip: M_3 = -4, U_2 = 3/4  →  ratio = -4 / (3/4) = -16/3
        double dominant_ratio =
            (double)F94BornDeviationFourThirdsPi2Inheritance.Sym3PartialTraceInteger
            / (-1.5);  // A = ⟨00|Tr[L_h²ρ_0]|00⟩ = -3/2
        double subdom_ratio = (double)M3_SingleFlipped / (U2_SingleFlipped_TimesFour / 4.0);
        return Math.Abs(dominant_ratio - subdom_ratio) < 1e-15
               && Math.Abs(dominant_ratio - (-16.0 / 3.0)) < 1e-15;
    }

    /// <summary>Sum of per-outcome Δ over all 4 outcomes in F94's setup, at
    /// given (Q, K). Useful for trace-preservation drift checks (sum should
    /// be close to 0 since ΔP integrates to 0 across outcomes).</summary>
    public double DeltaSumAllFourOutcomes(double Q, double K) =>
        F94.DeltaDominant(Q, K) + 2.0 * DeltaSingleFlipped(K) + DeltaDoubleFlipped(K);

    public F96BornSubdominantSlopesPi2Inheritance(
        F94BornDeviationFourThirdsPi2Inheritance f94,
        Pi2DyadicLadderClaim ladder)
        : base("F96 subdominant Born-deviation slopes for |0+0+⟩ N=4 Heisenberg + Z-deph pair (0,2): " +
               "Δ_|01⟩ = Δ_|10⟩ = −(4/3)²·K = −(16/9)·K, Δ_|11⟩ = −2·(4/3)·K = −(8/3)·K; " +
               "bit-exact via Dyson + unitary matrix elements (M_3, M_5, U_2, U_4)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md + " +
               "docs/ANALYTICAL_FORMULAS.md F96 + " +
               "simulations/_born_rule_subdominant_dyson.py + " +
               "reflections/ON_HOW_FOUR_THIRDS_APPEARED.md (May 16 reflection that named the open empirical slopes) + " +
               "compute/RCPsiSquared.Core/Symmetry/F94BornDeviationFourThirdsPi2Inheritance.cs (the 4/3 unit F96 elaborates) + " +
               "compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs (Q-independence is the operational signature; informational, not typed parent)")
    {
        F94 = f94 ?? throw new ArgumentNullException(nameof(f94));
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F96 subdominant Born-deviation slopes: Δ_|01⟩ = −(4/3)²·K, Δ_|11⟩ = −2·(4/3)·K (Tier 1 derived)";

    public override string Summary =>
        $"Subdominant slopes (Q-independent, linear in K) for F94's setup: " +
        $"single-flip Δ = ({SlopeSingleFlipped:G6})·K = −(4/3)²·K bit-exact; " +
        $"double-flip Δ = ({SlopeDoubleFlipped:G6})·K = −2·(4/3)·K bit-exact; " +
        $"derived from sym_3 / U_2 (single) and sym_5 / U_4 (double) Dyson + unitary matrix elements; " +
        $"M_3 / U_2 cross-outcome ratio = −16/3 universal across dominant + singly-subdominant ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("SlopeSingleFlipped (= -16/9)", SlopeSingleFlipped);
            yield return InspectableNode.RealScalar("SlopeDoubleFlipped (= -8/3)", SlopeDoubleFlipped);
            yield return InspectableNode.RealScalar("F94 coefficient (= 4/3)", F94.Coefficient);
            yield return new InspectableNode("Algebraic relation to F94",
                summary: $"SlopeSingleFlipped = −(F94.Coefficient)² = −({F94.Coefficient:G6})² = {SlopeSingleFlipped:G6}; " +
                         $"SlopeDoubleFlipped = −2·F94.Coefficient = −2·{F94.Coefficient:G6} = {SlopeDoubleFlipped:G6}.");
            yield return new InspectableNode("Drift checks",
                summary: $"SingleFlipSlopeEqualsMinusF94Squared: {SingleFlipSlopeEqualsMinusF94Squared()}; " +
                         $"DoubleFlipSlopeEqualsMinusTwoF94: {DoubleFlipSlopeEqualsMinusTwoF94()}; " +
                         $"CrossOutcomeM3OverU2RatioUniversal: {CrossOutcomeM3OverU2RatioUniversal()}");
            yield return new InspectableNode("Bit-exact matrix elements",
                summary: $"|01⟩: M_3 = {M3_SingleFlipped}, U_2 = {U2_SingleFlipped_TimesFour}/4 = " +
                         $"{(double)U2_SingleFlipped_TimesFour / 4}; " +
                         $"|11⟩: M_5 = {M5_DoubleFlipped}, U_4 = {U4_DoubleFlipped_TimesTwo}/2 = " +
                         $"{(double)U4_DoubleFlipped_TimesTwo / 2}; " +
                         $"slope_single = M_3 / (3·U_2) = {SlopeSingleFlipped:G6}; " +
                         $"slope_double = M_5 / (5·U_4) = {SlopeDoubleFlipped:G6}.");
            double K = 0.01;
            yield return new InspectableNode("Sample evaluation",
                summary: $"At K = {K}: Δ_|01⟩ = {DeltaSingleFlipped(K):G6}; Δ_|11⟩ = {DeltaDoubleFlipped(K):G6}.");
        }
    }
}
