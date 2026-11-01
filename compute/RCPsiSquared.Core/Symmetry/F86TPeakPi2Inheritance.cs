using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F86a's <c>t_peak = 1/(4γ₀)</c> has its "4" denominator sitting
/// on the Pi2-Foundation: its value is <c>a_{−1}</c> on the dyadic halving ladder.
///
/// <para>The "4" is not a free parameter: the 4γ₀ is the rung-2 rate 2γ₀·2 of the
/// Absorption Theorem (the slowest k = 1 rate at the EP, λ = −4γ₀), a rate four and not
/// the operator-space count d² (docs/proofs/PROOF_ABSORPTION_THEOREM.md, the remark on the
/// recurring four).</para>
///
/// <para>Mirror reading: 1/4 = <c>a_3</c> on the lower (memory) side of the ladder
/// (= <see cref="QuarterAsBilinearMaxvalClaim"/>), so <c>t_peak = a_3 · (1/γ₀)</c>
/// is the alternate composition. The two agree as a value identity on the ladder,
/// a_{−1}·a_3 = 1; the Quarter's ¼ is a different object from the rate four.</para>
///
/// <para>Tier1Derived: pure composition. F86a is Tier1Derived in
/// <c>docs/proofs/PROOF_F86A_EP_MECHANISM.md</c> (Statement 1), the theorem the
/// 2026-05-14 split carved out of the monolith whose name
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> now carries as the hub;
/// <see cref="TPeakLaw"/> wraps it as a typed
/// claim parameterised by γ₀.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86A_EP_MECHANISM.md</c> Statement 1 +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F86 +
/// <c>compute/RCPsiSquared.Core/F86/TPeakLaw.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F86TPeakPi2Inheritance : Claim, IF99AnchorBearing, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    /// <inheritdoc />
    /// <remarks>Parent role: feeds the F99 inheritance graph structurally
    /// but the claim's own value does not sit on the F86b α-axis.
    /// 1/4 mirror reading on inversion identity a_{-1}·a_3 = 4·(1/4) = 1 via
    /// Quarter, not direct F86b α.</remarks>
    public F99AnchorRole F99Role => F99AnchorRole.Parent;

    /// <inheritdoc />
    public IReadOnlyList<double> F99AnchorValues { get; } = Array.Empty<double>();

    public Pi2DyadicLadderClaim Ladder { get; }
    public TPeakLaw TPeak { get; }
    /// <summary>The "4" denominator in F86's <c>t_peak = 1/(4γ₀)</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> = 4, the rung-2
    /// rate 2·2.</summary>
    public double FourFactor => Ladder.Term(-1);

    /// <summary>The mirror partner: <c>1/4 = a_3</c> on the memory side. Equivalent
    /// reading via the inversion symmetry <c>a_n · a_{2−n} = 1</c>.</summary>
    public double OneOverFourFactor => Ladder.Term(3);

    /// <summary>The γ₀ used in the parent <see cref="TPeakLaw"/> instance.</summary>
    public double GammaZero => TPeak.GammaZero;

    /// <summary>The composed t_peak value: <c>1 / (FourFactor · GammaZero)</c>;
    /// bit-exactly equal to <see cref="TPeakLaw.Value"/>.</summary>
    public double LiveTPeak => 1.0 / (FourFactor * GammaZero);

    /// <summary>1/4 = (1/2)² bilinear-apex maxval — the typed mirror-partner
    /// parent. F86's <c>t_peak = 1/(4γ₀) = a_3 / γ₀ = Quarter / γ₀</c> on the
    /// memory side of the inversion identity <c>a_{-1} · a_3 = 4 · (1/4) = 1</c>.
    /// The class's prose docstring already says "1/4 = a_3 ... QuarterAsBilinearMaxvalClaim";
    /// added 2026-05-16 (Wave 6b) as a typed ctor parent to make the mirror-partner
    /// edge participate in the graph.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public F86TPeakPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        TPeakLaw tPeak,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("F86 t_peak's 4 denominator inherits from Pi2-Foundation (4 = a_{-1}, the rung-2 rate 2·2); mirror partner = 1/4 = a_3 = Quarter",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86A_EP_MECHANISM.md (Statement 1) + " +
               "docs/ANALYTICAL_FORMULAS.md F86 + " +
               "compute/RCPsiSquared.Core/F86/TPeakLaw.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (QuarterAsBilinearMaxval, typed mirror partner)")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        TPeak = tPeak ?? throw new ArgumentNullException(nameof(tPeak));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        "F86 t_peak's 4 denominator as Pi2-Foundation a_{-1}";

    public override string Summary =>
        $"t_peak = 1/(4γ₀): 4 = a_{{-1}} = {FourFactor} (the rung-2 rate 2·2); mirror 1/4 = a_3 = {OneOverFourFactor} (memory side); live t_peak = {LiveTPeak:G6} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F86a",
                summary: "t_peak = 1/(4γ₀) (Tier1Derived in PROOF_F86A_EP_MECHANISM Statement 1, hub PROOF_F86_QPEAK; universal across c, N, n, bond)");
            yield return InspectableNode.RealScalar("FourFactor (= a_{-1}, the rung-2 rate 2·2)", FourFactor);
            yield return InspectableNode.RealScalar("OneOverFourFactor (= a_3 mirror partner)", OneOverFourFactor);
            yield return InspectableNode.RealScalar("γ₀", GammaZero);
            yield return InspectableNode.RealScalar("LiveTPeak", LiveTPeak);
            yield return InspectableNode.RealScalar("TPeakLaw.Value (parent's pinned)", TPeak.Value);
            yield return new InspectableNode("inheritance",
                summary: "the 4 in the denominator is not a free parameter; it is the rung-2 rate 2γ₀·2 (the slowest k = 1 rate at the EP), valued a_{-1} on the upper side of the Pi2 ladder");
            yield return new InspectableNode("mirror reading",
                summary: "via inversion a_n · a_{2-n} = 1: t_peak = a_3 / γ₀ = (1/4)/γ₀; a value identity on the ladder, the Quarter's ¼ being a different object from the rate four");
        }
    }
}
