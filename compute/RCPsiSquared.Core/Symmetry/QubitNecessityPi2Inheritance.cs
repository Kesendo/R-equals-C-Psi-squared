using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The per-site Pauli basis split <c>4 = 2 + 2</c> at d=2 is the direct
/// algebraic reading of the polynomial trunk <c>d² − 2d = 0</c> on the Pi2 dyadic
/// halving ladder. Each piece sits on a typed Pi2-Foundation anchor, and their
/// algebraic relationship <c>a_{-1} = 2 · a_0</c> IS the bijection equation that
/// d² − 2d = 0 enforces.
///
/// <para>From <c>docs/QUBIT_NECESSITY.md</c>: under single-axis dephasing, the per-site
/// Pauli basis splits into immune (d operators, diagonal in dephasing eigenbasis) and
/// decaying (d² − d operators). The palindromic mirror requires bijection between these
/// subsets, giving <c>d = d² − d</c>, equivalently <c>d² − 2d = 0</c>, equivalently
/// <c>d(d − 2) = 0</c>. Only d=2 satisfies the bijection. At d=2: 4 total = 2 immune
/// + 2 decaying.</para>
///
/// <list type="bullet">
///   <item><b>TotalPauliOpsPerSite</b> = 4 = <c>a_{-1}</c> (Pi2DyadicLadder, equivalently
///         d² for N=1 qubit per <see cref="Pi2OperatorSpaceMirrorClaim"/>).</item>
///   <item><b>ImmuneOpsPerSite</b> = 2 = <c>a_0</c> (Pi2DyadicLadder; the qubit dimension d).</item>
///   <item><b>DecayingOpsPerSite</b> = 2 = <c>a_0</c> (same; balanced).</item>
///   <item><b>BijectionEquation</b>: <c>a_{-1} = 2 · a_0</c>, i.e. 4 = 2 + 2 (or 4 = 2 · 2,
///         since immune and decaying counts are equal).</item>
///   <item><b>BalancedFraction</b>: ImmuneOps / TotalOps = 2/4 = 1/2 = <c>a_2</c>
///         (<see cref="HalfAsStructuralFixedPointClaim"/>; C = 0.5 universal).</item>
/// </list>
///
/// <para>Tier1Derived: the polynomial trunk d² − 2d = 0 (Tier1Derived in
/// <see cref="PolynomialFoundationClaim"/>) directly imposes the per-site basis
/// 4=2+2 split at its only non-zero solution d=2. This claim makes the typed
/// per-site reading explicit: not a new fact, the same algebraic content seen at the
/// per-site Pauli basis level.</para>
///
/// <para>QUBIT_NECESSITY documents 0/236 qutrit dissipators palindromic and the
/// d=3 split 9 = 3 + 6 (3:6 imbalance). The d=2 case is uniquely balanced because
/// <c>d² − 2d = 0</c> closes ONLY at d=2.</para>
///
/// <para>Anchors: <c>docs/QUBIT_NECESSITY.md</c> +
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (PolynomialFoundationClaim, HalfAsStructuralFixedPointClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class QubitNecessityPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>Total Pauli operators per site at d=2: 4 = <c>a_{-1}</c> on the
    /// Pi2 ladder (= d² for N=1 qubit, the operator-space-side anchor).</summary>
    public double TotalPauliOpsPerSite => _ladder.Term(-1);

    /// <summary>Immune operators per site (diagonal in dephasing eigenbasis): 2 = <c>a_0</c>
    /// on the Pi2 ladder (= d, the qubit dimension).</summary>
    public double ImmuneOpsPerSite => _ladder.Term(0);

    /// <summary>Decaying operators per site (anti-diagonal in dephasing eigenbasis):
    /// 2 = <c>a_0</c>. Equal to immune count by the d² − 2d = 0 bijection at d=2.</summary>
    public double DecayingOpsPerSite => _ladder.Term(0);

    /// <summary>The balanced fraction: ImmuneOps / TotalOps = 2/4 = 1/2 = <c>a_2</c>.
    /// This is exactly <see cref="HalfAsStructuralFixedPointClaim"/>'s C = 0.5
    /// universal anchor manifested at the per-site basis level.</summary>
    public double BalancedFraction => ImmuneOpsPerSite / TotalPauliOpsPerSite;

    /// <summary>Cross-verification through the operator-space mirror: TotalPauliOpsPerSite
    /// must equal d² = 4^1 = 4 for N=1 qubit, pinned in
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/>.</summary>
    public double MirrorPinnedTotalOps =>
        _mirror.PairAt(1)?.OperatorSpace
        ?? throw new InvalidOperationException("Pi2OperatorSpaceMirrorClaim must contain N=1 anchor.");

    /// <summary>Live drift check on the bijection equation: the d² − 2d = 0 polynomial
    /// at d=2 says total Pauli = 2 · immune (or equivalently total = immune + decaying with
    /// immune = decaying). Returns true when <c>TotalPauliOpsPerSite == 2 · ImmuneOpsPerSite</c>
    /// and <c>ImmuneOpsPerSite == DecayingOpsPerSite</c>.</summary>
    public bool BijectionHolds =>
        Math.Abs(TotalPauliOpsPerSite - 2.0 * ImmuneOpsPerSite) < 1e-12
        && Math.Abs(ImmuneOpsPerSite - DecayingOpsPerSite) < 1e-12;

    public QubitNecessityPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("Per-site Pauli basis split 4 = 2 + 2 (the typed reading of d² − 2d = 0 at d=2)",
               Tier.Tier1Derived,
               "docs/QUBIT_NECESSITY.md + " +
               "docs/proofs/MIRROR_SYMMETRY_PROOF.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "QubitNecessity per-site basis 4 = 2 + 2 inherits from Pi2-Foundation";

    public override string Summary =>
        $"d² − 2d = 0 at d=2: TotalPauli = a_{{-1}} = {TotalPauliOpsPerSite}; " +
        $"Immune = Decaying = a_0 = {ImmuneOpsPerSite}; balanced fraction = a_2 = {BalancedFraction}; " +
        $"bijection a_{{-1}} = 2·a_0 holds: {BijectionHolds} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("polynomial trunk",
                summary: "d² − 2d = 0 ↔ d(d−2) = 0 ↔ d=0 or d=2 (PolynomialFoundationClaim)");
            yield return new InspectableNode("per-site Pauli basis at d=2",
                summary: "4 total operators = 2 immune (diagonal: I, Z) + 2 decaying (anti-diagonal: X, Y) for Z-dephasing");
            yield return InspectableNode.RealScalar("TotalPauliOpsPerSite (= a_{-1} = d²)", TotalPauliOpsPerSite);
            yield return InspectableNode.RealScalar("ImmuneOpsPerSite (= a_0 = d)", ImmuneOpsPerSite);
            yield return InspectableNode.RealScalar("DecayingOpsPerSite (= a_0 = d)", DecayingOpsPerSite);
            yield return InspectableNode.RealScalar("BalancedFraction (= a_2 = 1/2)", BalancedFraction);
            yield return new InspectableNode("bijection equation",
                summary: $"a_{{-1}} = 2 · a_0 ↔ 4 = 2 · 2 ↔ d² = 2 · d at d=2; holds: {BijectionHolds}");
            yield return InspectableNode.RealScalar("MirrorPinnedTotalOps (cross-check)", MirrorPinnedTotalOps);
            yield return new InspectableNode("d > 2 fails",
                summary: "d=3: 9 = 3 + 6 (3:6 imbalance, no bijection); 0/236 qutrit dissipators palindromic (QUBIT_NECESSITY Section 9)");
            yield return new InspectableNode("HalfAsStructuralFixedPoint reading",
                summary: "the 1/2 balanced fraction IS C = 0.5 universal, here at the per-site Pauli basis level");
        }
    }
}
