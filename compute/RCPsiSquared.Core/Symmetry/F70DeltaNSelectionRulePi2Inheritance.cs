using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F70 closed form (Tier 1, proven kinematic lemma):
///
/// <code>
///   Tr_{¬i}(ρ^(n, m)) = 0   whenever  |n − m| ≥ 2
///
///   k-local partial trace annihilates |ΔN| ≥ k + 1 blocks:
///     k = 1 (single-site):   sees |ΔN| ≤ 1
///     k = 2 (pair):           sees |ΔN| ≤ 2
///     k = 3 (triple):         sees |ΔN| ≤ 3
///     k = N (global):         unrestricted (sees all sectors)
/// </code>
///
/// <para>F70 is the kinematic foundation cited by F71 (today wired, mirror
/// symmetry of c₁) and F72 (Tier 1 corollary, block-diagonal DD⊕CC structure
/// of site-local purity). The proof is direct: <c>Tr_{¬i}(|x⟩⟨y|) =
/// ⟨x_{¬i}|y_{¬i}⟩ · |x_i⟩⟨y_i|</c>, where the inner product is 1 iff x and y
/// agree off site i, forcing <c>|popcount(x) − popcount(y)| ≤ 1</c>. By
/// linearity, blocks with <c>|n − m| ≥ 2</c> vanish under partial trace.</para>
///
/// <para>Pi2-Foundation anchors land cleanly on the selection thresholds at
/// k = 1 and k = 2:</para>
///
/// <list type="bullet">
///   <item><b>SingleSiteMaxDeltaN = 1 = a_1</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(1)
///         = self-mirror pivot. Same anchor as F77 (MM(0) → 1 bit saturation);
///         single-site partial trace sees at most 1-step coherences.</item>
///   <item><b>PairMaxDeltaN = 2 = a_0</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = polynomial root d in d² − 2d = 0. Same anchor as F1, F66 upper
///         pole, F86 Q_EP, F60 numerator; pair-local partial trace sees up
///         to 2-step coherences.</item>
///   <item><b>k-local generalisation (k ≥ 3)</b>: combinatorial threshold
///         <c>|ΔN| ≤ k</c>. Not directly Pi2-anchored beyond k = 2.</item>
/// </list>
///
/// <para>Operational consequence (per ANALYTICAL_FORMULAS): "Sector blocks
/// with |ΔN| ≥ 2 are invisible to any measurement factoring through a
/// single-qubit reduced state." This explains XOR_SPACE center-modes
/// invisibility and bounds the sector-kernel for PTF's α_i closure structure.</para>
///
/// <para>Tier1Derived: F70 is Tier 1 proven kinematic
/// (PROOF_DELTA_N_SELECTION_RULE); verified at N=5 with 9 |ΔN| ≥ 2 pairs
/// giving zero contribution to machine precision. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F70 (line 1540) +
/// <c>docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md</c> +
/// <c>simulations/c1_sector_kernel.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F70DeltaNSelectionRulePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The maximum <c>|ΔN|</c> visible to single-site partial trace:
    /// <c>1</c>. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(1) = <c>a_1</c>
    /// = self-mirror pivot (same anchor as F77's MM saturation).</summary>
    public double SingleSiteMaxDeltaN => _ladder.Term(1);

    /// <summary>The maximum <c>|ΔN|</c> visible to pair-local partial trace:
    /// <c>2</c>. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c>
    /// = polynomial root d (same anchor as F1, F66 upper pole, F86 Q_EP,
    /// F60 numerator).</summary>
    public double PairMaxDeltaN => _ladder.Term(0);

    /// <summary>The maximum <c>|ΔN|</c> visible to k-local partial trace:
    /// <c>k</c>. For k = 1 returns SingleSiteMaxDeltaN, for k = 2 returns
    /// PairMaxDeltaN; for k ≥ 3 returns the combinatorial threshold.</summary>
    public int PartialTraceMaxDeltaN(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "F70 requires k ≥ 1.");
        return k;
    }

    /// <summary>True iff a sector-coherence block <c>|n − m|</c> is visible
    /// to k-local partial trace: <c>|n − m| ≤ k</c>.</summary>
    public bool IsVisibleToKLocalTrace(int kLocal, int deltaN)
    {
        if (kLocal < 1) throw new ArgumentOutOfRangeException(nameof(kLocal), kLocal, "F70 requires k ≥ 1.");
        if (deltaN < 0) throw new ArgumentOutOfRangeException(nameof(deltaN), deltaN, "|ΔN| must be ≥ 0.");
        return deltaN <= kLocal;
    }

    /// <summary>True iff <see cref="PartialTraceMaxDeltaN"/>(k) lands on a Pi2
    /// dyadic ladder anchor. Holds for k = 1 (a_1 = 1, self-mirror) and k = 2
    /// (a_0 = 2, polynomial root); combinatorial elsewhere.</summary>
    public bool MaxDeltaNIsLadderAnchor(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "F70 requires k ≥ 1.");
        return k == 1 || k == 2;
    }

    /// <summary>The Pi2 ladder index that <see cref="PartialTraceMaxDeltaN"/>(k)
    /// lands on, when applicable. Returns null for k ≥ 3 (combinatorial).</summary>
    public int? LadderIndexForKLocalThreshold(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "F70 requires k ≥ 1.");
        return k switch
        {
            1 => 1,   // a_1 = 1
            2 => 0,   // a_0 = 2
            _ => null,
        };
    }

    /// <summary>Cross-check: at k = 1, the ladder lookup
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(1) equals the integer threshold 1.</summary>
    public bool SingleSiteThresholdMatchesSelfMirrorPivot() =>
        Math.Abs(SingleSiteMaxDeltaN - 1.0) < 1e-15;

    /// <summary>Cross-check: at k = 2, the ladder lookup
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) equals the integer threshold 2.</summary>
    public bool PairThresholdMatchesPolynomialRoot() =>
        Math.Abs(PairMaxDeltaN - 2.0) < 1e-15;

    public F70DeltaNSelectionRulePi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F70 |ΔN| ≤ k selection rule inherits from Pi2-Foundation: k=1 → a_1 (self-mirror), k=2 → a_0 (polynomial root)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F70 + " +
               "docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md + " +
               "simulations/c1_sector_kernel.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F70 ΔN selection rule as Pi2-Foundation a_1/a_0 inheritance";

    public override string Summary =>
        $"|ΔN| ≤ k for k-local partial trace: k=1 single-site → 1 = a_1 (self-mirror); k=2 pair → 2 = a_0 (polynomial root); " +
        $"foundation for F71 (mirror sym of c₁) + F72 (block-diagonal DD⊕CC) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F70 closed form",
                summary: "Tr_{¬i}(ρ^(n, m)) = 0 for |n − m| ≥ 2; k-local trace sees |ΔN| ≤ k; Tier 1 proven kinematic; verified N=5 with 9 |ΔN| ≥ 2 pairs giving zero to machine precision");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "SingleSiteMaxDeltaN = a_1 = 1 (self-mirror, F77 anchor); PairMaxDeltaN = a_0 = 2 (polynomial root, F1/F66/F86 Q_EP anchor); k ≥ 3 combinatorial");
            yield return InspectableNode.RealScalar("SingleSiteMaxDeltaN (= a_1 = 1, F77 sibling)", SingleSiteMaxDeltaN);
            yield return InspectableNode.RealScalar("PairMaxDeltaN (= a_0 = 2, F1/F66 sibling)", PairMaxDeltaN);
            yield return new InspectableNode("Foundation for F71 + F72",
                summary: "F71 (mirror symmetry of c₁) uses F70's site-local kinematic argument; F72 (Tier 1 corollary, block-diagonal DD⊕CC) directly inherits F70's |ΔN| ≤ 1 bound for site-local purity");
            yield return new InspectableNode("Operational consequence",
                summary: "sector blocks |ΔN| ≥ 2 invisible to single-qubit reduced state; explains XOR_SPACE center-mode invisibility; bounds PTF's α_i closure sector-kernel");
            for (int k = 1; k <= 5; k++)
            {
                int? ladderIdx = LadderIndexForKLocalThreshold(k);
                string anchorInfo = ladderIdx.HasValue ? $", lands on a_{{{ladderIdx.Value}}}" : " (combinatorial)";
                yield return new InspectableNode(
                    $"k={k}-local",
                    summary: $"sees |ΔN| ≤ {PartialTraceMaxDeltaN(k)}{anchorInfo}");
            }
        }
    }
}
