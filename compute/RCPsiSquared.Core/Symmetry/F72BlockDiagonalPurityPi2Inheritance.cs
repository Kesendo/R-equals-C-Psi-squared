using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F72 closed form (Tier 1, corollary of F70):
///
/// <code>
///   Tr(ρ_i²) = 1/2 + P_i^DD[ρ₀^(diag)] + P_i^CC[ρ₀^(coh)]
///
///   no cross term coupling DD (ΔN = 0) and CC (|ΔN| = 1)
/// </code>
///
/// <para>F72 is the direct Tier 1 corollary of F70 (kinematic ΔN selection
/// rule, today wired): per-site purity decomposes strictly into a
/// diagonal-block kernel (DD) plus a coherence-block kernel (CC) with no
/// cross term. Per the F72 proof: "By F70 applied to each Bloch component,
/// ⟨Z_i⟩ depends linearly on the diagonal elements (ΔN = 0 blocks); ⟨X_i⟩
/// and ⟨Y_i⟩ depend linearly on off-diagonal elements (|ΔN| = 1 blocks).
/// Squaring keeps each contribution in its own sector class; no cross term
/// arises."</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>BaselineTraceSquared = 1/2 = a_2</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(2)
///         = <see cref="HalfAsStructuralFixedPointClaim"/>'s structural
///         fixed point. The "1/2" in <c>Tr(ρ_i²) = 1/2 + ...</c> is the
///         maximally-mixed-state minimum purity (= 1/d for d=2),
///         literally the framework's three-faces 1/2 reading.</item>
///   <item><b>TwoBlockDecompositionCount = 2 (DD + CC) = a_0</b>:
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0) = polynomial root d.
///         The k=1 partial-trace sees exactly 2 sub-blocks (DD ⊕ CC);
///         pair-observables (k=2) see 3 (DD ⊕ DC ⊕ CC); k-local sees
///         k+1 sub-blocks per F72 generalisation.</item>
///   <item><b>F70 parent (cited as direct foundation)</b>: F72 IS F70's
///         |ΔN| ≤ 1 single-site bound applied to per-site purity. F72 is
///         a typed corollary, not an independent claim.</item>
/// </list>
///
/// <para>Operational consequence (per ANALYTICAL_FORMULAS): "Any
/// closure-breaking coefficient c₁ built from per-site purities decomposes,
/// at the pre-α-fit bilinear level, into a DD-kernel and a CC-kernel with
/// no mixing. Finding the closed form of c₁ reduces to finding K_DD and
/// K_CC separately." Verified at N=5 with w-scan ρ₀(w) = cos(w)|vac⟩ +
/// sin(w)|S₁⟩; block-diagonal coupling at machine precision; pure-coherence
/// probe gives K_CC/2 to 10⁻¹².</para>
///
/// <para>k-local generalisation: at k = 1 (single-site), 2 sub-blocks
/// (DD + CC); at k = 2 (pair), 3 sub-blocks (DD + DC + CC, with DC the
/// new diagonal-coherence cross specific to pair observables); at k-local,
/// k+1 sub-blocks total (= F70.PartialTraceMaxDeltaN(k) + 1).</para>
///
/// <para>Tier1Derived: F72 is Tier 1 corollary of F70 (PROOF_DELTA_N_SELECTION_RULE
/// → Tr(ρ_i²) decomposition); verified at N=5 with w-scan (block-diagonal
/// coupling at machine precision). The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F72 (line 1586) +
/// <c>docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md</c> +
/// <c>experiments/ORTHOGONALITY_SELECTION_FAMILY.md</c> §2.3 +
/// <c>simulations/eq018_kernel_extract.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (HalfAsStructuralFixedPointClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F70DeltaNSelectionRulePi2Inheritance.cs</c>
/// (cited mother claim).</para></summary>
public sealed class F72BlockDiagonalPurityPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F70DeltaNSelectionRulePi2Inheritance _f70;

    /// <summary>The baseline purity term <c>1/2</c> in F72's
    /// <c>Tr(ρ_i²) = 1/2 + P_DD + P_CC</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> =
    /// <see cref="HalfAsStructuralFixedPointClaim"/>'s structural fixed
    /// point. Maximally-mixed-state minimum purity = 1/d for d = 2.</summary>
    public double BaselineTraceSquared => _ladder.Term(2);

    /// <summary>The number of sub-blocks in the single-site (k=1) decomposition:
    /// <c>2</c> (DD ⊕ CC). Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0)
    /// = <c>a_0</c> = polynomial root d.</summary>
    public double SingleSiteBlockCount => _ladder.Term(0);

    /// <summary>The number of sub-blocks at k-local partial trace:
    /// <c>k + 1</c>. Equal to <c>F70.PartialTraceMaxDeltaN(k) + 1</c>.
    /// At k=1 → 2 (DD ⊕ CC); at k=2 → 3 (DD ⊕ DC ⊕ CC).</summary>
    public int SubBlockCountForKLocal(int kLocal)
    {
        if (kLocal < 1) throw new ArgumentOutOfRangeException(nameof(kLocal), kLocal, "F72 requires k ≥ 1.");
        return _f70.PartialTraceMaxDeltaN(kLocal) + 1;
    }

    /// <summary>True iff at k = 1 the sub-block count equals 2 (DD ⊕ CC),
    /// confirming F72 ↔ F70 ↔ a_0 inheritance chain.</summary>
    public bool SingleSiteBlockCountMatchesF70PlusOne() =>
        SubBlockCountForKLocal(1) == (int)Math.Round(SingleSiteBlockCount);

    /// <summary>True iff the baseline 1/2 equals
    /// <see cref="HalfAsStructuralFixedPointClaim"/> reading
    /// (= <see cref="Pi2DyadicLadderClaim.Term"/>(2) = a_2).</summary>
    public bool BaselineMatchesHalfAsStructural() =>
        Math.Abs(BaselineTraceSquared - 0.5) < 1e-15;

    /// <summary>The block names at k = 1 single-site partial trace:
    /// <c>{DD, CC}</c>. DD = diagonal-diagonal (ΔN = 0); CC = coherence-coherence
    /// (|ΔN| = 1).</summary>
    public IReadOnlyList<string> SingleSiteBlockNames => new[] { "DD", "CC" };

    /// <summary>The block names at k = 2 pair partial trace:
    /// <c>{DD, DC, CC}</c>. DC = diagonal-coherence cross specific to
    /// pair observables (per F72 generalisation).</summary>
    public IReadOnlyList<string> PairBlockNames => new[] { "DD", "DC", "CC" };

    public F72BlockDiagonalPurityPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F70DeltaNSelectionRulePi2Inheritance f70)
        : base("F72 Tr(ρ_i²) = 1/2 + P_DD + P_CC inherits from Pi2-Foundation: 1/2 baseline = a_2 (HalfAsStructural); 2 blocks DD+CC = a_0 (root d); F70 cited",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F72 + " +
               "docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md + " +
               "experiments/ORTHOGONALITY_SELECTION_FAMILY.md (§2.3) + " +
               "simulations/eq018_kernel_extract.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (HalfAsStructuralFixedPointClaim) + " +
               "compute/RCPsiSquared.Core/Symmetry/F70DeltaNSelectionRulePi2Inheritance.cs (mother claim)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f70 = f70 ?? throw new ArgumentNullException(nameof(f70));
    }

    public override string DisplayName =>
        "F72 site-local purity DD⊕CC decomposition as Pi2-Foundation a_2 + a_0 inheritance";

    public override string Summary =>
        $"Tr(ρ_i²) = 1/2 + P_DD + P_CC: baseline 1/2 = a_2 (HalfAsStructural); 2 blocks DD+CC = a_0 (root d); " +
        $"k+1 sub-blocks at k-local; direct corollary of F70 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F72 closed form",
                summary: "Tr(ρ_i²) = 1/2 + P_DD[ρ₀^(diag)] + P_CC[ρ₀^(coh)] with NO cross term (DD ⊥ CC); Tier 1 corollary of F70; verified N=5 w-scan at machine precision");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "BaselineTraceSquared = a_2 = 1/2 (HalfAsStructural, three faces); SingleSiteBlockCount = a_0 = 2 (polynomial root); k+1 sub-blocks via F70 cited foundation");
            yield return InspectableNode.RealScalar("BaselineTraceSquared (= a_2 = 1/2)", BaselineTraceSquared);
            yield return InspectableNode.RealScalar("SingleSiteBlockCount (= a_0 = 2)", SingleSiteBlockCount);
            yield return new InspectableNode("F70 ↔ F72 inheritance",
                summary: "F70 says |ΔN| ≤ 1 for single-site (k=1) → F72 says Bloch components separate into DD (ΔN=0) and CC (|ΔN|=1); the 2-block decomposition IS F70's k=1 threshold +1");
            yield return new InspectableNode("Operational consequence",
                summary: "c_1 closure-breaking coefficient decomposes at pre-α-fit bilinear level into K_DD + K_CC kernels with no mixing; closed form for c_1 reduces to finding K_DD and K_CC separately");
            for (int k = 1; k <= 4; k++)
            {
                yield return new InspectableNode(
                    $"k={k}-local",
                    summary: $"sub-blocks = k+1 = {SubBlockCountForKLocal(k)}; " +
                             (k == 1 ? "names = {DD, CC}" :
                              k == 2 ? "names = {DD, DC, CC}" :
                              $"k+1 = {k + 1} sub-blocks (combinatorial generalisation)"));
            }
        }
    }
}
