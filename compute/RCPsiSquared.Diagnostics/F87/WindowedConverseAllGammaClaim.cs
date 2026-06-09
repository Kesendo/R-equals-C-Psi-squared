using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F87 windowed-converse residual lemma. Phase A typed it as the open all-γ closure of
/// the windowed converse; Phase B (2026-06-09) reduced it to a positive-monomial certificate and
/// proved the spine, leaving two sharp residuals. It stays Tier1Candidate.
///
/// <para>PROVEN (Tier1Derived, via <see cref="WindowedConverseThresholdClaim"/> and
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/>): the two-reflection spine (𝓕=F⊗F, R=I⊗F ⟹
/// #A_L,#A_R,#Q all odd ⟹ bipartite soft + non-bipartite #A≥2ℓ threshold), the monomial structure,
/// and deg-1 positivity. The first nonvanishing odd power-sum of M=A+γQ is, for every hard pair, a
/// positive monomial c·γ^deg (deg in {1,3}, m* = 2ℓ+deg), verified bit-exact cell-wide at N=4 and at
/// N=5/N=6 reps. A positive monomial has no positive real root, so hard at every γ>0.</para>
///
/// <para>OPEN (the two residuals this lemma now reduces to, proven modulo them): R-deg and R-sign
/// (see the two open-work nodes). Closing both makes the windowed converse non-bipartite ⟹ hard a
/// closed-form general-N theorem and promotes F110/F111 to Tier1Derived.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> +
/// <c>simulations/f87_windowed_monomial_converse.py</c> +
/// <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.5/§7.7</c>.</para></summary>
public sealed class WindowedConverseAllGammaClaim : Claim
{
    public WindowedConverseAllGammaClaim()
        : base("F87 windowed converse, all-γ residual: non-bipartite ⟹ hard at every γ>0, the first nonvanishing odd power-sum is a positive monomial; proven modulo R-deg + R-sign",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md + " +
               "simulations/f87_windowed_monomial_converse.py + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.5/§7.7")
    {
    }

    /// <summary>The closing certificate.</summary>
    public string Theorem =>
        "For a non-bipartite windowed diagonal-cell pair, the first nonvanishing odd power-sum of M = A + γQ is a " +
        "positive monomial c·γ^deg (deg in {1,3}, m* = 2ℓ + deg). No positive real root ⟹ hard at every γ>0.";

    /// <summary>The proven spine (Tier1Derived, via WindowedConverseThresholdClaim).</summary>
    public string ProvenSpine =>
        "DERIVED (WindowedConverseThresholdClaim): the two reflections 𝓕=F⊗F, R=I⊗F force #A_L,#A_R,#Q all odd in " +
        "every odd power-sum word, giving bipartite ⟹ soft and non-bipartite ⟹ #A ≥ 2ℓ; the monomial structure and " +
        "deg-1 positivity are closed-form.";

    /// <summary>Open residual 1.</summary>
    public string RDeg =>
        "R-deg (OPEN; verified bit-exact cell-wide N=4, reps N=5/N=6): Tr(Q·A^{2ℓ}) on the (#A_L=ℓ,#A_R=ℓ) " +
        "odd-cycle-traversal class sums to zero for pure off-diagonal H, so p_{2ℓ+1} ≡ 0 and deg rises to 3 for " +
        "genuine cycles (the monomial at m* = 2ℓ+3). No uniform closed-form involution yet.";

    /// <summary>Open residual 2.</summary>
    public string RSign =>
        "R-sign (OPEN; verified 16/16 pure cycles at N=4): P_{m*,3} > 0 for genuine cycles, exactly the §7.5 " +
        "+N-population-Perron top-skew (the −N reflection mode absent). Not yet a closed-form nonneg identity.";

    /// <summary>The Tier1Candidate claims gated on closing both residuals.</summary>
    public string Consumers =>
        "F110 (HardCellYInversionPattern) and F111 (HardCellPureDTemplate) are gated on this lemma; closing R-deg + " +
        "R-sign promotes both from Tier1Candidate to Tier1Derived.";

    public override string DisplayName =>
        "F87 windowed converse, all-γ residual (Tier1Candidate, proven modulo R-deg + R-sign)";

    public override string Summary =>
        "the first nonvanishing odd power-sum is a positive monomial ⟹ hard ∀γ>0; the two-reflection spine + " +
        "monomial structure + deg-1 positivity are PROVEN (WindowedConverseThresholdClaim); the full theorem is " +
        $"proven modulo two open residuals R-deg + R-sign ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem (positive-monomial certificate)", summary: Theorem);
            yield return new InspectableNode("Proven spine (Tier1Derived)", summary: ProvenSpine);
            yield return new InspectableNode("Open residual R-deg", summary: RDeg);
            yield return new InspectableNode("Open residual R-sign", summary: RSign);
            yield return new InspectableNode("Consumers (gated on closing both residuals)", summary: Consumers);
        }
    }
}
