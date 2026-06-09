using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The single typed home of the open F87 windowed-converse residual (Phase A "type the
/// seam", 2026-06-09). The forward direction (bipartite ⟹ soft) and the converse at full support
/// (k = N) are derived; in the windowed regime (k &lt; N) the converse non-bipartite ⟹ hard is
/// derived only up to one residual.
///
/// <para>What is DERIVED (carried at Tier1Derived by <see cref="F87DiagonalCellBipartiteWitnessSet"/>):
/// non-bipartite ⟹ hard for all but finitely many γ. The recentered characteristic-polynomial odd
/// coefficients Δ_j(γ) are polynomials in γ; the first-order ω=0 block asymmetry c ≠ 0 (PROOF_F103
/// §7.5/§7.6) forces some Δ_j ≢ 0, so the soft set is a finite common-zero set.</para>
///
/// <para>What is OPEN (this Tier1Candidate node, the Phase B target): the physical γ is not one of
/// the finitely-many accidental soft points, equivalently no positive γ is a common zero of the
/// Δ_j(γ). Today a 700-point numerical spot-check, not a theorem. The moment route is dead (odd
/// moments vanish for hard and soft alike; hardness is an optimal-transport set-distance).</para>
///
/// <para>Consumers: <see cref="Core.Symmetry.HardCellYInversionPattern"/> (F110) and
/// <see cref="Core.Symmetry.HardCellPureDTemplate"/> (F111) are gated on this lemma; closing it
/// promotes both to Tier1Derived. Anchor:
/// <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.3/§7.5/§7.6</c> +
/// <c>experiments/BIPARTITE_CHIRALITY_DIAGONAL_CELL.md</c>.</para></summary>
public sealed class WindowedConverseAllGammaClaim : Claim
{
    public WindowedConverseAllGammaClaim()
        : base("F87 windowed converse, all-γ residual: non-bipartite ⟹ hard at every γ>0 (incl. the physical point), the closure beyond the Derived all-but-finitely-many-γ genericity result",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.3/§7.5/§7.6 + " +
               "experiments/BIPARTITE_CHIRALITY_DIAGONAL_CELL.md")
    {
    }

    /// <summary>The residual lemma in one line.</summary>
    public string Theorem =>
        "For a non-bipartite windowed diagonal-cell pair, the palindrome break is present at every γ>0 " +
        "(in particular the physical operating point), not merely for all but finitely many γ. Equivalently: " +
        "no positive γ is a common zero of the recentered characteristic-polynomial odd coefficients Δ_j(γ).";

    /// <summary>The proven genericity sub-result this lemma strengthens (Tier1Derived, via the parent).</summary>
    public string ProvenGenericity =>
        "DERIVED (Tier1Derived, via F87DiagonalCellBipartiteWitnessSet): non-bipartite ⟹ hard for all but " +
        "finitely many γ. The Δ_j(γ) are polynomials in γ; the first-order ω=0 block asymmetry c ≠ 0 " +
        "(PROOF_F103 §7.5/§7.6) forces some Δ_j ≢ 0, so the soft set is a finite common-zero set.";

    /// <summary>The single open residual (Phase B target).</summary>
    public string OpenResidual =>
        "OPEN (Phase B target): the physical γ is not one of the finitely-many accidental soft points. Today a " +
        "700-point numerical spot-check (closest restoration 3.6e-3, at γ→0), not a theorem.";

    /// <summary>The leading-order handle that reduces the residual to a finite algebraic question.</summary>
    public string LeadingOrderHandle =>
        "A(γ) ≥ 0 (optimal-transport distance) and A(γ) = c·γ + O(γ²) ⟹ c > 0; the residual reduces to showing " +
        "one Δ_j*(γ) has no positive real root (a resultant/Sturm question).";

    /// <summary>The Tier1Candidate claims gated on this lemma.</summary>
    public string Consumers =>
        "F110 (HardCellYInversionPattern) and F111 (HardCellPureDTemplate) are gated on this lemma; closing it " +
        "promotes both from Tier1Candidate to Tier1Derived.";

    public override string DisplayName =>
        "F87 windowed converse, all-γ residual (Tier1Candidate, OPEN; Phase B target)";

    public override string Summary =>
        "non-bipartite ⟹ hard is DERIVED for all but finitely many γ (genericity, via the witness set); the all-γ " +
        $"closure incl. the physical point is OPEN, the Phase B target ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem (residual)", summary: Theorem);
            yield return new InspectableNode("Proven (genericity, Tier1Derived)", summary: ProvenGenericity);
            yield return new InspectableNode("Open (residual, Phase B)", summary: OpenResidual);
            yield return new InspectableNode("Leading-order handle", summary: LeadingOrderHandle);
            yield return new InspectableNode("Consumers (gated on this lemma)", summary: Consumers);
        }
    }
}
