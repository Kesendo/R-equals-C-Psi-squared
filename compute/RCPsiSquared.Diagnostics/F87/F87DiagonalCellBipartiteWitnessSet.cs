using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 surface for the F103 §7 bipartite-chirality (diagonal-K) criterion as a
/// single registry-queryable claim. Wraps <see cref="F87DiagonalCellBipartiteWitness.StandardSet"/>
/// (4 canonical witnesses: soft, hard-odd-cycle, hard-pure-D-template, and the X-deph Z↔X mirror)
/// with typed parent edges to <see cref="F87TrichotomyClassification"/> (the F87 verdict the
/// criterion is checked against) and <see cref="Core.Symmetry.ChiralKClaim"/> (the chiral K that
/// certifies bipartite ⟹ soft). The individual parameterised
/// <see cref="F87DiagonalCellBipartiteWitness"/> stays a deferred type by design; this set is the
/// registered surface, mirroring how <see cref="F87StandardWitnessSet"/> surfaces the canonical
/// witnesses.
///
/// <para>The witnesses inside still classify lazily via <see cref="BipartiteChirality.Classify"/>;
/// this wrapper does NOT trigger classification at construction. Consumers that touch a witness's
/// <c>ActualClass</c> / <c>CriterionAgrees</c> pay the H/L-build cost at first access only.</para>
///
/// <para>Tier: Tier1Derived, matching the inner witnesses. <c>bipartite ⟹ soft</c> is derived
/// (chiral K); the converse <c>non-bipartite ⟹ hard</c> is derived modulo standard perturbation
/// theory (the first-order-block premise closed via degenerate-PT + analyticity, PROOF_F103
/// §7.5/§7.6). Promoted from Tier1Candidate in the formal promotion pass (2026-06-08).</para>
///
/// <para>"Derived modulo standard perturbation theory" means: derived at the genericity level
/// (non-bipartite ⟹ hard for all but finitely many γ). The remaining all-γ closure incl. the physical
/// point is isolated as <see cref="WindowedConverseAllGammaClaim"/>, closed 2026-06-10 with no residual
/// (Pascal-Gram positivity; the two-reflection spine is WindowedConverseThresholdClaim).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7 +
/// <see cref="Core.Symmetry.ChiralKClaim"/> + <see cref="PauliPairTrichotomy"/>; per-witness
/// anchors live on each inner <see cref="F87DiagonalCellBipartiteWitness.Anchor"/>.</para></summary>
public sealed class F87DiagonalCellBipartiteWitnessSet : Claim
{
    public ChainSystem Chain { get; }
    public IReadOnlyList<F87DiagonalCellBipartiteWitness> Witnesses { get; }

    public F87DiagonalCellBipartiteWitnessSet(ChainSystem chain)
        : base("F103 §7 bipartite-chirality witness set",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7 + " +
               "compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs + " +
               "compute/RCPsiSquared.Diagnostics/F87/PauliPairTrichotomy.cs + " +
               "simulations/f87_42_8_bipartite_fullcell.py + " +
               "simulations/f87_bipartite_chiral_witness.py")
    {
        Chain = chain ?? throw new System.ArgumentNullException(nameof(chain));
        Witnesses = F87DiagonalCellBipartiteWitness.StandardSet(chain);
    }

    public override string DisplayName =>
        $"F103 §7 diagonal-K (bipartite-chirality) witness set (N={Chain.N}, {Witnesses.Count} witnesses)";

    public override string Summary =>
        $"{Witnesses.Count} diagonal-cell witnesses (soft / hard-odd-cycle / hard-template / X-deph mirror); " +
        $"{PassCount()}/{Witnesses.Count} PASS ({Tier.Label()})";

    private int PassCount() => Witnesses.Count(w => w.Matches);

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("chain",
                summary: $"N={Chain.N}, {Chain.HType}, {Chain.Topology}");
            foreach (var w in Witnesses)
                yield return w;   // F87DiagonalCellBipartiteWitness : Claim : IInspectable
        }
    }
}
