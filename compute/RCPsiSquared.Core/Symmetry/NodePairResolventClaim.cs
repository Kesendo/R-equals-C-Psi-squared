using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The symmetry-free node-pair resolvent theorem and its bond consequences.
/// For a zero-free real symmetric open chain, if an eigenvector has nodes at both x and y,
/// then the corresponding reduced-resolvent entry vanishes. If a mode is blind at the watched
/// seat and has a node at either endpoint of a non-incident moved bond, that mode keeps its
/// energy and a blind representative for every knob value that keeps the chain zero-free.
/// On the uniform odd chain watched at its centre, the exact half-chain determinant closes the
/// pointwise converse for a non-incident bond at every ratio <c>r != 0,+1,-1</c>: a baseline blind
/// energy remains blind if and only if its baseline eigenvector has a node at a bond endpoint.
/// The off-centre and nonuniform pointwise converses remain open; the all-knob fixed-energy
/// necessity is a separate quantifier.</summary>
public sealed class NodePairResolventClaim : Claim
{
    /// <summary>Parent: F157 owns blindness as the Krylov complement and supplies the zero-free
    /// Jacobi node lemma used by both the resolvent theorem and the bond corollary.</summary>
    public SeatCutBlindnessClaim SeatBlindness { get; }

    public string Scope =>
        "A zero-free real symmetric open chain. The resolvent statement requires both nodes. " +
        "The bond statement is a sufficient direction for a non-incident bond whose moved " +
        "coupling stays nonzero (epsilon != -1): a node at either endpoint preserves the energy " +
        "and a blind eigenvector. In the uniform centre-watched odd-chain family, for " +
        "r != 0,+1,-1, this is also an if and only if for each baseline blind energy. The " +
        "off-centre and nonuniform pointwise converses remain open; equality-branch " +
        "and straddling mechanisms are not claimed.";

    public NodePairResolventClaim(SeatCutBlindnessClaim seatBlindness)
        : base(
            "Node-pair resolvent theorem: on a zero-free real symmetric open chain, two nodes " +
            "make the reduced-resolvent entry vanish; consequently, a node at either endpoint " +
            "is a sufficient condition for a non-incident nonzero bond detuning to preserve a " +
            "blind mode and its energy; on the uniform centre-watched odd chain, outside " +
            "r=0,+1,-1, that endpoint-node condition is also necessary pointwise",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_NODE_PAIR_RESOLVENT.md (primary: Theorem 1, Corollary A, and " +
            "Corollary B, and Corollary C) + simulations/node_pair_resolvent.py (G1, G1b, G6, G7, " +
            "G9, G9b-d, M8)")
    {
        SeatBlindness = seatBlindness ?? throw new ArgumentNullException(nameof(seatBlindness));
    }

    public override string DisplayName => "Node-pair resolvent and the scoped bond criterion";

    public override string Summary =>
        "Two nodes force the reduced-resolvent entry to vanish on a zero-free open chain; a node " +
        "at either endpoint is sufficient for a non-incident bond detuning with epsilon != -1 " +
        $"to preserve a blind mode and its energy. On the uniform centre-watched odd chain " +
        $"with r != 0,+1,-1 it is an if and only if; the off-centre and nonuniform pointwise " +
        $"converses remain open " +
        $"({Tier.Label()}).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return SeatBlindness;
            yield return new InspectableNode("live exact witness",
                summary: "Diagnostics recomputes the canonical resolvent zero, determinant " +
                         "factorization, uniform-centre iff cases, nonzero control and " +
                         "exceptional-set fences: inspect --root nodepair");
            yield return new InspectableNode("scope and remaining debt", summary: Scope);
        }
    }
}
