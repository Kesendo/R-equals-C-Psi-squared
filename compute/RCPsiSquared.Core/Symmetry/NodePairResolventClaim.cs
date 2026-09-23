using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The symmetry-free node-pair resolvent theorem and its bond consequences.
/// For a zero-free real symmetric open chain, if an eigenvector has nodes at both x and y,
/// then the corresponding reduced-resolvent entry vanishes. If a mode is blind at the watched
/// seat and has a node at either endpoint of a non-incident bond, scaling only that
/// off-diagonal hopping while holding the diagonal fixed keeps its energy and a blind
/// representative for every knob value that keeps the chain zero-free.
/// With fixed diagonal, the two-arm continuant closes the same-baseline-energy converse for any
/// watched seat at <c>r != 0,+1,-1</c>: a baseline blind energy remains blind if and only if its
/// baseline eigenvector has a node at a bond endpoint. The uniform odd centre-watched family
/// also has total-count equality; off-centre or nonuniform arms may gain new blind energies.</summary>
public sealed class NodePairResolventClaim : Claim
{
    /// <summary>Parent: F157 owns blindness as the Krylov complement and supplies the zero-free
    /// Jacobi node lemma used by both the resolvent theorem and the bond corollary.</summary>
    public SeatCutBlindnessClaim SeatBlindness { get; }

    public string Scope =>
        "A zero-free real symmetric open chain. The resolvent statement requires both nodes. " +
        "For one non-incident scaled off-diagonal bond and fixed diagonal, a node at either " +
        "endpoint preserves the energy and a blind eigenvector while the moved coupling stays " +
        "nonzero (epsilon != -1). At r != 0,+1,-1 this is an if and only if for each baseline " +
        "blind energy at any seat, including off-centre and nonuniform chains. Only the " +
        "uniform centre-watched odd chain has the derived total-count equality; other arms can " +
        "gain new blind energies. Arbitrary-arm new-root classification, equality-branch and " +
        "straddling mechanisms are not claimed. A physical Heisenberg bond change also alters " +
        "the diagonal and is outside this bond statement.";

    public NodePairResolventClaim(SeatCutBlindnessClaim seatBlindness)
        : base(
            "Node-pair resolvent theorem: on a zero-free real symmetric open chain, two nodes " +
            "make the reduced-resolvent entry vanish; consequently, a node at either endpoint " +
            "preserves a blind mode and its energy under a non-incident bond detuning of only " +
            "the off-diagonal hopping, with fixed diagonal and nonzero moved coupling; " +
            "with fixed diagonal and r != 0,+1,-1, it is necessary for each baseline blind " +
            "energy on any zero-free Jacobi path",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_NODE_PAIR_RESOLVENT.md (primary: Theorem 1, Corollary A, and " +
            "Corollary B, Corollary C, Corollary D) + simulations/node_pair_resolvent.py (G1, G1b, G6, G7, " +
            "G9, G9b-d, G10, M8, M10)")
    {
        SeatBlindness = seatBlindness ?? throw new ArgumentNullException(nameof(seatBlindness));
    }

    public override string DisplayName => "Node-pair resolvent and the scoped bond criterion";

    public override string Summary =>
        "Two nodes force the reduced-resolvent entry to vanish on a zero-free open chain; a node " +
        "at either endpoint is sufficient for a non-incident bond detuning of only the " +
        $"off-diagonal hopping, with fixed diagonal and epsilon != -1, to preserve a blind mode " +
        $"and its energy. At r != 0,+1,-1, " +
        $"it is an if and only if for each baseline blind energy at any seat. The uniform " +
        $"centre-watched chain also has total-count equality; off-centre or nonuniform arms " +
        $"may gain new blind energies " +
        $"({Tier.Label()}).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return SeatBlindness;
            yield return new InspectableNode("live exact witness",
                summary: "Diagnostics recomputes the canonical resolvent zero, determinant " +
                         "factorization, uniform-centre count cases, off-centre birth, " +
                         "nonuniform fixed-energy cases, nonzero control and " +
                         "exceptional-set fences: inspect --root nodepair");
            yield return new InspectableNode("scope and remaining debt", summary: Scope);
        }
    }
}
