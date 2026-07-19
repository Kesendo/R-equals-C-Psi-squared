using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F136, typed (Tier 1 derived): <b>the record letter law</b>. WHICH operator a witness
/// records is the parity of its shared neighborhood, derived channel-by-channel from F135's
/// Proposition 1 (<see cref="RecordParityLawClaim"/>). Partition the sites around the pair
/// (S, j) into <b>D</b> (shared dressers, m = |D|), <b>P</b> (private watchers of j), <b>Q</b>
/// (private S-neighbors); ratios r_k = Δ_jk/Δ_S. At t* = π/(4Δ_S), with the V_S = 0 hinge
/// (S keeps a neighbor besides j), exactly two mutually exclusive perfect records exist:
///
/// <code>
///   pointer record:  write bond present, deg(S) ≥ 2, EVERY watcher of j even
///                    → 1 bit of Z_S in j's equator (the ZY channel)
///   Bell record:     Q = ∅, m ≥ 1, all D odd, all P even
///                    → ρ_Sj = ½|Φ^{σ₁}⟩⟨Φ^{σ₁}| + ½|Ψ^{σ₂}⟩⟨Ψ^{σ₂}|, ZERO pointer content,
///                      the LETTER is the dresser parity: m odd → Y⊗Y, m even → X⊗X
///                      (sign(c₁c₂) = (−1)^m)
/// </code>
///
/// Aliveness is watcher parity: private watchers must be even for either family; the shared
/// dressers' parity selects it (all even + write bond → pointer, all odd → Bell, mixed →
/// dark); every integer violation is an exact kill. Signs are closed forms
/// (σ₁ = Π_D(−1)^{(1+r)/2}·Π_P(−1)^{r/2}, σ₂ = Π_D(−1)^{(1−r)/2}·Π_P(−1)^{r/2}; the pointer
/// sign is Π_{D∪P}(−1)^{r/2}). Dephasing dressing: the Bell record pays BOTH sites
/// (I = 1 − h₂((1+κ)/2), κ = e^{−2(γ_S+γ_j)t*}), the pointer record only γ_j, and γ on any
/// traced site is exactly invisible. A pure pendant S ROLE-SWAPS into j's witness (1 separable
/// YZ bit of Z_j at an odd watcher). Fan-out corollary: Bell witnesses need not be neighbors,
/// so anti-pointer redundancy is NOT deg-bounded (K_{R+1,2}: a Bell-record clique at
/// deg(S) = 2). Graph-level corollaries 7 + 8: at uniform coupling a pair is luminous iff
/// neighborhoods match or one member is the other's leaf, so the FULLY-witnessed worlds are
/// exactly the stars and the complete graphs (census 38 + 728 + 26704 connected graphs at
/// N = 4/5/6, winners the N labeled stars + K_N); and girth ≥ 5 + leafless ⇒ every pair dark
/// (the heavy-hex bulk holds zero luminous pairs; the square-lattice bulk is dark too, only
/// the isolated C₄ weaves).
///
/// <para><b>Typed parent.</b> <see cref="RecordParityLawClaim"/> (F135): every channel above
/// is Proposition 1 read on a channel class; Law A is the D = ∅ column of the trichotomy.
/// At uniform coupling the statics reduce to Hein-Eisert-Briegel 2004 composed with the
/// degree rotation (credited in the proof); the ratio arithmetic, signs, exclusivity and γ
/// dressing are the arc's own. Gates: <c>simulations/qd_letter_gates.py</c> (87/87
/// pre-registered), <c>simulations/qd_witness_play.py</c> (the census),
/// <c>simulations/qd_heavyhex_map.py</c> (the dark map). Live: <c>inspect --root record</c>
/// (<c>RecordLawWitness</c>). The MirrorWorld adoption (<c>Witness.cs</c>, run mode
/// <c>witness</c>) carries the same closed forms in the sober base.</para></summary>
public sealed class RecordLetterLawClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: every channel is F135's Proposition 1 read on a channel class.
    public RecordParityLawClaim ParityLaw { get; }

    public RecordLetterLawClaim(RecordParityLawClaim parityLaw)
        : base("The record letter law (F136): which operator a witness records is the parity of " +
               "its shared neighborhood -- private watchers must be even for any record; the " +
               "shared dressers' parity selects the family (all even + write bond -> pointer " +
               "record of Z_S, all odd -> Bell record with zero pointer content, mixed -> dark); " +
               "the Bell letter alternates Y/X with the dresser count (sign(c1c2) = (-1)^m); " +
               "signs are closed forms; the Bell record pays both sites, the pointer only the " +
               "witness, traced sites are free; a pendant S role-swaps; anti-pointer redundancy " +
               "is not deg-bounded; fully-witnessed worlds are exactly the stars and K_N and " +
               "girth >= 5 + leafless is dark",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_RECORD_LETTER_LAW.md + " +
               "experiments/QUANTUM_DARWINISM_POINTER_DOOR.md")
    {
        ParityLaw = parityLaw ?? throw new ArgumentNullException(nameof(parityLaw));
    }

    public override string DisplayName =>
        "Record letter law (F136): which operator a witness records is its shared-neighborhood parity";

    public override string Summary =>
        "two mutually exclusive perfect records at t*: the pointer record (write bond, all watchers " +
        "even; 1 bit of Z_S in ZY) and the Bell record (Q = ∅, dressers odd, privates even; " +
        "½Φ + ½Ψ, letter = dresser parity, m odd → YY / m even → XX); mixed parities are exactly " +
        "dark; Bell pays both sites, pointer only the witness, watching the writers is free; " +
        "fully-witnessed worlds = stars + K_N, girth ≥ 5 + leafless = dark " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the two families and their exclusivity",
                summary: "a perfect pointer record needs all of j's watchers even; a perfect Bell " +
                         "record needs the shared ones odd (with Q = ∅, privates even); no watcher is " +
                         "both, so a witness never records both. Every integer violation (odd private, " +
                         "mixed D, nonempty Q for Bell, missing write bond for pointer) is an exact kill.");
            yield return new InspectableNode("the letter is the dresser parity",
                summary: "m odd → ⟨Y_SY_j⟩ = ±1, m even → ⟨X_SX_j⟩ = ±1; sign(c₁c₂) = (−1)^m: the m-th " +
                         "shared dresser flips which Bell pair carries the bit; the direct S–j bond " +
                         "never moves the letter. The signs walk the closed σ-formulas (K₂,₁ at r = 3 " +
                         "reads YY = −1); the stacked-parity shape the one-four thesis predicted.");
            yield return new InspectableNode("the γ dressing: two record species, two prices",
                summary: "the Bell record decays at the SUM of the pair's watching rates " +
                         "(κ = e^{−2(γ_S+γ_j)t*}), the pointer record only at the witness's own rate " +
                         "(γ_S exactly invisible to it); no third site's watching ever touches either: " +
                         "an anti-pointer testimony pays double, and watching the messengers is free.");
            yield return new InspectableNode("the pendant role-swap and the fan-out",
                summary: "a pure pendant S (hypothesis H excludes it) reads BACKWARDS: with an odd " +
                         "watcher on j the pair holds j's pointer Z_j in S's equator (YZ), S the " +
                         "better-protected witness. And Bell witnesses need not be neighbors: on " +
                         "K_{R+1,2} every corner holds a perfect X⊗X bit at deg(S) = 2, a Bell-record " +
                         "clique; pointer redundancy is bounded by the bonds S owns, anti-pointer " +
                         "redundancy only by who shares its dressers.");
            yield return new InspectableNode("the graph-level corollaries (7 + 8)",
                summary: "uniform coupling: luminous ⟺ matched neighborhoods or one member the other's " +
                         "leaf; hence fully-witnessed worlds = exactly the stars and the complete " +
                         "graphs (census 38 + 728 + 26704 connected graphs at N = 4/5/6, winners the " +
                         "N labeled stars + K_N; the star shows all three readings, K_N weaves with " +
                         "zero pointer content); and girth ≥ 5 + leafless ⇒ every pair dark: the " +
                         "heavy-hex bulk is the anti-star, all private rooms; the square-lattice bulk " +
                         "is dark too (only the isolated C₄ weaves); on a dark lattice light can be " +
                         "aimed but not kept private.");
            yield return new InspectableNode("typed parent",
                summary: $"RecordParityLawClaim ({ParityLaw.Tier.Label()}): every channel is " +
                         "Proposition 1 read on a channel class; Law A is the D = ∅ column. At uniform " +
                         "coupling the statics are HEB 2004 + the degree rotation (credited); the ratio " +
                         "arithmetic, signs, exclusivity and γ dressing are the arc's own.");
            yield return new InspectableNode("live witness (inspect --root record)",
                summary: "RecordLawWitness recomputes the pair battery at inspect time: pointer / Bell " +
                         "/ role-swap / dark verdicts with letters, signs and prices from the closed " +
                         "forms, each checked against an independent full-state path (2^N closed-form " +
                         "density matrix, partial trace, measured I and correlators).");
        }
    }
}
