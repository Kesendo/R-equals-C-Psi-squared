using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3: the mirror splits at the Galois boundary (Tier 1 derived). How the palindrome
/// relates to the octic monodromy that generates Gal(F_8) = S_8.
///
/// <para>C-K (the base-space face passes): the q ↦ −q̄ reflection (from L(q)* = L(−q̄)) is an exact family
/// symmetry, so it automatically intertwines the monodromy; in the aligned strand labelling the induced
/// bijection comes out as the identity (σ_K = id), so each cluster EP carries the same braid as its
/// q ↦ −q̄ mirror: τ(−q̄*) = τ(q*). The branch-locus palindrome
/// (<see cref="F89BranchLocusPalindromeClaim"/>, the seams' POSITIONS mirror) lifted to the BRAIDS.</para>
///
/// <para>C-T (the fibre face does not): the spectral fold λ ↦ −λ̄ − 8 (mirror about Re λ = −4, exact at the
/// real base q where σ_T is built) induces a genuine strand involution σ_T (four fixed strands on the fold
/// Re λ = −4, plus two mirror-twin 2-cycles). It is NOT a loop-independent symmetry of the braiding: that
/// would require σ_T to commute with the full S_8 monodromy (<see cref="F89OcticMonodromyClaim"/>), forcing
/// σ_T into the centre of S_8, but Z(S_8) = 1 and σ_T ≠ id. So σ_T is NON-CENTRAL: it does not commute with
/// the monodromy, and the strand-pairing it encodes is braided away around the seams (which strand is whose
/// mirror is not fixed). The negative is a theorem, not a near-miss. (Conjugation by σ_T is still a
/// nontrivial inner automorphism of S_8, so the mirror does act on the tangle; the point is that it is
/// non-central, not that it induces no automorphism. The strictly stronger claim that σ_T permutes the
/// EP-transpositions among themselves is numerically false over the cluster set, but that is not what
/// Z(S_8) = 1 settles.) This is the from-below ground of "who watches whom has no fixed answer".</para>
///
/// <para>Live witness <c>inspect --root monodromymirror</c> (<c>MonodromyMirrorWitness</c>): the spectral
/// sanities, the per-EP q ↦ −q̄ intertwining (σ_K = id), and σ_T's fixed-and-twin structure. Anchors:
/// <c>experiments/F89_MONODROMY_MIRROR.md</c> + <c>reflections/ON_WHO_WATCHES_WHOM.md</c>.</para></summary>
public sealed class F89MonodromyMirrorClaim : Claim
{
    // Parent-edge marker (the S_8 monodromy this mirror analysis acts on).
    public F89OcticMonodromyClaim Monodromy { get; }
    // Parent-edge marker (the branch-locus palindrome, the seams' position mirror this lifts to braids).
    public F89BranchLocusPalindromeClaim BranchLocus { get; }

    public F89MonodromyMirrorClaim(F89OcticMonodromyClaim monodromy, F89BranchLocusPalindromeClaim branchLocus)
        : base("F89 path-3 mirror splits at the Galois boundary: the q↦−q̄ reflection (L(q)*=L(−q̄)) is an exact family symmetry so it intertwines the octic monodromy (σ_K = identity in the aligned labelling, τ(−q̄*) = τ(q*) for the cluster EPs), but the Re λ = −4 spectral fold induces a strand involution σ_T (four fixed on the fold + two mirror-twin 2-cycles) that is NON-CENTRAL: commuting with the full S_8 monodromy would force it central and Z(S_8) = 1, so σ_T does not commute with the braiding (not a loop-independent symmetry; conjugation by it is still a nontrivial inner automorphism of S_8)",
               Tier.Tier1Derived,
               "experiments/F89_MONODROMY_MIRROR.md + " +
               "reflections/ON_WHO_WATCHES_WHOM.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/MonodromyMirrorWitness.cs (inspect --root monodromymirror)")
    {
        Monodromy = monodromy ?? throw new ArgumentNullException(nameof(monodromy));
        BranchLocus = branchLocus ?? throw new ArgumentNullException(nameof(branchLocus));
    }

    public override string DisplayName =>
        "F89 mirror splits at the Galois boundary: q↦−q̄ intertwines the braiding (σ_K=id), the Re=−4 fold is non-central (Z(S_8)=1)";

    public override string Summary =>
        $"the mirror's base-space face q↦−q̄ intertwines the octic monodromy (σ_K = id; forced by L(q)*=L(−q̄)); the Re=−4 spectral fold induces a non-central involution σ_T which does not commute with the monodromy (Z(S_8)=1), so it is not a loop-independent symmetry of the braiding (conjugation by it is still an inner automorphism of S_8) ({Tier.Label()})";
}
