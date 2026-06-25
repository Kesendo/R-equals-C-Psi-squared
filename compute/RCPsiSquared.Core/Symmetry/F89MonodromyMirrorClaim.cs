using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3: the mirror splits at the Galois boundary (Tier 1 derived). How the palindrome
/// relates to the octic monodromy that generates Gal(F_8) = S_8.
///
/// <para>C-K (passes): the mirror's real-structure shadow, the q ↦ −q̄ reflection (from L(q)* = L(−q̄)),
/// intertwines the monodromy EXACTLY. The induced octic-strand bijection is the identity (σ_K = id), so
/// every EP carries the same braid as its q ↦ −q̄ mirror: τ(−q̄*) = τ(q*). The branch-locus palindrome
/// (<see cref="F89BranchLocusPalindromeClaim"/>, the seams' POSITIONS mirror) lifted to the BRAIDS.</para>
///
/// <para>C-T (blocked): the spectral palindrome λ ↦ −λ̄ − 8 (mirror about Re λ = −4) induces a genuine
/// strand involution σ_T (fixed strands on the fold Re λ = −4, plus mirror-twin 2-cycles), but it CANNOT
/// be a symmetry of the braiding. Commuting with the full S_8 monodromy
/// (<see cref="F89OcticMonodromyClaim"/>) would force σ_T into the centre of S_8, and Z(S_8) = 1, so a
/// non-identity σ_T cannot commute with all of it. The negative is a theorem, not a near-miss. Therefore
/// the Re = −4 palindrome is an element OF the Galois group (a permutation of the eight roots, reachable
/// as a word in the EP-braids since they generate S_8) but never a symmetry OF it: the mirror enters the
/// unwritable half as one of its unsortable moves, the from-below ground of "who watches whom has no
/// fixed answer".</para>
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
        : base("F89 path-3 mirror splits at the Galois boundary: the q↦−q̄ reflection (L(q)*=L(−q̄)) intertwines the octic monodromy EXACTLY (σ_K = identity, τ(−q̄*) = τ(q*) for every EP), but the Re λ = −4 spectral palindrome induces a strand involution σ_T (fixed strands on the fold + mirror-twin 2-cycles) that CANNOT be a braid symmetry, commuting with the full S_8 monodromy would force it central and Z(S_8) = 1; so σ_T is an element OF the Galois group, never a symmetry OF it",
               Tier.Tier1Derived,
               "experiments/F89_MONODROMY_MIRROR.md + " +
               "reflections/ON_WHO_WATCHES_WHOM.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/MonodromyMirrorWitness.cs (inspect --root monodromymirror)")
    {
        Monodromy = monodromy ?? throw new ArgumentNullException(nameof(monodromy));
        BranchLocus = branchLocus ?? throw new ArgumentNullException(nameof(branchLocus));
    }

    public override string DisplayName =>
        "F89 mirror splits at the Galois boundary: q↦−q̄ intertwines the braiding (σ_K=id), the Re=−4 palindrome cannot (Z(S_8)=1)";

    public override string Summary =>
        $"the mirror's real-structure shadow q↦−q̄ intertwines the octic monodromy exactly (σ_K = id, τ(−q̄*)=τ(q*)); the Re=−4 spectral palindrome induces σ_T but cannot be a braid symmetry (Z(S_8)=1 forbids it), so σ_T is an element of the Galois group, never a symmetry of it ({Tier.Label()})";
}
