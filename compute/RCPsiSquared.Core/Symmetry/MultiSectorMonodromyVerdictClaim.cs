using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Multi-Sector Monodromy verdict (Tier 1 candidate; computationally exact and gate-tested, the
/// free-fermion-additivity mechanism confirmed for the embedding, the λ-value construction open):
///
/// <para><b>Is the S₈ Galois braid the F89 (1,2) octic carries LOCALIZED to the (1,2) coherence sector, or SHARED
/// across the joint-popcount sectors of L(q)? The answer is N-DEPENDENT.</b></para>
///
/// <list type="bullet">
///   <item><b>N=4: CONFINED.</b> The braid-carrying defective √-EP is shared across exactly the D₄ orbit
///   {(1,2),(2,1),(2,3),(3,2)} of (1,2) under transpose + the bra/ket bit-flips; it does NOT reach the dense
///   half-filled core (2,2) (a Door-C echo: the braid lives on the dilute orbit, not the dense core).</item>
///   <item><b>N=5: SPREADS.</b> The braid reaches a symmetric 12-sector diamond
///   {(1,2),(1,3),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,2),(4,3)}, INCLUDING the dense cores (2,2),
///   (3,3), splitting into two cross-fold-conjugate families of 6, each carrying a BYTE-IDENTICAL shared
///   eigenvalue λ (same branch point, to 6 digits across the family, not merely the same gap).</item>
/// </list>
///
/// <para><b>Membership + the N=4 self-fold.</b> The braid set is the rule <c>{|bra−ket| = 1, both popcounts ∈
/// [1, N−1]} ∪ its q̃ ↦ N−q̃ cross-fold image</c>. The N-dependence is exactly the N=4 self-fold: the cross-fold
/// q̃ ↦ N−q̃ is a SELF-MAP on the |Δ|=1 orbit only at N=4 (where N−2 = 2, so (1,2) → (1,2)), collapsing the two
/// families onto one four-sector orbit (confined); at N ≥ 5 the self-fold lifts, Family B separates from Family A
/// and brings in the dense cores (spread). This is the same degenerate self-fold that runs through F89
/// (<see cref="F89CrossFoldSimilarityClaim"/> / the path-3 anchor), read on the braid.</para>
///
/// <para><b>Mechanism: free-fermion / AT additivity.</b> The elementary EP is a property of the |bra−ket|=1 SE-DE
/// coherence rung; a DIAGONAL mode-spectator (one more excitation on both bra and ket, the same single-particle
/// mode) leaves E_bra − E_ket and the whole EP (λ, character, gap) invariant, so (1,2) ≡ (2,3) ≡ (3,4)
/// byte-identically (verified to ~1e-14). Family B is the <see cref="F89CrossFoldSimilarityClaim"/> image
/// λ ↦ −λ̄ − 2N. So the SAME defective eigenvalue lives in every diamond sector: shared spectral content, a
/// symmetry broader than the naive Klein-four, the sector-resolved face of the free-fermion-additivity codim-1
/// collapse.</para>
///
/// <para><b>The λ value, from below.</b> The (1,2) block is the linear pencil L(q) = A + q·C, A = −2·diag(n_diff)
/// real (the AT rates, n_diff ∈ {1,3}) and C the ANTI-HERMITIAN free-fermion hopping. Since C is anti-Hermitian, at
/// real q the coherent term feeds only Im λ, so <b>Re λ = −2·⟨n_diff⟩_v exactly</b> (machine zero) for the
/// coalescing eigenvector: the λ VALUE is the AT Theorem-2 rate of a q-tuned MIXTURE of the rate-2 (overlap) and
/// rate-6 (no-overlap) eigenmodes (⟨n_diff⟩ = 2.31 at q=0.6209, 1.90 at q=1.078; the N=4 EP sits at the −4 midpoint,
/// ⟨n_diff⟩ = 2). Not off the AT theorem, only off the integer-quantized lines; the defectiveness is the
/// eigenvector coalescence (Jordan), separate from the rate.</para>
///
/// <para><b>What stays open (why this is Tier 1 CANDIDATE, not derived).</b> A CLOSED FORM for the mixture
/// ⟨n_diff⟩(q) (predicting it from the free-fermion mode geometry without diagonalizing the block) and the
/// codim-1-by-additivity theorem, plus general N: the verdict is verified only at N=4, 5, and the octic is S₈
/// (non-solvable in radicals), so no global radical form for the P₁₀ EP loci is expected.</para>
///
/// <para><b>Typed parents.</b> <see cref="F89OcticMonodromyClaim"/> (the S₈ braid this census spreads across
/// sectors) and <see cref="F89CrossFoldSimilarityClaim"/> (F89d, the exact antiunitary similarity λ ↦ −λ̄ − 2N
/// that pairs Family B with Family A and whose N=4 self-fold IS the confinement). Live:
/// <c>inspect --root sectorbraid</c>. Anchor: <c>experiments/F89_MULTI_SECTOR_MONODROMY.md</c>; the exact N=5
/// 12-set is the regression test <c>Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore</c>.</para></summary>
public sealed class MultiSectorMonodromyVerdictClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89OcticMonodromyClaim Octic { get; }
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89CrossFoldSimilarityClaim CrossFold { get; }

    public MultiSectorMonodromyVerdictClaim(F89OcticMonodromyClaim octic, F89CrossFoldSimilarityClaim crossFold)
        : base("The multi-sector monodromy verdict is N-dependent: the S₈ braid the F89 (1,2) octic carries is " +
               "CONFINED to the D₄ orbit {(1,2),(2,1),(2,3),(3,2)} at N=4 (the dense core (2,2) braid-free) but " +
               "SPREADS at N=5 to a symmetric 12-sector joint-popcount diamond (incl. the dense cores (2,2),(3,3)), " +
               "two cross-fold-conjugate families of 6 sharing a byte-identical eigenvalue λ. Membership = " +
               "{|bra−ket|=1, popcounts ∈ [1,N−1]} ∪ q̃↦N−q̃ cross-fold; the N-dependence is the N=4 self-fold " +
               "(the cross-fold is a self-map on the |Δ|=1 orbit only at N=4). Mechanism = free-fermion/AT " +
               "additivity (a diagonal mode-spectator leaves the EP invariant, so (1,2)≡(2,3)≡(3,4) byte-identical; " +
               "Family B = the F89d image λ↦−λ̄−2N). The λ VALUE is the AT Theorem-2 rate of the mixed defective " +
               "eigenvector: L(q)=A+qC with A=−2·diag(n_diff) real and C anti-Hermitian, so Re λ = −2·⟨n_diff⟩_v " +
               "exactly at real loci (a q-tuned mixture of the rate-2/rate-6 eigenmodes). A closed form for the " +
               "mixture ⟨n_diff⟩(q) and general N stay open.",
               Tier.Tier1Candidate,
               "experiments/F89_MULTI_SECTOR_MONODROMY.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/MultiSectorMonodromyCensus.cs + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SectorBraidWitness.cs (inspect --root sectorbraid)")
    {
        Octic = octic ?? throw new ArgumentNullException(nameof(octic));
        CrossFold = crossFold ?? throw new ArgumentNullException(nameof(crossFold));
    }

    public override string DisplayName =>
        "Multi-sector monodromy verdict: the S₈ braid is confined to the (1,2) orbit at N=4, spreads to a 12-sector diamond at N=5";

    public override string Summary =>
        $"N-dependent: N=4 CONFINED to the D₄ orbit (dense core braid-free), N=5 SPREADS to a 12-sector diamond " +
        $"(two cross-fold families sharing a byte-identical λ, incl. the dense core (2,2)); the N-dependence is the " +
        $"N=4 self-fold, the mechanism free-fermion/AT additivity, the λ-value construction open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N=4 verdict",
                summary: "CONFINED to the D₄ orbit {(1,2),(2,1),(2,3),(3,2)}; the dense half-filled core (2,2) is " +
                         "braid-free (Door-C echo). This is the degenerate self-fold (N−2 = 2, the cross-fold maps " +
                         "(1,2) onto itself).");
            yield return new InspectableNode("N=5 verdict",
                summary: "SPREADS to the 12-sector diamond {(1,2),(1,3),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3)," +
                         "(3,4),(4,2),(4,3)}, incl. the dense cores (2,2),(3,3); two cross-fold-conjugate families of 6 " +
                         "sharing a byte-identical λ.");
            yield return new InspectableNode("membership rule",
                summary: "{|bra−ket| = 1, both popcounts ∈ [1, N−1]} ∪ its q̃↦N−q̃ cross-fold image; the N-dependence " +
                         "is the N=4 self-fold (the cross-fold is a self-map on the |Δ|=1 orbit only at N=4).");
            yield return new InspectableNode("mechanism (free-fermion / AT additivity)",
                summary: "a diagonal mode-spectator (one excitation added to both bra & ket) leaves the EP invariant " +
                         "((1,2)≡(2,3)≡(3,4) to ~1e-14); Family B is the F89d cross-fold image λ↦−λ̄−2N.");
            yield return new InspectableNode("the λ value (from below)",
                summary: "L(q)=A+qC with A=−2·diag(n_diff) real and C anti-Hermitian, so at real q Re λ = −2·⟨n_diff⟩_v " +
                         "exactly (machine zero): the AT Theorem-2 rate of the coalescing eigenvector, a q-tuned mixture " +
                         "of the rate-2/rate-6 eigenmodes (⟨n_diff⟩ = 2.31 at q=0.6209; the N=4 EP at the −4 midpoint). " +
                         "Not off the AT theorem, only off the integer-quantized lines; defectiveness = the Jordan " +
                         "coalescence, separate from the rate.");
            yield return new InspectableNode("open",
                summary: "a CLOSED FORM for the mixture ⟨n_diff⟩(q) (from the free-fermion mode geometry, without " +
                         "diagonalizing) and the codim-1-by-additivity theorem; plus general N (verified at N=4, 5; " +
                         "the octic is S₈, so no global radical form for the P₁₀ loci).");
        }
    }
}
