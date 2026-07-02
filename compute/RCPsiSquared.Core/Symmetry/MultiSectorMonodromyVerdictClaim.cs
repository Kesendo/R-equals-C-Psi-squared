using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Multi-Sector Monodromy verdict (Tier 1 candidate; computationally exact and gate-tested, the
/// free-fermion-additivity mechanism confirmed for the embedding, the λ-value grounded and its closed-form
/// mixture resolved via the mode geometry, the codim-1-by-additivity theorem LANDED in
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> and general-N membership DERIVED in the CONTAINMENT
/// direction; the exclusion half stays census-evidence):
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
/// <para><b>Mechanism: free-fermion / AT additivity, now a theorem.</b> The elementary EP is a property of the
/// |bra−ket|=1 SE-DE coherence rung; the SITE-SUMMED spectator W(ρ) = Σ_l c_l†ρc_l (one excitation added to bra
/// AND ket at the same site, summed over sites; JW strings included) is an EXACT part-by-part intertwiner of the
/// full Liouvillian (<see cref="SpectatorIntertwinerClaim"/>, Theorem B of PROOF_CODIM1_BY_ADDITIVITY), so it
/// transports the whole EP (λ, character, chain) up the ladder and (1,2) ≡ (2,3) ≡ (3,4) byte-identically
/// (verified to ~1e-14; the earlier single-MODE spectator prose is the refuted reading: c_k†ρc_k provably fails
/// the dissipator part). Family B is the <see cref="F89CrossFoldSimilarityClaim"/> image λ ↦ −λ̄ − 2N. So the
/// SAME defective eigenvalue lives in every diamond sector: shared spectral content, a symmetry broader than the
/// naive Klein-four, the sector-resolved face of free-fermion additivity.</para>
///
/// <para><b>The λ value, from below.</b> The (1,2) block is the linear pencil L(q) = A + q·C, A = −2·diag(n_diff)
/// real (the AT rates, n_diff ∈ {1,3}) and C the ANTI-HERMITIAN free-fermion hopping. Since C is anti-Hermitian, at
/// real q the coherent term feeds only Im λ, so <b>Re λ = −2·⟨n_diff⟩_v exactly</b> (machine zero) for the
/// coalescing eigenvector: the λ VALUE is the AT Theorem-2 rate of a q-tuned MIXTURE of the rate-2 (overlap) and
/// rate-6 (no-overlap) eigenmodes (⟨n_diff⟩ = 2.31 at q=0.6209, 1.90 at q=1.078; the N=4 EP sits at the −4 midpoint,
/// ⟨n_diff⟩ = 2). Not off the AT theorem, only off the integer-quantized lines; the defectiveness is the
/// eigenvector coalescence (Jordan), separate from the rate.</para>
///
/// <para><b>The closed-form mixture, from the mode geometry.</b> Rewrite the pencil as L(q) = −6·I + 4·Ô + q·Ĥ:
/// a 1-ket + 2-bra free-fermion CONTACT problem, ⟨n_diff⟩ = 3 − 2·⟨Ô⟩ with Ô the ket-bra coincidence. In the
/// mode-product basis |k⟩⟨k₁,k₂| (<see cref="RCPsiSquared.Core.F86.JordanWigner.JwBlockBasis"/>) the diagonal
/// coincidence is ⟨Ô⟩_diag = I(k,k₁)+I(k,k₂), the mode-density overlap I(a,b) = Σ_l ψ_a(l)²ψ_b(l)², which is
/// EXACTLY QUANTIZED (Dirichlet sum): I = 1/(N+1) generic, 3/(2(N+1)) for a=b or chiral a+b=N+1. So
/// ⟨Ô⟩ = Σ|c_α|²·⟨Ô⟩_diag(α) + off-diagonal δ-multiplet mixing (δ = ε_k − ε_{k₁} − ε_{k₂}). The off-diagonal
/// mixing is the closed-form boundary: at N=4 it is EXACTLY ZERO (the self-fold cancels the chiral I-pairs), so
/// ⟨n_diff⟩ = N/2 = 2 (a second proof beside the Re λ=−N axis); at N≥5 it is nonzero, and the function ⟨n_diff⟩(q)
/// along the defective branch carries a √-branch point AT each EP (smooth+algebraic, quantized ⟨n_diff⟩=3 at q=0,
/// low-q 3−c·q²) so the per-locus value is the non-radical Re-part of an S₈ octic root. The S₈ "wall" is the
/// branch-point structure of ⟨n_diff⟩(q) itself. Instruments:
/// <see cref="RCPsiSquared.Core.F86.JordanWigner.JwBlockBasis"/>,
/// <see cref="RCPsiSquared.Core.F86.JordanWigner.XyJordanWignerModes"/> (Core);
/// SectorBraidModeGeometry (Diagnostics, the decomposition + branch sweep).</para>
///
/// <para><b>The theorem landed; what stays open (why this is still Tier 1 CANDIDATE).</b> The
/// codim-1-by-additivity theorem is LANDED (<c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>: the W intertwiner
/// Theorem B + the containment orbit corollary + the two-regime Theorem A), and general-N membership is now
/// DERIVED in the CONTAINMENT direction only: the orbit of (1,2) under {one climbing W-step, transpose, Klein
/// full flip, F89d cross-fold} reproduces the census braid sets exactly (parity law: diagonal core (p,p) iff
/// |2p−N| = 1, i.e. iff N odd; sizes 4N−8 odd / 4N−12 even), so the N=6 spread is a corollary, no census needed.
/// What keeps the claim at candidate: the EXCLUSION half of membership (no braid outside the orbit) is proven
/// only at the outer edge (n_diff ≡ 1 ⟹ normal pencil) and is otherwise census-evidence (exact N=4, 5; one
/// targeted probe N=6); and the GAP byte-identity across sectors is observed, not implied by the intertwiner
/// (which transports eigenvalue and Jordan depth, not the near-defective metric geometry). The closed-form
/// mixture is resolved (the quantized-overlap contraction above); the octic is S₈ (non-solvable in radicals),
/// so no global radical form for the P₁₀ EP loci exists, and this is derived (the √-branch-point structure of
/// the function), not merely expected.</para>
///
/// <para><b>Typed parents.</b> <see cref="F89OcticMonodromyClaim"/> (the S₈ braid this census spreads across
/// sectors), <see cref="F89CrossFoldSimilarityClaim"/> (F89d, the exact antiunitary similarity λ ↦ −λ̄ − 2N
/// that pairs Family B with Family A and whose N=4 self-fold IS the confinement), and
/// <see cref="SpectatorIntertwinerClaim"/> (Theorem B, the exact W intertwiner that carries the shared λ up the
/// diamond and derives the containment half of membership). Live: <c>inspect --root sectorbraid</c>. Anchor:
/// <c>experiments/F89_MULTI_SECTOR_MONODROMY.md</c>; the exact N=5 12-set is the regression test
/// <c>Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore</c>.</para></summary>
public sealed class MultiSectorMonodromyVerdictClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89OcticMonodromyClaim Octic { get; }
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89CrossFoldSimilarityClaim CrossFold { get; }
    // Parent-edge marker for Schicht-1 wiring: Theorem B, the exact W intertwiner behind the byte-identity
    // and the containment half of general-N membership (PROOF_CODIM1_BY_ADDITIVITY).
    public SpectatorIntertwinerClaim Intertwiner { get; }

    public MultiSectorMonodromyVerdictClaim(F89OcticMonodromyClaim octic, F89CrossFoldSimilarityClaim crossFold,
        SpectatorIntertwinerClaim intertwiner)
        : base("The multi-sector monodromy verdict is N-dependent: the S₈ braid the F89 (1,2) octic carries is " +
               "CONFINED to the D₄ orbit {(1,2),(2,1),(2,3),(3,2)} at N=4 (the dense core (2,2) braid-free) but " +
               "SPREADS at N=5 to a symmetric 12-sector joint-popcount diamond (incl. the dense cores (2,2),(3,3)), " +
               "two cross-fold-conjugate families of 6 sharing a byte-identical eigenvalue λ. Membership = " +
               "{|bra−ket|=1, popcounts ∈ [1,N−1]} ∪ q̃↦N−q̃ cross-fold; the N-dependence is the N=4 self-fold " +
               "(the cross-fold is a self-map on the |Δ|=1 orbit only at N=4). Mechanism = free-fermion/AT " +
               "additivity, now a THEOREM (the site-summed spectator W=Σc_l†ρc_l is an exact intertwiner and " +
               "transports the EP, so (1,2)≡(2,3)≡(3,4) byte-identical, SpectatorIntertwinerClaim; " +
               "Family B = the F89d image λ↦−λ̄−2N). The λ VALUE is the AT Theorem-2 rate of the mixed defective " +
               "eigenvector: L(q)=A+qC with A=−2·diag(n_diff) real and C anti-Hermitian, so Re λ = −2·⟨n_diff⟩_v " +
               "exactly at real loci (a q-tuned mixture of the rate-2/rate-6 eigenmodes). The mixture is RESOLVED via " +
               "the free-fermion mode geometry: ⟨n_diff⟩=3−2⟨Ô⟩ with the mode-density overlap I(a,b) exactly " +
               "quantized (1/(N+1) or 3/(2(N+1))); N=4 pins ⟨n_diff⟩=N/2 by a vanishing off-diagonal δ-multiplet " +
               "mixing (a 2nd proof beside the Re λ=−N self-fold axis), and ⟨n_diff⟩(q) has √-branch points at the " +
               "S₈ EPs so N≥5 is non-radical (the S₈ wall IS the function's branch-point structure). The " +
               "codim-1-by-additivity theorem is LANDED (docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) and general-N " +
               "membership is DERIVED in the CONTAINMENT direction only (the W-orbit corollary: cores iff |2p−N|=1, " +
               "sizes 4N−8 odd / 4N−12 even); the EXCLUSION half stays census-evidence (edge blocks excepted: " +
               "normal pencil) and the gap byte-identity is observed, not proven.",
               Tier.Tier1Candidate,
               "experiments/F89_MULTI_SECTOR_MONODROMY.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/MultiSectorMonodromyCensus.cs + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SectorBraidWitness.cs (inspect --root sectorbraid)")
    {
        Octic = octic ?? throw new ArgumentNullException(nameof(octic));
        CrossFold = crossFold ?? throw new ArgumentNullException(nameof(crossFold));
        Intertwiner = intertwiner ?? throw new ArgumentNullException(nameof(intertwiner));
    }

    public override string DisplayName =>
        "Multi-sector monodromy verdict: the S₈ braid is confined to the (1,2) orbit at N=4, spreads to a 12-sector diamond at N=5";

    public override string Summary =>
        $"N-dependent: N=4 CONFINED to the D₄ orbit (dense core braid-free), N=5 SPREADS to a 12-sector diamond " +
        $"(two cross-fold families sharing a byte-identical λ, incl. the dense core (2,2)); the N-dependence is the " +
        $"N=4 self-fold, the mechanism the W intertwiner theorem (PROOF_CODIM1_BY_ADDITIVITY), the mixture " +
        $"⟨n_diff⟩(q) resolved via the mode geometry (quantized overlap; N=4 ⟹ N/2, N≥5 √-branch/non-radical); " +
        $"general-N membership derived in the CONTAINMENT direction (W-orbit corollary), the exclusion half " +
        $"census-evidence, the gap byte-identity observed not proven ({Tier.Label()})";

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
            yield return new InspectableNode("mechanism (free-fermion / AT additivity, now Theorem B)",
                summary: "the site-summed spectator W(ρ)=Σ_l c_l†ρc_l (one excitation added to bra AND ket at the same " +
                         "site, summed over sites, JW strings included) is an EXACT intertwiner and transports the EP " +
                         "((1,2)≡(2,3)≡(3,4) to ~1e-14; SpectatorIntertwinerClaim / PROOF_CODIM1_BY_ADDITIVITY; the " +
                         "single-MODE spectator is the refuted reading: it fails the dissipator part); Family B is the " +
                         "F89d cross-fold image λ↦−λ̄−2N.");
            yield return new InspectableNode("the λ value (from below)",
                summary: "L(q)=A+qC with A=−2·diag(n_diff) real and C anti-Hermitian, so at real q Re λ = −2·⟨n_diff⟩_v " +
                         "exactly (machine zero): the AT Theorem-2 rate of the coalescing eigenvector, a q-tuned mixture " +
                         "of the rate-2/rate-6 eigenmodes (⟨n_diff⟩ = 2.31 at q=0.6209; the N=4 EP at the −4 midpoint). " +
                         "Not off the AT theorem, only off the integer-quantized lines; defectiveness = the Jordan " +
                         "coalescence, separate from the rate.");
            yield return new InspectableNode("the closed-form mixture (mode geometry)",
                summary: "⟨n_diff⟩ = 3 − 2⟨Ô⟩ (a 1-ket+2-bra free-fermion contact problem). The mode-density overlap " +
                         "I(a,b) = Σ_l ψ_a(l)²ψ_b(l)² is exactly quantized (1/(N+1) or 3/(2(N+1))); ⟨Ô⟩ = " +
                         "Σ|c_α|²·⟨Ô⟩_diag(α) + off-diagonal δ-multiplet mixing. N=4: mixing = 0 (self-fold) ⟹ " +
                         "⟨n_diff⟩ = N/2. N≥5: mixing ≠ 0, and ⟨n_diff⟩(q) has √-branch points at the S₈ EPs ⟹ " +
                         "non-radical. The S₈ wall IS the branch-point structure of the function.");
            yield return new InspectableNode("landed + open",
                summary: "LANDED: the codim-1-by-additivity theorem (docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md: W " +
                         "intertwiner + containment orbit corollary + two-regime Theorem A); general-N membership is " +
                         "DERIVED in the CONTAINMENT direction (cores iff |2p−N|=1, sizes 4N−8 odd / 4N−12 even; the " +
                         "N=6 spread is a corollary, no census needed). OPEN: the exclusion half of membership " +
                         "(census-evidence except the normal edge), the gap byte-identity (observed, not implied by " +
                         "the intertwiner), the interior-core kernel death (CLOSED at real loci by §6's " +
                         "rate-window lemma, its strictness N-uniform via the window-edge lemma; only the complex-q loci " +
                         "remain, the arc's ONE open item). Theorem A's D-half is CLOSED at N=5 (gate TwinScalarDHalfTests): " +
                         "twin-scalar at every genuinely-complex-q diabolic (additivity's codim-1 route extends to complex q), " +
                         "the pure-imaginary-q ones semisimple by Hermiticity. The closed-form mixture is resolved (the quantized-overlap contraction).");
        }
    }
}
