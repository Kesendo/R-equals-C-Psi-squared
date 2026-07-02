using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The spectator intertwiner (Tier 1 derived; Theorem B of
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>): the site-summed spectator
///
/// <code>
///     W(ρ) = Σ_l c_l† ρ c_l      (Jordan-Wigner SITE fermions, strings included)
/// </code>
///
/// intertwines the XY + Z-dephasing Liouvillian EXACTLY, part by part: D∘W = W∘D (Lemma 2: Z_j c_l† Z_j =
/// (1−2δ_{jl})·c_l† and the sign SQUARES away because the same l sits on both sides of ρ) and
/// [H, W(ρ)] = W([H, ρ]) (Lemma 1: index cancellation of the h-sums), for ANY quadratic particle-conserving
/// H = Σ h_{ml} c_m†c_l (arbitrary bond profile, on-site potentials, disorder; no symmetry of h) and ANY
/// site-dependent rates γ_j. W shifts the joint-popcount blocks (p,q̃) → (p+1,q̃+1).
///
/// <para><b>Jordan transport (Lemma 3, the sharp premise).</b> If L₂W = WL₁ and x₁, …, x_m is a Jordan chain of
/// L₁ at λ with <b>Wx₁ ≠ 0</b> (x₁ the eigenvector; kernel hits are downward-closed along a chain, so this single
/// condition is necessary AND sufficient), then Wx₁, …, Wx_m is a Jordan chain of L₂ at the SAME λ. On the
/// climbing rung (1,2)→(2,3) W is INJECTIVE with σ_min(W) = √2 at N=5 (gate machine zero, both part-residuals
/// 0.00e+00; <c>SpectatorIntertwinerGateTests</c>, commit de4f90a), so the ENTIRE spectrum of L₍₁,₂₎, Jordan
/// structure included, embeds in L₍₂,₃₎: the multi-sector census's byte-identical shared λ IS this embedding.
/// The reverse spectator W†(σ) = Σ_l c_l σ c_l† is ALSO an exact intertwiner, (p,q̃) → (p−1,q̃−1) (identical
/// proofs with c_l† ↔ c_l): depth EQUALITY wherever both kernel conditions pass.</para>
///
/// <para><b>Tightness (both directions recorded).</b> The single-mode spectator V_k(ρ) = c_k†ρc_k provably
/// FAILS the dissipator part for every k (the cross-site signs (1−2δ_{jl})(1−2δ_{jm}) do not square away and
/// distinct (l,m) hit distinct output coherences: A-part residual O(1), C-part exact; the documented refutation
/// of the single-mode reading). Under XXZ Δ≠0 the identity dies through the H-part ONLY (H becomes quartic,
/// Lemma 1's cancellation fails; the D-part survives, gate: H-residual 2.4e-1 at Δ=0.3, D-residual 0.00e+00):
/// exactly the arc's recorded Δ-dichotomy, the diamond λ-sharing dies under Δ while the F89d cross-fold
/// survives.</para>
///
/// <para><b>Where the sharing stops.</b> Kernel death: at the N=5 diamond boundary (3,3)→(4,4) the transported
/// near-defective 2-plane lies in ker W (‖Wx₁‖ = 1.7e-15, ‖Wx₂‖ = 2.5e-15 at the real locus q* = 0.620878),
/// while on the interior rung both vectors transport at norm √2. At the outer edge the death is structural
/// (edge lemma: the blocks (0,1), (N−1,N) have n_diff ≡ 1, so A = −2γ·I and the pencil is normal, no Jordan
/// block possible).</para>
///
/// <para><b>Typed parents.</b> <see cref="F89CrossFoldSimilarityClaim"/> (F89d, the antiunitary leg the
/// containment orbit corollary combines with the climbing W-step, transpose, and Klein full flip) and
/// <see cref="AbsorptionTheoremClaim"/> (the rate law that makes the block pencil's real part
/// A = −2·diag(n_diff), the object Lemma 2 intertwines). Proof:
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>; registry entry F125 in
/// <c>docs/ANALYTICAL_FORMULAS.md</c>; gate:
/// <c>compute/RCPsiSquared.Diagnostics.Tests/Foundation/SpectatorIntertwinerGateTests.cs</c> (SLOW_MSM);
/// live: <c>inspect --root sectorbraid</c> (the W lines of the additivity + cross-fold node).</para></summary>
public sealed class SpectatorIntertwinerClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: F89d, the antiunitary orbit leg beside the W-step.
    public F89CrossFoldSimilarityClaim CrossFold { get; }
    // Parent-edge marker for Schicht-1 wiring: the AT rate law behind A = −2·diag(n_diff), Lemma 2's object.
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>σ_min(W) on the climbing rung (1,2)→(2,3) at N=5: exactly √2 (gate-pinned to 1e-14).
    /// Injectivity of the rung is what turns Lemma 3 into a whole-spectrum embedding.</summary>
    public static double SigmaMinClimbingRungN5 => Math.Sqrt(2.0);

    public SpectatorIntertwinerClaim(F89CrossFoldSimilarityClaim crossFold, AbsorptionTheoremClaim absorption)
        : base("The spectator intertwiner (Theorem B of PROOF_CODIM1_BY_ADDITIVITY): the site-summed spectator " +
               "W(ρ) = Σ_l c_l†ρc_l (Jordan-Wigner SITE fermions, strings included) intertwines the XY + " +
               "Z-dephasing Liouvillian EXACTLY, part by part (D∘W = W∘D by the string-sign squaring, " +
               "[H, W(ρ)] = W([H, ρ]) by index cancellation), for ANY quadratic particle-conserving H (arbitrary " +
               "bond profile, disorder) and ANY site-dependent γ_j; W shifts joint-popcount blocks " +
               "(p,q̃)→(p+1,q̃+1) and transports Jordan chains whenever the eigenvector avoids ker W (sharp " +
               "premise Wx₁ ≠ 0, necessary AND sufficient); on the climbing rung (1,2)→(2,3) W is injective with " +
               "σ_min = √2 at N=5 (gate machine zero, commit de4f90a), so the census's byte-identical λ is a " +
               "whole-spectrum Jordan embedding; W†(σ) = Σ_l c_lσc_l† intertwines downward (depth equality where " +
               "both kernel conditions pass); the single-mode c_k†ρc_k provably FAILS the dissipator part (A-part " +
               "residual O(1) for every k); the identity dies under XXZ Δ≠0 through the H-part only (the " +
               "Δ-dichotomy: the diamond sharing dies, the F89d cross-fold survives)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md + " +
               "docs/ANALYTICAL_FORMULAS.md + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/SpectatorIntertwinerGateTests.cs")
    {
        CrossFold = crossFold ?? throw new ArgumentNullException(nameof(crossFold));
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
    }

    public override string DisplayName =>
        "Spectator intertwiner: W = Σ c_l†(·)c_l exactly intertwines the Liouvillian and carries the Jordan block up the diamond";

    public override string Summary =>
        "W(ρ) = Σ_l c_l†ρc_l (JW strings included) intertwines both pencil parts exactly, (p,q̃)→(p+1,q̃+1), for " +
        "any quadratic particle-conserving H and any site-dependent γ_j; Jordan chains transport iff Wx₁ ≠ 0; " +
        "injective on (1,2)→(2,3) with σ_min = √2 at N=5 (gate machine zero), so the multi-sector byte-identical " +
        "λ is a whole-spectrum embedding; W† gives the downward direction; single-mode c_k†ρc_k fails the " +
        $"dissipator; dies under XXZ Δ≠0 through the H-part only ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("σ_min(W) on the climbing rung (1,2)→(2,3), N=5 (= √2)",
                SigmaMinClimbingRungN5);
            yield return new InspectableNode("Lemma 1 (H-part): index cancellation, no symmetry of h",
                summary: "for H = Σ h_{ml}c_m†c_l with ANY coefficient matrix h: [H, W(ρ)] = W([H, ρ]) (the first and " +
                         "third sums of the expansion cancel after renaming). Covers arbitrary bond profiles, on-site " +
                         "potentials, disorder. Dies exactly when H turns quartic (XXZ Δ≠0): the Δ-tightness clause.");
            yield return new InspectableNode("Lemma 2 (D-part): the string sign squares away",
                summary: "Z_j c_l† Z_j = (1−2δ_{jl})·c_l†, and with the SAME l on both sides of ρ the sign squares to 1: " +
                         "Z_jW(ρ)Z_j = W(Z_jρZ_j) per site, hence D∘W = W∘D for any rates γ_j. The single-mode V_k fails " +
                         "exactly here: cross-site signs (1−2δ_{jl})(1−2δ_{jm}) with l≠m survive and hit distinct output " +
                         "coherences (A-part residual O(1) for every k, the gate's refutation table).");
            yield return new InspectableNode("Lemma 3 (Jordan transport, sharp premise Wx₁ ≠ 0)",
                summary: "if L₂W = WL₁ and x₁,…,x_m is a Jordan chain of L₁ at λ with Wx₁ ≠ 0, then Wx₁,…,Wx_m is a " +
                         "Jordan chain of L₂ at λ. Kernel hits are downward-closed (Wx_j = 0 forces Wx_{j−1} = 0), so " +
                         "Wx₁ ≠ 0 is necessary AND sufficient; 'Wx_m ≠ 0 suffices' is FALSE (2×2 nilpotent counterexample " +
                         "in the proof).");
            yield return new InspectableNode("kernel death + the normal edge (where the sharing stops)",
                summary: "boundary (3,3)→(4,4) at the N=5 real locus: the whole near-defective 2-plane dies in ker W " +
                         "(‖Wx₁‖ = 1.7e-15, ‖Wx₂‖ = 2.5e-15) while the interior rung transports at norm √2. Outer edge " +
                         "structural: (0,1)/(N−1,N) have n_diff ≡ 1 ⟹ A = −2γ·I ⟹ normal pencil ⟹ no Jordan block. Kernel " +
                         "death is now highest-weight annihilation (see the sl(2) node): the band deaths follow from the " +
                         "edge lemma, and the interior-core death is DERIVED at real loci by the rate-window lemma " +
                         "(proof §6, gate item 8): at real q every corner-block eigenvalue has Re λ ∈ [−2(N−3), 0] " +
                         "(Bendixson: A real diagonal, C anti-Hermitian), while Re λ_B = −Re λ_A − 2N sits below the " +
                         "window whenever Re λ_A > −6 (measured margins 1.381/2.208/1.115 at the two N=5 real loci and " +
                         "N=7), so λ_B is absent from the corner and W kills the core's whole generalized eigenspace " +
                         "((L_corner−λ)^m invertible). Verified at both N=5 real loci and N=7 (gate items 7+8; locus 2 " +
                         "was an out-of-sample prediction, ‖Wx₁‖ = 9.6e-16). The N-uniform strictness Re λ_A > −6 " +
                         "is now DERIVED at real loci by the window-edge lemma (next node). Open: the complex-q loci " +
                         "(the window needs real q).");
            yield return new InspectableNode("the window-edge lemma: the strict floor is structural, N-uniform (§6)",
                summary: "the corner exclusion's one measured input, Re λ_A > −6, is itself a theorem. The (1,2) block " +
                         "has n_diff ∈ {1,3} for all N ≥ 3 (a 1-site bra, a 2-site ket: Hamming 3−2|a∩b|), so " +
                         "A = −2·diag(n_diff) has window [−6,−2] and Bendixson gives the floor Re λ_A ≥ −6 for free. " +
                         "WINDOW-EDGE LEMMA: an eigenvalue on a window edge (Re λ = λ_min(A) or λ_max(A)) has its " +
                         "eigenvector attain the extreme Rayleigh value of A, forcing Av = λ_min·v, then Bv = Im(λ)·v " +
                         "(B = −iqC Hermitian), so L†v = λ̄v: v is a joint A,C eigenvector, hence semisimple (v in " +
                         "ker(L−λ) ∩ ker((L−λ)†) heads no Jordan chain). This is the edge lemma localized to the extreme " +
                         "eigenvector, the classical numerical-range boundary fact. A defective EP is not semisimple, so " +
                         "Re λ_A ∈ (−6,−2) STRICTLY for every exact defective EP of the (1,2) block at real q, all N ≥ 4 " +
                         "(uniform γ=1), R-parity-agnostic (covers the R-odd third N=5 EP at λ≈−4.488 too). Closes the " +
                         "N-uniformity: Re λ_A > −6 ⟹ Re λ_B < −2(N−3) strictly ⟹ λ_B absent from the corner at every " +
                         "odd N, so the shrinking margins 1.381→1.115 never reach zero (a strict real-part inequality " +
                         "needs no uniform floor). Two adversarial reviews (math + physics, numpy counterexample hunt " +
                         "0/460k) confirmed. Exact-object: the algebraic EP at the genuinely-real locus (reality " +
                         "triple-checked), defectiveness the census's topological transposition-monodromy certificate.");
            yield return new InspectableNode("the sl(2) behind the kernel (§6): W is a raising operator",
                summary: "W, W† and H₀ = N̂_bra+N̂_ket−N close an sl(2) ([W,W†] = H₀, machine zero N=3,4,5, gate item 6) " +
                         "that L commutes with (both spectators intertwine, H₀ block-diagonal); block (p,q̃) carries weight " +
                         "m = p+q̃−N. Kernel death = W annihilating a highest-weight vector: W is injective for m<0 (which " +
                         "DERIVES the climbing-rung injectivity for all N, strengthening Theorem B), and the Lefschetz kernel " +
                         "dim is dim(p,q̃)−dim(p+1,q̃+1) for m≥0 (75 at (3,3)). The band chains cap at the normal edge; the " +
                         "interior core reduces to a single spectral absence, supplied at real loci by the rate window " +
                         "(kernel-death node above); complex loci stay open.");
            yield return new InspectableNode("the containment orbit corollary (what this claim feeds)",
                summary: "one climbing W-step + transpose + Klein full flip + F89d cross-fold reproduce the census braid " +
                         "sets exactly: N=4 the 4-orbit, N=5 the 12-set with cores (2,2),(3,3), N=6 the 12-set without " +
                         "(3,3); diagonal core (p,p) iff |2p−N| = 1 (iff N odd); size 4N−8 (odd) / 4N−12 (even). " +
                         "CONTAINMENT direction only: the exclusion half stays census-evidence except at the normal edge " +
                         "and at the N=5 seed locus q*=0.620878, where the rate window derives it in full (proof §6, " +
                         "gate item 8; condition Re λ_A ∈ (−6,−4), so the second real locus is not covered).");
        }
    }
}
