# PROOF: the site-summed spectator intertwiner, and what additivity really buys the coalescences

**Status:** Theorem B (the intertwiner, Lemmas 1-3): Tier 1 derived, exact operator identities, gate-verified at machine zero (N=5, residuals 0.00e+00, commit `de4f90a`). Corollary (diamond containment, parity law): Tier 1 derived, CONTAINMENT DIRECTION ONLY; the exclusion half stays census-evidence. Theorem A (semisimplicity): the AT-locked half Tier 1 derived; the residual half a conditional twin-scalar statement, mechanism Tier 2, proven at the N=4 point. The sl(2) structure (§6): W, W†, H₀ = N̂_bra+N̂_ket−N close an sl(2) that L commutes with (machine zero N=3,4,5), which derives the climbing-rung injectivity structurally (strengthening Theorem B) and reduces the interior-core kernel death to a single spectral absence; the band deaths follow from the Edge lemma.
**Date:** 2026-07-02
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [PROOF_ABSORPTION_THEOREM.md](PROOF_ABSORPTION_THEOREM.md): the rate law Re λ = −2γ·n_diff and the light-content classification (§4.5) that makes the block pencil's real part a diagonal.
- The F89d cross-fold entry in [ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md) (~§F89d) / [`F89CrossFoldSimilarityClaim`](../../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs): the antiunitary similarity L₍₁,N−₂₎(q̄) = −P·conj(L₍₁,₂₎(q))·Pᵀ − 2N·I, machine-zero N=4..9, and its Klein legs including the unitary full flip X^⊗N.
- [PROOF_PI_FACTORS_AS_R_TIMES_D.md](PROOF_PI_FACTORS_AS_R_TIMES_D.md): the mirror group D₄ whose transpose leg this proof borrows.
- [hypotheses/DIABOLIC_BY_INTEGRABILITY.md](../../hypotheses/DIABOLIC_BY_INTEGRABILITY.md): the N=4 twin-scalar mechanism (H-scalar by Slater additivity, D-scalar at the AT midpoint) and the "what it is NOT" eliminations this proof inherits.
- [experiments/F89_PATH_K_DIABOLIC.md](../../experiments/F89_PATH_K_DIABOLIC.md): the diabolic census (11 complex-q diabolics at N=5; the N=6 flood of AT-locked crossings) and the codimension priors.
- [experiments/F89_MULTI_SECTOR_MONODROMY.md](../../experiments/F89_MULTI_SECTOR_MONODROMY.md): the sector census this proof upgrades: the 12-sector diamond sharing one defective eigenvalue byte-identically.

## Names imported from the arc (one line each, so this document parses cold)

The **census** = the numerical sector survey of [F89_MULTI_SECTOR_MONODROMY](../../experiments/F89_MULTI_SECTOR_MONODROMY.md). The **octic** = the degree-8 residual factor of the (1,2) block's characteristic polynomial after removing its rate-locked part (Galois group S₈). The **braid** = the monodromy of the octic's roots around its defective exceptional points (EPs); their loci are the simple zeros of the discriminant factor **P₁₀** (degree 10 in q²); a sector is "braid-free" when it carries no such defective coalescence. The **diamond** = the 12-sector set of blocks the census found sharing the coalescence. **AT** = the Absorption Theorem (Re λ = −2γ·n_diff, with n_diff the bra-ket Hamming distance); an **AT-locked** strand is an eigenvector family pinned to an integer rate rung, and **residual / H_B-mixed** names the complementary part of a block's spectrum, the part the hopping H_B genuinely mixes (the octic lives there); note that "residual" in the Status line and §10 separately means a numerical error norm. A **δ-multiplet** = the set of free-fermion mode triples sharing one frequency difference δ. **The 11** = the eleven complex-q diabolics recorded at N=5 (Builds-on). **R** = the site-reflection parity. **q̃** carries a tilde only to keep the ket popcount distinct from the coupling q (J = 2q, γ = 1 units); **q\*** = 0.620878 is the N=5 real defective locus the gate uses. **Klein full flip** = the unitary X^⊗N, one element of the F89d Klein four-group {I, bra flip, ket flip, full flip}. **H_B** = the block's hopping matrix (the Hamiltonian part of the pencil, the arc's symbol). **The gate** = the committed test class `SpectatorIntertwinerGateTests`, the source of every numerical input quoted here. The **coalescence gap** = |λ₁ − λ₂| of the closest eigenvalue pair at the probe q. **Tier 1 derived / Tier 2** = the repo's evidence grades (exact derivation vs supported mechanism).

## Preface

The multi-sector census found something that looked like a coincidence repeated six times: the same defective eigenvalue, to the last printed digit, living in six different coherence sectors of the same Liouvillian, with its mirror partner (λ ↦ −λ̄ − 2N, the cross-fold) in six more. Coincidences repeated six times are not coincidences; they are an operator nobody has written down yet. This note writes the operator down. It is three lines long: add one excitation to the bra and the ket at the same site, and sum over sites. Everything the census observed, the sharing, the outer boundary where it stops, the death under interaction, falls out of that one map and two four-line lemmas (the one interior-core boundary is reduced to a single spectral statement, §6).

## Abstract

Work on the operator space of the N-site XY chain (Δ=0, open ends, arbitrary bond profile J_b) under local Z-dephasing (arbitrary rates γ_j). The Liouvillian L = −i[H,·] + D preserves the joint-popcount blocks (p,q̃) (bra weight p, ket weight q̃); on each block it is the linear pencil L(q) = A + q·C with A real diagonal, A = −2·diag(n_diff) in uniform-γ=1 units (general γ_j form in §1; [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md)) and C anti-Hermitian at real q. Define, with the Jordan-Wigner site fermions c_l (strings included),

  **W(ρ) = Σ_l c_l† ρ c_l**  (the site-summed spectator: one excitation added to bra AND ket at the same site, summed over sites).

**Theorem B (the intertwiner).** W intertwines the full Liouvillian exactly, part by part: D∘W = W∘D and [H, W(·)] = W([H, ·]) as operator identities, for ANY quadratic particle-conserving H (no translation symmetry, no uniformity of J_b needed) and ANY site-dependent γ_j. W shifts blocks (p,q̃) → (p+1,q̃+1). Consequently (Lemma 3) every Jordan chain of L on block (p,q̃) whose eigenvector survives W transports to a Jordan chain of the same length at the SAME eigenvalue on block (p+1,q̃+1). On the climbing rung (1,2)→(2,3) W is injective (σ_min = √2 at N=5, gate-pinned), so the ENTIRE spectrum of L₍₁,₂₎, Jordan structure included, embeds in L₍₂,₃₎: the census's byte-identical λ within a cross-fold family IS this embedding (the other family carries the fold partner −λ̄ − 2N). The reverse map W†(σ) = Σ_l c_l σ c_l† is also an exact intertwiner, giving depth equality wherever both kernel conditions pass.

**Corollary (the diamond, containment).** The orbit of (1,2) under {one climbing W-step, the transpose, the Klein full flip, the F89d cross-fold} reproduces exactly the census's braid set: the 4-orbit at N=4, the 12-set with diagonal cores (2,2),(3,3) at N=5, the 12-set without (3,3) at N=6; in general the |p−q̃|=1 band united with its cross-fold image, size 4N−8 (N odd) / 4N−12 (N even), with a diagonal core (p,p) present iff |2p−N|=1, i.e. iff N is odd (the general-N statement is conditional on one per-N census input, the existence of the seed defective EP on block (1,2); the transport itself, W-injectivity by weight together with the transpose, full flip, and cross-fold, is derived, §5-§6, gate-checked at N=4, 5). This is the containment half: every band sector carries the shared defective λ, every fold-image sector its cross-fold partner −λ̄ − 2N (the two families of six the census recorded). The exclusion half (no braid outside) is proven only at the outer edge (§6); interior absence remains census evidence.

**Theorem A (what additivity buys the silent crossings).** Two regimes. (i) AT-locked crossings: two AT-locked eigenvector families are each analytic through a crossing and stay independent; their eigenvalue coincidences are codim-1 and automatically semisimple; this is the abundant flood (528 in one N=6 box). (ii) Residual coalescences: semisimplicity there is equivalent to the twin-scalar restriction of the pencil on the coalescing 2-plane; free-fermion additivity supplies the H-scalar half identically in q (the degenerate-multiplet descent), which is the codimension reduction; the D-scalar half is a genuine extra condition, proven at the N=4 point (the AT midpoint) and open in general. Under XXZ Δ≠0 the H-half dies (H becomes quartic): the recorded diabolic→defective flip and the death of the diamond sharing, while the D-half and the cross-fold survive.

Verification: [`SpectatorIntertwinerGateTests`](../../compute/RCPsiSquared.Diagnostics.Tests/Foundation/SpectatorIntertwinerGateTests.cs) (SLOW_MSM), all numbers in §10.

## §1 Setting and the one definition

Site basis states are bit strings a ∈ {0,1}^N; the (p,q̃) block is spanned by the coherences |a⟩⟨b| with popcount(a) = p, popcount(b) = q̃. On a block, Z-dephasing acts diagonally, D(|a⟩⟨b|) = −2(Σ_j γ_j [a_j ≠ b_j])·|a⟩⟨b| (uniform γ: −2γ·n_diff), and the XY Hamiltonian acts by hopping on each side. With J = 2q (the repo convention H = XYChain(N, 2q)) the block pencil is L(q) = A + qC.

The Jordan-Wigner site fermions are c_l† = (∏_{m<l} Z_m)·σ_l⁺ in the same site ordering the block basis uses. On a bit string, c_l†|a⟩ = s_l(a)·|a ∪ e_l⟩ for l ∉ a (zero otherwise), with the string sign s_l(a) = (−1)^{#occupied sites before l}. Define

  W(ρ) = Σ_{l=1}^{N} c_l† ρ c_l,  so on a coherence  W(|a⟩⟨b|) = Σ_{l ∉ a, l ∉ b} s_l(a)·s_l(b)·|a ∪ e_l⟩⟨b ∪ e_l|.

W adds one excitation at the SAME site to bra and ket, summed over sites; it maps block (p,q̃) into (p+1,q̃+1). Two traps recorded once: (i) the strings matter: Σ_l σ_l⁺ρσ_l⁻ is a different map and does NOT intertwine the Hamiltonian part; (ii) the site sum matters: the single-mode spectator V_k(ρ) = c_k†ρc_k (one Jordan-Wigner MODE added to both sides) intertwines the Hamiltonian part but provably NOT the dissipator (§3, remark), although the census's earlier prose described the mechanism in single-mode words.

## §2 Lemma 1: W intertwines any quadratic particle-conserving Hamiltonian

**Lemma 1.** Let H = Σ_{m,l} h_{ml} c_m†c_l with any coefficient matrix h. Then [H, W(ρ)] = W([H, ρ]) for all ρ.

*Proof.* From the canonical anticommutation relations, [H, c_a†] = Σ_m h_{ma} c_m† and [H, c_a] = −Σ_l h_{al} c_l. Expand:

  [H, Σ_a c_a†ρc_a] = Σ_a ( [H,c_a†]ρc_a + c_a†[H,ρ]c_a + c_a†ρ[H,c_a] )
                    = Σ_{m,a} h_{ma} c_m†ρc_a + W([H,ρ]) − Σ_{a,l} h_{al} c_a†ρc_l.

Renaming (a,l) → (m,a) in the third sum makes it identical to the first; they cancel. ∎

No property of h is used: the lemma covers arbitrary bond profiles, on-site potentials, and disorder. What it does NOT cover is a quartic H: an XXZ Δ-term makes [H, c_a†] more than linear in the fermions and the cancellation fails; §8 reads the consequence.

## §3 Lemma 2: W intertwines Z-dephasing, site by site

**Lemma 2.** For each site j, Z_j W(ρ) Z_j = W(Z_j ρ Z_j). Hence D∘W = W∘D for any rates γ_j.

*Proof.* Z_j = 1 − 2n_j is string-free, and Z_j c_l† Z_j = (1 − 2δ_{jl})·c_l† (the string of c_l† consists of Z's and commutes with Z_j; the local factor picks up −1 only at l = j). Therefore

  Z_j W(ρ) Z_j = Σ_l (Z_j c_l† Z_j)(Z_j ρ Z_j)(Z_j c_l Z_j) = Σ_l (1−2δ_{jl})² c_l† (Z_jρZ_j) c_l = W(Z_jρZ_j),

the sign squaring away because the SAME l sits on both sides of ρ. Summing γ_j(Z_j·Z_j − id) over j gives D∘W = W∘D. ∎

*Remark (why the single mode fails).* For V_k the two sides carry c_k† and c_k built from DIFFERENT sites after the mode expansion: the cross terms l ≠ m acquire (1−2δ_{jl})(1−2δ_{jm}), which does not square away, and each surviving term hits a distinct output coherence |a∪e_l⟩⟨b∪e_m|, so nothing can cancel. Gate table §10: the V_k dissipator residual is O(1) for every k while its Hamiltonian residual is machine zero. That is the refutation of the single-mode reading, kept as data.

## §4 Lemma 3: Jordan transport, with the sharp premise

**Lemma 3.** Let L₂W = WL₁ and let x₁, …, x_m be a Jordan chain of L₁ at λ ((L₁−λ)x_j = x_{j−1}, x₀ = 0; x₁ the eigenvector). If **Wx₁ ≠ 0**, then Wx₁, …, Wx_m is a Jordan chain of L₂ at λ; in particular L₂ has λ in its spectrum with Jordan depth ≥ m.

*Proof.* (L₂−λ)Wx_j = W(L₁−λ)x_j = Wx_{j−1}. Kernel hits are downward-closed along a chain (Wx_j = 0 forces Wx_{j−1} = (L₂−λ)Wx_j = 0), so Wx₁ ≠ 0 forces every Wx_j ≠ 0, and (L₂−λ)^{j−1}Wx_j = Wx₁ ≠ 0 makes the images linearly independent. ∎

The premise is sharp in both directions: Wx₁ ≠ 0 is necessary (a dead eigenvector transports nothing), and the tempting weakening "Wx_m ≠ 0 suffices" is FALSE: kernel hits truncate chains from the bottom (take L₁ the 2×2 nilpotent Jordan block, W = (0 1), L₂ = (0): WL₁ = L₂W, Wx₂ ≠ 0, yet L₂ is 1×1 with no chain).

## §5 Theorem B and the climbing rung

**Theorem B.** W is an exact intertwiner of the full pencil, L₍p+1,q̃+1₎(q)·W = W·L₍p,q̃₎(q) for all q (Lemmas 1+2 applied to the two parts; the relation holds part-by-part since L is linear in q). On the rung (1,2)→(2,3) at N=5, W is injective with σ_min(W) = √2 (gate-pinned; at N=4 the gate reads σ_min = 1). The FACT of injectivity is not a numerical input: §6 derives it from the sl(2) weight, the rung (1,2) sitting below the half-filling anti-diagonal (p+q̃ = 3 < N for N ≥ 4) where the raising operator W has no highest-weight vectors and is injective; only the norm value σ_min = √2 is measured. So by Lemma 3 the entire spectrum of L₍₁,₂₎, eigenvalues, multiplicities, and Jordan blocks, embeds in L₍₂,₃₎ at the same q. The reverse spectator W†(σ) = Σ_l c_l σ c_l† satisfies the same two lemmas (identical proofs with c_l† ↔ c_l), intertwining (p,q̃) → (p−1,q̃−1); where both an upward chain and its downward return avoid the kernels, the Jordan depth is EQUAL, not merely bounded.

This upgrades the census's central observation from numerics to theorem: the defective EP of the (1,2) octic, its eigenvalue, √-branch character, and chain, is not "repeated" in (2,3); it is the SAME chain, pushed through W. The byte-identical λ across a cross-fold family of six (recorded to 1.2e-14 in the census; a non-committed review-round cross-check saw the whole-spectrum match at 1.3e-11, and the committed gate pins the load-bearing residuals at machine zero; the other family of six carries the fold partner −λ̄ − 2N) is the visible face of the embedding.

**Scope.** The two lemmas never use uniformity: site-dependent γ_j and arbitrary J_b (disorder included) are covered. The one thing used is Δ=0 (H quadratic). What is NOT claimed: equality of the coalescence GAP across sectors (observed byte-identical in the census; the intertwiner gives eigenvalue and depth transport, not metric equality of the near-defective geometry).

## §6 Where the sharing stops: the sl(2) behind the kernel, and the normal edge

The intertwining identity is exact on every block, so the sharing can only stop where Lemma 3's premise fails: the transported eigenvector must die in ker W. The census (braid-free (4,4) at N=5) is therefore consistent only with Wx₁ = 0 for the (3,3) defective eigenvector; the gate measures exactly that, with a sharpening: the entire near-defective 2-plane dies, ‖Wx₁‖ = 1.7e-15 AND ‖Wx₂‖ = 2.5e-15 ((3,3)→(4,4) at the real locus), while on the interior rung both transport at norm √2 = σ_min. The diamond boundary is a kernel phenomenon, not a failure of the identity. What follows shows the kernel is not accidental: W is the raising operator of an sl(2), and every death is a highest-weight annihilation.

**The raising operator closes an sl(2).** In the doubled (bra, ket) picture W is one creation on each side at the same site, summed, and W† = Σ_l c_l(·)c_l† is its Hilbert-Schmidt adjoint (one annihilation on each side). Their commutator collapses: the cross terms (l ≠ m) cancel because c_l c_m† = −c_m† c_l flips the sign on each side, and the diagonal terms (l = m) sum to

  **[W, W†] = N̂_bra + N̂_ket − N =: H₀,**  with [H₀, W] = 2W and [H₀, W†] = −2W†,

where N̂_bra, N̂_ket are the bra- and ket-occupation numbers. So {W, W†, H₀} is an sl(2) triple: W raises, W† lowers, H₀ is the Cartan, and block (p,q̃) carries weight m = p + q̃ − N (from −N at the vacuum coherence to +N at the full one). Because W and W† are both exact intertwiners (Theorem B and its reverse, §5) while H₀ is block-diagonal, L commutes with the whole sl(2): it is an sl(2)-module endomorphism of the doubled operator space. The four identities hold at machine zero for N = 3, 4, 5, and the closure is N-independent.

**Injectivity is by weight, and Theorem B strengthens.** A raising operator annihilates exactly the highest-weight vectors, which live at weight ≥ 0. On any block with m < 0, that is p + q̃ < N (strictly below the half-filling anti-diagonal), W has no highest-weight vectors and is injective. The climbing rung (1,2)→(2,3) sits at m = −2 (at N=5; m = 3−N in general, negative for N ≥ 4), so its injectivity, taken in §5 as the one numerical input, is structural; more generally W is injective on every below-anti-diagonal rung, for every N. The whole-spectrum embedding of Theorem B thus needs no per-N measurement of injectivity (only the norm value σ_min = √2 is measured, not the fact of it). The Lefschetz count is exact per block: dim ker W|₍p,q̃₎ = 0 for m < 0, and dim(p,q̃) − dim(p+1,q̃+1) for m ≥ 0 (verified at every N=5 block; 100 − 25 = 75 at (3,3)).

**Kernel death is highest-weight annihilation.** The generalized λ-eigenspace M_λ is an sl(2)-submodule (W, W†, H₀ map it to itself, since L commutes with them), and on it L = λ + Nilp with Nilp an sl(2)-endomorphism. By Schur, M_λ = ⊕_j U_j ⊗ V_j (V_j the spin-j irrep, U_j its multiplicity space), and Nilp acts as ⊕_j nilp_j ⊗ I. A defective Jordan chain is a nilpotent block inside one U_j; through the tensor with V_j it appears at every weight of V_j, across the whole block range p + q̃ − N ∈ {−2j, …, +2j} at fixed d = p − q̃. Its top block carries the highest weight, where W kills it. So the death is set by the spin of the chain's component, not by the locus.

**The bands: capped by the normal edge.** For the |p−q̃| = 1 band, whose defective chain the §7 orbit reaches from the seed EP on (1,2) by transpose and full flip, the chain's next block up is a normal edge, (3,4)→(4,5) going up and (1,2)→(0,1) going down. Were W(x₁) ≠ 0 there, Lemma 3 would place a Jordan chain on the edge block; the Edge lemma below forbids it; so W(x₁) = 0. In sl(2) language the band chain is spin 1 (weights −2, 0, +2), topped at (3,4) and bottomed at (1,2), and W annihilates its two extreme weights. The outer-edge death is exactly this cap, now read as a weight statement.

At the outer edge the death is structural:

**Edge lemma.** On the blocks (0,1), (1,0), (N−1,N), (N,N−1), every coherence has n_diff = 1, so A = −2γ·I is scalar and L(q) = −2γ·I + q·C is a normal pencil at real q (C anti-Hermitian). A normal operator has no Jordan blocks: no defective EP can live on an edge block; whatever arrives there must arrive semisimple or die. (Gate: A₍₄,₅₎ = −2I to 0.00e+00.)

**The interior core: the one open piece, made sharp.** The diagonal cores (2,2),(3,3) are the d = 0 part of the diamond, at weights −1 and +1. Their kernel death (‖Wx₁‖ = 1.7e-15 at (3,3)→(4,4), and its full-flip image (2,2)→(1,1)) says the d = 0 λ-isotypic component is spin ½: top (3,3), bottom (2,2), reaching neither (4,4) at weight +3 nor (1,1) at weight −3. What is left to derive is precisely that this component is spin ½ and not spin 3/2, equivalently that the diamond eigenvalue λ carries no DEFECTIVE Jordan block in the interior blocks (4,4) and (1,1). Two natural shortcuts provably fail, which is why this piece is separate from the bands: (i) the Edge lemma does not reach it, since (4,4) has n_diff ∈ {0,2} (A = diag(0,−4), not scalar), so it is not normal and a Jordan block there is not forbidden a priori; (ii) the second, d-raising ladder V = Σ_l c_l†(·)c_l† does close its own sl(2) (Cartan N̂_bra − N̂_ket) but does NOT commute with L: it intertwines the dephasing ([D,V] = 0, the same string-sign squaring as W) yet not the Hamiltonian ([H, V(·)] ≠ V([H, ·]), the mirror of the single-mode V_k of §3, which failed on D and not H), so [L,V] ≠ 0 and the cross-fold cannot carry the band's cap onto the core. The interior death stays measured; the sl(2) has reduced it from "why does a 100-dimensional eigenvector fall into a 75-dimensional kernel" to the single spectral statement "λ has no Jordan block in the single-hole-on-both-sides cone (4,4) ≅ (1,1)." A targeted argument on that cone, or the multiplet descent of remainder 2 (§8), is the route.

## §7 The corollary: the diamond as an orbit, containment direction

Combine four exact maps, all Tier 1:

1. **the climbing W-step** (1,2)→(2,3): injective (σ_min = √2), transports everything (§5);
2. **the transpose** (p,q̃)↔(q̃,p): the mirror group's D leg ([PROOF_PI_FACTORS_AS_R_TIMES_D](PROOF_PI_FACTORS_AS_R_TIMES_D.md)), LINEAR (ρ↦ρᵀ, no conjugation). Its precise spectral action, derivable in two lines from Hᵀ = H and the dephasing's symmetry n_diff(a,b) = n_diff(b,a): T·L₍p,q̃₎(q) = L₍q̃,p₎(−q)·T, and since A is real and C purely imaginary entrywise, conj(L(q)) = L(−q̄); together spec₍q̃,p₎(q) = conj(spec₍p,q̃₎(q̄)). At a REAL locus with a REAL coalescing eigenvalue, and the gate reads Im λ_A ≈ 1e-16 there, the transposed sector therefore carries the SAME λ;
3. **the Klein full flip** X^⊗N: (p,q̃)→(N−p,N−q̃) at the SAME q, unitary ([F89d Klein section](../ANALYTICAL_FORMULAS.md); the F89d Klein four-group is {I, bra flip P, ket flip Q, full flip QP = X^⊗N});
4. **the F89d cross-fold** q̃↦N−q̃ at q̄: antiunitary, Jordan-preserving, an exact identity verified in exact arithmetic per N = 4..9, Tier 1 derived ([`F89CrossFoldSimilarityClaim`](../../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs)).

All "sharing at one q" summaries in this document implicitly use that the census loci are REAL (q̄ = q) and the coalescing eigenvalue is real there; at complex q the transpose and fold partners carry conjugated values at the conjugate parameter point, and the two-families picture must be read with that fine print. The fine print is itself census-confirmed: at N=4 the locus eigenvalue is COMPLEX (λ = −4 + 1.239i) and the census recorded the 4-orbit as conjugate-related, NOT byte-identical, exactly as the conjugation rule demands; at N=5 the locus eigenvalue is real and the families are byte-identical. At even N the four fold-fixed/fold-swapped band members around half filling are simultaneously band and fold image, so those sectors carry BOTH λ and −λ̄ − 2N (two coalescences in one block); the clean six-plus-six split is the odd-N (N=5) reading.

Starting from the (1,2) defective EP, these reach: the |p−q̃|=1 band {(p,p±1): both ∈ [1,N−1]} (one W-step + transpose + full flip suffice; the flip returns the upper rungs without ever climbing past half-filling, where W is rank-deficient by dimension count alone: e.g. (2,3)→(3,4) at N=5 maps dimension 100 into 50), united with the band's cross-fold image. Walked explicitly: N=4 gives the confined 4-orbit {(1,2),(2,1),(2,3),(3,2)} (the fold maps the band into itself); N=5 gives the recorded 12-set including the diagonal cores (2,2),(3,3); N=6 gives the 12-set {(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(4,5),(5,4),(1,4),(4,1),(2,5),(5,2)}, with (3,3) absent.

**Parity law.** A diagonal core (p,p) is reached iff it is the fold image of a band member: N−(p±1) = p ⟺ |2p−N| = 1 ⟺ N odd (the exactly-half-filled (N/2,N/2) is never reached; that is why even N misses its core). **Size.** |band ∪ fold(band)| = 4(N−2) − |overlap|, and the overlap is 4 at even N (the four members around half filling are fold-fixed or fold-swapped within the band) and 0 at odd N: **4N−12 (even), 4N−8 (odd)**. N=5 and N=6 both give 12, but they are different sets.

**Scoping, stated once and bluntly.** This corollary proves CONTAINMENT: every band sector in the orbit carries the shared defective eigenvalue and every fold-image sector its partner −λ̄ − 2N (each move transports the chain; the W-step's injectivity is derived by weight in §6, so the only per-N input is the existence of the seed defective EP on (1,2)). It does NOT prove exclusion: that no sector outside the orbit carries the braid is proven only at the outer edge (§6) and is otherwise the census's empirical result (exact at N=4, 5; one targeted probe at N=6). Any claim promotion cites the containment half only. The N=6 "no census needed" statement applies to the spread: the 12-set's membership follows from this corollary with zero new compute.

## §8 Theorem A: the two regimes of silence

The word "diabolic" in this arc covers two different objects, and the codimension story differs between them. This section separates them; the separation resolves a tension between two source documents (§9).

**Regime 1: AT-locked crossings, abundant and automatic.** An AT-locked strand is an exact eigenvector family of the full pencil, analytic in q, whose rate is pinned to an integer rung (Re λ = −2γ·n). Two such families cross where their frequencies meet: ONE analytic condition, codim-1-complex, dense curves of them (528 counted in one complex-q scan window at N=6, [F89_PATH_K_DIABOLIC](../../experiments/F89_PATH_K_DIABOLIC.md)). The crossing is semisimple for the classical block-diagonal reason: both eigenvectors exist analytically THROUGH the crossing and remain independent; nothing couples them at the crossing point because they are exact eigenvectors on both sides. For THIS regime the "two crossing modes come from different families, so their eigenvectors stay independent" voice is exactly right.

**Regime 2: residual coalescences, the 11, and the twin-scalar condition.** Inside the residual (H_B-mixed) part of a block, an eigenvalue coincidence is generically DEFECTIVE (one complex parameter, non-normal pencil: the branches meet in a √-branch point; codim-1 and loud). A SEMISIMPLE coincidence there needs the 2×2 restriction of the pencil on the coalescing plane to be scalar: two extra complex conditions, codim-3-complex in total, generically unreachable by one knob q. What free-fermion additivity buys is exactly ONE of the two extra conditions, identically in q:

**Twin-scalar lemma (conditional).** If the two coalescing directions descend from a single degenerate free-fermion multiplet (Slater additivity: E = Σε over each side's occupied modes, so equal mode multisets up to degeneracy give equal E), the H-part restriction is a multiple of the identity for ALL q: the H-half of the scalar condition costs nothing. The remaining D-half (the dephasing restriction scalar on the plane) is ONE genuine extra condition; where it is met, the coincidence is diabolic.

At the N=4 point both halves are established from below ([DIABOLIC_BY_INTEGRABILITY](../../hypotheses/DIABOLIC_BY_INTEGRABILITY.md): H|₂D = 2iJ·I by the 4-fold multiplet, D|₂D = −4γ·I at the AT midpoint, the rate midway between the 2γ and 6γ rungs, Re λ = −4γ, which is also the palindrome mirror line), and the discriminant records the collapse as the even-multiplicity factor (3q⁴+q²−1)²: coincidence and semisimplicity merging into one condition is the algebraic meaning of a squared factor. Whether the D-half is automatic (e.g. by chiral pairing within the δ-multiplet) at the 11 recorded complex-q diabolics of N=5 is NOT settled by the existing record (the census logs q, λ, loop character, exponent; not the multiplet descent); until that one check is run, the honest general statement is: **additivity collapses the semisimple coincidence from codim-3 to codim-≤2, and to codim-1 exactly where the D-half is supplied** (at N=4: by the AT midpoint on the mirror line).

**Δ-tightness (both regimes).** XXZ Δ≠0 makes H quartic: Lemma 1's cancellation dies (gate: W's H-residual jumps to 2.4e-1 at Δ=0.3 while the D-residual stays 0.00e+00), the multiplet descent breaks (E ≠ Σε), and the recorded phenomenology follows at once: diabolics defect instantly, the diamond sharing dies, while the F89d cross-fold (bit-flip parity, integrability-independent) survives. The theorem's name is honest: everything here is BY additivity, and only the transpose/flip/fold legs outlive its removal.

## §9 Bookkeeping: defective vs diabolic, and the two voices

- **The defective braid EPs (the P₁₀ simple zeros) need NO codimension collapse.** Codim-1 is the GENERIC count for a defective coalescence of a non-normal one-parameter family. What Theorem B explains about them is not their existence but their SHARING across sectors (the intertwiner) and where the sharing stops (kernel death, the normal edge).
- **The diabolics (the squared factor, and the complex-q census) are the codim story.** Regime 1 needs only analyticity of the AT-locked families; Regime 2 needs the twin-scalar lemma, of which additivity supplies half.
- **The two voices reconciled.** [F89_PATH_K_DIABOLIC](../../experiments/F89_PATH_K_DIABOLIC.md)'s synthesis ("the Hamiltonian and the diagonal dissipator do not mix the occupation-number sectors: two crossing modes come from different sectors") is the Regime-1 truth, over-applied in that sentence to Regime 2, where [DIABOLIC_BY_INTEGRABILITY](../../hypotheses/DIABOLIC_BY_INTEGRABILITY.md)'s "Not a commuting-symmetry separation" (both modes R=+1, same sector, twin-scalar instead) is the careful voice. Within one block no fine decomposition is invariant under the FULL L: the site sectors are dephasing-rigid but hopping-mixed, the mode sectors hopping-rigid but dephasing-mixed. That non-commuting pair is the whole game; where one of them goes scalar on a plane, the game goes silent.

## §10 Verification

[`SpectatorIntertwinerGateTests`](../../compute/RCPsiSquared.Diagnostics.Tests/Foundation/SpectatorIntertwinerGateTests.cs) (`SLOW_MSM`; run `dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "FullyQualifiedName~SpectatorIntertwinerGateTests"`), N=5, γ=1, H = XYChain(N, 2q), residuals normalized ‖X₂W−WX₁‖_F/(‖X₁‖₂·‖W‖_F), X ∈ {A, C} the pencil parts (the L₁, L₂ of §4 split part-by-part); λ_A = the (1,2) defective eigenvalue at q\*, λ_B its cross-fold partner:

| Measurement | Value | Reading |
|---|---|---|
| W (1,2)→(2,3): A-part / C-part residual | 0.00e+00 / 0.00e+00 (both q) | Lemmas 1+2, machine-exact |
| V_k, k=1..5: A-part residual | 3.8-4.1e-01 | single-mode dissipator failure (kept as refutation data) |
| V_k, k=1..5: C-part residual | ≤ 2e-16 | single-mode Hamiltonian half is exact |
| σ_min(W) on (1,2)→(2,3) | √2 to 1e-14 | injectivity of the climbing rung |
| interior transport ‖Wx₁‖, ‖Wx₂‖ at q*=0.620878 | 1.414, 1.414 | the chain crosses the rung at full norm |
| boundary (3,3)→(4,4): ‖Wx₁‖, ‖Wx₂‖ | 1.7e-15, 2.5e-15 | kernel death of the whole 2-plane (§6) |
| (3,3) pair at λ_B = −conj(λ_A)−2N | gap 4.771e-04 = (1,2) gap; \|mean−pred\| 5.5e-15 | cross-fold + embedding, byte-level |
| edge: A₍₄,₅₎ vs −2I; rank W (3,4)→(4,5) | 0.00e+00; 5/5 | edge lemma (§6) |
| Δ=0.3 (XXZ): W H-part / D-part residual | 2.4e-01 / 0.00e+00 | Δ-tightness (§8) |
| [W,W†] − (N̂_bra+N̂_ket−N), every block | 0.00e+00 (N=3,4,5) | sl(2) closure: W is the raising operator, Cartan H₀ (§6) |
| ker W per block vs Lefschetz dim(p,q̃)−dim(p+1,q̃+1) | exact, every block (75 at (3,3); 0 below the anti-diagonal) | injectivity by weight; kernel death = highest-weight (§6) |
| [L,W†] (2,3)→(1,2): A-part / C-part residual | 0.00e+00 / 0.00e+00 | L commutes with the whole sl(2), reverse spectator too |

The single-mode refutation and the W identity were additionally derived and verified from scratch (independent numpy build) by two independent reviewers on 2026-07-02; that review round is recorded in the local planning layer (deliberately not part of the repository), and its load-bearing numbers are reproduced by the committed gate above. The §6 sl(2) closure, its L-commutation, and the Lefschetz kernel dimensions were likewise independently re-derived the same day (convention-free, with a four-way stacking/sign probe and a disorder-general H), and are pinned by item 6 of the committed gate.

## What remains open

1. The interior-core kernel death ((3,3)→(4,4) and its full-flip image (2,2)→(1,1)) is reduced, not yet closed. §6 identifies W as an sl(2) raising operator (with W†, H₀ = N̂_bra+N̂_ket−N, all commuting with L), so kernel death is highest-weight annihilation: the band cases now follow from the Edge lemma, and the climbing-rung injectivity is derived. What remains is the interior core alone, reduced to the single spectral statement "the diamond eigenvalue λ carries no defective Jordan block in the single-hole-both-sides cone (4,4) ≅ (1,1)." The Edge lemma cannot reach it ((4,4) is not normal) and the second d-ladder does not commute with L, so the residue is genuinely separate; a targeted argument on the cone, or the multiplet descent of remainder 2, is the route.
2. Theorem A's D-half at the 11 complex-q diabolics of N=5 (one targeted eigenvector-descent check; until then codim-≤2 is the proven general statement).
3. Gap byte-identity across sectors: observed, not implied by the intertwiner.
4. The EXCLUSION half of membership beyond the outer edge (§7's scoping): census evidence by design; an intertwiner can never prove an absence.
