# The R₉₀ frozen divisor: watching profiles that pin an eigenvalue for every coupling

**Status:** Theorem (lower bound) Tier 1 derived; cofactor closed form + tightness criterion Tier 1 derived (Section 6); uniform-endpoint constants in closed form and tightness-for-generic-J at every N Tier 1 derived (Section 7, the two boundary clocks); the J-valuation lower bound, total and per pair, Tier 1 derived (Section 8, the distance ladder)
**Date:** 2026-07-25
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the anti-palindromic watching locus (every reflection pair of dephasing rates sums to the same value), the single-excitation corner block of the Liouvillian carries the eigenvalue λ = −4γ̄ with multiplicity at least ⌊N/2⌋, for every Hamiltonian coupling J (equality generic, verified N = 3..6). The mechanism is a rank bottleneck of a cell-level mirror, not an invariant subspace and not a spectral symmetry. The nonvanishing cofactor is a single N(N−1)/2 determinant in closed form, (−1)^N(4γ̄)^⌈N/2⌉·det((X P_{O₊} X)|_{V₋}), whose nonvanishing is exactly tightness (Section 6).
**Verification:** [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) (must print "R90 frozen divisor gate: ALL GREEN", ~7-10 min)
**Depends on:** [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md) (the R₉₀ reshuffle and its fixed locus), [GAMMA_FOLD_PAIR_OF_MIRRORS](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) (the X^N cross-dock used in the corollary), [PROOF_F139_SEAM_IDENTITY](PROOF_F139_SEAM_IDENTITY.md) (the sibling on the character side)

---

## What this means

Tune the watching so that mirrored sites share their rates evenly around the mean: the first and the last site together, the second and the second-to-last together, every reflection pair carrying the same total. On that one locus, and only there, the chain acquires eigenvalues that the Hamiltonian cannot move. However hard the sites talk to each other, at least ⌊N/2⌋ decay modes, generically exactly that many, sit frozen at one exact rate, one mode per balanced pair. They are not protected by a conserved quantity, and no subspace of states carries them; take any J and the eigenvectors have rearranged, yet the eigenvalue has not. What holds them is a bookkeeping fact: a mirror on the cells has more rooms it leaves alone than rooms it flips, and the surplus has nowhere to go.

There is a second reading, about how strongly they are held. Switch the coupling off entirely and twice as many modes sit at that rate; switch it on and half of them leave. How fast each one leaves is set by nothing but a distance: the mirrored pair whose sites are d apart needs the coupling raised to the power 2d before its mode stirs, because the chain has to walk the excitation from the one site over to the other first. The outermost pair, the two ends of the chain, is the slowest to let go. Section 8 proves the lower half of that law by counting rooms and steps.

F139 taught this arc that a wall can be a divisor instead of a symmetry. This is the same lesson on the home γ axis: a factor of the characteristic polynomial that divides out exactly on a locus, with no symmetry of the spectrum behind it.

## 1. Setting and definitions

Take the open chain of N qubits with an excitation-conserving Hamiltonian H(J) = J·H₁ whose single-excitation matrix h (the N×N matrix ⟨e_a|H₁|e_b⟩, where e_a is the computational state with the single excitation at site a) is real, symmetric, and invariant under the site reversal R: a ↦ N+1−a. The isotropic Heisenberg chain (h = hopping 2 on neighbours plus the ZZ diagonal) and the XY chain (hopping only) both qualify. Site-resolved Z-dephasing acts with rates γ_l (real; positivity is nowhere used in the proof, and the corollary in Section 4 exploits that), σ = Σ_l γ_l, γ̄ = σ/N.

The Liouvillian L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ) preserves the joint-popcount blocks (popcount of the bra index, popcount of the ket index). The **corner block** is the (1,1) block, spanned by the cells |e_a⟩⟨e_b|, written v_{(a,b)}. On this block

  L_block(J) = J·K − 2Γ,  K v_{(a,b)} = −i Σ_c (h_{ac} v_{(c,b)} − h_{cb} v_{(a,c)}),  Γ v_{(a,b)} = (γ_a + γ_b)·v_{(a,b)} for a ≠ b, 0 for a = b.

Two structures on the cells:

- **The locus.** The R₉₀ reshuffle of [F91](PROOF_F91_GAMMA_NINETY_DEGREES.md) acts on the rate profile; its fixed-point set is the anti-palindromic class **γ_l + γ_{R(l)} = 2γ̄ for every l** (at odd N this forces the middle rate to the mean). All statements below hold on this locus.
- **The mirror.** τQ is the linear involution of cells **(a,b) ↦ (Rb, Ra)** (transpose composed with reversal on both sides). It splits the diagonal cells D = {(a,a)} and the off-diagonal cells O = {(a,b), a≠b} into ±1 eigenspaces D±, O±. Throughout, P_S denotes the orthogonal projector onto the span of the cell set or subspace S (P_D, P_{D₊}, P_{D₋}, P_{O₊}).

Three dimension counts, by inspection: the τQ-fixed cells in O are the anti-diagonal (a, Ra) with a ≠ Ra, so **dim O₊ − dim O₋ = 2⌊N/2⌋**; the diagonal cells pair (a,a) ↔ (Ra,Ra), so **dim D₋ = ⌊N/2⌋**; and Γ vanishes on D.

## 2. The two oddness lemmas

**Lemma 1 (K is τQ-odd): τQ K τQ = −K.**

*Proof.* (τQ K τQ v)_{(a,b)} = (K τQ v)_{(Rb,Ra)} = −i Σ_c [h_{Rb,c}(τQv)_{(c,Ra)} − h_{c,Ra}(τQv)_{(Rb,c)}] = −i Σ_c [h_{Rb,c} v_{(a,Rc)} − h_{c,Ra} v_{(Rc,b)}]. Substituting c → Rc and using h_{Rx,Ry} = h_{xy} (R-invariance) and h symmetric turns this into +i Σ_c [h_{ac} v_{(c,b)} − h_{cb} v_{(a,c)}] = −(Kv)_{(a,b)}. ∎

Consequence: K maps O₊ into D₋ ⊕ O₋. The D₊ component of Kv is zero automatically for v ∈ O₊.

**Lemma 2 (the recentered rate operator is τQ-odd on O, exactly on the locus): τQ (2Γ − 4γ̄)|_O τQ = −(2Γ − 4γ̄)|_O.**

*Proof.* On the cell (a,b), a ≠ b, the conjugated operator reads 2(γ_{Rb} + γ_{Ra}) − 4γ̄, and the locus gives γ_{Ra} + γ_{Rb} = 4γ̄ − γ_a − γ_b, so the entry is −(2(γ_a+γ_b) − 4γ̄). ∎

On the diagonal cells the recentered rate is the constant −4γ̄ on both sides, which is even, not odd; this is why the full-block identity τQ(L_block + 4γ̄)τQ = −(L_block + 4γ̄) + 8γ̄·P_D carries a defect exactly on D, and why the argument below works on O and treats D as a constraint. The defect is not small: the block spectrum is not palindromic about −4γ̄, and no multiset symmetry argument applies.

## 3. The theorem

**Theorem (frozen divisor).** On the R₉₀-fixed locus, for every J,

  **dim ker(L_block(J) + 4γ̄) ≥ ⌊N/2⌋.**

Equivalently: (λ + 4γ̄) raised to the ⌊N/2⌋-th power divides det(λ − L_block(J)) as a polynomial identity in J on the locus. Generically the bound is an equality (verified N = 3..6).

*Proof (the pencil argument).* Let W := {v ∈ O₊ : P_{D₋} K v = 0}, a subspace of O₊ cut by dim D₋ = ⌊N/2⌋ linear conditions, so dim W ≥ dim O₊ − ⌊N/2⌋. Take v ∈ W:

1. Kv has no D-part: the D₊ part vanishes by Lemma 1 (oddness), the D₋ part by the definition of W. Hence Kv ∈ O₋.
2. (2Γ − 4γ̄)v ∈ O₋ by Lemma 2 (v is supported on O and τQ-even; an odd diagonal operator maps O₊ to O₋).
3. Therefore (L_block(J) + 4γ̄)v = J·Kv − (2Γ − 4γ̄)v lies in O₋, for every J.

So Φ_J: W → O₋, v ↦ (L_block(J) + 4γ̄)v is a linear map into a space of dimension dim O₋ = dim O₊ − 2⌊N/2⌋, and

  dim ker Φ_J ≥ dim W − dim O₋ ≥ (dim O₊ − ⌊N/2⌋) − (dim O₊ − 2⌊N/2⌋) = ⌊N/2⌋.

Every kernel vector is an exact eigenvector of L_block(J) at −4γ̄. ∎

The multiplicity is a dimension bottleneck of the mirror: the fixed cells of τQ (the anti-diagonal, one per balanced pair plus its transpose) give O₊ a surplus of 2⌊N/2⌋ rooms over O₋, the D₋ constraint taxes away half, and the remainder must freeze. Note what the proof does not use: no invariant subspace (the frozen eigenvectors move with J, only the eigenvalue stands still; Section 9's first bullet), no spectral palindromy (Section 2), no diagonalizability assumption.

Two structural by-products, both verified at machine precision in the gate: the frozen eigenvectors carry **zero weight on the diagonal cells** and have **τQ-even O-part**; the J-dependence of the frozen modes is only the ⌊N/2⌋-dimensional kernel line of Φ_J rotating inside the fixed subspace W.

## 4. Corollary: the antidiagonal corners, by the gamma fold

The block (1, N−1) (single excitation against N−1 excitations) is the image of the corner block under the one-sided X^N bridge, and the [gamma fold](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) turns that bridge into algebra: conjugating by right-multiplication with X^N sends L(γ⃗) to L(−γ⃗) − 2σ. The theorem applied at the profile −γ⃗ (the locus is preserved; the root becomes −4·(−γ̄) = +4γ̄) then gives, for every J,

  **the (1, N−1) block carries the eigenvalue 4γ̄ − 2σ with multiplicity ⌊N/2⌋.**

The two roots −4γ̄ and 4γ̄ − 2σ sum to −2σ: they are partners under the pair of mirrors, and recentered at the palindrome center (x = λ + σ) they sit at x = ±(σ − 4γ̄). At N = 4 the two coincide at the center x = 0. What once looked like two separate root families, including the N = 3 block (1,2) sighting at −2γ̄, is one theorem and its fold image: at N = 3, (1,2) is the antidiagonal corner, and 4γ̄ − 2σ = −2γ̄ there.

## 5. Why only the corners

The proof needs one affine recentering that makes the rate operator odd (Lemma 2). Under reversal a rate −2γ_S (S the disagreement set of the cell) goes to −2(2γ̄|S| − γ_S): the center depends on the size |S|. The corner block has a single off-diagonal size class (|S| = 2), so one center serves; the block (2,2) has classes |S| ∈ {0, 2, 4} and admits no single center, and indeed carries no eigenvalue at the corner's root −4γ̄ (gate control, N = 4, 5, 6). Away from the R₉₀ locus the whole structure disappears: partial balance (all but one condition satisfied) yields nothing, and the N = 3 closed form shows the defect as an explicit linear factor (M₍₁,₂₎ below is the (1,2) block of L, the N = 3 antidiagonal corner of Section 4),

  det(−2γ₂·I − M₍₁,₂₎) = 512·J⁴·(γ₁+γ₃)²·(γ₁+γ₃−2γ₂)·(4J² + (γ₁+γ₃)(γ₁+γ₃−2γ₂)),

so the balanced root exists exactly on the locus (given J ≠ 0), and the distance of the nearest eigenvalue grows linearly in the defect.

## 6. The cofactor: closed form, tightness, semisimplicity

The pencil argument of Section 3 bounds the multiplicity from below; this section computes what is left of the determinant after the frozen factor divides out, and the answer closes the tightness question along the way.

Write M̃ := L_block(J) + 4γ̄ for the recentered block, and split off its τQ-even part: M̃ = X + 4γ̄·P_D (this X is an operator on the corner block, not the many-body Pauli X^N of Section 4), where

  X = J·K − 2Δ,  Δ v_{(a,b)} = (δ_a + δ_b)·v_{(a,b)} on O, 0 on D,  δ_l := γ_l − γ̄.

X is τQ-odd (Lemmas 1 and 2; on the locus δ_{R(l)} = −δ_l), and X is **γ̄-free**: the mean rate enters M̃ only through the even defect 4γ̄·P_D. Let V₊ = D₊ ⊕ O₊ and V₋ = D₋ ⊕ O₋ be the parity eigenspaces of τQ, and note dim V₋ = ⌊N/2⌋ + dim O₋ = N(N−1)/2.

**Theorem (cofactor).** On the R₉₀-fixed locus, the characteristic polynomial p(ε) = det(εI − M̃) factors as ε^⌊N/2⌋·q(ε) with

  **q(0) = (−1)^N · (4γ̄)^⌈N/2⌉ · det( (X P_{O₊} X)|_{V₋} ),**

a determinant of size N(N−1)/2 whose entries are free of γ̄. Consequently the frozen multiplicity is exactly ⌊N/2⌋ if and only if this determinant is nonzero, and in that case the frozen eigenvalue is semisimple (algebraic = geometric = ⌊N/2⌋).

*Proof.* An odd operator exchanges the parity blocks, and an even one preserves them, so in the split V₊ ⊕ V₋ the matrix M̃ − εI has diagonal blocks Λ₊ = diag(4γ̄−ε on D₊, −ε on O₊) and Λ₋ = diag(4γ̄−ε on D₋, −ε on O₋), and off-diagonal blocks X₊₋, X₋₊ (the diagonal rate part −2Δ of X is odd, hence off-diagonal here). For ε ∉ {0, 4γ̄} the Schur complement gives

  det(M̃ − εI) = det(Λ₊) · det( Λ₋ − X Λ₊⁻¹ X |_{V₋} ),  Λ₊⁻¹ = (4γ̄−ε)⁻¹ P_{D₊} − ε⁻¹ P_{O₊}.

Pull the pole at ε = 0 out of the second factor: with A(ε) := Λ₋ − (4γ̄−ε)⁻¹·X P_{D₊} X |_{V₋} (regular at ε = 0),

  det(M̃ − εI) = (4γ̄−ε)^⌈N/2⌉ · (−1)^{dim O₊} · ε^{dim O₊ − dim V₋} · det( X P_{O₊} X |_{V₋} + ε·A(ε) ).

The exponent is dim O₊ − dim V₋ = ⌊N/2⌋ (Section 1's counts), det(XP_{O₊}X + εA(ε)) is regular at ε = 0, and both sides are rational functions equal off finitely many points, hence equal as rational functions; the left side is a polynomial, so the identity extends to ε = 0. This re-proves the divisor bound (order ≥ ⌊N/2⌋ at ε = 0, an independent second proof of the theorem of Section 3), and reading off the coefficient of ε^⌊N/2⌋, with p(ε) = (−1)^{N²}·det(M̃ − εI) (the block is N²×N²), gives q(0) = (−1)^{N² + dim O₊}(4γ̄)^⌈N/2⌉ det(XP_{O₊}X|_{V₋}). The sign: dim O₊ = N(N−1)/2 + ⌊N/2⌋, and N² + N(N−1)/2 + ⌊N/2⌋ ≡ N (mod 2) for every N (check both parities of N). Tightness: q(0) ≠ 0 says the algebraic multiplicity is exactly ⌊N/2⌋; the pencil gives geometric ≥ ⌊N/2⌋, and algebraic ≥ geometric always, so all three coincide. ∎

Three consequences, all pinned in the gate:

- **The γ̄-stratification is exact.** The cofactor is (4γ̄)^⌈N/2⌉ times a polynomial in (J, δ⃗) only: the whole γ̄-dependence of the residual spectrum at the frozen root is the even defect's ⌈N/2⌉ diagonal cells of D₊.
- **Tightness for generic J.** For a fixed locus profile, det(XP_{O₊}X|_{V₋}) is a polynomial in J of degree N(N−1) with leading coefficient det((K P_{O₊} K)|_{V₋}), which is nonzero for every N: Section 7 proves this by computing it in closed form, and the gate scans N = 3..10 on both the Heisenberg and the XY single-excitation matrix as a check. Hence the multiplicity is exactly ⌊N/2⌋ for all but finitely many J, at every N and every locus profile.
- **Small N in closed form** (symbolic, in the antisymmetric coordinates δ₁ = γ₁ − γ̄, δ₂ = γ₂ − γ̄):

    N = 3: q(0) = 2¹²·γ̄²·J⁴·(3J² − δ₁²)
    N = 4: q(0) = 2²⁰·γ̄²·J⁸·(8J⁴ − 4J²(3δ₁² + 2δ₁δ₂ + δ₂²) + (δ₁² − δ₂²)²)
    N = 5: q(0) = −2³⁰·γ̄³·J¹²·Q₅(J², δ₁, δ₂) with Q₅ of total degree 8 in (J, δ₁, δ₂), leading term 25J⁸

  (The XY chain differs only in the polynomial coefficients: its J-pure terms are 2J² and 5J⁴ at N = 3, 4.) The J-power in front is **2⌊N/2⌋⌈N/2⌉ = 2⌊N²/4⌋** (4, 8, 12, 18, 24 at N = 3..7), not the 4(N−2) the first three N suggest: the two expressions coincide exactly through N = 5 and split at N = 6 (exact-rational J → 0 valuations 18 and 24 at N = 6, 7; the N = 6 discriminator is gate-pinned). Mode-resolved, the valuation ladder is a distance ladder: the balanced pair (c, Rc) contributes one branch of the D₋ Schur complement with valuation J^{2(N+1−2c)}, twice the pair's site distance, because the Hamiltonian must walk the excitation from c to Rc before that frozen mode moves off its perch. [Section 8](#8-the-valuation-law-the-walk-to-the-anti-diagonal) proves the ladder as a lower bound, total and per pair; the first rung is also available directly, since the order-J² Schur complement on D₋ is −C_anti†C_anti with C_anti = P_anti K P_{D₋} (P_anti the projector onto the anti-diagonal cells) and rank C_anti = 1 at even N, 0 at odd N: only the distance-1 middle pair of even N is one hop from the anti-diagonal. From the theorem of this section only the total degree N(N−1) and the leading coefficient follow.

## 7. The uniform endpoint: the two boundary clocks

At the fully degenerate point of the locus, δ⃗ = 0 (uniform watching), X = J·K and the cofactor collapses to a J-monomial:

  q(0)|_{δ⃗=0} = (−1)^N·(4γ̄)^⌈N/2⌉·J^{N(N−1)}·D_N,  D_N := det( (K P_{O₊} K)|_{V₋} ).

This section computes D_N in closed form, for every N. Three lemmas reduce it to a Gram determinant, and the Gram determinant is a pure power of the chain's clock modulus.

**Lemma 3 (the pair basis diagonalizes K²).** Let λ_1 < ... < λ_N and u_1, ..., u_N be the eigenvalues and orthonormal eigenvectors of h, with R-parities π_i (R u_i = π_i u_i; the spectrum is nondegenerate, see Lemma 5). Writing w_{ij} := u_i ⊗ u_j for the cell-space product vectors, τQ w_{ij} = π_iπ_j·w_{ji}, so the vectors φ_{ij} := (w_{ij} − π_iπ_j w_{ji})/√2 over i < j form an orthonormal basis of V₋ (count: N(N−1)/2 ✓), and K² φ_{ij} = −(λ_i−λ_j)²·φ_{ij}.

**Lemma 4 (the Gram reduction).** K = −i·(a real symmetric matrix) is anti-Hermitian, K† = −K, which is what turns −K P_{D₊} K into +C†C below. K is τQ-odd, so P₊ K|_{V₋} = K|_{V₋} (P₊ the projector onto V₊), hence K P_{O₊} K = K² − K P_{D₊} K = K² + C†C on V₋, with C := P_{D₊} K|_{V₋}. By the Weinstein-Aronszajn identity,

  D_N = (−1)^{N(N−1)/2} · Π_{i<j} (λ_i−λ_j)² · det( I − W ),  W := C Λ⁻² C†,  Λ² := diag((λ_i−λ_j)²).

Moreover I − W = V†V with V_{i,k} := ⟨w_{ii}, d_k⟩ (d_k the D₊ basis): K φ_{ij} = −i(λ_i−λ_j)·ψ_{ij} with ψ_{ij} := (w_{ij} + π_iπ_j w_{ji})/√2, so the (λ_i−λ_j)² of Λ⁻² cancels and W_{kl} = Σ_{i<j} ⟨d_k, ψ_{ij}⟩⟨ψ_{ij}, d_l⟩, and {w_{ii}} ∪ {ψ_{ij}} is an orthonormal basis of V₊ containing d_k, so W_{kl} = δ_{kl} − Σ_i ⟨d_k, w_{ii}⟩⟨w_{ii}, d_l⟩. Since the w_{ii} are τQ-even (orthogonal to D₋), the N×N Gram matrix G_{ij} := ⟨w_{ii}, P_D w_{jj}⟩ = Σ_a u_i(a)²u_j(a)² equals VV†, and therefore

  det(I − W) = det(V†V) = pdet(G)  (the product of the nonzero eigenvalues of G, provided rank G = dim D₊ = ⌈N/2⌉, which Lemma 5 gives).

**Lemma 5 (the two boundary clocks).** The single-excitation eigenbasis of the open chain is a cosine angle lattice, with modulus M depending on the chain (not the running clock of [ClockHandLadder](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs), whose two hands tick on ω and γ; here "clock" is the timeless lattice of angles the excitation lives on):

- **Heisenberg (M = N).** λ_k = 4cos(kπ/N) + N − 5 and u_k(a) ∝ cos((2a−1)kπ/(2N)), k = 0..N−1 (the DCT-II basis). The check is two lines: the interior rows force the dispersion, the end rows hold identically for the half-integer cosine (the chain's +2 boundary defect from the ZZ diagonal is exactly what the product-to-sum identity 2cosθ·cos(θ/2) = cos(3θ/2) + cos(θ/2) absorbs), and matching the right end quantizes sin(Nθ) = 0.
- **XY (M = N+1).** λ_k = 4cos(kπ/(N+1)) and u_k(a) ∝ sin(akπ/(N+1)), k = 1..N (the DST-I basis, classical).

Both spectra are nondegenerate (distinct cosines), so Lemma 3 applies. For B_{a,k} := u_k(a)² everything reduces to geometric sums: in the Heisenberg case Σ_{k=0}^{N−1} cos(jkπ/N) = 1 for odd j and N·[N | j] for even j, and the products 2cos·cos split into sums with j = 2(a−b) and j = 2(a+b−1), whose resonances are exactly a = b and b = Ra; the XY case runs the same sums at modulus N+1. The result is one law for both chains:

  **B Bᵀ = (1 − 1/M)·𝟙𝟙ᵀ/N + (I + R)/(2M).**

B is doubly stochastic (B𝟙 = Bᵀ𝟙 = 𝟙), the anti-symmetric space is killed, and on the R-symmetric space (I+R)/2 acts as the identity: the spectrum of BBᵀ, hence of G = BᵀB, is exactly

  {1 (the flat vector), 1/M with multiplicity ⌈N/2⌉ − 1, 0 with multiplicity ⌊N/2⌋},

so rank G = ⌈N/2⌉ and **pdet(G) = M^{−⌊(N−1)/2⌋}**. ∎

**Theorem (uniform-endpoint constants).** For every N ≥ 3, on the Heisenberg (M = N) and XY (M = N+1) open chains,

  **D_N = (−1)^{N(N−1)/2} · Π_{i<j}(λ_i−λ_j)² · M^{−⌊(N−1)/2⌋} ≠ 0.**

Two consequences:

- **Tightness is now a theorem at every N** (upgrading the N = 3..10 scan of Section 6): D_N ≠ 0 makes det((XP_{O₊}X)|_{V₋}) a nonzero polynomial in J for every locus profile (leading coefficient J^{N(N−1)}·D_N), so the frozen multiplicity is exactly ⌊N/2⌋ for all but finitely many J; at the uniform point itself the determinant is the monomial J^{N(N−1)}·D_N, so there the multiplicity is exactly ⌊N/2⌋ for **every** J ≠ 0.
- **The J-pure constants of Section 6 are clock numbers.** Π(λ_i−λ_j)² is the discriminant of the (integer) characteristic polynomial of h, and the correction is a pure power of the clock modulus. N = 3 Heisenberg: disc = 2304, M^{−1} = 1/3, D₃ = −768 = −2⁸·3, reproducing the 3J² of Section 6; the gate pins the exact rational assembly at N = 4, 5 (at the gate's point γ̄ = 9/100, J = 4/3, the N = 5 uniform cofactor is exactly −2⁶⁴/2989355625).

The two moduli are the two boundary clocks of the open chain: the XY excitation lives on the Dirichlet lattice sin(akπ/(N+1)), which is the committed SE cosine lattice of modulus N+1 ([NivenRationalityRoot](../../compute/RCPsiSquared.Core/Symmetry/NivenRationalityRootClaim.cs) and the [F65 registry entry](../ANALYTICAL_FORMULAS.md) live there, as do the F129 collision combs and the F139 wall at modulus 11); the Heisenberg excitation lives on the Neumann half-integer lattice cos((2a−1)kπ/(2N)) of modulus N, which had no committed anchor before this proof. The frozen divisor's residual constant reads off which boundary clock the chain carries.

## 8. The valuation law: the walk to the anti-diagonal

Section 6 read the J → 0 valuation of the cofactor determinant off the data and called the per-pair ladder an observation. This section derives it. The single geometric input is that h is tridiagonal (nearest-neighbour hopping plus a diagonal; both chains of Section 1 qualify), and the quantity the answer is written in is the site distance of the balanced pair,

  d_c := N + 1 − 2c,  c = 1, ..., ⌊N/2⌋,  with Σ_c d_c = ⌊N/2⌋·⌈N/2⌉ = ⌊N²/4⌋.

**Theorem (valuation law).** On the R₉₀-fixed locus, for every profile and every N,

  **ord_J det( (X P_{O₊} X)|_{V₋} ) ≥ 2·Σ_c d_c = 2⌊N²/4⌋,**

and if the profile is generic (δ_a + δ_b ≠ 0 for every off-diagonal cell off the anti-diagonal), the ⌊N/2⌋ × ⌊N/2⌋ Schur complement S(J) of that matrix on D₋ satisfies the sharper per-pair form

  **ord_J S_{c,c'} ≥ d_c + d_{c'},  that is  S(J) = Λ·S̃(J)·Λ with Λ := diag(J^{d_c}) and S̃ regular at J = 0.**

Both hold with equality wherever they have been measured: the total at N = 3..10 for both chains, the per-pair orders at N = 3..7, all integer-exact (Section 10).

The proof is a Gram form, a grading, and a counting argument.

**Lemma 6 (the cofactor determinant is a Gram determinant).** In the cell basis h real symmetric makes K complex symmetric (K^T = K), Δ is diagonal, hence X^T = X, and the parity projectors are real. Writing Y := P_{O₊}X|_{V₋} as a matrix from an orthonormal basis of V₋ to one of O₊, of size dim O₊ × dim V₋ with dim O₊ − dim V₋ = ⌊N/2⌋,

  (X P_{O₊} X)|_{V₋} = Y^T Y,  so by Cauchy-Binet  det( (X P_{O₊} X)|_{V₋} ) = Σ_I det(Y_I)²,

the sum running over the maximal square minors of Y (row sets I of size dim V₋). Any other basis of O₊ serves as well: an unnormalized row basis turns the identity into Y^T D Y with D the diagonal of inverse squared norms, and Cauchy-Binet then carries those positive weights along, which is the form Section 11 uses. In particular ord_J of the left side is at least 2·min_I ord_J det Y_I. Cancellation between minors can only raise the left side, so it cannot damage the bound. This is not a positivity statement: Y is complex and the squares are squares of complex numbers. The Gram shape is used only to halve the bookkeeping, so that a count made once on Y pays twice.

**Lemma 7 (the level grading).** For a cell (a,b) put s(a,b) := a + b − (N+1) and

  **ℓ(a,b) := |s(a,b)|,**

the distance of the cell to the anti-diagonal {b = R(a)}. Then ℓ ∘ τQ = ℓ, so every parity basis vector is supported on a single level and ℓ is defined on the basis vectors of D±, O±, V±; and the entries of Y obey

  Y_{y,x} = 0 whenever |ℓ(y) − ℓ(x)| ≥ 2,  ord_J Y_{y,x} ≥ 1 whenever ℓ(y) ≠ ℓ(x).

*Proof.* τQ sends (a,b) to (R(b), R(a)) and therefore s to −s, which gives the invariance of ℓ. For the entries: Δ is diagonal in the cell basis, so it preserves every cell and hence s; the J-linear part of X is J·K, which couples (a,b) to (c,b) and to (a,c) with amplitudes h_{ac} and h_{cb}, zero unless the site indices are neighbours or equal, so one power of J moves s by at most one. Since X is affine in J, cells whose s differ by two or more are uncoupled at every order, and cells with s ≠ s' are uncoupled at order J⁰. Passing from s to ℓ costs nothing: if s and s' have opposite signs then |s − s'| = |s| + |s'|, so |s − s'| ≤ 1 already forces |ℓ − ℓ'| ≤ 1. ∎

The Heisenberg ZZ diagonal of h contributes same-cell terms only, so it is level-diagonal and transports nothing. That is why the two chains, whose h differ exactly by that diagonal, carry the same valuation.

**Lemma 8 (the level census).** Write n_j := dim(V₋ ∩ level j) and M_j := dim(O₊ ∩ level j). Then

  n_0 = 0,  M_0 = 2⌊N/2⌋,  and for j ≥ 1:  n_j = N − j,  M_j = N − j − ε_j,

where ε_j := 1 if j = d_c for some balanced pair c, and 0 otherwise.

*Proof.* Level 0 is exactly the anti-diagonal {(a, R(a))}, N cells, each of them τQ-fixed and hence in V₊: n_0 = 0, and of the N cells the odd-N centre cell is diagonal, leaving M_0 = 2⌊N/2⌋ in O₊. For j ≥ 1 the level consists of the N − j cells with s = +j and the N − j cells with s = −j; τQ maps the first set bijectively onto the second with no fixed cell, so V₊ and V₋ each receive N − j dimensions. A diagonal cell (a,a) has s = 2a − N − 1, so it sits at level j ≥ 1 exactly when j = |2a − N − 1|, that is j = d_c with a = c or a = R(c); those two cells are exchanged by τQ and contribute one dimension to D₊ and one to D₋. Subtracting the D₊ dimension from V₊ ∩ level j leaves M_j. ∎

The census is the whole mechanism in two lines. At each of the ⌊N/2⌋ distances d_c the columns outnumber the rows by exactly one, because the balanced pair puts a diagonal cell there. All the spare rows sit at level 0, on the anti-diagonal, where there are 2⌊N/2⌋ of them and no columns at all. These are the same τQ-fixed rooms whose surplus froze the eigenvalue in Section 3, seen from the other side: there they made the kernel, here they are the only place the excess can go.

**Lemma 9 (the transport bound).** Let A be a set of columns and I a set of rows with |I| = |A|, write Y_{I,A} for the square submatrix they cut out, and put

  F_j(I, A) := #{columns of A at level ≥ j} − #{rows of I at level ≥ j}.

Then ord_J det Y_{I,A} ≥ Σ_{j≥1} max(0, F_j(I, A)), and this pays as follows:

- for A the full column set, ord_J det Y_{I,A} ≥ Σ_{j≥1} Σ_{i≥j} ε_i = Σ_i i·ε_i = Σ_c d_c = ⌊N²/4⌋;
- for A the O₋ columns together with the single D₋ column of one pair c, ord_J det Y_{I,A} ≥ d_c.

*Proof.* Expand det Y_{I,A} over the bijections π from A onto I (σ is the total rate throughout this document). By Lemma 7 a term is nonzero only if |ℓ(π(x)) − ℓ(x)| ≤ 1 for every column x, and its J-order is then at least the number of columns with ℓ(π(x)) ≠ ℓ(x). Fix j ≥ 1 and let E_j := {x ∈ A : ℓ(x) ≥ j, ℓ(π(x)) ≤ j−1}. Every column of level ≥ j outside E_j goes to a row of level ≥ j, so #{columns of A at level ≥ j} ≤ #{rows of I at level ≥ j} + |E_j|, that is |E_j| ≥ F_j(I, A). The one-step constraint forces ℓ(x) = j exactly for x ∈ E_j, so the sets E_1, E_2, ... are pairwise disjoint, and each of their elements is one mismatch: the order of the term is at least Σ_j |E_j| ≥ Σ_j max(0, F_j(I, A)).

For the two specializations, bound the rows of I at level ≥ j by all the rows there, Σ_{i≥j} M_i, and use Lemma 8. With A the full column set the level counts are the n_i, so F_j ≥ Σ_{i≥j}(n_i − M_i) = Σ_{i≥j} ε_i ≥ 0, and summing over j counts each diagonal level i once for every 1 ≤ j ≤ i. With A the O₋ columns plus the D₋ column of the pair c, every diagonal column except that one is dropped, so the count at level i ≥ 1 falls to n_i − ε_i + [i = d_c] and F_j ≥ #{i ≥ j : i = d_c}, which sums over j to d_c. ∎

*Proof of the theorem.* The first display is Lemmas 6 and 9. For the second, let A be the O₋ columns together with the D₋ column of the pair c, and B the same with c'. In its rectangular form Cauchy-Binet reads det((Y^TY)_{A,B}) = Σ_I det(Y_{I,A})·det(Y_{I,B}), and Lemma 9's second specialization bounds the two factors by d_c and d_{c'} for every I. Hence the bordered determinant of T := (X P_{O₊} X)|_{V₋} has ord_J ≥ d_c + d_{c'}. For a generic profile the unbordered block T_{O₋,O₋} is a unit at J = 0: its value there is the Gram matrix of the vectors Δw over the O₋ basis, diagonal with entries proportional to (δ_a + δ_b)², nonzero exactly under the stated genericity. So the Schur complement entries S_{c,c'} = det T_{A,B} / det T_{O₋,O₋} inherit the order. Expanding det S over permutations then re-derives the total, every term having order at least Σ_c (d_c + d_{π(c)}) = 2Σ_c d_c. ∎

**The reading.** The frozen mode of the balanced pair (c, R(c)) is carried by a diagonal cell, and a diagonal cell sits at level d_c: as far from the anti-diagonal as the two sites are from each other. The census says that level owes one room, and that the only spare rooms in the whole block are on the anti-diagonal itself. So the debt has to be walked down to level 0, and the walk is the Hamiltonian's: one power of J is one hop, and no hop moves more than one level. The cost of the pair is its site distance; the Gram form charges it twice, once on each side of the determinant. This is the distance-is-t frame inside a determinant: what the excitation would need in time to cross from c to R(c) is what the determinant needs in powers of J before that frozen mode can move.

**Corollary (what the order counts).** At J = 0 the recentered block M̃ is diagonal, carrying −2(δ_a + δ_b) on the cell (a,b) and 4γ̄ on the diagonal cells, so on a generic locus profile it vanishes exactly on the 2⌊N/2⌋ anti-diagonal cells off the centre: there the frozen root has twice its generic multiplicity (gate, N = 3..6). Exactly ⌊N/2⌋ modes therefore leave the root as the coupling turns on, and since q(0) collects, up to sign, every eigenvalue of M̃ that is not frozen (Section 6), while only the departing ones among them vanish as J → 0, the theorem says the J-orders of the departing modes sum to at least 2Σ_c d_c. Under the measured equality the split is one departing mode per balanced pair, at order J^{2d_c}: the outermost pair, whose two sites are the ends of the chain, is the one that clings to the frozen root longest.

Two limits of the statement are worth naming. The theorem says only that the mode cannot move sooner; that it does move exactly then is the measured half, equivalent to S̃(0) being nonsingular, and it is what the equalities of Section 10 record. And the only place the chain geometry enters is Lemma 7, through the tridiagonal h: one power of J is one hop. Longer-range hopping crosses several levels per power and the ladder degrades accordingly, so the distance ladder is a statement about the chain, not about the mirror.

## 9. What it is not (the placement)

- **Not a decoherence-free structure.** The frozen eigenvectors are not in ker(ad_H): the kernel of K intersected with the anti-diagonal span is at most one-dimensional (one at even N, none at odd N), far below ⌊N/2⌋, and the frozen eigenvectors move with J. Only the eigenvalue stands still.
- **Not the uniform-γ commutant story.** At uniform γ (which sits on the locus as its fully degenerate point) the J-independent spectrum is the committed d_real ladder of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md), explained by weight-sector kernels ([absorption theorem](PROOF_ABSORPTION_THEOREM.md), [F50](PROOF_WEIGHT1_DEGENERACY.md)); those modes have J-independent eigenvectors. The frozen divisor is the site-resolved layer that survives when the profile is generic on the locus, and its mechanism is disjoint from the commutant.
- **Not a defective seed.** The frozen modes are semisimple (healthy left and right eigenvectors, overlap of order one); the Seed count of MirrorWorld concerns defective points, this concerns an eigenvalue pinned across a family. Siblings, not the same object.
- **The F139 kinship is the design lesson.** There the wall factor S₁₀ divides a character polynomial exactly, with no symmetry realizing the reflection; here (λ + 4γ̄)^⌊N/2⌋ divides the corner characteristic polynomial exactly on a locus, with no symmetry of the spectrum. Both walls are divisors.

## 10. Verification

The committed gate [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) checks, and must print "R90 frozen divisor gate: ALL GREEN":

- the block builder against the framework Liouvillian (sub-block equality, exact);
- the mirror identity with defect 8γ̄·P_D at machine zero for N = 4, 5, 6, and that the antilinear variant fails (the mirror is linear);
- the census: corner multiplicity ⌊N/2⌋ at −4γ̄ for N = 3..6 across J ∈ {0.6, 1, 2.3}, the (2,2) control empty (N = 4..6), the XY variant (h without the ZZ diagonal) and the antidiagonal corner at 4γ̄ − 2σ (the corollary), both at N = 4, 5;
- the pencil kernel dimensions (= ⌊N/2⌋) and the two eigenvector by-products;
- the partial-balance nulls;
- the N = 3 closed form, symbolically exact;
- at N = 6, exact Gaussian-rational arithmetic: the on-locus 36×36 corner determinant is exactly zero, and the transverse vanishing order in the defect is exactly 3 = ⌊N/2⌋;
- the cofactor theorem (Section 6): the closed form against the interpolated exact cofactor in Gaussian-rational arithmetic (N = 4, 5, Heisenberg; float cross-check XY N = 4), the symbolic N = 3 corner cofactor 2¹²γ̄²J⁴(3J² − δ₁²), and the nonvanishing of the leading coefficient det((K P_{O₊} K)|_{V₋}) for N = 3..10, Heisenberg and XY;
- the two boundary clocks (Section 7): the DCT-II / DST-I identification of the SE eigenbasis (machine zero, N = 3..10 both chains), the BBᵀ law with its {1, 1/M, 0} spectrum and pdet(G) = M^{−⌊(N−1)/2⌋}, the D_N closed form against the direct determinant, the exact rational assembly of the uniform cofactor at N = 4, 5 (Heisenberg, sympy discriminant), and uniform tightness at every sampled J (N = 3..6);
- the J-valuation discriminators (Section 6, last consequence): exact-rational cofactor ratios giving ord_J ≈ 18 at N = 6 (against 16) and ≈ 24 at N = 7 (against 20), a two-point estimate on exact arithmetic, since the orders themselves are pinned exactly by the next bullet; and the exact second-order rung rank(P_anti K P_{D₋}) = 1 at even N / 0 at odd N (N = 4..7);
- the valuation law (Section 8), in integer and rational arithmetic: the two grading properties of Lemma 7 and the level census of Lemma 8 against their closed forms (N = 3..10, both chains); the transport bound of Lemma 9 tested exhaustively over **all** maximal row sets (N = 3..6 Heisenberg, up to 816 minors), both as the inequality ord_J det Y_I ≥ its transport bound and as the two sharpness statements (the bound never falls below ⌊N²/4⌋, and the measured minimum over row sets attains it); the total ord_J det((X P_{O₊} X)|_{V₋}) = 2⌊N²/4⌋ for N = 3..10 on both chains, by exact interpolation rather than by ratio; the per-pair orders ord_J S_{c,c'} = d_c + d_{c'} for every pair, N = 3..7 Heisenberg; Lemma 6 entry by entry against an independent cell-basis build of the left side, and the same grading and valuation for a third h with non-uniform R-symmetric bonds and an arbitrary R-symmetric diagonal (N = 4..6), which is the generality the section claims; and the positive-weight Gram form of the leading Schur matrix with its positive definiteness (N = 4..6, the reduction recorded in Section 11). The one floating-point item in this group is the departure census of the Corollary (N = 3..6).

## 11. Open

- The upper half of the valuation law: Section 8 proves that the frozen modes cannot move before J^{d_c} and the measurements show that they move exactly then (N = 3..10 for the total, N = 3..7 per pair), which is the nonsingularity of the reduced matrix S̃(0). The open half now has a shape. Reading Cauchy-Binet at leading order gives

    S̃(0)·det T_{O₋,O₋}(0) = Σ_I w_I·ĝ_I(c)·ĝ_I(c'),  ĝ_I(c) := the coefficient of z^{d_c} in det Y_{I,A_c},

  where A_c is the column set of the proof above (the O₋ columns plus the D₋ column of the pair c), z := −iJ, and the w_I are strictly positive. In the variable z every entry of Y is real, since X = z·(iK) − 2Δ with iK and Δ real, and z is a unit times J so no valuation moves. The leading matrix is therefore a **real positive semidefinite Gram matrix**, and its determinant cannot cancel: nonsingularity is exactly the linear independence of the ⌊N/2⌋ vectors ĝ(c). The gate confirms the identity and positive definiteness at N = 4..6. What is missing is an argument that the ⌊N/2⌋ cheapest transports, one per distance, are independent as vectors over the row sets.
- The uniform-γ endpoint: how the frozen divisor's ⌊N/2⌋ modes embed into the enhanced d_real counts when all rate classes collapse. Section 7 supplies the frozen side in closed form; the d_real side of the ledger is still the committed open problem of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md).
- Adoption into MirrorWorld (the statement is finite linear algebra; the Sections 6-7 closed forms are entry-wise checkable up to one eigendecomposition of the N×N matrix h; candidate genre neighbour of Seed).
