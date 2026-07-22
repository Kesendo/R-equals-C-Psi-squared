# The R₉₀ frozen divisor: watching profiles that pin an eigenvalue for every coupling

**Status:** Theorem (lower bound) Tier 1 derived; cofactor closed form + tightness criterion Tier 1 derived (Section 6); uniform-endpoint constants in closed form and tightness-for-generic-J at every N Tier 1 derived (Section 7, the two boundary clocks)
**Date:** 2026-07-22
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the anti-palindromic watching locus (every reflection pair of dephasing rates sums to the same value), the single-excitation corner block of the Liouvillian carries the eigenvalue λ = −4γ̄ with multiplicity at least ⌊N/2⌋, for every Hamiltonian coupling J (equality generic, verified N = 3..6). The mechanism is a rank bottleneck of a cell-level mirror, not an invariant subspace and not a spectral symmetry. The nonvanishing cofactor is a single N(N−1)/2 determinant in closed form, (−1)^N(4γ̄)^⌈N/2⌉·det((X P_{O₊} X)|_{V₋}), whose nonvanishing is exactly tightness (Section 6).
**Verification:** [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) (must print "R90 frozen divisor gate: ALL GREEN", ~2-4 min)
**Depends on:** [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md) (the R₉₀ reshuffle and its fixed locus), [GAMMA_FOLD_PAIR_OF_MIRRORS](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) (the X^N cross-dock used in the corollary), [PROOF_F139_SEAM_IDENTITY](PROOF_F139_SEAM_IDENTITY.md) (the sibling on the character side)

---

## What this means

Tune the watching so that mirrored sites share their rates evenly around the mean: the first and the last site together, the second and the second-to-last together, every reflection pair carrying the same total. On that one locus, and only there, the chain acquires eigenvalues that the Hamiltonian cannot move. However hard the sites talk to each other, ⌊N/2⌋ decay modes sit frozen at one exact rate, one mode per balanced pair. They are not protected by a conserved quantity, and no subspace of states carries them; take any J and the eigenvectors have rearranged, yet the eigenvalue has not. What holds them is a bookkeeping fact: a mirror on the cells has more even rooms than odd rooms, and the surplus has nowhere to go.

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

The multiplicity is a dimension bottleneck of the mirror: the fixed cells of τQ (the anti-diagonal, one per balanced pair plus its transpose) give O₊ a surplus of 2⌊N/2⌋ rooms over O₋, the D₋ constraint taxes away half, and the remainder must freeze. Note what the proof does not use: no invariant subspace (the frozen eigenvectors move with J, only the eigenvalue stands still; Section 8's first bullet), no spectral palindromy (Section 2), no diagonalizability assumption.

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

Write M̃ := L_block(J) + 4γ̄ for the recentered block, and split off its τQ-even part: M̃ = X + 4γ̄·P_D, where

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
- **Tightness for generic J.** For a fixed locus profile, det(XP_{O₊}X|_{V₋}) is a polynomial in J of degree N(N−1) with leading coefficient det((K P_{O₊} K)|_{V₋}), which is nonzero for every N: Section 7 computes it in closed form (gate scan N = 3..10, both the Heisenberg and the XY single-excitation matrix). Hence the multiplicity is exactly ⌊N/2⌋ for all but finitely many J, at every N and every locus profile.
- **Small N in closed form** (symbolic, in the antisymmetric coordinates δ₁ = γ₁ − γ̄, δ₂ = γ₂ − γ̄):

    N = 3: q(0) = 2¹²·γ̄²·J⁴·(3J² − δ₁²)
    N = 4: q(0) = 2²⁰·γ̄²·J⁸·(8J⁴ − 4J²(3δ₁² + 2δ₁δ₂ + δ₂²) + (δ₁² − δ₂²)²)
    N = 5: q(0) = −2³⁰·γ̄³·J¹²·Q₅(J², δ₁, δ₂) with Q₅ of total degree 8 in (J, δ₁, δ₂), leading term 25J⁸

  (The XY chain differs only in the polynomial coefficients: its J-pure terms are 2J² and 5J⁴ at N = 3, 4.) The J-powers 4(N−2) in front are an observation, not yet derived; only the total degree N(N−1) and the leading coefficient come from the theorem.

## 7. The uniform endpoint: the two boundary clocks

At the fully degenerate point of the locus, δ⃗ = 0 (uniform watching), X = J·K and the cofactor collapses to a J-monomial:

  q(0)|_{δ⃗=0} = (−1)^N·(4γ̄)^⌈N/2⌉·J^{N(N−1)}·D_N,  D_N := det( (K P_{O₊} K)|_{V₋} ).

This section computes D_N in closed form, for every N. Three lemmas reduce it to a Gram determinant, and the Gram determinant is a pure power of the chain's clock modulus.

**Lemma 3 (the pair basis diagonalizes K²).** Let λ_1 < ... < λ_N and u_1, ..., u_N be the eigenvalues and orthonormal eigenvectors of h, with R-parities π_i (R u_i = π_i u_i; the spectrum is nondegenerate, see Lemma 5). Writing w_{ij} := u_i ⊗ u_j for the cell-space product vectors, τQ w_{ij} = π_iπ_j·w_{ji}, so the vectors φ_{ij} := (w_{ij} − π_iπ_j w_{ji})/√2 over i < j form an orthonormal basis of V₋ (count: N(N−1)/2 ✓), and K² φ_{ij} = −(λ_i−λ_j)²·φ_{ij}.

**Lemma 4 (the Gram reduction).** K is τQ-odd, so P₊ K|_{V₋} = K|_{V₋} (P₊ the projector onto V₊), hence K P_{O₊} K = K² − K P_{D₊} K = K² + C†C on V₋, with C := P_{D₊} K|_{V₋}. By the Weinstein–Aronszajn identity,

  D_N = (−1)^{N(N−1)/2} · Π_{i<j} (λ_i−λ_j)² · det( I − W ),  W := C Λ⁻² C†,  Λ² := diag((λ_i−λ_j)²).

Moreover I − W = V†V with V_{i,k} := ⟨w_{ii}, d_k⟩ (d_k the D₊ basis): K φ_{ij} = −i(λ_i−λ_j)·ψ_{ij} with ψ_{ij} := (w_{ij} + π_iπ_j w_{ji})/√2, so the (λ_i−λ_j)² of Λ⁻² cancels and W_{kl} = Σ_{i<j} ⟨d_k, ψ_{ij}⟩⟨ψ_{ij}, d_l⟩, and {w_{ii}} ∪ {ψ_{ij}} is an orthonormal basis of V₊ containing d_k, so W_{kl} = δ_{kl} − Σ_i ⟨d_k, w_{ii}⟩⟨w_{ii}, d_l⟩. Since the w_{ii} are τQ-even (orthogonal to D₋), the N×N Gram matrix G_{ij} := ⟨w_{ii}, P_D w_{jj}⟩ = Σ_a u_i(a)²u_j(a)² equals VV†, and therefore

  det(I − W) = det(V†V) = pdet(G)  (the product of the nonzero eigenvalues of G, provided rank G = dim D₊ = ⌈N/2⌉, which Lemma 5 gives).

**Lemma 5 (the two boundary clocks).** The single-excitation eigenbasis of the open chain is a cosine angle lattice, with modulus M depending on the chain (not the running clock of [ClockHandLadder](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs), whose two hands tick on ω and γ; here "clock" is the timeless lattice of angles the excitation lives on):

- **Heisenberg (M = N).** λ_k = 4cos(kπ/N) + N − 5 and u_k(a) ∝ cos((2a−1)kπ/(2N)), k = 0..N−1 (the DCT-II basis). The check is two lines: the interior rows force the dispersion, the end rows hold identically for the half-integer cosine (the chain's +2 boundary defect from the ZZ diagonal is exactly what the product-to-sum identity 2cosθ·cos(θ/2) = cos(3θ/2) + cos(θ/2) absorbs), and matching the right end quantizes sin(Nθ) = 0.
- **XY (M = N+1).** λ_k = 4cos(kπ/(N+1)) and u_k(a) ∝ sin(akπ/(N+1)), k = 1..N (the DST-I basis, classical).

Both spectra are nondegenerate (distinct cosines), so Lemma 3 applies. For B_{a,k} := u_k(a)² the three standard cosine sums (single sums = 0 past the flat mode; double products leave exactly the a = b and b = Ra resonances) give one law for both chains:

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

## 8. What it is not (the placement)

- **Not a decoherence-free structure.** The frozen eigenvectors are not in ker(ad_H): the kernel of K intersected with the anti-diagonal span is at most one-dimensional (one at even N, none at odd N), far below ⌊N/2⌋, and the frozen eigenvectors move with J. Only the eigenvalue stands still.
- **Not the uniform-γ commutant story.** At uniform γ (which sits on the locus as its fully degenerate point) the J-independent spectrum is the committed d_real ladder of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md), explained by weight-sector kernels ([absorption theorem](PROOF_ABSORPTION_THEOREM.md), [F50](PROOF_WEIGHT1_DEGENERACY.md)); those modes have J-independent eigenvectors. The frozen divisor is the site-resolved layer that survives when the profile is generic on the locus, and its mechanism is disjoint from the commutant.
- **Not a defective seed.** The frozen modes are semisimple (healthy left and right eigenvectors, overlap of order one); the Seed count of MirrorWorld concerns defective points, this concerns an eigenvalue pinned across a family. Siblings, not the same object.
- **The F139 kinship is the design lesson.** There the wall factor S₁₀ divides a character polynomial exactly, with no symmetry realizing the reflection; here (λ + 4γ̄)^⌊N/2⌋ divides the corner characteristic polynomial exactly on a locus, with no symmetry of the spectrum. Both walls are divisors.

## 9. Verification

The committed gate [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) checks, and must print "R90 frozen divisor gate: ALL GREEN":

- the block builder against the framework Liouvillian (sub-block equality, exact);
- the mirror identity with defect 8γ̄·P_D at machine zero for N = 4, 5, 6, and that the antilinear variant fails (the mirror is linear);
- the census: corner multiplicity ⌊N/2⌋ at −4γ̄ for N = 3..6 across J ∈ {0.6, 1, 2.3}, the (2,2) control empty (N = 4..6), the XY variant (h without the ZZ diagonal) and the antidiagonal corner at 4γ̄ − 2σ (the corollary), both at N = 4, 5;
- the pencil kernel dimensions (= ⌊N/2⌋) and the two eigenvector by-products;
- the partial-balance nulls;
- the N = 3 closed form, symbolically exact;
- at N = 6, exact Gaussian-rational arithmetic: the on-locus 36×36 corner determinant is exactly zero, and the transverse vanishing order in the defect is exactly 3 = ⌊N/2⌋;
- the cofactor theorem (Section 6): the closed form against the interpolated exact cofactor in Gaussian-rational arithmetic (N = 4, 5, Heisenberg; float cross-check XY N = 4), the symbolic N = 3 corner cofactor 2¹²γ̄²J⁴(3J² − δ₁²), and the nonvanishing of the leading coefficient det((K P_{O₊} K)|_{V₋}) for N = 3..10, Heisenberg and XY;
- the two boundary clocks (Section 7): the DCT-II / DST-I identification of the SE eigenbasis (machine zero, N = 3..10 both chains), the BBᵀ law with its {1, 1/M, 0} spectrum and pdet(G) = M^{−⌊(N−1)/2⌋}, the D_N closed form against the direct determinant, the exact rational assembly of the uniform cofactor at N = 4, 5 (Heisenberg, sympy discriminant), and uniform tightness at every sampled J (N = 3..6).

## 10. Open

- The inner structure of the residual polynomial Q_N away from the uniform point: the observed J^{4(N−2)} prefactor is not derived (the J-pure constants are now Section 7's clock numbers).
- The uniform-γ endpoint: how the frozen divisor's ⌊N/2⌋ modes embed into the enhanced d_real counts when all rate classes collapse. Section 7 supplies the frozen side in closed form; the d_real side of the ledger is still the committed open problem of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md).
- Adoption into MirrorWorld (the statement is finite linear algebra; the Sections 6-7 closed forms are entry-wise checkable up to one eigendecomposition of the N×N matrix h; candidate genre neighbour of Seed).
