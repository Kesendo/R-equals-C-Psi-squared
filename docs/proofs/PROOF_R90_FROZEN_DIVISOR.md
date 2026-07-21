# The R₉₀ frozen divisor: watching profiles that pin an eigenvalue for every coupling

**Status:** Theorem (lower bound) Tier 1 derived; tightness Tier 2 verified (N = 3..6, exact arithmetic at N = 6, symbolic closed form at N = 3)
**Date:** 2026-07-22
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the anti-palindromic watching locus (every reflection pair of dephasing rates sums to the same value), the single-excitation corner block of the Liouvillian carries the eigenvalue λ = −4γ̄ with multiplicity at least ⌊N/2⌋, for every Hamiltonian coupling J (equality generic, verified N = 3..6). The mechanism is a rank bottleneck of a cell-level mirror, not an invariant subspace and not a spectral symmetry.
**Verification:** [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) (must print "R90 frozen divisor gate: ALL GREEN", ~2 min)
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
- **The mirror.** τQ is the linear involution of cells **(a,b) ↦ (Rb, Ra)** (transpose composed with reversal on both sides). It splits the diagonal cells D = {(a,a)} and the off-diagonal cells O = {(a,b), a≠b} into ±1 eigenspaces D±, O±.

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

The multiplicity is a dimension bottleneck of the mirror: the fixed cells of τQ (the anti-diagonal, one per balanced pair plus its transpose) give O₊ a surplus of 2⌊N/2⌋ rooms over O₋, the D₋ constraint taxes away half, and the remainder must freeze. Note what the proof does not use: no invariant subspace (the per-pattern kernels are trivial off the trivial patterns), no spectral palindromy (Section 2), no diagonalizability assumption.

Two structural by-products, both verified at machine precision in the gate: the frozen eigenvectors carry **zero weight on the diagonal cells** and have **τQ-even O-part**; the J-dependence of the frozen modes is only the ⌊N/2⌋-dimensional kernel line of Φ_J rotating inside the fixed subspace W.

## 4. Corollary: the antidiagonal corners, by the gamma fold

The block (1, N−1) (single excitation against N−1 excitations) is the image of the corner block under the one-sided X^N bridge, and the [gamma fold](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) turns that bridge into algebra: conjugating by right-multiplication with X^N sends L(γ⃗) to L(−γ⃗) − 2σ. The theorem applied at the profile −γ⃗ (the locus is preserved; the root becomes −4·(−γ̄) = +4γ̄) then gives, for every J,

  **the (1, N−1) block carries the eigenvalue 4γ̄ − 2σ with multiplicity ⌊N/2⌋.**

The two roots −4γ̄ and 4γ̄ − 2σ sum to −2σ: they are partners under the pair of mirrors, and recentered at the palindrome center (x = λ + σ) they sit at x = ±(σ − 4γ̄). At N = 4 the two coincide at the center x = 0. What once looked like two separate root families, including the N = 3 block (1,2) sighting at −2γ̄, is one theorem and its fold image: at N = 3, (1,2) is the antidiagonal corner, and 4γ̄ − 2σ = −2γ̄ there.

## 5. Why only the corners

The proof needs one affine recentering that makes the rate operator odd (Lemma 2). Under reversal a rate −2γ_S (S the disagreement set of the cell) goes to −2(2γ̄|S| − γ_S): the center depends on the size |S|. The corner block has a single off-diagonal size class (|S| = 2), so one center serves; the block (2,2) has classes |S| ∈ {0, 2, 4} and admits no single center, and indeed carries no frozen extras (gate control, N = 4, 5, 6). Away from the R₉₀ locus the whole structure disappears: partial balance (all but one condition satisfied) yields nothing, and the N = 3 closed form shows the defect as an explicit linear factor (M₍₁,₂₎ below is the (1,2) block of L, the N = 3 antidiagonal corner of Section 4),

  det(−2γ₂·I − M₍₁,₂₎) = 512·J⁴·(γ₁+γ₃)²·(γ₁+γ₃−2γ₂)·(4J² + (γ₁+γ₃)(γ₁+γ₃−2γ₂)),

so the balanced root exists exactly on the locus (given J ≠ 0), and the distance of the nearest eigenvalue grows linearly in the defect.

## 6. What it is not (the placement)

- **Not a decoherence-free structure.** The frozen eigenvectors are not in ker(ad_H): the kernel of K intersected with the anti-diagonal span is at most one-dimensional (one at even N, none at odd N), far below ⌊N/2⌋, and the frozen eigenvectors move with J. Only the eigenvalue stands still.
- **Not the uniform-γ commutant story.** At uniform γ (which sits on the locus as its fully degenerate point) the J-independent spectrum is the committed d_real ladder of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md), explained by weight-sector kernels ([absorption theorem](PROOF_ABSORPTION_THEOREM.md), [F50](PROOF_WEIGHT1_DEGENERACY.md)); those modes have J-independent eigenvectors. The frozen divisor is the site-resolved layer that survives when the profile is generic on the locus, and its mechanism is disjoint from the commutant.
- **Not a defective seed.** The frozen modes are semisimple (healthy left and right eigenvectors, overlap of order one); the Seed count of MirrorWorld concerns defective points, this concerns an eigenvalue pinned across a family. Siblings, not the same object.
- **The F139 kinship is the design lesson.** There the wall factor S₁₀ divides a character polynomial exactly, with no symmetry realizing the reflection; here (λ + 4γ̄)^⌊N/2⌋ divides the corner characteristic polynomial exactly on a locus, with no symmetry of the spectrum. Both walls are divisors.

## 7. Verification

The committed gate [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) checks, and must print "R90 frozen divisor gate: ALL GREEN":

- the block builder against the framework Liouvillian (sub-block equality, exact);
- the mirror identity with defect 8γ̄·P_D at machine zero for N = 4, 5, 6, and that the antilinear variant fails (the mirror is linear);
- the census: corner multiplicity ⌊N/2⌋ at −4γ̄ for N = 3..6 across J ∈ {0.6, 1, 2.3}, the (2,2) control empty (N = 4..6), the XY variant (h without the ZZ diagonal) and the antidiagonal corner at 4γ̄ − 2σ (the corollary), both at N = 4, 5;
- the pencil kernel dimensions (= ⌊N/2⌋) and the two eigenvector by-products;
- the partial-balance nulls;
- the N = 3 closed form, symbolically exact;
- at N = 6, exact Gaussian-rational arithmetic: the on-locus 36×36 corner determinant is exactly zero, and the transverse vanishing order in the defect is exactly 3 = ⌊N/2⌋.

## 8. Open

- The nonvanishing cofactor in closed form for general N (the N = 3 instance is the part of the Section 5 formula beside the defect factor; its general shape is not derived).
- Tightness (multiplicity exactly ⌊N/2⌋, generically) as a statement rather than an observation.
- The uniform-γ endpoint: how the frozen divisor's ⌊N/2⌋ modes embed into the enhanced d_real counts when all rate classes collapse.
- Adoption into MirrorWorld (the statement is finite linear algebra, entry-wise checkable, eigensolver-free; candidate genre neighbour of Seed).
