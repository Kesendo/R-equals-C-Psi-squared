# PROOF: the partial palindrome at d>2 is the symmetric overlap of the disagreement count

**Status:** Tier 1 derived (one-line combinatorial identities, exact; verified machine-exact for the dissipator spectrum at d = 3, N = 2 and against brute enumeration on a (d, N) grid).
**Date:** 2026-06-11
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [QUBIT_NECESSITY](../QUBIT_NECESSITY.md): the per-site split d : (d²−d) and the polynomial trunk d² − 2d = 0 that closes only at d = 2 (typed as `QubitNecessityPi2Inheritance`).
- [PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md): the qubit rate law Re(λ) = −2γ·n_XY = −2γ·Hamming(i, j), the disagreement-count reading the dissipator gives the spectrum.
- [ON_THE_ONE_DIAGONAL](../../reflections/ON_THE_ONE_DIAGONAL.md): the dissipator IS the recentred disagreement-count diagonal; this proof is that diagonal read one dimension up.

## Abstract

The palindromic mirror is exact only for qubits (d = 2); [QUBIT_NECESSITY](../QUBIT_NECESSITY.md) proves this from the per-site balance d = d² − d, i.e. d² − 2d = 0, which closes only at d = 2. For d > 2 the spectrum is not random: N = 2 qutrits were observed to pair 36–52 of 81 eigenvalues, a residual structure no principle captured (OQ-002). This proof captures it. Under full-Cartan dephasing (the complete diagonal Cartan subalgebra as jump operators), the d levels are mutually **equidistant**, so the decay rate of a coherence |i⟩⟨j| is exactly −2γ·Hamming(i, j), the **same rate ladder as the qubit**. What differs is the multiplicity per rung: the number of coherences at Hamming distance k is

  **c_k = d^N · C(N, k) · (d−1)^k.**

The palindrome reflects rung k against rung N−k. For d = 2 the factor (d−1)^k = 1 and c_k = 2^N·C(N, k) is symmetric in k, so every rung finds its mirror: 100%. For d > 2 the (d−1)^k tilts the distribution toward large k, and only the symmetric overlap pairs. The dissipator's paired ceiling is

  **paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)},**

which equals d^{2N} (full) iff d = 2. For d = 3, N = 2 it is 54/81; the qutrit fraction erodes with N (66.7%, 66.7%, 51.9% for N = 1, 2, 3); d = 4 gives 50%, 50%, 31.2%. The (d−1) tilt base is exactly the per-site decaying-to-immune ratio (d²−d) : d = (d−1) : 1 from [QUBIT_NECESSITY](../QUBIT_NECESSITY.md), raised to the number of disagreeing sites. The full interacting Liouvillian is richer and is left open in §4.

## §1 The equidistant ladder

Take as jump operators the diagonal Cartan generators (for qutrits λ₃ = diag(1, −1, 0) and λ₈ = diag(1, 1, −2)/√3; generally a basis of the diagonal subalgebra). Each is diagonal, so a coherence |i⟩⟨j| is an eigenmode of the single-site dissipator with rate

  Σ_M [M_{ii}M̄_{jj} − ½(|M_{ii}|² + |M_{jj}|²)] = −½ Σ_M (M_{ii} − M_{jj})².

For the qutrit Cartan {λ₃, λ₈} this evaluates to −2 for every i ≠ j: (0,1) gives ½(4 + 0); (0,2) gives ½(1 + 3); (1,2) gives ½(1 + 3); all equal 2. The three levels are mutually equidistant. Hence per site the rate is −2γ if i ≠ j and 0 if i = j, and across N sites

  rate(|i⟩⟨j|) = −2γ · #{l : i_l ≠ j_l} = −2γ · Hamming(i, j),

the identical ladder to the qubit's Re(λ) = −2γ·n_XY ([PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md)). The diagonal generators are immune (rate 0); the d² − d off-diagonal generators decay. This is the per-site d : (d²−d) split of [QUBIT_NECESSITY](../QUBIT_NECESSITY.md), now read as a rung of the rate ladder.

## §2 The multiplicity tilt

Count the coherences at Hamming distance k. Choose the k disagreeing sites: C(N, k). On each, j_l ≠ i_l is one of d−1 values: (d−1)^k. The ket i ranges over all d^N states. Hence

  **c_k = d^N · C(N, k) · (d−1)^k,**   Σ_k c_k = d^N · Σ_k C(N,k)(d−1)^k = d^N · d^N = d^{2N}.

The factor (d−1)^k is the only place the dimension enters non-trivially. Its base d−1 is the per-site decaying : immune ratio (d² − d)/d, raised to the number of disagreements. At d = 2 it is 1^k = 1 and c_k = 2^N·C(N, k), the symmetric binomial.

## §3 The ceiling, and why only d = 2 is full

The dissipator spectrum is real, with rung k at rate −2γk and multiplicity c_k. The palindrome reflects λ ↦ −2·(Nγ) − λ about the center −Nγ, i.e. rung k against rung N−k. Modes at rung k find a partner iff there is a mode at rung N−k; the number that pair is, summing each two-rung pair and the self-mirrored middle (N even),

  paired(d, N) = Σ_{k < N/2} 2·min(c_k, c_{N−k}) + [N even]·c_{N/2}.

Because C(N, k) = C(N, N−k) and d − 1 ≥ 1, min(c_k, c_{N−k}) = d^N·C(N, k)·(d−1)^{min(k, N−k)}, so

  **paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)}.**

This is d^{2N} (everything pairs) iff (d−1)^{min(k,N−k)} = (d−1)^k for all k, i.e. iff d − 1 = 1, i.e. **d = 2**. For d > 2 the high rungs (k > N/2) carry more modes than their low-rung partners; the excess

  Σ_{k > N/2} d^N·C(N, k)·[(d−1)^k − (d−1)^{N−k}]

is unpaired. For d = 3, N = 2: c = [9, 36, 36], rung 0 (×9) pairs into rung 2 (×36) leaving 27, rung 1 (×36) self-mirrors; paired = 54, excess = 27. The qutrit fraction erodes with N (66.7%, 66.7%, 51.9%); d = 4 gives 50%, 50%, 31.2%. The d² − 2d = 0 uniqueness of [QUBIT_NECESSITY](../QUBIT_NECESSITY.md) reappears here as the unique fully-paired column of an N-indexed family.

## §4 The interacting spectrum: H degrades the palindrome

The ceiling above is the **dissipator's** palindrome, taken about the physical center −Nγ (the k ↔ N−k reflection, where the qubit palindrome is exact). Adding the Hamiltonian degrades it. The full Liouvillian L = L_H + L_D does not commute, and at **every fixed center H reduces the pairing**: about the physical center −Nγ = −2γ the dissipator's 54 drops to 48; about −3γ (where the two large rungs −2γ and −4γ are equinumerous, the dissipator's own best fit at 72) it drops to 60. There is no center at which adding H helps; the palindrome is fragile under H, the same fragility [QUBIT_NECESSITY](../QUBIT_NECESSITY.md) records for the full qubit mirror at N ≥ 3. (An earlier reading of this work mistook the full L's best-fit 60 at −3γ for "exceeding" the dissipator's 54 at −2γ: that compared two different centers. At equal center, H always reduces the count.)

What H does carry is a clean **real-part law in the symmetric case**. For the SU(3) Heisenberg the real parts sit exactly on Re(λ) = −2γ⟨Q⟩, where ⟨Q⟩ is the biorthogonal Hamming-distance expectation over each eigenmode (the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md)'s Rayleigh reading: the real part comes from the self-adjoint dissipator alone). The SU(3) symmetry quantizes ⟨Q⟩ into {0, 1, 1.5, 2} with multiplicities {6, 36, 12, 27}; the new −3γ rung is exactly ⟨Q⟩ = 1.5, a 50/50 mix of Hamming-1 and Hamming-2 coherences. This exactness is a symmetry effect, not a general law: a random Hermitian H breaks it (the real parts spread off the ladder, deviation ≈ 10⁻³).

So the interacting partial palindrome has **no H-independent closed form**: the paired count floats with the Hamiltonian (60 for the SU(3) Heisenberg, robust across all coupling ratios J/γ; 48 and 52 for less symmetric couplings; near zero for a generic H). The dissipator's closed form (54 about −Nγ, full iff d = 2) is the only invariant skeleton; the interacting count is a property of each H's symmetry. Verified in [`simulations/qutrit_interacting_palindrome.py`](../../simulations/qutrit_interacting_palindrome.py) (self-validating).

## §5 What is ours and what is the home

The d : (d²−d) per-site split, the binomial rung count, and the d = 2 necessity are catalogued (the necessity is typed as `QubitNecessityPi2Inheritance`; the qubit rate law is the Absorption Theorem). What this proof banks is the **synthesis**: the equidistant-ladder reading that carries the qubit rate law verbatim to d > 2, the multiplicity tilt (d−1)^k as the per-site ratio raised to the disagreement count, and the closed-form ceiling whose unique full column recovers d² − 2d = 0 as an N-family. The verification anchor is [`simulations/qutrit_partial_palindrome.py`](../../simulations/qutrit_partial_palindrome.py) (self-validating); the typed claim is `QuditPartialPalindromeCeiling` (`compute/RCPsiSquared.Core/Symmetry/`, parent `QubitNecessityPi2Inheritance`).

## §6 The operator realization: the product cap, Π's formula one dimension up, and the wreath family (added 2026-06-11, same day)

The sections above count the pairs; this section asks **which operator performs the pairing**, and the answer closes the review's open question with the trunk polynomial appearing a third time.

**Theorem (the product cap).** Any per-site mirror, W = ⊗_l q_l with arbitrary site-dependent q_l, one-sided or two-sided, linear or antilinear, that intertwines the dissipator palindrome W L_D = (−L_D − 2Nγ)·W on its support, pairs at most

  **(2d)^N of the d^{2N} coherences.**

*Proof.* The rate ladder is per-site additive, so the intertwining forces each q_l to be a strict class swap: q_l(dark) ⊆ lit and q_l(lit) ⊆ dark, where dark = {|i⟩⟨i|} (d dims, rate 0) and lit (d² − d dims, rate −2γ) per site. (If some q_l(e) mixed the characters, fixing the other sites' letters dark would put the image across two rungs, breaking the intertwining; product coefficients cannot cancel across independent site choices.) Then rank(q_l) = rank(q_l|dark) + rank(q_l|lit) ≤ min(d, d²−d) + min(d, d²−d) = 2d for d ≥ 2, and rank(W) = Π_l rank(q_l) ≤ (2d)^N. ∎

The cap is the full space iff (2d)^N = d^{2N}, i.e. iff **d² − 2d = 0**: the [QUBIT_NECESSITY](../QUBIT_NECESSITY.md) trunk, now as an operator bound (third appearance, after the per-site split and the ceiling column). And the cap is strictly below §3's combinatorial ceiling for every d ≥ 3, N ≥ 2 (at d = 3, N = 2: cap 36 < ceiling 54): since a global, non-product partial isometry reaches the ceiling exactly (greedy rung matching, exact intertwining on its support; verified), **the gap ceiling − (2d)^N is precisely the non-product part of the partial palindrome**. This is the inverse of the [golden-router story](PROOF_CEILING_GOLDEN_ROUTER.md): there the suspected non-locality dissolved because window sums gave per-site maps room the per-term test could not see; here the dissipator is strictly local, the rung count is rigid, and the non-locality is provable.

**The operator attaining the cap is the qubit palindromizer's own formula.** Define, exactly as in [F118](PROOF_PI_FACTORS_AS_R_TIMES_D.md) (Π_Z(ρ) = ρᵀ·X^⊗N),

  **Π_d(ρ) = ρᵀ · Shift^⊗N,**  Shift the cyclic clock shift |x⟩ ↦ |x+1⟩.

Per site Π_d sends the letter (i, j) to (j, i−1): dark (x, x) ↦ (x, x−1), lit-aligned (a, a−1) ↦ (a−1, a−1) dark. On the **shift-aligned subspace**, per-site span{(x,x)} ∪ {(a, a−1)} of dimension 2d per site, hence (2d)^N in all, Π_d is closed and the intertwining residual is **exactly zero** (verified at d = 3, N = 1..3 and d = 4, N = 1..2); on the complement it fails at O(γ), which is the provably unpaired part. There are two chiralities Π_d^± (the two shift directions, aligned subspaces {(a, a∓1)}); at d = 2 the two off-diagonals coincide, the chiralities merge, the aligned subspace is everything, and the mirror is full: **the qubit magic of §3 is, in operator language, the degeneracy of the two shift chiralities.**

**The mirror group becomes a wreath family.** With D the transpose (F118's reflection), the closure obeys, verified exactly at d = 2, 3, 4, 5:

  **ord(Π_d) = 2d,  |⟨Π_d, D⟩| = 2d²,  ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂**

(the elements are the index swaps-or-not composed with independent two-sided shifts (a, b) ∈ Z_d × Z_d; D-conjugation exchanges the two shift factors rather than inverting them, which is what makes the extension a wreath product and not a generalized dihedral group). At d = 2 this is exactly D₄: **the [F118 mirror group](PROOF_PI_FACTORS_AS_R_TIMES_D.md) is the d = 2 column of a d-indexed wreath family**, and for d ≥ 3 the reflection D swaps the Π_d^± chiralities instead of normalizing one mirror, the group-level face of the partiality.

The verification anchor for this section is [`simulations/qudit_product_mirror_cap.py`](../../simulations/qudit_product_mirror_cap.py) (self-validating, six blocks: cap + trunk arithmetic d = 2..5; exact alignment residuals; random class-swap rank compliance; the global ceiling-reacher; the group law; the d = 2 degeneracy). The typed claim is `QuditProductMirrorCap` (`compute/RCPsiSquared.Core/Symmetry/`, parents `QuditPartialPalindromeCeiling` + `QubitNecessityPi2Inheritance`). Open after this section: whether structured intermediates between product and global (translation-invariant non-product W's, the qudit analogue of the window-summed router) can exceed (2d)^N, and the SU(3)-Heisenberg interacting count 60 of §4 as a representation-theory exercise.
