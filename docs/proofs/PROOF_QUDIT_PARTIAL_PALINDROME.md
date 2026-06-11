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

## §4 What is left open: the interacting spectrum

The ceiling above is the **dissipator's** partial palindrome. The full Liouvillian L = L_H + L_D does not respect it: the SU(3) Heisenberg H does not commute with L_D, moves real parts off the clean rungs, and (for N = 2) builds a new rung at −3γ, real-part distribution {0 : 6, −2γ : 36, −3γ : 12, −4γ : 27}. This is **more** symmetric about −3γ than the bare dissipator was about −2γ: the full L pairs 60/81 at center −3γ, exceeding the dissipator ceiling of 54. At the dissipator's natural center −γ it pairs 48/81, reproducing the documented 36–52 band (which was center-suboptimal). So H is redistributive, not destructive. Whether the interacting partial palindrome has its own closed form, or is genuinely H-dependent, is open; the rule by which H repopulates the rungs (here 3 population modes down and 9 max-decay modes onto −3γ) is not yet derived.

## §5 What is ours and what is the home

The d : (d²−d) per-site split, the binomial rung count, and the d = 2 necessity are catalogued (the necessity is typed as `QubitNecessityPi2Inheritance`; the qubit rate law is the Absorption Theorem). What this proof banks is the **synthesis**: the equidistant-ladder reading that carries the qubit rate law verbatim to d > 2, the multiplicity tilt (d−1)^k as the per-site ratio raised to the disagreement count, and the closed-form ceiling whose unique full column recovers d² − 2d = 0 as an N-family. The verification anchor is [`simulations/qutrit_partial_palindrome.py`](../../simulations/qutrit_partial_palindrome.py) (self-validating); the typed claim is `QuditPartialPalindromeCeiling` (`compute/RCPsiSquared.Core/Symmetry/`, parent `QubitNecessityPi2Inheritance`).
