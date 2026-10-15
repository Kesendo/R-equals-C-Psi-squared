# Degeneracy Hunt: Where Does the High Multiplicity Come From?

**Status:** Finite-N census complete, and the maximal multiplicity is accounted for at N = 4..7. It belongs to a real eigenvalue on the grid, −2γk, whose eigenvectors are weight-k operators commuting with H: the commutant counts c_k of [Degeneracy Palindrome](DEGENERACY_PALINDROME.md#a-lower-bound-at-every-k-and-the-chains-commutant-counts), plus at N = 4 one weight-mixing mode. SU(2) is excluded as a Liouvillian symmetry; what it supplies is commutant elements, one orbit sum per letter multiset.
**Last refreshed:** 2026-10-15 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track A)
**Data:** `simulations/results/values_investigations/three_values_results.json`

---

## Motivation

[Symmetry Census](SYMMETRY_CENSUS.md) flagged max eigenvalue multiplicity = 14 at N=5 uniform chain. The known symmetries (U(1), spin-flip, reflection) predict at most 4× degeneracy (2 from flip × 2 from reflection). The gap between 4 and 14 was unexplained.

---

## 1. Multiplicity table (N=3-7, uniform chain, γ = 0.1)

| N | d² | Distinct eigenvalues | Max multiplicity | Count at max |
|---|-----|---------------------|-----------------|-------------|
| 3 | 64 | 26 | 6 | 2 |
| 4 | 256 | 127 | 14 | 1 |
| 5 | 1,024 | 488 | 14 | 2 |
| 6 | 4,096 | 2,207 | 19 | 2 |
| 7 | 16,384 | 8,136 | 22 | 2 |

**The measured sequence {6, 14, 14, 19, 22} is monotonically non-decreasing over N=3-7.** N=5 is not isolated within this finite range; no all-N monotonicity or multiplicity law follows from the table.

## 2. Eigenvector inspection (N=5, eigenvalue Re = -0.400)

The 14 degenerate eigenvectors at Re(λ) = −0.400 spread across multiple sectors:

| Sector (w_bra, w_ket) | Total weight |
|------------------------|-------------|
| (1,1) | 2.17 |
| (4,4) | 1.65 |
| (2,2) | 1.36 |
| (1,3) | 1.31 |
| (3,5) | 1.19 |
| (3,3) | 1.17 |

The table reports sector weights of one numerical eigenbasis. At a degeneracy,
an eigensolver may mix eigenvectors from invariant joint-popcount blocks, so
those cross-sector weights are basis-dependent. The basis-invariant statement
is that the same eigenvalue occurs in several blocks: projecting the degenerate
eigenspace with a commuting block projector produces block-supported
eigenvectors. The [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
gives Re(λ)=−2Σ_kγ_k⟨1_XY(k)⟩_v and, at uniform γ,
Re(λ)=−2γ⟨n_XY⟩_v. It accounts for the common real part −0.400 through
⟨n_XY⟩_v=2. The eigenvalue is real, −0.4 + 0i: with the script's H and γ the maximal
multiplicity sits at −0.4 and −0.6 for N = 5 (14-fold each), at −0.4 for
N = 4 (14) and at −0.4 and −0.8 for N = 6 (19 each), the real on-grid counts
d_real(k) of [Degeneracy Palindrome](DEGENERACY_PALINDROME.md). There are no
imaginary parts to match.
At N = 5 the 14 eigenvectors can be taken as weight-2 operators that commute
with H: the 12 orbit sums of the letter multisets with two X/Y letters, each
of which commutes with every SU(2)-invariant H and decays at exactly 4γ, and 2
more from the class-2 corner rows (1,1) and (4,4), which is the count
c_2 = 3N − 5 + 2⌊N/2⌋ = 14 of [Degeneracy Palindrome](DEGENERACY_PALINDROME.md#a-lower-bound-at-every-k-and-the-chains-commutant-counts).
The measured row {6, 14, 14, 19, 22} is the maximum of that document's d_real
rows at N = 3..7.

## 3. SU(2) Casimir check

| Test | Result |
|------|--------|
| [S², Z_k] | norm = 16.0 for all k (NOT zero) |
| [S², H] | norm = 0.0 (SU(2) invariant Hamiltonian) |
| [C_{S²}, L] | Frobenius norm = 28.6 (does NOT commute) |

**Conclusion:** SU(2) total spin is a symmetry of the Heisenberg Hamiltonian but is broken by Z-dephasing. The dephasing jump operators Z_k do not commute with S² (because Z only detects the z-component, not the total spin). SU(2) is therefore NOT a hidden symmetry of the Liouvillian.

The SU(2) candidate is excluded as a symmetry of L. Its trace survives in
the commutant: every letter-multiset orbit sum commutes with the
SU(2)-invariant H and is an exact eigenvector of uniform Z-dephasing, which
is where the high multiplicity comes from (§2).

---

## Files

- `simulations/three_values.py` (Track A computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)
- [Symmetry Census](SYMMETRY_CENSUS.md) (flagged the degeneracy question)

---

*April 12, 2026 finite-N census; SU(2) excluded as the missing symmetry.*

---

## § Update 2026-05-12: structural explanation in today's typed layer

The April-12 conclusion is correct but at the time we lacked the typed primitives that make the "rate formula coincidence" operational. With today's BlockSpectrum + SymmetryFamily infrastructure, the same observations have a precise structural reading:

1. **The "U(1)" mentioned in §1 was per-side U(1)**, predicting at most 2N+1 sectors. Today's [`JointPopcountSectors`](../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs) types the joint U(1)×U(1) per-side popcount conservation, giving (N+1)² block-diagonal sectors. The §2 weights show how one arbitrary numerical basis spans the degenerate eigenspace; they do not show physical mixing between these invariant sectors.

2. **X⊗N explains a subset of the cross-block recurrence.** [`XGlobalChargeConjugationPairing`](../compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs) (Tier1Derived 2026-05-12) types that under chain XY+Z-deph, sector `(p_c,p_r)` and `(N−p_c,N−p_r)` share spectra exactly. This supplies exact pairings such as `(1,1)↔(4,4)` and `(2,2)↔(3,3)`. It neither makes a numerical eigenvector's cross-sector weights invariant nor explains every block carrying the eigenvalue.

3. **The rate formula `Re(λ) = −2γ·⟨n_XY⟩_v` for a full eigenmode v at uniform γ** is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (typed as [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs)). It explains why modes in different blocks can share a decay rate. Only a basis coherence has the integer count `n_diff` outright, and neither equal `⟨n_XY⟩_v` nor equal real part forces the imaginary parts to coincide.

**Useful as forward regression test:** the multiplicity sequence `{6, 14, 14, 19, 22}` for N=3..7 at uniform γ is what `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` should reproduce under the same numerical clustering convention. A typed witness asserting the tolerance and multiplicity rule could replace the ad-hoc table at §1.

**Current conclusion:** SU(2) is not the missing Liouvillian symmetry. The maximal multiplicity is a real grid eigenvalue whose eigenspace is the weight-k commutant of H, counted exactly by c_k at N = 5, 6, 7 and by c_2 + 1 = 13 + 1 at N = 4, where one real mode mixes weights (§2). X⊗N pairs the blocks that carry it.

---

*The finite-N census and SU(2) exclusion remain the measured result.*
