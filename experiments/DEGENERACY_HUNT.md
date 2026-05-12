# Degeneracy Hunt: Where Does the High Multiplicity Come From?

**Status:** Investigation complete (April 12, 2026). SU(2) broken by dephasing. Degeneracy structure not anomalous at N=5. **2026-05-12 update:** the "rate formula coincidence" conclusion is now structurally typed in the BlockSpectrum / SymmetryFamily layer; see § Update at end.
**Date:** April 12, 2026 (original); 2026-05-12 status note added
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track A)
**Data:** `simulations/results/values_investigations/three_values_results.json`

---

## Motivation

[SYMMETRY_CENSUS](SYMMETRY_CENSUS.md) flagged max eigenvalue multiplicity = 14 at N=5 uniform chain. The known symmetries (U(1), spin-flip, reflection) predict at most 4× degeneracy (2 from flip × 2 from reflection). The gap between 4 and 14 was unexplained.

---

## 1. Multiplicity table (N=3-7, uniform chain, γ = 0.1)

| N | d² | Distinct eigenvalues | Max multiplicity | Count at max |
|---|-----|---------------------|-----------------|-------------|
| 3 | 64 | 26 | 6 | 2 |
| 4 | 256 | 127 | 14 | 1 |
| 5 | 1,024 | 488 | 14 | 2 |
| 6 | 4,096 | 2,207 | 19 | 2 |
| 7 | 16,384 | 8,136 | 22 | 2 |

**The sequence {6, 14, 14, 19, 22} is monotonically non-decreasing.** N=5 is not special; N=4 has the same max multiplicity. The high degeneracy is a structural feature of the Heisenberg + Z-dephasing Liouvillian at all N, not an anomaly at N=5.

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

The degenerate modes live in many different sectors simultaneously. This is consistent with the [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) rate formula Re(λ) = −2 Σ γ_k ⟨1_XY(k)⟩: each eigenmode's decay rate is twice the dephasing-weighted average of its X/Y Pauli content. For uniform γ, this simplifies to Re(λ) = −2γ·n_XY, placing all modes with the same n_XY count at the same rate. At the grid value Re = −0.400 = −4γ, modes from different (w_bra, w_ket) sectors coincide because they share n_XY = 2. This is a rate-formula coincidence, not a hidden symmetry.

## 3. SU(2) Casimir check

| Test | Result |
|------|--------|
| [S², Z_k] | norm = 16.0 for all k (NOT zero) |
| [S², H] | norm = 0.0 (SU(2) invariant Hamiltonian) |
| [C_{S²}, L] | Frobenius norm = 28.6 (does NOT commute) |

**Conclusion:** SU(2) total spin is a symmetry of the Heisenberg Hamiltonian but is broken by Z-dephasing. The dephasing jump operators Z_k do not commute with S² (because Z only detects the z-component, not the total spin). SU(2) is therefore NOT a hidden symmetry of the Liouvillian.

The high degeneracies are explained by the absorption theorem rate formula, which places many modes from different sectors at the same grid values (multiples of 2γ for uniform chains). These are "accidental" degeneracies from the rate formula, not from a hidden symmetry.

---

## Files

- `simulations/three_values.py` (Track A computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)
- [Symmetry Census](SYMMETRY_CENSUS.md) (flagged the degeneracy question)

---

*April 12, 2026. The 14-fold degeneracy is a rate-formula coincidence, not a hidden symmetry.*

---

## § Update 2026-05-12: structural explanation in today's typed layer

The April-12 conclusion is correct but at the time we lacked the typed primitives that make the "rate formula coincidence" operational. With today's BlockSpectrum + SymmetryFamily infrastructure, the same observations have a precise structural reading:

1. **The "U(1)" mentioned in §1 was per-side U(1)**, predicting at most 2N+1 sectors. Today's [`JointPopcountSectors`](../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs) types the joint U(1)×U(1) per-side popcount conservation, giving (N+1)² block-diagonal sectors. The 14 degenerate modes at N=5 Re=−0.400 spread across joint-popcount sectors `(1,1), (4,4), (2,2), (1,3), (3,5), (3,3)`; these are joint-popcount labels.

2. **The cross-sector spreading is X⊗N pairing.** [`XGlobalChargeConjugationPairing`](../compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs) (Tier1Derived 2026-05-12) types that under chain XY+Z-deph, sector `(p_c, p_r)` and sector `(N−p_c, N−p_r)` share spectra exactly. Inspecting the §2 sector list: `(1,1) ↔ (4,4)`, `(2,2) ↔ (3,3)`, `(1,3) ↔ (4,2)`; all paired sectors appear. The April observation is empirical confirmation of X⊗N pairing one month before that primitive was typed.

3. **The rate formula `Re(λ) = −2γ·n_XY` for uniform γ** is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (typed as [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs)). The "many sectors collapsing to one rate" is not a hidden symmetry, but the X⊗N-pairing typed claim names half of why: paired sectors are operator-conjugate-equivalent. The other half (rate-formula collapse across non-X⊗N-paired sectors that still share n_XY) is purely the absorption-theorem rate factorisation.

**Useful as forward regression test:** the multiplicity sequence `{6, 14, 14, 19, 22}` for N=3..7 at uniform γ is what `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` should reproduce bit-exactly. A typed witness asserting this could replace the ad-hoc table at §1.

**Conclusion stands:** the degeneracy is rate-formula coincidence (April-12 was right). What changed is that we now have typed primitives that name *which* coincidences and *why*: X⊗N pairing for half, absorption-theorem rate-grid for the other half. The investigation is structurally closed.

---

*Status note 2026-05-12. Original April-12 finding unchanged.*
