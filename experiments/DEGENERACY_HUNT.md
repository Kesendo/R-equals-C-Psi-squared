# Degeneracy Hunt: Where Does the High Multiplicity Come From?

**Status:** Investigation complete. SU(2) broken by dephasing. Degeneracy structure not anomalous at N=5.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track A)
**Data:** `simulations/results/values_investigations/three_values_results.json`

---

## Motivation

SYMMETRY_CENSUS.md flagged max eigenvalue multiplicity = 14 at N=5 uniform chain. The known symmetries (U(1), spin-flip, reflection) predict at most 4x degeneracy (2 from flip x 2 from reflection). The gap between 4 and 14 was unexplained.

---

## 1. Multiplicity table (N=3-7, uniform chain, gamma=0.1)

| N | d^2 | Distinct eigenvalues | Max multiplicity | Count at max |
|---|-----|---------------------|-----------------|-------------|
| 3 | 64 | 26 | 6 | 2 |
| 4 | 256 | 127 | 14 | 1 |
| 5 | 1,024 | 488 | 14 | 2 |
| 6 | 4,096 | 2,207 | 19 | 2 |
| 7 | 16,384 | 8,136 | 22 | 2 |

**The sequence {6, 14, 14, 19, 22} is monotonically non-decreasing.** N=5 is not special; N=4 has the same max multiplicity. The high degeneracy is a structural feature of the Heisenberg + Z-dephasing Liouvillian at all N, not an anomaly at N=5.

## 2. Eigenvector inspection (N=5, eigenvalue Re = -0.400)

The 14 degenerate eigenvectors at Re(lambda) = -0.400 spread across multiple sectors:

| Sector (w_bra, w_ket) | Total weight |
|------------------------|-------------|
| (1,1) | 2.17 |
| (4,4) | 1.65 |
| (2,2) | 1.36 |
| (1,3) | 1.31 |
| (3,5) | 1.19 |
| (3,3) | 1.17 |

The degenerate modes live in many different sectors simultaneously. This is consistent with the [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) rate formula Re(lambda) = -2 sum gamma_k <1_XY(k)>: each eigenmode's decay rate is twice the dephasing-weighted average of its X/Y Pauli content. For uniform gamma, this simplifies to Re(lambda) = -2 gamma n_XY, placing all modes with the same n_XY count at the same rate. At the grid value Re = -0.400 = -4 gamma, modes from different (w_bra, w_ket) sectors coincide because they share n_XY = 2. This is a rate-formula coincidence, not a hidden symmetry.

## 3. SU(2) Casimir check

| Test | Result |
|------|--------|
| [S^2, Z_k] | norm = 16.0 for all k (NOT zero) |
| [S^2, H] | norm = 0.0 (SU(2) invariant Hamiltonian) |
| [C_{S^2}, L] | Frobenius norm = 28.6 (does NOT commute) |

**Conclusion:** SU(2) total spin is a symmetry of the Heisenberg Hamiltonian but is broken by Z-dephasing. The dephasing jump operators Z_k do not commute with S^2 (because Z only detects the z-component, not the total spin). SU(2) is therefore NOT a hidden symmetry of the Liouvillian.

The high degeneracies are explained by the absorption theorem rate formula, which places many modes from different sectors at the same grid values (multiples of 2*gamma for uniform chains). These are "accidental" degeneracies from the rate formula, not from a hidden symmetry.

---

## Files

- `simulations/three_values.py` (Track A computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (rate formula)
- [Symmetry Census](SYMMETRY_CENSUS.md) (flagged the degeneracy question)

---

*April 12, 2026. The 14-fold degeneracy is a rate-formula coincidence, not a hidden symmetry.*
