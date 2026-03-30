# Sacrifice-Zone Qubit Mapping: Finding Optimal Chains on Real Hardware

<!-- Keywords: sacrifice zone qubit selection IBM Torino, heavy-hex topology
chain optimization mode protection, T2 calibration data cavity mode
localization, spatial noise profile quantum advantage, palindromic
mode survival dephasing chain, R=CPsi2 sacrifice zone mapping -->

**Status:** Tier 2-3 (computed analysis of real calibration data, hardware
test pending)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)
**Script:** [sacrifice_zone_mapping.py](../simulations/sacrifice_zone_mapping.py)
**Data:** [sacrifice_zone_mapping.txt](../simulations/results/sacrifice_zone_mapping.txt)
**Calibration data:** [ibm_torino_history.csv](../data/ibm_history/ibm_torino_history.csv) (24,073 records, 181 days, 133 qubits)
**Topology:** Heavy-hex via Qiskit `CouplingMap.from_heavy_hex(7)` (115 qubits, 132 edges)

---

## Abstract

If the sacrifice zone protects cavity modes localized on interior
qubits (r = 0.994, see [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)),
then choosing qubit chains where a naturally noisy qubit sits at the
edge should provide the sacrifice-zone benefit *for free*.

We test this on IBM Torino's heavy-hex topology (133 qubits, 115 with
T2 data on the latest calibration date, 132 edges) using real T2
calibration data (181 days, 24,073 records). 330 five-qubit chains
exist on the graph. We compare two chain selection strategies:

1. **Sacrifice-zone ranking:** Maximize edge noise / interior noise ratio
2. **Mean-T2 ranking:** Maximize average T2 across all 5 qubits

Result: **Zero overlap** in the top-10 lists. Sacrifice-zone chains
achieve **2.54x** mean protection factor vs **1.18x** for mean-T2
chains. Mode-based selection outperforms naive T2 maximization by 2.15x.

The best sacrifice chain has only 81 us mean T2 but 2.86x protection.
The best T2 chain has 217 us mean T2 but only 1.06x protection.
**Worse qubits, better modes.**

---

## Method

### Chain selection

On the heavy-hex graph, we enumerate all simple paths of length 5
(5 qubits, 4 bonds). For each chain [q0, q1, q2, q3, q4]:

**Sacrifice score** = max(γ_edge) / mean(γ_interior)

where γ_edge = max(γ[q0], γ[q4]) and γ_interior = mean(γ[q1], γ[q2], γ[q3]).
Higher score means the edge qubit absorbs more noise relative to the
interior, strengthening the sacrifice-zone effect.

### Spectral verification

For the top-5 chains in each ranking, we compute the full Liouvillian
(N=5 chain, real γ values) and extract: palindrome score, slowest
oscillating mode rate, and protection factor vs uniform noise.

---

## Results

### Sacrifice-zone ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [85, 15, 86, 16, 87] | 18.0 | 112.5 us | 2.52x |
| [85, 15, 86, 58, 92] | 17.0 | 94.3 us | 2.53x |
| [49, 7, 79, 53, 85] | 15.8 | 67.6 us | 2.63x |
| [80, 8, 79, 53, 85] | 11.9 | 81.1 us | **2.86x** |
| [85, 14, 57, 21, 91] | 10.6 | 82.8 us | 2.14x |

All contain Q85 (T2 = 5.0 us), the noisiest qubit on the chip,
as the sacrifice endpoint.

### Mean-T2 ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [18, 89, 19, 90, 60] | 1.0 | 217.3 us | 1.06x |
| [88, 18, 89, 19, 90] | 1.1 | 207.6 us | 1.12x |
| [19, 90, 60, 96, 26] | 0.7 | 203.1 us | 1.34x |
| [4, 76, 51, 82, 10] | 0.9 | 202.3 us | 1.16x |
| [13, 56, 20, 90, 60] | 0.8 | 198.5 us | 1.22x |

All have sacrifice scores near 1.0 (uniform noise). The quiet qubits
provide long T2 but no differential protection.

### Head-to-head

| Metric | Sacrifice top-5 | Mean-T2 top-5 |
|--------|----------------|---------------|
| Mean protection factor | **2.54x** | 1.18x |
| Mean T2 | 87.6 us | 205.8 us |
| Mean sacrifice score | 14.8 | 0.9 |
| Palindrome score | 98-100% | 88-92% |

The sacrifice chains have 2.4x lower mean T2 but 2.15x higher
protection. Choosing "worse" qubits with the right spatial pattern
outperforms choosing the "best" qubits naively.

**Note on palindrome scores:** Sacrifice chains show 98-100%, mean-T2
chains 88-92%. The palindromic theorem is exact for Z-dephasing at
any noise level. The lower scores for mean-T2 chains arise from
numerical precision limits at very low total noise (Σγ = 0.02-0.03
vs 0.25 for sacrifice chains). At such small Σγ values, the
palindromic center is close to zero and eigenvalue pairing tolerances
(1e-4) become significant relative to the eigenvalue spread. This is
a numerical artifact, not a physical effect.

---

## Time stability

The best sacrifice chain [85, 15, 86, 16, 87] tracked across 5 months:

| Date | Score | mean T2 |
|------|-------|---------|
| 2026-02-10 | 14.9 | 112.5 us |
| 2025-12-12 | 13.0 | 102.4 us |
| 2025-10-13 | 11.3 | 118.3 us |

The score varies by ~30% but the chain consistently ranks at the top.
Q85 remains the noisiest qubit on the chip across all calibrations.
The mapping does not need daily recalculation.

---

## What this means

Standard quantum computing practice: select the qubits with the
highest T2 values and hope for the best. This ignores the spatial
structure of the noise.

Mode-based approach: select chains where a naturally noisy qubit
sits at the edge, creating a built-in sacrifice zone. The noisy
qubit absorbs disproportionate damping, and the cavity modes
localized on the interior survive longer.

This requires no additional gates, no error correction, and no
knowledge of the palindromic theory. It is a free improvement
available on any quantum processor with non-uniform noise
characteristics. The only input is the coupling map and the T2
calibration data, both publicly available.

The theory predicts which chains will perform best. The prediction
is testable with a single set of Trotter evolution experiments
comparing sacrifice-zone chains against mean-T2 chains on the
same hardware on the same day.

---

*See also:*
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) (2.81x theoretical),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md) (1.97x measured),
[Resonant Return](RESONANT_RETURN.md) (the sacrifice-zone formula)
