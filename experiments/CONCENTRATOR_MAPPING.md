# Concentrator Qubit Mapping: Finding Optimal Chains on Real Hardware

**Naming note (2026-07-05):** renamed from "Sacrifice-Zone Qubit Mapping". The
noisy edge qubit sacrifices nothing; it concentrates the noise (the misnomer
was resolved 2026-03-28). The frozen `sacrifice_zone_mapping.*` script and data
keep their original names as the provenance of the run.

<!-- Keywords: concentrator qubit selection IBM Torino, heavy-hex topology
chain optimization mode protection, T2 calibration data cavity mode
localization, spatial noise profile quantum advantage, palindromic
mode survival dephasing chain, R=CPsi2 concentrator mapping -->

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
**Topology:** Heavy-hex (IBM's standard qubit layout where each node connects to 2 or 3 neighbors in a hexagonal pattern with extra "bridge" qubits on each edge) via Qiskit `CouplingMap.from_heavy_hex(7)` (115 qubits, 132 edges)

---

## What this document is about

This document shows how to find optimal qubit chains on real IBM hardware
by exploiting naturally noisy qubits as concentrators. Instead of
picking the qubits with the best T2 times (the standard approach), we
select chains where a noisy qubit sits at the edge, providing the
concentrator benefit for free. On IBM Torino's 133-qubit chip, this
mode-based selection outperforms naive T2 maximization by 2.15× in
protection factor, despite using qubits with 2.3× lower average T2.

---

## Abstract

If the concentrator protects cavity modes localized on interior
qubits (r = 0.994, see [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)),
then choosing qubit chains where a naturally noisy qubit sits at the
edge should provide the concentrator benefit *for free*.

We test this on IBM Torino's heavy-hex topology using real T2
calibration data (181 days, 24,073 records). On the latest calibration
date 133 qubits carry T2 data; the `from_heavy_hex(7)` coupling map used
here covers 115 of them, over 132 edges. 330 five-qubit chains
exist on the graph. We compare two chain selection strategies:

1. **Concentrator ranking:** Maximize edge noise / interior noise ratio
2. **Mean-T2 ranking:** Maximize average T2 across all 5 qubits

Result: **Zero overlap** in the top-10 lists. Concentrator chains
achieve **2.53x** mean protection factor vs **1.18x** for mean-T2
chains. Mode-based selection outperforms naive T2 maximization by 2.15x.

The best concentrator chain has only 81 us mean T2 but 2.86x protection.
The best T2 chain has 217 us mean T2 but only 1.06x protection.
**Worse qubits, better modes.**

---

## Method

### Chain selection

On the heavy-hex graph, we enumerate all simple paths of length 5
(5 qubits, 4 bonds). For each chain [q0, q1, q2, q3, q4]:

**Concentrator score** = max(γ_edge) / mean(γ_interior)

where γ_edge = max(γ[q0], γ[q4]) and γ_interior = mean(γ[q1], γ[q2], γ[q3]).
Higher score means the edge qubit absorbs more noise relative to the
interior, strengthening the concentrator effect.

### The rate convention, and what it rests on

**γ[q] = 1/(2·T₂[q])**. The model here is dephasing-only, its only jump
operators being √γ_q·Z_q, and a D[Z] channel at rate γ decays coherences at 2γ,
so that is the rate reproducing the measured T₂. See
[the glossary section](../docs/GLOSSARY.md) for the conversion and for when the
T₁-aware form is the right one instead.

The score is a ratio of γ values, so it is unchanged by a global factor on all
rates. That is **not** the same as being independent of the model, and the
difference matters here: the T₁-aware form is not a global factor, because the
T₁ share of the decay is wildly uneven across this chip. On Q85, the sacrifice
qubit in four of the five headline chains, T₁ = 2.9 µs against T₂ = 5.0 µs, so
the 1/(2T₁) term carries most of the coherence decay, while the interior qubits
sit near 40%. Under that model Q85's D[Z] rate would fall several-fold more than
the interior's and the ranking would not survive. A σ⁻ channel also breaks the Π
mirror (F82/F84), so the palindrome column would not survive either. Everything
below is a statement about the dephasing-only model these scripts implement.

Related, and unhandled by either script: Q53, which appears in two of the five
headline chains, reports T₂ = 62.4 > 2·T₁ = 44.8 on this calibration date. That
is a broken record (the typed layer clamps it, `IbmCalibration.cs`), and the
dephasing-only form is silent about it only because it never reads T₁.

### Spectral verification

For the top-5 chains in each ranking, we compute the full Liouvillian
(N=5 chain, real γ values) and extract: palindrome score, slowest
oscillating mode rate, and protection factor vs uniform noise.

---

## Results

### Concentrator ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [85, 15, 86, 16, 87] | 18.0 | 112.5 us | 2.51x |
| [85, 15, 86, 58, 92] | 17.0 | 94.3 us | 2.53x |
| [49, 7, 79, 53, 85] | 15.8 | 67.6 us | 2.63x |
| [80, 8, 79, 53, 85] | 11.9 | 81.1 us | **2.86x** |
| [85, 14, 57, 21, 91] | 10.6 | 82.8 us | 2.14x |

All contain Q85 (T2 = 5.0 us), the noisiest qubit on the chip,
as the concentrator endpoint.

### Mean-T2 ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [18, 89, 19, 90, 60] | 1.0 | 217.3 us | 1.06x |
| [88, 18, 89, 19, 90] | 1.1 | 207.6 us | 1.12x |
| [19, 90, 60, 96, 26] | 0.7 | 203.1 us | 1.34x |
| [4, 76, 51, 82, 10] | 0.9 | 202.3 us | 1.16x |
| [13, 56, 20, 90, 60] | 0.8 | 198.5 us | 1.22x |

All have concentrator scores near 1.0 (uniform noise). The quiet qubits
provide long T2 but no differential protection.

### Head-to-head

| Metric | Concentrator top-5 | Mean-T2 top-5 |
|--------|----------------|---------------|
| Mean protection factor | **2.53x** | 1.18x |
| Mean T2 | 87.6 us | 205.8 us |
| Mean concentrator score | 14.6 | 0.9 |
| Palindrome score | 96-98% | 85-92% |

The concentrator chains have 2.3x lower mean T2 but 2.15x higher
protection. Choosing "worse" qubits with the right spatial pattern
outperforms choosing the "best" qubits naively.

**The palindrome column is void, in both conventions.** It reports 96-98% for
concentrator chains and 85-92% for mean-T2 chains, and neither number measures
the palindrome. Measured on the mean-T2 top chain [18, 89, 19, 90, 60], the
palindromic symmetry is **exact**: comparing the 960 oscillatory rates against
their mirror image 2·Σγ − rate as sorted multisets gives a residual of
**1.8e-14**. The theorem holds to machine precision, as it must.

What the percentage measures is the matcher. `spectral_analysis` pairs each rate
with the **first** partner it finds within an **absolute** 1e-4, and the spectrum
is far denser than that: **927 of the 959** nearest-neighbour gaps are below
1e-4. So a rate almost never grabs its true partner, and it orphans that
partner's mate. Of the 410 matches the scorer accepts on that chain, **362 are
wrong by more than 1e-12**. The column measures how a greedy first-fit scrambles
inside the level clustering, which is why it can move in either direction: across
the ten chains the 2026-08-05 γ-book repair moved it by 0, 1, 2, 3 and 4 points,
and **upward** on two of them.

The repair this wants is therefore **not** a scale-relative tolerance. An exact
route exists, so the residual should be read rather than gated (the repo's "a
deviation is a sign" rule, case 1): the sorted-multiset comparison above gives
the answer with no threshold in it. Any tolerance below the minimum level spacing
returns exactly 100% and hides the same clustering, and 1e-6, 1e-8 and a
scale-relative 1e-4·Σγ all do. Until the scorer is replaced, treat these
percentages as carrying no information, ordinal or otherwise.

**What the γ-book repair changed here.** The 2026-08-05 repair halved every rate
in this run (γ = 1/T₂ → the D[Z] rate 1/(2T₂); see
[the glossary section](../docs/GLOSSARY.md)). This is **not** a change of units:
J is held at 1, so halving γ doubles Q = J/γ and moves the physical point. What
survived exactly are the quantities of γ-degree zero, and only those: the
concentrator scores, both rankings, the zero overlap, mean and min T2, the
crossing counts, and the 2.15x verdict. What moved: every absolute rate (halved),
and, slightly, the protection factors, which are a ratio of two rates taken at
**fixed H** and so are invariant only in the limit γ → 0. Two moved at the
printed precision (2.52x → 2.51x, and the top-5 mean 2.54x → 2.53x); the mean-T2
chains sit close enough to that limit to show no movement.

---

## Time stability

The best concentrator chain [85, 15, 86, 16, 87] tracked across 5 months:

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
sits at the edge, creating a built-in concentrator. The noisy
qubit absorbs disproportionate damping, and the cavity modes
localized on the interior survive longer.

This requires no additional gates, no error correction, and no
knowledge of the palindromic theory. It is a free improvement
available on any quantum processor with non-uniform noise
characteristics. The only input is the coupling map and the T2
calibration data, both publicly available.

The theory predicts which chains will perform best. The prediction
is testable with a single set of Trotter evolution experiments
comparing concentrator chains against mean-T2 chains on the
same hardware on the same day.

---

*See also:*
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) (2.81x theoretical),
[IBM Concentrator](IBM_CONCENTRATOR.md) (1.97x measured),
[Resonant Return](RESONANT_RETURN.md) (the concentrator formula)
