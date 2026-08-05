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
**Read on from here:** [Chain Selection Test](CHAIN_SELECTION_TEST.md), which
takes the two chains ranked here to a head-to-head and shares this run's
palindrome check, hence the same backward-error column
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
(N=5 chain, real γ values) and extract: palindrome backward error, slowest
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
| Palindrome backward error (see below) | 53.3-73.1 ε | 60.2-76.5 ε |

The concentrator chains have 2.3x lower mean T2 but 2.15x higher
protection. Choosing "worse" qubits with the right spatial pattern
outperforms choosing the "best" qubits naively.

**The palindrome column was void and has been replaced (2026-08-05).** It used
to report 96-98% for concentrator chains and 85-92% for mean-T2 chains, and
neither number measured the palindrome. Measured on the mean-T2 top chain
[18, 89, 19, 90, 60], the
palindromic symmetry holds **to the eigensolver's own accuracy**: comparing the
960 oscillatory rates against their mirror image 2·Σγ − rate as sorted multisets
gives a residual of **1.8e-14**. The theorem is proven analytically; what this
measures is that nothing in the numerics contradicts it at that scale.

What the percentage measures is the matcher. `spectral_analysis` pairs each rate
with the **first** partner it finds within an **absolute** 1e-4, and the spectrum
is far denser than that: **927 of the 959** nearest-neighbour gaps are below
1e-4. So a rate almost never grabs its true partner, and it orphans that
partner's mate. Of the 410 PAIRS the scorer accepts on that chain, **362 are
wrong by more than 1e-12**. (410 pairs cover 820 of the 960 rates, and
820/960 = 85.4%: that is where the published percentage comes from.) The column measures how a greedy first-fit scrambles
inside the level clustering, which is why it can move in either direction: across
the ten chains the 2026-08-05 γ-book repair moved it by 0, 1, 2, 3 and 4 points,
and **upward** on two of them.

The repair was therefore **not** a retuned tolerance, and the column now reports
the residual itself. Any tolerance below the minimum level spacing returns
exactly 100% and hides the same clustering, and 1e-6, 1e-8 and a scale-relative
1e-4·Σγ all do, so no threshold could have been the answer.

**What the column reports now, and whose it is.** The check is
`F1SpectrumStatistics.MaxF1PairingDistance`
(`compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs`), which the C# summary
calls the canonical F1 check, which stands behind
`MultisetAssert.NearestNeighbourEqual`, and which is already a live witness as
`BlockSpectrumWitness.PalindromePairingDistance` with committed values in
[`f1_n8_n9_metrics/`](../simulations/results/f1_n8_n9_metrics/). **The repo owned
this metric and these scripts were not using it**; they are now, via a port
carrying that provenance in its docstring. It is the max greedy
nearest-neighbour distance, WITH REMOVAL, between the eigenvalue multiset and its
F1 reflection {−2σ − λ}, on the **full complex** spectrum, so it is
multiplicity-aware and a dropped or duplicated eigenvalue cannot hide in it.

A rates-only version stood here for part of 2026-08-05: comparing only −Re(λ) as
a sorted multiset, which for a sorted list is exactly r[k] + r[n−1−k] = 2·Σγ.
That identity is correct and is the optimal pairing in one dimension, but it is a
strictly **weaker** test than F1, because it discards the imaginary parts: halving
one eigenvalue's Im leaves it bit-identical. It is recorded here because inventing
a metric next to one the repo already owns is the failure this column exists to
illustrate, twice in one day.

**The error model, since there is no exact route** (a non-Hermitian eigensolver).
Per the repo's no-rounding rule the number is published with its model rather than
against a threshold: the standard eigenvalue backward error is O(ε·‖L‖), so the
column is the distance in units of ε · spectral radius. Normalising instead by
max|rate| spreads the same ten chains by **15.0×**, because max|rate| is only the
real part while the spectrum is dominated by |Im| (ρ / max|rate| is 47× to 600×
here).

Measured: **53.3 to 76.5 ε** across the ten chains, a **1.44×** band (1.4353). That band
is at its floor, not merely small: **one sensitivity of the measurement alone
exceeds it, and a second is of the same size**.

Permuting only the ORDER in which the five jump operators are summed into L, at
identical physics and identical γ, moves this same number across all 120 orders
from 51.2 to 90.9, a **1.77×** spread. (Measured on the mean-T2 top chain
[18, 89, 19, 90, 60]; one chain, all 120 orders, and identical under both
in-place and out-of-place accumulation.) Ten physically different chains vary
less than one re-associated sum on a single chain. That is this repo's documented case-3 residual, a
deterministic function of an input the physics does not contain.

The second is a property of the matcher itself, and it is stated rather than left
implicit because this column has now been caught twice on measuring its matcher:
greedy nearest-neighbour is **order-dependent by construction**. Permuting only
the ARRAY ORDER of the eigenvalues, same L and same spectrum, gives 58.9 against
75.8 on one chain, a **1.29×** spread with no physics in it at all. So even the
canonical check partly measures its own matcher. It is used anyway: it is the
repo's canonical check, it is multiplicity-aware where the alternatives are not,
and both sensitivities push the same conclusion.

The reading is therefore the same for every chain: the palindrome holds to the
eigensolver's own accuracy, and **no chain ranks above another on this column**.
Ranking it would be ranking arithmetic order, which is why the sibling table in
[Chain Selection Test](CHAIN_SELECTION_TEST.md) prints "tie" there by
construction. What the number does **not** license is a claim of exactness: a
genuine violation below roughly 1e-13 absolute would sit inside the same band. F1
is proven analytically; this column is the numerical check, not the theorem.

**The defective scorer is a family of at least seventeen. This repair reached
two of them; the fifteen listed below still stand.** Named here rather than left quiet, because a substitution that
converts some sites and not others is worse than consistent wrongness: the
inconsistency reads as deliberate. The count grew from "two" to "five" to this
over three review rounds on 2026-08-05, so treat it as a floor and not as a
census; it was produced by grepping two code shapes, and a third shape would not
have been seen.

**Fixed here:** `sacrifice_zone_mapping.py` and `chain_selection_test.py` (this
document and [Chain Selection Test](CHAIN_SELECTION_TEST.md)).

**Shape A**, greedy first-fit inside a tolerance, reported as a percentage,
`simulations/`: `ibm_cavity_analysis.py` (1e-6, and its 100% / 100% / 100% row in
[IBM Cavity Spectral Analysis](IBM_CAVITY_SPECTRAL_ANALYSIS.md) is the saturation
artefact described above), `optimal_chain_search.py` (1e-4, the old code
verbatim, same function name), `combined_optimization.py`
(`max(1e-4, 1e-3·center)`, precisely the scale-relative retune ruled out above,
publishing 89%, 85%, 92%, 96%, 94%, 95% for **these same chains**),
`v_effect_thermal.py` (1e-3, feeding [Thermal Breaking](THERMAL_BREAKING.md),
whose 91% already carries a footnote calling it a tolerance artefact: the right
instinct attached to the wrong mechanism, since the tolerance saturates rather
than degrades).

**Shape B**, nearest-partner **without removal** (so several rates may claim the
same partner and it is not multiplicity-aware), with a `999` sentinel and only
the below-centre half scored: `analytical_spectrum_verify.py`,
`deep_band_structure.py`, `deep_computation.py`, `frequency_test.py`,
`mirror_symmetry_deep.py`, `mirror_transition.py`, `n5_optimal_cavity_size.py`,
`nested_mirror_asymptote.py`, `overnight_computation.py`.

**In C#**, and this is where it matters most: `MirrorAnalysis.CheckSymmetry`
(`compute/RCPsiSquared.Compute/MirrorAnalysis.cs`) is Shape B at a tolerance of
0.005, and `FillingThresholdCsr.ConjugationMatchFraction`
(`compute/RCPsiSquared.Diagnostics/`) is a first-fit variant sitting in the
**live** Diagnostics layer with a witness on top of it.

All of them should call `F1SpectrumStatistics.MaxF1PairingDistance`, which the
same solution already contains.

**The tolerance table below looks impossible and is not, and the reason it is
not is the defect itself.** Tightening an acceptance window should admit fewer
pairs; here it admits more. Measured on this chain:

| tolerance | pairs accepted | rates paired | rates orphaned | score |
|:---|---:|---:|---:|---:|
| 1e-4 (the retired default) | 410 | 820 | **140** | 85.4% |
| 1e-6 | 480 | 960 | 0 | 100% |

The cause is the greedy first-fit. A loose window lets a rate seize a partner
that is not its mirror; that partner is then consumed, and its own true mate can
be left with nothing available. Tighten the window and every rate finds *some*
partner, so the score reaches 100% by construction, which is exactly why 100% is
not evidence of anything either.

**It is not, however, true that a tight window makes every rate find its own**,
and that sentence stood here until 2026-08-05. Measured on this chain, at
tolerance 1e-6 the scorer accepts 480 pairs and prints 100%, yet **56 of those
480 are still wrong by more than 1e-12**; at 1e-8, 8 of 480 still are. The
mispairing rate falls from 71% to 12% to 2% while the printed score sits at 100%
throughout. Saturation, not correctness, is what the tighter window buys.

The two counts are **not** derivable from one another, and it is worth saying so
because the obvious guess is wrong: 362 mispairings do not imply 724 orphans.
Only 140 rates end unpaired, because a mispairing can itself join two rates that
would otherwise both have been orphaned. Both numbers are measured, not counted
from each other. These percentages carry no information, ordinal or otherwise,
and are kept only as the record of what the retired scorer produced.

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
