# Concentrator Qubit Mapping: Finding Optimal Chains on Real Hardware

**Naming note (2026-07-05):** renamed from "Sacrifice-Zone Qubit Mapping". The
noisy edge qubit sacrifices nothing; it concentrates the noise (the misnomer
was resolved 2026-03-28). The frozen `sacrifice_zone_mapping.*` script and data
keep their original names, which is where this work started: the engineering
rule was "protect the receiver side, sacrifice the sender side", and the noun
outlived the rule. The names are the lineage, not a freeze; the script was last
changed on 2026-09-06, when its coupling map was corrected.

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
**Topology:** IBM Torino's own coupling map via `FakeTorino` (133 qubits, 150 edges), the heavy-hex family where each node connects to 2 or 3 neighbours in a hexagonal pattern with extra "bridge" qubits on each edge

---

## What this document is about

This document shows how to find optimal qubit chains on real IBM hardware
by exploiting naturally noisy qubits as concentrators. Instead of
picking the qubits with the best T2 times (the standard approach), we
select chains where a noisy qubit sits at the edge, providing the
concentrator benefit for free. On IBM Torino's 133-qubit chip, this
mode-based selection outperforms naive T2 maximization by 2.6× in
protection factor, despite using qubits with 1.9× lower average T2.

---

## Abstract

If the concentrator protects cavity modes localized on interior
qubits (r = 0.994, see [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)),
then choosing qubit chains where a naturally noisy qubit sits at the
edge should provide the concentrator benefit *for free*.

We test this on IBM Torino's heavy-hex topology using real T2
calibration data (181 days, 24,073 records). On the latest calibration
date 133 qubits carry T2 data, and we search the device's own coupling map,
150 edges over those 133 qubits. 359 five-qubit chains exist on it. We compare
two chain selection strategies:

1. **Concentrator ranking:** Maximize edge noise / interior noise ratio
2. **Mean-T2 ranking:** Maximize average T2 across all 5 qubits

Result: **Zero overlap** in the top-10 lists. Concentrator chains
achieve **3.12x** mean protection factor vs **1.19x** for mean-T2
chains. Mode-based selection outperforms naive T2 maximization by 2.6x.

The best concentrator chain has only 138 us mean T2 but 3.42x protection.
The best T2 chain has 234 us mean T2 but only 1.12x protection.
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
T₁ share of the decay is wildly uneven across this chip. On Q85, the concentrator
qubit in all five headline chains, T₁ = 2.9 µs against T₂ = 5.0 µs, so
the 1/(2T₁) term carries most of the coherence decay, while the interior qubits
sit near 40%. Under that model Q85's D[Z] rate would fall several-fold more than
the interior's and the ranking would not survive. A σ⁻ channel also breaks the Π
operator (F82/F84), and beside the co-axial Z-dephasing these chains carry it
breaks the spectral palindrome too, so the column below would not survive
either. (σ⁻ *alone* would not break it: it would move the centre to
−Σγ/2, [F137](../docs/ANALYTICAL_FORMULAS.md#f137). It is the co-axial pair
that breaks.) Everything below is a statement about the dephasing-only model
these scripts implement.

Related, and unhandled by either script: Q53 reports T₂ = 62.4 > 2·T₁ = 44.8 on
this calibration date. That is a broken record (the typed layer clamps it,
`IbmCalibration.cs`), and the dephasing-only form is silent about it only
because it never reads T₁. Q53 sits in none of the ten chains the current search
returns; the broken record is a property of the calibration file rather than of
any chain, which is why it is kept here.

### Spectral verification

For the top-5 chains in each ranking, we compute the full Liouvillian
(N=5 chain, real γ values) and extract: palindrome backward error, slowest
oscillating mode rate, and protection factor vs uniform noise.

---

## Results

### Concentrator ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [85, 86, 87, 88, 89] | 32.2 | 138.1 us | 3.42x |
| [66, 67, 74, 86, 85] | 29.7 | 124.0 us | 3.11x |
| [85, 86, 87, 88, 94] | 29.0 | 130.2 us | **3.44x** |
| [68, 67, 74, 86, 85] | 22.3 | 98.8 us | 3.16x |
| [81, 82, 83, 84, 85] | 17.2 | 75.2 us | 2.49x |

All contain Q85 (T2 = 5.0 us), the noisiest qubit on the chip,
as the concentrator endpoint.

### Mean-T2 ranking (top 5)

| Chain | Score | mean T2 | Protection |
|-------|-------|---------|-----------|
| [10, 11, 12, 18, 31] | 1.1 | 234.5 us | 1.12x |
| [8, 9, 10, 11, 12] | 0.9 | 217.9 us | 1.31x |
| [9, 10, 11, 12, 18] | 2.0 | 216.5 us | 1.21x |
| [9, 10, 11, 12, 13] | 2.0 | 216.3 us | 1.21x |
| [10, 11, 12, 13, 14] | 2.1 | 214.6 us | 1.12x |

All have concentrator scores between 1 and 2 (nearly uniform noise). The quiet
qubits provide long T2 but no differential protection. None of these five
contains Q85: a chain selected for mean T2 does not run through the chip's
noisiest qubit, which is the whole of the contrast below.

### Head-to-head

| Metric | Concentrator top-5 | Mean-T2 top-5 |
|--------|----------------|---------------|
| Mean protection factor | **3.12x** | 1.19x |
| Mean T2 | 113.3 us | 220.0 us |
| Mean concentrator score | 26.1 | 1.6 |
| Palindrome backward error (see below) | 52.7-72.1 ε | 52.1-68.4 ε |

The concentrator chains have 1.9x lower mean T2 but 2.6x higher
protection. Choosing "worse" qubits with the right spatial pattern
outperforms choosing the "best" qubits naively.

**The palindrome column was void and has been replaced (2026-08-05).** It used
to report 96-98% for concentrator chains and 85-92% for mean-T2 chains, and
neither number measured the palindrome. Measured on a 5-site chain carrying the
γ of qubits [18, 89, 19, 90, 60] (those five were the mean-T2 top row of an
earlier search and are not adjacent on Torino, so this is a Liouvillian with
real measured rates rather than a layout; the palindrome does not know about
layout), the palindromic symmetry holds **to the eigensolver's own accuracy**:
comparing the
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

Measured: **52.1 to 72.1 ε** across the ten chains, a **1.384×** band. That band
is at its floor, not merely small: **one sensitivity of the measurement alone
exceeds it, and a second is of the same size**.

Permuting only the ORDER in which the five jump operators are summed into L, at
identical physics and identical γ, moves this same number across all 120 orders
from 51.2 to 90.9, a **1.77×** spread. (Measured on one 5-site chain carrying the T₂-derived γ of qubits
[18, 89, 19, 90, 60]; all 120 orders, identical under both in-place and
out-of-place accumulation. Those five qubits do not form a path on Torino, so
the run is a 5-site Heisenberg Liouvillian carrying real measured rates rather
than a realizable chain; nothing in this paragraph depends on the layout.) Ten physically different chains vary
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

**The defective scorer is a family of at least seventeen. Eight of them are
converted; the nine listed below still stand.** Named here rather than left quiet, because a substitution that
converts some sites and not others is worse than consistent wrongness: the
inconsistency reads as deliberate. The count grew from "two" to "five" to this
over three review rounds on 2026-08-05, so treat it as a floor and not as a
census; it was produced by grepping two code shapes, and a third shape would not
have been seen.

**Converted:** `sacrifice_zone_mapping.py` and `chain_selection_test.py` (this
document and [Chain Selection Test](CHAIN_SELECTION_TEST.md));
`ibm_cavity_analysis.py`, whose 100% / 100% / 100% row in
[IBM Cavity Spectral Analysis](IBM_CAVITY_SPECTRAL_ANALYSIS.md) was the
saturation artefact described above and now reads 67.5 / 72.4 / 57.5 ε·ρ;
`combined_optimization.py`, whose 89% / 85% / 92% / 96% / 94% / 95% row now
reads 60.8 / 66.9 / 78.0 / 72.0 / 62.2 / 68.0 ε·ρ, all six at the floor; and
`optimal_chain_search.py`, which is the one entry that turned out to publish
**nothing**. Its score was computed and never printed: deleting the fifteen
lines leaves the results file bit-identical, which is both the repair and the
proof it was dead. That is a caution about this inventory, not just about that
script, since the list was built by grepping for the code shape and a site's
presence here does not establish that a table depends on it. `v_effect_thermal.py`
is the sixth, under Shape A below; the two C# sites are the seventh and eighth,
and what they turned out to be is set out further down.

The five that still compute the check call one implementation;
`optimal_chain_search.py` has no check to call. The port had reached three hand-copies,
which is the cockpit signal that a primitive was missing, so it lives in
`simulations/framework/symmetry.py` as `fw.max_f1_pairing_distance` /
`fw.f1_distance_in_eps`, with its blind spots pinned in
`framework/tests/primitives/test_f1_pairing_distance.py`. The remaining Python
sites should import it rather than copy it again.

**Shape A** is empty. Its last member, `v_effect_thermal.py` (1e-3, feeding
[Thermal Breaking](THERMAL_BREAKING.md)), was converted on 2026-08-05, and it
carried a second defect the tolerance had hidden: three of its four call sites
scored a spectrum with an amplitude-damping channel against the *Z-dephasing*
centre. Converting it turned a symmetry that reads as decaying with temperature
(100% / 93% / 93% / 98%) into one that is exact at every temperature, which is
the largest thing this scorer family has hidden so far.

**Shape B**, nearest-partner **without removal** (so several rates may claim the
same partner and it is not multiplicity-aware), with a `999` sentinel and only
the below-centre half scored: `analytical_spectrum_verify.py`,
`deep_band_structure.py`, `deep_computation.py`, `frequency_test.py`,
`mirror_symmetry_deep.py`, `mirror_transition.py`, `n5_optimal_cavity_size.py`,
`nested_mirror_asymptote.py`, `overnight_computation.py`.

**The two C# sites were dealt with on 2026-08-06, and only one of them was what
this list said it was.**

`MirrorAnalysis.CheckSymmetry` (`compute/RCPsiSquared.Compute/MirrorAnalysis.cs`)
was Shape B at a tolerance of 0.005 and genuinely an F1 scorer: its reflection
d ↦ 2σ − d is F1's real part, centred on σ with partners summing to 2σ.

The embedding that carries rates into the canonical complex check, λ = −d, for
which −2σ − λ = −(2σ − d) is exactly the mirrored rate, did not stay in that
class. It lives in Core as `F1SpectrumStatistics.MaxF1RatePairingDistance`,
beside the canon it delegates to, and `MirrorAnalysis` is a thin caller adding
only the eps ratio and its empty-input policy. Six tests in
`Core.Tests/F1/F1SpectrumStatisticsTests.cs` pin it, and four of them were
checked against mutants rather than trusted: flipping the embedding sign kills
four, disabling the matcher's removal kills three. The producer
[`f1_rate_embedding_check.py`](../simulations/f1_rate_embedding_check.py) runs
the same identity on chain and star at N=2..5 and reports EXACT equality on all
eight rows, as the algebra requires, since d ↦ −d is an isometry and preserves
the index order the greedy matcher walks.

Two things fell out of doing it. The canonical complex check was itself
returning 0.0, a perfect score, for an EMPTY spectrum, which is the same quiet
failure the rate sibling's guard had just been written to forbid; both throw
now. And what the projection costs is not a number: the producer measures 1.74
to 9.56 across its eight rows, while the live witness at N=6, γ=0.5 reads 1.0.
The projection's weakness is structural, not measured, so no factor should be
quoted for it at all. `BlockSpectrumWitness` carries the two readings side by
side now (`inspect --root blockspectrum`), with the reason it must not be read
as a ratio. `MirrorAnalysis.Analyze` went with it,
as dead code; its output string had never been emitted anywhere in the repo.

`FillingThresholdCsr.ConjugationMatchFraction`
(`compute/RCPsiSquared.Diagnostics/`) is **not** an F1 scorer, and the earlier
instruction to point it at `MaxF1PairingDistance` was wrong. Its involution is
λ ↦ λ*, complex conjugation, not λ ↦ −2σ − λ; substituting the F1 check there
would have computed a different mathematical object. What it shares with the
family is only the matcher, and that is what was repaired: the match is
multiset matching with removal now. The repair is inert, which is the honest
reading and was measured rather than assumed. Against the pre-repair first-fit
matcher the fraction is unchanged to the printed digit on all three
configurations probed: 100.00% on the clean block at Δ=0, 1.00% clean at Δ=1,
and 0.00% under a random field. Only the last of those three is a configuration
the code actually runs, since both call sites pass a disorder field; the other
two were built for the probe. The defect never fired on the live input, and it
could not have: without removal the fraction can only be inflated, while both
consumers, the `FillingThresholdWitness` conjugation-match diagnostic used
alongside its GinUE comparison and the `frac < 0.1` test, read it against ≈ 0.
That diagnostic does not assign a symmetry class.

The nine remaining Python sites should import `fw.max_f1_pairing_distance`,
which is that same check ported once.

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

**Which of these numbers depend on the γ book.** The dephasing rate here is the
D[Z] rate γ = 1/(2T₂), not 1/T₂ (see [the glossary section](../docs/GLOSSARY.md)).
Choosing the other book is not a change of units: J is held at 1, so it moves Q =
J/γ and with it the physical point. The quantities of γ-degree zero do not care
which book is used: the concentrator scores, both rankings, the zero overlap,
mean and min T2, and the crossing counts. The protection factors do care a
little, being a ratio of two rates taken at **fixed H**, hence book-independent
only in the limit γ → 0; the mean-T2 chains sit close enough to that limit to
show no movement at the printed precision, the concentrator chains move in the
last digit.

---

## Time stability

The best concentrator chain [85, 86, 87, 88, 89] tracked across 5 months:

| Date | Score | mean T2 |
|------|-------|---------|
| 2026-02-10 | 33.75 | 138.1 us |
| 2025-12-12 | 33.70 | 166.1 us |
| 2025-10-13 | 30.09 | 159.5 us |

The score varies by ~11% but the chain consistently ranks at the top.
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
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) (2.80x theoretical),
[IBM Concentrator](IBM_CONCENTRATOR.md) (1.97x measured),
[Resonant Return](RESONANT_RETURN.md) (the concentrator formula)
