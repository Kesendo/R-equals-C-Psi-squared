# IBM Cavity Spectral Analysis: Why the Sacrifice Zone Works

<!-- Keywords: sacrifice zone cavity mode protection, IBM Torino spectral
analysis palindromic eigenvalues, Liouvillian zero noise cavity modes
hardware, spatial dephasing profile mode survival, Clebsch-Gordan
formula hardware verification, R=CPsi2 IBM cavity spectral -->

**Status:** Tier 2 (computed spectral analysis of real IBM hardware data)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md),
[IBM Concentrator](IBM_CONCENTRATOR.md)
**Script:** [ibm_cavity_analysis.py](../simulations/ibm_cavity_analysis.py)
**Data:** [ibm_cavity_analysis.txt](../simulations/results/ibm_cavity_analysis.txt)
**Raw IBM data:** [data/ibm_sacrifice_zone_march2026/](../data/ibm_sacrifice_zone_march2026/) (4 JSON files, March 24, 2026, IBM Torino)

---

## What this document is about

This document explains WHY the sacrifice-zone formula works, using the
full Liouvillian spectrum computed from real IBM Torino hardware data.
The answer: the formula does not protect individual qubits; it protects
cavity modes (the collective standing-wave oscillations of the entire
chain). Concentrating noise on one edge qubit creates a spatial gradient
that shields modes localized on the quiet qubits, giving the slowest
oscillating modes a 2.80× lifetime improvement.

---

## Abstract

The sacrifice-zone formula concentrates noise on one edge qubit while
protecting the rest. On IBM Torino hardware, this produces 1.97x
improvement in coherence survival. But WHY does concentrating noise help?

We answer this by computing the full Liouvillian spectrum of the N=5
chain using real IBM T2* data (Q85-Q94). The same 43 cavity mode
frequencies exist under all noise profiles. Only the damping changes.
The sacrifice zone does not protect qubits. It protects **cavity modes**.

Key results:
- Palindromic to the eigensolver's floor under strongly asymmetric IBM
  gammas (F1 pairing distance 72.4 eps x spectral radius; Q85 has 26x
  more dephasing than Q87)
- Slowest oscillating mode: 2.80x longer lifetime under sacrifice vs
  uniform (0.0229 vs 0.0640)
- 12 protected modes (rate < 0.05) under sacrifice, 0 under uniform
- Max decay rate = 2 x Σγ, exact by the palindrome rather than by the
  digits: the stationary mode at λ = 0 has an F1 partner at −2Σγ

---

## The connection

At zero noise (Σγ = 0), the Liouvillian has 120 stationary
modes and 904 oscillating modes across 43 distinct frequencies. These
are the **eigenfrequencies of the resonator** (see
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md)). The stationary count
matches the Clebsch-Gordan formula (the angular-momentum addition rule that counts how many spin multiplets exist for N coupled qubits) exactly (120 = Sum_J m(J,5)*(2J+1)^2).

When noise is turned on, the frequency comb persists and what changes
most is the damping: each mode acquires a decay rate that depends on how
much noise it "sees." The frequencies themselves are not fixed, and this
page runs a deliberately non-uniform profile, which is the case that moves
them (see [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md); the grouping
below bins at freq_tol = 0.1, and the worst shift it has to swallow under
this profile is 3.2e-2, so the margin is about threefold, not the order of
magnitude the median shift would suggest). Modes localized on quiet qubits (Q86-Q94) see less
noise. Modes touching the sacrifice qubit (Q85) see more.

The sacrifice zone works because it creates a **spatial gradient** in
the noise. Modes localized away from the sacrifice qubit are shielded.
The protection is not about qubits. It is about modes.

---

## IBM hardware data

Qubits: Q85 (sacrifice), Q86, Q87, Q88, Q94 (chain topology)

| Qubit | T2* (us) | 1/T2* (1/us) | γ = ½ · (1/T2*) (1/us) | Role |
|-------|---------|-------------|-------------------|------|
| Q85 | 3.73 | 0.2681 | 0.13405 | Sacrifice (26x the total rate of Q87) |
| Q86 | 61.35 | 0.0163 | 0.00815 | Interior |
| Q87 | 97.54 | 0.0103 | 0.00515 | Interior (quietest) |
| Q88 | 67.99 | 0.0147 | 0.00735 | Interior |
| Q94 | 95.03 | 0.0105 | 0.00525 | Edge (quiet) |

Σγ = 0.15995. Q85 carries 84% of the total budget of this column.

The γ column is half of the printed 1/T2* column, not half of 1/T2* recomputed
from the T2* column: the pipeline rounds the rate first and the script is fed
those rounded values, so Q87 reads 0.00515 where 1/(2·97.541) would give
0.00513. The digits above are the ones that actually enter the Liouvillian.

**Which column enters the Liouvillian, and why it is the γ column.** The
calibration reports a coherence rate 1/T2*. A `D[Z]` channel at rate γ decays
coherences at 2γ, so for a Z-dephasing-only model the D[Z] coefficient is
γ = 1/(2·T2*): row 1 of the T2 → γ table in
[GLOSSARY](../docs/GLOSSARY.md). This model builds `D[Z]` and nothing else, so
row 1 is what row 1 is for and no further argument is needed. Every rate on this
page is computed from the γ column.

Row 2 of that table, the T1-aware γ_Z = (1/T2 − 1/(2T1))/2, is not an
alternative here, and the reason is narrower than it first looks. Row 2's
coefficient reproduces the measured T2 only when a σ⁻ channel is built beside
it; drop the σ⁻ and the model no longer matches the calibration it came from.
Build the σ⁻ and the page changes subject, because σ⁻ beside co-axial
Z-dephasing is what breaks the SPECTRAL palindrome (GLOSSARY, the paragraph
after the table). Note the distinction that paragraph exists to protect: σ⁻
breaks Π, and Π breaking is not the palindrome breaking. By F137 a σ⁻ channel
alone leaves the spectrum palindromic at the halved centre. So the reason to
stay on row 1 is the model, not the symmetry.

**What that choice moves, and what it does not.** Two things scale exactly with
γ, and both are printed below: the palindrome centre Σγ, and the maximum decay
rate 2Σγ, which is exact by the theorem rather than by the digits. The ratios
among the γ's themselves are scale-free too. What is NOT scale-free is anything
that puts an interior rate against another: H is held at J = 1 while γ scales,
so the Liouvillian is affine in γ rather than homogeneous, the fixed H competes
with the scaled γ, and interior ratios drift (the protection factor reads 2.80x
here and would read 2.81x on the doubled column). The protected-mode count moves
much more than that, and not because of the physics: the 0.05 threshold is
absolute and does not scale with the profile. It was 16% of Σγ on the doubled
column and is 31% here, which is why the count reads 12 rather than 4. That count
is a statement about a fixed 0.05, not a property of the profile, and it should
not be read as three times the protection.

**One idealization the γ column carries, and what it costs the headline.** On
Q85 the same calibration reports T1 = 2.836 us, so 1/(2·T1) is about 66% of
1/T2*: the dephasing-only model books as dephasing a coherence loss that is
mostly relaxation. The consequence lands on the asymmetry figure this page
quotes. The 26x and the 84% are ratios of TOTAL coherence rates, and they are
untouched by the halving, which is why they are stated unchanged. They are not
the ratio of pure DEPHASING. Taking the same calibration's T1 into account
(row 2's γ_Z, computed only to price this, not to build a model) gives
Q85/Q87 = 10.9x and a Q85 share of 68%. So the sacrifice zone's dephasing
asymmetry is roughly 11x; the 26x is the asymmetry of the total watching, most
of Q85's share of it being T1. Read every "26x" on this page in that sense.

### Source data (all in repository)

| File | Content |
|------|---------|
| [`sacrifice_zone_hardware_20260324_191713.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hardware_20260324_191713.json) | Calibration data, T1/T2/T2* values, chain topology, experimental parameters |
| [`sacrifice_zone_hw_selective_dd_20260324_191523.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_selective_dd_20260324_191523.json) | Selective DD bitstring counts (4000 shots per time point) |
| [`sacrifice_zone_hw_uniform_dd_20260324_191614.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_uniform_dd_20260324_191614.json) | Uniform DD bitstring counts |
| [`sacrifice_zone_hw_no_dd_20260324_191713.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_no_dd_20260324_191713.json) | No DD bitstring counts |

Backend: IBM Torino (ibm_torino), March 24, 2026. J_coupling = 1.0.
Trotter steps: [2, 4, 6, 8, 10] at dt = 0.5 us.

---

## Results

### Three profiles compared

| Property | Zero noise | IBM sacrifice | Uniform |
|----------|-----------|--------------|---------|
| Stationary modes | 120 | 6 | 6 |
| Oscillating modes | 904 | 1018 | 1018 |
| Distinct frequencies | 43 | 120 | 112 |
| F1 pairing distance (ε · ρ), do not rank | 67.5 | 72.4 | 57.5 |
| Palindrome center | 0.0000 | 0.1599 | 0.1599 |
| Min decay rate (osc.) | 0 | 0.0229 | 0.0640 |
| Max decay rate (osc.) | 0 | 0.2971 | 0.2559 |
| Max decay rate (all) | 0 | 0.3199 | 0.3199 |
| Protected (rate < 0.05) | 904 | 12 | 0 |

Note: Max decay rate (all) = 2 x Σγ = 0.3199 applies to
non-oscillating modes (freq = 0) that represent pure decay. The
oscillating modes have lower maximum rates (0.2971 sacrifice,
0.2559 uniform).

Note on the F1 row: the "do not rank" in the label is there because the three
entries are not comparable quantities. Each is one draw from a permutation
orbit, and two of them happen to be the orbit's ceiling while the middle one is
not. See [the reading of that number](#the-palindrome-under-asymmetric-noise)
below.

### The 2.80x protection factor

The four slowest oscillating modes under the IBM sacrifice profile all
have frequency 7.2355, next to the zero-noise 7.2361, and decay rate 0.0229.
Under uniform noise, the slowest modes decay at 0.0640. Ratio: **2.80x**.

IBM hardware measured 1.97x improvement at early times. The computed
2.80x is the theoretical maximum. The hardware measurement is lower
because gate errors, crosstalk, and finite-time effects reduce the
effective protection.

### Mode survival comparison

| Rank | IBM rate | Uniform rate | Ratio | Frequency |
|------|---------|-------------|-------|-----------|
| 1 | 0.0229 | 0.0640 | 2.80x | 7.236 |
| 2 | 0.0229 | 0.0640 | 2.80x | 7.236 |
| 3 | 0.0229 | 0.0640 | 2.80x | 7.236 |
| 4 | 0.0229 | 0.0640 | 2.80x | 7.236 |
| 5 | 0.0467 | 0.0640 | 1.37x | 0.000 |
| 6 | 0.0467 | 0.0640 | 1.37x | 0.000 |
| 7 | 0.0491 | 0.0640 | 1.30x | 0.000 |
| 8 | 0.0491 | 0.0640 | 1.30x | 0.000 |
| 9 | 0.0496 | 0.0640 | 1.29x | 5.234 |
| 10 | 0.0496 | 0.0640 | 1.29x | 5.234 |

The top 4 modes (freq 7.236) are all 2.80x protected. These are modes
that oscillate at approximately 7.2 times the coupling strength. They
are the modes that the sacrifice zone was built to protect.

Modes 5-8 (non-oscillating, freq = 0) have only 1.3-1.4x protection.
These are pure decay modes. The sacrifice zone preferentially protects
the oscillating modes over the decaying ones.

---

## The palindrome under asymmetric noise

The IBM sacrifice profile has Q85 at 26x more dephasing than Q87.
Despite this extreme asymmetry, the palindrome survives to the
eigensolver's floor: every eigenvalue finds a partner summing to
-2 x Σγ = -0.3199, with a worst pairing distance of 72.4 eps x spectral
radius. This confirms the analytical proof: the palindrome depends on
the SUM of gammas, not their distribution.

**How that number should be read.** It is the repo's canonical F1
check, `F1SpectrumStatistics.MaxF1PairingDistance` (greedy
nearest-neighbour WITH removal, on the full complex spectrum, so it is
multiplicity-aware and not blind to the imaginary parts), reported
against its backward-error model: an eigensolver on a non-normal matrix
has no exact route, so the honest unit is eps x spectral radius rather
than a threshold. At N=5 the floor is O(10-100); the band is
N-dependent, so that range does not transfer (see below). It does **not**
license the word "exact": a genuine violation below about 1e-13
absolute would sit inside the same band. F1 is proven analytically; this
is the numerical check, not the theorem.

The three profiles must not be ranked by it, and the reason needs no appeal to
noise: the printed cells are not the same kind of number. The F1 distance comes
from a GREEDY matcher, so its value depends on the order the eigenvalues reach
it, and each cell is one draw from an ORBIT of values that one spectrum can
produce. Permuting the eigenvalue array never raises the result above the
matcher's best available pairing, so an orbit has a ceiling. Over 200 array
permutations at each of three RNG seeds, the IBM cell 72.4 and the uniform cell 57.5
are never once exceeded: both printed values ARE their ceiling. The zero-noise
cell 67.5 is not. Roughly three quarters of orderings beat it, its ceiling is
83.5 and its median between 70 and 73. Two ceilings and one interior draw are not a
ranking, and that is a defect in the row rather than a property of the systems.

The second orbit settles it. Permuting the summation order of the jump
operators moves the IBM number from 50.2 to 108.7 over all 120 orders, 76
distinct values, and the printed 72.4 sits in the MIDDLE of that range. So the
same cell is the ceiling of one orbit and a mid-orbit draw of another, and the
two orbits disagree about where it stands. No summary of a single number
survives that. Note also which nuisance is the larger: the IBM column moves by
1.45x to 1.48x under array order, depending on the RNG seed, and by 2.17x
under jump order, which draws no random numbers at all because it is a complete
enumeration; the array-order sweep
is the narrower exposure, not the more general one. The jump-order sweep speaks
only for the IBM column: on uniform every γ is equal and on zero noise every
jump operator is the zero matrix, so all 120 orders assemble a bit-identical L
(max absolute difference exactly 0.0 in both, against 1.7e-16 on IBM) and the
number does not move at all: one distinct value each, against 76 on IBM.

Every figure in these three paragraphs is emitted by
`python simulations/ibm_cavity_analysis.py --sweeps`, into the same results
file as the tables above. Read them there rather than reconstructing the
harness: the array-order counts are sampled (200 permutations at each of the
RNG seeds 0, 12345 and 7) and the paired test quotes RNG seed 0. "Seed" here
is numpy's, not the repo's [Seed](../compute/MirrorWorld/README.md) object.

What is NOT true, and an earlier version of this page claimed it, is that the
difference between the columns is bookkeeping noise that a larger sample would
absorb. It is not. Applying the SAME permutation to all three spectra, which is
how a shared nuisance has to be tested, the IBM entry exceeds the uniform entry
in 192 of 200 draws by a mean of 14.4, about 33 standard errors: the
order-dependence is almost entirely common mode and cancels in the difference.
The columns do differ systematically. What that difference is not is a
difference in palindromicity. All three sit at the eigensolver's floor, and what
separates them is how well-conditioned each eigenproblem is, the IBM profile
being the least normal of the three (the commutator norms below). Conditioning
is a fact about the arithmetic, not about whether F1 holds.

The same two sensitivities were found on a different chain in
[CONCENTRATOR_MAPPING](CONCENTRATOR_MAPPING.md).

A separate reason not to read the zero-noise column against the other two, and
it is about what that column measures rather than about the spread. At Σγ = 0
the Liouvillian is exactly anti-Hermitian, hence normal (spectral norm of the
commutator ‖[L, L†]‖₂ = 0.0 against 2.52 under the IBM profile; in the
Frobenius norm the contrast is 0.0 against 34.6), so it is not the non-normal
case the unit was chosen for. And the F1 reflection degenerates to λ ↦ −λ,
which holds for **any** Hermitian H: the Liouvillian spectrum is
{−i(E_a − E_b)}, so swapping the two indices negates every eigenvalue. The
column therefore carries no information about dephasing at all. The numerical
check agrees: random real-symmetric H that are not Heisenberg chains land at
the same floor at Σγ = 0. No band is quoted for them, because the theorem above
is the claim and any band would be an artifact of how many H were drawn.

The band is also not an N-independent constant. The same measurement over two
decades of J gives 2.1 to 4.0 at N=2 and 69.0 to 72.4 at N=5, so a larger chain
must be graded against its own N. Comparing at matched J the growth from N=2 to
N=5 is 28.5x, 35.1x and 17.5x at J = 0.1, 1 and 10; the point is the growth, not
any one of those ratios.

This column previously read "100% / 100% / 100%", from a greedy
first-fit inside an absolute tolerance of 1e-6. That score was measuring
its own matcher: tightening the tolerance saturates the printed
percentage long before it fixes the pairing, so a 100% entry was weaker
evidence than a 91% one, not stronger.

The max decay rate for non-oscillating modes equals 2 x Σγ (0.3199).
That one **is** exact, but by the theorem rather than by the printed
digits: the Liouvillian has a stationary mode at λ = 0, so F1 puts a
mode at λ = −2Σγ. And none decays faster, because every Lindbladian
eigenvalue has Re λ ≤ 0, so its F1 partner −2Σγ − λ does too, which is
exactly Re λ ≥ −2Σγ. Read off the run instead, the
same eigensolver floor applies as everywhere else on this page: the
computed maximum is 0.3199000000000025 against 2Σγ = 0.3199000000000000,
a gap of 2.5e-15, about 0.96 ε · spectral radius. These are the pure decay modes at maximum Pauli
weight (XOR drain, the modes where every qubit carries an X or Y operator and dephasing is strongest). The oscillating modes reach at most 0.2971
(sacrifice) and 0.2559 (uniform), always below this ceiling.

---

## What this means

The sacrifice-zone formula is not an engineering hack. It is a precise
intervention in the mode structure of a quantum resonator:

1. The cavity has 43 distinct frequencies (at zero noise)
2. Noise damps these modes without changing their frequencies
3. Concentrating noise on one edge qubit creates a spatial gradient
4. Modes localized away from the sacrifice see less damping
5. The slowest oscillating modes survive 2.80x longer
6. These are the modes that carry quantum information across the chain

The sacrifice zone tunes the resonator. What it changes most is which
notes ring longest; it detunes them a little as well, by less than the
binning used here can see.

---

*See also:*
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md) (the eigenfrequencies),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md) (24,073 records, r* threshold),
[IBM Concentrator](IBM_CONCENTRATOR.md) (hardware test, 2-3x measured),
[Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md) (the resonator paradigm),
[Energy Partition](../hypotheses/ENERGY_PARTITION.md) (2x decay law)
