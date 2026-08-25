# The Mediator Noise Gate at N=11

**Status:** Tier 1 (measured, reproduced from the committed record).
**Date:** 2026-08-23.
**Arc:** `the_gate_that_does_not_gate`, NextStep (1), in
[`OpenArcsRegistry.cs`](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)
(the entry begins at line 7896).

## The system, in full

Everything below runs on one object. Sites 0 to 10 in an open line, one bond
between each neighbouring pair, every bond Heisenberg with the same strength:

    H = Σ_bonds J (X_a X_b + Y_a Y_b + Z_a Z_b),   J = 1 for every bond

    0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 - 10

Site 5 is the one called the meta mediator. Blocks A = {0,1,2,3,4} and
B = {6,7,8,9,10} are the two halves it separates. Every site carries Z-dephasing
at γ = 0.05 except site 5, whose rate γ_M is the swept parameter. The initial
state is |Φ+⟩ on sites (0,1) tensored with |0⟩ on the other nine. The observable
is the peak of I(A:B) over a window t ∈ [0, 20], entropies in log base 2. The
integrator is RK4 at dt = 0.05.

**The word "bridge" describes the intent, not this graph.** The March code calls
this configuration Level 3 of a recursive construction,
`Topology.MediatorBridge(level, jInternal, jBridge, jMeta)`, in which two
five-site sub-bridges are joined through a meta mediator. All three couplings
default to 1.0, and at those defaults the recursion emits a plain uniform path
graph, which is what the March run printed as its bond list
(`simulations/results/mediator_bridge_scale.txt` lines 34 to 43, ten bonds, all
J=1). The same holds one level down: Level 2 is `Chain(5, [1,1,1,1])`,
a uniform five-site line with the mediator at site 2
(`compute/RCPsiSquared.Propagate/Topology.cs` lines 47 to 79). So both objects on
this page are uniform open chains with one distinguished interior site, and the
distinguished site is the exact centre in both. No bridge structure is present in
either. The arc records this independently.

## The question

On 2026-03-21 this configuration was run and γ_M was swept over
{0, 0.01, 0.05, 0.1, 0.2, 0.5}. Peak I(A:B) moved from 0.699531 to 0.646889,
about 7.5 % over the whole sweep. Twelve lines above in the same output file, the
meta coupling J_meta was swept over {0.25 … 3.0} and moved the same observable
from 0.438644 to 0.891614, a doubling. Read literally, the coupling is the strong
knob and the mediator's own noise is the weak one, which inverts the ranking in
[MEDIATOR_AS_QUANTUM_TRANSISTOR](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md),
where γ_M is "Knob 1: Decoherence Rate γ_M (Primary Gate)" (that document's
line 153).

Both sweeps were read as a maximum over a coarse grid of eight measurement TIMES,
{0, 2, 4, 6, 8, 10, 15, 20} (`compute/RCPsiSquared.Propagate/Program.cs` lines 850
and 879). One configuration of this Hamiltonian is committed three times over,
each on a different measurement grid, and the three disagree:

| measurement times | committed value | source |
|---|---|---|
| eight points, spacing 2 to 5 | 0.685503 | `mediator_bridge_scale.txt` line 90 |
| integers 0…20 | 0.733686 | `simulations/results/pull_principle.txt` line 33, published as 0.734 in [SCALING_CURVE](SCALING_CURVE.md) line 107, in a column headed MI(Bridge A:B) |
| 41 points, spacing 0.5 | 0.777324 | `mediator_bridge_scale.txt` line 62 |

(All three are γ_M = 0.05, J_meta = 1.0, i.e. the fully uniform chain. On the
coarse grid that configuration is simultaneously the t = 4.0 curve point, the
γ_M = 0.05 row and the J_meta = 1.00 row.) The spread across the three, 13 % of
the smallest, is larger than the 7.5 % effect being claimed from the first of
them.

So: is the flat γ_M response real insensitivity, or an artefact of how it was
measured?

## What the repo held before this page

**A positive law for exactly this quantity, which the repo files under a name
that does not contain the word mediator.** F64 (`docs/ANALYTICAL_FORMULAS.md`
lines 1849 to 1880) gives a single site's contribution: γ_eff = γ_B·|a_B|²,
scope line "Z-dephasing on any single site B", verified across chain, ring,
star hub and leaf, Y-junction and K₅. It has two readings: the first-order one
above, on the Hamiltonian's eigenvectors, and an exact one,
−Re(λ_k) = 2γ_B·|v_k(B)|², mode by mode on the coherence-sector Liouvillian at
any γ. Its own summary line reads: "γ_B appears as a constant prefactor.
It is not diminished by intervening sites." It names the quantity it says
matters, the site's own amplitude weight, and Reading 3 below tests the
first-order reading of it against this observable and refutes it. F66
(line 1925) records that interior placement of the dephased site is open and
carries one measurement showing the endpoint results do not carry over.
`docs/proofs/PROOF_ABSORPTION_THEOREM.md` Theorem 2 (line 354) is the same
content without F64's sector restriction: Re(λ) = −2 Σ_l γ_l·light_l(v), the
per-site light share, for any Hermitian H and any right eigenvector of the full
Liouvillian.

**One γ_M sweep at N=11 exists and it is the one this page audits.** The March
run itself is a clean one-factor sweep, all eleven rates at 0.05 with `gammas[5]`
overridden alone (`Program.cs` lines 874 to 880). What does not exist is such a
sweep on a measurement grid fine enough to read, which is the whole of Reading 1.
One other committed N=11 run touches the mediator's rate:
[RELAY_PROTOCOL](RELAY_PROTOCOL.md) lines 115 to 141 lowers site 5 from 0.05 to
0.005, but as stage 3 of a six-stage protocol that lowers a different site group
at each stage, and its three headline numbers are not one observable: 0.734 is a
maximum over the integer grid, while 0.759 and 0.723 are single readings at
t = 4.68, and 0.723 is a different Hamiltonian (2:1 couplings, not uniform J).
It is not a γ_M contrast and its numbers are not comparable to this page's.

**The precedent for the artefact is committed and is the same shape.** The
correction note of 2026-04-23 in
[RECEIVER_VS_GAMMA_SACRIFICE](RECEIVER_VS_GAMMA_SACRIFICE.md) line 25 records a
~0.38 spaced grid that missed the real peak, put uniform-J baselines "about
factor 2 too low" and thereby "inflated apparent J-modulation boost ratios"; two
claims were retracted there (lines 97 and 98), both of them N-trends.

**The mediator's noise was already measured at N=5, and it harms.**
[GAMMA_CONTROL](GAMMA_CONTROL.md) rows 125 to 127 hold a one-factor pair with
only γ₂ moving, and the sentence at line 136 states the reading: 0.416 at
γ_M = 0.01, 0.338 at 0.05, 0.192 at 0.20, falling monotonically. Those are
readings at the fixed time t = 5.0, not peaks, which matters here because
Reading 4 below finds the fixed-time and peak functionals differ by a factor of
about eight at N=11 and reverse direction between the two chain lengths. That
document was repaired on the same day as this page, in a parallel session; its
earlier headline read the opposite way because it compared profiles whose total
noise differed, and the correction is the 2026-08-23 entry of
`docs/CAUGHT_ERRORS.md` headed "a profile comparison whose arms carried
different Σγ".
It confirms the N=5 polarity datum and agrees with THE_OTHER_SIDE §22; what was
open is the N=11 behaviour, which is this page.

**Grid conventions exist but drift.** `experiments/RESONANT_RETURN.md` lines 61 to
67 is a labelled metric note naming both conventions in use ("Sum_MI@5", measured
at fixed t = 5.0, against "Peak Sum_MI", the maximum over all t > 0, with "Peak
values are always higher"), and RECEIVER_VS_GAMMA_SACRIFICE line 27 states its
grid in prose. What does not exist is a pinned convention: `docs/GLOSSARY.md`
defines Sum-MI as an observable without a reading time, and the default run mode of
the propagate engine alone uses four different measurement grids
(`Program.cs` lines 702, 762, 850 and 931).
`docs/CAUGHT_ERRORS.md` carries at least four grid-artefact findings, the nearest
in kind at line 467, where TEMPORAL_SACRIFICE claimed an "exact timestep" on a
0.5-spaced grid its own pending list calls too coarse.

**The March session retired the mediator for a different reason, and was right.**
[SCALING_CURVE](SCALING_CURVE.md) lines 131 to 134 lists "Special mediator nodes"
among its falsified items, on the ground that every qubit in a Heisenberg chain
mediates between its neighbours. That is a topological argument, and on a uniform
chain it is simply correct: site 5 is an ordinary interior site.

**Nothing else returned.** `fw.Confirmations` (24 entries) has no mediator, γ_M or
bridge-MI entry. `docs/GLOSSARY.md` has no row for mediator, γ_M or peak MI.

## The instrument

Both the Heisenberg bond term and the Z-dephasing dissipator conserve popcount,
and the initial state has support only on popcount 0 and 2, so the whole
trajectory lives in a space of dimension 1 + C(N,2): 56 at N=11 rather than 2048,
with ρ of size 56×56 rather than 2048×2048.

`simulations/bridge_sector.py` carries the C# engine's conventions verbatim
(qubit 0 as the most significant bit, the factor 2 on XX+YY, the mask
−2 Σ_{k: bit differs} γ_k, the |t − t_meas| < dt/2 measurement rule, Hermiticity
enforced after each full step, log base 2). It reproduces the committed
two committed columns of the March run in full: for I(A:B), the ten printed
Test-1 curve points from t = 2 to t = 20 (the trivial t = 0 row is not checked),
the six γ_M rows and the seven J_meta rows; and for I(A:D), the same ten curve
points. Thirty-three numbers, worst deviation 4.9e-7 against a file printed to
six decimals. What the check does NOT cover is the rest of that run: Test 0 at
N=5, and the sixty correlator values of Test 3. Test 1 took 587 s in March and
takes under a second here, on this machine; the two timings are not from the same
hardware and the ratio is not a benchmark.

## Reading 1: the grid halves the effect

Same Hamiltonian, same γ profiles, same RK4 step. Only the measurement grid
changes.

| measurement grid | span of peak I(A:B) over γ_M = 0 → 0.5 |
|---|---|
| March, 8 points | 7.53 % |
| dt = 0.5 | 15.93 % |
| dt = 0.1 | 15.56 % |
| dt = 0.05 (every RK4 step) | 15.57 % |

The converged answer is 15.6 %, just over twice the published 7.5 %. It is
converged in three separate senses: in the measurement grid (the table), in the
integrator (refining dt from 0.05 to 0.00625 moves it to 15.61 %), and in the
window (t_max = 20, 40 and 80 give the same 15.566 %).

The mechanism is the peak time. The maximum sits at t\* = 3.45 at γ_M = 0 and
drifts to 3.60 by γ_M = 0.5, so the coarse grid's nearest sample, t = 4.0, reads
a flank rather than a peak, and reads it at a different distance from the peak
for each γ_M.

The same correction has to be applied to the arm it is compared against, or a
corrected number would be set against an uncorrected one. On the converged grid
the J_meta sweep runs 0.447193 → 0.988352 over J_meta = 0.25 → 3.0, a rise of
121 % where the coarse grid said 103 %, and it is monotone: the dip at
J_meta = 1.0 in the committed table, which made the published curve
non-monotonic in between (the arc's registry entry recorded that phrasing and,
since 2026-08-23, records the dip as a sampling artefact), is a sampling
artefact of the same kind.

So the arc's qualitative reading survives its own correction, and the correction
runs against it: both arms grow, the mediator's more than the coupling's in
relative terms, and the coupling still moves the observable several times
further. Putting a number on "several" requires a convention, because the two
sweeps have incommensurate ranges: γ_M runs from 0, a multiplicatively infinite
range, while J_meta runs over a factor of 12. Two constructions, both stated
rather than chosen silently:

- raw spans over the ranges as swept: 121 % against 15.6 %, about 8:1
- dimensionless response at the operating point (J_meta = 1, γ_M = 0.05), central difference in log space with h = 0.05: ∂lnMI/∂lnJ_meta = +0.248 against ∂lnMI/∂lnγ_M = −0.026, about 9.5:1

The direction is robust and the order of magnitude is stable; the single digit is
not a property of the system.

## Reading 2: the closure is real, and this observable cannot locate it

| γ_M | 0 | 0.5 | 1 | 2 | 5 | 10 | 50 | 200 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| peak I(A:B), t_max = 20 | 0.801 | 0.677 | 0.613 | 0.536 | 0.455 | 0.426 | 0.194 | 0.039 | 0.003 |

Read on its own this table says the mediator's noise closes the correlation
between the halves, monotonically and eventually completely. The direction is
right and the rest of it is a property of the window, not of the system.

From γ_M around 5 upward the reported maximum sits at the window edge, so it
is not a peak but the last value sampled, and it rises when the window is
lengthened:

| γ_M | 5 | 10 | 20 | 50 | 100 | 200 |
|---|---|---|---|---|---|---|
| t_max = 20 | 0.455 | 0.426 | 0.348 | 0.194 | 0.096 | 0.039 |
| t_max = 80 | 0.461 | 0.461 | 0.460 | 0.425 | 0.336 | 0.212 |

At γ_M = 50 the value more than doubles. The γ_M = 0 reference is
window-stable at 0.801345, so half of it is 0.4007, and on the twenty-unit
window that is crossed between γ_M = 10 and 50 while on the eighty-unit
window it is crossed between 50 and 100. Every peak time in the lower row is
still at the window edge (80.05, 80.02, 80.01), so the eighty-unit answer is not
converged either.

**No closure threshold can be quoted from this observable.** The peak-over-window
functional does not converge at large γ_M: as the mediator is driven towards
the Zeno limit the correlation between the halves builds more slowly, so a longer
window always finds more of it, and "the maximum" tracks the window rather than
the physics. What survives is the ordering, at every fixed window the response is
monotone in γ_M and the mediator eventually dominates. Locating where needs a
window-stable functional, and that is an open item below rather than a number
here.

(The reduced RK4 steps in the table are set by max|mask|, which in the
popcount-{0,2} sector is 2γ_M + 6γ rather than 2γ_M plus twice the
sum of the others: at most four bits differ between bra and ket labels, so at most
three non-mediator sites can contribute at once. At γ_M = 50 that is 100.3
and the step is 0.014955, which is what the run prints.)

## Reading 3: the mediator is not a distinguished site

The obvious mechanism for the flatness is that the mediator's seat is special.
It has an exact structural peculiarity, so this is worth testing rather than
assuming. The chain is reflection-symmetric about its centre, so every
single-excitation eigenvector is symmetric or antisymmetric under that
reflection and every antisymmetric one vanishes at the centre site exactly.
Computed from the single-excitation block of this Heisenberg chain (the ZZ term
shifts the weights but leaves the nodes exact): ⌊N/2⌋ of the N modes have an
exact node at the centre, 5 of 11 here and 2 of 5 at N=5.

**The peculiarity is real and it is measurably inert.** Putting the noise on
each site in turn, everything else held fixed. Every row moves the same rate by
the same amount, so Σγ changes identically down the column and the comparison
across sites is total-noise matched by construction; only the seat differs:

| noisy site | exact nodes there | max \|a\|² | span of peak I(A:B) over γ = 0 → 0.5 |
|---|---|---|---|
| 0 | 0 | 0.1781 | 39.67 % |
| 1 | 0 | 0.1781 | 36.33 % |
| 2 | 0 | 0.1781 | 28.04 % |
| 3 | 0 | 0.1781 | 19.47 % |
| 4 | 0 | 0.1781 | 12.30 % |
| **5 (the mediator)** | **5** | **0.1818** | **15.57 %** |
| 6 | 0 | 0.1781 | 23.96 % |
| 7 | 0 | 0.1781 | 18.04 % |
| 8 | 0 | 0.1781 | 9.97 % |
| 9 | 0 | 0.1781 | 4.59 % |
| 10 | 0 | 0.1781 | 5.02 % |

Two things follow, and the second is the answer to the page's question.

**The static amplitude weight explains none of this.** The largest single-mode
weight is 0.1781 at every site and 0.1818 at the centre, a spread of 2 %, while
the measured response spans a factor of 8.6 across the chain. F64's first-order
reading is built on the H-eigenvector amplitude, and that amplitude is flat here
to within 2 %, so it cannot carry a factor of 8.6. (The table reports the largest
single-mode weight per site; F64 selects a particular mode rather than the
largest, so this rules out the reading as applied here, not every possible
mode-selection rule.) Nor does the nodal structure show up: the one site that has nodes
sits at 15.57 %, between its two neighbours and *above* site 4's 12.30 %. So the
nodal mechanism this reading was written to test is refuted by its own table.
What the response does track is where the excitation actually is: it falls away
from the prepared pair on sites 0 and 1, which is the light share of Theorem 2
read on the state rather than on the Hamiltonian's eigenvectors.

**The mediator's 15.6 % is unremarkable.** Among the nine interior sites the
response runs from 4.59 % to 36.33 %, and site 5 ranks fourth of the nine, just
below their median of 18.04 %.
Nothing about the flat γ_M response is a property of being the mediator. It is
what an interior site at that distance from the prepared pair gives, and site 5
is an interior site at that distance.

This is where the March session's own conclusion returns, reached from the
direction their data did not carry. [SCALING_CURVE](SCALING_CURVE.md) lines 131
to 134 retired "Special mediator nodes" on the topological ground that every
qubit in a Heisenberg chain mediates between its neighbours. That argument is
correct and this table is its measurement: on a uniform chain the mediator is an
ordinary site, and the transistor reading needs the site to be distinguished
before either of its knobs can mean what it says.

A caveat this table does not remove. Moving the noisy site changes two things at
once, its distance to the prepared pair and its role relative to the observable
(sites 0 to 4 lie inside block A, site 5 is traced out, sites 6 to 10 lie inside
block B). The table is therefore decisive against the nodal reading, which
predicts an anomaly at exactly one site and shows none, and it is not a clean
measurement of distance alone.

## Reading 4: the two chains are not measured the same way

THE_OTHER_SIDE §22, Test 3 reports the N=5 mediator noise sweep as 0.44 → 0.12
at γ_M = 0.5, and §23 extends the claim to the meta mediator by assertion. Set
against the N=11 numbers this looks like a factor of ten in sensitivity.

They are not the same measurement. Three things differ besides the chain length.
The N=5 number is I at the FIXED time t = 5.0 (`simulations/mediator_bridge.py`
lines 571 and 621) while N=11 takes a MAXIMUM over a grid. The N=5 partition is
the 2-qubit end pairs {0,1}:{3,4} while N=11 uses the 5-qubit halves. And the N=5
initial state is Bell(0,1) ⊗ |0⟩ ⊗ |+⟩|+⟩ while N=11 is Bell on vacuum. A fourth
difference is physical rather than observational and is not controlled anywhere
below: the source pair sits one bond from the mediator at N=5 and four bonds from
it at N=11.

The N=5 arm runs dense through the committed primitives of
`simulations/mediator_bridge.py`, so its conventions are the repo's own. §22
quotes two significant figures, but the March output behind it is committed and
carries four: `simulations/results/mediator_bridge.txt` lines 217 and 231 give
MI = 0.4405 at γ_M = 0 and 0.1155 at γ_M = 0.5. The re-derivation gives
0.440471 and 0.115494, which is agreement to every digit that file prints.

(A trap found while checking this, and it is in the reproduction recipes on this
page. Until 2026-08-23 `simulations/mediator_bridge.py` opened that output file
for writing at module level, so merely importing the module, which is what the
N=5 arm and every reader following these instructions does, truncated a tracked
result file to zero bytes. The open is now lazy. Anyone who ran an earlier
version of this page's recipe should check `git status` before committing.)

| configuration | drop over γ_M = 0 → 0.5 |
|---|---|
| N=5 as published (fixed t = 5, end pairs, \|++⟩ on B) | 73.8 % |
| N=5, fine-grid peak instead of the fixed t | 34.9 % |
| N=5, vacuum instead of \|++⟩ on B, fixed t | 49.2 % |
| **N=5, both changes** | **40.8 %** |
| **N=11, matched functional and state, 5-qubit halves** | **15.6 %** |
| **N=11, matched functional and state, 2-qubit end pairs** | **10.6 %** |
| N=11, fixed t = 5, 5-qubit halves | 1.7 % |
| N=11, fixed t = 5, 2-qubit end pairs | 13.0 % |
| N=11 as published (coarse max, halves) | 7.5 % |

**The partition is a choice and the page makes it visible rather than silently.**
At N=5 the end pairs ARE the halves flanking the mediator, so "halves" is the
structurally matched partition and "end pairs" is the size-matched one. Both are
defensible and they give different answers: the ratio is 40.8 / 15.6 = 2.6× on
the structural match and 40.8 / 10.6 = 3.9× on the size match. The axis is worth
a factor of 7.9 in its own right (1.66 % against 13.02 % at fixed t), so it cannot
be left implicit.

**No factor can be called the largest contributor.** The three changes do not
decompose: from 73.8 %, the functional alone gives −38.9 points and the state
alone −24.6 points, but both together give −33.0, not −63.5. Worse, the state
factor changes sign with the order in which the changes are applied: at fixed t
the vacuum swap lowers the drop (73.8 → 49.2), at the fine-grid peak it raises it
(34.9 → 40.8). Only the joint change is meaningful. The direction of the
functional's effect also reverses between the chains: at N=5 the fixed-t reading
is far more γ_M-sensitive than the peak (73.8 against 34.9), at N=11 it is far
less (1.7 against 15.6).

What survives is the comparison itself: matched, the two chains differ by a
factor between 2.6 and 3.9, not by the factor of ten the published pair suggests.
Roughly half to two thirds of that published gap is the observable rather than
the system, and Reading 3 offers a mechanism for the part that remains.

## What this leaves standing

On a uniform eleven-site Heisenberg chain with one noisy centre site, the centre
site's own dephasing is a weak knob within the operating range, weaker than the
coupling by roughly an order of magnitude. It is weak for no reason peculiar to
the mediator: the response to γ on any site falls away from the prepared pair,
and site 5 ranks fourth of the nine interior sites. The one structural
peculiarity it does have, being an exact node of ⌊N/2⌋ single-excitation modes,
is measurably inert. It is also not insensitive: at every fixed window the
response is monotone and the mediator eventually dominates, though this
observable cannot say where.

Two published numbers are corrected. The 7.5 % understates the γ_M response by a
factor of two, and that error is a coarse measurement grid, the failure mode the
repo's own 2026-04-23 correction note had already recorded in another experiment.
The comparison between the two chain lengths overstates the difference by a
factor of 2.5 to 3.7, and that error is NOT a grid: the N=5 number is a
fixed-time reading with no grid to be coarse, and Reading 4 shows the three
differences do not decompose, so no single one of them can be named the cause.

**This contradicts a live committed claim.**
MEDIATOR_AS_QUANTUM_TRANSISTOR's falsification of its hierarchy explicitly
exempts the gate: "The transistor properties (threshold, gate control,
directional bias) are real; the hierarchical scaling advantage is not"
(line 398). The measurement here says gate control by the mediator's own
dephasing is the weaker of the two knobs at N=11.

The scoping is narrower than it first looks, and it does not rescue the ranking.
That document's declared object is the same five-site chain measured here: its
section 1.2 sets the topology "A = {0, 1} - M(2) - B = {3, 4}", the N=5 mediator
bridge of THE_OTHER_SIDE §22. Its Appendix-A runs are a Heisenberg ring at N=3
to 5 (lines 512 to 537), and at N=3 a ring is a triangle where no site is a
mediator at all. Its Knob-1 threshold formula is imported from
STAR_TOPOLOGY_OBSERVERS, where the swept quantity is the UNIFORM rate on all
three sites and the measured quantity is a J threshold. So the mediator's own
rate was never isolated anywhere that ranking was established, and where it has
been isolated, at N=5 in GAMMA_CONTROL and at N=11 here, it is the weaker knob
both times.

## What stays open

**Find the law the site profile of Reading 3 obeys.** The response to γ on site
j runs 39.7, 36.3, 28.0, 19.5, 12.3, 15.6, 24.0, 18.0, 10.0, 4.6, 5.0 across the
chain, and this page offers no closed form for it. F64's static amplitude
reading is ruled out here: it is flat to 2 % across the chain while the response
spans a factor of 8.6. The candidate that fits the shape is the state's own
light share, Theorem 2's light_l(v) evaluated on the propagated state rather
than on an H-eigenvector, and the test is cheap: the excitation occupancy
⟨n_j(t)⟩ at t\* is a by-product of runs this page already makes. If that
reproduces the profile, the page's finding sharpens from "the mediator is
ordinary" to a statement about which sites can gate at all.

**The exact readings were not used and the general one already spans this
sector.** Reading 3 tests only F64's first-order H-eigenvector reading. F64 also
carries an exact form, −Re(λ_k) = 2γ_B·|v_k(B)|² mode by mode on the
coherence-sector Liouvillian at any γ, and Theorem 2 of the Absorption Theorem is
more general still: any Hermitian H, the full Liouville space, no sector
restriction. So the gap this page leaves is not a missing law but an unperformed
calculation, and it is a 56×56 one needing no new engine.

**Decide the partition convention, or record that there is none.** Reading 4
leaves 2.6× and 3.9× both standing. The repo has no convention for how to match
subsystem partitions across chain lengths, and this is the second document to
need one. Either a convention lands in the glossary or every cross-length MI
comparison carries both numbers.

**Locate the closure at all.** Reading 2 shows the peak-over-window functional
does not converge at large γ_M, so no threshold can be quoted: the twenty-unit
window brackets the half-value between γ_M = 10 and 50, the eighty-unit window
between 50 and 100, and neither is converged. This needs a window-stable
functional (steady-state MI, or MI integrated over time) before the question is
even well posed. The same functional is what the background item below waits
on.

**Both candidates have since been tried, and both fail, 2026-08-23/24.**
[The Blind Site](THE_BLIND_SITE.md) §7 rules out MI integrated over time: past
the transient every γ profile decays toward the same limit, so a longer window
dilutes every span alike and only the ratio between seats survives.
[The Seat That Cuts](THE_SEAT_THAT_CUTS.md) §1 rules out steady-state MI more
sharply. Wherever dephasing at a SINGLE seat leaves no blind subspace the
stationary state is the maximally mixed state of the excitation sector, at every
seat and every rate; for a support of several seats that needs connectivity as
well, and that page's opening carries the counterexample, so what licenses that
page's §1 kernels is the exact rank taken there profile by profile. The
steady-state mutual information is then a function of N and of the sector alone,
with no γ in it, and the span it reports is zero at every seat. Which
number it is depends on the sector, and for this page's own state that is not
the number quoted below. It is window-stable because it has stopped
looking. The premise was already proven next door:
`PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md` Consequence 2 makes the sector
populations a complete invariant for ρ(∞), and no γ enters it.

**What carries over to this page is the PREMISE, and neither the number nor the
blindness.** That page's closed form, log₂N − ((N+1)/N)·log₂((N+1)/2) bits, is
read in the **single-excitation** sector on halves of size (N−1)/2 with the
centre in neither. The state this page propagates has support on popcount 0 and
2, so it lands on a different limit and the closed form does not apply to it. The
blindness does not reach this page either, for the reason spelled out two
paragraphs down: every profile swept here puts γ > 0 on at least ten of the
eleven seats, so nothing is blind at any of them. What does apply, because it
belongs to the proof rather than to either, is that the limit here is fixed by
the sector populations alone and no γ profile can move it.

**And the obvious replacement does not work either, on this comparison.** The
**dimension** of the stationary manifold is window-free and rate-free and carries
the blindness count exactly, dim ker L_SE = 1 + deg gcd(χ(H_left), χ(H_right))
for dephasing at one seat of an open chain with NO ZERO BOND, read in the
single-excitation sector, which on the uniform HEISENBERG chain is
1 + (gcd(2j+1, N) − 1)/2 and on the uniform XY chain is gcd(j+1, N+1)
([The Seat That Cuts](THE_SEAT_THAT_CUTS.md) §§2-3; the count and the criterion
are **F157** since 2026-08-24, and the kernel identity displayed here is proved
since 2026-08-25 in
[the span and node-lemma proof](../docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md),
as the simple-spectrum case of dim ker L_SE = 1 + dim commutant(H restricted to
the seat's Krylov complement); the XY clause is written here as the bare gcd because
1 + (gcd(j+1, N+1) − 1) collapses to it, which is arithmetic and not a second
shape of law). The fence is not a
well-formedness condition: at a zero bond the two-halves form goes WRONG
rather than undefined, totally so on the Heisenberg book and not totally on
XY. Nothing on this page sits near it, the chain here being uniform and
zero-free throughout. But every profile in the
γ-sweep above puts γ > 0 on at least ten of the eleven seats (the γ = 0
configuration discussed later on this page is not one of the swept arms). That
by itself does not settle it, since a one-seat law says nothing about a ten-seat
support; what settles it is connectivity, this page's bridge being connected and
every arm's support reaching every seat or all but one. On a connected graph such
a support forces a one-dimensional kernel outright: a surviving off-diagonal
entry would need both its sites outside the support, of which there is at most
one, so the operator is site-diagonal, and commuting with a connected H then
leaves only the identity. The Blind Site §7's 23 uniform-chain profiles carry
the same conclusion a second way, their one-dimensional kernels certified by
the exact rank profile by profile. So the kernel is one-dimensional throughout,
and the dimension reports a span of zero here exactly as the mutual information
does. It discriminates a different
comparison, single-seat support, which this page never posed. The closure
question therefore stays open with two ruled-out candidates and no instrument
yet.

**A quiet-background comparison needs a different observable, and this one
cannot give it.** The natural question, how much of the weak response is the
mediator competing against the other ten sites' noise, cannot be asked with a
peak-over-window functional: with every rate at zero the dynamics is unitary,
I(A:B) never settles, and the "peak" is only the largest sample so far, climbing
to t\* = 79.3 in a window of 80 and still rising. The check is in
`bridge_node_weights.py`. Any answer here has to come from a window-stable
functional first.

**Distance and chain length are not separated.** The source pair sits one bond
from the mediator at N=5 and four bonds from it at N=11, and Reading 3 shows
distance from the prepared pair is the strongest thing varying along the chain.
Reading 4's residual factor is therefore not attributable to chain length. A
five-site chain with the source moved, or an eleven-site chain with the Bell pair
prepared adjacent to the centre, separates them at no cost.

## Reproduction

Run from the repository root, no arguments. Every table above appears in this
output except the three committed values near the top, which are read out of
`mediator_bridge_scale.txt` and `pull_principle.txt` and are not produced here.

```bash
python simulations/bridge_gate_calibrate.py   # the calibration check: 33 committed March numbers
python simulations/bridge_flat_gate.py        # Readings 1 and 2, and the J_meta correction
python simulations/bridge_level_decomp.py     # Reading 4, N=5 dense + N=11 sector
python simulations/bridge_node_weights.py     # Reading 3, Reading 2's window table,
#                                             and the convergence checks
#                                             all four import simulations/bridge_sector.py,
#                                             which is the propagator itself
```

Defaults in all of them: J = 1 on every bond, γ = 0.05 on every site but the
mediator, t_max = 20, RK4 dt = 0.05 reduced only where stability requires it (the
reduced step is printed per row). The N=5 arm of the decomposition propagates by
matrix exponential rather than RK4 and is not covered by the calibration check; against the
committed C# Test-0 value at t = 5 it agrees to 6e-5, which is the two
integrators disagreeing, not either of them being wrong.

The calibration check is the thing that would fail: any change to the propagator has to keep
reproducing those thirty-three numbers to the six decimals the March file
prints. It currently passes at a worst deviation of 4.9e-7.

## See also

- [THE_OTHER_SIDE](../hypotheses/THE_OTHER_SIDE.md), §22 and §23, the source of the N=5 pair and of the meta-mediator assertion
- [MEDIATOR_AS_QUANTUM_TRANSISTOR](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md), the knob ranking this page contradicts at N=11
- [GAMMA_CONTROL](GAMMA_CONTROL.md), the N=5 one-factor sweep, repaired the same day
- [RECEIVER_VS_GAMMA_SACRIFICE](RECEIVER_VS_GAMMA_SACRIFICE.md), the 2026-04-23 correction note, the repo's own precedent for a coarse grid manufacturing a trend
- [SCALING_CURVE](SCALING_CURVE.md), the third committed reading of this Hamiltonian and the March topological retirement of the mediator
- [RELAY_PROTOCOL](RELAY_PROTOCOL.md), the one committed N=11 run that varies the mediator's rate, in company
- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), Theorem 2, and F64 in [ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md), the per-site weighting Reading 3 applies
- [ON_FIVE_PAGES_THAT_NEVER_MET](../reflections/ON_FIVE_PAGES_THAT_NEVER_MET.md), the outward reading of the archaeology this page belongs to
