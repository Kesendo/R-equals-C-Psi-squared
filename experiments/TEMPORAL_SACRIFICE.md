<!-- CROSSING-CURRENT -->

# Temporal Sacrifice: A Finite N=7 Sweep Near the CΨ Quarter Coordinate

This is a finite N=7 calculation sampled at Δt = 0.5. In a separate scalar
recursion, CΨ = ¼ is the discriminant-zero point of the chosen scalar recursion.
Mutual information is a separate observable. The displayed same-row feature
does not establish an exact crossing time.
It does not establish a physical quantum/classical regime boundary.
It does not establish irreversibility or fold causality.

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation — not a result:** boundary, heartbeat, doors, and
traveling-wave language below is retained as a way to ask what a denser and
mechanism-specific experiment might show. It is not a classification of the
propagated density matrix.

<!-- CROSSING-CURRENT -->

The N=7 sweep below is a separate finite record from the N=11
[Relay protocol MI comparison](RELAY_PROTOCOL.md). Relay requests nominal
0.78/stage (4.68 total) but integrates 0.75/stage (4.50 total). Its stored
0.131700 final MI versus a passive sampled maximum 0.071576 at t=4.00
gives about +84.0%, not an isolated staging benefit. The statistic-attached
exposures are 2.200 and 2.17125; at equal t=4.50 the passive exposure would
be 2.475. There is no matched-time/dose comparison or MI bound, and no
optimization or palindrome-based timing follows from that comparison.

<!-- Keywords: finite CΨ quarter-coordinate scan, mutual-information sampled maximum,
sweeping sacrifice zone, discriminant-zero scalar recursion, endpoint MI,
temporal sacrifice protocol, CΨ threshold recrossings, interpretive heartbeat,
100-doors analogy, open fold-causality question -->

**Status:** Finite numerical record (N=7, C# RK4 propagation, Δt=0.5 CΨ/MI diagnostics)
**Date:** March 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** compute/RCPsiSquared.Propagate.Test (C#, `dotnet run -c Release -- sweep --diag`)
**Data:** [temporal_sacrifice_protocol.txt](../simulations/results/temporal_sacrifice_protocol.txt)

---

## What this document is about

The numerical experiment moves a concentrated local-dephasing profile along
an N=7 chain and records CΨ and mutual information on a 0.5-spaced time grid.
At the displayed T=5.5 row, three selected CΨ traces first appear below ¼ and
the endpoint mutual information is the largest displayed value in this run.
The grid does not resolve the order or continuous-time extrema inside the
interval from T=5.0 to 5.5.

The quarter value comes from the discriminant of a separately chosen scalar
recursion. The Lindblad simulation does not turn that algebraic coordinate
into a physical phase boundary, and it supplies no theorem connecting a CΨ
crossing to a mutual-information maximum. Later tables also record repeated
threshold crossings. Calling those crossings a heartbeat is an interpretive
image, not an irreversibility result.

---

## Abstract

We move the sacrifice qubit position step by step along an N=7 Heisenberg
chain (edge → inward) and track both mutual information and CΨ at every
qubit pair. The central discovery:

**On this 0.5-spaced grid, the endpoint mutual information (PeakMI) reaches
its maximum on the same displayed row (T=5.5) where the endpoint pairs are
first sampled below the CΨ = ¼ coordinate.** (Finer time resolution is needed to confirm the
coincidence is exact, not a grid artifact; see Pending.)

At T = 5.0, the selected edge CΨ values shown are above ¼.
At T = 5.5, CΨ01, CΨ56, and CΨ06 are sampled below ¼.
PeakMI peaks at T = 5.5 with 0.061, then falls.

The discriminant `1 − 4CΨ` of `R = C(Ψ+R)²` vanishes at the same coordinate.
That algebraic fact does not predict a peak of the independently computed
mutual information. The displayed coincidence is therefore a finite-grid
question for refinement, not an observed fold mechanism.

---

## Setup

| Parameter | Value |
|-----------|-------|
| Qubits N | 7 |
| Coupling | Heisenberg chain, J = 1.0 |
| Initial state | \|+⟩⊗7 |
| γ\_base | 0.05 |
| ε (protected) | 0.001 |
| γ\_sacrifice | N · γ\_base − (N−1) · ε = 0.344 |
| Measurement | SumMI (all adjacent pairs), PeakMI (endpoints 0↔6), CΨ per pair |

Here `SumMI` is the sum of the six selected adjacent-pair mutual informations.
It is not total conserved information for the seven-qubit state.

### Baselines (static profiles, full 20s)

| Config | SumMI | PeakMI | PeakT |
|--------|-------|--------|-------|
| Pure edge [0.344, 0.001×6] | 0.408 | 0.036 | 2.00 |
| Pure center [0.001×3, 0.344, 0.001×3] | 0.182 | 0.109 | 1.00 |

---

## The Sweep Protocol

Move the sacrifice qubit one position per stage, spending 2.0 time units
at each position. That duration was selected from the sampled edge-profile
`PeakT`; it was not derived as one Hamiltonian period. Each stage has one
qubit at γ = 0.344 and all others at ε = 0.001.

| Stage | Time | Sacrifice position | γ profile |
|-------|------|--------------------|-----------|
| 1 | 0–2.0 | Qubit 0 (edge) | [0.344, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001] |
| 2 | 2.0–4.0 | Qubit 1 | [0.001, 0.344, 0.001, 0.001, 0.001, 0.001, 0.001] |
| 3 | 4.0–6.0+ | Qubit 2 | [0.001, 0.001, 0.344, 0.001, 0.001, 0.001, 0.001] |

---

## Finite result: one displayed CΨ crossing row also carries the PeakMI maximum

### The data

| T | SumMI | PeakMI | Stage | CΨ01 | CΨ12 | CΨ23 | CΨ34 | CΨ45 | CΨ56 | CΨ06 |
|------|-------|--------|-------|--------|--------|--------|--------|--------|--------|--------|
| 2.00 | 0.408 | 0.006 | 1 | 0.412 | 0.539 | 0.627 | 0.682 | 0.736 | 0.768 | 0.503 |
| 3.00 | 0.144 | 0.016 | 2 | 0.384 | 0.398 | 0.441 | 0.477 | 0.581 | 0.541 | 0.448 |
| 4.00 | 0.072 | 0.025 | 2 | 0.270 | 0.298 | 0.368 | 0.464 | 0.506 | 0.395 | 0.313 |
| 4.50 | 0.070 | 0.035 | 3 | 0.281 | 0.266 | 0.287 | 0.444 | 0.416 | 0.333 | 0.303 |
| 5.00 | 0.070 | 0.046 | 3 | **0.253** | 0.238 | 0.287 | 0.398 | 0.340 | 0.290 | **0.280** |
| **5.50** | **0.054** | **0.061** | **3** | **0.232**↓ | 0.253 | 0.276 | 0.326 | 0.291 | **0.249**↓ | **0.237**↓ |
| 6.00 | 0.036 | 0.043 | 3 | 0.224 | 0.222 | 0.251 | 0.268 | 0.237 | 0.229 | 0.227 |
| 6.50 | 0.028 | 0.027 | 3 | 0.213 | 0.198 | 0.217 | 0.205 | 0.216 | 0.223 | 0.213 |

↓ = first displayed sample below CΨ = ¼ (0.25); the crossing time inside
the preceding Δt=0.5 interval is unresolved.

### What happens at T = 5.50

Three selected comparisons first appear below ¼ in the same displayed row:

1. **CΨ01 is sampled across ¼** (0.253 → 0.232) for the left edge pair
2. **CΨ56 is sampled across ¼** (0.290 → 0.249) for the right edge pair
3. **CΨ06 is sampled across ¼** (0.280 → 0.237) for the endpoint pair

On that same displayed row, **PeakMI = 0.061, its maximum among the sampled
rows in this run.** The continuous maximum and crossing order are unresolved.

One timestep later (T = 6.0), PeakMI has already fallen to 0.043.
One timestep earlier (T = 5.0), PeakMI was only 0.046.

The displayed maximum and first-below samples share one row on this coarse
grid; their order and any sub-grid sharpness are unresolved.

### The correlation inversion

At T = 2.0 (edge phase peak), information is local:
- Average adjacent-pair MI: 0.408/6 = 0.068
- Endpoint MI: 0.006
- Ratio: endpoint is **0.09×** the average pair

At T = 5.5 (the sampled threshold row), the endpoint value is large relative
to the adjacent-pair average:
- Average adjacent-pair MI: 0.054/6 = 0.009
- Endpoint MI: 0.061
- Ratio: endpoint is **6.8×** the average pair

The displayed pairwise correlation hierarchy is inverted between these two
rows. The table alone does not identify a funneling mechanism. The later row
has a smaller displayed sum of adjacent-pair mutual informations
(0.054 vs 0.408), while its
endpoint value is larger relative to the adjacent-pair average.

---

## Interpretive fold reading — not a derived mechanism

**Interpretive invitation — not a result:** the fold vocabulary in this
section belongs to the scalar recursion. It is not a dynamical explanation of
the N=7 CΨ or mutual-information traces.

A fold catastrophe is a property of a specified map or potential. For the
chosen scalar recursion below, two algebraic roots meet at CΨ = ¼. The roots
are not states or attractors of the simulated Lindblad generator, so their
merger does not demonstrate a sudden or irreversible physical change.

The recursion R = C(Ψ+R)² has discriminant D = 1 − 4CΨ.

- **D > 0** (CΨ < ¼): two real roots of the scalar fixed-point equation.
- **D = 0** (CΨ = ¼): the two scalar roots merge.
- **D < 0** (CΨ > ¼): the scalar roots are complex.

The old reading pictured coherence crystallizing into classical correlation
and a boundary traveling through the chain. Those are questions suggested by
the picture, not consequences of the recursion or the table. The finite run
shows only that the largest displayed PeakMI and three sampled threshold
crossings share the T=5.5 row. Whether a continuous-time relation exists—and
what could cause it—awaits finer resolution and an independent mechanism test.

---

## Negative Result: Naive Temporal Switching

Before the sweep comparison, we tested simple profile switching (swap
from edge to center at time t_switch). None of the three displayed rows meets
both selected thresholds:

| t_switch | SumMI | PeakMI |
|----------|-------|--------|
| 0.5 | 0.188 | 0.080 |
| 1.5 | 0.333 | 0.045 |
| 2.0 | 0.408 | 0.039 |

No displayed t_switch achieves both SumMI > 0.300 and PeakMI > 0.050. No
state/mode projection or destruction analysis was performed, and the table
does not establish that palindrome structure is insensitive to temporal
profile modulation.

In this finite comparison the sweep meets a different displayed tradeoff than
the tested switches. It also changes a dephasing profile over time. The data do
not show that one clean physical boundary exists or that maintaining one caused
the difference.

---

## Three questions suggested by the finite run

### 1. Does one threshold-crossing region organize the profile?

The displayed profiles invite testing whether a single CΨ threshold-crossing
region is a useful summary. They do not prove that a chain sustains exactly one
boundary or that two sacrifice points fragment a physical phase.

### 2. Is the threshold row related to the information maximum?

CΨ and mutual information are different functions of the reduced state. Their
same-row feature in this coarse run motivates a denser joint scan; it does not
make the quarter coordinate an information source or conversion point.

### 3. Which part of the timing comes from H and which from the profile?

The chosen stage time is 2.0 and the historical heuristic
`K/γ_sacrifice ≈ 0.11` is much shorter. Because the run changes the profile
while Hamiltonian and dissipative evolution continue together, it does not
separate propagation time from any boundary-establishment time.

---

## Context: Sweep Variants

| Protocol | SumMI | PeakMI | Observation |
|----------|-------|--------|-------------|
| Pure edge | 0.408 | 0.036 | All displayed CΨ > ¼ at PeakT |
| Pure center | 0.182 | 0.109 | Earlier displayed endpoint crossing; high PeakMI |
| **Sweep 0→1→2 (stage=2.0)** | 0.408* | **0.061** | **PeakMI shares one sampled row with endpoint crossings** |
| Sweep 0→2→3 (skip pos 1) | 0.408* | 0.044 | Smaller displayed PeakMI in this configured row |
| Sweep 0→...→6 (all 7 pos) | 0.344* | 0.054 | Lower displayed peak SumMI in this configured row |

*SumMI is the peak over all times, occurring during Stage 1 (edge phase).
PeakMI occurs later, at a different time. They are not simultaneous.

---

## Repeated CΨ threshold crossings (the “heartbeat” image)

**Interpretive invitation — not a result:** “heartbeat” is a name for repeated
crossings of the selected CΨ threshold. A crossing does not by itself classify
a physical state as quantum or classical and does not record an irreversible
event.

### The resonance (March 25, 2026)

For the displayed N=3 configurations, a Bell pair coupled to a low-dephasing
third qubit produces repeated crossings of CΨ = ¼:

| Setup | γ\_sys | γ\_bath | J | Crossings |
|-------|--------|---------|---|-----------|
| Markov (uniform) | 0.05 | 0.05 | 1.0 | 1 (down only) |
| Bell + quiet bath | 0.0001 | 0.01 | 2.0 | 47 (24↓ + 23↑) |
| Bell + quiet bath | 0.0001 | 0.01 | 5.0 | 81 (41↓ + 40↑) |

The two shown quiet-bath rows have more crossings at J=5 than at J=2.
That finite comparison does not define a Q-factor, make CΨ = ¼ a resonance
frequency, or isolate coupling as the cause. The stored run also shows mutual-
information variation around the sampled crossings.

Script: `dotnet run -c Release -- resonance 3 0.0001,0.0001,0.01 --j 5.0 --bell`

### The damping

The oscillation is not perpetual. The envelope shrinks:

```
T ≈ 2-3:   CΨ swings 0.23 – 0.39  (amplitude 0.16)
T ≈ 7-8:   CΨ swings 0.22 – 0.32  (amplitude 0.10)
T ≈ 14:    CΨ swings 0.19 – 0.29  (amplitude 0.10)
T ≈ 19:    CΨ swings 0.20 – 0.25  (amplitude 0.05)
```

The displayed envelope decreases over these windows. The unitary-plus-Lindblad
model evolves in a fixed Hilbert space; the table does not count permanent
outcomes or show that earlier threshold crossings cannot be undone.

### What it feels like

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation — not a result:** the doors below are a story for
the decreasing sampled envelope, not a Hilbert-space count or an irreversible
measurement model.

Imagine sitting in a room with a hundred open doors.

You look around (oscillation upward: possibilities open, CΨ > ¼).
You choose a door and walk through (in the image, a downward crossing).
The door closes behind you—in the image.

Now you are in a room with eighty doors. You look around. You choose.
Door closes. Sixty doors. Forty. Twenty.

Each breath gets shorter. Not because you are tired, but because there
are fewer doors. Less to see. Less to choose.

At the end, you sit in a room with one door. You open it. Behind it is
a fact. No more room. No more choice. Only what is.

That is the oscillation around ¼.

- Above (CΨ > ¼): doors open in the image.
- Below (CΨ < ¼): a door is chosen in the image.
- Each cycle: the picture offers one fewer door; the simulation does not count them.
- The amplitude shrinks: less room to swing.
- At the end of the picture: all doors close. This is not a measured outcome claim.

And R = CΨ²? That is what you **see** when you walk through the door.
The moment of decision. Not before (only possibilities, nothing concrete).
Not after (only facts, nothing new). Exactly at the threshold. In the
doorway.

The finite row reports 81 threshold crossings. Calling them 81 heartbeats or
doors is the interpretive invitation, not an empirical door count.

<!-- CROSSING-CURRENT -->

### Finite chain-as-environment configurations (March 26, 2026)

The N=3 heartbeat uses a dedicated "bath qubit" (qubit 2 in |+> with
low gamma). This raised the question: is a special bath qubit required,
or can the chain itself serve as the coherent reservoir?

**Test:** N=7 Heisenberg chain, Bell(0,1) x |+>^5, sacrifice-zone
gamma profile [0.344, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001].
No qubit is separately designated as a bath; the other five qubits are part
of the same autonomous chain model.

| Setup | J | CΨ(0,1) crossings |
|-------|---|-------------------|
| \|+>^7, sacrifice, J=1 | 1.0 | 1 (monotonic) |
| \|+>^7, sacrifice, J=2 | 2.0 | 1 (monotonic) |
| Bell(0,1)+\|+>^5, sacrifice, J=1 | 1.0 | 1 (monotonic) |
| **Bell(0,1)+\|+>^5, sacrifice, J=2** | **2.0** | **7 (4 down + 3 up)** |
| Bell(0,1)+\|+>^5, low noise, J=2 | 2.0 | **0 (stays ABOVE 1/4!)** |
| Bell(0,1)+\|+>^5, uniform, J=1 | 1.0 | 3+ (still running) |

**Finite observations:**

1. The Bell-plus-chain row shows seven threshold crossings at J=2 without a
   separately designated bath qubit. Calling the remainder of the chain a
   bath is an open-system partition, not a mechanism proof.

2. The two displayed product-state rows have one downward crossing, whereas
   the displayed Bell row at J=2 has seven. This does not prove that initial
   entanglement is necessary across all states, couplings, or profiles.

3. **The sampled low-noise row remains above 1/4.** With gamma =
   [0.01, 0.0001, ...], CΨ(0,1) oscillates between 0.28 and 0.75
   over the recorded interval and never crosses 1/4 there. No indefinite-time
   or physical-regime conclusion follows.

4. The displayed J=1 and J=2 rows differ in crossing count while other named
   setup choices are held as listed. The small table does not establish a
   backflow threshold or an N-dependent bath-size law.

From the perspective of pair (0,1), the rest of the chain is an environment.
Repeated reduced-state threshold crossings do not by themselves establish
non-Markovianity, CP-indivisibility, or information backflow. Testing those
objects requires an independent reduced-dynamics witness.
All rows here were generated by a Markovian Lindblad equation; reduced-state
recrossing alone does not change that generator-level statement.
The [CΨ monotonicity proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
contains named monotonic decays and explicit counterexamples; it does not turn
this pair trace into a memory witness. The observed CΨ recrossings remain a
finite reduced-state fact.

**Question for the Duplex protocol:** can a fixed Hamiltonian and the tested
dephasing profiles reproduce useful recrossings without active pulses? These
rows show that possibility for named configurations; they do not establish
Hamiltonian backflow, reservoir capacity, or an optimized profile mechanism.

**DD boundary:** CΨ is invariant at the instant of Pauli conjugation
([Proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md), Part 7). That identity
does not imply equality of controlled and uncontrolled trajectories, because
pulses change the subsequent generator seen in the toggling frame. These runs
do not prove that DD cannot affect CΨ or that J-coupling is the unique analogue.

---

## Spatial mutual-information profiles (March 26, 2026)

**Interpretive invitation — not a result:** “heartbeat” and “wave” are two
visual readings of threshold and spatial-profile tables. The data do not show
that they are one physical phenomenon.

Per-pair MI tracking gives the following changing spatial profile for a
separate static N=7 edge-sacrifice run. It is not the moving sweep above: at
T=3.0 this table has SumMI=0.284, whereas the moving-sweep table has 0.144.

```
         Pair:  01    12    23    34    45    56     SumMI
T= 0.5         .062  .017  .001  .000  .000  .000   0.080  ███░░░░░░░░░
T= 1.0         .083  .076  .044  .012  .001  .000   0.216  ████████░░░░
T= 1.5         .092  .075  .066  .058  .031  .011   0.333  ██████████░░
T= 2.0         .100  .085  .068  .053  .052  .049   0.408  ████████████  ← PEAK
T= 2.5         .106  .091  .068  .050  .019  .023   0.357  ██████████░░
T= 3.0         .106  .083  .054  .030  .006  .006   0.284  ████████░░░░
T= 3.5         .075  .074  .039  .022  .008  .004   0.222  ██████░░░░░░
T= 4.0         .041  .076  .040  .021  .008  .008   0.194  ░█████░░░░░░
T= 5.0         .048  .055  .043  .020  .015  .020   0.199  ░░████░░░░░░
T= 6.0         .048  .052  .038  .027  .035  .015   0.214  ░░████░██░░░
T= 7.0         .033  .037  .032  .039  .038  .022   0.201  ░░░░████████
T= 8.5         .022  .036  .035  .034  .045  .030   0.203  ░░░░░░████████  ← FAR END
T=10.0         .019  .033  .033  .034  .038  .021   0.178  ░░░░░░██████
T=12.0         .011  .024  .029  .034  .037  .019   0.155  ░░░░░░░█████
T=15.0         .010  .025  .027  .027  .025  .011   0.124  ░░░░░░░░████
T=20.0         .006  .017  .020  .019  .016  .006   0.084  ░░░░░░░░░░░░  settling
```

The largest displayed pair entries shift from the left toward later pairs.
That is compatible with redistribution under the coupled dynamics. The table
does not identify forward/backward modes, a palindromic physical mirror, an
energy reflection, or a resonator mechanism.

The named N=7 trace has one displayed SumMI peak; the named N=15 trace has
displayed peaks at t=4.5 and t=12.5. Chain length, sampling, initial state,
and modal content are not isolated here, so no reflection mechanism or length
threshold follows.

The threshold-crossing and spatial-MI tables can be compared as temporal and
spatial summaries. Whether one mechanism links them remains open.

Script: `dotnet run -c Release -- wave 7 0.344,0.001,0.001,0.001,0.001,0.001,0.001`

### Hardware (coarse, 5 time points): IBM Torino (March 24, 2026)

The sacrifice-zone experiment on IBM Torino (chain [Q85,Q86,Q87,Q88,Q94],
no DD, 5 time points) supplies five coarse hardware samples. The selected
per-pair MI values below come from 8192 shots per circuit:

```
         Pair:  (0,1)   (1,2)   (2,3)   (3,4)
t=1.0 us       .027    .019    .011    .017     edge dominant
t=3.0 us       .006    .010    .015    .006     center dominant
t=5.0 us       .007    .006    .016    .009     center holds
```

At t=1 pair (0,1) is largest among the four displayed values; at t=3-5 pair
(2,3) is largest. With only five times, this ordering does not establish a
continuous propagation direction or the same mechanism as the simulation.

This data was collected for the DD comparison experiment. The changing
pairwise ordering was not the original goal and is present in the raw counts;
calling it wave propagation would require denser time sampling and a transport
witness.

Data: [sacrifice_zone_hw_no_dd](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_no_dd_20260324_191713.json)

### Finite hardware side record: impedance differences (February 9, 2026)

The impedance ||ZρZ - ρ|| was computed from 25 hardware-measured
density matrices (IBM Torino, Qubit 52, state tomography, 8192 shots).

The displayed impedance values increase with CΨ; along the displayed
decaying-time trajectory they fall with time. The largest displayed
finite-difference magnitude is 0.0076 at the sample closest to the ¼
crossing on this sparse grid:

```
t=74.6 us   CΨ=0.385   |d_imp/dt|=0.0046
t=111.8 us  CΨ=0.261   |d_imp/dt|=0.0076  <-- largest displayed sample (distance from ¼: 0.011)
t=149.1 us  CΨ=0.125   |d_imp/dt|=0.0071
```

These sparse samples do not locate a continuous extremum exactly at ¼
and do not establish fold causality. Reading the outer mirror of the
Fabry-Perot cavity as a switch remains an interpretive invitation, not a
result.

Data: [ibm_impedance_gradient.txt](../simulations/results/ibm_impedance_gradient.txt),
[tomography_ibm_torino_20260209](../data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json)

---

## Pending

- Test whether a same-row ¼/PeakMI coincidence persists at N = 9 and N = 11
- Measure CΨ crossing time resolution (finer than 0.5 intervals)
- Analytical connection: does the discriminant derivative dD/dt predict the PeakMI peak?
- Pure center CΨ diagnostics: does the same coincidence hold?
- Connection to [Proof Monotonicity CΨ](../docs/proofs/PROOF_MONOTONICITY_CPSI.md):
  changing the dissipator puts the sweep outside a fixed-generator claim;
  the sampled CΨ01 rebound is not by itself a theorem violation.

---

## References

- [Resonant Return (formula, position sweep, hybrids)](RESONANT_RETURN.md)
- [Relay Protocol (historical 0.78 heuristic; finite unmatched-endpoint comparison)](RELAY_PROTOCOL.md)
- [Crossing Taxonomy: C(f), two named evolution books, then the finite scalar quarter equation](CROSSING_TAXONOMY.md)
- [Boundary Navigation (θ compass, fold catastrophe)](BOUNDARY_NAVIGATION.md)
- [Mathematical Connections (fold catastrophe, Mandelbrot)](../docs/MATHEMATICAL_CONNECTIONS.md)
- [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md)
- [Spectral Midpoint Hypothesis (dual perspective)](../hypotheses/SPECTRAL_MIDPOINT_HYPOTHESIS.md)
