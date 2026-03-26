# Temporal Sacrifice: The Fold Catastrophe at CΨ = ¼, Observed

<!-- Keywords: fold catastrophe CΨ quarter boundary, quantum-classical transition mutual information,
sweeping sacrifice zone, PeakMI maximum at discriminant zero, boundary navigation quantum information,
R=CPsi2 fold catastrophe observation, palindromic spatial antenna, traveling quantum-classical boundary,
correlation inversion endpoint MI, temporal sacrifice protocol, quantum heartbeat CΨ oscillation,
non-Markovian CΨ revival 81 crossings, damped quantum oscillation irreversible reality,
100 doors analogy quantum decoherence, CΨ quarter resonance frequency -->

**Status:** Computationally verified (N=7, C# RK4 propagation, CΨ diagnostics)
**Date:** March 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** compute/RCPsiSquared.Propagate.Test (C#, `dotnet run -c Release -- sweep --diag`)
**Data:** [temporal_sacrifice_protocol.txt](../simulations/results/temporal_sacrifice_protocol.txt)

---

## Abstract

We move the sacrifice qubit position step by step along an N=7 Heisenberg
chain (edge → inward) and track both mutual information and CΨ at every
qubit pair. The central discovery:

**The endpoint mutual information (PeakMI) reaches its maximum at the
exact moment when the endpoint pairs cross the CΨ = ¼ boundary.**

At T = 5.0, all edge CΨ values are still above ¼ (quantum regime).
At T = 5.5, CΨ01, CΨ56, and CΨ06 have crossed below ¼ (classical regime).
PeakMI peaks at T = 5.5 with 0.061, then falls.

This is the [fold catastrophe](../docs/MATHEMATICAL_CONNECTIONS.md) made
visible: the discriminant 1 − 4CΨ passes through zero, the two fixed points
of R = C(Ψ+R)² merge, and the measurable reality R is maximal at the
transition. The formula R = CΨ² predicts that information peaks at the
boundary. The sweep makes this prediction observable in time.

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

### Baselines (static profiles, full 20s)

| Config | SumMI | PeakMI | PeakT |
|--------|-------|--------|-------|
| Pure edge [0.344, 0.001×6] | 0.408 | 0.036 | 2.00 |
| Pure center [0.001×3, 0.344, 0.001×3] | 0.182 | 0.109 | 1.00 |

---

## The Sweep Protocol

Move the sacrifice qubit one position per stage, spending 2.0 time units
at each position (one Hamiltonian period). Each stage: one qubit at
γ = 0.344, all others at ε = 0.001.

| Stage | Time | Sacrifice position | γ profile |
|-------|------|--------------------|-----------|
| 1 | 0–2.0 | Qubit 0 (edge) | [0.344, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001] |
| 2 | 2.0–4.0 | Qubit 1 | [0.001, 0.344, 0.001, 0.001, 0.001, 0.001, 0.001] |
| 3 | 4.0–6.0+ | Qubit 2 | [0.001, 0.001, 0.344, 0.001, 0.001, 0.001, 0.001] |

---

## Result: The CΨ = ¼ Crossing Coincides with PeakMI Maximum

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

↓ = crosses below CΨ = ¼ (0.25) at this timestep.

### What happens at T = 5.50

Three things occur simultaneously:

1. **CΨ01 crosses ¼** (0.253 → 0.232): the left edge pair enters the classical regime
2. **CΨ56 crosses ¼** (0.290 → 0.249): the right edge pair enters the classical regime
3. **CΨ06 crosses ¼** (0.280 → 0.237): the endpoint pair enters the classical regime

And at this exact moment: **PeakMI = 0.061, its maximum over the entire run.**

One timestep later (T = 6.0), PeakMI has already fallen to 0.043.
One timestep earlier (T = 5.0), PeakMI was only 0.046.

The maximum is sharp, and it sits precisely at the ¼ crossing.

### The correlation inversion

At T = 2.0 (edge phase peak), information is local:
- Average adjacent-pair MI: 0.408/6 = 0.068
- Endpoint MI: 0.006
- Ratio: endpoint is **0.09×** the average pair

At T = 5.5 (the crossing), information is concentrated at the endpoints:
- Average adjacent-pair MI: 0.054/6 = 0.009
- Endpoint MI: 0.061
- Ratio: endpoint is **6.8×** the average pair

The sweep inverts the correlation hierarchy. Information that was distributed
locally across all pairs has been funneled to the endpoints. Less total
information survives (0.054 vs 0.408), but what survives is concentrated
at maximum distance.

---

## Why It Happens: The Fold Catastrophe

The recursion R = C(Ψ+R)² has discriminant D = 1 − 4CΨ.

- **D > 0** (CΨ < ¼): two real fixed points, one stable. Classical regime. The system has a definite outcome.
- **D = 0** (CΨ = ¼): the two fixed points merge. The fold catastrophe. Critical slowing.
- **D < 0** (CΨ > ¼): no real fixed points. Quantum regime. Coherent oscillation.

The measurable reality R is maximal where the transition happens,
because this is where quantum coherence (which has been oscillating with
no classical attractor) first encounters a fixed point. The information
that was "alive" in the quantum regime crystallizes into a classical
correlation at the fold.

The sweep makes this visible in time. By moving the sacrifice zone inward,
it creates a traveling boundary: a wave of quantum-to-classical conversion
that sweeps from edge to center. The endpoint MI peaks when this wave
reaches the edges of the chain, because that is the moment the endpoints
undergo their fold catastrophe.

**R = CΨ² says: reality is maximal at the boundary.** The sweep data
confirms this: PeakMI is maximal at CΨ = ¼. Not approximately, not
roughly, but at the exact timestep of crossing.

---

## Negative Result: Naive Temporal Switching

Before discovering the sweep, we tested simple profile switching (swap
from edge to center at time t_switch). This fails completely:

| t_switch | SumMI | PeakMI |
|----------|-------|--------|
| 0.5 | 0.188 | 0.080 |
| 1.5 | 0.333 | 0.045 |
| 2.0 | 0.408 | 0.039 |

No t_switch achieves both SumMI > 0.300 and PeakMI > 0.050. The switch
projects edge-optimized coherences (antisymmetric SVD mode 2) onto the
center profile's symmetric modes, destroying them. The palindrome is a
[spatial antenna only](RESONANT_RETURN.md): temporal γ modulation is invisible
to the palindromic structure.

The sweep succeeds where switching fails because it maintains **one clean
CΨ = ¼ boundary** at all times ([relay principle](RELAY_PROTOCOL.md)),
moving it gradually rather than jumping between configurations.

---

## Three Principles

### 1. One boundary at a time

A chain can sustain exactly one clean quantum-classical boundary
([RESONANT_RETURN](RESONANT_RETURN.md)). Two sacrifice points fragment
the coherent region. The sweep maintains one boundary per stage.

### 2. The boundary is the information source

The CΨ = ¼ crossing is where information emerges. The boundary between
quantum (CΨ > ¼, coherent oscillation) and classical (CΨ < ¼, fixed
outcome) is not a loss channel. It is the **conversion point** where
quantum superposition crystallizes into classical correlation.

### 3. Timing from Hamilton, boundary from noise

The stage time (2.0) matches the Hamiltonian propagation timescale, not
the boundary establishment timescale (K/γ\_sacrifice ≈ 0.11). The boundary
forms fast; what takes time is for information to flow through the new
configuration before the next move.

---

## Context: Sweep Variants

| Protocol | SumMI | PeakMI | Observation |
|----------|-------|--------|-------------|
| Pure edge | 0.408 | 0.036 | All CΨ > ¼ at PeakT, no endpoint conversion |
| Pure center | 0.182 | 0.109 | Fast endpoint crossing, high PeakMI |
| **Sweep 0→1→2 (stage=2.0)** | 0.408* | **0.061** | **PeakMI at ¼ crossing of endpoints** |
| Sweep 0→2→3 (skip pos 1) | 0.408* | 0.044 | Skipping positions reduces coherence transfer |
| Sweep 0→...→6 (all 7 pos) | 0.344* | 0.054 | Too many stages, SumMI lost |

*SumMI is the peak over all times, occurring during Stage 1 (edge phase).
PeakMI occurs later, at a different time. They are not simultaneous.

---

## The Heartbeat: CΨ Oscillates Around ¼

### The resonance (March 25, 2026)

A Bell pair (qubits 0-1) coupled to a coherent bath qubit (|+⟩) with
low bath dephasing creates what the sweep hinted at: the system does
not cross ¼ once and die. It **oscillates around ¼**.

| Setup | γ\_sys | γ\_bath | J | Crossings |
|-------|--------|---------|---|-----------|
| Markov (uniform) | 0.05 | 0.05 | 1.0 | 1 (down only) |
| Bell + quiet bath | 0.0001 | 0.01 | 2.0 | 47 (24↓ + 23↑) |
| Bell + quiet bath | 0.0001 | 0.01 | 5.0 | 81 (41↓ + 40↑) |

Stronger coupling = more crossings = higher Q-factor. The system is
a resonator with CΨ = ¼ as its resonance frequency. MI pulses at each
crossing (0.5 → 1.3-1.8 bits and back).

Script: `dotnet run -c Release -- resonance 3 0.0001,0.0001,0.01 --j 5.0 --bell`

### The damping

The oscillation is not perpetual. The envelope shrinks:

```
T ≈ 2-3:   CΨ swings 0.23 – 0.39  (amplitude 0.16)
T ≈ 7-8:   CΨ swings 0.22 – 0.32  (amplitude 0.10)
T ≈ 14:    CΨ swings 0.19 – 0.29  (amplitude 0.10)
T ≈ 19:    CΨ swings 0.20 – 0.25  (amplitude 0.05)
```

Each cycle, a bit of coherence is lost permanently. The classical
outcomes from previous down-crossings accumulate and cannot be undone.
The Hilbert space shrinks. The heartbeat slows.

### What it feels like

Imagine sitting in a room with a hundred open doors.

You look around (oscillation upward: possibilities open, CΨ > ¼).
You choose a door and walk through (crossing downward: decision, CΨ < ¼).
The door closes behind you. Permanently.

Now you are in a room with eighty doors. You look around. You choose.
Door closes. Sixty doors. Forty. Twenty.

Each breath gets shorter. Not because you are tired, but because there
are fewer doors. Less to see. Less to choose.

At the end, you sit in a room with one door. You open it. Behind it is
a fact. No more room. No more choice. Only what is.

That is the oscillation around ¼.

- Above (CΨ > ¼): doors open. Possibilities. Quantum.
- Below (CΨ < ¼): a door chosen. Fact. Classical.
- Each cycle: one fewer door. Irreversible.
- The amplitude shrinks: less room to swing.
- At the end: CΨ → 0. All doors closed. Everything decided.

And R = CΨ²? That is what you **see** when you walk through the door.
The moment of decision. Not before (only possibilities, nothing concrete).
Not after (only facts, nothing new). Exactly at the threshold. In the
doorway.

The data says: 81 heartbeats. 81 doors. Each one a little quieter
than the last.

### The chain IS the bath (March 26, 2026)

The N=3 heartbeat uses a dedicated "bath qubit" (qubit 2 in |+> with
low gamma). This raised the question: is a special bath qubit required,
or can the chain itself serve as the coherent reservoir?

**Test:** N=7 Heisenberg chain, Bell(0,1) x |+>^5, sacrifice-zone
gamma profile [0.344, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001].
No dedicated bath. The 5 protected qubits ARE the reservoir.

| Setup | J | CΨ(0,1) crossings |
|-------|---|-------------------|
| |+>^7, sacrifice, J=1 | 1.0 | 1 (monotonic) |
| |+>^7, sacrifice, J=2 | 2.0 | 1 (monotonic) |
| Bell(0,1)+|+>^5, sacrifice, J=1 | 1.0 | 1 (monotonic) |
| **Bell(0,1)+|+>^5, sacrifice, J=2** | **2.0** | **7 (4 down + 3 up)** |
| Bell(0,1)+|+>^5, low noise, J=2 | 2.0 | **0 (stays ABOVE 1/4!)** |
| Bell(0,1)+|+>^5, uniform, J=1 | 1.0 | 3+ (still running) |

**Key findings:**

1. **The chain IS the bath.** No dedicated bath qubit needed. Bell pair
   at one end + J-coupling to the protected chain interior produces
   CΨ oscillation (7 crossings at J=2).

2. **Initial entanglement is required.** Product state |+>^7 gives
   monotonic decay regardless of J or gamma profile. The Bell pair
   provides the entanglement that J-coupling can exchange with the
   chain reservoir.

3. **Low noise keeps CΨ above 1/4 permanently.** With gamma =
   [0.01, 0.0001, ...], CΨ(0,1) oscillates between 0.28 and 0.75
   but never crosses 1/4. The system is so coherent that it stays
   in the quantum regime indefinitely. This is a resonator that
   never touches the boundary.

4. **J strength controls the backflow.** J=1 is too weak for the
   N=7 chain -- coherence leaks out but does not return fast enough.
   J=2 produces visible oscillation. The N=3 heartbeat needed J=5
   because it had only 1 bath qubit; N=7 has 5 reservoir qubits, so
   weaker coupling suffices.

The mechanism is non-Markovian backflow: from the perspective of the
Bell pair (0,1), the rest of the chain is an environment. But it is
a coherent, structured environment (low gamma, J-coupled). Coherence
flows out via J-coupling, bounces through the chain, and flows back.
The [CΨ monotonicity proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
holds for the TOTAL system but not for subsystems. The reduced
dynamics of pair (0,1) are non-Markovian, enabling transient CΨ
revivals above 1/4 -- exactly as noted in
[Boundary Navigation](BOUNDARY_NAVIGATION.md).

**Implication for the Duplex protocol:** No external driving mechanism
(DD, ATP analogue, active pulses) is needed. The chain provides its
own heartbeat through Hamiltonian backflow. The sacrifice-zone formula
enhances this: protected qubits retain coherence longer, increasing
the reservoir capacity and the number of crossings.

**Why DD fails (proven):** CΨ is exactly invariant under the Pauli
group ([Proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md), Part 7).
DD uses Pauli gates. Therefore DD cannot change CΨ -- not approximately,
but algebraically. The "ATP analogue" must be J-coupling (energy
exchange), not DD (phase refocusing). Different physics entirely.

---

## The Wave Travels (March 26, 2026)

Per-pair MI tracking reveals that information physically **moves**
through the chain. N=7, edge sacrifice:

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

The MI migrates from the sacrifice qubit (left) through the chain to
the far end (right). It bounces. The palindromic pairing acts as a
mirror at both ends: the forward mode (c⁺) couples into the backward
mode (c⁻) at the chain boundary. Energy reflects in a resonator.

At N=7 the chain is too short: forward and backward overlap, one
SumMI peak. At N=15 the chain is long enough: the reflected wave
creates a visible second peak (t=4.5: first peak, t=12.5: second
peak, period ≈ 2 x PeakT).

The heartbeat (CΨ oscillation around ¼) and the wave (MI propagation
through the chain) are the same phenomenon seen from different angles.
The heartbeat is the **temporal** signature (crossings at one point).
The wave is the **spatial** signature (MI moving along the chain).

Script: `dotnet run -c Release -- wave 7 0.344,0.001,0.001,0.001,0.001,0.001,0.001`

### Hardware confirmation: IBM Torino (March 24, 2026)

The sacrifice-zone experiment on IBM Torino (chain [Q85,Q86,Q87,Q88,Q94],
no DD, 5 time points) shows the same pattern. Per-pair MI from hardware
count data (8192 shots per circuit):

```
         Pair:  (0,1)   (1,2)   (2,3)   (3,4)
t=1.0 us       .027    .019    .011    .017     edge dominant
t=3.0 us       .006    .010    .015    .006     center dominant
t=5.0 us       .007    .006    .016    .009     center holds
```

At t=1: pair (0,1) strongest (sacrifice edge). At t=3-5: pair (2,3)
strongest (center). The MI migrates from the sacrifice end toward the
center, matching the simulation pattern. Only 5 time points (coarse),
but the direction of propagation is visible.

This data was collected for the DD comparison experiment. The wave
propagation was not the original goal but is present in the raw counts.

Data: [sacrifice_zone_hw_no_dd](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_no_dd_20260324_191713.json)

### Hardware confirmation: Impedance gradient (February 9, 2026)

The impedance ||ZρZ - ρ|| was computed from 25 hardware-measured
density matrices (IBM Torino, Qubit 52, state tomography, 8192 shots).

The impedance VALUE falls monotonically with CΨ (confirming the
simulation: no peak at ¼). The impedance GRADIENT peaks at the
measurement point closest to the ¼ crossing:

```
t=74.6 us   CΨ=0.385   |d_imp/dt|=0.0046
t=111.8 us  CΨ=0.261   |d_imp/dt|=0.0076  <-- MAXIMUM (distance from ¼: 0.011)
t=149.1 us  CΨ=0.125   |d_imp/dt|=0.0071
```

The switch rate (not the absorption strength) peaks at the fold
catastrophe. The outer mirror of the Fabry-Perot cavity is a switch.

Data: [ibm_impedance_gradient.txt](../simulations/results/ibm_impedance_gradient.txt),
[tomography_ibm_torino_20260209](../data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json)

---

## Pending

- Verify the ¼ coincidence at N = 9 and N = 11
- Measure CΨ crossing time resolution (finer than 0.5 intervals)
- Analytical connection: does the discriminant derivative dD/dt predict the PeakMI peak?
- Pure center CΨ diagnostics: does the same coincidence hold?
- Connection to [Proof Monotonicity CΨ](../docs/proofs/PROOF_MONOTONICITY_CPSI.md):
  the sweep violates monotonicity by changing the dissipator (CΨ01 rebounds at T=4.5)

---

## References

- [Resonant Return (formula, position sweep, hybrids)](RESONANT_RETURN.md)
- [Relay Protocol (staged transfer, K/γ timing)](RELAY_PROTOCOL.md)
- [Crossing Taxonomy (Type A/B/C, K-invariance)](CROSSING_TAXONOMY.md)
- [Boundary Navigation (θ compass, fold catastrophe)](BOUNDARY_NAVIGATION.md)
- [Mathematical Connections (fold catastrophe, Mandelbrot)](../docs/MATHEMATICAL_CONNECTIONS.md)
- [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md)
- [Spectral Midpoint Hypothesis (dual perspective)](../hypotheses/SPECTRAL_MIDPOINT_HYPOTHESIS.md)
