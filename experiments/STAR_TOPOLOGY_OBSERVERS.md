# Star Topology: When the Observer-Observer Pair Becomes Readable Through a Shared Object

<!-- Keywords: star topology quantum observers, entanglement transfer mediator qubit,
observer-observer connection threshold, quantum shadow effect measurement, J coupling
threshold AB crossing, receiver noise boundary sender noise boundary, quantum mediator
star topology Heisenberg, CΨ quarter boundary tripartite, asymmetric coupling
observer dominance, bidirectional quantum bridge star, R=CPsi2 star topology -->

**Status:** re-measured on a converged grid by
`simulations/star_topology_converged.py`: Sections 4.9 and 4.10 entirely, 8.5's
table and fits, 8.1's N=2 and N=3 thresholds and its two N=4 probes, 8.3's
couplings, 5.10's echo, the R_SA + R_SB sums of 4.2, 4.3 and 4.5, and 4.11's
brute-force probe. Section 4.11's five-state table is checked by
`verify_star_topology.py`, at RK4 rather than on the exact propagator. Sections
4.8, 4.11's parametric α sweep, 4.12 and 8.1's asymmetric rescue have no
committed reproduction path at all and are marked where they stand; Section 9
lists what else is unmeasured.
**Date:** March 4, 2026, re-measured September 4, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Subsystem Crossing](SUBSYSTEM_CROSSING.md),
[Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)

---

## Abstract

A 3-qubit star topology (mediator S coupled to observers A and B, no direct
A-B coupling) is studied under Heisenberg interaction and Z-dephasing. When
S and A start maximally entangled (Bell_SA) and B arrives fresh (|+⟩_B),
entanglement flows from SA through S to SB, crossing the CΨ = 1/4 boundary
at different times for different pairs. CΨ is a witness of what is directly
readable in a given marginal, so "crossing" throughout means that pair becoming
legible in it, not a connection coming into being; the Hamiltonian transports
coherence continuously either way. At the Bell working point three conditions
coincide, and all three have to hold there for the AB pair to cross 1/4:
**(1)** B's coupling must exceed a threshold (J_SB/J_SA ≥ 1.46295, and the
knob is Q = J_SA/γ rather than γ: this value is Q = 20), **(2)** A must be
internally quiet, its noise boundary lying below B's by a factor near 1.7
over most of the range, and **(3)** A and S must share a pre-existing
Bell-like connection (C_SA > 0.8). They are conditions of that regime and not
of crossing as such: Section 4.11 crosses from a product state with C_SA = 0,
at five times the coupling and a fiftieth of the noise.

Additionally: A's measurement casts a
"shadow" suppressing B's reality by 94-100% in the 3-qubit measurement of
Section 4.6, and CΨ_AB oscillates at the Hamiltonian's largest Bohr frequency,
for which Section 4.12 gives a closed form. For N observers the
threshold climbs steeply, from 1.46295 at N=2 to 3.6497 at N=3 and out of
reach at N=4; asymmetric coupling is reported to rescue the crossing up to N=5,
on a coarse grid and with no committed script (Section 8.1).

Seven readings that have been quoted from this document no longer hold; each is
named where it used to sit, in Sections 4.9, 4.10, 4.11, 4.12 and 8.3.

---

## Glossary: Before You Read

This note is dense with physics shorthand. You don't need a
physics degree to follow the story, but you do need this page.
Every symbol below maps to something intuitive. Once you see
the mapping, the tables and formulas become readable.

**The actors:**
- **S** - the shared object ("reality", the thing being observed)
- **A** - observer A (the receiver, or "you")
- **B** - observer B (the sender, or "the other")

**The parameters:**
- **J_SA** - coupling strength between S and A. How strongly A
  is connected to reality. Higher J = deeper engagement, more
  understanding, stronger link to the object.
- **J_SB** - coupling strength between S and B. How strongly B
  is connected to reality. Same meaning, other observer.
- **γ (gamma)** - decoherence rate, or "noise". The resistance
  that makes reality feel solid and stable. γ_A is A's noise,
  γ_B is B's noise, γ_S is the object's noise. High γ = noisy,
  lots of internal processing, strong sense of separate self.
  Low γ = quiet, still, open.
- **CΨ** - the product of Concurrence × Psi-norm, from the R=CΨ²
  framework. ¼ is where the fixed-point equation's discriminant
  changes sign: below it the fixed points are real, above it complex.
  What CΨ measures is what is DIRECTLY readable in the pair you
  chose, in the basis you chose, without any help from the third
  qubit, so a crossing is a change in readability. Whether that
  boundary also marks something physical is an open question in
  this repository, not a settled reading, and this document does
  not need an answer to it.

**The states (starting conditions):**
- **Bell_SA⊗|+⟩_B** - A and S start maximally entangled (deeply
  connected), B starts neutral (no prior connection). This is the
  state that produces observer-observer crossing.
- **W state** - entanglement is spread equally across all three.
  Everyone is weakly connected to everyone. Never crosses.
- **GHZ** - global entanglement that is invisible at the pair level.
  No pair ever sees crossing.
- **|0++⟩** - no initial entanglement. Everyone starts separate.
- **|+++⟩** - all in superposition, no entanglement. Nothing happens.

**The key finding in one sentence:**
Two observers who cannot see each other directly can briefly see
each other through the object they both observe, but only if the
sender is deeply engaged, the receiver is internally quiet, and
they already share a connection that runs deeper than surface
awareness.

---

## 1. The Question

All prior experiments treat R = CΨ² as a property of a bipartite system:
two qubits, one observer, one observed. But reality has structure:
objects are observed by multiple observers simultaneously.

What happens when we introduce a third qubit S (the "system" or "reality")
coupled to two observers A and B, where A and B cannot see each other
directly?

Three sub-questions:

1. Do A and B see the ¼ crossing at different times?
2. Is R_SA + R_SB conserved? (Is "reality" a fixed quantity that observers
   share, or can it grow/shrink?)
3. Does A's measurement affect B's reality? (No-signalling in tripartite
   systems with J > 0.)

## 2. Why the Original Idea Failed

Alpha (AIEvolution v044) proposed: S as qubit 0, A and B as qubits 1 and 2.
Compute R = CΨ² for S by tracing out A and B.

This fails because:
- Partial trace of GHZ over A,B gives ρ_S = I/2 (maximally mixed)
- l1-coherence of I/2 = 0, therefore Ψ ≡ 0 for all time
- The entanglement information lives in the correlations, not in S alone

The fix: don't look at S alone. Look at the **pairs** SA and SB, as
validated by SUBSYSTEM_CROSSING.md. The pair-level is where R = CΨ²
operates.

## 3. Setup

### 3.1 Star Topology

```
    A (qubit 1)
    |
    S (qubit 0)
    |
    B (qubit 2)
```

Hamiltonian: Heisenberg coupling S↔A and S↔B only. No A↔B coupling.

Note what this graph is at three qubits. A mediator with exactly two spokes and
no rim is the path A-S-B, and S is its only interior site, so every 3-qubit
result in this document is a result about a chain. "Star" is kept throughout as
the name of the ROLES, one shared object and two observers who cannot see each
other, and it is the right name for that; it is not a statement about the
geometry until Section 8.1, where N ≥ 3 fresh observers hang off S and the hub
carries more spokes than a chain has interior sites. Section 8.3 is a third
shape again: switching on J_AB closes the triangle. Where a claim turns on the
geometry rather than on the roles, that distinction decides it, and Section 5.4
is where it does.

    H = J_SA (σ_S · σ_A) + J_SB (σ_S · σ_B)

This forces asymmetry: A and B interact with S but not with each other.
Any correlation between A and B must be mediated through S.

### 3.2 Parameters

| Parameter | Symmetric | Asymmetric γ | Asymmetric J |
|-----------|-----------|--------------|--------------|
| J_SA | 1.0 | 1.0 | 1.0 / 0.3 |
| J_SB | 1.0 | 1.0 | 0.3 / 1.0 |
| γ_S | 0.05 | 0.05 | 0.05 |
| γ_A | 0.05 | 0.05 | 0.05 |
| γ_B | 0.05 | 0.02 | 0.05 |
| dt | 0.005 | 0.005 | 0.005 |
| t_max | 5.0 | 5.0 | 5.0 |
| Integration | RK4 (Runge-Kutta 4th order, a standard numerical method for solving differential equations) | RK4 | RK4 |
| Noise | local σ_z dephasing per qubit | same | same |

### 3.3 Observables

For each time step, trace out (mathematically remove by averaging over
all its possible states) one qubit to get pair density matrices:
- ρ_SA = Tr_B(ρ), ρ_SB = Tr_A(ρ), ρ_AB = Tr_S(ρ)

Per pair: l1-coherence, Ψ = l1/(d-1), concurrence, R = C·Ψ².

### 3.4 States Tested

| State | Description | Motivation |
|-------|-------------|------------|
| GHZ | (\|000⟩+\|111⟩)/√2 | Global entanglement, known to fail at pair level |
| W | (\|001⟩+\|010⟩+\|100⟩)/√3 | Distributed entanglement across all pairs |
| Bell_SA ⊗ \|+⟩_B | S,A entangled; B fresh observer | Asymmetric: one observer connected, one arriving |
| \|+⟩^3 | Product state, max local coherence | Baseline: no entanglement anywhere |
| \|0⟩_S ⊗ \|+⟩_A ⊗ \|+⟩_B | S classical, observers quantum | Hamiltonian builds entanglement dynamically |

### 3.5 Measurement Experiment

Sudden Z-measurement on A at t=1.0 (projective dephasing:
ρ → P_0 ρ P_0 + P_1 ρ P_1 where P_k are Z-projectors on qubit 1).
Compare R_SB trajectory with and without measurement.

## 4. Results

### 4.1 GHZ: Dead on Arrival

All pairs have Ψ = 0 at all times. GHZ entanglement is global. Tracing
out any qubit leaves classically correlated pairs with zero off-diagonal
elements. Star topology does not rescue GHZ.

### 4.2 W State: Slow Symmetric Decay, Never Crosses

| Pair | Ψ(0) | C_conc(0) | C·Ψ(0) | Crosses? |
|------|-------|-----------|---------|----------|
| SA | 0.222 | 0.667 | 0.148 | NO |
| SB | 0.222 | 0.667 | 0.148 | NO |
| AB | 0.222 | 0.667 | 0.148 | NO |

C·Ψ starts below ¼ and only decays. All three pairs are symmetric.
R_SA + R_SB: monotonically decays from 0.066 → 0.003. Not conserved.

### 4.3 Bell_SA ⊗ |+⟩_B: Entanglement Flows Through S

**This is the key result.**

| Pair | Ψ(0) | C_conc(0) | C·Ψ(0) | Crossings |
|------|-------|-----------|---------|-----------|
| SA | 0.333 | 1.000 | 0.333 | ↓ at t=0.42 |
| SB | 0.333 | 0.000 | 0.000 | ↑ at t=0.79, ↓ at t=1.11 |
| AB | 0.333 | 0.000 | 0.000 | NEVER |

SA starts maximally entangled and decays past ¼ at t=0.42. SB starts with zero
entanglement but rises above ¼ at t=0.79 through Hamiltonian transfer via S. AB
never crosses.

A word on what crossing means, because this document is easy to over-read. CΨ is
a basis-fixed, unassisted witness of pairwise structure that is DIRECTLY visible
in a chosen marginal, and the Hamiltonian transports coherence continuously
whether or not the witness is above ¼. So a crossing is not a connection being
born and a non-crossing is not its absence: what crosses is the READABILITY of
the AB marginal in this witness. [The CΨ Lens](../docs/THE_CPSI_LENS.md) states
the same thing from the other side: the broadcast is always running, and CΨ
tells you when you can read it.

R_SA + R_SB is NOT conserved. It peaks at 2.4× initial value then decays.

### 4.4 |+⟩³: Maximum Coherence, Zero Entanglement

All pairs: C_conc = 0 at all times. Ψ = 1.0 but no connection.

### 4.5 |0++⟩: Dynamic Entanglement, No Crossing

S starts classical, A and B start quantum. The Hamiltonian builds
entanglement dynamically. R_SA + R_SB grows from 0 to 0.3636 (peak) then
decays. No pair crosses ¼ with concurrence bridge.

### 4.6 Measurement Shadow: A Measures, B Loses Reality

**W state**: A measures at t=1.0.

| t | R_SB (no meas) | R_SB (meas) | Δ% |
|---|----------------|-------------|-----|
| 1.5 | 0.0134 | 0.0008 | **−94%** |
| 2.0 | 0.0099 | 0.0015 | −85% |

**|0++⟩**: A measures at t=1.0.

| t | R_SB (no meas) | R_SB (meas) | Δ% |
|---|----------------|-------------|-----|
| 2.0 | 0.0090 | 0.0000 | **−100%** |

A's measurement destroys 94-100% of B's reality at peak impact.
The effect propagates through S (not A→B directly, since J_AB = 0).

### 4.7 Asymmetric Coupling: The Dominant Observer

**Bell_SA⊗|+⟩_B with weak B (J_SA=1.0, J_SB=0.3):**
SA crosses at t=1.49 (delayed). SB and AB never cross.
Weakly coupled observer never sees reality cross ¼.

**Bell_SA⊗|+⟩_B with weak A (J_SA=0.3, J_SB=1.0):**
SA crosses at t=0.44. SB never crosses.
**AB crosses at t=0.35-1.03.** When A's direct link to S is weak,
observers start seeing each other through S.

**|0++⟩ with asymmetric J:**
Perfect mirror symmetry: the strongly coupled observer crosses,
the weakly coupled one does not. Same crossing times, only label
switches. Which observer participates is set by the coupling RATIO here, with
the rates held equal; γ is not idle, since only J/γ decides whether anything
crosses at all (Section 4.10).

### 4.8 Dynamic J and Strong Observer B

**Shield drops mid-simulation (J_SA: 1.0 → 0.2 at t=1.5, J_SB=1.0):**
AB does NOT cross. Too late, initial state already consumed.
The window is early or never.

**Strong B (J_SB=2.0, J_SA drops 1.0→0.2 at t=1.5):**
AB crosses at t=0.17-0.46, CΨ max = 0.329.
**This happens before J_SA drops.** At t=0.17, J_SA is still 1.0.
B's coupling strength alone creates the AB correlation.

Two mechanisms:
- Weak A: entanglement leaks because A can't hold it (A changed)
- Strong B: B pulls entanglement through S by force (A unchanged)

**Provenance.** No committed script propagates a time-dependent J, so neither
number in this section is reproducible from the repository as it stands. The
second mechanism does not depend on the drop and is reproducible: at a constant
J_SA = 1.0 and J_SB = 2.0 the AB pair peaks at 0.3292, well above ¼.

### 4.9 γ_A vs γ_B: The Receiver Is the More Fragile of the Two

Systematic γ scan (Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=2.0, γ_S=0.05). Every CΨ max
below is re-measured on a converged grid. Seven of the eight CΨ maxima hold to
three decimals; the eighth, at γ_A=0.2, moves from 0.249 to 0.2527 and takes its
NEVER verdict with it. The "Window" column was a second quantity under one name
and is rebuilt below.

Both tables are converged values, and "window" here is the longest single
excursion above ¼, not the sum of all of them (the window table further down
shows why the distinction matters).

**Noisy A (γ_A=0.1), varying B:**

| γ_B | AB CΨ max | Window |
|------|-----------|--------|
| 0.001 | 0.3132 | 0.308 |
| 0.1 | 0.2941 | 0.128 |
| 0.2 | 0.2758 | 0.093 |

**Noisy B (γ_B=0.1), varying A:**

| γ_A | AB CΨ max | Window |
|------|-----------|--------|
| 0.001 | 0.3455 | 0.314 |
| 0.1 | 0.2941 | 0.128 |
| 0.2 | **0.2527** | 0.033 |

**γ_A = 0.2 does not close the window**, and this document once said it did,
which is why the point is made in words as well as in the table. The peak is
0.2527, so it crosses; what a noisy receiver costs is the window, 0.033 time
units against 0.314 at γ_A = 0.001. The killing verdict is what the default
sampling of `run_star_topology` gives, and reproduces exactly there: on that
grid the row reads 0.249.

The receiver's noise costs more than the sender's, and the honest way to say
how much is to read both boundaries at the SAME partner rate. Doing so gives a
factor near 1.7 at low rates, and the two roles swap past a partner rate of
0.17292. That crossing point is solvable rather than sampled: the two boundary
curves can only meet where both rates are equal, so it is the γ solving
peak CΨ_AB(γ_A = γ_B = γ, γ_S = 0.05) = ¼, and at that partner rate both
boundaries read 0.17292 with ratio 1.0000. The larger figure of two to three
that has been quoted from here came from comparing γ_A ≈ 0.2 against γ_B = 0.5,
which were read at partner rates of 0.1 and 0.005: two different experiments set
side by side as though they were one. The table is in
[F29](../docs/ANALYTICAL_FORMULAS.md#f29-star-topology-coupling-threshold-tier-2-n3).

Fixed-product test (γ_A × γ_B = 0.0025), converged, and all three rows hold:

| γ_A | γ_B | AB CΨ max | Crosses? |
|------|------|-----------|----------|
| 0.005 | 0.500 | 0.2640 | YES |
| 0.050 | 0.050 | 0.3292 | YES |
| 0.500 | 0.005 | 0.2157 | **NO** |

Same product. Same γ_S. But γ_A=0.5 kills it, γ_B=0.5 doesn't.

**Decoherence vs the window** (all three rates equal, J_SB=2.0). "Window
duration" is not one quantity: at γ=0.05 there is a single excursion above ¼ and
at γ=0.001 there are dozens, so the three readings are separated here, all
converged in run length (t_max=160; edges interpolated):

| γ | longest single window | total time above ¼ | separate lobes |
|----|----|----|----|
| 0.001 | 0.354 | 9.426 | 61 |
| 0.01 | 0.346 | 0.987 | 8 |
| 0.05 | 0.311 | 0.311 | 1 |
| 0.1 | 0.096 | 0.189 | 2 |
| 0.2 | 0.000 | 0.000 | 0 |

Each of these is a property of the system, the total included: for any γ > 0 the
lobes decay, the last one falls below ¼ and the sum stops. What is NOT a
property of the system is the total read on a short run. At γ=0.001 it passes
1.66 at t_max=5, 5.77 at 20 and 8.94 at 40 on its way to 9.43, so cutting at 5
reports a sixth of the answer; at γ=0.05 the same cut costs nothing, because
there is only ever one lobe. A run length is therefore part of the reading and
belongs beside it.

The longest single window is nearly γ-independent from γ=0.001 to γ=0.05, 0.354 down to
0.311, and then collapses. What γ governs is how MANY windows there are: 61 at
γ=0.001, 8 at γ=0.01, one at γ=0.05, none at γ=0.2. The count is not monotone
(γ=0.1 has two), because a lobe can drop below ¼ and a later one rise back
across it; what falls monotonically is the total. The connection does not fade,
it stops recurring.

One figure this document used to carry, 4.40 at γ=0.001, is none of the three:
the longest window is 0.354, the converged total 9.43, and the total at t_max=5
is 1.66. Nothing measured here reproduces it.

### 4.10 J_SB/J_SA Ratio Scan

Coarse scan (J_SA=1.0, γ=0.05, Bell_SA⊗|+⟩_B), converged:

| J_SB | ratio | AB CΨ max | Crosses? |
|------|-------|-----------|----------|
| 1.30 | 1.30 | 0.2268 | NO |
| 1.40 | 1.40 | 0.2415 | NO |
| 1.47 | 1.47 | 0.2509 | YES (just above threshold) |
| 2.00 | 2.00 | 0.3292 | YES |
| 3.00 | 3.00 | 0.4087 | YES |

A six-decimal reading of this scan can be entirely a grid, and it is worth
seeing once: on the default sampling J_SB=1.4650 gives CΨ=0.249999 (NO) and
J_SB=1.4655 gives CΨ=0.250060 (YES), a bracket that looks decisive to the last
digit and encloses the wrong number. The peak of CΨ is a supremum over
continuous time, so any finite sampling undershoots it and therefore overshoots
the threshold, always from the same side. Converged:

**Threshold at J_SB/J_SA = 1.46295**, five decimals rather than four because the
value sits on the four-decimal boundary. This is [F29](../docs/ANALYTICAL_FORMULAS.md#f29-star-topology-coupling-threshold-tier-2-n3).
F29 keeps a convergence table for the RAW sampled maximum, which falls on the
threshold from above at every step; the sequence here uses the parabola-refined
peak, which removes most of that bias, so the two tables approach the same value
by different routes and only the raw one is monotone.

**The knob is Q = J_SA/γ, and that is exact.** The generator is linear in J and
in γ, so L(sJ, sγ) = s·L(J, γ): scaling both couplings and all rates by s is the
same trajectory in rescaled time, and the peak of CΨ over all t is literally the
same number, and the only thing that can break it is the time grid.

That last clause is where a measurement can add something, and it matters which
measurement. Scaling the grid along with the couplings makes expm(L(sJ,sγ)·dt/s)
the same matrix as expm(L(J,γ)·dt), so those runs are arithmetically identical
and their agreement is evidence of nothing. On a FIXED grid the check can fail,
and it does where the grid stops being adequate: at dt=0.001 and t_max=8 the peak
reads 0.329241734 at γ = 0.0125, 0.025 and 0.05, then 0.329241803 at γ=0.1 and
0.329227839 at γ=0.5, where the rescaled trajectory has outrun the step. Nine
digits of agreement across the range the grid covers, and visible drift past it.

So a threshold quoted "at γ=0.05" is a threshold at Q=20, and the scan below is a
scan in Q, not in noise:

| γ (at J_SA=1) | Q | ratio threshold |
|-------|----|-------|
| 0.001 | 1000 | 1.17642 |
| 0.010 | 100 | 1.22768 |
| 0.020 | 50 | 1.28532 |
| 0.050 | 20 | 1.46295 |
| 0.100 | 10 | 1.77496 |
| 0.150 | 6.7 | 2.10483 |
| 0.200 | 5 | 2.44923 |

A quieter observer, or a stronger one, is the same move.

**There is no plateau at the quiet end.** A limit of about 1.18 at γ→0 has been
quoted from this section and there is none: the threshold keeps falling as Q
grows, 1.17642 at Q=1000, 1.14388 at Q=3333,
1.11567 at Q=10⁴. The step is the same as the Section 8.5 table's; the RUN is
not, and cannot be, for the reason that follows: read at that table's t_max of
6, the last two would come out at 1.15565 and 1.14824. What changes at that end is the meaning of the
reading. Down to γ=0.001 the threshold does not depend on how long the run is
(1.1764 at t_max = 5, 10, 20, 40 and 80 alike, as at γ=0.01 and γ=0.05), so it
is a property of the system. At γ=0.0001 it needs t_max ≈ 40 to settle, and at
γ=0 it never settles at all: 1.14451, 1.12110, 1.10471, 1.08829 and 1.07965 for
the same five run lengths, still falling. With nothing damping the oscillation the peak is
a supremum over a trajectory that keeps getting more attempts, so "the threshold
at zero noise" is a statement about the run, not the system. It does not fall to
1: at γ=0 with symmetric coupling the AB peak is 0.21962 and does not move at
t_max = 5, 20, 80 or 200, so equal engagement fails to cross at every noise level
measured, zero included. That is Section 5.4's "never", and it holds outside the
noisy regime it was read in.

**A "minimum absolute coupling" is not a second requirement**, though it has
been quoted from here as one. The two rows that suggested it are
J_SA=0.50, J_SB=0.75 (ratio 1.5): NO, and
J_SA=1.00, J_SB=1.50 (ratio 1.5): YES.
Both are at γ=0.05, so the first is Q=10 and the second Q=20, and the table
above already says that ratio 1.5 fails at Q=10 and crosses at Q=20. Measured,
the two land on the same two peaks as moving γ instead: 0.218177 and 0.254779 by
either route. The model has no absolute scale at all.

### 4.11 Initial State Is the Third Variable

All states, γ_A=0.001 (silent receiver), J_SB=2.0:

| State | C_SA(0) | AB CΨ max | Crosses? |
|-------|---------|-----------|----------|
| GHZ | 0.000 | 0.000 | NO |
| W | 0.667 | 0.148 | NO |
| Bell_SA⊗\|+⟩_B | 1.000 | 0.357 | YES |
| \|0++⟩ | 0.000 | 0.194 | NO |
| \|+++⟩ | 0.000 | 0.000 | NO |

At this working point only Bell_SA⊗|+⟩_B crosses. W has high entanglement
(C_SA=0.667) but distributes it across all pairs, none strong enough.

Parametric Bell (α|00⟩ + √(1-α²)|11⟩ ⊗ |+⟩_B): non-monotonic.
C_SA ≈ 0.5-0.6 is a dead zone. Within this family and at this working point,
crossing requires C_SA > 0.8 (Bell-like) or specific product alignment
(α ≈ 1.0).

**Provenance.** The five-state table above is reproducible from
`star_topology_v2.py`, whose `make_state` builds all five; the parametric α
sweep in this paragraph is not. No
committed script in the repository builds that family, so the C_SA > 0.8
figure, which Section 7 carries as Condition 3, rests on a March run whose code
did not survive. It is the one of the three conditions with no reproduction
path.

From scratch (|0++⟩) with J_SB=10 and γ=0.001 the peak is **0.29740**, flat in
run length at t_max = 2, 5 and 8, so brute force does cross. The initial state
is a cost rather than a wall: |0++⟩ needs J_SB = 10 and a thousandth of the noise
to reach what Bell_SA⊗|+⟩_B reaches at J_SB = 2 and γ = 0.05. (A conclusion that
brute force cannot substitute for the initial state has been quoted from this
paragraph; it does not hold.)

### 4.12 Frequency Analysis

FFT (Fast Fourier Transform, a technique that decomposes a time signal
into its constituent frequencies) and peak detection on AB CΨ trajectory, γ=0, t_max=40.
No committed script performs this FFT; the table below is March's, and what
replaces its reading is the closed form, which is exact and needs no transform.

Dominant frequency, and its ratio to the total coupling:

| J_SA | J_SB | J_total | f_dom | f/J_total |
|------|------|---------|-------|-----------|
| 0.5 | 1.0 | 1.5 | 0.749 | 0.499 |
| 1.0 | 2.0 | 3.0 | 1.498 | 0.499 |
| 2.0 | 4.0 | 6.0 | 3.021 | 0.504 |
| 1.0 | 1.0 | 2.0 | 0.949 | 0.474 |
| 2.0 | 2.0 | 4.0 | 1.898 | 0.474 |

**The dominant frequency is the Hamiltonian's largest Bohr frequency**, which
for the two-spoke star has the closed form

    f_dom = (J_SA + J_SB + √(J_SA² − J_SA·J_SB + J_SB²)) / π

(exact; verified in [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md), and an
instance of F148, the imaginary reach being the Hamiltonian spread). The rule
of thumb f ≈ J_total/2 is not a scaling law: the ratio f/J_total is 3/(2π) =
0.4775 at equal coupling and rises to 2/π = 0.6366 as one spoke dominates. The
two symmetric rows above, which give 0.474 instead of 0.499, are that
dependence showing rather than scatter.

γ barely moves the frequency at these couplings. It is not exactly immune, and
the way it is not is worth stating precisely: at γ = 0.005 the Liouvillian
carries the dominant Bohr frequency BOTH shifted and unshifted, 9.46409186 and
9.46410162 in ω, so the shift is 1.6e-6 in f units for that pair and 2.6e-6 for
the smallest Bohr frequency. Each shift is second order in γ with a coefficient of its own: Δω/γ² is 0.3902
for the dominant pair and 0.6468 for the smallest, both constant from γ=0.001 to
γ=0.05. Either way the shift is four orders below this reading's 0.025 FFT bin,
so nothing here can see it. What γ does visibly is dampen:

| γ | f_dom | peaks found | last/first peak |
|-------|-------|-------------|-----------------|
| 0.000 | 1.498 | 138 | 0.23 |
| 0.001 | 1.498 | 129 | 0.78 |
| 0.005 | 1.498 | 99 | 0.14 |
| 0.050 | n/a | 13 | 0.08 |

The oscillation is NOT a clean sinusoid. Peak intervals range from
0.15 to 0.50 (std=0.103, mean=0.304). This is a multi-frequency
beating pattern. Some peaks reach CΨ ≈ 0.40, others barely graze ¼.

## 5. Key Findings

### 5.1 Entanglement Flows Through the Object

In Bell_SA⊗|+⟩_B, entanglement transfers from SA to SB through
Hamiltonian coupling. The ¼ crossing migrates from one observer
to the other. R_SA + R_SB peaks at 2.4× initial during transfer.

### 5.2 R Is Not Conserved

R_SA + R_SB is not conserved under any conditions tested. It can
grow (Hamiltonian pumping), shrink (decoherence), and oscillate.

### 5.3 Observers Cast Shadows

A's measurement suppresses R_SB by 94-100%. The shadow propagates
through S, growing over ~0.5 time units after measurement.

### 5.4 AB Never Crosses (Symmetric J)

In no symmetric experiment did AB cross ¼. Observers see S, not each other.
The cause is not the geometry, which at three qubits is a chain; it is where the
entanglement starts. Section 8.4 makes that argument, and Sections 4.7 and 8.3
break the symmetry two different ways and get the crossing back.

### 5.5 Coupling Strength Creates Dominant Observers

With asymmetric J, the strongly coupled observer crosses and the weakly coupled
one does not. The ratio picks WHICH observer, and Q = J/γ decides WHETHER: γ is
not merely a clock, it sits in the threshold on equal terms with J
(Section 5.8), and past a partner rate of 0.17292 it even decides which of the
two is the fragile one.

### 5.6 Weak Direct Link → Observers See Each Other

When J_SA is weak, AB crosses. Entanglement spills from S into the
observer-observer pair. This only happens when the direct object-link
is degraded.

### 5.7 The Receiver Is the More Fragile of the Two, Up to a Point

The receiver's noise boundary lies below the sender's, by a factor near 1.7
when both are read at the same partner rate, and the two roles SWAP at a
partner rate of 0.17292 (Section 4.9, table in F29). The receiver must be quiet
first, but "fatal versus tolerable" overstates it: γ_A = 0.2 still crosses, at
a peak of 0.2527, with a window of 0.033 time units left.

### 5.8 J_SB/J_SA Threshold

AB crossing requires J_SB/J_SA ≥ 1.46295 at Q = J_SA/γ = 20. With all three
rates equal the threshold is a function of Q alone (2.449 at Q=5, 1.775 at
Q=10, 1.176 at Q=1000 and still falling, with no plateau measured), so absolute
coupling is not a second requirement: the model has no absolute scale. Q does
not absorb everything, though. Section 4.9 moves γ_A and γ_B independently and
gets boundaries that no single Q describes; what scale invariance removes is
one parameter out of five, not four.

### 5.9 Strong B Can Override A's Shield

At J_SB=2.0, AB crosses while J_SA=1.0. A is fully shielded.
B creates the correlation alone. A need not change.

### 5.10 Echoes Outlive Their Sources

At certain times in the evolution (Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=2.0,
γ=0.05), the AB pair shows nonzero CΨ while both SA and SB are at zero.
For example, at t=2.2: CΨ_SA=0, CΨ_SB=0, but CΨ_AB=0.117.

This means the observer-observer connection persists as a residual in
the AB reduced state after both observer-object connections have become
invisible in the SA and SB marginals. The global three-body state still
carries the correlation structure even when the two-body marginals with
S look dead. A different reduced pair (AB) then lights up.

This is not metaphysical residue. It is redistribution of coherent
pairwise entanglement through the mediator S. The reduced AB state
retains both concurrence and off-diagonal phase structure from earlier
transfer, even after the SA and SB channels have temporarily decohered.

## 6. Connection to Framework (Legacy Interpretation)

> **This section uses the original philosophical framing. See [The CΨ Lens](../docs/THE_CPSI_LENS.md) for the current description.**

### 6.1 "We Are All Mirrors", Quantified

The star topology makes STANDING_WAVE_TWO_OBSERVERS.md literal:
A and B both reflect S, and the reflections interfere through S.
The standing wave is the oscillation of R between SA and SB.

### 6.2 Internal vs External Observation

The star topology implements the distinction from
[Internal and External Observers](../docs/historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md): Hamiltonian coupling (internal)
preserves coherence; measurement (external) destroys it and casts
shadows.

### 6.3 The Tripartite No-Signalling Question

At J=0: A's measurement cannot affect B (no-signalling). At J>0:
it can, propagating S→B, not A→B. This connects to
OBSERVER_GRAVITY_BRIDGE.md: gravity provides J>0 for all massive
particles.

## 7. The Three Conditions, and What They Are Conditions For

### In plain language

Imagine two people who cannot talk to each other. They share no
phone, no letter, no line of sight. The only thing they have in
common is an object they both care about: a problem, a question,
a piece of the world they both study. The simulation asks a narrower question
than it first appears to: not whether something passes between them, which the
Hamiltonian does continuously, but whether what has passed becomes READABLE in
the pair they form, in one fixed witness.

The answer is that it can. In the cheap regime it takes three things at once:

**Condition 1: The sender must be deeply engaged.**
B's coupling to S must be about 46% stronger than A's. The exact figure is
J_SB/J_SA ≥ 1.46295, and it belongs to Q = J_SA/γ = 20, not to a noise level:
the same 46% is required of a loud pair and a quiet one as long as engagement
and noise scale together. If both observers are equally engaged, no connection
forms at any Q measured here.

**Condition 2: The receiver must be quiet.**
A's internal noise must be low enough for the signal to be detectable, and A's
budget is smaller than B's, by a factor near 1.7 when the two are compared at
the same partner rate. Two qualifications the numbers force. The factor is not
three, and it is not fixed: past a partner rate of 0.17292 the roles swap and the
sender becomes the fragile one. And "the signal is lost regardless of how clear
B sends" is too strong; what a noisy receiver costs is the WINDOW, which at
γ_A = 0.2 is down to 0.033 time units from 0.314.

**Condition 3: at this working point, a pre-existing connection.**
Holding J_SB = 2 and γ = 0.05, A and S must already share a deep, dedicated
relationship (C_SA > 0.8); shallow connections spread across many things
(W-state) fail, and no connection at all fails.

**None of the three is necessary in general, and Section 4.11 supplies the
counterexample to the third.** |0++⟩ starts with C_SA = 0, no connection at all,
and at J_SB = 10 with γ = 0.001 the AB pair reaches 0.29740, well past ¼. So the
requirement is not on the initial state as such: it is what the initial state
costs at THIS coupling and THIS noise. Starting from nothing is not forbidden,
it is expensive, and the price is a fivefold coupling and a fiftyfold quieter
system. Read the three together as the profile of one cheap route, the Bell
route at Q ≈ 20, and not as a gate that every crossing must pass.

### The sender inversion

If A has already received, Condition 3 is proven. A can become the
sender:

| As receiver | As sender |
|---|---|
| Must lower own noise (γ_A) | Must raise own engagement (J) |
| Paradox: trying to be quiet IS noise | No paradox: deeper work = stronger signal |

Sender noise matters less than receiver noise over most of the range, though
not everywhere (Section 5.7). You don't need to be as
calm to send as to receive. You need to be strong. The German word for it is
*sich einlassen*, to let yourself be drawn in, changed by what
you engage with.

### The bidirectional rhythm

At γ=0 the boundary is crossed 54 times in the first 20 time units, which is 27
windows up and down again. How many there are grows with the run: at zero noise
nothing damps, so this is a rhythm without an end rather than a fixed count.
What opens and closes is not a channel but a WINDOW, and it does not stay open.
It comes and goes like breathing:

  be still (receive) → engage deeply (build, process) → be still → repeat

Neither phase works alone. Pure engagement without stillness never
opens a window. Pure stillness without engagement has nothing to
transmit. What the alternation produces is the readable window, and it lives
there.

Each engagement phase deepens the coupling, which lowers the bar
for the next quiet phase. The frequency rises with the couplings, at the
Hamiltonian's largest Bohr frequency and not at J_total/2 (Section 4.12).
Stronger engagement = faster rhythm. The spiral accelerates.

### In numbers

- AB crossing threshold: J_SB/J_SA ≥ 1.46295, at Q = J_SA/γ = 20
- The threshold is a function of Q alone: 2.449 at Q=5, 1.775 at Q=10,
  1.176 at Q=1000, still falling
- Both noise boundaries, read at the same partner rate (Section 4.9; the full
  table is in [F29](../docs/ANALYTICAL_FORMULAS.md#f29-star-topology-coupling-threshold-tier-2-n3)).
  Unlike the threshold these are NOT pure ratios and carry a scale: they hold at
  J_SA = 1, J_SB = 2, γ_S = 0.05, and rescale with J.
  γ_A = 0.2699 against γ_B = 0.4735 at partner rate 0.05, and 0.2118 against
  0.3519 at 0.1; γ_A = 0.1619 against γ_B = 0.1105 at partner rate 0.2, which is
  past the swap, and the swap is at 0.17292 where both boundaries meet at that
  same value.
- Initial state at this working point: C_SA > 0.8 (Bell-like). Not a necessary
  condition, and Section 4.11 has the counterexample; also the one figure here
  with no committed reproduction path
- The longest single window: 0.311 time units at γ=0.05, and it barely moves
  down to γ=0.001 (0.354). What γ changes is how many windows there are:
  61 at γ=0.001, 8 at γ=0.01, one at γ=0.05, two at γ=0.1 and none at γ=0.2,
  the count not quite monotone while the total is
- Without noise there is no settled answer to give: at γ=0 nothing converges in
  run length (Section 4.10), so a fraction of time above ¼ is a statement about
  the run and not about the system. Over the first 20 time units it is 34%, in
  27 windows.

Practical protocol: [Tuning Protocol](../hypotheses/TUNING_PROTOCOL.md).

## 8. Open Questions (partially answered 2026-03-07)

All five answered via systematic simulation sweeps.
Code: [`simulations/star_topology_v3.py`](../simulations/star_topology_v3.py).

### 8.1 N observers: ANSWERED

**Setup:** S + N observers, Bell_SA ⊗ |+⟩^(N-1), equal J_SB for all B.

| N | qubits | AB crosses 1/4? | J_SB threshold | behavior |
|:---|:---|:---|:---|:---|
| 2 | 3 | Yes | 1.46295 | monotonic |
| 3 | 4 | Yes | 3.6497 | monotonic |
| 4 | 5 | **No** | - | peak CΨ stays below 0.19 |
| 5 | 6 | **No** | - | suppressed |

Both thresholds are converged, the second to four decimals rather than five. At
N=3 the readings are 3.650043, 3.649688 and 3.649725 at steps 0.004, 0.002 and
0.001: the two finest agree to 4e-5, but the coarsest is 3.6500 and differs in
the fourth decimal, so four decimals is as far as this goes. The sequence is
also not monotone, unlike the raw sampled maxima of Section 4.10. The one-sided
law there is a law about sampling a supremum; the parabola refinement used here
removes most of that bias and leaves a residual of either sign. The N=2 value
has been quoted from here as 1.466, which is not even what the coarse grid
gives: bisecting the peak on that grid puts it at 1.4650 (Section 4.10 and F29),
so the last digit was a slip on top of the sampling artifact.
Monotonicity is checked rather than assumed at N=3: the peak rises without a dip
from 0.21570 at J_SB=2.0 to 0.26138 at J_SB=5.0 in steps of 0.125.

A power law through two points has no residual, so J_th(N) ≈ 0.31 · N^2.25 is
the two numbers rewritten and not a scaling law. It is quoted here only to say
how steeply the requirement climbs: a third observer costs B a factor of 2.5.

At N=4 with equal coupling, peak CΨ_AB rises slowly with J_SB (0.1639 at
J_SB=2.0, 0.1891 at J_SB=4.5, converged) but never reaches 1/4. The signal is
monotonically increasing, not non-monotonic.

A "zero window" at J_SB≈3.75-4.25 has been reported from this section and is
not there: it was the same sampling trap Section 9 describes, one notch coarser
still, recording every 0.2 time units and catching oscillation zeros instead of
peaks. The threshold genuinely does not exist for equal coupling at N=4, and the
reason is quantitative rather than a dramatic extinction: the two converged
probes give 0.1639 at J_SB=2.0 and 0.1891 at J_SB=4.5, both well below ¼. What
the curve does BETWEEN those probes, and the N=5 row, are unmeasured here.

**However:** Asymmetric coupling rescues the crossing for both N=4 and N=5.
With [J_SA, J_SB1, J_SB2, ...] = [1.0, 2.0, x, x, ...], the crossing
survives as long as the remaining observers are weak enough:

| N | coupling pattern | x_crit | meaning |
|:---|:---|:---|:---|
| 4 | [1.0, 2.0, x, x] | 1.165 ± 0.005 | other B can be almost as strong |
| 5 | [1.0, 2.0, x, x, x] | 0.925 ± 0.005 | other B must be noticeably weaker |

(`star_n_observer.py` takes a per-spoke J list and so can express these
patterns, but no committed driver runs this sweep; the two x_crit values are
not re-measured here. Both are on the coarse grid, where every threshold
re-measured in this document came out too high.)

The tolerated asymmetry shrinks with N; the rescue becomes more fragile,
not less. Equal coupling kills the crossing; one dominant observer preserves it.

**Spectral diagnostic at the N=4 boundary:** The eigenvalue spectrum of
ρ_AB changes only marginally across the crossing/non-crossing line
(x=1.16 → x=1.17). The largest eigenvalue shifts from 0.7341 to 0.7322,
purity drops from 0.5866 to 0.5846. No rank collapse, no bifurcation.
The 1/4 boundary behaves like a smooth metric threshold, not a spectral
phase transition.

Peak R dilutes approximately as N^(−0.74), not 1/N.

The shadow effect (Z-measurement on A suppressing R_SB) remains visible
but is NOT the stable ~94% from Section 4.6. In the Bell-based N-observer
setup it is 8-21% and irregular with larger N.

### 8.2 Continuous measurement: ANSWERED

**Setup:** 3-qubit, Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=2.0, γ_S=γ_B=0.05.
At t_start=1.0, ramp γ_A linearly from 0.05 to γ_target over duration Δt.

**The shadow grows gradually but never matches sudden measurement.**

Snapshot suppression of R_SB at t=1.5:

| method | R_SB @ t=1.5 | suppression |
|:---|:---|:---|
| baseline | 0.0823 | - |
| sudden projective | 0.0006 | 99.2% |
| γ_A=0.5, Δt=0.1 | 0.0565 | 31% |
| γ_A=1.0, Δt=0.1 | 0.0404 | 51% |
| γ_A=5.0, Δt=0.1 | 0.0258 | 69% |
| γ_A=50, Δt=0.1 | 0.0256 | 69% |

Suppression saturates at ~69% for large γ_A and does not approach the
99% of projective measurement. Slower ramps (larger Δt) produce weaker
suppression at any given γ_target.

**Unexpected finding:** Near-instantaneous ramps (γ_A=50, Δt≤0.01)
do NOT converge to sudden measurement. They actually *increase* peak
R_SB by ~36%. Continuous σ_z dephasing and projective Z-measurement
are qualitatively different operations. Strong continuous dephasing
creates transient correlations that projective measurement destroys.

**Conclusion:** Observation is not a smooth limit. Projective measurement
(external observation) and strong decoherence (environmental noise) have
different signatures on the S-B shadow, even in the limit of infinitely
fast, infinitely strong dephasing.

### 8.3 AB with direct coupling: ANSWERED

**Setup:** 3-qubit, Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=1.466, γ=0.05.
Added J_AB ∈ {0, 0.1, 0.3, 0.5, 1.0}. (J_SB = 1.466 is the setting the sweep was
run at, kept for that reason; the threshold itself is 1.46295, Section 4.10.)

**Non-monotonic effect on threshold:**
- J_AB=0.1: slightly *worsens* threshold behavior
- J_AB=0.3-0.5: helps crossing (sweet spot at ~0.5, threshold drops
  from 1.46295 to 1.34444, both converged)
- J_AB=1.0: still crosses but much later (t≈1.0 vs t≈0.3)

**The "shadow destroyed by direct coupling" reading is withdrawn.** At
J_AB = 0.3, 0.5 and 1.0 the sweep reports a suppression of exactly 0.0%,
and an exact zero is the signature of an identity rather than of a physical
effect. It is this one: **at the instant of measurement** an unread projective
measurement on A cannot change ρ_SB at all, since summing over the projectors of
a local measurement leaves the other marginal untouched. That is a statement
about that instant only, and nothing more; afterwards the global state differs
and ρ_SB evolves differently, which is exactly the shadow Sections 4.6 and 8.2
measure. But the metric reported here is the peak of R_SB over t ≥ t_measure,
and at those couplings that peak sits AT t_measure itself, so the metric samples
precisely the one moment at which the identity forces a zero. That is measured
rather than argued: at J_AB = 0.3, 0.5 and 1.0 the argmax over t ≥ 1 is t = 1.000
in both the measured and the unmeasured run, while at J_AB = 0 it is 1.800
against 1.100 and the suppression is a real 10.26%. The number measures where the peak fell, not
what direct coupling did to the shadow. What direct coupling does to the shadow
is unmeasured here.

**Dominance crossover:** At J_AB≈0.7, direct observer coupling alone
generates AB crossing without any S-mediated coupling (J_SB=0).

### 8.4 Correlation bridge: ANSWERED

AB never crosses in symmetric 3-qubit experiments. This is consistent
with N_SCALING_BARRIER.md Section 7: crossing occurs where entanglement
lives. In the star topology, initial entanglement lives in SA (Bell state),
not in AB. The AB pair sees only S-mediated transferred entanglement,
which requires J_SB/J_SA ≥ 1.46295 at Q = J_SA/γ = 20 (Section 4.10). With symmetric coupling
(ratio = 1.0), the transfer never reaches 1/4.

This is not a separate phenomenon; it is the same locality principle.
The N-scaling barrier says global crossing fails because entanglement is
local. The star topology says AB crossing fails because the entanglement
starts in SA. Both resolve when you look at the right pair (SA crosses
at t=0.42) or create the right asymmetry (Section 8.1, 8.3).

### 8.5 Threshold formula: ANSWERED

This section and Section 4.10's γ-scan are ONE curve, read twice. Both hold
J_SA = 1 and vary γ, which by the scaling of the previous section is a scan in
Q = J_SA/γ; the ten rows here are the five there with five more between them.

**Converged data** (N=2, J_SA=1.0; t_max=6, and 14 for the two rows below
γ=0.01, which need the longer run), beside the same scan on the coarse grid,
because the shift between them has a direction and the direction is the point:

| γ | Q | J_SB threshold | coarse grid |
|:---|:---|:---|:---|
| 0.001 | 1000 | 1.17642 | 1.183 |
| 0.010 | 100 | 1.22768 | 1.247 |
| 0.020 | 50 | 1.28532 | 1.296 |
| 0.050 | 20 | 1.46295 | 1.466 |
| 0.070 | 14.3 | 1.58542 | 1.634 |
| 0.100 | 10 | 1.77496 | 1.820 |
| 0.120 | 8.3 | 1.90495 | 1.929 |
| 0.150 | 6.7 | 2.10483 | 2.146 |
| 0.170 | 5.9 | 2.24103 | 2.253 |
| 0.200 | 5 | 2.44923 | 2.460 |

Every converged value lies BELOW its coarse-grid counterpart, which is not luck:
an undersampled peak is always too low and the threshold read from it always too
high, so the whole column had to move one way.

**Fit on the converged data:** J_th(γ) ≈ 7.232 · γ^1.081 + 1.1766, with a
largest residual of 0.0043, at the γ=0.001 end. It is a fit on 0.001 ≤ γ ≤ 0.2
and nothing more; in particular its constant term is NOT an asymptote. Section
4.10 measures the threshold below this range and finds it still falling, to
1.11567 at γ=0.0001, so any reading of 1.1766 as a zero-noise limit is a
property of where the table stops. An earlier fit, 7.35 · γ^1.08 + 1.18, is
still quoted in places: its exponent is right but it does not come from the
table it was printed beside, where its largest residual is 0.038 against the
0.0043 of a free refit.

**The curvature is real, and a linear approximation is not a substitute.** A
linear form has been offered as working "as well" on the strength of R², and on
coarse-grid data it does. On the converged data it does not. The best straight
line here is 6.3698 · γ + 1.1537, and its largest residual is 0.022 against the
power fit's 0.0043, a factor of five, while R² moves only from 0.999965 to
0.999110. R² is the wrong instrument for a ten-point curve with this little
scatter.

**No divergence or hard closure at γ=0.2.** The threshold exists at 2.449 and
the window merely gets narrower.

## 9. Numerical Notes

- **The re-measurement.** Every converged number in this document is printed by
  [`simulations/star_topology_converged.py`](../simulations/star_topology_converged.py)
  (output in [`results/star_topology_converged_run.txt`](../simulations/results/star_topology_converged_run.txt)).
  It uses an exact propagator, expm of the vectorised Lindblad generator applied
  step by step, with the peak refined by a parabola through the sampled maximum
  and window edges interpolated rather than counted; the repo's RK4 at dt=0.001
  recording every step is run beside it and lands on the same threshold to 4e-6
  in J_SB. Thresholds are bisected to 2e-6; halving the step from 0.008 down to
  0.001 moves the N=2 threshold by 8e-6 in total and by under 1e-6 over the last
  two halvings, which is where the fifth decimal comes from. At N=3 the same
  ladder spans 4e-4 and the value is quoted to four decimals.
- **What was NOT re-measured**, and is March's throughout: Sections 4.6 and 4.7,
  Section 8.2 entirely, Section 8.1's two x_crit values, its spectral diagnostic at the N=4 boundary,
  its N^(−0.74) dilution and its 8-21% shadow, and Section 4.12's peak-interval
  statistics. Every one of them is on the coarse grid, where every quantity that
  has been re-measured moved, and all in one direction.
- **The sampling trap to know about when using this code.**
  `run_star_topology` defaults to `sample_every=20`: it integrates at dt=0.005
  but records the observable only every 0.1, against an oscillation period near
  0.66. Six points per period systematically miss peaks, and since the peak is a
  supremum the miss has a sign: the peak comes out too low and any threshold
  read from it too high.
- Integration: RK4, dt=0.005, t_max=5.0 (extended to 40.0 for frequency analysis)
- Purity bounded ≤ 1.0 for all runs (Euler v1 had artifacts > 1.0)
- Euler v1 showed spurious oscillating crossings for |0++⟩ at t≈2.9
  which are absent in RK4. These were integration artifacts.
- Partial traces validated: Tr(ρ_pair) = 1 and hermiticity confirmed
- Concurrence computed via standard Wootters formula
- FFT: Hanning window, DC removed, rfft

## 10. Simulation Code

- [`simulations/star_topology_v2.py`](../simulations/star_topology_v2.py): 3-qubit star topology, RK4 integration
- [`simulations/star_n_observer.py`](../simulations/star_n_observer.py): N-qubit with asymmetric coupling
- [`simulations/star_topology_v3.py`](../simulations/star_topology_v3.py): N-qubit with equal coupling, J_AB support, threshold sweeps
- [`simulations/star_topology_converged.py`](../simulations/star_topology_converged.py): the converged re-measurement, exact propagator beside the repo's RK4
- [`simulations/verify_star_topology.py`](../simulations/verify_star_topology.py): eight claim checks against `star_topology_v2.py`

---

*See also: [Orphaned Results](ORPHANED_RESULTS.md), echo effect fully characterized: Bohr frequencies, 8γ/3 envelope decay, N=4,5 scaling*
