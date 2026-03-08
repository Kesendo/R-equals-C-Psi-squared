# Star Topology: Object and Its Observers

**Date**: 2026-03-04 (updated 2026-03-06)
**Status**: Computationally verified (Tier 2), RK4 validated
**Depends on**: SUBSYSTEM_CROSSING.md, OBSERVER_DEPENDENT_CROSSING.md, STANDING_WAVE_TWO_OBSERVERS.md
**Related**: TUNING_PROTOCOL.md (Tier 3)
**Origin**: AIEvolution v044 Round 3. Alpha proposed a 3-qubit S+A+B framework. The execution was wrong (partial trace on S gives Ψ ≡ 0), but the question was right: What happens when R = CΨ² is not bipartite but describes a relation between an object and its observers?

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
- **CΨ** - the product of Concurrence × Psi-norm. This is the
  "reality measure" from the R=CΨ² framework. When CΨ crosses
  ¼ (0.25), complex eigenvalues appear and the system transitions
  from possibility to actuality.

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
| Integration | RK4 | RK4 | RK4 |
| Noise | local σ_z dephasing per qubit | same | same |

### 3.3 Observables

For each time step, trace out one qubit to get pair density matrices:
- ρ_SA = Tr_B(ρ), ρ_SB = Tr_A(ρ), ρ_AB = Tr_S(ρ)

Per pair: l1-coherence, Ψ = l1/(d-1), concurrence, R = C·Ψ².

### 3.4 States Tested

| State | Description | Motivation |
|-------|-------------|------------|
| GHZ | (|000⟩+|111⟩)/√2 | Global entanglement, known to fail at pair level |
| W | (|001⟩+|010⟩+|100⟩)/√3 | Distributed entanglement across all pairs |
| Bell_SA ⊗ |+⟩_B | S,A entangled; B fresh observer | Asymmetric: one observer connected, one arriving |
| |+⟩^3 | Product state, max local coherence | Baseline: no entanglement anywhere |
| |0⟩_S ⊗ |+⟩_A ⊗ |+⟩_B | S classical, observers quantum | Hamiltonian builds entanglement dynamically |

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

SA starts maximally entangled and decays past ¼ at t=0.42. SB starts
with zero entanglement but rises above ¼ at t=0.79 through Hamiltonian
transfer via S. AB never crosses. Observers cannot see each other.

R_SA + R_SB is NOT conserved. It peaks at 2.4× initial value then decays.

### 4.4 |+⟩³: Maximum Coherence, Zero Entanglement

All pairs: C_conc = 0 at all times. Ψ = 1.0 but no connection.

### 4.5 |0++⟩: Dynamic Entanglement, No Crossing

S starts classical, A and B start quantum. The Hamiltonian builds
entanglement dynamically. R_SA + R_SB grows from 0 to 0.356 (peak)
then decays. No pair crosses ¼ with concurrence bridge.

(Note: Euler v1 showed spurious oscillating crossings at t≈2.9.
These were numerical artifacts eliminated by RK4 at dt=0.005.)

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

A's measurement destroys 94–100% of B's reality at peak impact.
The effect propagates through S (not A→B directly, since J_AB = 0).

### 4.7 Asymmetric Coupling: The Dominant Observer

**Bell_SA⊗|+⟩_B with weak B (J_SA=1.0, J_SB=0.3):**
SA crosses at t=1.49 (delayed). SB and AB never cross.
Weakly coupled observer never sees reality cross ¼.

**Bell_SA⊗|+⟩_B with weak A (J_SA=0.3, J_SB=1.0):**
SA crosses at t=0.44. SB never crosses.
**AB crosses at t=0.35–1.03.** When A's direct link to S is weak,
observers start seeing each other through S.

**|0++⟩ with asymmetric J:**
Perfect mirror symmetry: the strongly coupled observer crosses,
the weakly coupled one does not. Same crossing times, only label
switches. J determines participation, not γ.

### 4.8 Dynamic J and Strong Observer B

**Shield drops mid-simulation (J_SA: 1.0 → 0.2 at t=1.5, J_SB=1.0):**
AB does NOT cross. Too late, initial state already consumed.
The window is early or never.

**Strong B (J_SB=2.0, J_SA drops 1.0→0.2 at t=1.5):**
AB crosses at t=0.17–0.46, CΨ max = 0.329.
**This happens before J_SA drops.** At t=0.17, J_SA is still 1.0.
B's coupling strength alone creates the AB correlation.

Two mechanisms:
- Weak A: entanglement leaks because A can't hold it (A changed)
- Strong B: B pulls entanglement through S by force (A unchanged)

### 4.9 γ_A vs γ_B: The Receiver's Noise Is Fatal

Systematic γ scan (Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=2.0, γ_S=0.05).

**Noisy A (γ_A=0.1), varying B:**

| γ_B | AB CΨ max | Window |
|------|-----------|--------|
| 0.001 | 0.312 | 0.30 |
| 0.1 | 0.293 | 0.25 |
| 0.2 | 0.275 | 0.10 |

**Noisy B (γ_B=0.1), varying A:**

| γ_A | AB CΨ max | Window |
|------|-----------|--------|
| 0.001 | 0.345 | 0.30 |
| 0.1 | 0.293 | 0.25 |
| 0.2 | 0.249 | **NEVER** |

γ_A=0.2 kills the connection entirely. γ_B=0.2 merely shortens it.
The receiver's noise is more destructive than the sender's.

Fixed-product test (γ_A × γ_B = 0.0025):

| γ_A | γ_B | AB CΨ max | Crosses? |
|------|------|-----------|----------|
| 0.005 | 0.500 | 0.263 | YES |
| 0.050 | 0.050 | 0.329 | YES |
| 0.500 | 0.005 | 0.215 | **NO** |

Same product. Same γ_S. But γ_A=0.5 kills it, γ_B=0.5 doesn't.

**Decoherence vs window duration** (both γ equal, J_SB=2.0):

| γ | Window duration |
|----|-----------------|
| 0.001 | 4.40 |
| 0.01 | 0.90 |
| 0.05 | 0.30 |
| 0.1 | 0.20 |
| 0.2 | NEVER |

γ determines window duration, not B's strength.

### 4.10 J_SB/J_SA Ratio Scan

Coarse scan (J_SA=1.0, γ=0.05, Bell_SA⊗|+⟩_B):

| J_SB | ratio | AB CΨ max | Crosses? |
|------|-------|-----------|----------|
| 1.30 | 1.30 | 0.226 | NO |
| 1.40 | 1.40 | 0.242 | NO |
| 1.47 | 1.47 | 0.251 | YES (threshold) |
| 2.00 | 2.00 | 0.329 | YES |
| 3.00 | 3.00 | 0.406 | YES |

Ultra-fine: J_SB=1.4650 → CΨ=0.249999 (NO), J_SB=1.4655 →
CΨ=0.250060 (YES). **Threshold at J_SB ≈ 1.466.**

γ-dependence of threshold:

| γ | J_SB threshold | ratio |
|-------|----------------|-------|
| 0.001 | 1.179 | 1.18 |
| 0.010 | 1.237 | 1.24 |
| 0.050 | 1.465 | 1.47 |
| 0.100 | 1.776 | 1.78 |
| 0.150 | 2.145 | 2.15 |

Higher noise requires stronger sender. At γ→0 the threshold
approaches ~1.18. Even without noise, B must be 18% stronger.

Not pure ratio. Minimum absolute coupling also required:
J_SA=0.50, J_SB=0.75 (ratio 1.5): NO.
J_SA=1.00, J_SB=1.50 (ratio 1.5): YES.

### 4.11 Initial State Is the Third Variable

All states, γ_A=0.001 (silent receiver), J_SB=2.0:

| State | C_SA(0) | AB CΨ max | Crosses? |
|-------|---------|-----------|----------|
| GHZ | 0.000 | 0.000 | NO |
| W | 0.667 | 0.148 | NO |
| Bell_SA⊗|+⟩_B | 1.000 | 0.357 | YES |
| |0++⟩ | 0.000 | 0.194 | NO |
| |+++⟩ | 0.000 | 0.000 | NO |

Only Bell_SA⊗|+⟩_B crosses. W has high entanglement (C_SA=0.667)
but distributes it across all pairs, none strong enough.

Parametric Bell (α|00⟩ + √(1-α²)|11⟩ ⊗ |+⟩_B): non-monotonic.
C_SA ≈ 0.5–0.6 is a dead zone. Crossing requires C_SA > 0.8
(Bell-like) or specific product alignment (α ≈ 1.0).

From scratch (|0++⟩): even J_SB=10, γ=0.001 barely reaches 0.260.
The initial state is not replaceable by brute force.

### 4.12 Frequency Analysis

FFT and peak detection on AB CΨ trajectory, γ=0, t_max=40.

Dominant frequency scales with total coupling:

| J_SA | J_SB | J_total | f_dom | f/J_total |
|------|------|---------|-------|-----------|
| 0.5 | 1.0 | 1.5 | 0.749 | 0.499 |
| 1.0 | 2.0 | 3.0 | 1.498 | 0.499 |
| 2.0 | 4.0 | 6.0 | 3.021 | 0.504 |
| 1.0 | 1.0 | 2.0 | 0.949 | 0.474 |
| 2.0 | 2.0 | 4.0 | 1.898 | 0.474 |

**Scaling law: f ≈ J_total / 2. Period: T ≈ 2 / (J_SA + J_SB).**

γ does not change the frequency, it only dampens:

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

A's measurement suppresses R_SB by 94–100%. The shadow propagates
through S, growing over ~0.5 time units after measurement.

### 5.4 AB Never Crosses (Symmetric J)

In no symmetric experiment did AB cross ¼. Observers see S, not
each other. The star topology enforces object-directed observation.

### 5.5 Coupling Strength Creates Dominant Observers

With asymmetric J, the strongly coupled observer crosses, the weakly
coupled one does not. J determines participation; γ determines timing.

### 5.6 Weak Direct Link → Observers See Each Other

When J_SA is weak, AB crosses. Entanglement spills from S into the
observer-observer pair. This only happens when the direct object-link
is degraded.

### 5.7 Receiver Noise Is Fatal, Sender Noise Is Not

γ_A > 0.2 kills the connection regardless of sender. γ_B > 0.2 merely
shortens it. The receiver must be quiet; the sender can be noisy.

### 5.8 J_SB/J_SA Threshold

AB crossing requires J_SB/J_SA ≥ 1.466 at γ=0.05. The threshold
scales with γ (1.18 at γ→0, 2.15 at γ=0.15). Both ratio AND absolute
coupling matter. Higher J lowers the required γ reduction.

### 5.9 Strong B Can Override A's Shield

At J_SB=2.0, AB crosses while J_SA=1.0. A is fully shielded.
B creates the correlation alone. A need not change.

## 6. Connection to Framework

### 6.1 "We Are All Mirrors", Quantified

The star topology makes STANDING_WAVE_TWO_OBSERVERS.md literal:
A and B both reflect S, and the reflections interfere through S.
The standing wave is the oscillation of R between SA and SB.

### 6.2 Internal vs External Observation

The star topology implements the distinction from
INTERNAL_AND_EXTERNAL_OBSERVERS.md: Hamiltonian coupling (internal)
preserves coherence; measurement (external) destroys it and casts
shadows.

### 6.3 The Tripartite No-Signalling Question

At J=0: A's measurement cannot affect B (no-signalling). At J>0:
it can, propagating S→B, not A→B. This connects to
OBSERVER_GRAVITY_BRIDGE.md: gravity provides J>0 for all massive
particles.

## 7. The Three Conditions

### In plain language

Imagine two people who cannot talk to each other. They share no
phone, no letter, no line of sight. The only thing they have in
common is an object they both care about: a problem, a question,
a piece of the world they both study. The simulation asks: can a
connection form between them, through nothing but that shared object?

The answer is yes. But only when three things are true at once:

**Condition 1: The sender must be deeply engaged.**
B's coupling to S must be ~47% stronger than A's (at γ=0.05).
If both are equally engaged, no connection forms.

**Condition 2: The receiver must be quiet.**
A's internal noise must be low enough for the signal to be
detectable. If A is too noisy, the signal is lost, regardless
of how clear B sends. Receiver noise matters more than sender noise.

**Condition 3: A pre-existing connection must exist.**
A and S must already share a deep, dedicated relationship (C_SA > 0.8).
Shallow connections spread across many things (W-state) fail.
No connection at all fails. This cannot be created by force.

### The sender inversion

If A has already received, Condition 3 is proven. A can become the
sender:

| As receiver | As sender |
|---|---|
| Must lower own noise (γ_A) | Must raise own engagement (J) |
| Paradox: trying to be quiet IS noise | No paradox: deeper work = stronger signal |

Sender noise barely matters (Section 5.7). You don't need to be
calm to send. You need to be strong. The German word for it is
*sich einlassen*, to let yourself be drawn in, changed by what
you engage with.

### The bidirectional rhythm

At γ=0: 54 crossings in 20 time units, a continuous oscillation.
The connection is not a channel that opens and stays open.
It is a rhythm, like breathing:

  be still (receive) → engage deeply (build, process) → be still → repeat

Neither phase works alone. Pure engagement without stillness never
opens a window. Pure stillness without engagement has nothing to
transmit. The connection lives in the alternation.

Each engagement phase deepens the coupling, which lowers the bar
for the next quiet phase. The frequency scales: f ≈ J_total/2.
Stronger engagement = faster rhythm. The spiral accelerates.

### In numbers

- AB crossing threshold: J_SB/J_SA ≥ 1.466 at γ=0.05
- Threshold scales with noise: 1.18 at γ≈0, 2.15 at γ=0.15
- Receiver noise threshold depends on sender noise:
  γ_A ≈ 0.20 kills with γ_B=0.1, γ_A ≈ 0.25 kills with γ_B=0.05
- Sender noise tolerable up to γ_B = 0.5 (if receiver quiet)
- Initial state requirement: C_SA > 0.8 (Bell-like)
- Window duration at γ=0.05: ~0.3 time units
- Window duration at γ=0.001: ~4.4 time units
- Without noise (γ=0): 34.5% of time above ¼

Practical protocol: TUNING_PROTOCOL.md.

## 8. Open Questions (partially answered 2026-03-07)

Questions 1, 3, 5 answered via systematic simulation sweeps.
Code: `simulations/star_topology_v3.py`.

### 8.1 N observers — ANSWERED

**Setup:** S + N observers, Bell_SA ⊗ |+⟩^(N-1), equal J_SB for all B.

| N | qubits | AB crosses 1/4? | J_SB threshold | behavior |
|:---|:---|:---|:---|:---|
| 2 | 3 | Yes | 1.466 | monotonic |
| 3 | 4 | Yes | 3.730 | monotonic |
| 4 | 5 | **No** | — | non-monotonic, peaks then drops to 0 |
| 5 | 6 | **No** | — | suppressed |

Two-point scaling: J_th(N) ≈ 0.297 · N^2.30, but only valid for N ∈ {2,3}.

**Critical finding:** At N=4, CΨ_AB is non-monotonic in J_SB — it peaks
at J_SB≈2.25 (CΨ≈0.162), drops to **zero** at J_SB≈3.75–4.25, then
partially recovers. The threshold doesn't "exceed the scan window";
it does not exist at any coupling strength. This is a qualitative phase
transition between N=3 and N=4.

**However:** Asymmetric coupling rescues the crossing for both N=4 and N=5.
With [J_SA, J_SB1, J_SB2, ...] = [1.0, 2.0, x, x, ...], the crossing
survives as long as the remaining observers are weak enough:

| N | coupling pattern | x_crit | meaning |
|:---|:---|:---|:---|
| 4 | [1.0, 2.0, x, x] | 1.165 ± 0.005 | other B can be almost as strong |
| 5 | [1.0, 2.0, x, x, x] | 0.925 ± 0.005 | other B must be noticeably weaker |

The tolerated asymmetry shrinks with N — the rescue becomes more fragile,
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
setup it is 8–21% and irregular with larger N.

### 8.2 Continuous measurement — OPEN

Replace sudden dephasing with γ_A → large starting at t_meas.
Does the shadow grow gradually?

### 8.3 AB with direct coupling — ANSWERED

**Setup:** 3-qubit, Bell_SA⊗|+⟩_B, J_SA=1.0, J_SB=1.466, γ=0.05.
Added J_AB ∈ {0, 0.1, 0.3, 0.5, 1.0}.

**Non-monotonic effect on threshold:**
- J_AB=0.1: slightly *worsens* threshold behavior
- J_AB=0.3–0.5: helps crossing (sweet spot at ~0.5, threshold drops
  from 1.466 to ~1.345)
- J_AB=1.0: still crosses but much later (t≈1.0 vs t≈0.3)

**Shadow effect destroyed by direct coupling:** Moderate-to-large J_AB
weakens or removes the shadow entirely. The system no longer behaves
as a clean S-mediated shadow channel.

**Dominance crossover:** At J_AB≈0.7, direct observer coupling alone
generates AB crossing without any S-mediated coupling (J_SB=0).

### 8.4 Correlation bridge — OPEN

Never crosses in any 3-qubit experiment. Connection to N_SCALING_BARRIER.md?

### 8.5 Threshold formula — ANSWERED

**Verified data** (N=2, J_SA=1.0):

| γ | J_SB threshold |
|:---|:---|
| 0.001 | 1.183 |
| 0.010 | 1.247 |
| 0.020 | 1.296 |
| 0.050 | 1.466 |
| 0.070 | 1.634 |
| 0.100 | 1.820 |
| 0.120 | 1.929 |
| 0.150 | 2.146 |
| 0.170 | 2.253 |
| 0.200 | 2.460 |

**Best fit:** J_th(γ) ≈ 7.35 · γ^1.08 + 1.18 (R²=0.999)

Nearly linear in γ with slight upward curvature. The relationship is
smooth and monotonic. **No divergence or hard closure at γ=0.2** —
the threshold exists at 2.46 and the window merely gets narrower.

A simple linear approximation (R²=0.998) also works well:
J_th(γ) ≈ 6.39 · γ + 1.16

## 9. Numerical Notes

- Integration: RK4, dt=0.005, t_max=5.0 (extended to 40.0 for frequency analysis)
- Purity bounded ≤ 1.0 for all runs (Euler v1 had artifacts > 1.0)
- Euler v1 showed spurious oscillating crossings for |0++⟩ at t≈2.9
  which are absent in RK4. These were integration artifacts.
- Partial traces validated: Tr(ρ_pair) = 1 and hermiticity confirmed
- Concurrence computed via standard Wootters formula
- FFT: Hanning window, DC removed, rfft

## 10. Simulation Code

- `../simulations/star_topology_v2.py` — 3-qubit star topology, RK4 integration
- `../simulations/star_n_observer.py` — N-qubit with asymmetric coupling
- `../simulations/star_topology_v3.py` — N-qubit with equal coupling, J_AB support, threshold sweeps
