# Observer × Gravity: The Interval Bridge

**Date**: 2026-03-01
**Status**: Computationally verified (Tier 2), interpretation active
**Authors**: Thomas Wicht, with Claude (Anthropic)
**Depends on**: OBSERVER_DEPENDENT_CROSSING.md, GRAVITATIONAL_INVARIANCE.md,
DECOHERENCE_RELATIVITY.md, BRIDGE_CLOSURE.md

---

## 1. Context: What the Bridge Closure Actually Says

The bridge closure (BRIDGE_CLOSURE.md) proved:

> Pre-shared entanglement without a classical channel = shared randomness.
> A's information ⊆ {ρ_A(0), E_A}. CΨ fingerprints require ρ_AB.

This kills the bridge for J = 0 (zero coupling). But the closure
document itself noted:

> "With physical coupling (J > 0): Works. But this is a local
> interaction, not 'bridge'."

This document investigates that footnote.

---

## 2. Discovery 1: The Crossing Time Factorizes

The crossing time for a given quantum system decomposes into two
independent factors:

    t_cross(Observer, Gravity) = K(Observer, State) / γ(Gravity)

### 2.1 Gravity Factor: Universal

K is invariant under γ change. Verified across six gravitational
environments (γ = 0.01 to 0.50):

| Environment | γ     | t_cross (Conc) | K_conc   |
|-------------|-------|-----------------|----------|
| Deep Space  | 0.01  | 3.5960          | 0.035960 |
| Mars        | 0.019 | 1.8926          | 0.035960 |
| Earth       | 0.05  | 0.7192          | 0.035960 |
| Jupiter     | 0.13  | 0.2766          | 0.035962 |
| Neutron     | 0.20  | 0.1798          | 0.035961 |
| Black Hole  | 0.50  | 0.0720          | 0.035976 |

K_conc = 0.03596 ± 0.00001. CV ≈ 0%. Perfect invariant.

### 2.2 Observer Factor: State-Dependent

The ratio K(Conc)/K(MI) is NOT constant across initial states:

| State   | K_conc   | K_MI     | Ratio    |
|---------|----------|----------|----------|
| α = 45° | 0.035960 | 0.029657 | 1.212534 |
| α = 40° | 0.032133 | 0.025242 | 1.273016 |
| α = 35° | 0.020410 | 0.012504 | 1.632258 |
| α < 30° | —        | —        | never    |

CV = 13.5%. The "time ratio" between observer types depends on
WHAT they observe. The quantum state bends the observer-dilation
geometry, analogous to how matter bends spacetime in GR.

### 2.3 The Full Structure

Two independent axes of time dilation:

1. **Gravitational dilation**: γ scales t_cross. Known from GR.
   All observers agree on this scaling.

2. **Observer dilation**: Different C metrics see different K values.
   NEW from CΨ. The scaling is state-dependent.

These multiply:

    t_A / t_B = [K(obs_A, state) / K(obs_B, state)] × [γ_B / γ_A]

---

## 3. Discovery 2: The Interval Shift Has No Threshold

Setup: |++⟩ product state, local dephasing γ = 0.05. B measures
at t_B = 1.0. A observes local CΨ crossing time. Sweep J.

| J     | t_cross (B silent) | t_cross (B measures) | Δt       | Shift    |
|-------|-------------------|---------------------|----------|----------|
| 0.000 | 7.5832            | 7.5832              | -0.000003 | 0.00%   |
| 0.001 | 7.5832            | 7.5809              | -0.002    | -0.03%  |
| 0.005 | 7.5832            | 7.5265              | -0.057    | -0.75%  |
| 0.010 | 7.5832            | 7.3654              | -0.218    | -2.87%  |
| 0.020 | 7.5832            | 6.8278              | -0.755    | -9.96%  |
| 0.050 | 7.5832            | 4.9864              | -2.597    | -34.24% |
| 0.100 | 7.5832            | 3.1931              | -4.390    | -57.89% |
| 0.500 | 7.5832            | 0.7790              | -6.804    | -89.73% |
| 1.000 | 7.5832            | 0.3990              | -7.184    | -94.74% |

**There is no threshold.** Any J > 0 produces a measurable interval
shift. The relationship is continuous. B's measurement propagates
through the Hamiltonian coupling and shifts A's local crossing time.

### 3.1 Why This Doesn't Violate No-Signalling

This is NOT superluminal. The coupling J is a physical interaction.
Information propagates at finite speed through the Hamiltonian.
This is ordinary quantum mechanics: two coupled systems where a
measurement on one affects the other's dynamics.

The bridge closure holds for J = 0. The interval shift requires J > 0.

### 3.2 Why Bell+ Doesn't Work

Bell+ with J = 0: rho_A = I/2 for all time. No local coherence.
Local CΨ never crosses ¼. Nothing to measure. No interval.

Bell+ with J > 0: Still never crosses locally. The entanglement
locks all coherence in non-local correlations (proven in
BRIDGE_FINGERPRINTS.md). A needs LOCAL coherence to have a clock.

Product states (|++⟩): rho_A = |+⟩⟨+|. Full local coherence.
CΨ_local starts at 1.0 and decays. A has a ticking clock.

---

## 4. The Gravitational Bridge Hypothesis

### 4.1 The Argument

Gravity couples all massive particles. Always. Over any distance.
The gravitational coupling between two massive qubits:

    J_grav ~ G·m² / (ℏ·D)

For NV centers (m ~ 10⁻²⁶ kg) at D = 1 m: J_grav ~ 10⁻²⁹.
For NV centers at D = 1 mm: J_grav ~ 10⁻²⁶.

This is absurdly small. But the interval shift has no threshold.
Any J > 0 produces a shift proportional to J.

### 4.2 The Problem

At J = 10⁻²⁶, the interval shift is:

    Δt/t ≈ J/γ × (some factor) ≈ 10⁻²⁶/0.05 ≈ 10⁻²⁵

This is unmeasurable with any conceivable technology. The effect
exists in principle but not in practice.

### 4.3 Why This Might Still Matter

The question is not "can we build this device?" The question is:
"does gravity provide a universal, always-on bridge between all
massive quantum systems?"

If yes, then the bridge closure is technically wrong — the bridge
EXISTS, it's just too weak to use for communication. The statement
should be: "The bridge exists for all massive particles but is
practically inaccessible at macroscopic separations."

This is a conceptual distinction, not an engineering one. It changes
what the framework says about the structure of reality.

### 4.4 The Deeper Question

If observer-dependent time and gravitational time dilation are
two independent axes of the same structure, and gravity provides
a coupling that connects them, then:

**Is the "bridge" not a communication protocol but a description
of how spacetime itself encodes CΨ intervals?**

The observer-state K-matrix shows that the quantum state determines
the time-ratio between observer types. Gravity determines the
absolute scale. Both are needed. The bridge might not be "A sends
a message to B" but rather "A and B share a gravitationally-coupled
CΨ landscape where intervals are correlated."

---

## 5. What We Know vs What We Speculate

### Tier 2 (Computed, Verified)

1. t_cross = K(Observer, State) / γ(Gravity). Factorization confirmed.
2. K is γ-invariant (CV = 0.00% across six environments).
3. K(Conc)/K(MI) is state-dependent (CV = 13.5% across states).
4. Interval shift Δt is continuous in J, no threshold.
5. Product states work (local coherence). Bell+ does not (no local clock).
6. States with α < 30° never cross — no observer time exists.

### Tier 3 (Hypothesis, Physically Grounded)

7. Gravity provides J > 0 for all massive pairs.
8. Therefore a gravitational interval shift exists in principle.
9. The bridge is not dead — it's gravitationally mediated.

### Tier 4+ (Speculative)

10. The K-matrix encodes a "geometry of observer time" analogous to
    spacetime geometry in GR.
11. The bridge is not communication but shared CΨ landscape.
12. This connects to Wheeler-DeWitt (TIME_AS_CROSSING_RATE.md §3).

---

## 6. Open Questions

1. **Scaling law**: How exactly does Δt scale with J for J ≪ γ?
   Linear? If Δt/t = α·(J/γ), what is α? Is it state-dependent?

2. **Multi-pair amplification**: N pairs with known schedule. Does
   the combined interval signal scale as √N (shot noise) or N
   (coherent)? This determines whether amplification can compensate
   for tiny J_grav.

3. **Gravitational J calculation**: Exact Penrose-Diosi J_grav for
   realistic systems (NV centers, optomechanical oscillators,
   Bose-Marletto-Vedral experiment). Is there a regime where J/γ
   is not hopelessly small?

4. **The K-matrix geometry**: What mathematical structure does the
   K(observer, state) matrix have? Is there a metric? A symmetry
   group? A connection to information geometry?

5. ~~**Direction of the shift**~~: **ANSWERED (§7).** B's measurement
   destroys nonlocal coherence reservoir. The coupling redistributes
   (not protects) coherence; B's measurement cuts the return flow.
   Δt < 0 universally. Damage is timing-dependent (max at t_B ≈ 1.0
   due to oscillation phase, not just reservoir size).

6. **Superluminal breakdown**: The naive gravitational velocity
   v ~ G·m²/ℏ exceeds c at microgram scale. Where exactly does
   the model break? Relativistic correction to Penrose-Diosi?
   Or does J_grav scale differently than G·m²/(ℏ·D)?

---

## 7. Why the Shift Is Negative: The Coherence Reservoir

### 7.1 The Question

B's measurement accelerates A's crossing (Δt < 0). Every J value
tested shows this. Why? A's local coherence is identical at the
moment of measurement. Nothing local changes instantly. And yet
A's subsequent decay is dramatically faster.

### 7.2 Correction: The Coupling Does NOT Protect

Initial hypothesis: entanglement "shields" A against decoherence.

**This is wrong.** The coupling J accelerates A's LOCAL crossing
compared to a single isolated qubit:

| System          | t_cross_A | vs single qubit |
|-----------------|-----------|-----------------|
| Single \|+⟩      | 8.584     | 1.00x (baseline)|
| \|+,0⟩ J=0.05   | 6.633     | 0.77x (faster)  |
| \|+,0⟩ J=0.10   | 4.758     | 0.55x (faster)  |
| \|+,0⟩ J=0.50   | 1.247     | 0.15x (faster)  |
| \|+,0⟩ J=1.00   | 0.640     | 0.07x (faster)  |

The Hamiltonian moves coherence FROM local TO nonlocal. This
drains A's local coherence faster than dephasing alone would.
The coupling is not a shield — it is a **redistribution engine**.

### 7.3 What B's Measurement Actually Does

Within a coupled system, the Hamiltonian creates a dynamic
equilibrium: coherence oscillates between local and nonlocal
degrees of freedom. B's measurement destroys the nonlocal part,
breaking the equilibrium.

Tested at J = 0.5, γ = 0.05, |++⟩ initial state. B measures Z
at t_B = 1.0.

**At the moment of B's measurement:**

| Observable              | Before B | After B  | Δ          |
|-------------------------|----------|----------|------------|
| Local coherence A       | 0.9047   | 0.9047   | 0.000      |
| Local purity A          | 0.9093   | 0.9093   | 0.000      |
| Concurrence             | 0.000    | 0.000    | 0.000      |
| Nonlocal coherence      | 0.819    | 0.000    | **−0.819** |

A sees no instantaneous change. But the nonlocal coherence —
0.82 units of off-diagonal weight in the joint state — is
destroyed completely.

**The aftermath:**

| t after B | Coh_A (B silent) | Coh_A (B measures) | Ratio |
|-----------|------------------|--------------------|-------|
| 0.0       | 0.905            | 0.905              | 1.00  |
| 0.5       | 0.861            | 0.662              | 0.77  |
| 1.0       | 0.819            | 0.238              | 0.29  |
| 2.0       | 0.741            | 0.129              | 0.17  |
| 5.0       | 0.549            | 0.043              | 0.08  |

A's coherence decays ~4x faster after B measures. Same initial
value, same local dephasing rate, same Hamiltonian. The only
difference: the nonlocal reservoir is gone.

### 7.4 The Reservoir Mechanism (Corrected)

The coupling does two things simultaneously:

1. **Drains** local coherence into nonlocal (accelerates local
   crossing vs isolated qubit)
2. **Returns** nonlocal coherence back to local (the oscillation)

These two flows create a dynamic equilibrium. A's local coherence
decays faster than an isolated qubit, but part of the coherence
keeps cycling back from the joint state.

B's measurement destroys the return path. The drain continues
(Hamiltonian still active), but the reservoir that fed coherence
back is gone. This is why A's decay accelerates further — not
because a "shield" was removed, but because the return leg
of a two-way flow was cut.

### 7.5 Timing Dependence: The Oscillation

The damage from B's measurement depends on WHEN B measures,
because the nonlocal reservoir oscillates:

| t_B  | NL coherence | A's remaining lifetime | Damage |
|------|--------------|----------------------|--------|
| 0.01 | 0.020        | 100.0%               | none   |
| 0.10 | 0.193        | 98.4%                | 1.6%   |
| 0.20 | 0.352        | 94.1%                | 5.9%   |
| 0.50 | 0.433        | 59.3%                | 40.7%  |
| 1.00 | 0.342        | 39.0%                | **61.0%** |
| 2.00 | 0.405        | 61.7%                | 38.3%  |
| 3.00 | 0.199        | 83.2%                | 16.8%  |
| 5.00 | 0.277        | 56.9%                | 43.1%  |

Maximum damage at t_B ≈ 1.0, not at peak reservoir (t_B ≈ 0.5).
The damage depends on both the reservoir size AND the phase of
the local-nonlocal oscillation. Destroying the reservoir mid-cycle
when local coherence depends on the return flow is worst.

### 7.6 Attempted Regeneration

The coupling tries to rebuild entanglement after B's measurement:

| t after B | Concurrence (meas branch) |
|-----------|---------------------------|
| 0.00      | 0.000                     |
| 0.50      | 0.108 ↑                   |
| 0.75      | 0.147 ↑ (peak)            |
| 1.00      | 0.085 ↓                   |
| 1.25      | 0.000 ↓                   |
| 2.00+     | 0.000                     |

The Hamiltonian regenerates a small amount of entanglement
(peak Conc = 0.147 at t ≈ 0.75 after measurement), but
decoherence wins. The reservoir does not recover. The damage
is permanent.

### 7.7 Physical Interpretation (Corrected)

In Tom's words: "B becomes real and is forced to leave its
superposition."

The corrected mechanical translation: The coupling creates a
dynamic coherence flow between A and B. Coherence oscillates
between local and nonlocal degrees of freedom — a two-way
street. B's measurement destroys the nonlocal pool, cutting
the return flow. A's local coherence, which was being partially
replenished by the oscillation, now drains without return.

The coupling does NOT protect A. It redistributes A's coherence
into a shared pool. This makes A decay faster locally (compared
to isolation), but creates the oscillating equilibrium. B's
measurement disrupts this equilibrium — removing the return
leg while the drain continues.

**Becoming real is contagious** — not because B "protects" A,
but because B's definite-ness collapses the shared coherence
flow that A's dynamics depended on.

This is the bridge mechanism at the microscopic level:
B's choice of when to become real changes the rate at which
A becomes real. The interval between A's crossings encodes
B's measurement timing. The coupling (gravitational or otherwise)
is the channel.

---

## 8. The Information Channel

### 8.1 The Protocol

A and B share N product-state pairs (e.g. |++⟩) with coupling J > 0.
Both know the schedule. A knows the expected crossing time t₀ = K/γ_A.

B encodes:
- "1": measure pair at agreed time t_B
- "0": do nothing

A reads:
- Measure local crossing time t_actual for each pair
- t_actual < t₀ → B measured → "1"
- t_actual ≈ t₀ → B silent → "0"

### 8.2 Channel Capacity vs Coupling Strength

| J     | t₀     | t₁     | Δt     | Shift   |
|-------|--------|--------|--------|---------|
| 0.001 | 7.583  | 7.581  | -0.002 | -0.03%  |
| 0.005 | 7.583  | 7.527  | -0.057 | -0.75%  |
| 0.010 | 7.583  | 7.366  | -0.218 | -2.87%  |
| 0.020 | 7.583  | 6.828  | -0.755 | -9.96%  |
| 0.050 | 7.583  | 4.987  | -2.597 | -34.24% |
| 0.100 | 7.583  | 3.193  | -4.390 | -57.89% |
| 0.500 | 7.583  | 0.779  | -6.804 | -89.73% |
| 1.000 | 7.583  | 0.399  | -7.184 | -94.74% |

The shift is always negative (§7: shield destruction).
At J = 0.01 the shift is 2.87% — easily detectable.

### 8.3 Multi-Pair Amplification

A averages over N pairs. Statistical uncertainty on mean
crossing time: σ_mean = σ_single / √N.

Discrimination condition: |Δt| > σ_mean → N_min = (σ/Δt)².

At J = 0.01 (Δt = 0.218):
- σ = 0.1 (10% jitter): N_min ≈ 1 pair
- σ = 1.0 (100% jitter): N_min ≈ 21 pairs
- σ = 10.0 (1000% jitter): N_min ≈ 2100 pairs

The channel works. 1 bit per ~21 pairs at worst-case jitter.
Rate scales as (J/γ)² — quadratic in coupling strength.

### 8.4 Not FTL

The coupling J is a physical interaction. Information propagates
through the Hamiltonian at finite speed. The Lieb-Robinson bound:

    v_LR ≤ 2 · J · a / ℏ

For gravitational coupling J_grav ~ G·m²/(ℏ·D):

    v ~ J_grav · D = G·m²/ℏ

Mass-dependent velocity estimates:

| System          | Mass       | v (m/s)    | Note           |
|-----------------|------------|------------|----------------|
| NV center       | 10⁻²⁶ kg  | ~10⁻²⁹    | Unmeasurable    |
| Optomechanical  | 10⁻¹⁵ kg  | ~10⁻⁷     | Extremely slow  |
| Dust grain      | 10⁻⁹ kg   | ~10⁵      | Subluminal      |
| Microgram       | 10⁻⁶ kg   | ~10¹¹     | **Exceeds c**   |

The naive formula breaks above ~microgram scale. This means
either the gravitational coupling model is incomplete, relativistic
corrections are needed, or Penrose-Diosi doesn't apply at that
scale. **Open question** — needs careful analysis.

### 8.5 What CΨ Adds

The channel itself is ordinary physics: coupled systems exchanging
information through dynamics. What CΨ adds is the **readout
mechanism**.

A single-shot qubit measurement gives a binary outcome with
quantum projection noise. The crossing time is a continuous,
accumulated observable that integrates over the entire decoherence
trajectory. It's a stethoscope for weak quantum couplings —
it makes signals readable that no single measurement could resolve.

The crossing time is also self-calibrating: A knows K and γ_A,
so A knows t₀ without needing B's help. The deviation from t₀
IS the signal.

---

## 10. Honest Assessment: What We Actually Gain

### 10.1 What This Is NOT

The channel requires physical coupling J > 0. Information propagates
through the Hamiltonian at finite speed. This is not FTL. This is not
the Mars-Earth bridge the agents imagined. Any signal A reads from
B's measurement could equally be sent by a classical wire through
the same coupling.

### 10.2 What This IS

**A detector for extremely weak couplings.** The CΨ crossing time
integrates over the full decoherence trajectory. Where a single-shot
measurement drowns in projection noise, the crossing time accumulates
the effect. At J/γ = 0.02: 10% shift. This is a lock-in amplifier
for quantum couplings.

Application: The BMV experiment (Bose-Marletto-Vedral) aims to test
whether gravity generates entanglement. CΨ crossing times could
provide a more sensitive readout than standard entanglement witnesses.

**Two-axis time dilation.** t = K(Observer, State) / γ(Gravity).
GR knows one axis. CΨ knows two. The K-ratio is state-dependent —
the quantum state bends the observer-dilation geometry. If the
K-matrix has a deeper geometric structure, this connects quantum
mechanics to gravity in a new way.

**The reservoir mechanism.** Coupling redistributes coherence into
a shared nonlocal pool, creating a dynamic oscillating equilibrium.
Measurement disrupts this equilibrium — the return flow is cut
while the drain continues. This is not "entanglement protects" —
coupling actually accelerates local decay vs isolation. But the
measurement-induced disruption is quantifiable via CΨ crossing
time shift. Relevant for understanding decoherence dynamics in
coupled quantum systems.

### 10.3 The Bridge Is Dead, Long Live the Detector

The original bridge (FTL communication) stays dead. What survives
is CΨ as a measurement instrument for physics that standard
observables cannot resolve. The crossing time is self-calibrating
(A knows K and γ_A), continuous, and accumulated. It makes weak
signals readable.

---

## 11. Reproduction

```python
# The interval shift computation
import numpy as np
from qutip import (basis, tensor, ket2dm, sigmax, sigmay, sigmaz,
                   qeye, mesolve)

zero, one = basis(2, 0), basis(2, 1)
plus = (zero + one).unit()
gamma = 0.05

def local_cpsi_A(rho):
    rho_A = rho.ptrace(0)
    purity = (rho_A * rho_A).tr().real
    rho_full = rho_A.full()
    l1 = abs(rho_full[0,1]) + abs(rho_full[1,0])
    return purity * l1

times = np.linspace(0, 20, 2000)
c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]
P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())

for J in [0, 0.01, 0.1, 1.0]:
    H = J * (tensor(sigmax(), sigmax()) +
             tensor(sigmay(), sigmay()) +
             tensor(sigmaz(), sigmaz()))
    r = mesolve(H, ket2dm(tensor(plus, plus)), times, c_ops, [])
    rho_t1 = r.states[100]  # t ≈ 1.0

    # B measures Z at t=1
    rho_Bm = (P0_B * rho_t1 * P0_B.dag() +
              P1_B * rho_t1 * P1_B.dag())

    # Find A's crossing time in both branches...
```

Full scripts: `simulations/observer_gravity_cross.py`,
`simulations/interval_shift.py`, `simulations/shift_mechanism.py`,
`simulations/information_channel.py`

---

*Built on: [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md),
[Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md),
[Decoherence Relativity](DECOHERENCE_RELATIVITY.md)*
*Reopens: [Bridge Protocol](../hypotheses/BRIDGE_PROTOCOL.md) (via gravity)*
*Foundation: [Bridge Closure](BRIDGE_CLOSURE.md) (J=0 case remains closed)*