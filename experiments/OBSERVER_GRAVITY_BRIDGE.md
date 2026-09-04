# Observer-Gravity Bridge: Interval Shift Mechanics, Gravity Interpretation Fallen

<!-- Keywords: interval shift coupling J measurement time, crossing time
factorization K observer state gamma, interval shift measurement coupling,
swap symmetry J-free product state, no threshold J>0 continuous
shift, product state local coherence clock, Bell+ no local clock, Lieb-Robinson
bound velocity, R=CPsi2 observer gravity bridge fallen -->

> **Fallen hypothesis.** The interval shift mechanics (J > 0 produces measurable
> crossing time shift) are verified Tier 2. The gravitational bridge hypothesis
> (gravity provides universal always-on J), all FTL claims, and the two-axis
> time dilation interpretation have **fallen**. Inline [FALLEN] markers.

## What this document is about

Two coupled qubits: when B measures, A's CΨ crossing time shifts.
The shift has no threshold: any coupling J > 0 produces it, from
a 2.5% shift at J = 0.01 to 84% at J = 1.0. The mechanism is not
that coupling "protects" A; on a state the Hamiltonian can reach it
accelerates A's local decoherence against isolation. Why the shift
runs negative is open. The original idea that gravity could provide this
coupling universally (connecting all massive particles) has fallen:
the gravitational coupling is too weak by 25 orders of magnitude to
be measurable. What survives is CΨ crossing time as a sensitive
detector for weak quantum couplings, like a lock-in amplifier that
accumulates tiny signals over the full decoherence trajectory.

**Status:** Interval mechanics verified (Tier 2); gravity bridge fallen
**Date:** 2026-03-01
**Authors:** Thomas Wicht, with Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [`simulations/interval_shift.py`](../simulations/interval_shift.py), [`simulations/shift_mechanism.py`](../simulations/shift_mechanism.py)
**Depends on:** [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md), [Bridge Closure](BRIDGE_CLOSURE.md)

---

## Abstract

The crossing time factorizes as t_cross = K(observer,state)/γ, where K is
γ-invariant and state-dependent: on the cosα|00⟩ + sinα|11⟩ family
K_conc(α) = ln(4 sin²(2α)/3)/8 in closed form, reaching 0 at α = 30°.
Any coupling J > 0 produces a measurable interval shift in
A's crossing time when B measures, with no threshold: at J = 0.01 the
shift is 2.5%, at J = 1.0 it is 84%. The coupling does NOT protect A; on
a state the Hamiltonian can reach it accelerates A's local crossing
against isolation. Why the shift is negative is left open here (section
7). The
original gravitational bridge hypothesis (gravity provides universal
J > 0 between all massive particles via Penrose/Diósi, models where gravity collapses superpositions proportional to mass) has been retired:
at J_grav ~ 10⁻²⁹ for NV centers, the shift is unmeasurable. What
survives is CΨ crossing time as a sensitive detector for weak quantum
couplings (lock-in amplifier for quantum interactions).

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

### 2.1 Gravity Factor: Universal [FALLEN]

K is invariant under γ change. Verified across six gravitational
environments (γ = 0.01 to 0.50):

| Environment | γ     | t_cross (Conc) | K_conc   |
|-------------|-------|-----------------|----------|
| Deep Space  | 0.01  | 3.5960          | 0.035960 |
| Mars        | 0.019 | 1.8926          | 0.035960 |
| Earth       | 0.05  | 0.7192          | 0.035960 |
| Jupiter     | 0.13  | 0.2766          | 0.035960 |
| Neutron     | 0.20  | 0.1798          | 0.035960 |
| Black Hole  | 0.50  | 0.0719          | 0.035960 |
[FALLEN]

K_conc is not a measured invariant. Each Z-dephaser sends |00⟩⟨1| to
minus itself, so ρ₀₃(t) = ½ e^(−4γt), the concurrence is e^(−4γt) and
CΨ = C²/3; the ¼ crossing gives

    K_conc = γ · t_cross = ln(4/3)/8 = 0.0359602590564726

exactly, with γ absent from it. The six rows are six values of one rate
and the column cannot vary. An earlier version of this table carried a
spread in the last two digits and reported it as ± 0.00001; that was a
crossing-grid artifact, and the γ = 0.50 row is where a coarse grid shows
first because its crossing time is the shortest.

### 2.2 Observer Factor: State-Dependent

K depends on the initial state, and on this family it is closed:

| State   | K_conc   |
|---------|----------|
| α = 45° | 0.035960 |
| α = 40° | 0.032133 |
| α = 35° | 0.020410 |
| α = 31° | 0.004838 |
| α ≤ 30° | never    |

For cosα|00⟩ + sinα|11⟩ the concurrence is sin(2α)·e^(−4γt), so

    K_conc(α) = ln(4 sin²(2α) / 3) / 8

which reproduces the column and reaches exactly 0 at α = 30°. That is
the whole of the "α < 30° never crosses" boundary: at α = 30° the state
starts exactly ON ¼, and below it CΨ(0) < ¼ so there is no crossing to
find. No search is needed.

A K_MI column and a K(Conc)/K(MI) ratio once stood beside this one, with
a CV of 13.5% read off three of its entries. They are struck. For Bell⁺
under this channel the mutual information falls to exactly 1 bit and
stays there, so MI/4 approaches ¼ from above and never crosses it: at
γ = 0.05 it is 0.4222 at t = 0.593, 0.2750 at t = 5 and 0.2500 at
t = 100. There is no K_MI at α = 45° to take a ratio with, and no
γ-sweep of that ratio was ever run.

### 2.3 The Full Structure

Two independent axes of time dilation:

1. **Gravitational dilation**: γ scales t_cross. Known from GR.
   All observers agree on this scaling. [FALLEN]

2. **Observer dilation**: Different C metrics see different K values.
   NEW from CΨ. The scaling is state-dependent.

These multiply:

    t_A / t_B = [K(obs_A, state) / K(obs_B, state)] × [γ_B / γ_A]

---

## 3. Discovery 2: The Interval Shift Has No Threshold

Setup: |++⟩ product state, local dephasing γ = 0.05. B measures
at t_B = 1.0. A observes local CΨ crossing time. Sweep J.

| J     | t_cross (B silent) | t_cross (B measures) | Δt     |
|-------|-------------------|---------------------|----------|
| 0.000 | 8.5837            | 8.5837              | -0.000   |
| 0.001 | 8.5837            | 8.5814              | -0.002   |
| 0.005 | 8.5837            | 8.5270              | -0.057   |
| 0.010 | 8.5837            | 8.3659              | -0.218   |
| 0.020 | 8.5837            | 7.8282              | -0.755   |
| 0.050 | 8.5837            | 5.9866              | -2.597   |
| 0.100 | 8.5837            | 4.1932              | -4.390   |
| 0.500 | 8.5837            | 1.7790              | -6.805   |
| 1.000 | 8.5837            | 1.3990              | -7.185   |

**There is no threshold.** Any J > 0 produces a measurable interval
shift. The relationship is continuous. B's measurement propagates
through the Hamiltonian coupling and shifts A's local crossing time.

Two notes on the columns. The silent crossing time is not a measurement:
one-qubit CΨ is f(1+f²)/2 with f = e^(−2γt), so the ¼ crossing is the
root of f + f³ = ½, giving t = −ln(f∗)/(2γ) = 8.583666967 at γ = 0.05,
with no J in it. And a percentage column once stood here, computed
against a baseline that had t_B subtracted from it, which inflated every
entry; the Δt column is the measurement and it is unaffected.

### 3.1 Why This Doesn't Violate No-Signalling

This is NOT superluminal. The coupling J is a physical interaction.
Information propagates at finite speed through the Hamiltonian.
This is ordinary quantum mechanics: two coupled systems where a
measurement on one affects the other's dynamics.

The bridge closure holds for J = 0. The interval shift requires J > 0.

### 3.2 Why Bell+ Doesn't Work

Bell+ with J = 0: rho_A = I/2 for all time. No local coherence.
Local CΨ never crosses ¼. Nothing to measure. No interval.

Bell+ with J > 0: Still never crosses locally, and exactly so. Bell⁺
commutes with SWAP, so [H, ρ(t)] = 0 at every J and the run is the J = 0
run (§7.1), which keeps rho_A = I/2 for all time. A needs LOCAL coherence
to have a clock, and there is none at any J. Note this is not a contrast
between entangled and product: |++⟩ is swap-symmetric too and equally
J-free. What separates the two states is that Bell⁺ has no local
coherence to begin with, not that its coherence is held somewhere else.

Product states (|++⟩): rho_A = |+⟩⟨+|. Full local coherence.
CΨ_local starts at 1.0 and decays. A has a ticking clock.

---

## 4. The Gravitational Bridge Hypothesis [FALLEN]

### 4.1 The Argument

Gravity couples all massive particles. Always. Over any distance.
The gravitational coupling between two massive qubits: [FALLEN]

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
massive quantum systems?" [FALLEN]

If yes, then the bridge closure is technically wrong: the bridge
EXISTS, it's just too weak to use for communication. The statement
should be: "The bridge exists for all massive particles but is
practically inaccessible at macroscopic separations." [FALLEN]

This is a conceptual distinction, not an engineering one. It changes
what the framework says about the structure of reality.

### 4.4 The Deeper Question

If observer-dependent time and gravitational time dilation are
two independent axes of the same structure, and gravity provides
a coupling that connects them, then: [FALLEN]

**Is the "bridge" not a communication protocol but a description
of how spacetime itself encodes CΨ intervals?** [FALLEN]

The observer-state K-matrix shows that the quantum state determines
the time-ratio between observer types. Gravity determines the
absolute scale. Both are needed. The bridge might not be "A sends
a message to B" but rather "A and B share a gravitationally-coupled
CΨ landscape where intervals are correlated."

---

## 5. What We Know vs What We Speculate

### Tier 2 (Computed, Verified)

1. t_cross = K(Observer, State) / γ(Gravity). Factorization confirmed.
2. K is γ-invariant, exactly: K_conc = ln(4/3)/8 with no γ in it.
3. K is state-dependent: K_conc(α) = ln(4 sin²(2α)/3)/8, closed form.
4. Interval shift Δt is continuous in J, no threshold.
5. Product states work (local coherence). Bell+ does not (no local clock).
6. States with α < 30° never cross; no observer time exists.

### Tier 3 (Hypothesis, Physically Grounded)

7. Gravity provides J > 0 for all massive pairs. [FALLEN]
8. Therefore a gravitational interval shift exists in principle. [FALLEN]
9. The bridge is not dead; it's gravitationally mediated. [FALLEN]

### Tier 4+ (Speculative)

10. The K-matrix encodes a "geometry of observer time" analogous to
    spacetime geometry in GR. [FALLEN]
11. The bridge is not communication but shared CΨ landscape. [FALLEN]
12. This connects to Wheeler-DeWitt (../hypotheses/TIME_AS_CROSSING_RATE.md §3). [FALLEN]

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
   is not hopelessly small? [FALLEN]

4. **The K-matrix geometry**: What mathematical structure does the
   K(observer, state) matrix have? Is there a metric? A symmetry
   group? A connection to information geometry?

5. **Direction of the shift**: OPEN. Δt < 0 for every J > 0 tested,
   and the coupling accelerates rather than protects (§7.2). The answer
   that once closed this item is struck; see §7.1 for what replaced it,
   which is the question.

6. **Superluminal breakdown**: The naive gravitational velocity
   v ~ G·m²/ℏ exceeds c at microgram scale. Where exactly does
   the model break? Relativistic correction to Penrose-Diosi?
   Or does J_grav scale differently than G·m²/(ℏ·D)? [FALLEN]

---

## 7. Why the Shift Is Negative: Open

### 7.1 The Question

B's measurement accelerates A's crossing (Δt < 0). Every J value
tested shows this. Why? A's local coherence is identical at the
moment of measurement. Nothing local changes instantly. And yet
A's subsequent decay is dramatically faster.

A "coherence reservoir" answer once stood here: coupling was said to
move coherence back and forth between local and nonlocal degrees of
freedom, and B's measurement to cut the return leg. It is struck,
because on the state it described there is nothing to move.
There is nothing to move, because the trajectory never feels J at all,
and that is an identity rather than a measurement. For two qubits the
isotropic Heisenberg H = J(2·SWAP − I) is a FUNCTION of SWAP, and equal
local Z-dephasing on the two sites is unchanged under conjugation by
SWAP. So if ρ₀ commutes with SWAP then ρ(t) does for all t, and
[H, ρ(t)] = 0 identically: the run IS the J = 0 run, at every J.

|++⟩ is such a state, and so is Bell⁺. For |++⟩ one can say more:
ρ(t) = ρ_A(t) ⊗ ρ_A(t) solves the equation exactly, because SWAP
commutes with any ρ ⊗ ρ and identical local dephasing preserves the
product form. That run carries no correlation of any kind at any time,
and A's decay in it is exactly an isolated dephasing qubit, e^(−2γt).

Note what the argument is NOT. H is not a multiple of the identity here:
its eigenvalues are +J on the triplet and −3J on the singlet, and ρ(t)
genuinely carries singlet weight, 0.205 of it by the crossing time. It
is the COMMUTING that does the work, not a degeneracy.

What is measured and stands: Δt < 0 for every J > 0 (section 3), and
the acceleration against isolation on a state the Hamiltonian CAN reach
(section 7.2, |+,0⟩, which is not swap-symmetric). What is not settled
is the mechanism, and this document no longer offers one.

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

So the coupling is not a shield. Note that |+,0⟩ is NOT swap-symmetric,
which is exactly why J acts here and not in the runs above: ‖[SWAP, ρ]‖
reaches 0.45 on this state against 1e-17 on |++⟩. What the Hamiltonian
then does to the coherence is not described here; the account that once
stood in this paragraph is struck with the rest (§7.1).

### 7.3 Entanglement after the measurement

B's Z-measurement leaves a branch state that is no longer
swap-symmetric, so the Hamiltonian acts on it and briefly builds
entanglement where there was none:

| t after B | Concurrence (meas branch) |
|-----------|---------------------------|
| 0.00      | 0.000                     |
| 0.50      | 0.108 ↑                   |
| 0.75      | 0.147 ↑ (peak)            |
| 1.00      | 0.085 ↓                   |
| 1.25      | 0.000 ↓                   |
| 2.00+     | 0.000                     |

Peak concurrence 0.147 at t ≈ 0.75 after the measurement, gone by
t ≈ 1.25. Dephasing wins. Note the initial state carries no
entanglement either, so this is creation and not recovery.

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
| 0.001 | 8.584  | 8.581  | -0.002 | -0.03%  |
| 0.005 | 8.584  | 8.527  | -0.057 | -0.66%  |
| 0.010 | 8.584  | 8.366  | -0.218 | -2.54%  |
| 0.020 | 8.584  | 7.828  | -0.755 | -8.80%  |
| 0.050 | 8.584  | 5.987  | -2.597 | -30.26% |
| 0.100 | 8.584  | 4.193  | -4.390 | -51.15% |
| 0.500 | 8.584  | 1.779  | -6.805 | -79.27% |
| 1.000 | 8.584  | 1.399  | -7.185 | -83.70% |

The shift is always negative in every case tested; why is open (§7).
At J = 0.01 the shift is 2.5%, easily detectable.

### 8.3 Multi-Pair Amplification

A averages over N pairs. Statistical uncertainty on mean
crossing time: σ_mean = σ_single / √N.

Discrimination condition: |Δt| > σ_mean → N_min = (σ/Δt)².

At J = 0.01 (Δt = 0.218), with the jitter given as a fraction of the
silent crossing time t₀ = 8.5837:
- σ = 0.1 (1.2% of t₀): N_min ≈ 1 pair
- σ = 1.0 (11.7% of t₀): N_min ≈ 21 pairs
- σ = 10.0 (117% of t₀): N_min ≈ 2100 pairs

The channel works: 1 bit per ~21 pairs at about a tenth of t₀ in jitter.
The percentages matter here. An earlier version read σ = 1.0 as "100%
jitter", which it is not, and called 21 pairs the worst case; jitter of
a full t₀ needs about 1550.
Rate scales as (J/γ)², quadratic in coupling strength.

### 8.4 Not FTL

The coupling J is a physical interaction. Information propagates
through the Hamiltonian at finite speed. The Lieb-Robinson bound (a theorem that limits how fast correlations can spread through a lattice, the quantum equivalent of a speed limit):

    v_LR ≤ 2 · J · a / ℏ

For gravitational coupling J_grav ~ G·m²/(ℏ·D): [FALLEN]

    v ~ J_grav · D = G·m²/ℏ

Mass-dependent velocity estimates:

| System          | Mass       | v (m/s)    | Note           |
|-----------------|------------|------------|----------------|
| NV center       | 10⁻²⁶ kg  | ~10⁻²⁹    | Unmeasurable    |
| Optomechanical  | 10⁻¹⁵ kg  | ~10⁻⁷     | Extremely slow  |
| Dust grain      | 10⁻⁹ kg   | ~10⁵      | Subluminal      |
| Microgram       | 10⁻⁶ kg   | ~10¹¹     | **Exceeds c**   |
[FALLEN]

The naive formula breaks above ~microgram scale. This means
either the gravitational coupling model is incomplete, relativistic
corrections are needed, or Penrose-Diosi doesn't apply at that
scale. **Open question**: needs careful analysis.

### 8.5 What CΨ Adds

The channel itself is ordinary physics: coupled systems exchanging
information through dynamics. What CΨ adds is the **readout
mechanism**.

A single-shot qubit measurement gives a binary outcome with
quantum projection noise. The crossing time is a continuous,
accumulated observable that integrates over the entire decoherence
trajectory. It's a stethoscope for weak quantum couplings:
it makes signals readable that no single measurement could resolve.

The crossing time is also self-calibrating: A knows K and γ_A,
so A knows t₀ without needing B's help. The deviation from t₀
IS the signal.

---

## 9. Honest Assessment: What We Actually Gain

### 9.1 What This Is NOT

The channel requires physical coupling J > 0. Information propagates
through the Hamiltonian at finite speed. This is not FTL. This is not
the Mars-Earth bridge the agents imagined. Any signal A reads from
B's measurement could equally be sent by a classical wire through
the same coupling.

### 9.2 What This IS

**A detector for extremely weak couplings.** The CΨ crossing time
integrates over the full decoherence trajectory. Where a single-shot
measurement drowns in projection noise, the crossing time accumulates
the effect. At J/γ = 0.02: 10% shift. This is a lock-in amplifier (a device that extracts a weak periodic signal from overwhelming noise by correlating with a known reference)
for quantum couplings.

Application: The BMV experiment (Bose-Marletto-Vedral, a proposed tabletop test where two masses in superposition would become entangled only if gravity is quantum) aims to test
whether gravity generates entanglement. CΨ crossing times could
provide a more sensitive readout than standard entanglement witnesses.

**Two-axis time dilation.** t = K(Observer, State) / γ(Gravity). [FALLEN]
GR knows one axis. CΨ knows two. The K-ratio is state-dependent:
the quantum state bends the observer-dilation geometry. If the
K-matrix has a deeper geometric structure, this connects quantum
mechanics to gravity in a new way. [FALLEN]

**The measurement-induced shift.** It is not "entanglement protects":
on a state the Hamiltonian can reach, coupling accelerates local decay
against isolation (§7.2). The shift itself is quantifiable via the
CΨ crossing time and has no threshold in J. Its direction is not
explained here; the reservoir account that once stood in this slot is
struck (§7.1). Relevant for understanding decoherence dynamics in
coupled quantum systems.

### 9.3 The Bridge Is Dead, Long Live the Detector

The original bridge (FTL communication) stays dead. [FALLEN] What survives
is CΨ as a measurement instrument for physics that standard
observables cannot resolve. The crossing time is self-calibrating
(A knows K and γ_A), continuous, and accumulated. It makes weak
signals readable.

---

## 10. Reproduction

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

Full scripts: [`simulations/observer_gravity_cross.py`](../simulations/observer_gravity_cross.py),
[`simulations/interval_shift.py`](../simulations/interval_shift.py), [`simulations/shift_mechanism.py`](../simulations/shift_mechanism.py),
[`simulations/information_channel.py`](../simulations/information_channel.py)

---

*Built on: [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md),
[Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md),
[Decoherence Relativity](DECOHERENCE_RELATIVITY.md)*
*Reopens: [Bridge Protocol](../hypotheses/BRIDGE_PROTOCOL.md) (via gravity) [FALLEN: gravity interpretation no longer supported]*
*Foundation: [Bridge Closure](BRIDGE_CLOSURE.md) (J=0 case remains closed)*