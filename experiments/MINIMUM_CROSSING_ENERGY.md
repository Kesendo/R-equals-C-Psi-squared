# Minimum Energy for CΨ Crossing

**Date**: 2026-03-01
**Status**: Computationally verified (Tier 2)
**Authors**: Thomas Wicht, with Claude (Anthropic)
**Depends on**: OBSERVER_GRAVITY_BRIDGE.md, DYNAMIC_ENTANGLEMENT.md

---

## 1. The Question

What is the minimum energy a quantum system needs to undergo
a CΨ = ¼ crossing? Does the crossing have an energy threshold?

## 2. The Answer

**There is no energy threshold.** The crossing condition is not
about energy but about CΨ_max — the peak CΨ reachable during
the system's evolution. CΨ_max depends on the competition
between Hamiltonian dynamics (J) and decoherence (γ).

---

## 3. Evidence

### 3.1 Same Energy, Different Crossing Behavior

Family: cos(α)|00⟩ + sin(α)|11⟩ under Heisenberg H.

Both |00⟩ and |11⟩ are triplet states with E = +J. Therefore
⟨H⟩ = +J for ALL values of α. The energy is constant.

| α     | ⟨H⟩  | CΨ(0)  | Crosses? |
|-------|------|--------|----------|
| 45°   | J    | 0.3333 | YES      |
| 35°   | J    | 0.2943 | YES      |
| 31°   | J    | 0.2599 | YES (barely) |
| 30°   | J    | 0.2500 | NO       |
| 25°   | J    | 0.1956 | NO       |
| 15°   | J    | 0.0833 | NO       |

Same energy. Different crossing behavior. The boundary is at
α_critical = exactly 30°, where CΨ(0) = exactly 0.25 = ¼.

### 3.2 Zero Initial CΨ, Still Crosses

Family: cos(α)|01⟩ + sin(α)|10⟩ under Heisenberg H.

| α     | CΨ(0)  | CΨ_max | t(max) | Crosses? |
|-------|--------|--------|--------|----------|
| 45°   | 0.3333 | 0.3333 | 0.00   | YES      |
| 25°   | 0.1956 | 0.2956 | 0.37   | YES      |
| 15°   | 0.0833 | 0.3028 | 0.38   | YES      |
| 5°    | 0.0101 | 0.3078 | 0.39   | YES      |

ALL states in this family cross, even at CΨ(0) = 0.01.
The Hamiltonian pumps CΨ upward first. The |01⟩↔|10⟩
exchange interaction generates entanglement rapidly.

### 3.3 Product States: Hamiltonian Creates the Crossing

Starting from CΨ(0) = 0 (zero entanglement):

| State   | CΨ(0) | CΨ_max | Crosses? |
|---------|--------|--------|----------|
| \|0,1⟩  | 0.000  | 0.309  | YES      |
| \|1,0⟩  | 0.000  | 0.309  | YES      |
| \|+,0⟩  | 0.000  | 0.295  | YES      |
| \|0,+⟩  | 0.000  | 0.295  | YES      |
| \|+,1⟩  | 0.000  | 0.295  | YES      |
| \|+,+⟩  | 0.000  | 0.000  | NO       |
| \|0,0⟩  | 0.000  | 0.000  | NO       |
| \|1,1⟩  | 0.000  | 0.000  | NO       |

|0,1⟩ and |1,0⟩ cross despite starting with zero entanglement.
The Heisenberg coupling creates entanglement from the spin
exchange interaction.

|+,+⟩, |0,0⟩, |1,1⟩ never cross because they are eigenstates
(or within a single symmetry sector) of the Hamiltonian. No
dynamics means no entanglement generation.

### 3.4 The Competition: J/γ Ratio

For |0,1⟩ product state, sweeping J/γ:

| J/γ   | CΨ_max | Crosses? |
|-------|--------|----------|
| 0.1   | 0.003  | NO       |
| 0.5   | 0.045  | NO       |
| 1.0   | 0.100  | NO       |
| 2.0   | 0.169  | NO       |
| 5.0   | 0.248  | NO       |
| 10.0  | 0.286  | YES      |
| 20.0  | 0.309  | YES      |
| 50.0  | 0.323  | YES      |
| 100.0 | 0.328  | YES      |

Critical J/γ ≈ 5–10 for this state. Below this ratio,
decoherence wins before the Hamiltonian can build enough
entanglement. Above it, the Hamiltonian wins.

CΨ_max saturates near 1/3 for large J/γ. This is the
maximum CΨ reachable from |0,1⟩ under Heisenberg dynamics
(limited by the state geometry, not by energy).

---

## 4. The Three Regimes

The crossing condition defines three regimes:

**Regime 1: CΨ(0) > ¼**
Starts above threshold. Decoherence drives CΨ downward through ¼.
Crossing always occurs. Time scale: t_cross = K/γ.

**Regime 2: CΨ(0) < ¼ but CΨ_max > ¼**
Hamiltonian pumps CΨ upward first. CΨ rises above ¼, then
decoherence drives it back down through ¼. Crossing occurs
on the downward pass. Requires J/γ above a state-dependent
critical ratio.

**Regime 3: CΨ_max < ¼**
CΨ never reaches the threshold. No crossing. No observer time.
Either: (a) J/γ too small (decoherence wins), or
(b) state is an eigenstate of H (no dynamics at all).

### 4.1 Physical Interpretation

The crossing is not an energy barrier. It is a **coherence
barrier**. The system must accumulate enough entanglement
(measured by CΨ) to reach ¼. The Hamiltonian is the pump.
Decoherence is the drain. The crossing happens when the pump
can fill the pool faster than the drain empties it.

For eigenstates of H: the pump is off. No dynamics, no
entanglement, no crossing, no time. This connects to the
Wheeler-DeWitt ground state (Ĥ|Ψ⟩ = 0, no ticks, no events).

---

## 5. What Determines CΨ_max

CΨ_max depends on four factors:

1. **Initial state geometry**: How much "raw material" for
   entanglement exists. |0,1⟩ has maximal exchange potential.
   |0,0⟩ has none.

2. **Hamiltonian structure**: Which interactions generate
   entanglement. Heisenberg σ·σ is efficient for exchange.
   σ_z⊗σ_z alone (Ising) would behave differently.

3. **J/γ ratio**: Competition between generation and destruction.
   Critical ratio is state-dependent (≈5–10 for |0,1⟩).

4. **Symmetry sector**: States within a single eigenspace of H
   have CΨ_max = CΨ(0). No dynamics to pump.

---

## 6. Reproduction

```python
# Critical angle for cos(α)|00⟩ + sin(α)|11⟩
# Binary search confirms: α_critical = 30.000000°
# At α = 30°: CΨ(0) = 0.250000, sin²(α) = 0.250000

# J/γ sweep for |0,1⟩ product state
# Critical J/γ ≈ 5-10 (CΨ_max crosses 0.25)
```

Full script: `simulations/minimum_energy.py`

---

*Built on: [Observer × Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md),
[Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md)*
*Connects to: [Time as Crossing Rate](../hypotheses/TIME_AS_CROSSING_RATE.md)
(Wheeler-DeWitt: no dynamics → no time)*