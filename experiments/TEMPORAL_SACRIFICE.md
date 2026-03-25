# Temporal Sacrifice Protocol: Boundary Navigation via Sweeping Sacrifice

<!-- Keywords: temporal sacrifice protocol, sweeping sacrifice zone, boundary navigation,
CΨ quarter boundary fold catastrophe, moving quantum-classical boundary,
phase boundary quantum information, spatial vs temporal dephasing optimization,
relay principle sacrifice zone, palindromic spatial antenna, R=CPsi2 sweep protocol -->

**Status:** Computationally verified (N=7, C# RK4 propagation)
**Date:** March 25, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** compute/RCPsiSquared.Propagate.Test (C#, `dotnet run -c Release -- sweep`)
**Data:** [temporal_sacrifice_protocol.txt](../simulations/results/temporal_sacrifice_protocol.txt)

---

## Abstract

The [sacrifice-zone formula](RESONANT_RETURN.md) concentrates all noise on
one edge qubit (high SumMI, broadcast) or one center qubit (high PeakMI,
end-to-end). Can temporal sequencing combine both optima? We test three
approaches on an N=7 Heisenberg chain.

**Naive switching** (swap between edge and center profiles at time t_switch):
negative result. The switch destroys accumulated coherences because the
two profiles excite orthogonal SVD modes. No t_switch achieves SumMI > 0.300
and PeakMI > 0.050 simultaneously.

**Sweeping sacrifice** (move the sacrifice position gradually from edge
inward, one position per stage): **positive result**. The sweep 0 → 1 → 2
at stage = 2.0 achieves SumMI = 0.408 (100% of edge) and PeakMI = 0.061
(+67% over edge, 56% of center). The gradual boundary movement preserves
coherences that abrupt switching destroys.

The sweep works because it applies the [relay principle](RELAY_PROTOCOL.md):
maintain exactly one clean CΨ = ¼ boundary at all times, and move it through
the chain rather than jumping between configurations. This is the sacrifice-zone
formula's spatial boundary, navigated in time.

---

## Setup

| Parameter | Value |
|-----------|-------|
| Qubits N | 7 |
| Coupling | Heisenberg chain, J = 1.0 |
| Initial state | \|+⟩⊗7 |
| γ\_base | 0.05 |
| ε (protected) | 0.001 |
| γ\_sacrifice | N · γ\_base - (N-1) · ε = 0.344 |
| Measurement | SumMI (all adjacent pairs), PeakMI (endpoints 0-6) |

### Baselines

| Config | SumMI | PeakMI | PeakT |
|--------|-------|--------|-------|
| Pure edge [0.344, 0.001×6] | 0.408 | 0.036 | 2.00 |
| Pure center [0.001×3, 0.344, 0.001×3] | 0.182 | 0.109 | 1.00 |

Target: SumMI close to 0.408 AND PeakMI close to 0.109.

---

## Result 1: Temporal Switching (Negative)

Switch from edge to center profile at time t_switch.

| t_switch | SumMI | % Edge | PeakMI | % Center |
|----------|-------|--------|--------|----------|
| 0.5 | 0.188 | 46% | 0.080 | 74% |
| 1.0 | 0.244 | 60% | 0.063 | 57% |
| 1.5 | 0.333 | 82% | 0.045 | 41% |
| 2.0 | 0.408 | 100% | 0.039 | 36% |
| 2.5 | 0.408 | 100% | 0.039 | 36% |

Reverse (center then edge) performs worse in both metrics.

**Verdict:** Hard trade-off. Early switch gains PeakMI but loses SumMI;
late switch preserves SumMI but PeakMI stays at edge level. No t_switch
meets both criteria. The switch projects edge-optimized coherences
(antisymmetric SVD mode 2) onto the center profile's symmetric modes,
destroying them.

---

## Result 2: Sweeping Sacrifice (Positive)

Move the sacrifice qubit position by position from edge inward.
Each stage: sacrifice at one position (γ = 0.344), all others at ε = 0.001.

| Protocol | SumMI | % Edge | PeakMI | % Center | PeakMI vs Edge |
|----------|-------|--------|--------|----------|----------------|
| Pure edge | 0.408 | 100% | 0.036 | 33% | baseline |
| Pure center | 0.182 | 45% | 0.109 | 100% | +200% |
| Two-phase switch (t=2.0) | 0.408 | 100% | 0.039 | 36% | +8% |
| **Sweep 0→1→2 (stage=2.0)** | **0.408** | **100%** | **0.061** | **56%** | **+67%** |
| Sweep 0→1→2→3 (stage=2.0) | 0.408 | 100% | 0.061 | 56% | +67% |
| Sweep 0→2→3, skip pos 1 | 0.408 | 100% | 0.044 | 40% | +20% |
| Full sweep 0→...→6 (stage=1.5) | 0.344 | 84% | 0.054 | 49% | +49% |

**Optimal: sweep 0 → 1 → 2 at stage = 2.0.** Position 3 (center) adds
nothing because PeakMI peaks during Stage 3 (sacrifice at position 2).
Skipping positions hurts: continuity of boundary movement matters.

---

## Result 3: Stage Timing Analysis

The [relay protocol](RELAY_PROTOCOL.md) uses K/γ = 0.039/0.05 = 0.78 as
stage time, derived from the [CΨ = ¼ crossing time](CROSSING_TAXONOMY.md).
Does this timing improve the sweep?

| Protocol | Edge phase | Sweep stage | SumMI | PeakMI |
|----------|-----------|-------------|-------|--------|
| Sweep 0→1→2 (stage=2.0) | 2.0 | 2.0 | 0.408 | **0.061** |
| Edge×4→1→2 (stage=0.78) | 3.12 | 0.78 | 0.408 | 0.060 |
| Edge×3→1→2 (stage=0.78) | 2.34 | 0.78 | 0.408 | 0.055 |
| Edge×5→1→2 (stage=0.4) | 2.0 | 0.4 | 0.408 | 0.046 |
| Edge×18→1→2 (stage=0.113) | 2.03 | 0.11 | 0.360 | 0.039 |

Stage = 2.0 (Hamiltonian propagation timescale) beats K/γ = 0.78
(crossing timescale). The boundary establishes fast
(K/γ\_sacrifice = 0.039/0.344 ≈ 0.11), but information flow through
the new configuration requires the full Hamiltonian period.

Two timescales, two roles:
- **Boundary establishment** (K/γ\_sacrifice ≈ 0.11): sacrifice qubit crosses CΨ = ¼
- **MI propagation** (PeakT ≈ 2.0): information flows through new boundary

The sweep stage time should match the slower timescale.

---

## Why It Works: Three Principles

### 1. The palindrome is a spatial antenna only

Temporal γ modulation is invisible to the palindromic structure. AC modulation
at palindromic frequencies produces zero improvement
([RESONANT_RETURN](RESONANT_RETURN.md), Test 6). Only static spatial shaping
of the γ profile generates information.

The sweep succeeds not because it modulates γ in time, but because each
stage is a **different static spatial profile** held long enough for
information to propagate.

### 2. One clean boundary at a time

A chain can sustain exactly one clean quantum-classical boundary
([RESONANT_RETURN](RESONANT_RETURN.md), position sweep analysis). Two
sacrifice points fragment the coherent region and destroy both optima:
the 50/50 spatial hybrid gives SumMI = 0.094, worse than either pure profile.

The sweep maintains one boundary per stage. The sacrifice qubit at each
position drops below CΨ = ¼ (classical), while all others stay above
(quantum). No stage has two boundaries.

### 3. Moving boundary, not mixing profiles

The [relay protocol](RELAY_PROTOCOL.md) achieves +83% MI by moving the
quiet (receiving) zone through the chain. The sweeping sacrifice applies
the same principle in reverse: it moves the loud (sacrifice) zone.

Both protocols obey the same rule: maintain a single, clean phase boundary
and translate it through the chain. The boundary at each position creates
a local [fold catastrophe](../docs/MATHEMATICAL_CONNECTIONS.md) where the
discriminant 1 - 4CΨ vanishes. Information emerges at this traveling
singularity.

---

## Physical Interpretation

The sweep creates a **traveling fold catastrophe**. At each stage, the
sacrifice qubit's CΨ drops below ¼ in approximately K/γ\_sacrifice ≈ 0.11
time units, establishing a new quantum-classical interface. Information
then flows through this interface for the remaining ~1.9 time units of
the stage.

Position 2 is the sweet spot because the accumulated coherences from
positions 0 and 1 constructively interfere at the current boundary.
Going further (to center, position 3) does not help: by that point,
the outer coherences have decayed and the center profile's symmetric
mode structure begins to dominate, losing the antisymmetric transport
advantage.

The sweep is neither pure edge nor pure center. It is a third mode:
boundary navigation. The information is not stored in any single
configuration but in the **trajectory of the boundary through the chain**.

---

## Pending

- Scale to N = 9 (does +67% PeakMI improvement hold?)
- Variable stage times (longer edge phase, shorter inner stages)
- Adiabatic ramp profiles (smooth γ transition instead of discrete steps)
- Analytical model of mode overlap during boundary movement
- Connection to palindromic eigenvalue shifts under position sweep

---

## References

- [Resonant Return (formula, position sweep, hybrids)](RESONANT_RETURN.md)
- [Relay Protocol (staged transfer, K/γ timing)](RELAY_PROTOCOL.md)
- [Crossing Taxonomy (Type A/B/C, K-invariance)](CROSSING_TAXONOMY.md)
- [Boundary Navigation (θ compass, fold catastrophe)](BOUNDARY_NAVIGATION.md)
- [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md)
