# Combined Optimization: Mapping + Selective DD

<!-- Keywords: sacrifice zone selective dynamic decoupling combined,
chain selection mode protection IBM Torino, cavity mode localization
synergistic optimization, April 9 hardware test experiment plan,
R=CPsi2 combined optimization -->

**Status:** Tier 2 (simulation of 6 scenarios with real IBM data).
Hardware test planned April 9, 2026.
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
**Script:** [combined_optimization.py](../simulations/combined_optimization.py)
**Data:** [combined_optimization.txt](../simulations/results/combined_optimization.txt)

---

## Abstract

We have three independent optimizations that have never been combined:

1. **Chain selection (mapping):** Choose a chain where a noisy qubit
   sits at the edge. Cost: zero. Benefit: 2.81x mode protection.
2. **Selective DD:** Apply dynamical decoupling only to inner qubits.
   Cost: DD pulses on 3-4 qubits. Benefit: 1.97x (measured on hardware).
3. **Combined:** Both at once.

We simulate all six combinations of chain type (sacrifice vs mean-T2)
and DD strategy (none, uniform, selective) using real IBM Torino data.

Result: **Selective DD adds 1.46x on top of mapping** (synergistic).
The combination achieves 3.72x mode protection efficiency, the highest
of any strategy tested.

---

## The six scenarios

| # | Chain | DD | Σγ | Slowest rate | vs uniform |
|---|-------|----|----|-------------|------------|
| 1 | Mean-T2 | None | 0.048 | 0.0183 | 1.06x |
| 2 | Mean-T2 | All | 0.024 | 0.0091 | 1.06x |
| 3 | Mean-T2 | Selective | 0.033 | 0.0106 | 1.25x |
| 4 | Sacrifice | None | 0.320 | 0.0455 | 2.81x |
| 5 | Sacrifice | All DD | 0.214 | 0.0255 | 3.35x |
| 6 | Sacrifice | **Selective** | 0.290 | 0.0312 | **3.72x** |

**vs uniform** = slowest rate under uniform noise (same Σγ) divided by
slowest rate under this profile. Measures how efficiently the noise
distribution protects modes, normalized for total noise budget.

---

## Key findings

### Selective DD is synergistic with mapping

| Comparison | Factor |
|-----------|--------|
| Mapping only (4) vs uniform (1) | 2.81x |
| Mapping + Selective DD (6) vs Mapping only (4) | 1.46x additional |
| Mapping + Selective DD (6) vs uniform (1) | 3.72x total |

DD on the inner qubits reduces their effective γ (T2echo/T2* = 2.0-2.5x
improvement). Since the protected cavity modes are localized on the
interior (r = 0.994), this directly reduces the noise they see. The
sacrifice qubit keeps its high noise (no DD), maintaining the gradient.

### Mean-T2 chain gains almost nothing from DD strategy

Scenarios 1-3 all show vs_uniform near 1.0 (1.06x to 1.25x). With
uniform noise across the chain, there is no spatial gradient to exploit.
DD reduces the absolute rates (lower Σγ), but the mode protection
efficiency stays flat.

### Absolute rates vs efficiency

The sacrifice chain has higher total noise (Σγ = 0.29-0.32) and
therefore higher absolute decay rates than the mean-T2 chain (Σγ = 0.02-0.05).
The advantage is in **efficiency**: per unit of total noise, the sacrifice
chain extracts 3.72x more protection for the slowest modes.

On hardware, this translates to: the sacrifice chain's protected modes
survive longer PER UNIT OF NOISE than the mean-T2 chain's modes. Whether
this overcomes the higher total noise depends on the observable (SumMI,
fidelity, etc.) and the measurement time window.

---

## IBM experiment plan (April 9, 2026)

### Objective

Verify that sacrifice-zone chain selection + selective DD outperforms
standard practice on IBM Torino hardware.

### Setup

| Parameter | Value |
|-----------|-------|
| Chain A | [85, 86, 87, 88, 94] (sacrifice, Q85 = 3.7 us T2*) |
| Chain B | [18, 89, 19, 90, 60] (mean-T2, all > 140 us T2) |
| DD configs | None, Uniform (all 5), Selective (inner 3-4) |
| Trotter | dt = 0.5 us, steps = [2, 4, 6, 8, 10] |
| Shots | 4000 per config per time point |
| Total | 30 circuits, 120,000 shots |

### Expected results (from simulation)

| Config | vs_uniform |
|--------|-----------|
| Chain A + Selective DD | 3.72x |
| Chain A + No DD | 2.81x |
| Chain A + Uniform DD | 3.35x |
| Chain B + Uniform DD | 1.06x |

### Success criteria

**Primary:** Chain A + Selective DD shows measurably higher SumMI
ratio (vs its own uniform) than Chain B + Uniform DD.

**Secondary:** Selective DD on Chain A outperforms Uniform DD on
Chain A (the 1.46x synergy).

**Fallback:** If DD effects are washed out by gate errors, Chain A
(no DD) vs Chain B (no DD) still tests the mapping prediction.

---

## Predicted time evolution

Full Lindblad simulation (rho(t) = exp(Lt) rho(0), |+>^5 initial state)
reveals the dominant effect is not mode protection but **coupling dynamics**.

### SumMI(t) for all 6 scenarios

| t (us) | 1 T2 | 2 T2+DD | 3 T2+Sel | 4 Sac | 5 Sac+DD | 6 Sac+Sel |
|--------|------|---------|----------|-------|----------|-----------|
| 0.5 | 0.000 | 0.000 | 0.001 | 0.047 | 0.039 | 0.056 |
| 1.0 | 0.000 | 0.000 | 0.001 | 0.121 | 0.106 | 0.149 |
| 1.5 | 0.001 | 0.000 | 0.001 | 0.158 | 0.146 | **0.201** |
| 2.0 | 0.001 | 0.000 | 0.001 | 0.115 | 0.108 | 0.145 |
| 2.5 | 0.000 | 0.000 | 0.001 | 0.071 | 0.066 | 0.089 |
| 3.0 | 0.000 | 0.000 | 0.001 | 0.058 | 0.054 | 0.074 |

### The dominant effect

The mean-T2 chains show **near-zero SumMI** (< 0.001) at all times.
The initial state |+>^5 is a Heisenberg eigenstate: it does not evolve
under H. With very low noise, the chain stays frozen. High purity
(0.79-0.89 at t=2.5) but no information transfer.

The sacrifice chains show **rich oscillatory dynamics** with a peak at
t = 1.5 us. The high noise on Q85 breaks the eigenstate symmetry and
drives Hamiltonian dynamics. The noise is not just damping; it is the
**engine** that creates the dynamics. Without sufficient noise, nothing
happens.

### Selective DD adds ~27% across all times

Sac+Sel (scenario 6) outperforms Sac only (scenario 4) by 1.24-1.34x
at every time point. This is consistent with the eigenvalue prediction
(1.46x spectral protection translates to ~1.27x in observable SumMI).

### What this means for the IBM experiment

The comparison on hardware will not be "sacrifice vs mean-T2 at equal
conditions." It will be "a chain with dynamics vs a chain without."
The sacrifice chain wins not because it has better mode protection
per se, but because the noise gradient creates dynamics that quiet
chains cannot match.

**Script:** [time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py)
**Data:** [time_evolution_6scenarios.txt](../simulations/results/time_evolution_6scenarios.txt),
[time_evolution_plotdata.csv](../simulations/results/time_evolution_plotdata.csv)

---

*See also:*
[Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md) (chain selection),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md) (1.97x measured March 24),
[Resonant Return](RESONANT_RETURN.md) (the sacrifice-zone formula)
