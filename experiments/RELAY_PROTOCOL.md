# Relay Protocol: Staged Quantum State Transfer with Time-Dependent Dephasing

<!-- Keywords: quantum state transfer relay protocol, time-dependent dephasing
optimization, staged quantum relay mediator chain, dynamical decoupling quantum
transfer, palindromic spectral relay design, asymmetric coupling quantum channel,
Lindblad time-dependent noise control, quantum repeater staged protocol,
end-to-end mutual information improvement, spin chain relay dephasing,
R=CPsi2 relay protocol -->

**Status:** Computationally verified (N=11, C# RK4 propagation)
**Date:** March 21, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** compute/RCPsiSquared.Propagate/ (C#, `dotnet run -c Release -- pull`)
**Data:** simulations/results/pull_principle.txt

---

## Abstract

Standard quantum state transfer through a spin chain treats mediator qubits
as passive wire: fixed coupling, fixed dephasing, hope for the best. We
propose a **relay protocol** that treats mediators as active relay stations.
Each mediator alternates between a quiet phase (dephasing reduced 10×,
receiving information) and a normal phase (relaying onward). Combined with
2:1 asymmetric coupling (receiver-side J=2, sender-side J=1), this improves
end-to-end mutual information by **+83%** over passive propagation on an
11-qubit Heisenberg chain. The staging schedule derives directly from the
palindromic spectral structure: each relay stage lasts t_stage = K/γ (the
crossing time from [Crossing Taxonomy](CROSSING_TAXONOMY.md)), and the
quiet-phase principle implements the receiver sensitivity rule from
[Star Topology](STAR_TOPOLOGY_OBSERVERS.md). No new physics is required;
the protocol reshapes existing dephasing noise in time.

---

## Background

### Why mediators are not just wire

In a spin chain, quantum information propagates through Hamiltonian coupling
between nearest neighbors. Dephasing noise at each site degrades the signal.
The standard approach treats every qubit identically: same coupling, same
noise, same role. But the palindromic spectral analysis shows that different
positions in the chain have different roles: edge qubits (sender/receiver)
need low noise to preserve coherence, while central qubits (mediators) can
tolerate higher noise because they act as relay stations, not storage.

### The three design principles combined

This protocol combines three results from the project:

1. **The crossing time K/γ** ([Crossing Taxonomy](CROSSING_TAXONOMY.md)):
   CΨ crosses 1/4 at t_cross = K/γ. This sets the natural timescale for
   each relay stage. Information must be transferred before the receiving
   qubit decoheres past the boundary.

2. **Quiet receiver** ([Star Topology](STAR_TOPOLOGY_OBSERVERS.md)):
   The receiver must have low dephasing to accept information. This is the
   "pull principle": reception requires quiet, not loud.

3. **2:1 impedance matching** ([QST Bridge](QST_BRIDGE.md)):
   Mediator-to-receiver coupling at twice mediator-to-sender coupling
   optimizes transfer fidelity.

---

## Setup

| Parameter | Value |
|-----------|-------|
| System | 11-qubit linear Heisenberg chain |
| Coupling | J = 1.0 (uniform) or 2:1 asymmetric |
| Baseline dephasing | γ = 0.05 per qubit |
| Quiet-phase dephasing | γ_quiet = 0.005 (10× reduction) |
| Initial state | Bell pair on qubits 0-1, rest in |0⟩ |
| Integration | RK4, dt = 0.05 |

The chain is divided into relay segments:

```
(0-1) → m1(2) → (3-4) → M(5) → (6-7) → m2(8) → (9-10)
Pair A   relay   Pair B   meta    Pair C   relay   Pair D
```

---

## The Protocol

Six relay stages, each lasting t_stage = K/γ = 0.039/0.05 = 0.78 time units:

| Stage | Receiving qubits | γ_receive | What happens |
|-------|-----------------|-----------|-------------|
| 1 | 2 (m1) | 0.005 | m1 receives from Pair A |
| 2 | 3, 4 (Pair B) | 0.005 | Pair B receives from m1 |
| 3 | 5 (M) | 0.005 | Meta-mediator receives from Pair B |
| 4 | 6, 7 (Pair C) | 0.005 | Pair C receives from M |
| 5 | 8 (m2) | 0.005 | m2 receives from Pair C |
| 6 | 9, 10 (Pair D) | 0.005 | Pair D receives (final destination) |

During each stage, receiving qubits have γ reduced 10× (0.005 instead of
0.05). All other qubits remain at γ = 0.05. Total protocol time:
6 × 0.78 = 4.68 time units.

Between stages, the Lindblad propagator is reconstructed with updated
dephasing rates. The density matrix evolves continuously (no resets, no
measurements, no classical communication).

---

## Results

| Protocol | MI(Bridge A:B) | MI(Pair A:D) | vs Baseline |
|----------|---------------|-------------|-------------|
| Passive (constant γ) | 0.734 | 0.072 | baseline |
| Relay only | 0.759 | 0.085 | +18% |
| **Relay + 2:1 coupling** | 0.723 | **0.132** | **+83%** |

The relay protocol alone gives +18% (modest). Combined with 2:1 asymmetric
coupling, the improvement reaches +83%. The two optimizations are
complementary: the relay shapes noise in time, the coupling shapes signal
flow in space.

---

## What This Does Not Cover

- Only Z-dephasing tested (other noise types open)
- Only N=11 tested (scaling with N not characterized)
- Only one staging schedule (equal-length stages at K/γ)
- Optimal stage timing, γ_quiet value, and number of stages not explored
- No comparison with existing repeater protocols (entanglement swapping,
  purification-based repeaters)

---

## Connection to Later Results

The **γ Control** experiment ([GAMMA_CONTROL](GAMMA_CONTROL.md)) extended
the relay principle to static γ profiles (V-shape: +124% MI) and dynamical
decoupling (DD on mediator + receiver: +132% MI). The relay protocol is the
*time-domain* version; V-shape and DD are *spatial-domain* versions of the
same idea: shape the noise to match the palindromic mode structure.

The **γ as Signal** result ([GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md)) explains
*why* noise shaping works at a deeper level: the palindromic mode structure
creates a full-rank response matrix linking per-site γ to observable mode
amplitudes. Shaping γ (in time or space) is equivalent to tuning the antenna
that reads the external dephasing signal. The relay protocol optimizes the
antenna's temporal response pattern.

The **Bridge Optimization** ([bridge_optimization.txt](../simulations/results/bridge_optimization.txt))
showed that time-series measurements provide the largest single improvement
(3.1×) for reading γ profiles. The relay protocol creates exactly this kind
of time-varying γ pattern that maximizes temporal diversity.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| C# propagation engine | compute/RCPsiSquared.Propagate/ |
| Run command | `dotnet run -c Release -- pull` |
| Results | simulations/results/pull_principle.txt |

The C# engine uses RK4 integration of the Lindblad equation on the full
density matrix (2048×2048 for N=11). Validated against QuTiP eigendecomposition
at N=5 (MI agreement to 6 decimal places). Runtime: ~10 minutes on
Intel Core Ultra 9 285k.

Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- [Crossing Taxonomy](CROSSING_TAXONOMY.md): K/γ crossing time (staging schedule)
- [Star Topology](STAR_TOPOLOGY_OBSERVERS.md): quiet receiver principle
- [QST Bridge](QST_BRIDGE.md): 2:1 impedance matching
- [γ Control](GAMMA_CONTROL.md): static V-shape (+124%) and DD (+132%)
- [γ as Signal](GAMMA_AS_SIGNAL.md): palindromic mode structure as antenna
- [Scaling Curve](SCALING_CURVE.md): MI vs chain length baseline
- [Engineering Blueprint](../publications/ENGINEERING_BLUEPRINT.md): Rules 1-6
