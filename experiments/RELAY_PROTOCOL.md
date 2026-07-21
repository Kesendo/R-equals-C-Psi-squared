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
**Data:** [`simulations/results/pull_principle.txt`](../simulations/results/pull_principle.txt)

---

## What this document is about

Sending quantum information through a chain of qubits is like passing
a whisper down a line of people in a noisy room. The standard approach
treats every person the same: same volume, same timing. This document
asks: what if we told each person in turn to cup their ear (reduce
noise) right when the whisper reaches them, then go back to normal?

That is the relay protocol. Instead of treating qubits as passive
wire, we treat them as active relay stations that take turns listening.
Combined with a second trick (making the receiving end pull harder than
the sending end pushes, like a 2:1 gear ratio), this improves
information transfer by 83%. No new physics, no extra hardware: just
reshaping when and where the noise is quiet.

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
   Impedance matching is a concept from electrical engineering: for
   maximum power transfer, the receiver's "pull" must match the
   sender's "push." Here, setting mediator-to-receiver coupling at
   twice mediator-to-sender coupling optimizes transfer fidelity.

---

## Setup

| Parameter | Value |
|-----------|-------|
| System | 11-qubit linear Heisenberg chain |
| Coupling | J = 1.0 (uniform) or 2:1 asymmetric |
| Baseline dephasing | γ = 0.05 per qubit |
| Quiet-phase dephasing | γ_quiet = 0.005 (10× reduction) |
| Initial state | Bell pair on qubits 0-1, rest in \|0⟩ |
| Integration | RK4 (Runge-Kutta 4th order, a standard numerical method for evolving differential equations step by step), dt = 0.05 |

The chain is divided into relay segments:

```
(0-1) → m1(2) → (3-4) → M(5) → (6-7) → m2(8) → (9-10)
Pair A   relay   Pair B   meta    Pair C   relay   Pair D
```

---

## The Protocol

Six relay stages, each lasting t_stage = K/γ = 0.039/0.05 = 0.78 time
units (K = 0.039 is the February tool's feedback-model concurrence
value; the exact standard-Lindblad book gives K = 0.036, t_stage = 0.72.
The protocol below was designed and run with 0.78; its +83% result
stands as measured, and the stage length is a design choice, not a
physics constant):

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

Between stages, the Lindblad propagator (the mathematical engine that
evolves the quantum state forward in time, including both coherent
dynamics and noise) is reconstructed with updated dephasing rates. The density matrix evolves continuously (no resets, no
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

The **γ Control** experiment ([γ Control](GAMMA_CONTROL.md)) extended
the relay principle to static γ profiles (V-shape: +124% MI) and dynamical
decoupling (DD on mediator + receiver: +132% MI). The relay protocol is the
*time-domain* version; V-shape and DD are *spatial-domain* versions of the
same idea: shape the noise to match the palindromic mode structure.

The **γ as Signal** result ([γ as Signal](GAMMA_AS_SIGNAL.md)) explains
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

## Seen again through the Q axis (2026-05-31)

We came back to this protocol months later, with the dimensionless ratio
**Q = J/γ₀** in hand, and the whole thing reads differently, sharper, than
it could in March. Three things we could not see then.

**The baseline γ was a guess that turned out to be a measurement.** The setup
table just picks γ = 0.05 as a working value. On 2026-05-29 we extracted the
local dephasing rate of a real chip from an [IBM Kingston](GAMMA0_IS_ALWAYS_THERE.md)
J-scan and found γ₀ ≈ 0.05. The number we wrote into the setup line was, in
hindsight, a prediction of the hardware's actual γ₀.

**The three tricks were one knob.** Impedance matching (J=2 vs J=1), the
quiet phase (γ → γ/10), and the staging clock K/γ are all moves along the
single [Q axis](../docs/Q_REGIME_ANCHORS.md). Dividing each configuration by
γ₀ places the whole protocol at Q ∈ {20, 40, 200, 400}: passive and the
sender side at 20, the 2:1 receiver side at 40, the quiet phase at 200, the
quiet 2:1 receiver at 400. The anchor map's own caveat names J=1 → Q=20 as
"the original deep-quantum baseline from before the framework Q-band structure
was identified." This protocol ran on exactly that default. Every value sits a
full order of magnitude or two above everything that structures the Q axis:
the exceptional point at Q_EP ≈ 1.5–2, the transfer-resonance peak band at
1.2–1.8, the onset at 0.2–0.35.

**That deep-Q placement is not overshoot; it is the only band where transport
exists.** This is the part [The Flow Between Two Singularities](THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)
later made plain. Below the rotation onset (small Q) a single excitation just
diffuses and decays in place, forgetting; it never reaches the far end. Above
it, the excitation sloshes site to site as a wave that propagates and reflects,
remembering. Transport across a chain *is* that slosh, and the slosh lives only
above the EP, seen on Kingston hardware as an excitation that crossed a 3-site
chain to the far end and revived. The single-bond resonance peak at Q ≈ 1.5 is
a local quantity (the birth of the slosh, the EP itself), not arrival at the
destination. To carry an excitation across eleven sites and have it *arrive*
rather than equipartition mid-chain, you must sit deep in the memory regime. So
the relay was twice right before it knew the axis: it guessed γ₀, and it sat in
the only transport band, the deep plateau, before the EP that defines that band
was found.

This reframes the two measured numbers. The quiet phase (Q: 20 → 200) crosses
no regime boundary; it moves from plateau to deeper plateau, extending the
*lifetime* of the carrying memory against decay rather than switching transport
on. That is why it helps only modestly (+18%): the slosh already arrives at
Q=20, and more Q only stretches the coherence, it does not change whether the
wave reaches the end. The +83% lives on a different, orthogonal axis: the 2:1
asymmetry breaks the chain's left-right symmetry and biases *where* the wave
flows, toward the receiver. The time axis (Q) is saturated inside the transport
band; the whole remaining lever is spatial.

One honest new constraint the simulation could not show. The March run was pure
RK4, where the decay envelope is set by T2. On the real chip the envelope is
**Trotterization-limited** (~9 μs from two-qubit gate error at ~0.5% per gate),
not T2-limited (~200 μs). Any hardware version of this protocol must fit its six
stages (6 × 0.78 = 4.68 time units) inside that gate budget, not the far more
generous coherence budget the simulation assumes.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| C# propagation engine | compute/RCPsiSquared.Propagate/ |
| Run command | `dotnet run -c Release -- pull` |
| Results | [`simulations/results/pull_principle.txt`](../simulations/results/pull_principle.txt) |

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
- [The Flow Between Two Singularities](THE_FLOW_BETWEEN_TWO_SINGULARITIES.md): transport is the slosh above the EP (2026-05-31 revisit)
- [γ₀ Is Always There](GAMMA0_IS_ALWAYS_THERE.md): the measured γ₀ ≈ 0.05 behind the setup's guessed baseline
- [Q-Regime Anchor Map](../docs/Q_REGIME_ANCHORS.md): the Q = J/γ₀ axis the protocol sits on
- Main README Section 6: eight engineering consequences
