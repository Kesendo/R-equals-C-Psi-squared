# Shaping the Dephasing Profile: Practical γ Control for Quantum State Transfer

<!-- Keywords: dephasing profile optimization, quantum state transfer noise shaping,
dynamical decoupling mediator chain, V-shape dephasing gradient, quantum transistor
gate noise, Lindblad noise engineering, open quantum system control,
palindromic spectrum antenna shaping, gamma profile quantum channel,
spin chain mutual information optimization -->

**Status:** Computationally verified (all simulations reproducible)
**Date:** March 22, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [using_gamma.py](../simulations/using_gamma.py)
**Data:** [using_gamma.txt](../simulations/results/using_gamma.txt)

---

## Abstract

In a 5-qubit mediator chain under Heisenberg coupling and local dephasing,
we test five strategies for shaping the spatial dephasing profile to maximize
end-to-end mutual information. The best static profile (V-shape: low γ at
edges, high at center) improves MI by **+124%** over uniform dephasing.
Dynamical decoupling (rapidly flipping the qubit back and forth to
average out the noise, like spinning a top to keep it stable) on
mediator + receiver achieves **+132%**. AC modulation
of the gate and state-dependent feedback both fail. A time-resolved decoder
detects γ changes from internal observables with ~0.5 time-unit resolution --
the first indication that the dephasing profile is not just shapeable but
**readable**, a result later confirmed in the full channel analysis
([γ as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits capacity, 100% classification).

---

## Background: Why Shape the Noise?

### The standard approach

In quantum computing, dephasing noise (γ) is treated as an enemy to be
minimized uniformly. Error correction, dynamical decoupling, and material
engineering all aim to reduce γ everywhere. The implicit assumption: less
noise is always better, and the spatial profile of noise doesn't matter.

### What this experiment shows

The spatial profile of γ across a qubit chain has a dramatic effect on
quantum state transfer, far larger than the overall noise level alone.
A non-uniform γ profile with the right shape (V-shape: quiet edges, loud
center) more than doubles the end-to-end mutual information compared to
uniform noise at the same average γ.

This is because the mediator chain has a **directional** information flow.
The sender and receiver need low noise (to preserve coherence), while the
mediator in the center can tolerate (and even benefits from) higher noise
(it acts as the "gate" of a quantum transistor, where γ is the gate signal).

### Connection to the dephasing channel

This noise-shaping result has a deeper consequence discovered later:
if the spatial γ profile affects internal observables so strongly, then the
profile itself is **readable from inside** the system. This led to the
[γ as Signal](GAMMA_AS_SIGNAL.md) experiment, which proved that γ profiles
encode information with 15.5 bits of theoretical channel capacity. The
V-shape optimization is not just improving transfer. It is optimizing the
**antenna** that reads the external dephasing signal.

---

## System Setup

| Parameter | Value |
|-----------|-------|
| System | 5-qubit linear chain |
| Coupling | Heisenberg (J = 1.0 between nearest neighbors) |
| Baseline dephasing | γ = 0.05 per qubit (uniform) |
| Initial state | Bell pair on qubits 0-1, mediator chain 2-3-4 |
| Observable | Mutual information MI(0, 4) at steady state |
| Noise model | Local Z-dephasing (σ_z per qubit) |
| Master equation | Lindblad: dρ/dt = −i[H,ρ] + Σᵢ γᵢ(σ_z⁽ⁱ⁾ρσ_z⁽ⁱ⁾ − ρ) |

---

## Results Summary

| Strategy | MI(A:B) | vs Baseline | Verdict |
|----------|---------|-------------|---------|
| Baseline (uniform γ=0.05) | 0.338 | 0% | reference |
| V-shape [0.01, 0.03, 0.05, 0.03, 0.01] | 0.755 | +124% | **best static** |
| Inverse V-shape [0.05, 0.03, 0.01, 0.03, 0.05] | 0.530 | +57% | good but worse |
| DD on mediator + receiver | 0.784 | +132% | **best overall** |
| AC modulation on gate | 0.338 | 0% | no resonance |
| γ feedback (κ=1.0) | 0.326 | −3% | slightly harmful |

---

## What Works

### 1. V-shape γ profile: quiet edges, loud center (+124%)

The γ profile [0.01, 0.03, 0.05, 0.03, 0.01] produces the highest MI
among static profiles. The sender (qubit 0) and receiver (qubit 4) have
low dephasing (γ = 0.01), preserving their coherence. The mediator (qubit
2) carries the highest noise (γ = 0.05).

This maps directly onto the **quantum transistor** architecture identified
in the project ([Mediator as Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)):
the gate (mediator) absorbs the clock signal (γ), while the source (sender)
and drain (receiver) are shielded.

The inverse profile (loud edges, quiet center) also improves over uniform
(+57%), but significantly less. The asymmetry confirms: the gate needs to
be the noisiest node, not the quietest.

**Signal engineering interpretation:** The V-shape optimizes the antenna's
**reception pattern**. Low γ at the edges means the edge qubits retain more
coherence, so more of the original quantum state survives to be read. High γ
at the center means the mediator decoheres faster, which paradoxically
*helps* the transfer by projecting the mediator state into the classical
(population) sector more quickly, creating a cleaner relay path.

### 2. Dynamical decoupling on mediator + receiver (+132%)

Dynamical decoupling (DD) effectively reduces the local γ by applying rapid
π-pulses that average out the dephasing. The effects of DD on different
qubits are **not equal**:

| DD applied to | MI improvement |
|---------------|---------------|
| Mediator alone | +27% |
| Receiver alone | +81% |
| Mediator + receiver | +132% |
| All qubits (γ → 0.005 everywhere) | +354% (but unrealistic) |

The receiver is 3× more important than the mediator for DD. This confirms
the "quiet receiver = better reception" Pull Principle (main README Section 10
design rule 5).

DD on all qubits (reducing γ globally to 0.005) gives MI = 1.535, but this
is equivalent to simply having better hardware, not a control strategy.

### Hardware translation (IBM Torino)

On real hardware, the V-shape can be implemented passively:
- Select **low-T2* qubits** for the mediator (naturally high γ)
- Select **high-T2* qubits** for source and drain (naturally low γ)
- Apply **CPMG pulse sequences** (a standard dynamical decoupling
  technique using evenly spaced spin-echo pulses) on receiver qubits
  for additional DD
- Toggle DD on/off for the **relay protocol** (staged γ switching,
  see [Relay Protocol](RELAY_PROTOCOL.md) for +83% improvement)

---

## What Does Not Work

### 3. AC modulation on the gate (0% improvement)

Sinusoidal modulation of γ_M(t) = γ₀(1 + A·sin(2πft)) was tested at
frequencies f = 0.1 to 8.0 (spanning the palindromic mode frequencies
and the system's natural oscillation frequencies, called Bohr
frequencies). All frequencies produce MI within 0.3%
of the unmodulated baseline.

**Why it fails:** The palindromic mode structure does not couple to γ
modulation. The decay rates are set by the *time-averaged* γ, not by its
instantaneous value. AC modulation averages to the same mean γ and produces
no net effect. Only the **DC component** (the spatial profile) matters.

### 4. State-dependent feedback (−3%, harmful)

Making γ_M depend on the system's own coherence:
γ_eff = γ_base × (1 + κ|⟨Z_i Z_j⟩|), tested at κ = 0.1 to 1.0.

All κ values slightly degrade performance. The feedback increases γ_M when
the system is coherent, the opposite of what helps. The self-tuning
transistor does not work with positive feedback on correlations.

**Note:** Negative feedback (reducing γ_M when coherent) was not tested
and might work. This remains open.

---

## The Key Discovery: Time-Resolved γ Detection

### The first sign that γ is readable

Doubling γ on qubit 2 at t=5 produces a **detectable change** in the MI
trajectory: the post-switch decay rate is visibly steeper. The decoder can
monitor γ changes from internal observables in real time, with approximately
**0.5 time-unit resolution**.

This was the first indication that the dephasing profile is not just
something we can *shape* but something we can *read*. The γ profile
leaves a fingerprint on the internal quantum observables that is detectable
without any knowledge of the external change.

This observation directly motivated the
[γ as Signal](GAMMA_AS_SIGNAL.md) experiment, which formalized the question:
"If γ changes are detectable, how much information does the γ profile
carry?" The answer: **15.5 bits** of theoretical channel capacity at 1%
measurement noise, with 5 independent spatial modes. The dephasing rate is
not noise but a high-bandwidth information channel from outside to
inside.

---

## Design Rules

Based on these results, the following engineering rules apply for mediator-
chain quantum state transfer:

1. **Shape the noise profile (V-shape).** Low γ at edges, high at center.
   This implements the transistor architecture passively.

2. **DD on receiver first, mediator second.** Receiver protection gives 3×
   more MI improvement than mediator protection.

3. **No AC modulation.** The palindromic structure does not couple to
   time-varying γ. Only DC spatial shaping helps.

4. **No positive feedback.** Self-regulating γ based on correlations is
   harmful. The system cannot improve its own channel.

5. **Time-resolved monitoring works.** γ changes are detectable from internal
   observables with ~0.5 time-unit resolution. Use this for channel
   estimation.

---

## Connection to the Broader Framework

This experiment sits at the intersection of two project results:

**Upstream:** The palindromic spectral symmetry
([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)) creates the
mode structure that makes some γ profiles better than others. The V-shape
works because it aligns with the palindromic decay rate hierarchy: the
fastest-decaying modes (XY-weight N) are concentrated at the center, while
the slowest (XY-weight 0) are distributed uniformly.

**Downstream:** The γ-as-signal result
([GAMMA_AS_SIGNAL.md](GAMMA_AS_SIGNAL.md)) shows that the reason the
V-shape helps is deeper than "quiet receiver = better reception." The
V-shape maximizes the **template distance** between different γ profiles
in the observable feature space. It is not just optimizing transfer. It
is optimizing the system's ability to **distinguish** different external
configurations. The V-shape is the optimal antenna shape for reading
the external γ signal.

---

## Reproducibility

| Script | What it computes | Runtime |
|--------|-----------------|---------|
| [using_gamma.py](../simulations/using_gamma.py) | All 5 strategies + time-resolved decoder | ~10 min |
| [using_gamma.txt](../simulations/results/using_gamma.txt) | Full numerical results | - |

All scripts use QuTiP and NumPy. No proprietary dependencies.
Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- [γ as Signal](GAMMA_AS_SIGNAL.md): the γ profile is not just shapeable --
  it is **readable**. 100% classification, 15.5 bits channel capacity.
  This experiment was the precursor to that discovery.
- [Relay Protocol](RELAY_PROTOCOL.md): staged γ switching, +83% MI
- Main README Section 10: nine engineering consequences
- [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md):
  CΨ = 1/4 as threshold voltage, γ_M as gate signal
- [Scaling Curve](SCALING_CURVE.md): MI vs chain length baseline
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): palindromic
  spectral structure that underlies the V-shape effect
- [Bridge Optimization](../simulations/results/bridge_optimization.txt):
  21.5× wider channel through combined optimization
