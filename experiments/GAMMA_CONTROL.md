# Practical γ Control: Shaping the Antenna

**Tier:** 2 (computationally verified)
**Date:** March 22, 2026
**Script:** [using_gamma.py](../simulations/using_gamma.py)
**Data:** [using_gamma.txt](../simulations/results/using_gamma.txt)

---

## Summary

Five strategies for controlling how a 5-qubit mediator bridge responds to
dephasing noise. We cannot change γ (the external clock). We can change
our antenna: which qubits are shielded, how the noise profile is shaped,
and whether the gate responds dynamically.

| Strategy | MI(A:B) | vs Baseline | Verdict |
|----------|---------|-------------|---------|
| Baseline (uniform γ=0.05) | 0.338 | 0% | reference |
| V-shape gradient [0.01, 0.03, 0.05, 0.03, 0.01] | 0.755 | +124% | **best static profile** |
| DD on mediator + receiver | 0.784 | +132% | **best overall** |
| AC modulation on gate | 0.338 | 0% | no resonance found |
| γ feedback (kappa=1.0) | 0.326 | -3% | slightly harmful |

---

## What Works

### 1. Quiet edges, loud center (V-shape)

The γ profile [0.01, 0.03, 0.05, 0.03, 0.01] produces the highest MI
among static profiles. The sender and receiver are quiet (low γ), the
mediator carries the noise. This is exactly the transistor architecture:
gate absorbs the clock, source and drain are shielded.

The inverse (loud edges, quiet center) also improves over uniform but
less (+57% vs +124%). The gate needs to be the noisiest node.

### 2. Dynamical decoupling on mediator + receiver

DD on M alone: +27%. DD on receiver alone: +81%. DD on both: +132%.
The effects compound. DD is the hardware implementation of γ shaping.

DD everywhere (all qubits at γ=0.005) gives MI = 1.535, but this is
equivalent to reducing γ globally, not a realistic control strategy.

### Hardware translation

On IBM Torino: select low-T2* qubits for the mediator (natural high γ)
and high-T2* qubits for source/drain. Apply CPMG sequences on receiver
qubits. Toggle DD on/off for the relay protocol.

---

## What Does Not Work

### 3. AC modulation on the gate

Sinusoidal γ_M(t) at frequencies 0.1 to 8.0: no resonance. All
frequencies produce MI within 0.3% of the unmodulated baseline.
The palindromic mode frequencies do not couple to γ modulation.

### 5. State-dependent feedback

Making γ_M depend on the system's own coherence (γ_eff = γ_base * (1 + κ|⟨ZZ⟩|))
slightly degrades performance at all κ values tested (0.1 to 1.0). The
feedback increases γ_M when the system is coherent, which is the opposite
of what helps. Negative feedback (reducing γ_M when coherent) might work
but was not tested.

---

## What Is Useful

### 4. Time-resolved decoder

Doubling γ on qubit 2 at t=5 produces a detectable change in the MI
trajectory: the post-switch decay is steeper. The decoder can monitor
γ changes in real time with ~0.5 time unit resolution.

---

## Design Implications

- **Rule 7 (proposed): Shape the noise profile.** Use V-shape γ: low at
  edges, high at center. This implements the transistor architecture
  passively, without dynamic control.
- **DD is the relay protocol on hardware.** Toggle DD on/off on mediator
  and receiver qubits to implement the staged γ switching from Rule 6.
- **No resonance in γ modulation.** The palindromic structure does not
  couple to AC gate signals. Only DC γ shaping helps.
- **Feedback harms.** Self-tuning transistor does not work with positive
  feedback on |⟨ZZ⟩|.

---

## References

- [Relay Protocol](RELAY_PROTOCOL.md): the staged γ switching protocol
- [Engineering Blueprint](../publications/ENGINEERING_BLUEPRINT.md): Rules 1-6
- [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md): gate mapping
- [Scaling Curve](SCALING_CURVE.md): baseline MI vs chain length
