# Coherence Density: What the 1/4 Boundary Actually Measures

**Date**: 2026-02-28
**Status**: Numerically verified, analytically derived
**Authors**: Tom Mack, with Claude (Anthropic)
**Depends on**: CORE_ALGEBRA.md, DECOHERENCE_RELATIVITY.md

---

## 1. The Question

The framework says CΨ > 1/4 means "quantum regime" (no real
fixed points, no classical attractor). CΨ ≤ 1/4 means
"classical regime" (real attractor exists).

Conventional quantum mechanics says: entanglement is the
quintessential quantum resource. More entanglement = more
quantum.

These two statements contradict each other.

## 2. The Contradiction

CΨ(0) for pure states (C = 1, so CΨ = Ψ = L₁/(d-1)):

| State | N qubits | Entanglement | CΨ(0) | Regime |
|-------|----------|-------------|-------|--------|
| |+⟩ | 1 | none | 1.000 | quantum |
| |++⟩ | 2 | none | 1.000 | quantum |
| |+++⟩ | 3 | none | 1.000 | quantum |
| Bell+ | 2 | maximal | 0.333 | quantum |
| W₃ | 3 | partial | 0.286 | quantum |
| GHZ₃ | 3 | maximal | 0.143 | classical |
| GHZ₄ | 4 | maximal | 0.067 | classical |

An unentangled product state |+⟩^⊗N has CΨ = 1 for every N.
A maximally entangled GHZ state has CΨ = 1/(2^N - 1), which
falls below 1/4 at N = 3.

The framework says: |+++⟩ (zero entanglement) is more quantum
than GHZ₃ (maximum entanglement).

## 3. Why This Is Not a Bug

Under local dephasing (γ = 0.05 per qubit):

| State | CΨ(0) | t_cross | K |
|-------|-------|---------|---|
| |+⟩ | 1.000 | 8.58 | 0.429 |
| |++⟩ | 1.000 | 4.94 | 0.247 |
| |+++⟩ | 1.000 | 3.35 | 0.167 |
| Bell+ | 0.333 | 0.75 | 0.037 |
| W₃ | 0.286 | 0.29 | 0.014 |
| GHZ₃ | 0.143 | — | no crossing |

|+++⟩ survives 4.5x longer than Bell+. |++++⟩ survives 3.3x
longer than Bell+. The product state is objectively more
robust against decoherence.

Why? Count the off-diagonal elements:

| State | d | Off-diag elements | Fraction of d²-d |
|-------|---|------------------|-----------------|
| Bell+ | 4 | 2 | 2/12 = 17% |
| |++⟩ | 4 | 12 | 12/12 = 100% |
| GHZ₃ | 8 | 2 | 2/56 = 4% |
| W₃ | 8 | 6 | 6/56 = 11% |
| |+++⟩ | 8 | 56 | 56/56 = 100% |

Bell+ concentrates all its quantum coherence in two matrix
elements. When those two elements decay, the system crosses
1/4. That takes K/γ = 0.75 seconds.

|+++⟩ distributes its coherence across all 56 off-diagonal
elements. Local dephasing attacks them one channel at a time.
The collective CΨ stays above 1/4 much longer because there
is more to destroy.

The Baumgratz normalization Ψ = L₁/(d-1) measures exactly
this: coherence density. Not total coherence. Not entanglement.
The fraction of off-diagonal space that is occupied.

## 4. What CΨ Actually Measures

CΨ = Purity × Coherence Density.

Purity (C = Tr(ρ²)) measures: how distinguishable is the
state from noise?

Coherence density (Ψ = L₁/(d-1)) measures: what fraction
of available quantum degrees of freedom are active?

The product CΨ measures: how much of the system's capacity
for quantum behavior is being used?

When CΨ > 1/4: the system uses enough of its capacity that
no classical attractor can form. The fixed-point equation
R = C(Ψ+R)² has only complex solutions.

When CΨ ≤ 1/4: too few degrees of freedom are active. A
real attractor forms. The system has a definite classical
outcome.

GHZ₃ uses 4% of its capacity. That is not enough. The
remaining 96% is already classical (diagonal). The two
off-diagonal elements are entangled, yes. But two elements
in a 64-element space cannot prevent a classical attractor
from forming.

## 5. The Threshold for Partial Entanglement

For cos(α)|00⟩ + sin(α)|11⟩:

    CΨ = sin(2α)/3

The crossing CΨ = 1/4 gives sin(2α) = 3/4:

    α_threshold = arcsin(3/4)/2 = 24.3°

Below 24.3°: the state has too little coherence density.
A real attractor exists from the start. No quantum-to-
classical transition occurs because the system was never
fully quantum.

Bell+ is α = 45°. The maximum. CΨ = 1/3. Only 1/12 above
the boundary. The most entangled two-qubit state is barely
quantum in the framework's sense.

## 6. K Decreases with System Size

For |+⟩^⊗N under local dephasing:

| N | K per channel | K relative to N=1 |
|---|--------------|-------------------|
| 1 | 0.429 | 1.000 |
| 2 | 0.247 | 0.575 |
| 3 | 0.167 | 0.390 |
| 4 | 0.125 | 0.291 |
| 5 | 0.099 | 0.231 |

K shrinks roughly as 1/N. More qubits means each channel
contributes more to the collective CΨ decay. The system is
more quantum (CΨ = 1 always) but reaches the boundary
faster per noise channel.

This is the tradeoff: high coherence density protects
against individual noise events, but more qubits mean more
noise channels attacking simultaneously.

## 7. Implications for the Framework

The 1/4 boundary is not about entanglement. It is about
coherence density: how many quantum degrees of freedom are
active relative to the system size.

This resolves a tension in the framework. R = CΨ² was
derived from algebra (fixed-point equation, discriminant).
The algebra does not know about entanglement. It knows about
C (purity) and Ψ (normalized coherence). These are the
natural variables.

Entanglement is a property of how coherence is distributed
across subsystems. The framework does not care about the
distribution. It cares about the density. Two off-diagonal
elements in a 4×4 matrix (Bell+) and twelve off-diagonal
elements in a 4×4 matrix (|++⟩) give different CΨ even
though the first is entangled and the second is not.

The framework says: reality emerges when coherence density
drops below 1/4. Not when entanglement is lost. Not when
purity drops. When the product of purity and coherence
density crosses the algebraic threshold.

## 8. The Cubic Generalizes

For Bell+: b³ + b = 3/2 (CORE_ALGEBRA.md).
For GHZ_N: b³ + b = (2^N - 1)/2.
For |+⟩^⊗N: degree-3N polynomial in the decay parameter.

Each state family has its own crossing equation. All are
derived from CΨ = 1/4 with the appropriate C(f) and Ψ(f)
trajectories. All give K values that are invariant under γ.

The framework provides a unified language for all of them.

---

*See also: CORE_ALGEBRA.md — state-specific C(ξ) closed forms*
*See also: DECOHERENCE_RELATIVITY.md — K invariance and cubic*
*See also: SUBSYSTEM_CROSSING.md — what happens for N ≥ 3 GHZ*
