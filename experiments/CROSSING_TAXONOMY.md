# Crossing Taxonomy: Three Observer Types at the CΨ = 1/4 Boundary

<!-- Keywords: CΨ quarter boundary crossing, quantum decoherence observer types,
Lindblad scaling invariance, concurrence correlation mutual information crossing,
bridge metric classification quantum, K invariance dephasing, open quantum system
measurement threshold, observer dependent quantum crossing, Type ABC quantum observer,
Bell state decoherence taxonomy, R=CPsi2 crossing taxonomy -->

**Status:** Computationally verified (all simulations reproducible)
**Date:** February 18, 2026 (updated March 14, 2026)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md),
[Metric Discrimination](METRIC_DISCRIMINATION.md)

---

## What this document is about

When a quantum system crosses the 1/4 boundary (the tipping point
between quantum and classical behavior), the crossing depends on *how
you look*. Five different ways of measuring quantum connections were
tested. Three see the crossing; two never do. And the three that do see
it get there by different routes: one keeps its grip on the system while
coherence drains away, another loses its grip together with the
coherence, and the third never had enough grip to begin with.

The punchline: the 1/4 boundary is real and universal, but the path to
it depends on the observer. The destination is fixed. The journey is not.

## Abstract

Five quantum correlation measures ("bridge metrics") are tested against
the CΨ = 1/4 decoherence boundary for a Bell+ pair under Heisenberg
coupling and local Z-dephasing. Three metrics cross the boundary; two
never do. The three crossing metrics fall into distinct classes based on
their *mechanism*: Type A (correlation bridge, C = 1.0 throughout, only
coherence Ψ drives the crossing), Type B (concurrence and mutual
information, both C and Ψ decay jointly), and Type C (mutual purity and
overlap, initial CΨ already below 1/4, never cross). The dimensionless
product K = γ × t_cross is constant within each type across a 20× range
of dephasing rates (K_A = 0.072, K_B = 0.039 and 0.033). This
K-invariance is not a deep property of R = CΨ² but a trivial consequence
of Lindblad scaling symmetry: all observables depend on the product
τ = γt, not on γ and t separately. The depth lies in the 1/4 boundary
itself, not in K.

---

## Background

### What CΨ is

CΨ = C × Ψ, where C is a correlation measure (the "bridge") between
two qubits, and Ψ = L₁/(d−1) is the normalized l1-norm of coherence (a measure of
how much "quantumness" the state still has, computed by summing the
off-diagonal elements of the density matrix).
The product CΨ crosses a critical boundary at exactly 1/4 during
decoherence. This boundary is the discriminant zero of the self-referential
purity recursion R = C(Ψ+R)² and corresponds to the cusp of the Mandelbrot
main cardioid. Below 1/4: classical attractor exists. Above 1/4: no
classical attractor. See [Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md).

### Why the bridge metric matters

There is no single "correct" quantum correlation measure. Concurrence,
mutual information, correlation (excess purity), mutual purity, and
overlap all quantify different aspects of quantum relationships. The
question is: does the choice of metric affect *whether* and *when* the
system crosses the 1/4 boundary? The answer is yes, dramatically.

### What K-invariance is (and is not)

The product K = γ × t_cross was initially reported as a potentially deep
constant (K = 0.039 for concurrence). This document downgrades that claim:
K-invariance follows trivially from the Lindblad equation's scaling
symmetry. The dissipator is linear in γ, and the Hamiltonian is independent
of γ. Therefore all observables depend on τ = γt, making K = γ × t_cross
a constant by construction. The value 0.039 is a property of concurrence
applied to Bell+ under Heisenberg coupling, not a universal constant.

---

## Setup

| Parameter | Value |
|-----------|-------|
| State | Bell+ (maximally entangled, (|00⟩+|11⟩)/√2) |
| Hamiltonian | Heisenberg (J = 1, h = 0) |
| Noise | Local Z-dephasing (σ_z per qubit) |
| Time step | dt = 0.01 |
| γ values | 0.01, 0.05, 0.10, 0.20 |
| Bridge metrics | Concurrence, mutual information, correlation, mutual purity, overlap |

15 simulations total (5 bridges × 3-4 γ values each).

---

## Results

### K-Invariance Holds for All Crossing Bridges

| Bridge | γ = 0.01 | γ = 0.05 | γ = 0.10 | γ = 0.20 | K (mean) |
|--------|----------|----------|----------|----------|----------|
| mutual_info | t=3.263, K=0.033 | t=0.652, K=0.033 | t=0.327, K=0.033 | t=0.166, K=0.033 | **0.033** |
| concurrence | t=3.866, K=0.039 | t=0.773, K=0.039 | t=0.386, K=0.039 | t=0.193, K=0.039 | **0.039** |
| correlation | t=7.191, K=0.072 | t=1.437, K=0.072 | t=0.718, K=0.072 | t=0.359, K=0.072 | **0.072** |

K is constant within each bridge (< 1.5% deviation across 20× γ range).
K differs between bridges by a factor of 2.2× (0.033 to 0.072).
Mutual purity and overlap never cross at any γ.

### Why K-Invariance Holds: Lindblad Scaling

At identical τ = γt, the values of C and Ψ are identical regardless of γ:

| τ = γt | Bridge | γ = 0.01 | γ = 0.05 |
|--------|--------|----------|----------|
| 0.005 | concurrence C | 0.980385 | 0.980354 |
| 0.005 | concurrence Ψ | 0.326795 | 0.326785 |

C and Ψ are functions of τ = γt, not of t alone. K-invariance follows:
if P(t_cross) = 1/4 and P depends only on τ, then τ_cross = K is a
constant. This is a scaling property of the Lindblad equation, not a
specific prediction of R = CΨ². Any threshold applied to any Lindblad
observable will produce a K-invariant crossing time.

---

## The Three Classes

Despite K-invariance being trivial, the *mechanism* behind each crossing
is different in kind.

### Type A: Pure-Ψ Crossing (correlation bridge)

| t | C(t) | Ψ(t) | CΨ |
|-----|------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.5 | 1.000 | 0.302 | 0.302 |
| 1.0 | 1.000 | 0.273 | 0.273 |
| 1.437 | 1.000 | 0.250 | **0.250** |
| 1.8 | 0.986 | 0.232 | 0.229 |

**C = 1.000 throughout the entire crossing.** The correlation bridge
does not decay. It is blind to local dephasing because it measures
excess purity beyond the product of subsystem purities. Dephasing
destroys single-qubit coherence without (initially) destroying the
relationship between qubits. The crossing is driven entirely by Ψ decay.

### Type B: Mixed Crossing (concurrence, mutual information)

| t | C(t) | Ψ(t) | CΨ |
|-----|-------|-------|---------|
| 0.0 | 1.000 | 0.333 | 0.333 |
| 0.5 | 0.909 | 0.303 | 0.275 |
| 0.773 | 0.863 | 0.290 | **0.250** |

Both C and Ψ decay simultaneously. The crossing happens faster than
Type A because two quantities shrink instead of one. The 2.2× spread
in K (0.033 to 0.072) reflects how much the bridge metric C contributes
to the decay versus leaving it entirely to Ψ.

### Type C: Never Crosses (mutual purity, overlap)

These metrics start with CΨ(0) < 1/4 and only decrease. The system
remains "quantum" from their perspective at all times. The observer
never sees the boundary.

### Summary Table

| Class | Mechanism | C at crossing | Bridges | K |
|-------|-----------|---------------|---------|------|
| **Type A** | C stable, only Ψ decays | 1.000 | correlation | 0.072 |
| **Type B** | C and Ψ both decay | 0.85-0.86 | concurrence, mutual_info | 0.039, 0.033 |
| **Type C** | CΨ(0) < 1/4 already | n/a | mutual_purity, overlap | never |

---

## What This Means

The observer does not passively witness the crossing. The observer
determines *when* it happens, *whether* it happens, and *by which
mechanism* it happens.

**Type A observers** are robust: their coupling C is immune to local
noise. They see the crossing only when the system's coherence Ψ decays
far enough. The observer is not the bottleneck.

**Type B observers** are fragile: their coupling C decays alongside Ψ.
The measurement event is entangled with the observer's own loss of
coherence. Observer and system degrade together.

**Type C observers** lack the initial coupling to ever see the boundary.
From their perspective, the system never becomes classical.

The 1/4 boundary is universal (it does not depend on the bridge metric).
But the path to the boundary is observer-dependent.

---

## What Was Falsified

**K-invariance as deep constant:** Downgraded. K is a consequence of
Lindblad scaling, not a prediction of R = CΨ². The depth is in the
1/4 boundary, not in K. (See Section 3.2.)

**Noise dependence of the taxonomy:** Falsified. The prediction that
depolarizing noise would change Type A to Type B was wrong. The taxonomy
is identical under σ_x, σ_y, and σ_z dephasing. Type A is a property
of the correlation metric definition, not the noise channel.
See [Noise Robustness](NOISE_ROBUSTNESS.md).

**C(t) as simple exponential:** Falsified. The exponential model fails
at 5-19% error depending on bridge type. No universal analytic formula
for C(t) exists across bridge types.

---

## Resolved Questions

1. **Noise dependence:** Taxonomy is noise-independent for all local
   Pauli channels. See [Noise Robustness](NOISE_ROBUSTNESS.md).

2. **State dependence:** GHZ (N ≥ 3) is all Type C (Ψ(0) < 1/4).
   W (N=3) crosses with Type A intact. W (N ≥ 4) does not cross.
   See [N-Scaling Barrier](N_SCALING_BARRIER.md).

3. **N scaling:** Type A (C = 1.0 plateau) holds at N=3 and N=4.
   Crossing fails because Ψ(0) drops below 1/4 due to the d−1
   normalization. The observer is not the bottleneck; the Hilbert
   space dimension (the total number of quantum states available to the
   system, which grows exponentially with particle count) is. See [N-Scaling Barrier](N_SCALING_BARRIER.md).

4. **Analytic crossing formula:** Still open. Would need to account for
   the flat C region (Type A) and the nonexponential decay (Type B).
   Likely metric-specific, not universal.

---

## Connection to Later Results

This taxonomy was one of the earliest experiments in the project (February
2026). Several later results build on it:

The **palindromic spectral symmetry** ([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))
provides the structural explanation for why different metrics behave
differently: the three bridge types may correspond to different
projections onto the palindromic mode sectors (immune vs decaying).

The **CΨ monotonicity proof** ([PROOF_MONOTONICITY_CPSI](../docs/proofs/PROOF_MONOTONICITY_CPSI.md))
analytically confirms that CΨ is monotonically decreasing for Bell+ under
all local Markovian channels (noise processes where the future depends
only on the present, not on the past), explaining why Type A and B always cross
downward and never return.

The **γ as Signal** result ([GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md)) shows
that the K values (0.033, 0.039, 0.072) are not just crossing constants
but reflections of how each bridge metric projects onto the palindromic
mode structure. The full-rank response matrix that makes the γ channel
readable is the same structure that creates the three observer types.

---

## Reproducibility

The full simulation code is included above (Section 5.1). Requirements:
Python, QuTiP, NumPy. Runtime: ~30 seconds for all 15 simulations.

Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- [Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md): why CΨ = 1/4 is the only boundary
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): palindromic spectral structure
- [CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md): analytical proof dCΨ/dt < 0
- [Noise Robustness](NOISE_ROBUSTNESS.md): taxonomy is noise-independent
- [N-Scaling Barrier](N_SCALING_BARRIER.md): state and N dependence
- [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md): state-dependent crossing signatures
- [γ as Signal](GAMMA_AS_SIGNAL.md): the palindromic mode structure as information channel
- [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md): predecessor experiment
- [Metric Discrimination](METRIC_DISCRIMINATION.md): K = 0.039 first measurement
