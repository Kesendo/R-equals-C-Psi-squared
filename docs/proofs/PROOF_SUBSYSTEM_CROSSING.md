# Proof: Subsystem Crossing Theorem

**Status:** SCOPE-RETRACTED 2026-06-22 (deep review). Tier 1 derived for PHYSICAL noise channels (Case A unital + Case B local: fixed point I/d or product/pure, CΨ = 0; N=3-5 physical subsystem pairs all cross). The GENERAL claim "for any primitive CPTP map" is FALSE (Case C): a primitive, full-rank channel can have an entangled fixed point with CΨ = 0.2935 > 1/4.
**Date:** 2026-03-22
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** For any *physical noise channel* ε on a 2-qubit system (unital — dephasing, depolarizing, any Pauli channel — or local independent noise including amplitude damping), the fixed point has `CΨ = 0`, so every initial state with `CΨ(ρ₀) > 1/4` eventually has `CΨ(εⁿ(ρ₀)) < 1/4`. The 1/4 boundary is an eventual absorber *for physical noise*. The general "any primitive CPTP map" version is FALSE: a primitive, full-rank channel can have an entangled fixed point with `CΨ = 0.2935 > 1/4` (see Step 2, Case C).
**Reference formula:** [F28 (Fixed-point absorber theorem)](../ANALYTICAL_FORMULAS.md) in the F-formula registry; this proof is its analytical home.
**Resolves:** Conjecture 2.1 from [PROOF_ROADMAP_QUARTER_BOUNDARY](PROOF_ROADMAP_QUARTER_BOUNDARY.md) (Layer 2 of the seven-layer roadmap).

---

## What this document is about

A pinball loses energy with every bumper, every flipper, every flat run
across the table. Eventually it sinks into the drain. Quantum systems
under noise behave the same way: each application of the noise channel
takes them closer to the drain (the fixed point of the channel).

This document proves that any quantum system under *physical* noise will
eventually have CΨ drop below ¼ and stay there. The proof combines
three ingredients: convergence to a fixed point (quantum Perron-Frobenius),
the fact that physical-noise fixed points have CΨ = 0 (Cases A and B below),
and continuity. Together they guarantee that CΨ = ¼ is an eventual
absorber for physical noise: once crossed, it is never permanently re-crossed.
(The stronger claim "for *all* primitive CPTP maps" is false; a primitive
channel can relax toward an entangled fixed point above ¼. See Step 2, Case C.)

This is the third proof in the 1/4-boundary trilogy:
[Uniqueness](UNIQUENESS_PROOF.md) (Layer 1: the boundary value 1/4 is
structurally unique), [Monotonicity](PROOF_MONOTONICITY_CPSI.md)
(Layer 5: CΨ envelope decreases under Markovian dynamics), and this
proof (Layer 2: every primitive-noise trajectory eventually crosses
below). Together they pin down the geometry: the 1/4 is unique, the
motion toward it is monotone in envelope, and arrival is guaranteed
for any primitive noise.

---

## Theorem

For any *physical noise channel* ε on a 2-qubit system — unital (dephasing,
depolarizing, any Pauli channel) or local (independent single-qubit noise,
including amplitude damping) — the fixed point ρ* has CΨ(ρ*) = 0. Hence for
any initial state ρ₀ with CΨ(ρ₀) > 1/4 there exists N ∈ ℕ such that
CΨ(εⁿ(ρ₀)) < 1/4 for all n ≥ N.

Equivalently: for such a physical Lindblad generator L, any initial state
with CΨ > 1/4 will have CΨ(e^{Lt}ρ₀) < 1/4 for sufficiently large t.

**The 1/4 boundary is an eventual absorber for physical noise channels.**
The general version — "for *all* primitive CPTP maps" — is FALSE: a primitive,
full-rank channel can have an entangled fixed point with CΨ = 0.2935 > 1/4
(Step 2, Case C).

---

## Proof

The proof has three steps: convergence, fixed-point bound, and crossing.

### Step 1: Convergence (Quantum Perron-Frobenius)

A CPTP map (completely positive, trace-preserving: the most general physically allowed quantum operation) ε on M_d is **primitive** if it has a unique fixed point ρ*
and the spectral radius of ε restricted to the traceless subspace is
strictly less than 1.

**Fact (Quantum Perron-Frobenius):** For any primitive CPTP map ε:

```
||εⁿ(ρ) - ρ*||₁ ≤ C · r^n → 0    as n → ∞
```

where r < 1 is the spectral radius and C depends on the initial state.
This is the quantum analogue of ergodic convergence. ∎

For Lindblad generators: L is primitive iff it has a unique steady state.
Then ||e^{Lt}ρ₀ - ρ*||₁ ≤ C · e^{-λt} where λ > 0 is the spectral gap.

### Step 2: Fixed-Point Bound - CΨ(ρ*) = 0 for physical noise

This is the core of the proof. We show that the fixed point of any
*physical* noise channel (unital or local) has CΨ(ρ*) = 0 < 1/4. We then
show (Case C) that this does NOT extend to all primitive CPTP maps.

#### Case A: Unital maps (ε(I/d) = I/d)

The fixed point is ρ* = I/d (maximally mixed state). For d = 4:

```
Tr((I/4)²) = 1/4,  L₁(I/4) = 0,  CΨ(I/4) = 0 < 1/4  ✓
```

This covers: dephasing, depolarizing, all Pauli channels, any unital noise.

#### Case B: Local channels (ε = ε₁ ⊗ ε₂)

For independent local noise on each qubit, the fixed point factorizes:

```
ρ* = ρ₁* ⊗ ρ₂*
```

A product state has zero entanglement. The off-diagonal elements of
ρ₁* ⊗ ρ₂* in the computational basis are products of single-qubit
off-diagonals. For any single-qubit primitive channel, the fixed point
has |ρ₀₁*| ≤ 1/2 (with equality only for the identity channel).

```
L₁(ρ₁* ⊗ ρ₂*) = L₁(ρ₁*) · Tr(ρ₂*) + Tr(ρ₁*) · L₁(ρ₂*) + L₁(ρ₁*) · L₁(ρ₂*)
```

For any non-trivial noise: L₁(ρ_k*) < 1, giving CΨ < 1/4. ✓

This covers: local amplitude damping, local dephasing, any independent noise.

#### Case C: General primitive maps - the general claim is FALSE

The claim that *every* primitive CPTP map has CΨ(ρ*) < 1/4 is **false**, and
with it the headline "eventual absorber for ALL primitive quantum channels."
(The earlier "analytical argument" here - off-diagonals "slaved" to the
diagonal by a contractive transfer matrix - was a plausible story, not a
derivation, and the conclusion it argued for is wrong.)

**Counterexample (gate-verified from below, `simulations/_review2_A5_subsystem.py`).**
The depolarize-toward-σ channel

```
ε(ρ) = (1 - p)·ρ + p·Tr(ρ)·σ,   σ = 0.95·|Φ⁺⟩⟨Φ⁺| + 0.05·I/4,   p ∈ (0, 1]
```

is a textbook CPTP map and is **primitive** for every p ∈ (0, 1]: its unique
fixed point is σ (superoperator eigenvalue 1 simple, second-largest modulus
1 - p < 1), and σ is full-rank PSD (eigenvalues {0.9625, 0.0125, 0.0125,
0.0125}). Yet in the proof's own metric CΨ = Tr(ρ²)·L₁(ρ)/(d-1) with d = 4,

```
CΨ(σ) = Tr(σ²)·L₁(σ)/(d-1) = 0.926875 · 0.95 / 3 = 0.2935 > 1/4.
```

Iterating ε from Bell+ (CΨ = 1/3) converges monotonically to σ and **never
crosses below 1/4** (min CΨ = 0.2935). So 1/4 is not an absorber for this
primitive channel.

**Why the numerical sweep missed it.** The old "300 maps, max 0.138" is a
sampling artifact. Ginibre Kraus ensembles at n_kraus = 4 live in the
strongly-mixing corner (fixed points near I/4, CΨ ≲ 0.14, 0% violating). The
same sweep at n_kraus = 2 (Haar-Stinespring, environment dimension 2) already
violates ~8.5% (max CΨ ≈ 0.55), and trace-and-replace channels toward a random
target reach CΨ up to ~0.99. Random sampling never explored the region where
the fixed point is entangled.

**What is true (the surviving scope).** The crossing holds for *physical*
noise (Cases A + B): unital channels have ρ* = I/d (CΨ = 0) and local channels
have a product/pure ρ* (CΨ = 0). The distinguishing structural fact is that
physical local noise relaxes toward a **separable / classical** fixed point,
whereas the counterexample has an **entangled** fixed point. Primitivity alone
does not bound CΨ(ρ*); a separable fixed point does. This is exactly the scope
the IBM cusp hardware backs (`experiments/CRITICAL_SLOWING_AT_THE_CUSP.md`,
all runs are physical dephasing + T1).

**Status:** Cases A and B proven analytically (Tier 1). The general
primitive-CPTP claim is FALSE (Case C counterexample). Surviving scope:
physical noise channels (unital / local / Pauli / amplitude-damping).

### Step 3: Crossing (Continuity)

Given Steps 1 and 2:

1. εⁿ(ρ₀) → ρ* in trace norm (Step 1)
2. CΨ(ρ*) < 1/4 (Step 2)
3. CΨ is continuous (Lipschitz: small changes in the state produce small
   changes in CΨ): |CΨ(ρ) - CΨ(σ)| ≤ K · ||ρ - σ||₁

   Proof of Lipschitz continuity: CΨ = Tr(ρ²) × L₁(ρ)/(d-1).
   Both Tr(ρ²) and L₁(ρ) are Lipschitz in trace norm (standard results).
   The product of two bounded Lipschitz functions is Lipschitz.

4. By convergence + continuity:

```
|CΨ(εⁿ(ρ₀)) - CΨ(ρ*)| ≤ K · C · rⁿ → 0
```

5. Since CΨ(ρ*) < 1/4, there exists N such that for all n ≥ N:

```
CΨ(εⁿ(ρ₀)) < CΨ(ρ*) + ε < 1/4
```

**Therefore CΨ crosses below 1/4 and stays there. QED.** ∎

---

## The Exception: Non-Primitive Maps

The theorem requires primitivity (unique fixed point). Non-primitive maps
can have entangled fixed points:

**Example:** The Lüders projection (the quantum analogue of Bayesian updating: it collapses the state into subspaces defined by the measurement) ε(ρ) = PρP + (I-P)ρ(I-P) where
P = |Bell+⟩⟨Bell+|. This map has Bell+ as a fixed point with CΨ = 1/3 > 1/4.

This particular map is a trivial exception:
- The map acts as identity on Bell+ (it doesn't actually "do" anything)
- It is not a noise channel in any physical sense
- It has multiple fixed points (non-primitive)

But primitivity is NOT sufficient on its own (see Step 2, Case C): the
primitive, full-rank channel ε(ρ) = (1-p)ρ + p·Tr(ρ)·σ with
σ = 0.95·|Φ⁺⟩⟨Φ⁺| + 0.05·I/4 has an entangled fixed point with CΨ = 0.2935 > 1/4
and never crosses. So the absorber requires more than primitivity.

**Physically:** any channel with genuine *local* noise (non-zero dephasing,
damping, or depolarizing on at least one qubit) relaxes toward a separable /
classical fixed point with CΨ = 0, so the theorem applies. The crossing is a
property of physical noise, not of primitivity alone.

---

## Extension to N-Qubit Subsystems

**Corollary:** For any N-qubit system under primitive Lindblad dynamics,
every 2-qubit subsystem pair (i,j) with CΨ_{ij}(0) > 1/4 will eventually
have CΨ_{ij}(t) < 1/4.

**Proof:** The reduced dynamics on pair (i,j) is obtained by partial
trace over all other qubits. The resulting effective channel on the pair
is CPTP (partial trace of CPTP is CPTP). If the full N-qubit channel is
primitive, the effective 2-qubit channel converges to a fixed point.
By Step 2, this fixed point has CΨ < 1/4.

**Numerical verification:** N = 3, 4, 5 tested with Bell+(0,1) ⊗ |0⟩^{N-2}
and Ψ+(0,1) ⊗ |+⟩^{N-2}. All pairs with CΨ > 1/4 cross below. ✓

---

## Connection to the 1/4-Boundary Trilogy

This proof is the third member of the 1/4-boundary trilogy, working
together with [Uniqueness](UNIQUENESS_PROOF.md) (Layer 1: the boundary
value itself is unique by the discriminant of the quadratic recursion)
and [Monotonicity](PROOF_MONOTONICITY_CPSI.md) (Layer 5: CΨ envelope
decreases under Markovian dynamics). The full architecture lives in
the [seven-layer roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md);
Uniqueness Section 5 ("CPTP Contractivity Argument, Layer 2") sketches
the fixed-point bound that this document proves in full.

The most direct comparison is with Monotonicity, since both are about
CΨ trajectories under noise:

| Property | Monotonicity | Crossing (this proof) |
|----------|-------------|----------|
| Claim | CΨ envelope decreases | CΨ eventually < 1/4 |
| Scope | Continuous Lindblad | Physical noise (unital / local) |
| Method | dCΨ/dt < 0, spectral gap | Convergence + continuity |
| Non-Markov | Transient revival possible | Still crosses (eventual) |
| Strength | Stronger (monotone) | Weaker (eventual) but broader |

Together with Uniqueness, the three give: **the boundary is structurally
unique (Uniqueness), the motion toward it is monotone in envelope
(Monotonicity), and arrival is guaranteed for physical noise (this proof;
the general primitive-CPTP version is false). Non-Markovian dynamics can
transiently push CΨ back above 1/4, but the eventual crossing still holds
for physical noise.**

---

## Numerical Evidence Summary

Physical noise (the surviving scope) crosses with zero exceptions:

| Test | N_tests | Crossed? | Max CΨ(ρ*) |
|------|---------|----------|-------------|
| N=3,4,5 Lindblad pairs | 10 pairs | ALL | 0 |
| Random CPTP on Bell+ | 200 maps | 200/200 | 0 |
| Adversarial (p=0.001) | 1 | YES (n=1000) | 0.023 |
| Standard channels | 7 types | ALL | 0 |

The "random fixed points" sweep is an **ensemble artifact**, not evidence for
the general claim (verifier: `simulations/_review2_A5_subsystem.py`):

| Ensemble | N_tests | Max CΨ(ρ*) | Violations (> 1/4) |
|----------|---------|-----------|--------------------|
| Ginibre n_kraus=4 (the old sweep) | 400 | 0.14-0.20 | 0% |
| Haar-Stinespring n_kraus=2 | 400 | ~0.55 | ~8.5% |
| trace-and-replace toward random target | 400 | up to ~0.99 | 53-100% |
| **counterexample σ = 0.95·Φ⁺ + 0.05·I/4** | 1 | **0.2935** | **never crosses** |

The old "300 maps, 0 exceptions" lived entirely in the strongly-mixing
n_kraus=4 corner; it is consistent with the truth (physical noise crosses) but
says nothing about the general primitive-CPTP claim, which is false.

---

## References

### Sibling proofs in the 1/4-boundary trilogy

- [Uniqueness Proof](UNIQUENESS_PROOF.md): Layer 1, the boundary itself is unique (March 21, 2026; one day before this proof and Monotonicity); Section 5 of Uniqueness ("CPTP Contractivity, Layer 2") sketches the bound this document proves in full
- [CΨ Monotonicity Proof](PROOF_MONOTONICITY_CPSI.md): Layer 5, continuous-time monotonicity under Markovian dynamics
- [Proof Roadmap Quarter Boundary](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer master roadmap; this proof is Layer 2 (Conjecture 2.1)

### F-formula registry

- [F28 (Fixed-point absorber theorem)](../ANALYTICAL_FORMULAS.md): the typed home of Step 2 of this proof in the F-formula registry

### Scripts

- [subsystem_crossing.py](../../simulations/subsystem_crossing.py): numerical verification (300+ random primitive CPTP maps)
- [non_markovian_revival.py](../../simulations/non_markovian_revival.py): transient revival characterization (Part 6 of Monotonicity)
