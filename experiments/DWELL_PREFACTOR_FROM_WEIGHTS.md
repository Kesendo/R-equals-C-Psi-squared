# Dwell-Time Prefactor from Pauli Sector Weights

**Status:** Partial success (April 5, 2026)
**Script:** [dwell_prefactor_from_weights.py](../simulations/dwell_prefactor_from_weights.py)
**Predecessors:**
[BOUNDARY_NAVIGATION](BOUNDARY_NAVIGATION.md) (dwell time),
[PRIMORDIAL_SUPERALGEBRA_CAVITY](PRIMORDIAL_SUPERALGEBRA_CAVITY.md) (face-swap algebra)

---

## What this document is about

A quantum state can be described in two very different ways. One way is static: you write down how much of the state lives in each type of Pauli component (the four basic quantum operators I, X, Y, Z and their combinations), and you get a list of numbers that describe the state at one moment in time, like a snapshot. The other way is dynamic: you let the state evolve under decoherence, watch the numbers change, and derive quantities from the trajectory, like a film.

For a long time these two views were used in parallel in this project. The static view produced things like the light-lens decomposition, which says every quantum mode is part structure (the lens face, immune to decoherence) and part signal (the light face, absorbed by decoherence). The dynamic view produced things like the dwell time at the CÎ¨ = 1/4 crossing: how long the state lingers near the fold before falling through.

The question was whether these two views carry the same information. Specifically: the dwell time near the fold has a prefactor, a fixed number that scales how the duration depends on how close you pass to the boundary. For the Bell+ state under Z-dephasing, that number is 1.080088. It was computed from the dynamic equations. The question is whether the same number can be read off directly from the static Pauli weights, without ever solving the dynamics.

The answer is half yes, half no, and the half-line is precise. For states whose Pauli content sits only in even-weight sectors (like Bell+, which has no single-qubit terms), the dwell prefactor is a pure algebraic function of one static quantity: the weight in the light face at the crossing moment. Static and dynamic are the same object. For states with odd-weight content (like the simple product state of two qubits in superposition), the dwell prefactor needs more information than the sector weights alone. It needs the actual sizes of individual Pauli coefficients, which the sector weights, being sums of squares, have already lost.

This is a precise result about where algebra replaces dynamics and where it does not. It is also a worked example of the project's broader claim that static and dynamic descriptions are often two views of the same structure, and that the useful work is figuring out exactly when they diverge.

---

## Abstract

For Bell+ under Z-dephasing, the dwell-time prefactor at the CΨ = 1/4
crossing is expressible as a **pure function of the k = 2 sector weight**:

    prefactor = (2 + 4W₂) / (1 + 6W₂)

where W₂ is the light-face ({X,Y}) sector weight at the crossing moment.
This gives 1.080088, matching the direct computation exactly. The static
face-swap algebra and the dynamic cusp passage are the same object for
Bell+: an algebraic identity, not an analogy.

For states with odd-weight Pauli content (like |+⟩^{⊗2}), the prefactor
additionally requires individual Pauli coefficient magnitudes (via √W₁),
not just sector weights. The extra structure is the sign pattern of Pauli
coefficients, which determines L1-coherence but is not reconstructible
from sector weights alone.

---

## 1. Bell+ Pauli Decomposition Under Z-Dephasing

Bell+ = (|00⟩ + |11⟩)/√2 has the Pauli expansion:

    ρ(0) = (1/4)(I⊗I + X⊗X - Y⊗Y + Z⊗Z)

Under Z-dephasing at rate γ per qubit, each Pauli coefficient decays as
e^{-2γ·w(P)·t} where w(P) is the XY-weight (number of X or Y factors):

    ρ(t) = (1/4)(I⊗I + f·X⊗X - f·Y⊗Y + Z⊗Z),    f = e^{-4γt}

The state has **only k = 0 and k = 2 content**: {I⊗I, Z⊗Z} are stationary
(lens-face), {X⊗X, Y⊗Y} decay together (light-face). No k = 1 terms.

Verified numerically: purity and coherence from the Pauli decomposition
match the closed forms C = (1+f²)/2 and Ψ = f/3 to machine precision.

---

## 2. CΨ(t) from Sector Weights

The sector weights are:

    W₀ = (1/4)(a²_{II} + a²_{ZZ}) = 1/2    (constant)
    W₂ = (1/4)(a²_{XX} + a²_{YY}) = f²/2    (decays as e^{-8γt})

Purity: C = W₀ + W₂ = 1/2 + f²/2
Coherence: Ψ = f/3 = √(2W₂)/3

The CΨ product involves both the weight W₂ (through purity) and the
square root √W₂ (through coherence). The state-specific relation is
W₂ = 9Ψ²/2, or equivalently Ψ = √(2W₂)/3.

---

## 3. The Weight-Based Derivative

Each sector weight decays at a rate proportional to its XY-weight:

    dW_k/dt = -4γk·W_k

For Bell+ (W₁ = 0):

    dC/dt = dW₂/dt = -8γW₂
    dΨ/dt = -4γΨ    (since Ψ = f/3 and df/dt = -4γf)

Applying the product rule:

    dCΨ/dt = dC/dt · Ψ + C · dΨ/dt
           = (-8γW₂)·Ψ + (1/2 + W₂)·(-4γΨ)
           = -4γΨ(2W₂ + 1/2 + W₂)
           = **-2γΨ(1 + 6W₂)**

This is the weight-based formula for dCΨ/dt, verified to match
Formula 25's derivative to 14 decimal places.

### At the crossing (CΨ = 1/4)

Since CΨ = C·Ψ = 1/4, we have Ψ = 1/(4C) = 1/(2 + 4W₂). Substituting:

    |dCΨ/dt|_cross = 2γ(1 + 6W₂)/(2 + 4W₂)

And the dwell-time prefactor becomes:

    **prefactor = (2 + 4W₂) / (1 + 6W₂)**

For Bell+: W₂ at the crossing is 0.3709 (from f_cross = 0.8612), giving
prefactor = 3.4834/3.2256 = **1.080088**, matching the direct computation
exactly.

---

## 4. Why This Works for Bell+ and Not Generally

### Bell+: success (algebraic identity)

Bell+ has no k = 1 Pauli content. The coherence Ψ = f/3 involves only
the k = 2 coefficient, and the relationship Ψ = √(2W₂)/3 connects Ψ
directly to the sector weight. The derivative dΨ/dt = -4γΨ follows from
the exponential decay of the k = 2 sector alone.

### |+⟩^{⊗2}: partial (needs coefficient magnitudes)

|+⟩^{⊗2} = |+⟩⊗|+⟩ has Pauli content across k = 0, 1, and 2:

    ρ(t) = (1/4)(I⊗I + g·I⊗X + g·X⊗I + g²·X⊗X),    g = e^{-2γt}

Sector weights: W₀ = 1/4, W₁ = g²/2, W₂ = g⁴/4.

The coherence is Ψ = (2g + g²)/3, which involves g = √(2W₁) linearly.
The derivative dΨ/dt = -4γg(1+g)/3 requires g, not W₁ = g²/2. The
square root is the obstruction: **L1-coherence sums absolute values of
Pauli coefficients, while sector weights sum their squares.**

|+⟩^{⊗2} crossing: g_cross = 0.6102, dwell prefactor = **1.7248**.
The Bell+ weight formula (2+4W₂)/(1+6W₂) applied to |+⟩^{⊗2} gives
1.7704, which differs by 2.6%. The discrepancy is the k = 1 contribution.

### The dividing line

States with **only even-weight Pauli content** (Bell+, GHZ, any state
that is a superposition of computational basis states differing in all
bits): the prefactor is a pure weight function. The face-swap algebra
and the cusp dynamics are algebraically identical.

States with **odd-weight content** (product states, W states, any state
with single-qubit coherence): the prefactor additionally requires the
Pauli coefficient magnitudes, not just sector weights. The connection is
structural but not a pure weight identity.

---

## 5. What the Prefactor Depends On

| Ingredient | Source | Weight-only? |
|-----------|--------|-------------|
| dC/dt     | -4γΣ(k·W_k) | Yes |
| C at crossing | Σ W_k | Yes |
| Ψ at crossing | L1-coherence / (d-1) | **No** (needs coefficient signs) |
| dΨ/dt     | d(L1)/dt / (d-1) | **No** (needs coefficient magnitudes) |
| CΨ = 1/4 condition | C·Ψ = 1/4 | Constrains Ψ given C |

The prefactor is fully determined by:
1. The sector weights W_k at the crossing (algebraic, weight-only)
2. The functional form Ψ(W_0, W_1, W_2, ...) relating L1-coherence to
   the weights (state-specific, depends on Pauli coefficient signs)

For Bell+, item 2 is trivial: Ψ = √(2W₂)/3. For general states, item 2
carries the sign-pattern information that sector weights do not encode.

---

## Verdict

**Partial success.** For Bell+ (and states with only even-weight Pauli
content), the dwell-time prefactor is a pure function of the light-face
sector weight W₂ at the crossing. The static face-swap algebra and the
dynamic cusp passage are algebraically identical: **prefactor = (2+4W₂)/(1+6W₂)**.

For general states, the prefactor requires the Pauli coefficient
magnitudes, not just sector weights. The missing information is the
sign pattern of the Pauli decomposition, which determines L1-coherence
(Ψ) but is not encoded in the sector weights W_k = (1/d)Σ|a_P|².

---

## Data Files

- [dwell_prefactor_from_weights.txt](../simulations/results/dwell_prefactor_from_weights.txt)
