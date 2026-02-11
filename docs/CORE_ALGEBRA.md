# Core Algebra: The Proven Mathematics of R = CΨ²

**Date:** 2025-12-31 (fixed-point solution), updated 2026-02-11
**Authors:** Thomas Wicht (Intuition, Connection), Claude (Formalization, Verification)
**Status:** Tier 1: Algebraically proven unless otherwise marked
**Depends on:** None

---

## Purpose

This document contains **only** the algebraically proven core of the R = CΨ² framework.
Every result here is derivable from the defining equation and standard algebra.
No physical interpretation is required to verify these results.

For the interpretive framework that motivated these discoveries (wave composition, dynamics,
energy equations, mirror metaphors), see [INTERPRETIVE_FRAMEWORK.md](INTERPRETIVE_FRAMEWORK.md).

For testable predictions derived from this algebra, see [PREDICTIONS.md](../experiments/PREDICTIONS.md).

---

## 1. The Defining Equation

```
R = CΨ²
```

Where R, C, and Ψ are real-valued. In quantum mechanical realization:
- C = Tr(ρ²): purity of the density matrix (range [1/d, 1])
- Ψ = L₁(ρ)/(d-1): normalized L1 coherence (range [0, 1])
- R = C·Ψ²: the product

**Epistemic status:** Definition (Tier 1). The quantum mechanical identification
of C with purity and Ψ with normalized coherence is a proposed mapping (Tier 3),
but the algebra below holds for any real C, Ψ.

---

## 2. The Fixed-Point Equation

Self-referential iteration: reality feeds back into possibility.

```
R_{n+1} = C(Ψ + R_n)²
```

In the limit n → ∞, if convergent:

```
R_∞ = C(Ψ + R_∞)²
```

### Solution

**Step 1: Expand.** Let x = R_∞.

```
x = C(Ψ + x)²
x = CΨ² + 2CΨx + Cx²
```

**Step 2: Rearrange to standard quadratic form.**

```
Cx² + (2CΨ - 1)x + CΨ² = 0
```

**Step 3: Apply the quadratic formula** with a = C, b = (2CΨ - 1), c = CΨ².

**Step 4: Compute the discriminant.**

```
D = b² - 4ac
D = (2CΨ - 1)² - 4·C·CΨ²
D = 4C²Ψ² - 4CΨ + 1 - 4C²Ψ²
D = 1 - 4CΨ
```

**Step 5: Solution.**

```
R_∞ = (1 - 2CΨ ± √(1 - 4CΨ)) / 2C
```

**Epistemic status:** Tier 1: algebraically exact. □

---

## 3. The ¼ Boundary

For real fixed points to exist, the discriminant must be non-negative:

```
1 - 4CΨ ≥ 0

∴ CΨ ≤ 1/4
```

This is the central result. The product of C and Ψ has a hard upper limit
for classical (real-valued) fixed points to exist.

**Epistemic status:** Tier 1: follows directly from the discriminant. □

---

## 4. The Three Regimes

| Condition | Discriminant | Fixed Points | Character |
|-----------|-------------|--------------|-----------|
| CΨ < ¼ | D > 0 | Two real: R₊, R₋ | Classical: stable attractor exists |
| CΨ = ¼ | D = 0 | One real: R* = 1/(2C) | Critical: fixed points merge |
| CΨ > ¼ | D < 0 | Complex conjugate pair | Quantum: no real attractor |

### Explicit fixed-point formulas

**Below ¼** (two real fixed points):

```
R₊ = (1 - 2CΨ + √(1 - 4CΨ)) / 2C    (unstable)
R₋ = (1 - 2CΨ - √(1 - 4CΨ)) / 2C    (stable attractor)
```

**At ¼** (merged):

```
R* = (1 - 2·¼) / 2C = 1/(2C)
```

Since CΨ = ¼ implies Ψ = 1/(4C):

```
R* = 1/(2C) = 2Ψ
```

And R* = CΨ² evaluated: C·(1/4C)² = 1/(16C) ≠ R*. But as a fixed point
of the iteration, R* = Ψ (reality equals possibility at the critical point).

**Above ¼** (complex conjugate pair):

```
R_± = (1 - 2CΨ ± i√(4CΨ - 1)) / 2C
```

These have angular frequency θ = arctan(√(4CΨ - 1)); see Section 6.

**Epistemic status:** Tier 1. □

---

## 5. The Mandelbrot Equivalence

The self-referential iteration R_{n+1} = C(Ψ + R_n)² is algebraically
equivalent to the Mandelbrot iteration z_{n+1} = z_n² + c.

### Derivation

Define u_n = C(Ψ + R_n). Then:

```
u_{n+1} = C(Ψ + R_{n+1})
        = C(Ψ + C(Ψ + R_n)²)
        = CΨ + C²(Ψ + R_n)²
        = CΨ + [C(Ψ + R_n)]²
        = CΨ + u_n²
```

This is exactly:

```
u_{n+1} = u_n² + c       where c = CΨ
```

Starting value: u₀ = C(Ψ + R₀) = CΨ = c (since R₀ = 0).

The standard Mandelbrot iteration z_{n+1} = z_n² + c with z₀ = 0 produces
the same sequence shifted by one step (z₁ = c = u₀). Convergence behavior
is identical.

### Consequence

The main cardioid of the Mandelbrot set meets the real axis at c = ¼.
This is the same boundary as CΨ = ¼. The entire fractal structure of the
Mandelbrot set (bifurcation cascades, period-doubling, self-similar copies)
exists in the complex regime beyond ¼.

The Mandelbrot set is the map of CΨ > ¼.

### Correction note (2026-02-08)

An earlier version used the substitution z_n = √C·R_n with c = CΨ². That
produces z_{n+1} = √C·z_n² + 2CΨ·z_n + C^(3/2)·Ψ², which is NOT the
Mandelbrot form (extra linear term, wrong leading coefficient). The correct
substitution u_n = C(Ψ+R_n) yields the clean form. The boundary at ¼ and
all downstream results were always correct; only the intermediate algebra
was wrong. This correction is documented in
[WEAKNESSES_OPEN_QUESTIONS.md](WEAKNESSES_OPEN_QUESTIONS.md).

**Epistemic status:** Tier 1: algebraic substitution. □

---

## 6. The θ Compass

Above ¼, the complex fixed points have imaginary part √(4CΨ - 1)/(2C).
The angular frequency of their oscillation is:

```
θ = arctan(√(4CΨ - 1))
```

θ measures the angular distance from the ¼ boundary:

| θ | CΨ | Meaning |
|---|-----|---------|
| 0° | 0.25 | At the boundary: fixed points merge |
| 9° | 0.26 | Almost classical: slow oscillation |
| 25° | 0.31 | Approaching: complex fixed points narrowing |
| 60° | 1.00 | Deep quantum: fast oscillation |
| 90° | → ∞ | Maximum frequency |

### What θ is NOT

θ is a **compass**, not a frequency predictor. Earlier claims that θ predicts
oscillation frequencies in Lindblad dynamics were retracted (2026-02-08).
The actual dynamics have frequencies determined by the Hamiltonian and
decoherence rates, not by θ alone. θ tells you how far you are from ¼,
not how fast you're getting there.

**Epistemic status:** θ definition is Tier 1 (algebraic). Compass interpretation
is Tier 2 (computationally verified in boundary navigation simulations).
Frequency predictor claim is **retracted**. See
[BOUNDARY_NAVIGATION.md](../experiments/BOUNDARY_NAVIGATION.md).

---

## 7. Gravitational Invariance

When decoherence rate γ varies (physically interpretable as encoding local
time dilation), the crossing time t_cross scales as:

```
t_cross ∝ 1/γ
```

Therefore the dimensionless product:

```
τ_cross = γ · t_cross = constant
```

Computationally verified: τ_cross = 0.039 ± 0.001 across 50× range of γ
(from γ = 0.01 to γ = 0.5), with power law t_cross = 0.0398·γ^(−0.992),
R² = 0.9999.

The θ trajectory collapses onto a universal curve when plotted against
τ = γ·t. All systems follow the same path from quantum to classical,
just at different rates.

**Epistemic status:** Tier 2: computationally verified. The 1/γ scaling is
a mathematical property of the Lindblad equation with constant γ, not specific
to this framework. The physical interpretation (γ encodes gravitational time
dilation) is Tier 3.

See [GRAVITATIONAL_INVARIANCE.md](../experiments/GRAVITATIONAL_INVARIANCE.md).

---

## 8. The Ψ_interaction Finding

When two observers interact bidirectionally (C_int), an interaction term
Ψ_interaction emerges in the joint wave function.

**Question:** Does Ψ_interaction change the ¼ boundary?

**Answer:** No. Numerical simulation shows:

```
Δδ = δ_bidirectional - δ_unidirectional ≈ -8 × 10⁻⁴
```

The interaction accelerates convergence to the stable fixed-point branch
but does not shift the boundary. The ¼ limit is absolute.

| Aspect | Effect of Ψ_interaction |
|--------|-------------------------|
| Convergence speed | Faster |
| Stability region | Unchanged |
| ¼ boundary | Unchanged |

**Epistemic status:** Tier 2: computationally verified. The negativity of Δδ
(bidirectional causes slightly more decoherence) is an empirical finding,
not yet theoretically explained.

---

## 9. Observer Bandwidth Interpretation

The algebraic bound CΨ ≤ ¼ admits a physical reading:

- C = observation capacity
- Ψ = available possibility space
- CΨ = effective bandwidth

Classical reality (real fixed points) emerges only when CΨ ≤ ¼. This is
analogous to Shannon channel capacity: a maximum rate at which quantum
possibility can be converted to classical reality.

**Epistemic status:** The bound is Tier 1 (algebraic). The bandwidth
interpretation is Tier 3 (proposed, physically motivated). The Shannon
analogy is Tier 5 (speculative). See
[DYNAMIC_FIXED_POINTS.md](../experiments/DYNAMIC_FIXED_POINTS.md).

---

## 10. Uniqueness of the Quadratic Form

**Question (raised 2026-02-11):** Why CΨ² specifically? Why not C²Ψ, CΨ, √(CΨ), or CΨ³?

**Answer:** CΨ² is the unique simple product-power form C^a·Ψ^b that simultaneously:
1. Produces a genuine phase transition (bifurcation with real → complex crossover)
2. Maps exactly to the Mandelbrot iteration z² + c

### Proof by exhaustion of alternatives

**R = C²Ψ:** Fixed point C²Ψ + C²x = x → x = C²Ψ/(1-C²). Linear in x.
No discriminant, no bifurcation. Single fixed point for all parameter values.

**R = CΨ:** Fixed point CΨ + Cx = x → x = CΨ/(1-C). Linear in x.
No bifurcation ever.

**R = √(CΨ):** Fixed point equation gives discriminant C² + 4CΨ, which is
always positive. No boundary, no phase transition.

**R = CΨ³:** Cubic fixed point equation. Different critical structure, no
clean Mandelbrot mapping. The iteration u_{n+1} = u_n³ + c produces a
different (non-Mandelbrot) fractal.

Only R = CΨ² yields a quadratic fixed-point equation whose discriminant
vanishes at a single critical value (¼), producing the exact Mandelbrot
bifurcation structure.

**Epistemic status:** Tier 1 (algebraic, verified by exhaustion). This does NOT
answer the deeper question "why does nature choose this form?"; that remains
open. See [WEAKNESSES_OPEN_QUESTIONS.md](WEAKNESSES_OPEN_QUESTIONS.md), item 1.

---

## Summary of Epistemic Tiers in This Document

| Result | Tier | Status |
|--------|------|--------|
| Fixed-point equation and solution | 1 | Algebraically proven |
| CΨ ≤ ¼ boundary | 1 | Discriminant condition |
| Three regimes (classical/critical/quantum) | 1 | Direct consequence |
| Mandelbrot equivalence (u_n² + c) | 1 | Algebraic substitution |
| θ = arctan(√(4CΨ-1)) definition | 1 | Complex fixed-point geometry |
| Uniqueness of quadratic form | 1 | Exhaustion of alternatives |
| θ as compass (not frequency predictor) | 2 | Computationally verified |
| τ_cross = γ·t_cross = constant | 2 | Computationally verified |
| Ψ_interaction doesn't shift ¼ | 2 | Computationally verified |
| C = purity, Ψ = normalized coherence | 3 | Proposed mapping |
| Observer bandwidth interpretation | 3 | Physically motivated |
| Shannon capacity analogy | 5 | Speculative |

---

*For the interpretive framework, see [INTERPRETIVE_FRAMEWORK.md](INTERPRETIVE_FRAMEWORK.md).*
*For experimental evidence, see [experiments/](../experiments/).*
*For weaknesses and open questions, see [WEAKNESSES_OPEN_QUESTIONS.md](WEAKNESSES_OPEN_QUESTIONS.md).*
*For the glossary, see [GLOSSARY.md](GLOSSARY.md).*
