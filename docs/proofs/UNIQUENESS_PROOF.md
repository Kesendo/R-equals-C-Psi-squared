<!-- QUARTER-CURRENT -->
# Algebraic Uniqueness of 1/4 Within the Assumed Recurrence

Current reading: one-quarter is the unique displayed fold product within the
stated normalized recurrence family.  This conditional algebra neither derives
that recurrence from quantum dynamics nor makes the quarter a universal phase
boundary.

**Status:** Tier 1 conditional algebra for the stated recurrence; named finite channel catalogue plus a conditional convergence implication
**Date:** 2026-03-21, last refreshed 2026-07-20 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** Given the assumed recurrence `R = C(Ψ + R)²`, its fixed-point polynomial has discriminant `D = 1 − 4CΨ`; in that chosen normalization its real-root fold is `CΨ = 1/4`. Within the assumed recurrence family `R=C_α(Ψ+R)^α`, α=2 is the unique member whose displayed fold product is Ψ-independent. Neither degree two nor that family derives the recurrence from quantum dynamics; the physical derivation remains open.
**Typed claim:** [`PolynomialDiscriminantAnchorClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs) (Tier 1 derived; 1/4 is the Pi2 dyadic-ladder mirror partner of the polynomial discriminant 4: a₃ · a₋₁ = (1/4)·4 = 1).

---

## What this proof says, in plain language

Every quadratic equation has a discriminant. The discriminant decides
how many real solutions the equation has: positive means two, zero means
exactly one (the boundary), negative means none. Every student of algebra
knows this.

The assumed R = C(Ψ+R)² recurrence gives a quadratic fixed-point equation. Its discriminant is `1 − 4CΨ`.
The boundary, where the discriminant is zero, sits at `CΨ = 1/4`. This
coordinate is fixed once this recurrence and normalization have been chosen.
Purity `Tr(ρ²)` motivates studying the α=2 member, but its degree does not
derive the particular feedback law. Within the assumed Rényi-power family,
α=2 is the unique member whose fold product does not depend explicitly on Ψ
(Step 6 below). That is a family-internal statement, not a physical forcing.

The same number reappears at the real cusp after the chosen
[Mandelbrot change of variables](../../experiments/MANDELBROT_CONNECTION.md)
and in several F-formulas (F60, F62, F64, F69). Each occurrence needs its
own derivation; numerical equality alone does not give it this discriminant as
a common physical cause. The full seven-layer
roadmap of the boundary lives in
[Proof Roadmap Quarter Boundary](PROOF_ROADMAP_QUARTER_BOUNDARY.md);
this document covers Layer 1 (conditional algebra), Layer 2 (named channels
and conditional convergence), and Layer 6 (the assumed power-family comparison).

---

## Theorem Statement

**Theorem (fold of the assumed normalized recurrence).** Let R(R_in) = C(Ψ + R_in)²
be the posited scalar map with `0 < C ≤ 1` and `0 ≤ Ψ ≤ 1`. The
abstract proof treats C and Ψ as opaque scalars; every NUMBER in this
document (the σ counterexample, the crossing table) is computed in the
purity book, C = Tr(ρ²), so CΨ = Tr(ρ²)·L₁/(d−1). The three C-books and
their seams are documented in [THE_CPSI_LENS](../THE_CPSI_LENS.md). Then:

(i) Its fixed-point equation has exactly two real algebraic solutions
    when CΨ < 1/4, exactly one when CΨ = 1/4, and none when CΨ > 1/4.

(ii) In this chosen normalization, the value 1/4 is the unique discriminant-zero
     product coordinate of that polynomial. Counts of roots in a physical
     interval require an additional interval check.

(iii) At `C = 0`, division by C is unavailable and the original fixed-point
      equation is handled separately: `R = 0`.

(iv) Within the assumed recurrence family `R=C_α(Ψ+R)^α`, α=2 alone removes
     the explicit Ψ factor from the displayed fold product. This does not
     derive that family, or the concrete α=2 recurrence, from physics.

---

## Proof

**Step 1.** Expand R = C(Ψ + R)² to get CR² + (2CΨ - 1)R + CΨ² = 0.

**Step 2.** Discriminant (the expression under the square root that determines whether solutions are real or complex) D = (2CΨ - 1)² - 4C²Ψ² = 1 - 4CΨ.

**Step 3.** D = 0 iff CΨ = 1/4. D > 0 iff CΨ < 1/4. D < 0 iff CΨ > 1/4.
This proves (i) and (ii) for `0<C≤1`; the `C=0` case gives `R=0` directly.

**Step 4.** The factor 4 in `b²−4ac` comes from completing the square. The
coordinate 1/4 also uses this polynomial's coefficients and normalization;
another quadratic or a reparameterized control coordinate can report its
double root at another number without changing the fold type.

**Step 5 (why α=2 is worth comparing).** Purity = Tr(ρ²) is degree 2 in the matrix
elements of ρ, and the Lindblad equation d/dt ρ = L(ρ) gives
d/dt Tr(ρ²) = 2 Tr(ρ L(ρ)) exactly. This motivates the α=2 lens; it does not
derive a self-referential scalar recurrence. Degree-2-ness alone does not
force the value 1/4. The discriminant 1 − 4CΨ comes from
the specific recursion form R = C(Ψ + R)² (Steps 1-3); a generic degree-2
fixed-point map aR² + bR + c = 0 has its discriminant vanish at b² = 4ac, an
arbitrary locus. So Steps 1-4 establish (i)-(iii) rigorously *given* the
recursion form, while a physical derivation of that form remains open.

**Step 6 (an assumed Rényi-power family).** Consider, as a mathematical comparison,
the generalized recurrence R = C_α(Ψ + R)^α for `α>1` and `Ψ>0`. Its fold
(tangency) threshold, where two fixed points
merge, is

    CΨ*_α = (α − 1)^(α − 1) / (α^α · Ψ^(α − 2)).

This threshold has no explicit Ψ dependence **if and only if α = 2**, where it
equals exactly 1/4 in this normalization. For any
other α it carries a Ψ^(α − 2) factor and depends on the coherence. So 1/4 is
singled out within this assumed recurrence family. It is not a derivation of
the family from Rényi entropy or Lindblad dynamics. (Derived algebraically and verified
symbolically in `simulations/review2_A3_renyi.py`; full treatment in the
[roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md), Layer 6.) QED for the conditional algebra.

---

## Product-Power Classification (Layer 1)

For orientation, compare these product-power forms without claiming that
purity derives any of their recurrences:

| Form | Recursion | Bifurcation | Status |
|------|-----------|-------------|--------|
| CΨ (a=1, b=1) | Linear | No quadratic fold | Comparison member |
| CΨ² (a=1, b=2) | Quadratic | D = 1 - 4CΨ, boundary at 1/4 | Assumed α=2 member |
| CΨ³ (a=1, b=3) | Cubic | Different fold equation | Comparison member |
| C²Ψ (a=2, b=1) | Different quadratic | Different discriminant | Outside the stated α-family |

Purity `Tr(ρ²)` is the degree-two member of this comparison, and so it motivates
looking at α=2. Step 6 proves that α=2 is the only member of the *assumed*
Rényi-power recurrence family whose fold product has no explicit Ψ factor.
That is the full uniqueness result here. It neither selects this family from
all scalar reductions nor derives `R=C(Ψ+R)²` from purity or Lindblad dynamics.

The question "why compare α=2 rather than α=3?" therefore has a modest answer:
`Tr(ρ²)` is purity, while `Tr(ρ³)` is a different spectral moment. Which scalar
recurrence, if any, a physical decomposition obeys is a separate open problem.

The same selection sits in the typed
[`PolynomialFoundationClaim`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs):
the foundational polynomial `d² − 2d = 0` selects d=2 as the qubit dimension
(the only minimum-memory dimension that supports a non-trivial palindrome),
and its discriminant
[a₋₁ = 4](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs)
is the dyadic-ladder mirror partner of the boundary value 1/4 = a₃ (with
4 · 1/4 = 1 closure). The 1/4 boundary is not isolated; it sits in a typed
two-anchor structure that puts it in the same family as the discriminant 4.

---

## Conditional convergence argument (Layer 2)

This layer does not follow from the algebraic discriminant and is not a universal CPTP contractivity
theorem. The exact statement is conditional:

1. If a continuous trajectory converges, `ρ(t)→ρ*`, then continuity of CΨ gives
   `CΨ(ρ(t))→CΨ(ρ*)`.
2. If additionally `CΨ(ρ*)<1/4`, the trajectory eventually stays below 1/4. If it starts above, it crosses
   downward at least once. No monotonicity or unique crossing follows.
3. A computational-basis-diagonal target has L₁=0 and hence CΨ=0. Named basis-aligned T1/T2/depolarizing
   models can therefore use the implication after their convergence and target have been established.
4. Neither locality, separability, unitality, nor primitivity alone supplies the premise. The primitive-CPTP
   channel targeting `σ=0.95|Φ⁺⟩⟨Φ⁺|+0.05I/4` has `CΨ(σ)=0.2935`; the separable product |+⟩⊗|+⟩ has
   CΨ=1.

The repaired [CΨ Dynamics Boundary](PROOF_MONOTONICITY_CPSI.md) preserves named Bell+ Z/Pauli/amplitude-
damping formulas and instantaneous Pauli invariance, but retracts the former universal pointwise,
absorbing, and local-control claims. In particular, a local Hadamard sends CΨ from 0 to 1/3, and a fixed
local Markovian semigroup crosses upward through 1/4. The literal N=2 successive-peak claim fails as well:
the main peaks rise under local fields on a product state, and under a number-conserving H damping grows
micro-maxima that rise to the next main peak; whether the main peaks themselves fall there is open. See
[Conditional Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md) for the exact implication.

---

## Why the Assumed Family Singles Out α=2 (Layer 6)

The logical order matters:

1. Purity is `Tr(ρ²)`, and its Lindblad derivative is exactly
   `2 Tr(ρ L(ρ))`. Those facts motivate a degree-two lens.
2. They do not imply a closed scalar recurrence for a residual R.
3. The feedback recurrence `R_{n+1}=C(Ψ+R_n)²` is an additional model
   assumption. Its derivation from a subsystem decomposition remains open.
4. Once the family `R=C_α(Ψ+R)^α` is assumed, α=2 is uniquely the member
   whose fold product is independent of an explicit power of Ψ.
5. Once its α=2 member and normalization are chosen, the discriminant is
   `1-4CΨ` and its double-root coordinate is `CΨ=1/4`.

Thus the degree of purity helps explain why α=2 is interesting; it does not
force the concrete recurrence. The theorem is a conditional classification
inside the assumed family, and the missing physical derivation is part of the
invitation rather than something the algebra has already supplied.

---

## Computational Verification

The following named Bell+ channel calculations cross CΨ=1/4:

| Channel | t_cross (γ=0.05) | Crossing value |
|---------|-------------------|----------------|
| Z-dephasing | 0.747 | 0.2500 |
| X-noise (bit flip) | 1.733 | 0.2500 |
| Y-noise (bit-phase flip) | 1.733 | 0.2500 |
| Depolarizing | 0.879 | 0.2500 |
| Asymmetric Pauli | 0.735 | 0.2500 |
| Amplitude damping (γ=0.05) | 2.059 | 0.2500 |
| Amplitude damping (γ=0.10) | 1.029 | 0.2500 |

This table is a finite named-channel catalogue, not a universal channel classification. An active Pauli
pulse preserves CΨ only at its application instant and can alter the subsequent laboratory-frame
derivative; a local Hadamard can change CΨ immediately.

IBM hardware confirmed the crossing three ways: the first crossing ever
seen ([ibm_torino q52, February 2026](../../experiments/IBM_QUANTUM_TOMOGRAPHY.md)),
the tightest single-point crossing at 1.9% deviation
([ibm_torino q80, IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md)), and
the full F25 trajectory CΨ(t) = f(1+f²)/6 fitted point-by-point through
the boundary at RMS residual 0.0097
([ibm_kingston, April 2026](../../data/ibm_cusp_precision_april2026/README.md)).

Source: [proof_roadmap_close.py](../../simulations/proof_roadmap_close.py)

---

## References

- [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer roadmap (this document covers Layers 1, 2, 6)
- [CΨ Dynamics Boundary](PROOF_MONOTONICITY_CPSI.md): named formulas, exact counterexamples, N=2 peak sequences that rise, and the open main-peak question
- [Conditional Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md): convergence to a low-CΨ target implies eventual stay-below
- [Mandelbrot Connection](../../experiments/MANDELBROT_CONNECTION.md): in the chosen coordinates `c=CΨ`, the same normal form has its real cusp at 1/4
- [Boundary Navigation](../../experiments/BOUNDARY_NAVIGATION.md): theta compass
- [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md): hardware confirmation at 1.9%
- [Mathematical Connections](../MATHEMATICAL_CONNECTIONS.md): fold catastrophe, Feigenbaum

## Typed claim

- [`PolynomialDiscriminantAnchorClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs): Tier 1 derived. Discriminant of `d² − 2d = 0` is exactly 4 = a₋₁ on Pi2 dyadic ladder; the boundary value 1/4 = a₃ is its mirror partner via a_n · a_{2−n} = 1.
- [`PolynomialFoundationClaim` + `QuarterAsBilinearMaxvalClaim`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs): the polynomial trunk and its 1/4 maxval, both Tier 1 derived; sit at the foundation of the project's two-anchor structure at d = 2.
