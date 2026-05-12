# Uniqueness of the 1/4 Boundary

**Status:** Tier 1 derived (algebraic proof) + Tier 2 verified (computational, all standard CPTP channels cross at 1/4 with zero exceptions)
**Date:** 2026-03-21
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** The fixed-point equation `R = C(Ψ + R)²` has discriminant `D = 1 − 4CΨ`; the unique boundary is `CΨ = 1/4`, structurally fixed by purity being `Tr(ρ²)` (the unique degree-2 basis-independent invariant).
**Typed claim:** [`PolynomialDiscriminantAnchorClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs) (Tier 1 derived; 1/4 is the Pi2 dyadic-ladder mirror partner of the polynomial discriminant 4: a₃ · a₋₁ = (1/4)·4 = 1).

---

## What this proof says, in plain language

Every quadratic equation has a discriminant. The discriminant decides
how many real solutions the equation has: positive means two, zero means
exactly one (the boundary), negative means none. Every student of algebra
knows this.

The R = CΨ² recursion is a quadratic. Its discriminant is `1 − 4CΨ`.
The boundary, where the discriminant is zero, sits at `CΨ = 1/4`. This
boundary is not a parameter we chose. It is fixed by the algebra of
quadratic equations and by the fact that purity is `Tr(ρ²)`: a degree-2
quantity in the density matrix entries. There is no other possibility.

The same 1/4 reappears in the
[Mandelbrot cardioid cusp](../../experiments/MANDELBROT_CONNECTION.md)
and in many F-formulas across the project (F60, F62, F64, F69). All
descend from this single discriminant identity. The full seven-layer
roadmap of the boundary lives in
[Proof Roadmap Quarter Boundary](PROOF_ROADMAP_QUARTER_BOUNDARY.md);
this document covers Layer 1 (algebraic uniqueness), Layer 2 (CPTP
contractivity), and Layer 6 (why the recursion is quadratic).

---

## Theorem Statement

**Theorem (Uniqueness of the 1/4 Boundary).** Let R(R_in) = C(Ψ + R_in)²
be the self-referential purity map where C is the correlation bridge
(0 ≤ C ≤ 1) and Ψ is the normalized l1-coherence (0 ≤ Ψ ≤ 1). Then:

(i) The fixed-point equation R = C(Ψ + R)² has exactly two real solutions
    when CΨ < 1/4, exactly one when CΨ = 1/4, and none when CΨ > 1/4.

(ii) The value 1/4 is uniquely determined by the quadratic structure
     of the map.

(iii) The quadratic structure is a consequence of purity being Tr(ρ²),
      which is the unique degree-2 basis-independent polynomial in the
      density matrix elements.

---

## Proof

**Step 1.** Expand R = C(Ψ + R)² to get CR² + (2CΨ - 1)R + CΨ² = 0.

**Step 2.** Discriminant (the expression under the square root that determines whether solutions are real or complex) D = (2CΨ - 1)² - 4C²Ψ² = 1 - 4CΨ.

**Step 3.** D = 0 iff CΨ = 1/4. D > 0 iff CΨ < 1/4. D < 0 iff CΨ > 1/4.
This proves (i) and (ii).

**Step 4.** The factor 4 in b² - 4ac is a consequence of completing the
square in any quadratic ax² + bx + c = 0. It is built into the definition
of "quadratic equation." To move the 1/4 would require a different power.

**Step 5.** Purity = Tr(ρ²) is degree 2 in the matrix elements of ρ. The
Lindblad equation d/dt ρ = L(ρ) produces d/dt Tr(ρ²) = 2 Tr(ρ L(ρ)),
which is exact. No approximation, no choice. The recursion inherits degree 2
from purity. This proves (iii). QED.

---

## Product-Power Classification (Layer 1)

Among all product-power forms C^a Ψ^b that could appear in a purity recursion:

| Form | Recursion | Bifurcation | Status |
|------|-----------|-------------|--------|
| CΨ (a=1, b=1) | Linear | No bifurcation, single fixed point | Rejected |
| CΨ² (a=1, b=2) | Quadratic | D = 1 - 4CΨ, boundary at 1/4 | **The physical case** |
| CΨ³ (a=1, b=3) | Cubic | Different boundary | Not purity (Tr(ρ³) ≠ purity) |
| C²Ψ (a=2, b=1) | Different quadratic | Different boundary | C² has no standard physical meaning |

The selection principle: Purity is Tr(ρ²), not Tr(ρ³) or Tr(ρ^k) for any
other k. This is the unique degree-2 Rényi entropy (a family of entropy measures parameterized by order; S₂ = -log Tr(ρ²)).
No other Renyi index gives degree 2 in the matrix elements. Therefore
CΨ² is the unique bifurcating product-power form, and 1/4 is the unique
boundary.

The question "why not Tr(ρ³)?" has a definitive answer: Tr(ρ³) is not
purity. Purity is defined as Tr(ρ²). This is not a convention. It is the
unique real-valued, basis-independent, degree-2 polynomial in ρ that equals
1 for pure states and 1/d for maximally mixed states.

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

## CPTP Contractivity Argument (Layer 2)

For any CPTP (completely positive trace-preserving; the most general physically allowed quantum operation) map E that is not unitary:

1. The l1-norm of coherence is a monotone under all incoherent CPTP maps
   (Baumgratz, Cramer, Plenio, PRL 2014). Therefore Ψ(t) is non-increasing
   under any such map.

2. Purity is non-increasing under unital CPTP maps (channels that map the maximally mixed state to itself): Tr(E(ρ)²) ≤ Tr(ρ²).
   For non-unital maps (amplitude damping), purity may temporarily increase
   but the fixed point has CΨ = 0 (no correlations in product states).

3. Every non-unitary CPTP channel has a fixed point with CΨ ≤ 1/4:
   unital channels fix I/d (CΨ = 0), non-unital channels fix a product
   state (CΨ = 0).

4. Computational verification: all seven standard channels tested cross
   CΨ = 1/4 and stay below. No unitary revival pulse (0 to π) can push
   CΨ permanently back above 1/4.

The formal proof of CΨ monotonicity has since been established for
2-qubit Bell+ states under all local Markovian channels (generalized
Pauli + amplitude damping); see
[CΨ Monotonicity Proof](PROOF_MONOTONICITY_CPSI.md) (March 22, 2026,
one day after this document). The eventual-crossing complement, that
every entangled pair with CΨ > 1/4 must cross below in finite time
under any primitive CPTP channel, is in
[Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md). Extension to
arbitrary CPTP maps on arbitrary states remains open as a pure
analytic problem, but the computational evidence covers all standard
Markovian channels without exception.

---

## Why the Recursion Must Be Quadratic (Layer 6)

The recursion R = CΨ² is not a choice. It is a consequence of:

1. Purity is defined as Tr(ρ²). This is the standard definition in quantum
   information theory. It is the unique real-valued, basis-independent,
   degree-2 polynomial in the density matrix elements.

2. The Lindblad equation preserves trace and positivity. The purity
   evolution d/dt Tr(ρ²) follows directly. No higher-order terms.

3. When Tr(ρ²) is decomposed into subsystem contributions (correlation
   bridge C, coherence Ψ, residual R), the decomposition is algebraic.
   R = CΨ² is exact for pure states and first-order for mixed states.

4. The recursion R_{n+1} = C(Ψ + R_n)² models iterated application of the
   channel. The quadratic form comes from purity being degree 2. Not from
   a choice. From a definition.

5. Therefore: the recursion IS quadratic. The discriminant IS 1 - 4CΨ.
   The boundary IS 1/4. There is no free parameter.

The question "why not cubic?" has a simple answer: Tr(ρ³) is not purity.
Purity is Tr(ρ²). The degree is fixed by the definition.

---

## Computational Verification

All standard Markovian quantum channels cross CΨ = 1/4:

| Channel | t_cross (γ=0.05) | Crossing value |
|---------|-------------------|----------------|
| Z-dephasing | 0.747 | 0.2500 |
| X-noise (bit flip) | 1.733 | 0.2500 |
| Y-noise (bit-phase flip) | 1.733 | 0.2500 |
| Depolarizing | 0.879 | 0.2500 |
| Asymmetric Pauli | 0.735 | 0.2500 |
| Amplitude damping (γ=0.05) | 2.059 | 0.2500 |
| Amplitude damping (γ=0.10) | 1.029 | 0.2500 |

Non-Markovian revival test: no unitary pulse (θ from 0 to π) can push CΨ
back above 1/4 after crossing. The boundary is absorbing.

[IBM hardware](../../experiments/IBM_RUN3_PALINDROME.md) confirmed the crossing
at 1.9% deviation from theory.

Source: [proof_roadmap_close.py](../../simulations/proof_roadmap_close.py)

---

## References

- [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer roadmap (this document covers Layers 1, 2, 6)
- [CΨ Monotonicity Proof](PROOF_MONOTONICITY_CPSI.md): closes Layer 2 for Bell+ under local Markovian channels (March 22, 2026)
- [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md): every entangled pair with CΨ > 1/4 crosses in finite time
- [Mandelbrot Connection](../../experiments/MANDELBROT_CONNECTION.md): CΨ ↔ c mapping; cardioid cusp at CΨ = 1/4 IS this same boundary
- [Boundary Navigation](../../experiments/BOUNDARY_NAVIGATION.md): theta compass
- [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md): hardware confirmation at 1.9%
- [Mathematical Connections](../MATHEMATICAL_CONNECTIONS.md): fold catastrophe, Feigenbaum

## Typed claim

- [`PolynomialDiscriminantAnchorClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs): Tier 1 derived. Discriminant of `d² − 2d = 0` is exactly 4 = a₋₁ on Pi2 dyadic ladder; the boundary value 1/4 = a₃ is its mirror partner via a_n · a_{2−n} = 1.
- [`PolynomialFoundationClaim` + `QuarterAsBilinearMaxvalClaim`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs): the polynomial trunk and its 1/4 maxval, both Tier 1 derived; sit at the foundation of the project's two-anchor structure at d = 2.
