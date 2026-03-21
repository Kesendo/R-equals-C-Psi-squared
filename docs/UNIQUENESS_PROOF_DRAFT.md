# Uniqueness of the 1/4 Boundary: Formal Theorem Draft

**Tier:** 1 (algebraic proof) + 2 (computational verification)
**Date:** March 21, 2026
**Status:** Draft. The algebra is complete. The physical justification
for the quadratic structure needs formalization.

---

## Theorem Statement

**Theorem (Uniqueness of the 1/4 Boundary).** Let R(R_in) = C(Ψ + R_in)²
be the self-referential purity map where C is the correlation bridge
(0 ≤ C ≤ 1) and Ψ is the normalized l1-coherence (0 ≤ Ψ ≤ 1). Then:

(i) The fixed-point equation R = C(Ψ + R)² has exactly two solutions
    when CΨ < 1/4 (one stable, one unstable), exactly one solution
    when CΨ = 1/4 (marginal), and no real solutions when CΨ > 1/4.

(ii) The value 1/4 is uniquely determined by the quadratic structure
     of the map: it is the discriminant zero of the quadratic
     CR² + (2CΨ - 1)R + CΨ² = 0.

(iii) No reparameterization of C and Ψ that preserves their
      operational definitions (purity-based and l1-based respectively)
      can change the boundary value.

---

## Proof

**Step 1.** Expand R = C(Ψ + R)² to get CR² + (2CΨ - 1)R + CΨ² = 0.

**Step 2.** Discriminant D = (2CΨ - 1)² - 4C²Ψ² = 4C²Ψ² - 4CΨ + 1 - 4C²Ψ² = 1 - 4CΨ.

**Step 3.** D = 0 iff CΨ = 1/4. D > 0 iff CΨ < 1/4. D < 0 iff CΨ > 1/4. QED for (i) and (ii).

**Step 4 (iii).** The factor 4 in the discriminant formula b² - 4ac is
a consequence of completing the square in the general quadratic ax² + bx + c = 0.
It is built into the definition of "quadratic equation." To change the 1/4
would require changing the power from 2 to something else.

**Step 5.** The power 2 comes from Purity = Tr(ρ²), which is inherently
degree 2 in the matrix elements of ρ. No physical redefinition changes this.

---

## Computational Verification

All standard Markovian quantum channels cross CΨ = 1/4 (not any other value):

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

Source: [proof_roadmap_close.py](../simulations/proof_roadmap_close.py)

---

## What Remains Open

The step from "purity is degree 2" to "the recursion must be quadratic"
requires showing that no higher-order corrections (Renyi entropies
S_alpha for alpha ≠ 2, or non-linear functions of ρ) modify the
effective fixed-point equation.

Specifically: if we replace Tr(ρ²) with Tr(ρ³) (Renyi-3), the
recursion becomes CUBIC, and the bifurcation boundary shifts.
The claim is that Tr(ρ²) is the physically correct choice because
it is the unique Renyi entropy that appears in the purity recursion
of the Lindblad equation. This needs to be made rigorous.

The connection to catastrophe theory (fold catastrophe with CΨ as
bifurcation parameter) should be made explicit.

---

## References

- [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer roadmap
- [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md): CΨ ↔ c mapping
- [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md): theta compass
- [IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md): hardware confirmation at 1.9%
