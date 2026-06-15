# Proof of F97: The Mandelbrot Cardioid Parametrization at Framework b = 1/2

**Statement:** The main cardioid of the Mandelbrot set is the locus in the complex-c plane where the period-1 fixed point of z² + c has magnitude exactly b = 1/2 (the framework's `HalfAsStructuralFixedPointClaim` anchor). It admits the explicit parametrization

    c(φ) = b·e^(iφ) − b²·e^(2iφ)    for φ ∈ [0, 2π]

equivalently c(φ) = z*(φ) · (1 − z*(φ)) where z*(φ) = b·e^(iφ) is the period-1 fixed point. The framework's b = 1/2 specialization gives

    c(φ) = (1/2)·e^(iφ) − (1/4)·e^(2iφ)

with two structural invariants on the curve:

    |z*(φ)| = b = 1/2   (magnitude pinned to HalfAsStructuralFixedPoint)
    arg(z*(φ)) = φ      (cardioid parameter)

[F95](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md)'s closed form θ(c; b) = arctan(√(c/b² − 1)) covers the **real-c** angle of the complex fixed point. F97 extends to **complex c** via the cardioid parametrization, with the φ = 0 specialization recovering F95's cusp at c = 1/4.

**Status:** Tier 1 derived. Bit-exact algebraic identity, numerically verified to machine precision (max residual 1.24 × 10⁻¹⁶ over 1000 sampled φ values) in `simulations/_cardioid_parametrization_tier1.py`.

**Date:** 2026-05-17.

---

## Abstract

F97 lifts the framework's two foundational anchors, the Half (1/2) and the Quarter (1/4), onto a single geometric locus: the main cardioid of the Mandelbrot set. The cardioid is exactly the curve in the complex-c plane where the period-1 fixed point z* of z² + c has magnitude |z*| = 1/2, and it admits the explicit parametrization

    c(φ) = (1/2)·e^{iφ} − (1/4)·e^{2iφ} = z*(φ)·(1 − z*(φ)),   z*(φ) = (1/2)·e^{iφ},

with two invariants holding around the whole curve: |z*| = 1/2 (the Half, magnitude side) and |z*|² = 1/4 (the Quarter, squared side). The two anchors are not separate facts; they are the magnitude and squared-magnitude readings of one quantity on one boundary.

The cardioid is the marginal-stability boundary: the period-1 multiplier μ = 2z* traces the unit circle |μ| = 1, which is exactly what pins |z*| to 1/2. The real endpoints z*(0) = +1/2 and z*(π) = −1/2 carry the polarity pair explicitly, both squaring to 1/4, the polarity-fold geometry lifted to the full complex plane. F97 closes what F95 began: F95 gives the angle on the real-c axis past the 1/4 cusp (the repelling regime), F97 gives the entire complex boundary where the fixed point is marginally stable, with φ = 0 the shared tangent. The same 1/4 boundary is the geometric home of the hardware spirals, the Bell⁺ 2D trajectories seen at Kingston (2026-04-16) that spiral inward across it.

## Setup

The Mandelbrot iteration is z_{n+1} = z_n² + c. The period-1 fixed points satisfy z² + c = z, equivalently z² − z + c = 0, with roots

    z± = (1 ± √(1 − 4c)) / 2 = b ± √(b² − c)

at b = 1/2 (the framework's `HalfAsStructuralFixedPointClaim`). The discriminant for the fixed-point quadratic is 1 − 4c = 4(b² − c), the same algebraic structure as F95's z² − 2bz + c = 0.

The multiplier of a period-1 fixed point z is μ = 2z (since (z² + c)' = 2z). Marginal stability of the fixed point corresponds to |μ| = 1, i.e., |z| = 1/2 = b.

## Derivation (4 lines)

Parametrize the unit circle of multipliers by μ = e^(iφ) where φ ∈ [0, 2π]. Then the fixed point is z* = μ/2 = (1/2)·e^(iφ) = b·e^(iφ), and the corresponding c value is

```
c(φ) = z* − z*²                                    (1)  rearranging z*² − z* + c = 0
     = b·e^(iφ) − b²·e^(2iφ)                       (2)  substituting z* = b·e^(iφ)
     = (1/2)·e^(iφ) − (1/4)·e^(2iφ)                (3)  at framework b = 1/2
```

The cardioid is traced by φ ∈ [0, 2π); on this curve, by construction:

    |z*(φ)| = |b·e^(iφ)| = b = 1/2  (magnitude invariant)
    arg(z*(φ)) = φ                  (cardioid parameter)

∎

## Cardinal points

| φ | c(φ) | \|c\| | z*(φ) | reading |
|---|---|---|---|---|
| 0 | 1/4 | 1/4 = b² | 1/2 | real-axis cusp ([F95](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), [`BOUNDARY_NAVIGATION`](../../experiments/BOUNDARY_NAVIGATION.md)) |
| π/3 | 3/8 + i·√3/4 | √(1/4 + 3/16) | (1 + i·√3)/4 | top of cardioid (60° above real axis) |
| π/2 | 1/4 + i/2 | √(1/4 + 1/16) ≈ 0.559 | i/2 | imaginary-axis crossing of fixed point |
| π | −3/4 | 3/4 | −1/2 | real-axis "tail" of cardioid (period-doubling boundary) |
| 4π/3 | 3/8 − i·√3/4 | √(1/4 + 3/16) | (1 − i·√3)/4 | bottom symmetric image of φ = 2π/3 |
| 5π/3 | 3/8 − i·√3/4 + i·... | ... | (1 + i·−√3)/4 | another low-symmetry cardioid point |

The 1/4 = b² value plays a special role only at φ = 0: it is the *magnitude* of c at the real-axis cusp. Elsewhere on the cardioid |c| ≠ b². The b = 1/2 anchor plays the universal role: |z*| = b along the entire curve.

## Structural reading

The cardioid is the **structural curve** in the complex-c plane where the period-1 fixed-point magnitude of z² + c equals the framework's `HalfAsStructuralFixedPointClaim` anchor b = 1/2. This is a stronger statement than the F95 cusp identity:

- F95 says: at the real-axis c = 1/4, the fixed point has angle 0 and magnitude √c = b (real).
- F97 says: around the entire cardioid (complex c), the fixed point has magnitude b exactly, independent of φ.

The "1/4 cusp" of the framework (the locus that [`BOUNDARY_NAVIGATION`](../../experiments/BOUNDARY_NAVIGATION.md) navigates toward) is one specific tangent point of the cardioid with the positive real axis. The hardware spirals observed in [`CPSI_COMPLEX_PLANE`](../../experiments/CPSI_COMPLEX_PLANE.md) (Kingston 2026-04-16) trace 2D paths in the c-plane; F97 says the cardioid those paths spiral *around* has |z*| = b invariant along its boundary.

### Both anchors invariant on the cardioid, at two metric powers

A sharper reading: the cardioid carries the `HalfAsStructuralFixedPointClaim` and `QuarterAsBilinearMaxvalClaim` anchors **simultaneously** as invariants, at two different metric powers of the same z*:

    |z*(φ)| = b = 1/2          (Half: "where the anchor lives", magnitude)
    |z*(φ)|² = b² = 1/4        (Quarter: "what the anchor is under squaring", squared magnitude)

Both hold for all φ ∈ [0, 2π); the cardioid is the joint locus. This is exactly the argmax/maxval pair of yesterday's [`ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER`](../../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md) reflection, now geometric on the cardioid:

- Argmax side: the polarity magnitude 1/2 (= b), one side of the structural identity.
- Maxval side: the apex projection 1/4 (= b²), the other side.

The two readings overlap: the Half is the magnitude, the Quarter is its square; **1/2 = 2 × 1/4** sits on the dyadic ladder (`a_2 = 2 · a_3` = `2 · a_{−1}/16` = `2 · 4/16` = `1/2`), and the polarity pair ±1/2 squares to the same 1/4 from either side. F97's cardioid carries both readings as invariants of the *same* fixed-point quantity, not two independent constants on different objects.

By contrast, |c(φ)|² is *not* invariant on the cardioid:

    |c(φ)|² = 5/16 − (1/4)·cos(φ)

ranging from 1/16 (= 1/4² at φ = 0, the cusp) to 9/16 (= (3/4)² at φ = π, the tail). The Quarter b² = 1/4 equals |c|² only at the cusp; elsewhere |c| varies but the |z*|, |z*|² invariants hold.

### Role table

The four typed Pi2 anchors enter the cardioid story:

| Anchor | Role on cardioid |
|---|---|
| `HalfAsStructuralFixedPointClaim` (b = 1/2) | Magnitude \|z*\| invariant around the curve (argmax side) |
| `QuarterAsBilinearMaxvalClaim` (b² = 1/4) | Squared magnitude \|z*\|² invariant around the curve (maxval side); also \|c\| at the real-axis cusp only |
| `NinetyDegreeMirrorMemoryClaim` (i, 90°) | Complex-parameter generator that lifts c from real-axis to full complex plane |
| `PolynomialFoundationClaim` (d² − 2d = 0) | The c = 0 case where z* = 0 (degenerate fixed point) |

## Relation to F95

F95 is the φ = 0 specialization: at c = 1/4 (real-axis cusp), F95's θ = arctan(√(4c − 1)) gives θ = arctan(0) = 0, matching F97's arg(z*) = 0 at φ = 0.

For c slightly above 1/4 on the real axis (still in F95's domain), z± = b ± i·√(c − b²), so |z+| = √(b² + (c − b²)) = √c. Here |z+| ≠ b in general; F95's regime is *off the cardioid* but on the real axis.

The cardioid is the boundary of the period-1 attracting region; F95's c > 1/4 real-axis regime is *outside* the cardioid where the period-1 fixed point is repelling (|μ| > 1).

So F95 and F97 cover complementary regions:
- F95 (real c > b² = 1/4): repelling fixed point on the positive real axis, angle θ(c; b) of the complex-conjugate pair (which are the symmetric complex roots, not the marginal-stability fixed point on the cardioid).
- F97 (complex c on the cardioid): marginally stable fixed point, magnitude pinned to b, angle equal to the cardioid parameter φ.

Both share the same b = 1/2 anchor and emerge from the same z² − 2bz + c = 0 algebra; they project different aspects of its geometry.

## Connection to hardware data

The [`CPSI_COMPLEX_PLANE`](../../experiments/CPSI_COMPLEX_PLANE.md) hardware run on `ibm_kingston` (2026-04-16) observed two Bell⁺ pairs tracing 2D logarithmic spirals in the complex-c plane around the cusp at c = 1/4. The radial decay is set by γ (F25 Bell⁺ closed form, F57 K_dwell), the angular rotation by Ω (the residual Z-detuning, later actively steered in the [`f95_angle_steering_kingston_may2026`](../../docs/ANALYTICAL_FORMULAS.md#f95) Confirmation).

F97 places this experimentally observed 2D spiral into the cardioid framing: the hardware trajectory is a c(t) path in the complex plane that crosses the cardioid boundary (|z*| = b transition) before spiraling into the stable interior. The framework's `HalfAsStructuralFixedPoint` anchor is operationally the "where the period-1 stability transition lives", with the cardioid as its geometric locus in the c-plane.

## Numerical verification

`simulations/_cardioid_parametrization_tier1.py` traces 1000 φ values across [0, 2π] and verifies:

- Magnitude invariance: max |z*(φ)| − min |z*(φ)| over 1000 samples is 0 to machine precision (std ≈ 2.5 × 10⁻¹⁷).
- Algebraic identity c(φ) = z*(1 − z*): max residual 1.24 × 10⁻¹⁶ over 1000 φ (bit-exact).
- F95 cusp recovery: c(0) = 0.25 + 0i exactly; Im(c) = 0.0, Re(c) = 0.25.

## Anchors

- Numerical + algebraic verification: [`simulations/_cardioid_parametrization_tier1.py`](../../simulations/_cardioid_parametrization_tier1.py)
- F95 (companion: real-c angle): [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), [`F95 ANALYTICAL_FORMULAS entry`](../ANALYTICAL_FORMULAS.md#f95)
- Hardware 2D extension precursor: [`experiments/CPSI_COMPLEX_PLANE.md`](../../experiments/CPSI_COMPLEX_PLANE.md) (Kingston 2026-04-16; complex-CΨ already observed as 2D spirals)
- Februar boundary precursor: [`experiments/BOUNDARY_NAVIGATION.md`](../../experiments/BOUNDARY_NAVIGATION.md) (real-c θ-compass)
- Quarter-boundary roadmap (Layer 7 next-move slot): [`docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md`](PROOF_ROADMAP_QUARTER_BOUNDARY.md)
- Reflection that named the cardioid magnitude reading (2026-05-16): [`reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md`](../../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md)
- Mandelbrot connection synthesis: [`experiments/MANDELBROT_CONNECTION.md`](../../experiments/MANDELBROT_CONNECTION.md)
- F97 typed claim: [`F97CardioidHalfFixedPointPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F97CardioidHalfFixedPointPi2Inheritance.cs)
