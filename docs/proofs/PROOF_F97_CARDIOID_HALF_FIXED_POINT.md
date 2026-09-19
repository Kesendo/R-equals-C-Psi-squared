<!-- QUARTER-CURRENT -->
# Proof of F97: period-one cardioid parametrization

Current reading: choosing the marginal fixed point
`z*=exp(i*phi)/2` gives `c=z*-z*^2=exp(i*phi)/2-exp(2i*phi)/4` and
`|2z*|=1`.  Only this selected root is guaranteed marginal; the other root is
not generally marginal.  The formula neither inherits F95 nor turns a measured
radial trace into a recurrence orbit.

In mathematical notation: `c = e^{iφ}/2 − e^{2iφ}/4`.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the longer derivation and its earlier hardware analogies
follow; the current period-one statement is given above.

# Proof of F97: The Mandelbrot Cardioid Parametrization at Framework b = 1/2

**Statement:** The main cardioid of the Mandelbrot set is the locus in the
complex-c plane where a selected period-1 fixed point of `z²+c` has multiplier
magnitude one, equivalently `|z*|=1/2`. It admits the explicit parametrization

    c(φ) = b·e^(iφ) − b²·e^(2iφ)    for φ ∈ [0, 2π]

equivalently c(φ) = z*(φ) · (1 − z*(φ)) where z*(φ) = b·e^(iφ) is the period-1 fixed point. The framework's b=1/2 case gives

    c(φ) = (1/2)·e^(iφ) − (1/4)·e^(2iφ)

with two structural invariants on the curve:

    |z*(φ)| = b = 1/2   (selected fixed-point magnitude)
    arg(z*(φ)) ≡ φ (mod 2π)      (cardioid parameter)

The Python producer reports φ wrapped to [0,2π), while the typed `FixedPointArgument` API reports the principal `Atan2` branch; at φ=3π/2 it reports −π/2. These are two representatives of the same argument class.

[F95](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md) is a real-c, finite-b>0 root-angle
identity. F97 instead selects a marginal period-one root and solves for complex
c. At φ=0 both calculations contain the number 1/4, but neither claim inherits
the other.

**Status:** Tier 1 derived. The algebraic identity is exact; a separate 1000-point floating reconstruction has maximum residual 1.24×10⁻¹⁶.

**Date:** 2026-05-17.

---

## Abstract

F97 parameterizes the main cardioid by choosing a period-one fixed point z* of
z²+c with |z*|=1/2. Solving c=z*−z*² gives

    c(φ) = (1/2)·e^{iφ} − (1/4)·e^{2iφ} = z*(φ)·(1 − z*(φ)),   z*(φ) = (1/2)·e^{iφ},

with |z*|=1/2 and therefore |z*|²=1/4 around the curve. These are elementary
properties of the selected root, not typed ancestry for F97 and not statements
about the varying parameter magnitude |c|.

The selected root is marginal because μ=2z* and |μ|=1. The other quadratic
root is not generally marginal: at φ=π/2 its multiplier magnitude is √5.
Measured Bell-pair CΨ trajectories are different complex-valued readouts and
are not c-orbits or cardioid-boundary crossings.

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
    arg(z*(φ)) ≡ φ (mod 2π)                  (cardioid parameter)

∎

## Cardinal points

| φ | c(φ) | \|c\| | z*(φ) | reading |
|---|---|---|---|---|
| 0 | 1/4 | 1/4 = b² | 1/2 | real-axis cusp ([F95](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), [`BOUNDARY_NAVIGATION`](../../experiments/BOUNDARY_NAVIGATION.md)) |
| π/3 | 3/8 + i·√3/8 | √(9/64 + 3/64) = √3/4 ≈ 0.433 | (1 + i·√3)/4 | upper-right point; multiplier parameter phi=60°, arg(c)=30° |
| π/2 | 1/4 + i/2 | √(1/4 + 1/16) ≈ 0.559 | i/2 | imaginary-axis crossing of fixed point |
| π | −3/4 | 3/4 | −1/2 | real-axis "tail" of cardioid (period-doubling boundary) |
| 4π/3 | −1/8 − i·3√3/8 | √(1/64 + 27/64) = √7/4 ≈ 0.661 | (−1 − i·√3)/4 | conjugate image of φ = 2π/3 (240°) |
| 5π/3 | 3/8 − i·√3/8 | √(9/64 + 3/64) = √3/4 ≈ 0.433 | (1 − i·√3)/4 | conjugate image of φ = π/3 (300°) |

The 1/4 = b² value plays a special role only at φ = 0: it is the *magnitude* of c at the real-axis cusp. Elsewhere on the cardioid |c| ≠ b². The b = 1/2 anchor plays the universal role: |z*| = b along the entire curve.

## Structural reading and scope

The cardioid is the curve in complex c where the **selected** period-one root
has magnitude 1/2. This is a recurrence statement, separate from F95:

- F95 gives the upper-root angle for a real-coefficient positive-b quadratic.
- F97 chooses a marginal root and maps it to the complex recurrence parameter.

The real point c=+1/4 is the cardioid cusp. A radial set |CΨ_com|=1/4
contains many other phases and is not this cusp or this curve.

### One selected root, two elementary powers

For the selected marginal root, two elementary equalities hold:

    |z*(φ)| = 1/2
    |z*(φ)|² = 1/4

Both hold for all φ. Their arithmetic relation does not make the cardioid a
joint typed locus with other claims that happen to use 1/2 or 1/4.

The relevant distinction is `|z*|²=1/4` versus the parameter magnitude
`|c|`; only the former is constant on the curve.

By contrast, |c(φ)|² is *not* invariant on the cardioid:

    |c(φ)|² = 5/16 − (1/4)·cos(φ)

The squared magnitude ranges from 1/16 (= 1/4² at φ = 0, the cusp) to 9/16 (= (3/4)² at φ = π, the tail). The Quarter b² = 1/4 equals |c| only at the cusp; elsewhere |c| varies but the |z*|, |z*|² invariants hold.

### Former role table, now separated

The current typed F97 claim is parentless. The following numbers or operations
may resemble other registry entries, but they are not parent roles:

| Anchor | Role on cardioid |
|---|---|
| 1/2 | selected-root magnitude on the curve |
| 1/4 | selected-root squared magnitude; also c at φ=0 |
| complex phase | parameter used to traverse the curve |
| c=0 | recurrence roots 0 and 1, distinct rather than degenerate |

## Comparison with F95

At c=1/4, F95's real-coefficient root angle is zero and F97's selected
cardioid root also has argument zero. This is a shared evaluation, not a
specialization or ancestry relation.

For c slightly above 1/4 on the real axis (still in F95's domain), z± = b ± i·√(c − b²), so |z+| = √(b² + (c − b²)) = √c. Here |z+| ≠ b in general; F95's regime is *off the cardioid* but on the real axis.

The cardioid is the boundary of the period-1 attracting region; F95's c > 1/4 real-axis regime is *outside* the cardioid where the period-1 fixed point is repelling (|μ| > 1).

The two claims answer different questions:
- F95 (real c > b² = 1/4): repelling fixed point on the positive real axis, angle θ(c; b) of the complex-conjugate pair (which are the symmetric complex roots, not the marginal-stability fixed point on the cardioid).
- F97 (complex c on the cardioid): marginally stable fixed point, magnitude pinned to b, argument congruent to φ modulo 2π.

They share a quadratic form at b=1/2 but select different domains and objects.

## Hardware non-connection

The Kingston files contain reconstructed complex CΨ values for Bell-like
pairs. Their radial decay and measured phase are finite hardware readouts. They
are not the recurrence parameter c unless an additional map is supplied.

Therefore F97 does not classify those hardware paths as interior/exterior and
does not claim they cross the cardioid. The data can motivate a picture, but it
does not test this recurrence theorem.

## Numerical verification

`simulations/cardioid_parametrization_tier1.py` traces 1000 φ values across [0, 2π] and verifies:

- Magnitude invariance: max |z*(φ)| − min |z*(φ)| over 1000 samples is 5.55 × 10⁻¹⁷ (std ≈ 2.5 × 10⁻¹⁷).
- Algebraic identity c(φ) = z*(1 − z*): max floating residual 1.24 × 10⁻¹⁶ over 1000 φ.
- Endpoint control: c(0) = 0.25 + 0i exactly in symbolic arithmetic.

## Anchors

- Numerical + algebraic verification: [`simulations/cardioid_parametrization_tier1.py`](../../simulations/cardioid_parametrization_tier1.py)
- F95 (companion: real-c angle): [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), [`F95 ANALYTICAL_FORMULAS entry`](../ANALYTICAL_FORMULAS.md#f95)
- Complex-CΨ comparison/motivation only: [`experiments/CPSI_COMPLEX_PLANE.md`](../../experiments/CPSI_COMPLEX_PLANE.md) (the document supplies no map from that hardware readout to recurrence c and no test of F97)
- Februar boundary precursor: [`experiments/BOUNDARY_NAVIGATION.md`](../../experiments/BOUNDARY_NAVIGATION.md) (real-c θ-compass)
- Quarter-boundary roadmap (Layer 7 next-move slot): [`docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md`](PROOF_ROADMAP_QUARTER_BOUNDARY.md)
- Reflection that named the cardioid magnitude reading (2026-05-16): [`reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md`](../../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md)
- Mandelbrot connection synthesis: [`experiments/MANDELBROT_CONNECTION.md`](../../experiments/MANDELBROT_CONNECTION.md)
- F97 typed claim: [`F97CardioidHalfFixedPointPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F97CardioidHalfFixedPointPi2Inheritance.cs)
