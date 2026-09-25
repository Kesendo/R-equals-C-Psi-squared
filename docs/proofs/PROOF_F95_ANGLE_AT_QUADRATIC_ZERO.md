# Proof of F95: Angle Emergence at Quadratic Discriminant Zero

**Statement:** For a monic quadratic z² − 2bz + c = 0 with real c and finite b > 0, the principal argument of its upper-half-plane root (when c > b², the discriminant-negative regime) is:

    θ(c; b) = arctan( √(c/b² − 1) )    for c > b²
    θ = 0                              for c = b²
    θ undefined                        for c < b²

With the numerical specialization b=1/2, the threshold is b²=1/4 and the
formula becomes

    θ(c) = arctan( √(4c − 1) )         for c > 1/4

The same arithmetic appeared in the February angle table, but that finite
state readout and the recurrence are separate applications of a quadratic.

**Status:** Tier 1 derived. 4-line polynomial calculation. Against the Februar θ-compass table it reproduces the printed angles to their 0.1° rounding at five of six points; the sixth (c = 0.256) is off by 0.3° because c itself is printed rounded (c = (1 + tan² 9.1°)/4 ≈ 0.2564).

**Date:** 2026-05-16 (evening).

---

## Abstract

F95 names the angle that appears when a positive-b quadratic crosses its discriminant zero. For monic z² − 2bz + c with real c and finite b > 0, the complex root pair has arguments

    θ(c; b) = arctan(√(c/b² − 1))   for c > b²,   θ = 0 at c = b²,   undefined for c < b²,

and at b=1/2 this becomes θ(c)=arctan(√(4c−1)) beyond c=1/4.
The identity is scale-invariant under `(b,c)→(sb,s²c)` for s>0. A finite
angle table may use it as a coordinate without becoming a recurrence orbit.

The proved reading is narrower: a discriminant-negative positive-b quadratic
has a complex root, and its principal argument needs this second coordinate.
That algebra does not force complex amplitudes in a physical theory. F94's
Dyson coefficient and F97's period-one cardioid are independent claims with no
typed ancestry through F95.

## Setup

Take the monic quadratic over the reals:

    p(z) = z² − 2bz + c,    c ∈ ℝ, finite b > 0

with leading coefficient 1, linear-term half b, constant term c. The discriminant is

    D = (−2b)² − 4·1·c = 4(b² − c)

and the two roots are

    z± = (2b ± √D) / 2 = b ± √(b² − c)

The three regimes:

| sign(D) | sign(c − b²) | roots | angle structure |
|---|---|---|---|
| D > 0 | c < b² | two distinct real roots z₊ = b + √(b²−c), z₋ = b − √(b²−c) | no angle (purely real) |
| D = 0 | c = b² | one degenerate real root z = b (double) | angle = 0 (the boundary) |
| D < 0 | c > b² | complex conjugate pair z± = b ± i·√(c − b²) | angle = arctan(Im/Re) |

The angle emergence is the D < 0 case. In that regime the imaginary part is √(c − b²) and the real part of each root is b.

## Derivation (4 lines)

For c > b²:

```
z± = b ± i·√(c − b²)               (1)  quadratic formula in the D < 0 regime

Re(z+) = b
Im(z+) = √(c − b²)                  (2)  real and imaginary parts

arg(z+) = arctan( Im(z+) / Re(z+) )
        = arctan( √(c − b²) / b )    (3)  definition of complex argument

        = arctan( √( (c − b²)/b² ) )
        = arctan( √( c/b² − 1 ) )     (4)  algebraic simplification ∎
```

The convention here picks the upper-half-plane root z₊ (positive imaginary part); because b > 0, its principal argument lies in [0, π/2), and the algebraic step from √(c−b²)/b to √((c−b²)/b²) is valid. The lower-half-plane root z₋ is the complex conjugate, with argument −arctan(√(c/b² − 1)), reflecting through the real axis. Both roots together carry a single magnitude (|z±| = √(b² + (c − b²)) = √c) and a paired ±θ angle. The b < 0 branch would require the second-quadrant principal argument π−θ, while b = 0 gives π/2 for c > 0; neither branch is part of the current F95 API.

## Numerical specialization (b = 1/2)

Substituting b=1/2 gives:

    threshold:  b² = 1/4
    formula:    θ(c) = arctan( √(c/(1/4) − 1) ) = arctan( √(4c − 1) )

which numerically reproduces the February table when its scalar is used as c.

## Numerical verification

The script [`simulations/angle_at_zero_tier1_candidate.py`](../../simulations/angle_at_zero_tier1_candidate.py) implements `angle_at_quadratic_crossing(c, b)` and checks against the Februar θ-compass table:

| c used in the table | Februar θ (deg) | F95 derived (deg) | match |
|---|---|---|---|
| 1/3 = 0.3333 | 30.0 | 30.0000 | ✓ exact |
| 0.308 | 25.8 | 25.7184 | ✓ within 0.1° |
| 0.286 | 20.7 | 20.7804 | ✓ within 0.1° |
| 0.266 | 14.1 | 14.1969 | ✓ within 0.1° |
| 0.256 | 9.1 | 8.8062 | ✗ 0.3° drift (the printed c is rounded, not a formula error) |
| 0.250 | 0.0 | 0.0000 | ✓ exact (threshold) |

The single drift at 0.256 comes from the rounded c the table prints: near
¼ the angle is steep in c, and 9.1° corresponds to c = (1 + tan² 9.1°)/4 ≈ 0.2564.

## Scope: a quadratic identity, not a genealogy

The identity applies to any finite real b>0 and real c in the stated domain.
One may substitute b=1 or b=1/2, but a shared polynomial shape does not make
another claim a parent or a physical realization of F95.

For the recurrence `z²-z+c=0`, b=1/2 is one genuine application. The
TransitionBridge decay-root quadratic is another. A CΨ angle table is only a
coordinate application unless its physical equation is separately supplied.

The angle is the principal argument of one selected quadratic root. No
conclusion about a "quadratic state" or physical amplitude follows from that
fact alone.

## Anchors

- Numerical verification: [`simulations/angle_at_zero_tier1_candidate.py`](../../simulations/angle_at_zero_tier1_candidate.py)
- F-registry entry: [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md) §F95
- Companion reflection: [`reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](../../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md)
- Februar Mandelbrot-specific precursor: [`experiments/BOUNDARY_NAVIGATION.md`](../../experiments/BOUNDARY_NAVIGATION.md)
- Typed claim: [`F95AngleAtQuadraticZeroPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs), currently parentless despite its historical class name
- Independent comparison: [F94](../ANALYTICAL_FORMULAS.md#f94), a named-ring Dyson coefficient rather than a magnitude-side sibling
