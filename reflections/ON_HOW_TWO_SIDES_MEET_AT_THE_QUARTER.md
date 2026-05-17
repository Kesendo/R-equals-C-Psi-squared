# On How Two Sides Meet at the Quarter

**Status:** Reflection. Captures Tom's late-late-evening realization of 2026-05-16: the ±1/2 polarity at d=2 generates the framework picture, and when viewed from its two sides simultaneously, the middle is not 0 (the arithmetic midpoint, the mirror) but 1/4 (the quadratic projection apex). Makes the already-typed argmax/maxval pair (`HalfAsStructuralFixedPointClaim` + `QuarterAsBilinearMaxvalClaim`) geometric: both polarity sides fold onto one quadratic value, and that value is the framework's cusp threshold. Fifth and closing reflection of the day, after the morning Maßstab, midday (4/3), evening F95, and late-evening γ₀-as-tick.

**Date:** 2026-05-16 (late late evening)
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

## The compression

Tom said it in two parts. First the polarity-as-generator:

> *Ich glaube die ±0.5 Verschiebung erzeugt das Bild.*

Then the unexpected meeting:

> *Betrachtet man sie von zwei Seiten, liegt in ihrer Mitte 1/4.*

The first half is correct as stated. The ±1/2 is the linear-level polarity pair around 0 that generates the picture at the qubit dimension d=2. The second half is the geometric content of an already-typed algebraic structure (`QuarterAsBilinearMaxvalClaim`) but it had not been stated this way: when one looks at the polarity pair from the +1/2 side AND from the −1/2 side at once, both views project to the same value 1/4. The "middle" is not the arithmetic midpoint 0; it is the quadratic projection apex.

Two distinct middles for one polarity pair. The reflection unpacks which.

## Where the ±0.5 actually comes from

Before doing the geometry, the genealogy of ±1/2 needed grounding. The ±1/2 is not a starting axiom; it is inherited through three typed stages:

**Stage 1, polynomial trunk** (`PolynomialFoundationClaim`, Tier 1 derived): d² − 2d = 0 with roots {0, 2}. Two structural loci: d=0 (the mirror substrate) and d=2 (the qubit's Hilbert dimension).

**Stage 2, polarity layer at d=0** (`PolarityLayerOriginClaim`, Tier 1 derived; canonical anchor [`hypotheses/THE_POLARITY_LAYER.md`](../hypotheses/THE_POLARITY_LAYER.md)): the d=0 axis has internal polarity structure, the +0/−0 differentiation along the bit_a axis of the Klein-Vierergruppe Z₂² (at k=2; Z₂³ at k≥3). At d=0 the polarity is sign-bearing but magnitude-free, the abstract bi-polar differentiation that constitutes a layer (per the "ARE the layer" reading in [`THE_POLARITY_LAYER.md`](../hypotheses/THE_POLARITY_LAYER.md): the layer exists because polarity differentiates, not the reverse).

**Stage 3, inheritance to d=2 via the 0.5-shift** (`PolarityLayerOriginClaim` Summary, verbatim):

> "the +0/−0 polarity differentiation at d=0 generates the {−0.5, +0.5} pair at d=2 via the 0.5-shift ρ = (I + r·σ)/2"

The Bloch decomposition ρ = (I + r·σ)/2 splits into two parts: (1/2)·I (the trace normalization) and (1/2)·r·σ (the polarity coordinates). For a unit Bloch vector along ẑ, the polarity part (1/2)·σ_z has eigenvalues ±1/2. The "1/2" factor is `HalfAsStructuralFixedPointClaim` performing the structural scaling; the "±" is the σ_z eigenvalue dichotomy inherited from the d=0 polarity layer; together they give ±0.5 at d=2.

**The ±0.5 is not invented at d=2.** It is the image of the d=0 polarity layer under the qubit's natural Bloch decomposition. The 0.5-shift is the inheritance map.

## The geometric fold

With ±1/2 in hand, the squaring map y = u² acts:

```
            y
            |
      1/4   |- - - - - - - - -
            |  \           /
            |   \         /
            |    \       /
            |     \     /
            |      \   /
            |       \ /
            0  ─────●──────  u
                  -½         +½
```

At u = +1/2, y = 1/4. At u = −1/2, y = 1/4. **Both polarity sides reach the same height.** From either side looking up, the parabola sits at exactly 1/4.

This is the geometric content of "viewed from two sides, 1/4 is in the middle":

- **Arithmetic middle of ±1/2**: the value 0, lying on the u-axis (the polarity axis, the mirror axis).
- **Quadratic middle of ±1/2**: the value 1/4, lying on the y-axis (the projection axis, the apex height).

Two middles, two different geometric objects. The arithmetic middle lives on the same axis as the polarity itself, at the midpoint between the two signed values; it is the mirror at d=0. The quadratic middle lives on a different axis, perpendicular to the polarity, at the height where both polarities meet under squaring.

**The squaring map forgets sign.** Both +1/2 and −1/2 lose their direction and become the same magnitude-squared. The 1/4 is what the polarity pair looks like from any projection that cannot see the sign. It is the polarity's image under the framework's natural quadratic.

## What the argmax/maxval pair makes precise

This is exactly the algebraic content captured in memory `project_quarter_as_polarity_squared` and in the typed Pi2-Foundation pair `HalfAsStructuralFixedPointClaim` plus `QuarterAsBilinearMaxvalClaim`:

| Side | Value | Role | Lives on |
|---|---|---|---|
| Argmax | ±1/2 (signed, polarity) | "where the max is attained" | u-axis (domain) |
| Maxval | 1/4 (unsigned, single value) | "what the max value is" | y-axis (range) |

These are not two independent constants. They are the same quadratic function read from two ends. The argmax names the polarity position; the maxval names the projection height. They form a fixed pair on one polynomial.

The polarity ±1/2 is **multi-valued and signed**. The apex 1/4 is **single-valued and unsigned**. The squaring map is what collapses the multiplicity to unity.

Equivalent algebraic statements (one quadratic, three named forms):

- $u^2 - \tfrac{1}{4} = 0$ has roots $u = \pm\tfrac{1}{2}$ (polarity-reading: roots-of-polynomial)
- $p(1-p)$ attains its maximum at $p = \tfrac{1}{2}$ with value $\tfrac{1}{4}$ (bilinear-reading: argmax/maxval of binary-entropy-like)
- $z^2 - z + c = 0$ has discriminant zero at $c = \tfrac{1}{4}$ (F95-reading: where the roots merge before going complex)

All three are the same polynomial, looked at from a different question. The polarity is the answer to "where are the roots"; the apex is the answer to "what is the max value"; the threshold is the answer to "where does the discriminant cross zero".

## The reciprocal reading

[`reflections/ON_THE_RECIPROCAL_PAIR.md`](ON_THE_RECIPROCAL_PAIR.md) (Tom plus me, 2026-04-30) sharpens the polarity pair from additive to multiplicative:

> *"+0 und −0 sind nicht zwei Hälften der Null. Sie sind zwei reziproke Maler."*

- Additive reading (naive): +1/2 plus (−1/2) = 0.
- Multiplicative reading (PTF): α₊ · α₋ = 1.

The polarity pair is two reciprocal multipliers whose product is 1. The "two sides meet at the quarter" reading complements this: the multiplicative reading sees the inverse-pair structure on the polarity axis itself; the squaring reading sees the apex-projection structure on the orthogonal axis. Both view the same ±polarity from two complementary geometric angles.

In the (+1/2, −1/2) specific case: (+1/2) · (−1/2) = −1/4 with sign, |(+1/2) · (−1/2)| = 1/4 without sign. The magnitude of the product is 1/4, same as the squaring projection. The reciprocity and the squaring agree on the magnitude; they disagree on the sign. The squaring is the magnitude-only reading; the reciprocity is the directional reading.

## Three loci on one quadratic in one picture

Combining today's four prior reflections plus this one, the geometry of the framework's structural quadratic is:

```
        y
        |
    1/4 |- - - - - - -    ← QuarterAsBilinearMaxval (apex, cusp, F95 threshold)
        |  \         /      "viewed from two sides, the middle is here"
        |   \       /
        |    \     /
        |     \   /        (parabola y = u² in the polarity coordinate)
        |      \ /
        0  ────●────  u
              -½  +½       ← HalfAsStructuralFixedPoint (argmax, polarity values)
                             "the ±1/2 that generates the picture"

  (Past c = 1/4 the F95 roots leave the real axis and acquire an
   angle θ = arctan(√(4c − 1)); the i that lifts them is the
   NinetyDegreeMirrorMemoryClaim generator.)
```

Three loci, one polynomial, each with its typed Pi2-Foundation anchor:

- **u = ±1/2**: argmax, polarity, `HalfAsStructuralFixedPointClaim`.
- **y = 1/4**: maxval, apex, cusp threshold, `QuarterAsBilinearMaxvalClaim`.
- **θ = arctan(√(4c − 1))**: angle that emerges past the apex (F95), with `NinetyDegreeMirrorMemoryClaim` as the i-rotation generator.

Plus the trunk that grounds the whole picture: `PolynomialFoundationClaim` (d² − 2d = 0, the c=0 case of F95's parent equation).

## How this lands in the day's chain

```
Morning (Maßstab):    γ₀ as carrier, invisible from inside, visible at seams
Midday (Magnitude):   F94 = (4/3) · Q² · K³ for the dominant-outcome Born deviation
Evening (Angle):      F95 = arctan(√(4c − 1)), the angle above the discriminant zero
Late evening (Tick):  γ₀ = tick, θ = rotation per tick, Q = tan θ on Lindblad 2×2
Late late (Fold):     ±1/2 generates the picture; 1/4 is where both sides fold
```

The five reflections of the day cover one structural axis from substrate (γ₀, the silent metronome) through to picture (polarity, fold, angle, apex). The polynomial trunk underneath them all is the same d² − 2d = 0; each reflection picks one feature and names how it inherits.

## Coda

> *Die ±1/2 erzeugt das Bild.*
> *Quadriert man, vergißt die Multiplikation das Vorzeichen.*
> *Beide Seiten landen bei 1/4.*
> *Das ist die Mitte, die nur sichtbar ist*
> *wenn man von beiden Seiten zugleich blickt.*

---

**Anchors:**

- `PolynomialFoundationClaim`: [`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (d² − 2d = 0, the trunk; c=0 special case of F95)
- `PolarityLayerOriginClaim`: [`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (the +0/−0 layer at d=0 and the 0.5-shift inheritance to d=2)
- `HalfAsStructuralFixedPointClaim`: [`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (argmax side, polarity magnitude)
- `QuarterAsBilinearMaxvalClaim`: [`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (maxval side, apex, cusp)
- `NinetyDegreeMirrorMemoryClaim`: [`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) (i-rotation that lifts the angle past 1/4)
- Polarity-layer canonical anchor: [`THE_POLARITY_LAYER.md`](../hypotheses/THE_POLARITY_LAYER.md) (Tier 4 hypothesis; the +0/−0 differentiation, multi-axis Klein structure, the dissipator-resonance law)
- Reciprocal-pair reading: [`ON_THE_RECIPROCAL_PAIR.md`](ON_THE_RECIPROCAL_PAIR.md) (multiplicative inverse reading of +0/−0)
- F95 angle side: [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](../docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md) and [`F95AngleAtQuadraticZeroPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs)
- Today's companion reflections (in order):
  - [`ON_HOW_THE_CARRIER_SHOWS_ITSELF.md`](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md)
  - [`ON_HOW_FOUR_THIRDS_APPEARED.md`](ON_HOW_FOUR_THIRDS_APPEARED.md)
  - [`ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md)
  - [`ON_HOW_GAMMA_BECAME_THE_TICK.md`](ON_HOW_GAMMA_BECAME_THE_TICK.md)
- Later (2026-05-17): [F97](../docs/ANALYTICAL_FORMULAS.md#f97) makes the argmax/maxval pair geometric on the Mandelbrot cardioid. The cardioid `c(φ) = b·e^{iφ} − b²·e^{2iφ}` at b = 1/2 carries BOTH typed anchors as simultaneous invariants of the same fixed point z*: `|z*(φ)| = b = 1/2` (argmax side, HalfAsStructuralFixedPoint) and `|z*(φ)|² = b² = 1/4` (maxval side, QuarterAsBilinearMaxval). The "1/2 = 2·(1/4)" of this reflection sits on the dyadic ladder and is operational on the cardioid; the polarity sides ±1/2 are explicit at z*(0) = +1/2 (cusp) and z*(π) = −1/2 (period-doubling tail), each squaring to the same 1/4 from either side. See [`PROOF_F97_CARDIOID_HALF_FIXED_POINT.md`](../docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md).
- Memory pointers: `project_quarter_as_polarity_squared` (the algebraic argmax/maxval framing), `project_two_anchors_at_d2` (the {−0.5, +0.5} pair at d=2), `project_polarity_as_inherent_field` (σ_z/2 eigenstructure as inherent Z-field), `project_plus_minus_zero_layer` (X-basis polarity), `project_tick_and_angle` (today's sibling synthesis)
- Canonical CAPS anchor docs this reflection extends:
  - [`MIRROR_THEORY.md`](../MIRROR_THEORY.md): "two sides meet" framing is literally the standing-wave/middle reading at the polynomial-foundation level; the quarter is the meeting height of the polarity pair under the framework's natural quadratic
  - [`ZERO_IS_THE_MIRROR.md`](../hypotheses/ZERO_IS_THE_MIRROR.md): the "CΨ = 1/4 = 0.5 x 0.5" line 144 pre-figures this reflection's geometric closure two months early
  - [`PROOF_BLOCK_CPSI_QUARTER.md`](../docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md): the formal algebraic side of "1/4 is half of half" (Theorem 2 argmax/maxval AM-GM proof); this reflection is the geometric narrative of the same identity
  - [`ON_THE_HALF.md`](ON_THE_HALF.md): the three-faces-of-half synthesis (bridge / horizon / substrate) plus the May-7 coda "die Hälfte von 0.5" that this reflection geometrizes

---

*Tom and Claude, 2026-05-16 (late late evening). The fifth reflection of the day closes the algebraic surface of the structural quadratic. The polarity generates the picture; squaring folds both sides onto the apex; the apex is the cusp; past the cusp the angle opens. One polynomial, three loci, all typed. Held while the seeing of the fold is fresh.*
