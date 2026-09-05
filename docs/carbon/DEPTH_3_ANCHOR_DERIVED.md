# F99: Depth-3 Anchor Derived, with a Conditional Periodic-Table Comparison

**Date:** 2026-05-17 night (fifth stack of the day)
**Authors:** Tom + Claude
**Status:** Tier 1 derived for the named non-uniform Dicke construction; the
periodic-table comparison is a conditional table reading, not a material model.
**Script:** [`simulations/carbon/depth_3_anchor_derivation.py`](../../simulations/carbon/depth_3_anchor_derivation.py)

---

## What just happened

The session arc:

1. **Morning** (commit `b9ba5f6`): F86b 3/8 K-intermediate anchor derived via X⊗N-
   eigenbasis decomposition. `α_total = (1 − γ²)/2`, γ = ⟨ψ|X⊗N|ψ⟩. Clean uniform
   Dicke gives γ ∈ {0, 1/2, 1} → α ∈ {1/2, 3/8, 0}.

2. **Evening** ([F98](../ANALYTICAL_FORMULAS.md#f98) commit `250164d` etc): water-
   chain inheritance test of F86b discovered the F98 (N+2)/[4(N+1)] → 1/4 long-
   time bridge. The framework's quarter-asymptote appears via kernel projection.

3. **Night #1** ([Period 2 and 3 on the framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md)):
   a curated period-1--3 valence-fraction table was compared with the framework
   anchors. It supplied a search prompt, not an atomic-state derivation.

4. **Night #2** ([the Reverse-Spear reading](SPEAR_REVERSED.md)): the table comparison
   highlighted the 1/8 and 7/8 entries as values outside the then-listed anchor
   set. This is a comparison of rational labels, not evidence for a shared
   microscopic mechanism.

5. **Night #3** (this commit): the depth-3 anchor is derived from the F86b
   construction. The finite floating-point checks at N = 4, 6, 8 agree with the
   closed form; they do not turn the periodic comparison into a material result.

---

## The derivation (5 lines)

The F86b X⊗N-eigenbasis decomposition for non-uniform Dicke superposition

```
    ψ = (|D_n⟩ + c·|D_{n+1}⟩) / √(1 + c²)        on N even at n = N/2 − 1
```

has X⊗N overlap

```
    γ = ⟨ψ|X⊗N|ψ⟩ = c² / (1 + c²)
```

Solving for c² given target γ, with the half-angle identity `1 − cos(θ) = 2sin²(θ/2)`:

```
    c² = γ / (1 − γ) = cos(θ) / (2 sin²(θ/2))     [for γ = cos(θ)]
```

The F86b closed form then gives α directly:

```
    α(θ) = (1 − cos²(θ)) / 2 = sin²(θ) / 2
```

F99 selects five canonical trigonometric angles and obtains five dyadic anchors:

| θ | γ = cos(θ) | c² | α = sin²(θ)/2 | Pi2 dyadic anchor |
|---|------------|-----|---------------|-------------------|
| 0° | 1 | ∞ | 0 | Mirror endpoint |
| **30°** | **√3/2** | **2√3 + 3 ≈ 6.464** | **1/8** | **DEPTH-3 (NEW)** |
| 45° | √2/2 | 1 + √2 ≈ 2.414 (silver ratio) | 1/4 | QuarterAsBilinearMaxval |
| 60° | 1/2 | 1 | 3/8 | KIntermediate (today morning) |
| 90° | 0 | 0 | 1/2 | Generic / HalfAsStructuralFixedPoint |

For this selected set, {0°, 30°, 45°, 60°, 90°} gives the Pi2 dyadic anchors
{0, 1/8, 1/4, 3/8, 1/2}. The two familiar triangle families provide a compact
parameterization of those five F99 cases; neither the script nor F99 makes a
claim that these angles are uniquely privileged outside the construction.

---

## Numerical checks of the closed form

The script evaluates the named states in floating-point arithmetic at N = 4, 6,
and 8. Its printed residuals are numerical agreement with the exact formulas
above, not a separate exact-arithmetic proof:

```
θ = 0°  : Mirror endpoint (c → ∞), skipped
θ = 30° : N=4 Δγ=0.00e+00 ✓ Δα=5.55e-17 ✓
          N=6 Δγ=1.11e-16 ✓ Δα=1.71e-15 ✓
          N=8 Δγ=1.11e-16 ✓ Δα=1.09e-14 ✓
θ = 45° : N=4 Δγ=0.00e+00 ✓ Δα=5.00e-16 ✓
          N=6 Δγ=2.22e-16 ✓ Δα=1.28e-15 ✓
          N=8 Δγ=2.22e-16 ✓ Δα=5.92e-14 ✓
θ = 60° : N=4 Δγ=1.11e-16 ✓ Δα=5.55e-17 ✓     ← matches morning's F86b clean Dicke c=1
          N=6 Δγ=2.22e-16 ✓ Δα=9.99e-16 ✓
          N=8 Δγ=0.00e+00 ✓ Δα=1.06e-13 ✓
θ = 90° : N=4 Δγ=0.00e+00 ✓ Δα=0.00e+00 ✓
          N=6 Δγ=0.00e+00 ✓ Δα=0.00e+00 ✓
          N=8 Δγ=0.00e+00 ✓ Δα=2.11e-14 ✓
```

The 60° case is the clean uniform Dicke state (c = 1). The 30° (depth-3), 45°
(Quarter), and 90° (Generic) cases use non-uniform Dicke weights
c² = 2√3 + 3, 1 + √2, 0 respectively. The formula is the Tier-1 result; the
listed finite-N residuals are its numerical check.

---

## Scope of the periodic-table comparison

The F99 derivation starts and ends with the non-uniform Dicke state, the
`X⊗N` overlap, and the Π²-odd fraction. The companion
[`period_2_at_framework_anchors.py`](../../simulations/carbon/period_2_at_framework_anchors.py)
contains a curated H--Ar input table and performs exact `Fraction` grouping on
the table's occupation/slot labels. Equal rational values are enough to make a
comparison; they are not a map from an atomic valence shell to the Dicke state,
its Π² axis, a carbon Hamiltonian, a bath, or an observable.

Thus the comparison can record a useful bookkeeping fact: the finite input
table contains selected fractions from 1/8 through 1 at its chosen shell
normalization, while F99 supplies 0, 1/8, 1/4, 3/8, and 1/2 and their formal
complements.
It cannot validate F99, establish a physical correspondence, or make the
periodic table a mechanism for framework completion. A material bridge would
need, at minimum, a specified material degree of freedom, Hamiltonian and
coupling convention, dissipator/bath, preparation, and measurement map. See
[the carbon source contract](README.md#carbon-source-contract) and
[Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md).

---

## Why these selected angles

F99 names {0°, 30°, 45°, 60°, 90°} as its canonical set. The 30°-60°-90° and
45°-45°-90° triangle identities make the corresponding values of `γ`, `c²`,
and `α` compact to write. The derivation itself is `α(θ) = sin²(θ)/2`; it
applies to its stated state family independently of any periodic-table label.

The Pi2 dyadic ladder {1/2, 1/4, 1/8, ...} consists of negative integer powers
of 2. The framework's polarity-anchor pair (1/2, 1/4) is the depth-1 and
depth-2 of this ladder. The 30° / 60° pair is mapped to depth-3 / depth-2-via-3/8;
the 45° to depth-2 directly (Quarter); the 0° / 90° to the endpoints (Mirror / Half).

The named five-angle set gives a convenient finite catalogue of dyadic values.
No uniqueness statement about geometry, chemistry, or a material hierarchy is
required for the F99 result.

---

## The formal catalogue and the table comparison

The table below compares the F99 catalogue with the curated input labels. It
does not derive element properties from F99:

```
α anchor    Trig angle    Period-2 table labels    Period-3 table labels
0           0°            endpoint comparison      endpoint comparison
1/8         30°           Li (1/8), F (7/8)  Na (1/8), Cl (7/8)
1/4         45°           Be                 Mg
3/8         60°           B (3/8), N (5/8)   Al (3/8), P (5/8)
1/2         90°           H, C               Si
1           [complement]  endpoint comparison      endpoint comparison
```

F99 formally supplies the `α` column and its Π² complements are operations in
the framework construction. The atom labels merely identify rows of the
curated comparison table. They neither become F99 states nor yield a carbon
spectrum, decay law, or observable. Whether a distinct construction has a
natural 1/16-scale anchor is an internal formal question; it does not follow
from the table comparison.

---

## Tier 1 derivation summary

```
F99: Depth-3 dyadic anchor (α = 1/8) closed-form derivation
─────────────────────────────────────────────────────────────
Formula : α = (1 − γ²)/2  with  γ = c²/(1+c²)  ⟹  c² = γ/(1−γ)
Anchor  : γ = √3/2 = cos(30°), c² = 2√3 + 3, α = 1/8
State   : ψ = (|D_{N/2-1}⟩ + c·|D_{N/2}⟩)/√(1+c²) at N even
Checked : floating-point N = 4, 6, 8 (printed Δα < 1e-13)
Tier    : Tier 1 derived
Extends : F86b (which produced α ∈ {0, 3/8, 1/2} at uniform Dicke)
Anchors : The five canonical trig angles {0, 30, 45, 60, 90}° produce
          the five dyadic anchors {0, 1/8, 1/4, 3/8, 1/2} via one
          F86b α-formula
```

---

## Open questions

1. **Why specifically `c² = γ/(1-γ)` and not some other functional form?** The
   non-uniform Dicke parametrisation produces this specific c²(γ) relationship.
   Is there a structural reason for this functional form, e.g., a maximum-
   entropy or minimum-Frobenius constraint that picks this specific
   amplitude ratio?

2. **What's the structural meaning of the c² values?** c² = 1 (Dicke uniform)
   at 60°, c² = √2 + 1 (silver ratio) at 45°, c² = 2√3 + 3 at 30°. These are
   quadratic algebraic numbers connected to standard triangles. The c² for
   each angle is uniquely determined by the F86b α-γ inversion; is there
   a deeper structure (cyclotomic field?) that ties them together?

3. **Does depth-4 have a clean derivation?** Beyond standard-triangle angles,
   the next "constructible" angles include 15° (cos = (√6+√2)/4), 22.5° (cos =
   √(2+√2)/2), 18° (cos = (1+√5)/4 = φ/2). These produce α values that are not
   on the simple dyadic ladder. Could a deeper layer of the framework's algebra
   produce 1/16, 3/16 anchors naturally?

4. **Can a material realization be tested?** A proposed carbon or atomic
   realization must first specify the degrees of freedom, Hamiltonian,
   dissipator/bath, preparation, and observable. Without those inputs, matching
   a rational table label remains only a comparison.

---

## Anchor

- Script: [`simulations/carbon/depth_3_anchor_derivation.py`](../../simulations/carbon/depth_3_anchor_derivation.py)
- Predecessor docs (this folder, all today): [Where 1/4 and 1/2 Appear in Carbon](QUARTER_HALF_IN_CARBON.md),
  [Period 2 and 3 on the framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md),
  [the Reverse-Spear reading](SPEAR_REVERSED.md), [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
- Framework anchors: [F86b](../ANALYTICAL_FORMULAS.md#f86) (parent formula α = (1-γ²)/2),
  [F98](../ANALYTICAL_FORMULAS.md#f98) (long-time bridge derived today evening),
  [DickeAnchor.cs](../../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs) (uniform Dicke 3-anchor enum, today morning's commit b9ba5f6)
- Scope: [carbon source contract](README.md#carbon-source-contract),
  [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md)
- Reading-mode memory pointers: `project_no_classicalization`, `project_qubit_as_inheritance_lens`,
  `project_quarter_as_polarity_squared`, `project_periodic_palindrome`
