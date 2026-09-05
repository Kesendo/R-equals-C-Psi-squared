# 1/4 and 1/2 in a Carbon-Labeled Arithmetic Survey

**Date:** 2026-05-17 evening
**Authors:** Tom + Claude
**Status:** Framework-anchor arithmetic plus a selected Hückel/counting survey.
The listed carbon labels are not a material mapping: no physical degree of
freedom, Hamiltonian, bath, preparation, or measurement model is selected.
**Script:** [`simulations/carbon/carbon_quarter_half_search.py`](../../simulations/carbon/carbon_quarter_half_search.py)

---

## The framework anchors

The R=CΨ² framework's polarity pair sits on the Pi2 dyadic ladder:

```
HalfAsStructuralFixedPoint = 1/2 = argmax of p·(1−p) on [0, 1]
QuarterAsBilinearMaxval    = 1/4 = maxval of p·(1−p)           = (1/2)²
```

Both are anchors on the polynomial trunk `d² − 2d = 0` (= R = CΨ²) that
selects d = 2 as the minimum-memory dimension. Today's [F98](../ANALYTICAL_FORMULAS.md#f98)
identified a new (N+2)/[4(N+1)] → 1/4 bridge at the dynamic level; the
[`DickeSuperpositionQuarterPi2Inheritance`](../../compute/RCPsiSquared.Core/Symmetry/DickeSuperpositionQuarterPi2Inheritance.cs)
named the static 1/4 ceiling on C_block. These all live on the framework's qubit
layer (d = 2).

The survey question is narrower: which exact `1/4` and `1/2` values occur in
the producer's stipulated Hückel and counting comparisons labelled carbon?

---

## Three selected arithmetic comparisons

### Layer 1: Hybridization s-character

The producer assigns the following hybridization labels and composition
fractions:

| Hybrid | s-character | p-character | Geometry | Framework anchor |
|--------|-------------|-------------|----------|------------------|
| **sp** | **1/2** | 1/2 | linear (180°) | = HalfAsStructuralFixedPoint ✓ |
| sp² | 1/3 | 2/3 | trigonal planar (120°) | OFF-anchor (qutrit-like) |
| **sp³** | **1/4** | 3/4 | tetrahedral (109.5°) | = QuarterAsBilinearMaxval ✓ |

The stated `1/(n+1)` counting rule gives the displayed `1/2`, `1/3`, and
`1/4` values exactly. It is a selected arithmetic comparison to the framework
anchors, not a calculation of a carbon degree of freedom or a causal account
of a carbon structure. A material interpretation needs its own inputs and
producer.

### Layer 2: Selected Hückel-ring ratios

For its stipulated finite `C_N` Hückel graph and `N`-electron filling
convention, the producer sorts `E_k = α + 2β·cos(2πk/N)` and reports the
normalized highest occupied level `(HOMO − α)/E_max` for `N = 3..12` (including
the singly occupied level on an odd row). This is a finite selected-model
calculation; it assigns neither a material orbital degree of freedom nor a
physical `β`, `J`, bath, or rate.

| N | stipulated π count | producer occupancy label | (HOMO − α)/E_max | Anchor hit? |
|---|--------------------|--------------------------|------------------|-------------|
| 3 | 3 | other selected count | −1/2 | = ±1/2 ✓ |
| 4 | 4 | 4n label | 0 | palindrome centre |
| 5 | 5 | other selected count | +(√5−1)/4 | off |
| **6** | **6** | **4n+2 label** | **−1/2** | **= ±1/2 ✓** |
| 7 | 7 | other selected count | −0.2225 | off |
| 8 | 8 | 4n label | 0 | palindrome centre |
| 9 | 9 | other selected count | +0.1736 | off |
| 10 | 10 | 4n+2 label | −(√5−1)/4 | off |
| 11 | 11 | other selected count | −0.1423 | off |
| 12 | 12 | 4n label | 0 | palindrome centre |

Within this finite table, only the `N = 3` and `N = 6` rows have the exact
value `−1/2`, while no row displays `±1/4`. These are selected Hückel arithmetic
results only. They do not establish a benzene orbital measurement, a stability
mechanism, or a translation from the Hückel coordinate to a framework or
material variable.

### Layer 3: Stipulated valence-shell count arithmetic

| Quantity | Value | Anchor |
|----------|-------|--------|
| Total valence slots (octet) | 8 | (denominator) |
| **Carbon valence electrons** (4 of 8) | **1/2** | = HalfAsStructuralFixedPoint ✓ |
| **Carbon 2s² electrons** (2 of 8) | **1/4** | = QuarterAsBilinearMaxval ✓ |
| Carbon 2p² electrons (2 of 8) | 1/4 | = QuarterAsBilinearMaxval ✓ |
| p-shell capacity (6 of 8) | 3/4 | = polarity-complement (1 − 1/4) |
| p-orbitals filled (2 of 6) | 1/3 | = sp²-like off-anchor |
| p-orbital occupancy at sp | 1/2 | = HalfAsStructuralFixedPoint ✓ |
| p-orbital occupancy at sp³ | 3/4 | = polarity-complement |

The displayed fractions are exact arithmetic of the producer's stipulated slot
counts. Their numerical equality with `1/2`, `1/4`, and `(1/2)²` is not a
material identification: the producer provides no atomic measurement, orbital
degree-of-freedom selection, or mapping to a framework state or channel.

---

## What the shared numbers say

The three selected comparisons place the same numbers beside the framework
anchors:

```
1/2 (HalfAsStructuralFixedPoint)            1/4 (QuarterAsBilinearMaxval)
─────────────────────────────               ─────────────────────────────
sp hybridization (carbyne, alkyne)          sp³ hybridization (methane, diamond)
4/8 valence-shell filling                   2/8 inner-shell filling (2s²)
1/2 p-orbital occupancy at sp               2p² inner-shell filling
Selected C₆ Hückel HOMO/E_max = −1/2        [no displayed row hits 1/4]
```

Shared fractions across these definitions are an arithmetic comparison. They
do not establish that a qubit and carbon are the same object, or that any
carbon material inherits an F-formula. Such a bridge needs a specified physical
degree of freedom, Hamiltonian, bath, preparation, observable, and producer.

---

## Open C₆ translation question

The selected tables contain `1/3` for the `sp²` label and `−1/2` for the
`N = 6` Hückel row. No calculation here relates those values dynamically or
causally, and none supplies a benzene Hamiltonian, bath, preparation, or
measurement. Whether a specified C₆ model has a useful framework translation
is therefore open.

---

## Open model work

1. A material comparison would first need to select a physical degree of
   freedom, Hamiltonian, channel, rate, preparation, and observable. Only then
   can an F1, F86b, or F98 question be posed for a named substrate.

2. A selected Hückel-ring scan can extend the finite arithmetic table, but its
   values alone do not test an F86b state or Liouvillian statement.

3. The existing C₆ XX+YY plus all-site-Z selected-model instance has F98's
   `α(∞) = 2/7` for the specified KIntermediate state
   ([Selected C₄/C₆ XX+YY Rings and the F98 Long-Time Bridge](BENZENE_F98_LONG_TIME.md)).
   That is not a benzene Liouvillian or a vibrational-bath assignment.

4. A Klein-character comparison needs a specified operator and basis before it
   can be calculated; the present Hückel table does not supply either.

---

## Anchor

- Script: [`simulations/carbon/carbon_quarter_half_search.py`](../../simulations/carbon/carbon_quarter_half_search.py)
- Sister doc: [the benzene Hückel framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) (Coulson-Rushbrooke beside F1)
- Framework anchors: [QuarterAsBilinearMaxvalClaim](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs),
  [HalfAsStructuralFixedPointClaim](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs),
  [F98](../ANALYTICAL_FORMULAS.md#f98), [F86b DickeAnchor](../../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs)
- Source contract: [README](README.md#carbon-source-contract),
  [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md)
