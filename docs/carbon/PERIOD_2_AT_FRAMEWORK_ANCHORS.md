# Period 2 + 3 Table Fractions Compared with Framework Polarity Anchors

**Date:** 2026-05-17 night
**Authors:** Tom + Claude
**Status:** exact `Fraction` grouping for the script's curated input table.
The comparison is not a material, biological, or carbon-model derivation.
**Script:** [`simulations/carbon/period_2_at_framework_anchors.py`](../../simulations/carbon/period_2_at_framework_anchors.py)
**Scope:** framework anchors retain their stated formal tiers. The table rows
are labels supplied to the script; no atomic state, bath, spectrum, dynamics,
or observational map is constructed here.

---

## What got tested

The R=CΨ² framework names four polarity-fraction anchors on the qubit dyadic
ladder, all proven (or freshly derived) in 2026:

| Fraction | Framework anchor | Source |
|----------|------------------|--------|
| **1/4** | `QuarterAsBilinearMaxval` | maxval of p·(1−p) = (1/2)²; Mandelbrot cardioid maxval; F97 |
| **3/8** | F86b KIntermediate Dicke α_total | **derived this morning (2026-05-17, commit b9ba5f6)** via X⊗N-eigenbasis decomp |
| **1/2** | `HalfAsStructuralFixedPoint` | argmax of p·(1−p); polarity-pair {−1/2, +1/2}; d=2 selector |
| **3/4** | polarity-complement of Quarter | 1 − 1/4; F86b α_total at γ² > 1 boundary; F88b popcount-mirror complement |

The script contains a curated H--Ar table with a chosen occupation/slot
normalization and computes each supplied ratio as a `fractions.Fraction`. The
question is strictly bookkeeping: which entries of that input table equal the
listed framework fractions? The script is not an independent source for the
atomic assignments and does not construct a material Hamiltonian or channel.

---

## Result

**Within that curated input table,** all four selected framework fractions
occur among the period-2 and period-3 labels. The equality is exact arithmetic
on the supplied fractions; it is not a physical realization claim.

| Framework fraction | Period-2 input labels | Period-3 input labels |
|--------------------|-----------------------|-----------------------|
| 1/4 (Quarter) | Be (2/8) | Mg (2/8) |
| 3/8 (KIntermediate) | B (3/8) | Al (3/8) |
| 1/2 (Half) | C (4/8) | Si (4/8) |
| 3/4 (polarity complement) | O (6/8) | S (6/8) |
| 1 (endpoint) | Ne (8/8) | Ar (8/8) |
| other input fractions | Li, N, F | Na, P, Cl |

The period-1 H input row is likewise 1/2 under its separately supplied
two-slot normalization. This does not identify that row with a qubit or a
framework state.

---

## The listed CHNOPS subset

The script's CHNOPS subset contains the following selected ratios:

| Input label | Fraction | Relation to selected set |
|-------------|----------|--------------------------|
| H | 1/2 | Half |
| C | 1/2 | Half |
| N | 5/8 | not in the initial four-fraction set |
| O | 3/4 | polarity complement |
| P | 5/8 | not in the initial four-fraction set |
| S | 3/4 | polarity complement |

The resulting count is four of six entries in this chosen subset. It does not
partition chemical or biological function, and it does not establish a
correlation with molecular structure. Any such claim requires an independently
specified molecular dataset, comparison rule, and a test that can fail.

---

## The 3/8 table coincidence

The 3/8 anchor was derived this morning (2026-05-17, commit `b9ba5f6`) via
X⊗N-eigenbasis decomposition for the Dicke superposition
`(|D_{N/2−1}⟩ + |D_{N/2}⟩)/√2`, formalised in
[`DickeAnchor.cs`](../../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs).
The closed form: α_total = (1 − γ²)/2 with γ = ⟨ψ|X⊗N|ψ⟩; the KIntermediate
case has γ = 1/2, giving α_total = 3/8.

The supplied B row has the rational label 3/8 under the script's selected
eight-slot normalization. That equals the KIntermediate value, but it does not
identify an atomic shell with the Dicke state or imply a common polarity axis.
The equality is a candidate for a future model comparison only after its
degrees of freedom, Hamiltonian, environment, preparation, and readout have
been named.

---

## Scope of the comparison

There are two distinct objects here: F86b/F99's qubit-state algebra and a
curated table of valence-shell fractions. The present calculation only compares
numbers from them. It does not supply a carbon or biological Hamiltonian,
dissipator, spectral prediction, dynamical law, or experimental observable.

The material mapping remains open. A bridge has to name the physical degree of
freedom, coupling convention, bath channel and rate, preparation, and
measurement map. In particular, neither an element symbol nor a matching
fraction fixes `Q = J/γ`; see [the carbon source contract](README.md#carbon-source-contract)
and [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md).

---

## Other fractions in the input table

The remaining period-1--3 input rows fall at 1/8, 5/8, and 7/8. Relative to
the initially selected anchor set, each differs by 1/8 from at least one listed
fraction. This is an exact statement about the finite input table and the
chosen set; it supplies no biological, chemical, spectral, or dynamical
classification.

---

## Anchor

- Script: [`simulations/carbon/period_2_at_framework_anchors.py`](../../simulations/carbon/period_2_at_framework_anchors.py)
- Sibling docs: [1/4 and 1/2 in carbon](QUARTER_HALF_IN_CARBON.md) (carbon-specific 1/2 + 1/4 layers),
  [the benzene Hückel framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) (Coulson-Rushbrooke beside F1)
- Framework anchors: [F98](../ANALYTICAL_FORMULAS.md#f98) (today's KIntermediate bridge),
  [DickeAnchor.cs](../../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs) (today's 3/8 derivation),
  [Pi2KnowledgeBaseClaims.cs](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)
  (HalfAsStructuralFixedPoint, QuarterAsBilinearMaxval)
- Scope: [carbon source contract](README.md#carbon-source-contract),
  [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md)
