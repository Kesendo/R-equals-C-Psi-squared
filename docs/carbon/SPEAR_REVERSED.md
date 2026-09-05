# Reverse-Spear: A Conditional Periodic-Fraction Comparison

**Date:** 2026-05-17 night (third stack of the day)
**Authors:** Tom + Claude
**Status:** F99 supplies the named framework result: the five canonical angles
give `{0, 1/8, 1/4, 3/8, 1/2}` for its specified non-uniform-Dicke
construction. The periodic labels in this page are curated producer inputs;
they do not establish a material mapping or a causal mechanism.
**Script:** [`simulations/carbon/spear_reversed_missing_anchors.py`](../../simulations/carbon/spear_reversed_missing_anchors.py)

## Scope

The producer places a curated period-2/3 fraction lookup next to framework
fractions. This is a comparison of numbers and labels, not a reverse inference
from a material system to framework completeness. No physical degree of
freedom, Hamiltonian, bath, preparation, observable, or measurement mapping
is specified for the periodic entries.

F99's non-uniform Dicke construction gives `α = 1/8` at `θ = 30°` with
`γ = √3/2`; its complement is `1 − α = 7/8`. These are formal framework
relations, independent of the curated periodic lookup.

---

## Framework arithmetic beside curated labels

| Fraction | Framework relation | Curated period-2 label | Curated period-3 label |
|----------|--------------------|------------------------|------------------------|
| 1/8 | F99: `α = sin²(30°)/2` | Li | Na |
| 1/4 | `QuarterAsBilinearMaxval` | Be | Mg |
| 3/8 | F86b KIntermediate Π²-odd | B | Al |
| 1/2 | `HalfAsStructuralFixedPoint` | C | Si |
| 5/8 | F86b Π²-even complement, `1 − 3/8` | N | P |
| 3/4 | `1 − 1/4` quarter complement | O | S |
| 7/8 | `1 −` F99's `1/8` | F | Cl |

The F99 and complement arithmetic is the framework result. The final two
columns reproduce the producer's curated labels only; matching fractions does
not identify an atomic state, chemical mechanism, or material realization of
the framework quantity.

---

## Dyadic arithmetic

Every `n/8` has a binary decomposition in `{1/2, 1/4, 1/8}`. F99 now supplies
the previously absent formal `1/8` row; `7/8` is its arithmetic complement.
The decomposition is a statement about fractions, not a periodic-table or
material theorem.

## Open work beyond the named result

F99 closes the `1/8` construction for its specified non-uniform Dicke state.
Whether another state family, decomposition, or Lindblad channel has a distinct
anchor formula is open and requires its own operator, state, channel, and
producer. The curated periodic labels do not select any of those inputs.

---

## Anchor

- Curated lookup producer: [`spear_reversed_missing_anchors.py`](../../simulations/carbon/spear_reversed_missing_anchors.py).
- Framework sources: [F99 Niven completeness](F99_NIVEN_COMPLETENESS.md),
  [the depth-3 anchor derivation](DEPTH_3_ANCHOR_DERIVED.md), and
  [F99](../ANALYTICAL_FORMULAS.md#f99).
- Source contract: [Carbon README](README.md#carbon-source-contract),
  [Q Belongs to No Substance](../Q_BELONGS_TO_NO_SUBSTANCE.md).
