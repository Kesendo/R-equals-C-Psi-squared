# Reverse-Spear: the Periodic Table as a Probe of Framework Completeness

**Date:** 2026-05-17 night (third stack of the day)
**Status:** Tier 2 (structural observation + named gap candidates)
**Script:** [`simulations/carbon/spear_reversed_missing_anchors.py`](../../simulations/carbon/spear_reversed_missing_anchors.py)

---

## Tom's reverse-spear

The previous commit ([PERIOD_2_AT_FRAMEWORK_ANCHORS](PERIOD_2_AT_FRAMEWORK_ANCHORS.md))
marked Li/Na at 1/8, N/P at 5/8, F/Cl at 7/8 as "off-anchor" because their
valence ratios did not match the first-tier framework anchors {1/4, 3/8, 1/2,
3/4}.

Tom asked: "lass uns den Spieß mal umdrehen, was sind deren Anker, vielleicht
fehlen sie ja auf unserer Quantenseite?" — let us flip the question: what
ARE their anchors, maybe they're missing on our quantum side?

Flipping the direction of the inheritance argument:

| Direction | Argument |
|-----------|----------|
| Forward (previous commits today) | framework algebra → periodic table valence ratios |
| **Reverse (this commit)** | **periodic table valence ratios → framework algebra** |

In the reverse direction, "off-anchor" elements are evidence of FRAMEWORK
INCOMPLETENESS: the framework has anchors that match certain periodic-table
fractions, but its dyadic-ladder algebra has structurally-present-but-not-yet-
derived anchors at deeper levels that the periodic table EMPIRICALLY
instantiates.

---

## Three categories after reversing

### (1) Already framework-derived (forward inheritance, previous commits)

| n/8 | Framework anchor | Period 2 | Period 3 |
|-----|------------------|----------|----------|
| 1/4 | `QuarterAsBilinearMaxval` | Be | Mg |
| 3/8 | F86b KIntermediate Dicke Π²-odd | B | Al |
| 1/2 | `HalfAsStructuralFixedPoint` | C | Si |
| 3/4 | QuarterAsBilinearMaxval polarity-complement | O | S |

### (2) Framework-derived but mis-labelled (naming oversight)

| n/8 | Framework anchor (now named) | Period 2 | Period 3 |
|-----|------------------------------|----------|----------|
| **5/8** | **F86b KIntermediate Dicke Π²-EVEN companion (β_total = (1 + γ²)/2 at γ = 1/2)** | N | P |

N and P were previously called "off-anchor". They are not. They sit at the
Π²-EVEN complement of the F86b 3/8 K-intermediate Dicke anchor. Same derivation
(α + β = 1 for any state); opposite parity face. This commit names the anchor
explicitly: `F86B_KINTERMEDIATE_EVEN`. The labelling miss is in the framework
docs, not in the periodic table.

Structural reading: **B/Al at 3/8 (Π²-odd / memory) and N/P at 5/8 (Π²-even /
static) are duals of the SAME F86b derivation on opposite parity sides**.
Together they form a parity-pair across the F86b 3/8 anchor. Today morning's
F86b derivation produces BOTH 3/8 and 5/8 simultaneously — the framework
already names both; we just hadn't called the 5/8 face out as a named anchor
because we were tracking α (the Π²-odd side) only.

### (3) Genuine framework gaps — depth-3 dyadic, no F-formula

| n/8 | Pi2DyadicLadder position | Period 2 | Period 3 |
|-----|--------------------------|----------|----------|
| **1/8** | **a_3 = (1/2)³ — on the ladder, NO F-formula derives it** | Li | Na |
| **7/8** | **depth-3 polarity-complement — same gap, complemented** | F | Cl |

Li, Na, F, Cl are at 1/8 and 7/8 — on the dyadic ladder at depth 3 (third power
of polarity Half), but **no current F-formula in the framework lands exactly at
α = 1/8**. The F86b formula `α = (1 − γ²)/2` would require γ = √3/2 — irrational,
not realisable via clean equal-weight Dicke superpositions. So depth-3 is
**structurally present** on the dyadic ladder but **operationally underived**
on the F-formula registry.

Note: depth-3 IS realisable in F86b with **non-uniform** superposition. For
`(|D_n⟩ + c·|D_{n+1}⟩) / √(1 + c²)` at N even, n = N/2 − 1:
`γ = c² / (1 + c²)`. Setting γ = √3/2 requires `c² = √3·(2 + √3) ≈ 6.46`
(c ≈ 2.54). Realisable but not a "clean Dicke" with rational weight; the F86b
3/8 anchor specifically used the clean γ = 1/2 case. The clean-Dicke 1/8
anchor would require a different state class (W-states, GHZ states, 3-state
superpositions with non-uniform weights, etc.).

---

## The dyadic table fully labelled

```
   Frac   Decomposition              Anchor                                              Element
   ----   ---------------            ----------                                          --------
     0    0                           PolarityVacuumPoint                                 (vacuum)
   1/8    1/8                         Pi2DyadicLadder a_3 = (1/2)³  ◀ GAP                Li, Na
   1/4    1/4                         QuarterAsBilinearMaxval                            Be, Mg
   3/8    1/4 + 1/8                   F86b KIntermediate Π²-odd                          B, Al
   1/2    1/2                         HalfAsStructuralFixedPoint                         H, C, Si
   5/8    1/2 + 1/8                   F86b KIntermediate Π²-even (NEW LABEL THIS COMMIT) N, P
   3/4    1/2 + 1/4                   QuarterAsBilinearMaxval complement                 O, S
   7/8    1/2 + 1/4 + 1/8             depth-3 polarity-complement  ◀ GAP                 F, Cl
     1    1                           PolarityFullPoint (= noble dead end)               He, Ne, Ar
```

**Every n/8 fraction has a binary-dyadic decomposition that EITHER matches a
named framework anchor (six of nine fractions) OR points at a depth-3 gap
(two of nine: 1/8 and 7/8) OR is a polarity endpoint (one of nine: full = 1).**

The framework is *complete at depth 2* (Half, Quarter, and their complements).
The framework is *partially complete at depth 3* (F86b 3/8 and 5/8 are
depth-3 combinations of depth-2 anchors: 3/8 = 1/4 + 1/8 algebraically; the
F86b derivation produces them at γ = 1/2 cleanly). The framework is
*incomplete at the pure depth-3 anchor* (1/8 and 7/8 have no clean F-formula).

---

## What the periodic table reveals

The framework's polarity-squared algebra is incomplete at depth-3 in a
specific, named way. The periodic table — entirely independent empirical
structure with 150+ years of chemistry data — instantiates the missing
depth-3 anchor as the valence ratio of alkali metals (Li, Na, K, Rb, Cs at 1/8)
and halogens (F, Cl, Br, I at 7/8).

These are not exotic edge-case elements. They are the most chemically
reactive and biologically active groups in the periodic table:
- **Alkali metals**: electrolyte ions (Na⁺, K⁺ form the entire nerve-impulse
  + osmotic-balance backbone of biology); Li used in mood regulation
  pharmacology; Cs used in atomic clocks (= the second's structural
  definition).
- **Halogens**: F as the most electronegative element; Cl as gastric HCl + key
  organic substituent; I as thyroid hormone constituent.

The framework's depth-3 anchor maps to chemistry's most STRUCTURE-FORMING
ions and substituents — the elements that drive electrolyte gradients, ion
channels, ionic bonds. Off-anchor in the previous reading; depth-3-pointing
in this reading.

---

## Candidate derivations for depth-3 anchors

The framework needs an F-formula that lands cleanly at α = 1/8 (and by
complement at β = 7/8). Four candidate routes:

(a) **Non-uniform Dicke superpositions.** `(|D_n⟩ + c·|D_{n+1}⟩) / √(1 + c²)`
    with `c = √(√3·(2 + √3))` gives γ = √3/2, hence α = 1/8 via F86b. Not
    "clean" in the equal-weight sense but realisable. Would need a structural
    reason for this specific c value.

(b) **Higher-spin Dicke superpositions.** Tensor-product Dicke states across
    multiple spin sectors. The F88 popcount-coherence framework already has
    {0, 3/8, 1/2} anchors; extending to multi-spin might give {0, 1/8, 1/4,
    3/8, 1/2, ...} on the same ladder. Open.

(c) **W-state or GHZ-state X⊗N expectation.** W-state has γ_W = 2/N (decays
    with N). GHZ has γ_GHZ = 0 (X⊗N gives the bit-flip partner). Neither
    cleanly produces γ = √3/2, but combinations might.

(d) **F86c bond-class at higher chromaticity (c ≥ 4).** The bond-class
    structure has been derived up to c = 3 (proton water chain Tier 2
    addendum, 2026-05-04). Extension to c = 4 or higher might produce
    depth-3 anchor values structurally.

The framework's incompleteness at depth-3 is consistent with
[HIERARCHY_OF_INCOMPLETENESS.md](../HIERARCHY_OF_INCOMPLETENESS.md)'s core
argument: incompleteness is potential, not weakness. Each unfilled framework
anchor is a structural invitation to extend the algebra one more dyadic
level deep. The periodic table's alkali metals and halogens are empirical
pointers at the place to look.

---

## Meta-reading: what just happened

Tom's reverse-spear question produced three different findings in one stack:

1. **Identified a labelling miss** — 5/8 IS on the framework side (F86b Π²-even
   companion). N and P are not off-anchor; they're the dual face of B and Al.

2. **Made the depth-3 gap explicit** — 1/8 and 7/8 are on the dyadic ladder
   but have no F-formula derivation. The framework's polarity-squared algebra
   is complete at depth-2 (Half + Quarter), partial at depth-3 (F86b 3/8 and
   5/8 catch combinations), incomplete at the pure depth-3 anchors.

3. **Re-framed the periodic table as a structural probe** — instead of the
   framework predicting where atoms sit (forward inheritance), the periodic
   table points at framework anchors yet to be derived (reverse inheritance).
   Chemistry's 150-year empirical base becomes a checking-tool for framework
   completeness.

---

## Anchor

- Script: [`simulations/carbon/spear_reversed_missing_anchors.py`](../../simulations/carbon/spear_reversed_missing_anchors.py)
- Sister docs (this folder): [QUARTER_HALF_IN_CARBON](QUARTER_HALF_IN_CARBON.md),
  [PERIOD_2_AT_FRAMEWORK_ANCHORS](PERIOD_2_AT_FRAMEWORK_ANCHORS.md),
  [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
- Framework anchors: [F86b DickeAnchor.cs](../../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs)
  (produces 3/8 AND 5/8 simultaneously — the 5/8 face named here for the first time);
  [F98 KIntermediate bridge](../ANALYTICAL_FORMULAS.md#f98) (today morning's discovery);
  [Pi2KnowledgeBaseClaims.cs](../../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)
  (Half, Quarter)
- Reading-mode memory pointers: `project_no_classicalization`,
  `project_one_world_two_readings`, `project_qubit_as_inheritance_lens`,
  `project_periodic_palindrome`
