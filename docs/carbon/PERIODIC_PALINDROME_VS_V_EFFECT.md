# The Periodic Palindrome and the V-Effect: A Curated-Input Null Test

**Date:** 2026-05-22
**Authors:** Tom + Claude
**Status:** Tier 1 computational result for the stated per-pair analysis of the
repository's curated periodic-property arrays, inside a Tier 5 hierarchy framing.
The hierarchy is an interpretive frame, not a structural or material theorem
established by this null test. The test reports that the selected Level-0
V-Effect boundary localization is not reproduced by these curated deviations.
**Script:** [`simulations/periodic_palindrome_veffect.py`](../../simulations/periodic_palindrome_veffect.py)
**Answers:** [README](README.md) open question 6.
**See also:** [the hardened periodic palindrome](PERIODIC_PALINDROME_HARDENED.md): a hardened sign-flip-null gate shows the palindrome's *presence* (not just its deviations) is mostly smoothness; what survives the ramp leans anti-F1 at the light elements and significantly-but-ambiguously F1-respecting at the heavy ones.

---

## The question

The selected Level-0 V-Effect calculation has its stated boundary-sector
localization. [The Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md)
offers a Tier 5 framing that compares that model shape with periodic-table
language. This page does not establish that comparison as a physical or
structural identity.

Question 6 asks a narrower, empirical question.
[`periodic_palindrome.py`](../../simulations/periodic_palindrome.py) finds a
*periodic palindrome*: across the repository's curated property arrays, paired
values have sums near a constant under its shuffle null. That null has the
known smooth-ramp limitation documented by the producer and hardened follow-up.
Does the *deviation pattern* of these arrays reproduce the Level-0 V-Effect's
boundary-sector localization? This page tests that limited question only.

## The test

The aggregate coefficient of variation that `periodic_palindrome.py` reports
hides where a curated row deviates. The per-pair statistic does not: it is each
pair sum's distance from the row mean, with `k = 0` the outermost pair and
increasing `k` moving inward. The test asks whether that index pattern resembles
the selected V-Effect boundary pattern. It does not make the index a shell,
sector, or material mapping.

## Result within the curated inputs: the Level-0 localization does not transfer

Across all 11 period × property combinations (IE periods 2-6, Pauling EN
periods 2-5, Allen EN periods 2-3), the inner-pair mean deviation and the
outer-pair mean stay the same order of magnitude. The inner/outer ratio ranges
from 0.69 (IE period 6, where the outer pairs deviate *more*) to 2.5 (EN
period 4), clustering near 1. There is nothing resembling the V-Effect's own
dichotomy, where boundary residuals (≥ 8) and extreme residuals (≤ 10⁻¹⁵) are
separated by ten orders of magnitude.

Within each curated first-ionization-energy array, the producer's per-pair
calculation reports these largest absolute deviations:

| Period | Largest-deviation pair | Deviation |
|--------|------------------------|-----------|
| 2 | (B, O) | 13.6% |
| 3 | (Al, S) | 14.5% |
| 4 | (Cr, Ga) | 20.6% |
| 5 | (Mo, In) | 15.4% |

No atomic configuration labels or atomic-mechanism account is produced by this
calculation. Any explanatory candidate for these curated pair/deviation data is
conditional and needs its own specified source and producer. The selected data
do show that a pair nearer the middle need not have the largest deviation, which
is sufficient for the present null test.

## The answer to Question 6

Two levels, two answers.

**The hierarchy remains a Tier 5 framing.** This analysis neither proves nor
refutes a material or structural relation between periodic-table language and
the V-Effect.

**The selected Level-0 mechanism does not transfer to these arrays.** The
per-pair deviations do not localize to the inner pairs in the way required by
the V-Effect boundary picture. The result is a null test on the producer's
curated inputs, not a theorem about the periodic table, atomic chemistry, or a
material mechanism.

The direct outcome is therefore narrow: the repository has no support here for
identifying the periodic deviations with the V-Effect's XY-weight boundary
break. Other explanatory candidates remain open and require their own sources
and producers.

## Anchor

- Script: [`simulations/periodic_palindrome_veffect.py`](../../simulations/periodic_palindrome_veffect.py)
  (per-pair deviation analysis), built on
  [`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py)
  (the palindrome test itself)
- Hierarchy framing: [the Hierarchy of Incompleteness](../HIERARCHY_OF_INCOMPLETENESS.md)
- The Level-0 V-Effect: [V-Effect boundary localization](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md),
  [the V-Effect palindrome](../../experiments/V_EFFECT_PALINDROME.md)
- Parent: [README.md](README.md)
