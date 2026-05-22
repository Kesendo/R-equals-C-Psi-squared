# Periodic Palindrome vs V-Effect: a Correspondence That Does Not Hold

**Date:** 2026-05-22
**Status:** Negative result (Tier 1 computational). The proposed correspondence
between the periodic-table palindrome's deviations and the V-Effect's
boundary-sector breaking does not survive a direct test. Documented per the
repo's negative-results convention.
**Script:** [`simulations/periodic_palindrome_veffect.py`](../../simulations/periodic_palindrome_veffect.py)
**Answers:** [README](README.md) open question 6.

---

## The question

Two framework-adjacent palindromes were proposed to be the same phenomenon.

[`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py)
establishes the **periodic palindrome**: across a period, per-element properties
(first ionization energy, Pauling and Allen electronegativity) form pairs whose
sums sit near a constant, value_k + value_{N+1−k} ≈ const, a pattern
periodic_palindrome.py finds statistically significant by a shuffle null test
(p-values down to 10⁻⁴).

[V_EFFECT_BOUNDARY_LOCALIZATION](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md)
establishes the **V-Effect boundary localization**: the framework's Liouvillian
palindrome break is confined to boundary sectors 0 < w < N (w = XY-weight); the
extreme sectors w = 0 and w = N are immune, exactly, to 10⁻¹⁵.

Question 6 asks whether these match: do the periodic palindrome's deviations sit
where the V-Effect predicts breaking? The proposed map sends closed shells (the
period ends) to the immune extreme sectors and partly-filled shells (the period
middle) to the breaking boundary sectors.

## The test

The aggregate coefficient of variation that periodic_palindrome.py reports hides
where in a period the palindrome deviates. The per-pair deviation does not: each
pair sum's distance from the period mean, with pair k = 0 the outermost pair (the
period's two ends) and increasing k moving toward the partly-filled middle. The
V-Effect correspondence makes a sharp prediction: the deviation should
concentrate in the inner pairs and leave the outer pairs near-palindromic.

## Result: the localization prediction fails

Across all 11 period × property combinations (IE periods 2-6, Pauling EN
periods 2-5, Allen EN periods 2-3), the inner-pair mean deviation and the
outer-pair mean stay the same order of magnitude. The inner/outer ratio ranges
from 0.69 (IE period 6, where the outer pairs deviate *more*) to 2.5 (EN
period 4), clustering near 1. There is nothing resembling the V-Effect's own
dichotomy, where boundary residuals (≥ 8) and extreme residuals (≤ 10⁻¹⁵) are
separated by ten orders of magnitude.

What the deviations *do* track is specific anomalous electron configurations. The
single most-deviating pair, per IE period:

| Period | Dominant pair | Configurations | Deviation |
|--------|---------------|----------------|-----------|
| 2 | (B, O) | 2p¹ and 2p⁴ | 13.6% |
| 3 | (Al, S) | 3p¹ and 3p⁴ | 14.5% |
| 4 | (Cr, Ga) | anomalous 3d⁵4s¹, and 4p¹ | 20.6% |
| 5 | (Mo, In) | anomalous 4d⁵5s¹, and 5p¹ | 15.4% |

These are the textbook ionization-energy anomalies: the s²p¹ effect (group 13,
the p electron above a filled s² is loosely held), the post-half-filled pairing
penalty (group 16, the fourth p electron pays exchange energy), and the
half-filled-d stability of chromium and molybdenum. The partly-filled pairs that
are *not* anomalous, for instance (C, N) = 2p²/2p³, are among the most
palindromic (1.7% deviation at period 2).

## The answer to Question 6

**No.** The periodic palindrome's deviations are not the V-Effect's boundary
break. The deviations are real, but they sit at specific anomalous configurations
and are accounted for by standard atomic physics (Hund's rule, exchange-pairing
energy, the s²p¹ effect). The V-Effect's boundary-sector mechanism cannot predict
them: its XY-weight sectors carry no information distinguishing an anomalous
configuration (p¹, p⁴, d⁵) from a regular one (p², p³), and its sharp
extreme-immune / boundary-break dichotomy has no analog in the periodic data,
where deviation is spread across inner and outer pairs within a factor of two.

The two palindromes are different objects. The V-Effect's w is a label on Pauli
operators in a 4^N Liouville space; a shell-filling count is a number of
electrons. Equating them is a category slip: both structures have "a palindrome
that weakens somewhere", but the *where* and the *mechanism* do not correspond.
The periodic palindrome remains a genuine and statistically significant pattern;
its deviations are atomic chemistry, not an inherited V-Effect signature.

## Anchor

- Script: [`simulations/periodic_palindrome_veffect.py`](../../simulations/periodic_palindrome_veffect.py)
  (per-pair deviation analysis), built on
  [`simulations/periodic_palindrome.py`](../../simulations/periodic_palindrome.py)
  (the palindrome test itself)
- Compared against: [V_EFFECT_BOUNDARY_LOCALIZATION](../../experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md),
  [V_EFFECT_PALINDROME](../../experiments/V_EFFECT_PALINDROME.md)
- Parent: [README.md](README.md)
