# N=5 Check: Special or Selection-Biased?

**Status:** Complete finite comparison. No global structural verdict follows.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track C)

---

## Motivation

N=5 appears repeatedly in this repo: IBM Torino had 5-qubit chains, CONCENTRATOR_GEOMETRY centers on N=5, the V-Effect analysis used N=5, the 14-fold degeneracy appears at N=5. Tom's observation: "N=5 hat sich schon oefters als perfekt gezeigt." Is this pattern recognition finding real structure, or selection bias from repeated use?

## The table (N=3-8, uniform chain, γ = 0.1)

| N | d² | Max mult | Frac distinct | Slow rate | Slow/Σγ | Max sector dim |
|---|-----|----------|--------------|-----------|-----------------|----------------|
| 3 | 64 | 6 | 0.406 | 0.266 | 0.887 | 9 |
| 4 | 256 | 14 | 0.4961 | 0.299 | 0.748 | 36 |
| **5** | **1,024** | **14** | **0.4766** | **0.319** | **0.637** | **100** |
| 6 | 4,096 | 19 | 0.5388 | 0.332 | 0.553 | 400 |
| 7 | 16,384 | 22 | 0.497 | 0.340 | 0.485 | 1,225 |
| 8 | 65,536 | n/a | n/a | 0.347 | 0.433 | 4,900 |

N=8 max multiplicity and frac distinct are not available (would require full 65536×65536 eigendecomposition, 73 GB RAM, C# engine only). The slow-mode rate is from the SE-restricted Liouvillian (64×64, trivial). A sixth metric (palindromic pair count) is omitted: it is 100% for all N by the [mirror symmetry theorem](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

## Analysis: is N=5 extremal on any axis?

**Max multiplicity {6, 14, 14, 19, 22}:** Monotonically non-decreasing. N=5 ties with N=4, not a peak.

**Fraction distinct {0.406, 0.4961, 0.4766, 0.5388, 0.497}:** N=5 is a
strict local minimum between N=4 and N=6 for this metric. It is not the minimum
of the displayed N=3..7 window, because N=3 is lower. This is a real local
feature and not a global size-selection result.

**Slow-mode rate {0.266, 0.299, 0.319, 0.332, 0.340, 0.347}:** Monotonically increasing, approaching an asymptotic limit. N=5 is not an inflection point (the rate of increase slows smoothly).

**Slow rate / Σγ {0.887, 0.748, 0.637, 0.553, 0.485, 0.433}:** Monotonically decreasing. The slow mode occupies a shrinking fraction of the total dephasing budget as N grows. No N=5 feature.

**Max sector dimension {9, 36, 100, 400, 1225, 4900}:** Grows as C(N, ⌊N/2⌋)², approximately exponential. No N=5 feature.

## Verdict

Four of the five displayed metrics do not select N=5 inside this window. The
distinct-frequency fraction does: `0.4961 -> 0.4766 -> 0.5388` makes N=5 a
strict local minimum. The table does not establish why that feature occurs,
whether it continues, or whether some other predeclared metric would select a
different N. Hardware history explains why N=5 was sampled often, but it is not
a proof that every N=5 feature is an artifact.

---

## Files

- `simulations/three_values.py` (Track C: N-scaling computation)
- `simulations/results/values_investigations/three_values_results.json` (raw data)
- [Concentrator Geometry](CONCENTRATOR_GEOMETRY.md) (SE slow-mode rate data)
- [Symmetry Census](SYMMETRY_CENSUS.md) (sector dimension formula)

---

*April 12, 2026; current reading September 15, 2026. The displayed window
contains one strict local N=5 minimum and no licensed global verdict.*
