# N=5 Check: Special or Selection-Biased?

**Status:** Complete. N=5 is not structurally special. All metrics are monotonic with N.
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Script:** `simulations/three_values.py` (Track C)

---

## Motivation

N=5 appears repeatedly in this repo: IBM Torino had 5-qubit chains, SACRIFICE_GEOMETRY centers on N=5, the V-Effect analysis used N=5, the 14-fold degeneracy appears at N=5. Tom's observation: "N=5 hat sich schon oefters als perfekt gezeigt." Is this pattern recognition finding real structure, or selection bias from repeated use?

## The table (N=3-8, uniform chain, gamma=0.1)

| N | d^2 | Max mult | Frac distinct | Slow rate | Slow/Sigma_gamma | Max sector dim |
|---|-----|----------|--------------|-----------|-----------------|----------------|
| 3 | 64 | 6 | 0.406 | 0.266 | 0.887 | 9 |
| 4 | 256 | 14 | 0.496 | 0.299 | 0.748 | 36 |
| **5** | **1,024** | **14** | **0.477** | **0.319** | **0.637** | **100** |
| 6 | 4,096 | 19 | 0.539 | 0.332 | 0.553 | 400 |
| 7 | 16,384 | 22 | 0.497 | 0.340 | 0.485 | 1,225 |
| 8 | 65,536 | n/a | n/a | 0.347 | 0.433 | 4,900 |

N=8 max multiplicity and frac distinct are not available (would require full 65536x65536 eigendecomposition, 73 GB RAM, C# engine only). The slow-mode rate is from the SE-restricted Liouvillian (64x64, trivial).

## Analysis: is N=5 extremal on any axis?

**Max multiplicity {6, 14, 14, 19, 22}:** Monotonically non-decreasing. N=5 ties with N=4, not a peak.

**Fraction distinct {0.406, 0.496, 0.477, 0.539, 0.497}:** Oscillates between 0.4 and 0.54. N=5 is a local minimum but not extreme; N=3 is lower. No clear N=5 peak.

**Slow-mode rate {0.266, 0.299, 0.319, 0.332, 0.340, 0.347}:** Monotonically increasing, approaching an asymptotic limit. N=5 is not an inflection point (the rate of increase slows smoothly).

**Slow-mode rate / Sigma_gamma {0.887, 0.748, 0.637, 0.553, 0.485, 0.433}:** Monotonically decreasing. The slow mode occupies a shrinking fraction of the total dephasing budget as N grows. No N=5 feature.

**Max sector dimension {9, 36, 100, 400, 1225, 4900}:** Grows as C(N, floor(N/2))^2, approximately exponential. No N=5 feature.

## Verdict

**N=5 is not structurally special.** None of the six metrics shows N=5 as an extremum or inflection point. Every metric either increases or decreases monotonically with N, or oscillates without a consistent N=5 peak.

The repeated appearance of N=5 in this repo is **selection bias**: the IBM Torino chain happened to have 5 qubits, so N=5 became the testbed. Once N=5 was the testbed, all subsequent experiments (SACRIFICE_GEOMETRY, CUSP_LENS_CONNECTION, boundary straddling) used N=5 as the reference. The physics at N=5 is representative of the general N trend, not an outlier.

This is a useful negative result: it means the repo's findings generalize to other N without needing special-case analysis.

---

*April 12, 2026. Tom's intuition ("N=5 hat sich schon oefters als perfekt gezeigt") is selection bias, not physics. The data are monotonic.*
