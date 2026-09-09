# Moving the preparation and the readout

**Date:** 2026-09-09

## What this is about

We had found a special meeting point of two fading quantum motions, but our
first way of looking at it needed hundreds of thousands of repetitions.
Perhaps we were listening in the wrong place. This note moves the prepared
phase and the measurement around a four-spin chain. Preparing one middle
spin and reading the other makes a small change in the end bonds much
easier to distinguish. The same improvement works nearby, away from the
special point. The useful discovery here is how to look at the chain.

## Abstract

The [first virtual experiment](ROUTE_B_N4_VIRTUAL_READOUT.md) reads a
coherence where it was prepared. A finite search finds a substantially
cheaper parameter contrast by preparing a phase on one middle spin and
reading the other. The benefit persists away from the EP.
At fixed q = 2 and γ = 1, a search over 24 product preparations, 72 local
readouts and 401 times reduces the ideal mean-SNR-5 cost at the EP from
598,838–721,419 to 88,863–89,947 shots per parameter setting for end biases
±0.05. Two off-EP centers retain comparable performance. These are finite
catalogue results for a difference-of-means criterion, not a general
measurement optimum or evidence of an EP-specific sensing advantage.

## Search and accounting

Keep the same N = 4 XY model, γ = 1, q = 2, equal end weights
1 + ε/2 and time τ = γt. Search 24 product states: choose one excited
spin, a distinct spin in |+⟩ or |+i⟩, and put the other two in |0⟩.
Search 72 readouts: X or Y on any spin, either unconditional (8 choices)
or restricted to one of the eight spectator bit patterns (64 choices).
The latter assigns zero to all other patterns, counting every shot.

For each preparation/readout, use the same time in both parameter settings
and M independent shots at each setting. Maximize the difference-of-means
SNR over 401 times in [0,2]. Compare ε_center with ε_center ± 0.05 at
three centers: ε_EP − 0.2, ε_EP, ε_EP + 0.2, where ε_EP = √2 − 2.
This is a numerical design search. An experimental search would have its
own additional cost. The criterion uses the measured mean and variance,
not the full outcome likelihood or quantum Fisher information.

## Best results in this catalogue and time grid

| Center relative to EP | M each, Δε = −0.05 | M each, Δε = +0.05 |
|---|---:|---:|
| −0.2 | 87,588 | 87,393 |
| 0 | 88,863 | 89,947 |
| +0.2 | 94,580 | 96,542 |

All entries target mean SNR 5, so each comparison costs 2M shots in total.
At the EP, the original protocol required 598,838 and 721,419 respectively.
The reduction is a factor 6.74 and 8.02. The best grid times move from
0.685/0.680 to 0.345/0.340. No experimental errors are included.

One representative winner prepares (|0100⟩ + i|0110⟩)/√2, with site zero
the rightmost bit, and measures X on site two. Retain its ±1 outcome only
when sites (zero, one, three) have bits (0,1,0), otherwise assign zero.
Thus the preparation writes ρ₄₆ while the readout measures 2 Re ρ₂₆.
The reflected protocol has the same performance; floating-point ordering
can select either representative. |+⟩ with the complementary Y quadrature
also supplies phase-related alternatives. This is a product preparation
and one local measurement basis setting, with spectator Z measurements.

For this representative, the selected λ = −4 + 2i Jordan contribution
still has nonzero linear coefficient. Numerical contour evaluation gives
a ≈ −0.0520833333333333i and b ≈ −0.0416666666666666 in
2 Re[exp(λτ)(a + bτ)]. Its |b| is smaller than the original protocol's
0.0833333333333333 despite the improved parameter contrast. The score
does not measure the magnitude of this Jordan coefficient alone.

Keeping the EP-selected preparation and readout fixed at the other centers
recovers their best catalogue scores after changing only the time. Even
keeping the EP time fixed needs only about 89,000–100,000 shots per setting
across these centers. The improvement therefore survives without retuning
the preparation/readout. There is no privileged EP advantage in these
comparisons; the lower-ε center is slightly cheaper when its time is chosen.
This finite comparison does not exclude a different result for a different
measurement family or estimation task.

## Reproduction and checks

Run `python simulations/route_b_n4_readout_search.py` with
OPENBLAS_NUM_THREADS=1. Results are in
`simulations/results/route_b_n4_readout_search.json`, including the selected
protocols, times, contrasts and transferred protocols at fixed EP time.

The full density generator propagates all 24 preparations. Repeated short
matrix exponentials are checked against direct exponentials at two times
for every preparation and parameter setting. All three-outcome probabilities
are checked across the grid. Initial |+i⟩ reads +1 under Y, guarding the
complex conjugation convention. Identical parameter settings have zero
contrast. The original baseline is reproduced against the separate
exponential-action propagation. The selected Jordan contribution is
evaluated independently of the contrast score by contour projection.

The [histogram-filter comparison](ROUTE_B_N4_HISTOGRAM_FILTER.md) tests
information discarded by this mean using optimally weighted outcome
scores under the same mean-SNR criterion. It reuses the same measurement
record and compares full, ternary and spectator-only processing.
