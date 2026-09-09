# Listening to the whole measurement record

**Date:** 2026-09-09

## What this is about

Our measurement already returns four bits, but we had been reducing them to
one number: +1, −1, or zero. Many different patterns became the same zero.
Signal processing suggests keeping the patterns and asking how much each
one helps distinguish the two settings of the chain. With the same
preparation and measurements, that change makes the distinction cheaper.
Most of the gain is already present in the three spectator spins' bit
patterns. We were throwing away useful information before analyzing it.

## Abstract

For the fixed product preparation and local measurement selected in the
[readout search](ROUTE_B_N4_READOUT_SEARCH.md), covariance-weighted outcome
scores reduce the ideal mean-SNR-5 cost at the N = 4 EP from 88,863–89,947
to 17,570–18,497 shots per parameter setting. Optimal weighting of the
three pooled outcomes already reaches 22,582–23,930; the eight spectator
patterns alone reach 18,243–19,211. The preparation and measurement basis
are unchanged. Times are optimized on the same 401-point grid, with a
fixed-time comparison also retained. This is an optimum over linear
scores for known outcome distributions, not finite-shot test power,
quantum Fisher information, or an EP-specific sensing advantage.

## The borrowed tool and its exact mapping

The signal-processing tool is the covariance-weighted matched filter:
weight a known signal difference against its fluctuations. The
[NPTEL generalized matched-filter lecture](https://archive.nptel.ac.in/content/storage2/courses/117103018/module10/lec36/1.html)
describes inverse-covariance weighting for Gaussian detection. Here the
data are categorical counts. We transfer the mean-SNR optimization and
derive it for their exact covariance; no Gaussian detection theorem is
assumed for finite-shot error probabilities.

Let p and r be the outcome probabilities at two known parameter settings,
d = r − p, s = p + r. Two independent empirical histograms, each from M
shots, have difference covariance C/M, where

    C = diag(s) − ppᵀ − rrᵀ.

A real outcome score w has squared mean SNR

    M (wᵀd)² / (wᵀCw).

Define g = Σ dᵢ²/sᵢ, omitting empty bins. The choice wᵢ = dᵢ/sᵢ obeys
sᵀw = 0 and Cw = (1 − g/2)d. Cauchy–Schwarz in the covariance metric
therefore gives the maximum squared SNR per M as g/(1 − g/2).
Adding a constant to all scores or rescaling them changes no SNR.
Identical distributions give zero; disjoint supports give zero-variance
separation. All reported comparisons have overlapping supports.

This is also a statistical information-loss test: merging bins restricts
the possible scores, so it cannot improve the optimum. Full outcomes must
do at least as well as either the ternary or spectator-only record.

## Protocol and results

Keep N = 4, XY, local Z dephasing γ = 1, q = 2, end weights 1 + ε/2.
Prepare (|0100⟩ + i|0110⟩)/√2, measure X on site two and Z on the others
(site zero is the rightmost bit). The 16 outcomes are eight spectator
patterns times the two X signs. The earlier readout keeps the signs for
spectator mask 2 and pools every other pattern into zero.

Compare ε_center with ε_center ±0.05. Each row uses two independent
datasets of M shots, costing 2M total. At ε_center = √2 − 2:

| Processing | M each, −0.05 | M each, +0.05 |
|---|---:|---:|
| Original +1/−1/0 mean | 88,863 | 89,947 |
| Best score of three outcomes | 22,582 | 23,930 |
| Best score of eight spectator patterns | 18,243 | 19,211 |
| Best score of all 16 outcomes | 17,570 | 18,497 |

Each method chooses its time in 0 ≤ γt ≤ 2, spacing 0.005. The full record
chooses 0.385/0.375. At the original mean's times, 0.345/0.340, it still
needs only 17,994/18,800 shots per setting. Thus the gain does not depend
on changing the physical measurement time. All records are obtainable
from the same measurement shots; different processing needs no new basis.

At centers shifted by −0.2 and +0.2 from the EP, the full-record costs
are respectively 14,138–14,953 and 21,407–22,408. The benefit persists
away from the EP. The small gap between spectator-only and full-record
costs shows that most of this contrast is already accessible from the
spectator populations. It does not identify a Jordan signature.

Weights use the two model distributions and are chosen before collecting
data. Calibration, learning weights from data, nuisance parameters, SPAM
and T1 errors are outside this ideal calculation. A likelihood-ratio test
would answer a different question about error probabilities; it is not
implemented or claimed here.

## Reproduction

Run `python simulations/route_b_n4_histogram_filter.py` with
OPENBLAS_NUM_THREADS=1. Its JSON result is
`simulations/results/route_b_n4_histogram_filter.json`.
Full-density exponential action supplies outcome probabilities, checked
independently by rotation into the measurement basis. The score formula is
checked against the covariance pseudoinverse and direct score moments.
The probability-normalization null direction is excluded from numerical
inversion. Every grid point checks the information hierarchy. Controls
include identical distributions, disjoint supports and contrast destroyed
by pooling outcomes. The original mean reproduces the prior experiment.
