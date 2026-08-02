# D7: Q-Factor Distribution (Arcsine)

**What this derivation is about:** The Q-factor (quality factor: how many oscillations a mode completes before it decays to half amplitude) of the modes F7 describes follows an arcsine distribution: most modes cluster at the extremes (very high Q or very low Q), with few in the middle. Those modes are the N-1 oscillating ones of the chain's (0,1) coherence block, the N-dimensional span of the |0><j| between the ferromagnet and the single excitations; this derivation says nothing about the rest of the Liouvillian's palindromic modes. It is derived analytically from the dispersion relation, and the distribution's finite-N statistics are checked against their closed forms for N=5 through N=50.

**Source formulas:** 2 ((0,1) coherence block dispersion), 7 (Q-factor spectrum)
**Tier:** 1-2 (algebraic, converges with N)
**Status:** VERIFIED (N=5-50) in the following sense, and it is worth being exact about which step each part rests on. Shape is arcsine; endpoints and mean are exact at all N; the finite-N variance is closed exactly as A^2(N-2)/(2(N-1)), with the arcsine value A^2 cos^2(pi/N)/2 its N -> infinity limit. All of that is algebra on top of F7's closed form, and the table below compares two closed forms to each other, so it cannot and does not test the SCOPE claim (which modes these are). What pins F7's Q_k to actual Liouvillian eigenvalues is D10 and F2: a full eigendecomposition at N=2-6 and block closure at N=2-10.
**Origin:** Cascade A (derivation_cascades.py)

## Derivation

From F7: Q_k = 2J/gamma * (1 - cos(pi*k/N)), k = 1, ..., N-1.

For large N, treat k/N as continuous variable q in (0,1):

    Q(q) = 2J/gamma * (1 - cos(pi*q))
    dQ/dq = 2J*pi/gamma * sin(pi*q)
    rho(Q) = |dq/dQ| = gamma / (2J*pi*sin(pi*q))

Express sin(pi*q) in terms of Q:

    cos(pi*q) = 1 - Q*gamma/(2J)
    sin(pi*q) = sqrt(1 - cos^2(pi*q))
              = sqrt(Q*gamma/(2J) * (2 - Q*gamma/(2J)))

    rho(Q) = 1 / (pi * sqrt((Q - Q_min) * (Q_max - Q)))

where Q_min = 2J/gamma*(1-cos(pi/N)), Q_max = 2J/gamma*(1+cos(pi/N)).

This is the **arcsine distribution** (Beta(1/2, 1/2) (the standard U-shaped probability distribution that peaks at both endpoints) scaled to [Q_min, Q_max]).

Properties (write A = 2J/gamma):
- U-shaped: modes cluster at Q_min and Q_max (band edges)
- Mean: Q_mean = (Q_min + Q_max)/2 = A (exact at all N, from sum_{k=1}^{N-1} cos(pi*k/N) = 0)
- Variance (continuum/arcsine): (Q_max - Q_min)^2 / 8 = A^2 cos^2(pi/N) / 2
- Variance (exact, finite-N discrete mode set): A^2 (N-2)/(2(N-1)). This is exact at every
  finite N, not a limit: it follows from Var(Q) = A^2 (E[cos^2] - E[cos]^2) with E[cos] = 0 and
  sum_{k=1}^{N-1} cos^2(pi*k/N) = (N-2)/2 (so E[cos^2] = (N-2)/(2(N-1))). Both variances converge
  to A^2/2 as N -> infinity, and they differ at finite N by a computable amount (see below).

## Numerical verification

| N  | Var (exact discrete) | Var (arcsine, N -> inf) | Gap |
|----|---------------------|------------------------|------------|
| 5  | 600.0               | 523.6                  | 14.6%      |
| 10 | 711.1               | 723.6                  | 1.7%       |
| 20 | 757.9               | 780.4                  | 2.9%       |
| 50 | 783.7               | 796.8                  | 1.7%       |

The "Var (exact discrete)" column is reproduced bit-for-bit by the closed form A^2(N-2)/(2(N-1))
(here A = 2J/gamma = 40, so A^2 = 1600): 1600*3/8 = 600.0, 1600*8/18 = 711.1, 1600*18/38 = 757.9,
1600*48/98 = 783.7. So the "Gap" column is NOT a fit residual; it is the exact difference between
two closed forms, the finite-N discrete variance A^2(N-2)/(2(N-1)) and its continuum-limit arcsine
value A^2 cos^2(pi/N)/2, both of which converge to A^2/2 as N -> infinity. The sign of the gap flips
between N=5 (arcsine under-predicts) and N >= 10 (arcsine over-predicts), a finite-N crossover, not a
model failure. Mean is exact at all N (cosine sum identity from D2); the arcsine is the exact shape of
the N -> infinity limit (KS-confirmed in the bulk at every tested N).

Script: [`simulations/derivation_cascades.py`](../../../simulations/derivation_cascades.py)

## Replaces

Numerical Q-factor histogram construction for the (0,1) coherence
block. The distribution is fully characterized by Q_min and Q_max (two
parameters from formulas 2 + 7). The U-shape means most modes are near
the extremes, not near the mean.

**Scope.** Uniform gamma, open Heisenberg chain. Both parent formulas
are uniform-gamma statements, and the reason is structural rather than
a missing measurement: under a gamma PROFILE the block stays exactly
closed, but diag(gamma) stops commuting with the block's Laplacian, so
Q_k is no longer (frequency)/(2 gamma) with a common denominator. The
frequencies move as well
([Concentrator Optics](../../../experiments/CONCENTRATOR_OPTICS.md)
Result 3). So the arcsine shape carries no prediction for
sacrifice-zone design, which is a profile by construction.
