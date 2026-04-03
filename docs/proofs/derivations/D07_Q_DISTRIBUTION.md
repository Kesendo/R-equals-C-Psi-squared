# D7: Q-Factor Distribution (Arcsine)

**What this derivation is about:** The Q-factor (quality factor: how many oscillations a mode completes before it decays to half amplitude) of all palindromic modes follows an arcsine distribution: most modes cluster at the extremes (very high Q or very low Q), with few in the middle. This is derived analytically from the dispersion relation and verified numerically for N=5 through N=50.

**Source formulas:** 2 (w=1 dispersion), 7 (Q-factor spectrum)
**Tier:** 1-2 (algebraic, converges with N)
**Status:** VERIFIED (N=5-50, variance converges to arcsine prediction)
**Origin:** Cascade A (derivation_cascades.py)

## Derivation

From formula 7: Q_k = 2J/gamma * (1 - cos(pi*k/N)), k = 1, ..., N-1.

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

Properties:
- U-shaped: modes cluster at Q_min and Q_max (band edges)
- Mean: Q_mean = (Q_min + Q_max)/2 = 2J/gamma (exact, from D2)
- Variance: (Q_max - Q_min)^2 / 8

## Numerical verification

| N  | Var (numerical) | Var (arcsine) | Rel. error |
|----|----------------|--------------|------------|
| 5  | 600.0          | 523.6        | 14.6%      |
| 10 | 711.1          | 723.6        | 1.7%       |
| 20 | 757.9          | 780.4        | 2.9%       |
| 50 | 783.7          | 796.8        | 1.7%       |

Convergence improves with N (continuous limit). Mean is exact
at all N (cosine sum identity from D2).

Script: simulations/derivation_cascades.py

## Replaces

Numerical Q-factor histogram construction. The distribution is
fully characterized by Q_min and Q_max (two parameters from
formulas 2 + 7). The U-shape means most modes are near the
extremes, not near the mean. This is relevant for sacrifice-zone
design (edge modes dominate both fast and slow Q regimes).
