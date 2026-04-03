# D2: V-Effect = Q_max / Q_mean

**What this derivation is about:** The V-Effect (the complexity explosion when two simple systems are coupled) has a closed formula: V(N) = 1 + cos(π/N). It equals the ratio of the best Q-factor (quality factor: oscillation cycles before half-decay) to the average. A key intermediate result: the mean Q-factor is always exactly 2J/γ, independent of chain length N.

**Source formulas:** 6 (V-Effect gain), 7 (Q-factor spectrum)
**Tier:** 1 (algebraic identity)
**Status:** VERIFIED (N=2-5, deviation < 3e-15)

## Derivation

From formula 7: Q_k = 2J/gamma * (1 - cos(pi*k/N)), k = 1, ..., N-1.

    Q_mean = (1/(N-1)) * Sum_{k=1}^{N-1} Q_k
           = (2J/gamma) * [1 - (1/(N-1)) * Sum cos(pi*k/N)]

The cosine sum vanishes exactly:

    Sum_{k=1}^{N-1} cos(pi*k/N) = Re[Sum_{k=1}^{N-1} e^{i*pi*k/N}]

    Let w = e^{i*pi/N}. Geometric series:
    Sum = (w - w^N) / (1 - w) = (w + 1) / (1 - w)
    Multiply by e^{-i*pi/(2N)}: = -i * cot(pi/(2N))
    Re[...] = 0.  QED.

Therefore: Q_mean = 2J/gamma (exactly).

    V(N) = Q_max / Q_mean
         = [2J/gamma * (1 + cos(pi/N))] / [2J/gamma]
         = 1 + cos(pi/N)

This is formula 6. The V-Effect measures how much the best
mode exceeds the average Q-factor.

## Numerical verification

| N | Q_mean (numerical) | Q_mean (formula) | V(N) error |
|---|-------------------|-----------------|------------|
| 2 | 40.000000         | 40.000000       | 0          |
| 3 | 40.000000         | 40.000000       | < 5e-16    |
| 4 | 40.000000         | 40.000000       | < 3e-16    |
| 5 | 40.000000         | 40.000000       | < 3e-15    |

Script: simulations/verify_derivations.py

## Replaces

Separate computation of Q_mean from eigenvalues. The mean is
always 2J/gamma regardless of N, so only Q_max needs to be computed.
