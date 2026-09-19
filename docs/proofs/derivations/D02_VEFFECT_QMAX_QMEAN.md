# D2: F6 Q-edge gain = Q_max / Q_mean

**What this derivation is about:** The **F6 Q-edge gain** has the closed formula V(N) = 1 + cos(π/N). It is the within-N ratio of the best Q-factor to the average Q-factor. “V-Effect gain” is retained only as a qualified historical alias; this exact ratio is separate from the finite **V-Effect census**. A key intermediate result is Q_mean = 2J/γ in the declared single-excitation model.

**Source formulas:** 6 (F6 Q-edge gain; historical alias: V-Effect gain), 7 (Q-factor spectrum)
**Tier:** 1 (algebraic identity)
**Status:** VERIFIED (N=2-5, deviation < 3e-15)

## Derivation

From F7: Q_k = 2J/gamma * (1 - cos(pi*k/N)), k = 1, ..., N-1.

    Q_mean = (1/(N-1)) * Sum_{k=1}^{N-1} Q_k
           = (2J/gamma) * [1 - (1/(N-1)) * Sum cos(pi*k/N)]

The cosine sum vanishes exactly:

    Sum_{k=1}^{N-1} cos(pi*k/N) = Re[Sum_{k=1}^{N-1} e^{i*pi*k/N}]

    Let w = e^{i*pi/N}. Geometric series:
    Sum = (w - w^N) / (1 - w) = (w + 1) / (1 - w)
    Multiply by e^{-i*pi/(2N)}: = +i * cot(pi/(2N))
    Re[...] = 0.  QED.

Therefore: Q_mean = 2J/gamma (exactly).

    V(N) = Q_max / Q_mean
         = [2J/gamma * (1 + cos(pi/N))] / [2J/gamma]
         = 1 + cos(pi/N)

This is F6. It measures how much the best Q-factor exceeds the average
within the same N and model; it is not the V-Effect census or a coupling-causality claim.

## Exact and numerical verification

For `w=cos(x)+i sin(x)`, the verifier simplifies the denominator-cleared
residual `(w+1)-i*cot(x/2)*(1-w)` to exact zero in SymPy. Flipping the sign
leaves `2+2i` at N=2. This sign-sensitive identity is independent of the
floating-point cosine-sum and ratio checks below.

| N | Q_mean (numerical) | Q_mean (formula) | V(N) error |
|---|-------------------|-----------------|------------|
| 2 | 40.000000         | 40.000000       | 0          |
| 3 | 40.000000         | 40.000000       | < 5e-16    |
| 4 | 40.000000         | 40.000000       | < 3e-16    |
| 5 | 40.000000         | 40.000000       | < 3e-15    |

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Separate computation of Q_mean from eigenvalues. The mean is
always 2J/gamma regardless of N, so only Q_max needs to be computed.
