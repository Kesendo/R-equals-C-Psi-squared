# D4: Dimensional Factor in Crossing

**Source formulas:** 12 (single-qubit crossing), 25 (Bell+ crossing)
**Tier:** 1 (algebraic from CPsi definition)
**Status:** VERIFIED (exact to machine precision)

## Derivation

CPsi = C * Psi = Tr(rho^2) * L1 / (d-1).

For states where C = (1+f^2)/2 and L1 = f (single coherence decay):

    CPsi = f * (1+f^2) / (2*(d-1))

At crossing CPsi = 1/4:

    f * (1+f^2) = (d-1) / 2

This gives:

| System       | d  | d-1 | Crossing condition |
|-------------|----|----|-------------------|
| Single qubit | 2  | 1   | f(1+f^2) = 1/2   |
| Bell+ 2-qubit| 4  | 3   | f(1+f^2) = 3/2   |
| 3-qubit      | 8  | 7   | f(1+f^2) = 7/2   |

The factor (d-1) = 3 between single-qubit and Bell+ is exact.
It comes from the normalization of the L1 coherence by (d-1).

## Numerical verification

| System | f* (analytical) | CPsi at f* | Target |
|--------|----------------|-----------|--------|
| d=2    | 0.42385        | 0.25000   | 0.25   |
| d=4    | 0.86122        | 0.25000   | 0.25   |

Ratio of crossing conditions: 3/2 / 1/2 = 3 = d-1. EXACT.

Script: simulations/verify_derivations.py

## Replaces

Understanding why larger systems have "easier" crossings
(higher f* at crossing = less decoherence needed).
