# D5: Dynamic Palindromic Mode Count

**Source formulas:** 4 (Stat(N)), 22 (XOR drain), 23 (XOR vanishing)
**Tier:** 1 (combinatorial)
**Status:** VERIFIED (N=2-5, exact match at gamma -> 0)

## What this derivation is about

Of the 4^N modes in the Liouvillian (the superoperator governing
open quantum evolution), how many are oscillating palindromic pairs,
how many are stationary, and how many drain through the XOR channel?
This derivation gives the exact count in each category and verifies
it numerically for N=2-5.

---

## Derivation

Total Liouvillian modes: 4^N (dimension of superoperator space).

Three non-overlapping categories at gamma = 0:
1. Stationary: Stat(N) = Sum_J m(J,N) * (2J+1)^2 (formula 4),
   where J is the total angular momentum quantum number and m(J,N) is
   the multiplicity (how many times spin J appears in N coupled qubits)
2. XOR drain: N+1 modes at rate 2·N·γ (formula 22)
3. Oscillating palindromic: everything else

    Oscillating = 4^N - (N+1) - Stat(N)

**Important caveat:** This decomposition is exact only at gamma -> 0.
At finite gamma, the Hamiltonian mixes weight sectors of the same
parity (w with w +/- 2). This mixing redistributes modes between
rate groups but does not change the total count of oscillating
vs stationary vs XOR modes; it only shifts their rates.

## Numerical verification (gamma = 1e-6)

| N | 4^N  | Stat(N) | XOR | Oscillating | Match? |
|---|------|---------|-----|-------------|--------|
| 2 | 16   | 10      | 3   | 3           | YES    |
| 3 | 64   | 24      | 4   | 36          | YES    |
| 4 | 256  | 54      | 5   | 197         | YES    |
| 5 | 1024 | 120     | 6   | 898         | YES    |

Script: simulations/verify_derivations.py

## Replaces

Eigenvalue counting and classification. The oscillating fraction
approaches 1 exponentially as N grows (from formula 23).
