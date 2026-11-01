# D5: Dynamic Palindromic Mode Count

**Source formulas:** 4 (Stat(N)), 22 (XOR drain), 23 (XOR vanishing)
**Tier:** 1 (combinatorial)
**Status:** VERIFIED (N=2-5, Heisenberg chain, the γ = 0 count)

## What this derivation is about

Of the 4^N modes in the Liouvillian (the superoperator governing
open quantum evolution), how many oscillate and how many carry no
frequency? This derivation gives the count at γ = 0 and verifies it
for N=2-5, together with what γ > 0 changes.

---

## Derivation

Total Liouvillian modes: 4^N (dimension of superoperator space).

Two categories at γ = 0:
1. Zero frequency: the commutant of H, Stat(N) = Sum_J m(J,N) * (2J+1)^2
   (F4), where J is the total angular momentum quantum number and m(J,N)
   is the multiplicity (how many times spin J appears in N coupled qubits)
2. Oscillating: everything else

    Oscillating = 4^N - Stat(N)

The N+1 XOR modes at rate 2·N·γ (F22) commute with H, so they are
counted inside Stat(N), not beside it.

At γ > 0 the count is not kept: the dissipator splits part of the
degenerate commutant, and that part acquires a frequency: at N=4, 8 modes, with |Im| growing as γ² (1.5e-5 at
γ = 0.01, 3.6e-4 at γ = 0.05); at N=5, 24 modes, some growing as γ² and
some as γ⁴, the smallest |Im| 1.9e-6 at γ = 0.05. So at γ = 0.05 and 0.1
the oscillating count is 210 at N=4 and 928 at N=5, while N=2 and N=3
keep 6 and 40. The count 4^N - Stat(N) belongs to γ = 0 itself, and it
is read there exactly.

## Numerical verification

Heisenberg chain, J = 1. The γ = 0 column is exact: in the Pauli-string
basis L_H is an integer matrix, and its nullity over a prime bounds the
rational nullity from above, while the SU(2) symmetry of H bounds the
commutant from below by Stat(N); at N = 2 to 5 the two meet, over two
primes. The γ = 0.05 column is read from the eigensolver: the largest
spurious |Im| on the real side is at most eps·‖L‖, frequencies are
counted above a threshold of 1e3·eps·‖L‖, and none may lie within a
factor 1e3 above it. The count is the same at γ = 0.1.

| N | 4^N  | Stat(N) = nullity of L_H | XOR (inside Stat) | Oscillating, γ = 0  | Oscillating, γ = 0.05 |
|---|------|--------------------------|-------------------|--------------------|-----------------------|
| 2 | 16   | 10                       | 3                 | 6                  | 6                     |
| 3 | 64   | 24                       | 4                 | 40                 | 40                    |
| 4 | 256  | 54                       | 5                 | 202                | 210                   |
| 5 | 1024 | 120                      | 6                 | 904                | 928                   |

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Eigenvalue counting and classification. The oscillating fraction
approaches 1 exponentially as N grows (from F23).
