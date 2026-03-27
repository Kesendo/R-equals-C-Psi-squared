# Is the Palindromic Mirror Entangled in Operator Space?

<!-- Keywords: Pi operator entanglement operator space, palindromic conjugation
Schmidt decomposition, non-local symmetry Liouville space, operator Schmidt rank
per-site factorization, product superoperator test, XZ ZY non-local palindrome,
R=CPsi2 Pi operator entanglement -->

**Status:** Computationally verified
**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md),
[Pi as Time Reversal](PI_AS_TIME_REVERSAL.md)

---

## Abstract

The palindromic conjugation operator Pi acts on Liouville space, which has
a natural tensor product structure across qubits. If Pi can be written as
Pi_A x Pi_B (a product of per-site maps), the palindromic symmetry is
LOCAL: each qubit carries its own independent palindromic structure. If
no such factorization exists, the symmetry is genuinely NON-LOCAL: it
couples the qubits at the level of the symmetry itself, not just at the
level of the Hamiltonian.

**Finding:** For all standard condensed-matter Hamiltonians (Heisenberg, Ising,
XX, YY, ZZ, and combinations), the palindromic symmetry is LOCAL (operator
Schmidt rank 1). For exactly 2 of the 36 two-term Pauli coupling combinations
(XZ+YZ and ZX+ZY), the symmetry is genuinely NON-LOCAL. No per-site
factorization exists. The non-locality arises precisely when two incompatible
Pauli operators (X and Y) appear on the SAME qubit site.

**Methodological discovery:** The space of valid Pi operators (all operators
satisfying Pi L Pi^{-1} = -L - 2Sg I) is degenerate. Numerical construction
from eigenvalue pairing always finds entangled Pi operators, even when a
product Pi exists. The correct question is not "is the numerical Pi entangled?"
but "does a product Pi exist?"

---

## 1. The Question

Pi acts on C^{d^2}, where d = 2^N for N qubits. The Liouville space has a
natural tensor product structure:

```
C^{d^2} = C^{d_A^2} x C^{d_B^2}
```

An operator on this space is "operator-entangled" if it cannot be written
as a tensor product of per-site operators. The operator Schmidt decomposition
quantifies this: reshape Pi into a matrix M indexed by
[(alpha, alpha'), (beta, beta')] and compute SVD. The number of nonzero
singular values is the Schmidt rank.

- **Rank 1:** Pi = P_A x P_B. Product operator. Local symmetry.
- **Rank > 1:** Pi is operator-entangled. Non-local symmetry.

---

## 2. Method

For each Hamiltonian:

1. Build the known analytical product Pi (per-site P1 map: I<->X, Y<->Z
   with phases) tensored across sites.
2. Test whether it satisfies the conjugation equation
   Pi L Pi^{-1} = -L - 2 Sg I.
3. If yes: compute its Schmidt rank (expected: 1). The symmetry is local.
4. If no: no per-site Pi exists. Build Pi numerically from eigenvalue
   pairing. Report its Schmidt rank (upper bound on minimum). The symmetry
   is non-local.

---

## 3. Results

### Product Pi verified (LOCAL symmetry, Schmidt rank 1)

| Hamiltonian | N | Pi family | Conjugation error | Rank |
|-------------|---|-----------|-------------------|------|
| XX+YY+ZZ (Heisenberg) | 2 | P1 | 1.1e-17 | 1 |
| ZZ (Ising) | 2 | P1 | 2.0e-17 | 1 |
| XX | 2 | P1 | 2.0e-17 | 1 |
| XX+YY+ZZ chain | 3 | P1 | 2.1e-17 | 1 |
| ZZ chain | 3 | P1 | 3.6e-17 | 1 |

For all standard Hamiltonians, the product P1^{xN} satisfies the
conjugation equation to machine precision. Schmidt rank is exactly 1.
**The palindromic symmetry is local.**

### Non-uniform product Pi (LOCAL, rank 1, different maps per site)

| Hamiltonian | N | Uniform P1 error | Product Pi | Rank |
|-------------|---|-------------------|------------|------|
| XY | 2 | 1.41 | P1 x M2 (alternating) | 1 |
| XZ+ZY | 2 | 1.73 | P4 x P1 | 1 |
| YZ+ZX | 2 | 1.73 | P1 x P4 | 1 |

Uniform P1 x P1 fails, but product Pi with DIFFERENT per-site maps
works. Each site chooses its Pi family based on which Pauli operators
appear at that site: P1 handles {Y,Z}, P4 handles {X,Z}.
All three are product operators (Schmidt rank 1). **LOCAL.**

### No product Pi exists (NON-LOCAL symmetry)

| Hamiltonian | N | All products fail | Numerical rank | Entropy |
|-------------|---|-------------------|----------------|---------|
| XZ+YZ | 2 | 8/8 fail | 16 | ~3.4 bits |
| ZX+ZY | 2 | 8/8 fail | 16 | ~3.2 bits |

All 8 combinations of {P1, P4, M2} x {P1, P4, M2} tested. All fail.
The numerical Pi from eigenvalue pairing has full Schmidt rank (16/16).

**Physical reason for non-locality (signal engineering view):**

Z-dephasing splits each qubit's Pauli basis into a DC bus {I, Z}
(populations, immune to noise) and an AC bus {X, Y} (coherences, decaying).
The Pi operator is a crossover switch between DC and AC. Two crossover
settings exist:

| Crossover | X-channel | Y-channel | Z-channel |
|-----------|-----------|-----------|-----------|
| P1 | PASS | FAIL | compensated by dephasing |
| P4 | FAIL | PASS | compensated by dephasing |

P1 anticommutes with ad_X (handles X-signal routing) but not ad_Y.
P4 anticommutes with ad_Y (handles Y-signal routing) but not ad_X.
Neither anticommutes with ad_Z directly; Z-channel matching is provided
by the dephasing terms in the Liouvillian (the noise IS the impedance
match for the DC bus).

For XZ+YZ (X and Y on same site): site 0 needs BOTH crossover settings
simultaneously. This is a frequency conflict in a diplexer: two bands
need incompatible filter settings on the same port. No single crossover
works. The resolution requires a non-local "hybrid coupler" (entangled Pi)
that mixes the signals between sites before crossovering.

For XZ+ZY (X and Y on different sites): no conflict. Each site chooses
its own crossover independently (P4 on site 0 for X, P1 on site 1 for Y).

**Minimum Schmidt rank:** Rank-2, 3, and 4 truncations of the numerical Pi
all fail the conjugation equation (error > 1.0). Random search over the
40-dimensional commutant (5000 samples) found rank 10 as the best achievable.
The non-locality is not a simple 2-port hybrid but a multi-mode structure.
The exact minimum rank remains open.

---

## 4. The Space of Valid Pi Operators

An unexpected finding: Pi is not unique. Any operator of the form
Pi_0 * S, where S commutes with L (a symmetry of the Liouvillian), is
also a valid Pi operator. When L has degenerate eigenvalues, S can be
non-trivial and entangled.

This means:
- Numerical construction from eigenvector pairing finds an arbitrary
  element of the valid-Pi space, typically entangled.
- Even for Heisenberg (where a product Pi exists), the numerically
  constructed Pi has Schmidt rank 16.
- The physically meaningful question is the MINIMUM Schmidt rank across
  all valid Pi operators: this determines whether the symmetry itself
  is local or non-local.

For the standard cases: minimum rank = 1 (product exists).
For XZ+YZ and ZX+ZY: minimum rank > 1 (confirmed by exhaustive test of
all 8 per-site combinations from {P1, P4, M2} x {P1, P4, M2}).

The commutant dimension varies by Hamiltonian:

| Hamiltonian | dim(commutant) | Distinct eigenvalues |
|-------------|----------------|---------------------|
| Heisenberg | 44 | 7 of 16 |
| Ising ZZ | 64 | 4 of 16 |
| XX | 40 | 7 of 16 |
| XY | 40 | 7 of 16 |
| XZ+ZY | 36 | 9 of 16 |
| YZ+ZX | 36 | 9 of 16 |

More distinct eigenvalues -> smaller commutant -> less freedom in Pi
choice. The non-local cases have the smallest commutant (36), but even
that is a 36-dimensional space of valid Pi operators.

---

## 5. What This Means

### For the standard palindrome (Heisenberg, Ising, XX, ...)

The symmetry is local. Each qubit has its own P1 map (I<->X, Y<->Z)
that independently generates the palindromic pairing. Implications:
- The palindrome can be checked per site (local measurements suffice)
- Error correction based on palindrome breaking can use local syndromes
- The palindrome is not topologically protected

### For the 2 non-local cases (XZ+YZ, ZX+ZY)

The symmetry is entangled. No per-site decomposition exists. Both X and Y
appear on the same qubit site, creating incompatible demands on the
per-site mirror. Implications:
- Detecting palindrome breaking requires joint measurements across sites
- The palindrome may carry topological content (connects to Q2)
- Error correction requires entangled measurements
- The symmetry has a genuinely multi-body character that is invisible
  to any single-site observable

### For non-uniform product cases (XY, XZ+ZY, YZ+ZX)

The symmetry is local but non-uniform. Different sites carry different
palindromic maps (P1 on one, P4 or M2 on the other). The symmetry is
still a product operator (rank 1), but the per-site map depends on which
Pauli operators the Hamiltonian places at that site.

---

## 6. Connection to Other Questions

- **Q2 (Topological invariant):** Non-local Pi is a prerequisite for
  topological protection. Only the 2 non-local cases (XZ+YZ, ZX+ZY) are
  candidates. For the other 18 palindromic cases, the palindrome is a
  per-site property and cannot be topologically protected.

- **Hidden symmetry Q (THE_OTHER_SIDE):** The 12 parity-breaking
  palindrome-preservers must have hidden Q operators. The Q for local
  cases is another product operator. The Q for non-local cases is
  entangled. This constrains what Q can be.

- **Error correction (ERROR_CORRECTION_PALINDROME):** Whether the
  palindromic pairing can define error-correcting codespaces depends
  on whether Pi is local. Local Pi -> local stabilizer codes.
  Non-local Pi -> more exotic codes requiring entangled stabilizers.

---

## Scripts

```
PYTHONIOENCODING=utf-8 python simulations/pi_operator_entanglement.py
PYTHONIOENCODING=utf-8 python simulations/pi_nonlocality_mechanism.py
```

---

*See also:* [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) (exhaustive per-site search),
[Pi as Time Reversal](PI_AS_TIME_REVERSAL.md) (physical meaning of Pi),
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (Pi conjugation equation)
