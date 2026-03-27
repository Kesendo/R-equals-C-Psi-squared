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
(XZ+ZY and YZ+ZX), the symmetry is genuinely NON-LOCAL. No per-site
factorization exists.

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

### No product Pi exists (NON-LOCAL symmetry)

| Hamiltonian | N | P1 error | Numerical rank | Entanglement entropy |
|-------------|---|----------|----------------|---------------------|
| XZ+ZY | 2 | 1.73 | 16 | 3.35 bits |
| YZ+ZX | 2 | 1.73 | 16 | 3.15 bits |

P1 fails with error O(1) (not a small correction, a complete failure).
No per-site map from either P1 or P4 family works. The numerical Pi
from eigenvalue pairing has full Schmidt rank (16/16), with operator
entanglement entropy near the maximum of 4.0 bits.

### Alternating family (XY coupling)

| Hamiltonian | N | P1 error | Status |
|-------------|---|----------|--------|
| XY | 2 | 1.41 | Alternating M1 x M2 (different per-site maps) |

P1 x P1 fails, but the alternating family uses different maps on different
sites: M1 x M2. This is still a product operator (rank 1), just not
uniform. The XY case is LOCAL with non-uniform per-site maps.

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
For XZ+ZY and YZ+ZX: minimum rank > 1 (proven by exhaustive per-site
search in [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md)).

---

## 5. What This Means

### For the standard palindrome (Heisenberg, Ising, XX, ...)

The symmetry is local. Each qubit has its own P1 map (I<->X, Y<->Z)
that independently generates the palindromic pairing. Implications:
- The palindrome can be checked per site (local measurements suffice)
- Error correction based on palindrome breaking can use local syndromes
- The palindrome is not topologically protected

### For the 2 non-local cases (XZ+ZY, YZ+ZX)

The symmetry is entangled. No per-site decomposition exists. Implications:
- Detecting palindrome breaking requires joint measurements across sites
- The palindrome may carry topological content (connects to Q2)
- Error correction requires entangled measurements
- The symmetry has a genuinely multi-body character that is invisible
  to any single-site observable

### For the alternating cases (XY, YX, ...)

Intermediate: the symmetry is local but non-uniform. Different sites
carry different palindromic maps. The symmetry is still a product
operator (rank 1), but the per-site structure varies.

---

## 6. Connection to Other Questions

- **Q2 (Topological invariant):** Non-local Pi is a prerequisite for
  topological protection. The 2 non-local cases are the candidates.
  For the 20 local/alternating cases, the palindrome cannot be
  topologically protected (it's a per-site property).

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
```

---

*See also:* [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) (exhaustive per-site search),
[Pi as Time Reversal](PI_AS_TIME_REVERSAL.md) (physical meaning of Pi),
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (Pi conjugation equation)
