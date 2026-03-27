# Is the Palindromic Mirror Entangled in Operator Space?

<!-- Keywords: Pi operator entanglement operator space, palindromic conjugation
Schmidt decomposition, non-local symmetry Liouville space, operator Schmidt rank
per-site factorization, product superoperator test, XZ YZ non-local palindrome,
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

The palindromic spectral symmetry of open quantum systems is generated
by a conjugation operator Pi. This operator acts on Liouville space,
which for multi-qubit systems has a natural tensor product structure
across qubit sites. We ask: can Pi be decomposed into independent
per-site operators (Pi = Pi_A ⊗ Pi_B), or is it inherently non-local?

**Result:** For all standard Hamiltonians (Heisenberg, Ising, XX, XY,
and combinations), Pi is a product operator: each qubit has its own
independent mirror. For exactly 2 of 36 two-term Pauli coupling
combinations, Pi is genuinely non-local (minimum operator Schmidt
rank = 9 out of 16). The non-locality arises when two incompatible
Pauli operators (X and Y) appear on the same qubit site, creating a
signal routing conflict that no single per-site crossover can resolve.

---

## 1. Why This Question Matters

The palindromic symmetry pairs every Liouvillian eigenmode with a
partner at -lambda - 2Sg. If this symmetry is a per-site property
(Pi = Pi_A ⊗ Pi_B), then:
- The palindrome can be checked and exploited locally
- Error correction based on palindrome breaking uses local syndromes
- The structure is not topologically protected

If Pi is non-local (entangled across sites), then:
- Detecting palindrome breaking requires joint measurements
- Error correction needs entangled stabilizers
- A topological invariant may protect the palindrome

Understanding which case applies determines what the palindrome is
useful for in practice.

---

## 2. Setup and Method

### Notation

"XZ" denotes the Hamiltonian coupling term sigma_X ⊗ sigma_Z (Pauli X
on site 0 tensored with Pauli Z on site 1). "XZ+YZ" means
H = J(sigma_X ⊗ sigma_Z + sigma_Y ⊗ sigma_Z), where both X and Y
act on site 0, and Z acts on site 1.

### Operator Schmidt decomposition

Pi is a matrix on Liouville space (dimension d^2 = 16 for 2 qubits).
This space has a tensor product structure: C^4_A ⊗ C^4_B, where
each C^4 represents one qubit's operator space (spanned by I, X, Y, Z).

To check whether Pi factorizes as Pi_A ⊗ Pi_B, we reshape it: group
the row index into (site-A-output, site-B-output) and the column index
into (site-A-input, site-B-input), then rearrange into a matrix M with
rows (A-output, A-input) and columns (B-output, B-input). The singular
values of M are the operator Schmidt coefficients. If only one is
nonzero (rank 1), Pi is a product. If more are nonzero, Pi is entangled.

### Testing strategy

Two per-site Pi families are known from
[Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md):
- **P1:** swaps I↔X and Y↔Z (with phases)
- **P4:** swaps I↔Y and X↔Z (with phases)

For each Hamiltonian, we test all combinations of per-site maps:
P1⊗P1, P4⊗P4, P1⊗P4, P4⊗P1, and the alternating map M2 (a variant
of P4 with different phases). If any combination satisfies the
conjugation equation Pi L Pi^{-1} = -L - 2Sg I to machine precision,
the symmetry is local. If all fail, we construct Pi numerically from
eigenvalue pairing and analyze its Schmidt rank.

**Key methodological insight:** Earlier work
([Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md)) tested
only UNIFORM products (same map on all sites). Testing NON-UNIFORM
products (different maps per site) resolved 2 cases previously
classified as non-local (XZ+ZY and YZ+ZX) as local.

---

## 3. Results

### Local cases (product Pi exists, Schmidt rank 1)

**Uniform product (same map on all sites):**

| Hamiltonian | N | Pi family | Conjugation error | Rank |
|-------------|---|-----------|-------------------|------|
| XX+YY+ZZ (Heisenberg) | 2 | P1⊗P1 | 1.1e-17 | 1 |
| ZZ (Ising) | 2 | P1⊗P1 | 2.0e-17 | 1 |
| XX | 2 | P1⊗P1 | 2.0e-17 | 1 |
| XX+YY+ZZ chain | 3 | P1⊗P1⊗P1 | 2.1e-17 | 1 |
| ZZ chain | 3 | P1⊗P1⊗P1 | 3.6e-17 | 1 |

**Non-uniform product (different maps per site):**

| Hamiltonian | N | Product Pi | Conjugation error | Rank |
|-------------|---|------------|-------------------|------|
| XY | 2 | P1⊗M2 | 2.0e-17 | 1 |
| XZ+ZY | 2 | P4⊗P1 | 1.4e-17 | 1 |
| YZ+ZX | 2 | P1⊗P4 | 1.4e-17 | 1 |

For XZ+ZY: X appears on site 0, Y on site 1. Site 0 uses P4
(handles X-routing), site 1 uses P1 (handles Y-routing).
No conflict because the incompatible requirements sit on different sites.

### Non-local cases (no product Pi exists)

| Hamiltonian | N | Products tested | All fail? | Min. rank |
|-------------|---|-----------------|-----------|-----------|
| XZ+YZ | 2 | 8 of 8 | yes | 9 |
| ZX+ZY | 2 | 8 of 8 | yes | 9 |

Both X and Y appear on the SAME site. All 8 combinations from
{P1, P4, M2} ⊗ {P1, P4, M2} fail. The minimum Schmidt rank is 9
(out of 16 maximum), found by systematic optimization over the
40-dimensional space of valid Pi operators.

---

## 4. The Signal Routing Mechanism

### DC/AC bus structure

Z-dephasing splits each qubit's Pauli basis into two buses:
- **DC bus** {I, Z}: immune to noise (populations, decided outcomes)
- **AC bus** {X, Y}: decaying under noise (coherences, undecided)

Pi is a crossover switch that swaps DC and AC. Two settings exist:

### Crossover compatibility

Each crossover anticommutes with one Pauli channel's adjoint action
and fails the other:

| Crossover | X-channel ({Pi, ad_X} = 0?) | Y-channel ({Pi, ad_Y} = 0?) |
|-----------|----------------------------|----------------------------|
| P1 | PASS | FAIL |
| P4 | FAIL | PASS |

Neither P1 nor P4 anticommutes with ad_Z. The Z-channel compatibility
is provided by the dephasing terms themselves: the noise IS the
impedance match for the DC bus.

### The diplexer conflict

For XZ+YZ (X and Y on same site): site 0 must route BOTH the X-channel
and the Y-channel simultaneously. P1 routes X but fails Y. P4 routes Y
but fails X. No single crossover setting handles both bands. This is
a frequency conflict in a diplexer: two bands need incompatible filter
settings on the same port.

For XZ+ZY (X and Y on different sites): no conflict. Site 0 uses P4
for X, site 1 uses P1 for Y. Each port handles one band.

### Minimum coupling depth

The non-local Pi that resolves the conflict requires 9 coupled modes
out of 16 (Schmidt rank 9). Systematic optimization (L-BFGS-B, 20
restarts per target rank) shows ranks 2 through 8 are unreachable.
This is not a simple 2-port hybrid coupler but a multi-mode MIMO
structure involving more than half of Liouville space.

---

## 5. The Space of Valid Pi Operators

Pi is not unique. Any operator Pi_0 S, where S commutes with L
(a symmetry of the Liouvillian), is also a valid Pi operator. The space
of valid operators is isomorphic to the commutant of L.

| Hamiltonian | Commutant dim. | Distinct evals | Local? |
|-------------|----------------|----------------|--------|
| Heisenberg | 44 | 7 of 16 | yes (P1⊗P1) |
| Ising ZZ | 64 | 4 of 16 | yes (P1⊗P1) |
| XX | 40 | 7 of 16 | yes (P1⊗P1) |
| XY | 40 | 7 of 16 | yes (P1⊗M2) |
| XZ+ZY | 36 | 9 of 16 | yes (P4⊗P1) |
| XZ+YZ | 40 | 9 of 16 | **no** (rank 9) |

More eigenvalue degeneracy means larger commutant and more freedom in
choosing Pi. Even for the smallest commutant (dimension 36), the space
of valid operators is vast. This is why numerical construction from
eigenvector pairing always produces entangled Pi, even when a product
solution exists. The physically meaningful question is not "is the
numerical Pi entangled?" but "does a product Pi exist anywhere in the
valid space?"

---

## 6. Implications

- **Error correction:** For 18 of 20 palindromic cases, palindrome
  breaking can be detected with local measurements (per-site syndromes).
  Only for XZ+YZ and ZX+ZY are entangled measurements needed.

- **Topological protection (Q2):** Only the 2 non-local cases could
  carry a topological invariant. For the other 18, the palindrome is a
  per-site property and not topologically protected.

- **Hidden symmetry Q:** The 12 parity-breaking palindrome-preservers
  ([THE_OTHER_SIDE](../hypotheses/THE_OTHER_SIDE.md)) must have hidden Q
  operators. For local palindromic cases, Q is also a product operator.
  For the 2 non-local cases, Q is entangled.

---

## Scripts

```
PYTHONIOENCODING=utf-8 python simulations/pi_operator_entanglement.py
PYTHONIOENCODING=utf-8 python simulations/pi_nonlocality_mechanism.py
```

---

*See also:* [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md),
[Pi as Time Reversal](PI_AS_TIME_REVERSAL.md),
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
