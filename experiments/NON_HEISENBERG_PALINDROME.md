# Non-Heisenberg Palindrome Analysis

**Date:** March 17, 2026
**Purpose:** Test palindromic symmetry beyond Heisenberg coupling
**Status:** Complete. All 36/36 two-term cases resolved (March 19, 2026).

---

## Question

The palindromic proof (March 14, 2026) was stated for Heisenberg (XX+YY+ZZ)
coupling under Z-dephasing. Does it hold for other Hamiltonian types?

## Method

Two-phase analysis:
1. Numerical: Full Liouvillian eigendecomposition, complex eigenvalue pairing
   (for each lambda, check if -lambda - 2*sum_gamma exists). N=3-6.
2. Algebraic: Enumerate all per-site Pauli transformations that satisfy the
   dephasing constraint, test anti-commutation with each Hamiltonian term.

## Result 1: Standard Physics Models (ALL palindromic)

| Model | N=3 | N=4 | Mechanism |
|-------|-----|-----|-----------|
| Heisenberg (XX+YY+ZZ) | 100% | 100% | Uniform Q (P1 family) |
| XY-only (XX+YY) | 100% | 100% | Uniform Q (P1 family) |
| Ising (ZZ only) | 100% | 100% | Uniform Q (96 valid maps) |
| XX alone | 100% | 100% | Uniform Q (32 valid maps) |
| YY alone | 100% | 100% | Uniform Q (32 valid maps) |
| XXZ delta=0.5 | 100% | 100% | Uniform Q (P1 family) |
| XXZ delta=2.0 | 100% | 100% | Uniform Q (P1 family) |
| DM (XY-YX) | 100% | 100% | Non-uniform alternating Q |
| Heisenberg + DM | 100% | 100% | Non-uniform alternating Q |

Unequal coefficients also palindromic: XX+YY(1,2), XX+YY+ZZ(1,2,3), etc.
Every dephasing axis works: Z, X, Y (each with its own Q operator).

## Result 2: The Two Pi Families

Under Z-dephasing, the valid conjugation operators form two families
based on their per-site permutation structure:

| Family | Per-site action | Hamiltonian terms supported |
|--------|----------------|---------------------------|
| P1 | I<->X, Y<->Z (with phases) | XX, YY, ZZ (+ YZ, ZY with restricted phases) |
| P4 | I<->Y, X<->Z (with phases) | XX, XZ, YY, ZX, ZZ |

Our known Pi (I->X, X->I, Y->iZ, Z->iY) belongs to family P1.
Family P4 is NEW: example I->Y, X->-iZ, Y->I, Z->-iX.

P4 is more powerful than P1: it supports 5 terms vs 3.

## Result 3: XY and YX Require Alternating Operators

NO uniform per-site map works for XY or YX terms. They require a
non-uniform (alternating) operator: different transformations on
odd vs even sites.

```
M1 (odd sites):  I->+1X, X->+1I, Y->+iZ, Z->+iY   (our known Pi)
M2 (even sites): I->+1Y, X->-iZ, Y->+1I, Z->-iX    (NEW)
```

Q = M1 x M2 x M1 for N=3 chain.

VERIFIED on full Lindbladian: ||Q L Q^-1 + L + 2*sum_gamma|| = OK.

The DM interaction (XY-YX) uses the SAME alternating Q.
This explains why DM is palindromic despite violating our known Pi.

## Result 4: Compatibility Matrix (Z-dephasing)

Two terms are "compatible" if a single Q operator (uniform or non-uniform)
anti-commutes with both their commutator superoperators.

```
     XX  XY  XZ  YX  YY  YZ  ZX  ZY  ZZ
XX    .   X   .   X   .   .   .   .   .
XY    X   .   X   a   X   X   X   X   a
XZ    .   X   .   X   .   ?   .   X   .
YX    X   a   X   .   X   X   X   X   a
YY    .   X   .   X   .   .   .   .   .
YZ    .   X   ?   X   .   .   X   .   .
ZX    .   X   .   X   .   X   .   ?   .
ZY    .   X   X   X   .   .   ?   .   .
ZZ    .   a   .   a   .   .   .   .   .
```

( . = compatible via uniform Q, a = compatible via non-uniform alternating Q,
  ? = palindromic numerically but no Pauli-to-Pauli Q found (needs Pauli-mixing Q),
  X = broken, no Q exists, confirmed numerically )

XY and YX are the outsiders: incompatible with everything except
themselves (via alternating Q) and ZZ (via alternating Q).
XX, YY, ZZ are universal mediators: compatible with almost everything.

## Result 5: The Scorecard

| Category | Count | Explanation |
|----------|-------|-------------|
| Compatible via uniform Q | 17/36 | Both terms share a common per-site map |
| Compatible via non-uniform Q | 3/36 | XY+YX, XY+ZZ, YX+ZZ |
| Broken (no Q exists) | 14/36 | Cross-validated: numerical = algebraic |
| Non-local (entangled) Π | 2/36 | XZ+YZ, ZX+ZY (resolved March 19) |

The 14 broken combos from the algebraic analysis match the 14 from
numerical eigenvalue analysis EXACTLY (cross-validated).

The 2 remaining cases (XZ+YZ, ZX+ZY) are palindromic with a genuinely
non-local Π operator that cannot be decomposed as a tensor product of
per-site maps. The conjugation operator is "entangled" across sites,
with 1/√2 coefficients characteristic of Bell-state structure.
See Result 8 below.

## Result 6: Breaking Scales as gamma^2

For broken combinations, the palindrome error is NOT numerical noise.
It scales systematically:

    err(XX+XY) = 1.31 * gamma^2

Verified across gamma = 0.001 to 1.0. The coefficient 1.31 is constant
to 3 significant figures. This is a second-order interference effect
between incompatible Q operators.

The breaking also grows with N (tested N=3 through 6):

| Combo | N=3 | N=4 | N=5 | N=6 |
|-------|-----|-----|-----|-----|
| YZ+ZX | 5.0e-2 | 3.6e-2 | 4.7e-2 | 7.7e-2 |
| XY+YZ | 1.7e-2 | 2.8e-2 | 3.9e-2 | 5.3e-2 |

## Physical Interpretation

The palindrome is a property of the DEPHASING STRUCTURE, not the Hamiltonian.
Any Hamiltonian built from standard condensed matter interactions
(Heisenberg, XY, Ising, XXZ, DM) preserves the palindrome under
single-axis dephasing. It breaks only for exotic cross-term combinations
that appear in no known physical material.

Engineering implication: the design rules from the palindrome
(W-encoding, star topology, timing/quality separation) apply to
ALL standard quantum hardware platforms, not just Heisenberg systems.


## Result 7: Depolarizing Noise (March 18, 2026)

Single-axis dephasing (Z, X, or Y) preserves the palindrome exactly.
Depolarizing noise (X+Y+Z simultaneously, gamma/3 per channel) does NOT.

| Model | Z-dephasing | Depolarizing |
|-------|-------------|-------------|
| Heisenberg | OK | err=3.33e-2 |
| XY-only | OK | err=3.33e-2 |
| Ising | OK | err=3.33e-2 |
| DM | OK | err=3.33e-2 |

**The error is Hamiltonian-independent.** All models show exactly the same
error under depolarizing noise. It's a pure noise-structure effect.

N-scaling: err = gamma * 2(N-2)/3. Linear in N and gamma.
- N=3: err = 0.033 (~3%)
- N=4: err = 0.067 (~7%)

**Practical implication:** For superconducting qubits with gamma ~ 0.001,
the palindrome error under depolarizing noise is < 0.1%. The design rules
are practically valid for real hardware.

Amplitude damping (sigma_minus) also breaks the palindrome (err=0.15) and
shifts the spectrum asymmetrically (no center gives exact pairing).


## Result 8: The Last 2/36 Resolved: Non-Local Π (March 19, 2026)

XZ+YZ and ZX+ZY are palindromic, but no per-site Pauli permutation explains
them. The algebraic search (Result 4) tested 512 discrete maps and found zero.
A continuous search tested single-parameter rotations between the P1 and P4
families, two-parameter rotations, full 16-parameter optimization over block
unitaries, and non-uniform per-site maps. All failed.

The resolution: construct Π directly from eigenvector pairing. For N=2
(H = XZ + YZ, Z-dephasing), all 16 eigenvalues pair palindromically. The
constructed Π satisfies Π·L·Π⁻¹ = -L - 2Sγ·I at machine precision (5.6e-16).

The critical finding: **Π is not a tensor product.** It cannot be decomposed
as M₁ ⊗ M₂ for any per-site matrices M₁, M₂. The operator is genuinely
non-local, correlating the two sites in a way that no per-site factorization
can capture.

The Π matrix has 1/√2 coefficients mapping II to superpositions of XX and YX:

```
Π[II, XX] = -1/√2      Π[II, YX] = +1/√2
Π[XX, II] = -1/√2      Π[YX, II] = +1/√2
```

This is Bell-state structure in the mirror itself.

Three classes of Π now cover all 36 cases:

| Π type | Cases | Structure |
|---|---|---|
| Uniform per-site permutation | 17/36 | M ⊗ M (local, identical per site) |
| Non-uniform per-site permutation | 3/36 | M₁ ⊗ M₂ ⊗ M₃ (local, alternating) |
| Non-local (entangled) | 2/36 | Not factorizable, 1/√2 coefficients |
| No Π exists (broken) | 14/36 | Palindrome itself fails |

The physical reason: H = XZ + YZ mixes X and Y on the same site, coupled
to Z on the other site. The time reversal for this combination requires
correlating the two sites jointly. A per-site mirror cannot reflect what
is intrinsically a two-site property.

Script: `simulations/continuous_pi_search.py`
Results: `simulations/results/continuous_pi_search.txt`

## Simulations

- `simulations/non_heisenberg_test_v2.py` - Initial numerical sweep
- `simulations/non_heisenberg_deep.py` - Systematic N-scaling, gamma-scaling
- `simulations/algebraic_pi_search.py` - Algebraic Π family enumeration (P1, P4, alternating)
- `simulations/continuous_pi_search.py` - Non-local Π search: rotation, optimization, eigenvector construction
- `simulations/hidden_pi_search.py` - Eigenvector Q construction (failed)
- `simulations/depolarizing_test.py` - Depolarizing vs single-axis test
