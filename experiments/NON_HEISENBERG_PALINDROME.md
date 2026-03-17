# Non-Heisenberg Palindrome Analysis

**Date:** March 17, 2026
**Purpose:** Test palindromic symmetry beyond Heisenberg coupling
**Status:** New experimental result

---

## Question

The palindromic proof (March 14, 2026) was stated for Heisenberg (XX+YY+ZZ)
coupling. Does it hold for other Hamiltonian types?

## Method

N=3 chain, Z-dephasing gamma=0.05. Full complex eigenvalue pairing:
for each eigenvalue lambda, check if -lambda - 2*sum_gamma exists.

## Results: Standard Models (ALL palindromic)

| Model | Palindromic | Error |
|-------|-------------|-------|
| Heisenberg (XX+YY+ZZ) | YES | < 1e-10 |
| XY-only (XX+YY) | YES | < 1e-10 |
| Ising (ZZ only) | YES | < 1e-10 |
| XX alone | YES | < 1e-10 |
| YY alone | YES | < 1e-10 |
| DM (XY-YX) | YES | < 1e-10 |
| XXZ delta=0.5 | YES | < 1e-10 |
| XXZ delta=2.0 | YES | < 1e-10 |
| Heisenberg + DM | YES | < 1e-10 |

**Confirmed at N=3 and N=4.** Every single Pauli pair (all 9: XX through ZY)
is palindromic individually. All standard physics models are palindromic.

## Results: Unequal Coefficients (palindromic)

| Model | Palindromic | Error |
|-------|-------------|-------|
| XX+YY (1,2) | YES | < 1e-10 |
| XX+YY+ZZ (1,2,3) | YES | < 1e-10 |
| XX+ZZ (1,1) | YES | < 1e-10 |
| XX+ZZ (1,2) | YES | < 1e-10 |

**Key finding:** The palindrome does NOT require equal coupling coefficients.
J_XX != J_YY != J_ZZ is fine. This is broader than the original proof suggested.

## Results: Dephasing Axis (all palindromic)

| Dephasing | Model | Palindromic |
|-----------|-------|-------------|
| Z | Heisenberg | YES |
| X | Heisenberg | YES |
| Y | Heisenberg | YES |

The palindrome holds for any single-axis dephasing, not just Z.
Each axis has its own Pi operator (related by rotation).

## Results: Cross-Term Combinations (mixed)

Under Z-dephasing, systematic test of all two-term combinations:

**Always palindromic:**
- ZZ + any other term (all 6 combinations)
- XX+XZ, XX+ZX, XX+YZ, XX+ZY
- YY+YZ, YY+XZ, YY+ZX, YY+ZY
- XZ+YZ

**Slightly broken (err ~ 3e-3 to 5e-2):**
- XX+XY (err=3.3e-3)
- XX+YX (err=3.3e-3)
- YY+XY (err=3.3e-3)
- YY+YX (err=3.3e-3)
- XY+XZ (err=8.9e-3)
- XY+YZ (err=1.7e-2)
- XY+ZX (err=1.7e-2)
- XY+ZY (err=8.9e-3)
- XZ+ZY (err=5.0e-2)

**Severely broken:**
- Random Hamiltonians (different coefficients per bond): 0% match at tol=1e-6

## Pattern

The rule for Z-dephasing appears to be:
1. Any term involving Z on at least one site pairs safely with anything
2. The "diagonal" pairs (XX, YY, ZZ) are safe with arbitrary coefficients
3. Cross-terms (XY, YX) break the palindrome when combined with
   terms sharing a non-Z Pauli on the same site (XX+XY, YY+XY)
4. The breaking is small (3e-3) for structured Hamiltonians but
   severe for random ones

The dephasing axis defines which Pauli is "safe". For Z-dephasing,
ZZ-type terms are always safe. XY-type cross-terms are the troublemakers.

## Physical Interpretation

All physically relevant spin models are palindromic:
- Heisenberg (magnetic ordering): palindromic
- XY model (superfluidity): palindromic
- Ising (classical limit): palindromic
- XXZ (anisotropic magnets): palindromic
- DM interaction (spin-orbit coupling): palindromic
- Heisenberg + DM (real materials): palindromic

The palindrome breaks only for Hamiltonians with cross-terms that
no standard condensed matter model uses. The symmetry is broader
than proven but not universal.

## Open Questions

1. Can the proof be extended to cover all "diagonal + antisymmetric" tensors?
2. What is the analytical mechanism by which XY+XX breaks the palindrome?
3. Is the small breaking (3e-3) a sign of an approximate symmetry, or
   a true symmetry breaking that happens to be numerically small at N=3?
4. Does the breaking grow or shrink with N?

## Simulations

- `simulations/non_heisenberg_test_v2.py` - Full systematic analysis
