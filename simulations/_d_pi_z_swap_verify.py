"""Verifies the strukturelle Identität:

    D · Π_Z · D = Π_Y

surfaced during the Welle 10d Task 1 spec-review audit (2026-05-27).

Where:
  - Π_Z, Π_Y are the F1 palindrome operators for Z- and Y-dephasing
    respectively (signed permutation matrices on the 4^N Pauli basis,
    constructed via the per-letter ActOnLetter rules).
  - D = diag((-1)^n_Y(k)) is the real diagonal unitary involution on
    4^N Pauli basis, where n_Y(k) counts the number of Y letters in
    the Pauli string indexed by k.

Significance:
  - The standard codebase pattern combines T (vec_F transform) with a
    vec_R-style commutator -i(H ⊗ I − I ⊗ H^T), producing Pauli-basis
    matrices that are D · L_natural · D rather than L_natural.
  - This identity says D-conjugation IS the Z↔Y dephasing-letter swap
    at the operator-space level. So all F1 residual norms / inner
    products / spectra computed via the "twisted" pipeline are
    automatically Z↔Y equivariant — a non-trivial structural symmetry.

Verification: numpy double precision, N = 1, 2, 3, 4. Bit-exact at every
N (residual matrix L_inf norm = 0 exactly, no machine-epsilon leakage).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.symmetry import build_pi_full

LETTERS = ["I", "X", "Z", "Y"]
BIT_B = {"I": 0, "X": 0, "Z": 1, "Y": 1}


def n_y_count(letters):
    return sum(1 for L in letters if L == "Y")


def build_D(N: int) -> np.ndarray:
    """D = diag((-1)^n_Y(k)) on 4^N Pauli basis (PauliIndex flat-encoding)."""
    from itertools import product
    strings = list(product(LETTERS, repeat=N))
    diag = np.array([(-1) ** n_y_count(s) for s in strings], dtype=complex)
    return np.diag(diag)


def main():
    print("Verifying D · Π_Z · D = Π_Y at N = 1, 2, 3, 4")
    print("=" * 72)
    print()

    all_pass = True
    for N in [1, 2, 3, 4]:
        # Build Π_Z and Π_Y via the framework's build_pi_full.
        pi_z = build_pi_full(N, dephase_letter="Z")
        pi_y = build_pi_full(N, dephase_letter="Y")
        D = build_D(N)

        # Compute D · Π_Z · D (note D is its own inverse since D² = I).
        lhs = D @ pi_z @ D
        rhs = pi_y

        diff = np.max(np.abs(lhs - rhs))
        status = "PASS" if diff == 0.0 else "FAIL"
        print(f"N={N}: 4^N = {4**N}, D shape = {D.shape}, Π_Z shape = {pi_z.shape}")
        print(f"       max|D·Π_Z·D − Π_Y| = {diff:.3e}  →  {status}")

        # Sanity: D is involution.
        d_squared = D @ D
        d_inv_diff = np.max(np.abs(d_squared - np.eye(4 ** N, dtype=complex)))
        print(f"       D · D = I (sanity): max|D² − I| = {d_inv_diff:.3e}")

        # Sanity: Π_Z and Π_Y are order 4.
        pi_z_4 = pi_z @ pi_z @ pi_z @ pi_z
        pi_y_4 = pi_y @ pi_y @ pi_y @ pi_y
        eye = np.eye(4 ** N, dtype=complex)
        print(f"       Π_Z^4 = I: max|Π_Z⁴ − I| = {np.max(np.abs(pi_z_4 - eye)):.3e}")
        print(f"       Π_Y^4 = I: max|Π_Y⁴ − I| = {np.max(np.abs(pi_y_4 - eye)):.3e}")
        print()

        if diff != 0.0:
            all_pass = False

    print("=" * 72)
    if all_pass:
        print("ALL PASS: D · Π_Z · D = Π_Y bit-exact at N = 1, 2, 3, 4")
        print("  → The vec_F vs vec_R convention twist in the standard codebase")
        print("    pipeline corresponds exactly to the Z↔Y dephasing-letter swap")
        print("    on operator space. F1 residual norms / inner products / spectra")
        print("    computed via the 'twisted' pipeline are automatically Z↔Y")
        print("    equivariant.")
    else:
        print("FAIL at one or more N — identity broken.")
        sys.exit(1)


if __name__ == "__main__":
    main()
