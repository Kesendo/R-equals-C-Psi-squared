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


def verify_per_site_identity_symbolic():
    """Verify d_l · π_Z_local · d_l = π_Y_local symbolically via sympy.

    The N-site identity D · Π_Z · D = Π_Y reduces, by the mixed-product
    property of the Kronecker product, to this single 4×4 check (Step 3
    of PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md). Verifying it symbolically
    in sympy with exact rationals + I gives zero machine-epsilon residual
    and closes the universal-N case in finite symbolic computation.

    Basis order: (I, X, Z, Y), matching PauliLetter packing a + 2·b.
    π_Z_local follows PiOperator.ActOnLetter(σ, PauliLetter.Z) rules:
        I → X · 1, X → I · 1, Z → Y · i, Y → Z · i.
    π_Y_local follows PiOperator.ActOnLetter(σ, PauliLetter.Y) rules:
        I → X · 1, X → I · 1, Z → Y · −i, Y → Z · −i.
    Column-as-input convention: pi[newK, k] = sign.
    """
    import sympy as sp

    I_s, mi = sp.I, -sp.I
    one, mone, zero = sp.Integer(1), sp.Integer(-1), sp.Integer(0)

    # π_Z_local on basis (I, X, Z, Y), column-as-input convention.
    pi_z = sp.Matrix([
        [zero, one,  zero, zero],
        [one,  zero, zero, zero],
        [zero, zero, zero, I_s ],
        [zero, zero, I_s,  zero],
    ])

    # π_Y_local on basis (I, X, Z, Y), column-as-input convention.
    pi_y = sp.Matrix([
        [zero, one,  zero, zero],
        [one,  zero, zero, zero],
        [zero, zero, zero, mi  ],
        [zero, zero, mi,   zero],
    ])

    # d_l = diag(1, 1, 1, -1): n_Y is 1 only for the Y entry (index 3).
    d = sp.diag(one, one, one, mone)

    lhs = d * pi_z * d
    rhs = pi_y
    diff = sp.simplify(lhs - rhs)

    is_zero = all(diff[i, j] == 0 for i in range(4) for j in range(4))
    return is_zero, diff


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

    print()
    print("=" * 72)
    print("Symbolic per-site identity check (sympy exact):")
    is_pass, diff_matrix = verify_per_site_identity_symbolic()
    status = "PASS" if is_pass else "FAIL"
    print(f"  d_l · π_Z_local · d_l == π_Y_local: {status}")
    if not is_pass:
        import sympy as sp
        print("  diff matrix:")
        sp.pprint(diff_matrix)
        sys.exit(1)
    print("  → Per-site 4×4 reduction closed symbolically; universal-N case")
    print("    follows by Kronecker product (mixed-product property).")
    print("    See docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md.")


if __name__ == "__main__":
    main()
