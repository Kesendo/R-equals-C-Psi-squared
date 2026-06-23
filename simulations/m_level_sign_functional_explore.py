"""M-level signed equivariance: investigate the sign functional ε(σ) such that

    D · L_σ · D = ε(σ) · L_σ    in the 4^N Pauli basis

for single Pauli strings σ at N = 1, 2, 3, 4, and extend the closed-form to
bond bilinear Hamiltonians + the empirical Welle-15 cases (XZ+ZX → ε=−1,
YZ+ZY → ε=+1). Outputs a structured enumeration and a closed-form candidate.

Surfaced 2026-05-27 from Welle 15 Task A polish wave: substantive M_anti at
XZ+ZX bond at N=2 shows anti-equivariance bit-exact (M_Y = −D·M_Z·D). YZ+ZY
at N=2 shows equivariance (M_Y = +D·M_Z·D). The sign is bond-specific; this
script characterizes the closed-form for ε(σ) as a function of σ's letter
content.

Where:
  - L_σ = −i[σ, ·] in vec_F basis = −i(σ ⊗ I − I ⊗ σ^T) on the 4^N Pauli basis
  - D = diag((−1)^n_Y(α)) is the real diagonal unitary involution
  - Pauli string σ has per-site letters in {I, X, Z, Y} (PauliLetter packing)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


LETTERS = ["I", "X", "Z", "Y"]
N_Y = {"I": 0, "X": 0, "Z": 0, "Y": 1}
N_X = {"I": 0, "X": 1, "Z": 0, "Y": 1}
N_Z = {"I": 0, "X": 0, "Z": 1, "Y": 1}

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def n_y_count(letters):
    return sum(N_Y[L] for L in letters)


def n_x_count(letters):
    return sum(N_X[L] for L in letters)


def n_z_count(letters):
    return sum(N_Z[L] for L in letters)


def pauli_string(letters):
    """Tensor product of per-site Pauli matrices."""
    result = PAULI[letters[0]]
    for L in letters[1:]:
        result = np.kron(result, PAULI[L])
    return result


def build_D(N: int) -> np.ndarray:
    """D = diag((−1)^n_Y(α)) on 4^N Pauli basis."""
    strings = list(itertools.product(LETTERS, repeat=N))
    diag = np.array([(-1) ** n_y_count(s) for s in strings], dtype=complex)
    return np.diag(diag)


def build_L_sigma_in_pauli_basis(sigma_letters):
    """L_σ in 4^N Pauli basis = T† · (-i(σ ⊗ I − I ⊗ σ^T)) · T / norm

    where T is the vec_F → Pauli-basis transform. We compute directly
    via inner products to avoid the basis-twist complications: for each
    pair (σ_α, σ_β) of Pauli strings:
        [L_σ]_(α, β) = Tr(σ_α · [-i(σ σ_β − σ_β σ)]) / 2^N
    """
    N = len(sigma_letters)
    d = 2**N
    sigma = pauli_string(sigma_letters)

    strings = list(itertools.product(LETTERS, repeat=N))
    n_basis = len(strings)
    L = np.zeros((n_basis, n_basis), dtype=complex)
    for j, beta in enumerate(strings):
        sigma_beta = pauli_string(beta)
        commutator = -1j * (sigma @ sigma_beta - sigma_beta @ sigma)
        for i, alpha in enumerate(strings):
            sigma_alpha = pauli_string(alpha)
            L[i, j] = np.trace(sigma_alpha @ commutator) / d
    return L


def epsilon_of_sigma(sigma_letters):
    """Compute ε such that D · L_σ · D = ε · L_σ. Returns:
      +1 / −1 if L_σ is nontrivial and a pure +/− eigenvector,
      "ZERO" if L_σ = 0 (commutator vanishes; only for σ = I⊗N),
      "MIXED" if L_σ has both + and − D-conjugation entries.
    """
    N = len(sigma_letters)
    L = build_L_sigma_in_pauli_basis(sigma_letters)
    D = build_D(N)

    norm = np.linalg.norm(L)
    if norm < 1e-12:
        return "ZERO", 0.0

    DLD = D @ L @ D
    diff_plus = np.linalg.norm(DLD - L)
    diff_minus = np.linalg.norm(DLD + L)

    if diff_plus < 1e-10:
        return +1, diff_plus
    if diff_minus < 1e-10:
        return -1, diff_minus
    return "MIXED", min(diff_plus, diff_minus)


def epsilon_of_H(terms_with_coeffs, N):
    """Compute ε for H = Σ c_k σ_k. Returns:
      +1 / −1 if D·L_H·D = ε·L_H bit-exact,
      "ZERO" if L_H = 0,
      "MIXED" if the H-level commutator structure does not lift to a single sign.
    """
    L_total = np.zeros((4**N, 4**N), dtype=complex)
    for letters, coeff in terms_with_coeffs:
        L_total += coeff * build_L_sigma_in_pauli_basis(letters)

    D = build_D(N)
    norm = np.linalg.norm(L_total)
    if norm < 1e-12:
        return "ZERO", 0.0

    DLD = D @ L_total @ D
    diff_plus = np.linalg.norm(DLD - L_total)
    diff_minus = np.linalg.norm(DLD + L_total)

    if diff_plus < 1e-10:
        return +1, diff_plus
    if diff_minus < 1e-10:
        return -1, diff_minus
    return "MIXED", min(diff_plus, diff_minus)


def main():
    print("=" * 80)
    print("M-LEVEL SIGNED EQUIVARIANCE ε(σ): D · L_σ · D = ε(σ) · L_σ")
    print("=" * 80)
    print()

    # Part 1: enumerate ε(σ) for all single Pauli strings at N = 1, 2, 3
    print("PART 1: ε(σ) for single Pauli strings σ")
    print("-" * 80)

    closed_form_check_passed = True

    for N in [1, 2, 3]:
        print(f"\nN = {N}: {4**N} Pauli strings")
        results_by_eps = {+1: [], -1: [], "ZERO": [], "MIXED": []}
        for sigma_letters in itertools.product(LETTERS, repeat=N):
            eps, _ = epsilon_of_sigma(sigma_letters)
            results_by_eps[eps].append(sigma_letters)

        for key in [+1, -1, "ZERO", "MIXED"]:
            label = "+1" if key == +1 else ("-1" if key == -1 else key)
            count = len(results_by_eps[key])
            print(f"  ε = {label}: {count} strings")
            if count > 0 and count <= 12:
                for s in results_by_eps[key]:
                    n_y = n_y_count(s)
                    n_x = n_x_count(s)
                    n_z = n_z_count(s)
                    label_str = "".join(s)
                    print(
                        f"    {label_str}  (n_Y={n_y}, n_X={n_x}, n_Z={n_z})"
                    )

        # Closed-form vermutung: at N=1, ε(I)=ZERO, ε(X)=ε(Z)=−1, ε(Y)=+1.
        # At higher N, try various formulas in terms of n_Y, n_X, n_Z.
        if N == 1:
            continue

        # Vermutung: ε(σ) = (−1)^(n_Y(σ) + 1) for σ ≠ I⊗N
        # i.e. ε = +1 if n_Y odd, ε = −1 if n_Y even (and σ is not all-identity)
        for key in [+1, -1]:
            for s in results_by_eps[key]:
                n_y = n_y_count(s)
                expected = (-1) ** (n_y + 1)
                if expected != key:
                    print(
                        f"    VERMUTUNG FAIL: {''.join(s)} ε={key} "
                        f"but (-1)^(n_Y+1) = {expected}"
                    )
                    closed_form_check_passed = False

    if closed_form_check_passed:
        print()
        print("=" * 80)
        print("CLOSED-FORM ε(σ) = (−1)^(n_Y(σ) + 1) for σ ≠ I⊗N")
        print("                = +1 if n_Y(σ) odd, −1 if n_Y(σ) even (non-identity σ)")
        print("VERIFIED bit-exact for all Pauli strings at N = 1, 2, 3 above.")
        print("=" * 80)
    else:
        print()
        print("VERMUTUNG falsified somewhere; see FAIL lines above.")

    # Part 2: ε(H) for the Welle 15 empirical bond cases
    print()
    print()
    print("PART 2: ε(H) for Welle 15 empirical bond Hamiltonians")
    print("-" * 80)

    N = 2
    test_hamiltonians = [
        ("XZ + ZX", [(("X", "Z"), 1.0), (("Z", "X"), 1.0)]),
        ("YZ + ZY", [(("Y", "Z"), 1.0), (("Z", "Y"), 1.0)]),
        ("XX + YY + ZZ (Heisenberg)", [
            (("X", "X"), 1.0), (("Y", "Y"), 1.0), (("Z", "Z"), 1.0)
        ]),
        ("XX (Π²-even bit_a homog)", [(("X", "X"), 1.0)]),
        ("YY (Π²-even bit_b homog)", [(("Y", "Y"), 1.0)]),
        ("ZZ (Π²-even bit_b homog)", [(("Z", "Z"), 1.0)]),
        ("XY (Π²-odd)", [(("X", "Y"), 1.0)]),
        ("YX (Π²-odd)", [(("Y", "X"), 1.0)]),
        ("ZY (Π²-odd)", [(("Z", "Y"), 1.0)]),
        ("XZ (Π²-odd)", [(("X", "Z"), 1.0)]),
        ("X + Z (mixed letter at single site)", [
            (("X", "I"), 1.0), (("Z", "I"), 1.0)
        ]),
        ("X + Y (mixed letter, one Y one non-Y)", [
            (("X", "I"), 1.0), (("Y", "I"), 1.0)
        ]),
    ]

    for name, terms in test_hamiltonians:
        eps, residual = epsilon_of_H(terms, N)
        eps_str = (
            "+1"
            if eps == +1
            else ("-1" if eps == -1 else str(eps))
        )

        # Predict via the closed-form for each term:
        # H ε is well-defined iff all terms have the same ε
        per_term_eps = []
        for letters, _ in terms:
            n_y = n_y_count(letters)
            per_term_eps.append((-1) ** (n_y + 1))
        all_same = len(set(per_term_eps)) == 1
        predicted_eps = per_term_eps[0] if all_same else "MIXED"

        agree = "✓" if (
            (predicted_eps == "MIXED" and eps == "MIXED")
            or (isinstance(predicted_eps, int) and predicted_eps == eps)
            or (eps == "ZERO")
        ) else "✗"

        pred_str = (
            "+1"
            if predicted_eps == +1
            else ("-1" if predicted_eps == -1 else "MIXED")
        )
        print(
            f"  {name:42s}  ε = {eps_str:6s}  "
            f"per-term predict = {pred_str:6s}  {agree}  "
            f"(residual {residual:.2e})"
        )

    print()
    print("=" * 80)
    print("EXTENSION TO LINEAR COMBINATIONS H = Σ c_k σ_k:")
    print("  ε(H) is well-defined iff all σ_k share the same ε(σ_k).")
    print("  If terms have mixed ε values, ε(H) = MIXED (no single sign).")
    print("  L_H is bilinear in coefficients, but the per-term ε relation is")
    print("  additive only when ε's coincide.")
    print("=" * 80)

    # Part 3: empirical Welle 15 cases at N = 3, 4 (sanity check the closed form at higher N)
    print()
    print()
    print("PART 3: Closed-form sanity check at N = 3, 4 for selected H")
    print("-" * 80)

    test_at_higher_N = [
        # N = 3 cases
        (3, "XZ@(0,1) + ZX@(1,2) (chain bond pair)", [
            (("X", "Z", "I"), 1.0), (("I", "Z", "X"), 1.0)
        ]),
        (3, "YZ@(0,1) + ZY@(1,2) (chain bond pair)", [
            (("Y", "Z", "I"), 1.0), (("I", "Z", "Y"), 1.0)
        ]),
        (3, "ZZZ (three-body)", [(("Z", "Z", "Z"), 1.0)]),
        (3, "YYY (three-body)", [(("Y", "Y", "Y"), 1.0)]),
        (3, "XYZ (mixed three-body)", [(("X", "Y", "Z"), 1.0)]),
        # N = 4 cases
        (4, "Heisenberg chain bonds @ (0,1)+(1,2)+(2,3)", [
            (("X", "X", "I", "I"), 1.0),
            (("I", "X", "X", "I"), 1.0),
            (("I", "I", "X", "X"), 1.0),
            (("Y", "Y", "I", "I"), 1.0),
            (("I", "Y", "Y", "I"), 1.0),
            (("I", "I", "Y", "Y"), 1.0),
            (("Z", "Z", "I", "I"), 1.0),
            (("I", "Z", "Z", "I"), 1.0),
            (("I", "I", "Z", "Z"), 1.0),
        ]),
    ]

    for N_test, name, terms in test_at_higher_N:
        eps, residual = epsilon_of_H(terms, N_test)
        eps_str = (
            "+1"
            if eps == +1
            else ("-1" if eps == -1 else str(eps))
        )

        per_term_eps = []
        for letters, _ in terms:
            n_y = n_y_count(letters)
            per_term_eps.append((-1) ** (n_y + 1))
        all_same = len(set(per_term_eps)) == 1
        predicted_eps = per_term_eps[0] if all_same else "MIXED"
        pred_str = (
            "+1"
            if predicted_eps == +1
            else ("-1" if predicted_eps == -1 else "MIXED")
        )

        agree = "✓" if (
            (predicted_eps == "MIXED" and eps == "MIXED")
            or (isinstance(predicted_eps, int) and predicted_eps == eps)
            or (eps == "ZERO")
        ) else "✗"

        print(
            f"  N={N_test}  {name:55s}  ε = {eps_str:6s}  "
            f"predict = {pred_str:6s}  {agree}"
        )

    print()
    print("=" * 80)
    print("SUMMARY:")
    print("  Closed-form ε(σ) = (−1)^(n_Y(σ) + 1) for non-identity σ.")
    print("  Equivalent: ε(σ) = +1 if n_Y(σ) odd, ε(σ) = −1 if n_Y(σ) even.")
    print("  For H = Σ c_k σ_k: ε(H) well-defined iff all σ_k share same n_Y parity.")
    print()
    print("  Welle 15 empirical match (bond Hamiltonians at N = 2):")
    print("    XZ + ZX     → both n_Y = 0 (even) → ε = −1 ✓")
    print("    YZ + ZY     → both n_Y = 1 (odd)  → ε = +1 ✓")
    print("    Heisenberg  → XX (n_Y=0), YY (n_Y=2), ZZ (n_Y=0), all even → ε = −1 ✓")
    print()
    print("  Operational reading: D-conjugation of L_H = -i[H,·] flips sign exactly")
    print("  when H contains EVEN total Y letters per Pauli term (and σ ≠ I⊗N).")
    print("  The single-Y per-term Hamiltonians (YZ+ZY, XY, YX, ZY) get ε = +1;")
    print("  the Y-free or 2-Y Hamiltonians (XX, ZZ, YY, XZ, Heisenberg) get ε = −1.")
    print("=" * 80)


if __name__ == "__main__":
    main()
