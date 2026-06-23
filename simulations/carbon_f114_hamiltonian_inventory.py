"""F114 ε-classification inventory over Carbon-relevant Hamiltonian terms.

The bridge from formula to understanding: F114's closed form
ε(σ) = (−1)^(n_Y(σ) + 1) for σ ≠ I⊗N classifies any Hamiltonian term
into one of three families:

  ε = −1:  n_Y even per term (real-symmetric operator content)
  ε = +1:  n_Y odd per term (TR-odd operator content under D-conjugation)
  Mixed:   terms split across both parity classes; ε(H) not well-defined

For a sum H = Σ c_k σ_k, ε(H) is well-defined iff all σ_k share the same
n_Y parity. Otherwise ε(H) = Mixed.

D = diag((−1)^n_Y(α)) on the 4^N Pauli basis is the operator-space lift of
complex conjugation K, which is time reversal for spinless systems. ε(σ)
is precisely the time-reversal-parity sign of L_σ = −i[σ, ·] under D-conj:
    T · L_σ · T⁻¹ = ε(σ) · L_σ    (TR-parity of the commutator superoperator)

We let F114 classify a realistic Carbon Hamiltonian inventory: π-electron
hopping after Jordan-Wigner, Hubbard density-density, Heisenberg spin
exchange, external B-field Zeeman terms (each direction separately),
spin-orbit type bilinears, magnetic ring current, Dzyaloshinskii-Moriya
antisymmetric exchange, and realistic combinations.

Each row reports:
  - Term (chemistry name + Pauli operator content)
  - n_Y per string, parity per string
  - Predicted ε via closed form
  - Verified ε via direct matrix computation D · L · D vs ε · L

The bridge: see which Carbon Hamiltonian terms map to which F114 class.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


LETTERS = ["I", "X", "Z", "Y"]
PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def n_y(letters):
    return sum(1 for L in letters if L == "Y")


def pauli_string_op(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def letters_from_flat(k, N):
    out = []
    ki = k
    for _ in range(N):
        out.append(LETTERS[ki & 3])
        ki >>= 2
    return out


def flat_from_letters(letters):
    k = 0
    for i, L in enumerate(letters):
        k += LETTERS.index(L) << (2 * i)
    return k


def vec_to_pauli_transform(N):
    d = 2**N
    n_basis = 4**N
    T = np.zeros((d * d, n_basis), dtype=complex)
    for k in range(n_basis):
        sigma = pauli_string_op(letters_from_flat(k, N))
        T[:, k] = sigma.flatten("F")
    return T


def commutator_superop_pauli_basis(H, N):
    """L_H = -i[H, ·] in the 4^N Pauli basis."""
    d = 2**N
    I = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, I) - np.kron(I, H.T))
    T = vec_to_pauli_transform(N)
    return T.conj().T @ L_vec @ T / d


def build_D(N):
    """D = diag((-1)^n_Y(α)) on 4^N Pauli basis."""
    n_basis = 4**N
    diag = np.array(
        [(-1) ** n_y(letters_from_flat(k, N)) for k in range(n_basis)],
        dtype=complex,
    )
    return np.diag(diag)


def predicted_epsilon_for_term(letters):
    """ε for a single Pauli string per F114 closed form."""
    if all(L == "I" for L in letters):
        return "Zero"
    return +1 if n_y(letters) % 2 == 1 else -1


def predicted_epsilon_for_H(term_list):
    """ε for H = Σ σ_k. Well-defined iff all non-identity terms share n_Y parity."""
    seen_plus = False
    seen_minus = False
    for letters in term_list:
        eps = predicted_epsilon_for_term(letters)
        if eps == "Zero":
            continue
        if eps == +1:
            seen_plus = True
        else:
            seen_minus = True
    if seen_plus and seen_minus:
        return "Mixed"
    if not seen_plus and not seen_minus:
        return "Zero"
    return +1 if seen_plus else -1


def verified_epsilon_matrix_level(term_list, N, tol=1e-10):
    """Verify ε by computing D · L_H · D and comparing to ±L_H bit-exact.
    Returns +1, -1, "Mixed" (with residual), or "Zero"."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for letters in term_list:
        H = H + pauli_string_op(letters)
    if np.linalg.norm(H) < 1e-12:
        return "Zero", 0.0
    L = commutator_superop_pauli_basis(H, N)
    if np.linalg.norm(L) < 1e-12:
        return "Zero", 0.0
    D = build_D(N)
    DLD = D @ L @ D
    diff_plus = np.linalg.norm(DLD - L)
    diff_minus = np.linalg.norm(DLD + L)
    if diff_plus < tol:
        return +1, diff_plus
    if diff_minus < tol:
        return -1, diff_minus
    return "Mixed", min(diff_plus, diff_minus)


# Carbon Hamiltonian inventory. Each entry:
#   (chemistry_label, term_list, physical_interpretation)
# term_list is a list of Pauli-letter tuples; each tuple = one Pauli string in H.
# For multi-site, work at N=2 for clarity (concepts generalize to N-ring).

N_TEST = 2

INVENTORY = [
    # === Hückel core: pure electronic, real-symmetric ===
    ("Hopping (X⊗X)",
     [("X", "X")],
     "real part of Hückel π-hopping after Jordan-Wigner"),
    ("Hopping (Y⊗Y)",
     [("Y", "Y")],
     "imag part of Hückel π-hopping after Jordan-Wigner"),
    ("Bond operator B = X⊗X + Y⊗Y",
     [("X", "X"), ("Y", "Y")],
     "full Hückel hopping per bond; the Hückel atom of the framework"),

    # === Density-density / Coulomb ===
    ("Density-density Z⊗Z",
     [("Z", "Z")],
     "on-site Coulomb (Hubbard U: n_a · n_b)"),
    ("Heisenberg (X⊗X + Y⊗Y + Z⊗Z)",
     [("X", "X"), ("Y", "Y"), ("Z", "Z")],
     "isotropic spin exchange (real-coefficient)"),

    # === External Zeeman field per direction ===
    ("Zeeman x: σ_x on site 0",
     [("X", "I")],
     "external B-field in x direction (TR-even operator)"),
    ("Zeeman y: σ_y on site 0",
     [("Y", "I")],
     "external B-field in y direction (TR-ODD operator, F114 sees it)"),
    ("Zeeman z: σ_z on site 0",
     [("Z", "I")],
     "external B-field in z direction (TR-even operator)"),

    # === Spin-orbit type terms ===
    ("Spin-orbit Y⊗Y",
     [("Y", "Y")],
     "two TR-odd factors → TR-even composite; same class as XX/ZZ"),
    ("Spin-orbit X⊗Y",
     [("X", "Y")],
     "single Y per term → TR-odd Pauli string"),
    ("Spin-orbit Y⊗X",
     [("Y", "X")],
     "single Y per term → TR-odd Pauli string"),

    # === Magnetic ring current (induced current operator on a bond) ===
    ("Ring current Y⊗Z + Z⊗Y",
     [("Y", "Z"), ("Z", "Y")],
     "induced magnetic current per bond (both terms n_Y = 1)"),
    ("Ring current Y⊗Z − Z⊗Y",
     [("Y", "Z"), ("Z", "Y")],
     "antisymmetric ring current (same n_Y parity; signs don't change ε)"),

    # === Dzyaloshinskii-Moriya ===
    ("DM (X⊗Y + Y⊗X)",
     [("X", "Y"), ("Y", "X")],
     "Dzyaloshinskii-Moriya antisymmetric exchange (both terms n_Y = 1)"),

    # === Mixed-parity Hamiltonians (realistic Carbon) ===
    ("Hückel + Zeeman_y",
     [("X", "X"), ("Y", "Y"), ("Y", "I")],
     "real Hückel + external y-field (TR-breaking perturbation)"),
    ("Heisenberg + Zeeman_y",
     [("X", "X"), ("Y", "Y"), ("Z", "Z"), ("Y", "I")],
     "isotropic exchange + y-field"),
    ("Hückel + DM",
     [("X", "X"), ("Y", "Y"), ("X", "Y")],
     "Hückel + Dzyaloshinskii-Moriya (TR-breaking via cross-Pauli term)"),
    ("Hückel + Ring current",
     [("X", "X"), ("Y", "Y"), ("Y", "Z"), ("Z", "Y")],
     "Hückel + induced ring current (magnetic-field-induced TR breaking)"),
]


def main():
    print("=" * 90)
    print("F114 ε-classification inventory over Carbon-relevant Hamiltonian terms (N = 2)")
    print("=" * 90)
    print()
    print("Layout per row:")
    print("  - Chemistry label")
    print("  - n_Y parity profile (string n_Y values across terms)")
    print("  - Closed-form predicted ε  (from F114 rule)")
    print("  - Matrix-verified ε       (direct D · L_H · D check)")
    print("  - Physical interpretation")
    print()

    print(f"{'Label':<35} {'n_Y per term':<18} {'predicted ε':<12} {'verified ε':<14} {'agree':<6}")
    print("-" * 90)

    classification_count = {-1: 0, +1: 0, "Mixed": 0, "Zero": 0}
    for label, term_list, interp in INVENTORY:
        n_ys = [n_y(letters) for letters in term_list]
        predicted = predicted_epsilon_for_H(term_list)
        verified, residual = verified_epsilon_matrix_level(term_list, N_TEST)
        agree = "✓" if predicted == verified else "✗"
        if isinstance(verified, int) or verified in ("Mixed", "Zero"):
            classification_count[verified] = classification_count.get(verified, 0) + 1
        pred_str = str(predicted) if predicted in ("Mixed", "Zero") else f"{predicted:+d}"
        ver_str = str(verified) if verified in ("Mixed", "Zero") else f"{verified:+d}"
        print(f"{label:<35} {str(n_ys):<18} {pred_str:<12} {ver_str:<14} {agree:<6}")

    print()
    print("Physical interpretation summary:")
    print()
    for label, term_list, interp in INVENTORY:
        predicted = predicted_epsilon_for_H(term_list)
        pred_str = str(predicted) if predicted in ("Mixed", "Zero") else f"ε={predicted:+d}"
        print(f"  {label:<35} {pred_str:<8}  {interp}")

    print()
    print("=" * 90)
    print("CLASSIFICATION COUNTS")
    print("=" * 90)
    print(f"  ε = −1 (TR-even / real-symmetric class):  {classification_count.get(-1, 0)} terms")
    print(f"  ε = +1 (TR-odd / single-Y per term):      {classification_count.get(+1, 0)} terms")
    print(f"  ε = Mixed (coexisting parities):          {classification_count.get('Mixed', 0)} terms")
    print(f"  ε = Zero (trivial / identity):            {classification_count.get('Zero', 0)} terms")

    print()
    print("=" * 90)
    print("READING")
    print("=" * 90)
    print()
    print("The ε = −1 family (n_Y even per term) groups all real-symmetric Pauli")
    print("operators on the 2-site Hückel-type structure: π-electron hopping (XX, YY,")
    print("B = XX+YY), density-density (ZZ), Heisenberg exchange (XX+YY+ZZ), Zeeman in")
    print("x or z direction. These transform identically under D-conjugation: a pure −1.")
    print()
    print("The ε = +1 family (n_Y odd per term) groups all single-Y Pauli content:")
    print("Zeeman in y direction, single-Y bilinears (XY, YX, YZ, ZY etc.), magnetic-")
    print("field couplings to the imaginary axis of spin, induced ring currents,")
    print("Dzyaloshinskii-Moriya antisymmetric exchange. These transform with a +1.")
    print()
    print("Realistic Carbon Hamiltonians (Hückel + Zeeman_y, Heisenberg + DM, etc.)")
    print("often live in ε = Mixed: terms split across both parity classes, so the")
    print("F114 sign does not collapse to a single number on H. D-conjugation gives a")
    print("per-term-signed superposition rather than a global rescaling.")
    print()
    print("F114 ε is a closed-form TR-parity diagnostic at the operator-algebra level.")
    print("The classification opens three families on the Carbon Hamiltonian space; a")
    print("single number per Hamiltonian (or 'Mixed') summarises a structural property")
    print("that is otherwise invisible without spelling out the Pauli decomposition.")


if __name__ == "__main__":
    main()
