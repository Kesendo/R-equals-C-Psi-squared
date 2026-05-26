"""F87 <-> F112 connection probe.

F87 trichotomy classifies (Pauli pair, dephase letter) as truly/soft/hard
based on M = Pi L Pi^-1 + L + 2 sigma I behavior:
  truly: M = 0
  soft:  M != 0 but spec(L) palindromic
  hard:  spec(L) breaks palindromy

F112 (just sketched) says: polarity_coordinates_from_L gives asymmetry = 0
when c (the dissipator's collapse operator) is bit_b-homogeneous.

In the standard F87 setup, c = single Pauli letter D (X, Y, or Z), which
is trivially bit_b-homogeneous. So F112 predicts asymmetry = 0 for ALL
F87-style L's, regardless of trichotomy class.

This probe verifies that prediction:
  - Pick representatives of each F87 class (truly, soft, hard) at N=3
  - Build L with H = the pair's bilinear chain, dephase = Z
  - Check polarity asymmetry

Hypothesis: asymmetry = 0 across all three classes.

If true: F112 does NOT subdivide F87 trichotomy. F87 lives in M's
Frobenius magnitude + palindrome structure; F112 lives in M_anti's
+i/-i Pi-eigenvalue split. Orthogonal structural axes.

If false (e.g., F87-hard sometimes breaks polarity balance):
F112 is a finer classification within F87-hard.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np

from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def ps(letters):
    op = PAULI[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI[l])
    return op


def L2P(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def site_op(N, site, letter):
    """Single Pauli on `site`, identity on other sites."""
    letters = ['I'] * N
    letters[site] = letter
    return ps(letters)


def build_L_pair(N, pair_letters, dephase_letter, gamma=0.1):
    """Build L with H = chain of pair bilinears, dephase = single letter.

    pair_letters: list of (a, b) tuples, e.g. [('X', 'X'), ('Y', 'Y')]
    Chain: H = sum over bonds b of sum over (a, b) in pair_letters of a_l b_{l+1}.
    """
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for a_letter, b_letter in pair_letters:
            letters = ['I'] * N
            letters[b] = a_letter
            letters[b + 1] = b_letter
            H = H + ps(letters)

    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    # Dephase on each site
    for l in range(N):
        Dl = site_op(N, l, dephase_letter)
        L = L + gamma * np.kron(Dl, Dl.conj())
    return L, H


def main():
    N = 3
    gamma = 0.1
    sigma = N * gamma

    # Representative pairs per F87 trichotomy class under Z-dephasing.
    # Truly: XX (Heisenberg subset, M = 0 by F1)
    # Soft (Pi^2-Z-even non-truly): YZ + ZY (F108 anomaly)
    # Hard (Pi^2-Z-odd): XY (the canonical Pi^2-Z-odd bilinear)
    cases = [
        ("TRULY: Heisenberg (XX+YY+ZZ)",
         [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')], 'Z'),
        ("TRULY: XX only",
         [('X', 'X')], 'Z'),
        ("SOFT (Pi^2-Z-even non-truly): YZ + ZY",
         [('Y', 'Z'), ('Z', 'Y')], 'Z'),
        ("HARD (Pi^2-Z-odd): XY pure",
         [('X', 'Y')], 'Z'),
        ("HARD (Pi^2-Z-odd): YX pure",
         [('Y', 'X')], 'Z'),
        ("HARD (Pi^2-Z-odd): XZ + ZX",
         [('X', 'Z'), ('Z', 'X')], 'Z'),
        ("MIXED truly + soft + hard",
         [('X', 'X'), ('Y', 'Z'), ('X', 'Y')], 'Z'),
    ]

    print(f"N = {N}, gamma = {gamma}, sigma = {sigma}")
    print(f"All cases use c = single Pauli letter (Z-deph), so F112 predicts asymmetry = 0.")
    print()
    print(f"{'Case':<60}  {'||M||^2':>12}  {'asym':>14}  {'rel_asym':>12}  F87 status")
    print("-" * 130)

    for label, pair_letters, dephase in cases:
        L_vec, H = build_L_pair(N, pair_letters, dephase, gamma)
        L_pauli = L2P(L_vec, N)
        result = polarity_coordinates_from_L(L_pauli, N, sigma)
        ns_M = result['norm_sq']['M']
        asym = result['asymmetry']
        rel = abs(asym) / max(ns_M, 1e-15)
        f112_marker = "BAL" if rel < 1e-10 else "BREAK"

        # Classify F87 status: M = 0 -> truly; M != 0 -> need to check palindrome
        if ns_M < 1e-10:
            f87_status = "TRULY (M=0)"
        else:
            # Check spectrum palindromicity (rough: spec(L) symmetric around -sigma?)
            # For demo, just report based on label
            if "TRULY" in label:
                f87_status = "TRULY"
            elif "SOFT" in label:
                f87_status = "SOFT"
            elif "HARD" in label:
                f87_status = "HARD"
            elif "MIXED" in label:
                f87_status = "MIXED"
            else:
                f87_status = "?"

        print(f"{label:<60}  {ns_M:>12.4f}  {asym:>+14.4e}  {rel:>12.4e}  {f112_marker} | {f87_status}")

    print()
    print("F87 <-> F112 connection summary:")
    print("=" * 80)
    print("If all asym = 0 across truly/soft/hard:")
    print("  F112 polarity balance is INSENSITIVE to F87 trichotomy.")
    print("  F87 lives in ||M||^2 magnitude + spectrum palindrome.")
    print("  F112 lives in M_anti's Pi +i/-i split.")
    print("  These are ORTHOGONAL structural axes that happen to share the bit_b grading.")
    print()
    print("If F87-hard cases break polarity balance:")
    print("  F112 sub-classifies F87-hard into bit_b-homogeneous-c-preserving vs not.")
    print()
    print("If F87-soft cases break polarity balance:")
    print("  F112 detects something invisible to F87 trichotomy.")


if __name__ == '__main__':
    main()
