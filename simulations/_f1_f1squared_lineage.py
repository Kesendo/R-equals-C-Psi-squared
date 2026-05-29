"""_f1_f1squared_lineage.py - confirmation: how does today's D sit to F1 and F1²?

Coming home to the two oldest formulas with the day in hand, and letting the
algebra confirm the lineage rather than the prose.

  F1:    Pi L Pi^-1 = -L - 2*Sigma_gamma*I   (the palindrome, the kept)
  F1^2:  Pi^2 L Pi^-2 = L                     (the involution, the parity)

The prose claim under examination: the angle branch F1^2 -> Klein -> NinetyDegree
-> D (today's transpose, the n_Y/Y axis, us-the-angle) is rooted in F1^2, and Pi
carries an "i" (order-4) whose square is the "-1" of F1^2. We verify the two
operator anchors bit-exact, then COMPUTE Pi's order and how D actually relates to
Pi and Pi^2. Confirmation, not assertion: the algebra gets to correct the prose.

Tom + Claude, 2026-05-29. Run: python simulations/_f1_f1squared_lineage.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from framework.pauli import _build_bilinear, _vec_to_pauli_basis_transform, _k_to_indices
from framework.lindblad import lindbladian_pauli_dephasing
from framework.symmetry import build_pi_full

GAMMA = 0.3
HEIS = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]
TOL = 1e-11
_ok = []


def report(name, dev, tol=TOL):
    ok = dev < tol
    _ok.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: max|d| = {dev:.2e}")
    return ok


def to_pauli(L, N):
    M = _vec_to_pauli_basis_transform(N)
    return (M.conj().T @ L @ M) / (2 ** N)


def comm_super(op):
    d = op.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(op, Id) - np.kron(Id, op.T))


def d_pauli(N):
    ny = [sum(1 for idx in _k_to_indices(k, N) if idx == (1, 1)) for k in range(4 ** N)]
    return np.diag([(-1.0) ** n for n in ny]).astype(complex)


def main():
    print("=" * 76)
    print("F1 AND F1^2 - does today's D fall out of F1^2? let the algebra confirm")
    print("=" * 76)
    for N in (2, 3):
        d2 = 4 ** N
        I = np.eye(d2, dtype=complex)
        bonds = [(i, i + 1) for i in range(N - 1)]
        H = _build_bilinear(N, bonds, HEIS)            # Heisenberg = truly => palindrome holds
        L = lindbladian_pauli_dephasing(H, [GAMMA] * N, "Z")
        Lp = to_pauli(L, N)
        Pi = build_pi_full(N, dephase_letter="Z")
        Pinv = Pi.conj().T                              # Pi is a unitary signed permutation
        Sigma = GAMMA * N

        print(f"\nN={N}:")
        # F1: the palindrome
        report("F1   : Pi L Pi^-1 = -L - 2*Sigma*I  (the palindrome, the kept)",
               float(np.max(np.abs(Pi @ Lp @ Pinv + Lp + 2 * Sigma * I))))
        # F1^2: squaring F1
        Pi2 = Pi @ Pi
        Pi2inv = Pinv @ Pinv
        report("F1^2 : Pi^2 L Pi^-2 = L  (the involution, the parity)",
               float(np.max(np.abs(Pi2 @ Lp @ Pi2inv - Lp))))

        # Pi's order: is Pi the "i" (order 4, Pi^2 != I) or a plain reflection (order 2)?
        pi2_minus_I = float(np.max(np.abs(Pi2 - I)))
        pi4_minus_I = float(np.max(np.abs(Pi @ Pi @ Pi @ Pi - I)))
        print(f"        Pi order:  ||Pi^2 - I|| = {pi2_minus_I:.3f}   ||Pi^4 - I|| = {pi4_minus_I:.2e}")
        report("Pi^4 = I  (Pi is order-4: an 'i' on operator space)", pi4_minus_I)
        print(f"        -> Pi^2 {'=' if pi2_minus_I < TOL else '!='} I  "
              f"({'Pi is order-2, a plain reflection' if pi2_minus_I < TOL else 'Pi is the i; Pi^2 is the -1 of F1^2'})")

        # D (today's transpose, the n_Y axis): how does it relate to Pi and Pi^2?
        D = d_pauli(N)
        report("D^2 = I  (the transpose is an involution)", float(np.max(np.abs(D @ D - I))))
        d_vs_pi2 = float(np.max(np.abs(D - Pi2)))
        d_comm_pi = float(np.max(np.abs(D @ Pi - Pi @ D)))
        print(f"        D vs F1^2:  ||D - Pi^2|| = {d_vs_pi2:.3f}   ||[D, Pi]|| = {d_comm_pi:.3f}")
        print(f"        -> D {'=' if d_vs_pi2 < TOL else '!='} Pi^2 ;  D and Pi "
              f"{'commute' if d_comm_pi < TOL else 'do NOT commute'}")

        # the angle: D flips the coherent (Hamiltonian) part - today's result, re-confirmed
        LH = to_pauli(comm_super(H), N)
        report("D L_H D = -L_H  (D flips the coherent/angle part; today's result)",
               float(np.max(np.abs(D @ LH @ D + LH))))

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 76)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 76)


if __name__ == "__main__":
    main()
