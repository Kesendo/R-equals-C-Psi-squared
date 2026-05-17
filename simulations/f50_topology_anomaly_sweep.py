#!/usr/bin/env python3
"""F50 topology sweep: which (graph, N) combinations satisfy d_real(-2γ) = 2N?

Background (2026-05-17): the F50 typed claim asserts
`d_real(Re = -2γ) = 2N for any connected graph` based on the
[`PROOF_WEIGHT1_DEGENERACY`](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md)
SWAP-invariant construction. This script does the empirical sweep that
surfaced the K_3 N=3 anomaly (2N+2 = 8 instead of 2N = 6) and verified that
no other tested connected graph violates the formula.

Builds the Heisenberg + Z-dephasing Liouvillian at J = γ = 1 for each
(graph, N), eigendecomposes, counts pure-real eigenvalues at Re = -2γ.

Findings (see PROOF_WEIGHT1_DEGENERACY § Appendix for the structural reading):
  - Chain N=2..5 → count = 2N ✓
  - Ring C_n at n=4, 5 → count = 2N ✓
  - Star K_{1,n-1} at n=3, 4, 5 → count = 2N ✓
  - Complete K_n at n=4, 5 → count = 2N ✓
  - Paw / bowtie / book (triangles inside larger graphs at N=4, 5) → count = 2N ✓
  - N=3 K_3 (= ring = triangle = complete on 3 vertices) → count = 8 (2N+2) ✗

The 2 K_3 extras are weight-1 operators commuting with H_K_3 (matrix-commutator
sum cancellation across the three bonds) but NOT with H_chain alone. They
correspond to the 2-dim standard irrep of S_3 = Aut(K_3) acting on the
weight-1 c=1 Pauli-string sector. Adding any external bond breaks the
S_3 symmetry and removes the extras.
"""
from __future__ import annotations

import sys
import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def two_site(N, i, j, a, b):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, a if k == i else (b if k == j else PAULI["I"]))
    return out


def single_site(N, i, op):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, op if k == i else PAULI["I"])
    return out


def heisenberg(N, bonds, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (two_site(N, a, b, PAULI["X"], PAULI["X"])
                        + two_site(N, a, b, PAULI["Y"], PAULI["Y"])
                        + two_site(N, a, b, PAULI["Z"], PAULI["Z"]))
    return H


def liouvillian(H, gamma, N):
    d = 2 ** N
    I = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for l in range(N):
        Z = single_site(N, l, PAULI["Z"])
        L += gamma * (np.kron(Z.T, Z) - np.kron(I, I))
    return L


def count_pure_real_at(L, target=-2.0, tol=1e-6):
    eigs = np.linalg.eigvals(L)
    return sum(1 for e in eigs if abs(e.real - target) < tol and abs(e.imag) < tol)


def sweep_section(title, graphs, N, expected_2N):
    print(f"--- {title} (expected {expected_2N}) ---")
    print(f"  {'graph':>30} {'bonds':>6} {'count':>6} {'diff':>6}")
    for name, bonds in graphs.items():
        H = heisenberg(N, bonds)
        L = liouvillian(H, 1.0, N)
        cnt = count_pure_real_at(L)
        diff = cnt - expected_2N
        mark = "OK" if diff == 0 else f"+{diff}"
        print(f"  {name:>30} {len(bonds):>6} {cnt:>6} {mark:>6}")
    print()


def main():
    print("=" * 70)
    print("F50 topology sweep (Heisenberg + Z-deph at J = γ = 1)")
    print("=" * 70)
    print()

    sweep_section("N = 3 connected graphs", {
        "P_3 (chain)": [(0, 1), (1, 2)],
        "K_3 (= ring = triangle = complete)": [(0, 1), (1, 2), (2, 0)],
    }, N=3, expected_2N=6)

    sweep_section("N = 4 connected graphs", {
        "P_4 (chain)":          [(0, 1), (1, 2), (2, 3)],
        "K_{1,3} (star)":        [(0, 1), (0, 2), (0, 3)],
        "C_4 (ring)":            [(0, 1), (1, 2), (2, 3), (3, 0)],
        "paw (triangle + leaf)": [(0, 1), (1, 2), (2, 0), (2, 3)],
        "K_4 - e (diamond)":     [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)],
        "K_4 (complete)":        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    }, N=4, expected_2N=8)

    sweep_section("N = 5 connected graphs", {
        "P_5 (chain)":               [(0, 1), (1, 2), (2, 3), (3, 4)],
        "K_{1,4} (star)":            [(0, 1), (0, 2), (0, 3), (0, 4)],
        "C_5 (ring)":                [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
        "K_5 (complete)":            [(i, j) for i in range(5) for j in range(i + 1, 5)],
        "bowtie (two triangles)":    [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)],
        "book (three triangles)":    [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (0, 4), (1, 4)],
    }, N=5, expected_2N=10)

    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print("  Anomaly: N = 3 K_3 only (count = 8 vs expected 2N = 6).")
    print("  All other tested (graph, N) combinations: count = 2N.")
    print("  See PROOF_WEIGHT1_DEGENERACY § Appendix for the structural reading.")


if __name__ == "__main__":
    main()
