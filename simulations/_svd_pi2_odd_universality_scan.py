"""Systematically scan Π²-odd 2-body universality.

Scope:
1. All Π²-odd 2-body bond bilinears at chain N=3, 4, 5.
2. All Π²-even non-truly 2-body bilinears (for contrast).
3. Different topologies: chain, star.
4. Look for unitary equivalence M_XY ↔ M_XZ.

Goal: characterize the universality precisely. Is it strict (all Π²-odd
identical) or does it hold within sub-classes? Is it topology-generic?
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np
from itertools import product as iproduct

import framework as fw
from framework.pauli import _build_bilinear, _k_to_indices, bit_b
from framework.symmetry import pi_squared_eigenvalue
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]
def star_bonds(N): return [(0, i) for i in range(1, N)]


def build_M(N, terms, bonds):
    H = _build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])
    L = lindbladian_z_dephasing(H, [1.0]*N)
    return palindrome_residual(L, N*1.0, N)


def cluster_svs(M, tol=1e-6):
    svs = np.linalg.svd(M, compute_uv=False)
    out = []
    for s in svs:
        placed = False
        for i, (v, c) in enumerate(out):
            if abs(s - v) < tol:
                out[i] = (v, c+1); placed = True; break
        if not placed:
            out.append((s, 1))
    return sorted(out, key=lambda x: -x[0])


def cluster_signature(clusters, tol=1e-4):
    """Make a hashable signature for cluster comparison."""
    return tuple((round(v, 4), m) for v, m in clusters)


def char_letter(idx):
    return ['I', 'X', 'Y', 'Z'][idx]


def pi2_parity(P, Q):
    return (bit_b(P) + bit_b(Q)) % 2


def is_truly(N, terms, bonds):
    """A Hamiltonian is truly if M = 0 exactly."""
    M = build_M(N, terms, bonds)
    return float(np.linalg.norm(M)) < 1e-9


def main():
    print("=" * 78)
    print("Scan 1: All Π²-odd 2-body single bilinears (chain)")
    print("=" * 78)
    print()

    paulis = ['I', 'X', 'Y', 'Z']

    print("Chain N=4:")
    bonds_c4 = chain_bonds(4)
    pi2_odd_terms = []
    for p in paulis:
        for q in paulis:
            if (p == 'I' and q == 'I'): continue
            # Pure 2-body: both non-identity
            if p == 'I' or q == 'I': continue
            if pi2_parity(p, q) == 1:
                pi2_odd_terms.append((p, q))

    sigs_at_n4 = {}
    for term in pi2_odd_terms:
        M = build_M(4, [term], bonds_c4)
        sig = cluster_signature(cluster_svs(M))
        sigs_at_n4.setdefault(sig, []).append(term)

    for sig, terms in sigs_at_n4.items():
        print(f"  Cluster signature {sig}:")
        for t in terms:
            print(f"    {t[0]}{t[1]}")
    print()

    print("Chain N=5:")
    bonds_c5 = chain_bonds(5)
    sigs_at_n5 = {}
    for term in pi2_odd_terms:
        M = build_M(5, [term], bonds_c5)
        sig = cluster_signature(cluster_svs(M))
        sigs_at_n5.setdefault(sig, []).append(term)

    for sig, terms in sigs_at_n5.items():
        print(f"  Cluster signature {sig}:")
        for t in terms:
            print(f"    {t[0]}{t[1]}")
    print()

    print("=" * 78)
    print("Scan 2: All Π²-even non-truly 2-body bilinears (chain)")
    print("=" * 78)
    print()

    pi2_even_terms = []
    for p in paulis:
        for q in paulis:
            if p == 'I' or q == 'I': continue
            if pi2_parity(p, q) == 0:
                pi2_even_terms.append((p, q))

    print("Chain N=4:")
    sigs_even_n4 = {}
    truly_terms_n4 = []
    for term in pi2_even_terms:
        truly = is_truly(4, [term], bonds_c4)
        if truly:
            truly_terms_n4.append(term)
        else:
            M = build_M(4, [term], bonds_c4)
            sig = cluster_signature(cluster_svs(M))
            sigs_even_n4.setdefault(sig, []).append(term)

    print(f"  Truly (M=0): {[t[0]+t[1] for t in truly_terms_n4]}")
    for sig, terms in sigs_even_n4.items():
        print(f"  Cluster signature {sig}:")
        for t in terms:
            print(f"    {t[0]}{t[1]}")
    print()

    print("=" * 78)
    print("Scan 3: Star topology Π²-odd at N=4")
    print("=" * 78)
    print()

    bonds_s4 = star_bonds(4)
    sigs_star_n4 = {}
    for term in pi2_odd_terms:
        M = build_M(4, [term], bonds_s4)
        sig = cluster_signature(cluster_svs(M))
        sigs_star_n4.setdefault(sig, []).append(term)

    for sig, terms in sigs_star_n4.items():
        print(f"  Cluster signature {sig}:")
        for t in terms:
            print(f"    {t[0]}{t[1]}")
    print()

    print("=" * 78)
    print("Scan 4: Mixed-parity Hamiltonians (chain N=4)")
    print("=" * 78)
    print()

    mixed_cases = [
        ("YZ + XY (even + odd)",  [('Y', 'Z'), ('X', 'Y')]),
        ("ZY + XZ (even + odd)",  [('Z', 'Y'), ('X', 'Z')]),
        ("YZ+ZY + XY+YX (mix)",   [('Y', 'Z'), ('Z', 'Y'), ('X', 'Y'), ('Y', 'X')]),
    ]
    for name, terms in mixed_cases:
        M = build_M(4, terms, bonds_c4)
        idx_p = np.array([k for k in range(4**4) if pi_squared_eigenvalue(_k_to_indices(k, 4)) == 1])
        idx_m = np.array([k for k in range(4**4) if pi_squared_eigenvalue(_k_to_indices(k, 4)) == -1])
        diag_pp = float(np.linalg.norm(M[np.ix_(idx_p, idx_p)]))
        diag_mm = float(np.linalg.norm(M[np.ix_(idx_m, idx_m)]))
        off_pm = float(np.linalg.norm(M[np.ix_(idx_p, idx_m)]))
        cls = cluster_svs(M)
        print(f"  {name}: ‖[+,+]‖={diag_pp:.3f}, ‖[-,-]‖={diag_mm:.3f}, ‖[+,-]‖={off_pm:.3f}")
        print(f"    clusters: {[(round(v,3), m) for v, m in cls[:8]]} ...")


if __name__ == '__main__':
    main()
