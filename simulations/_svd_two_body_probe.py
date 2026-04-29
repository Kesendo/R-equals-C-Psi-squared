"""Quick probe: does 2-body bond-bilinear M have analogous structure?

Single-body H = Σ_l c_l·P_l now has a closed-form additive decomposition for
M's SVD (see _svd_active_spectator.py). For 2-body bond-bilinears
H = Σ_bonds (A_b·B_b' + ...), L doesn't trivially decompose per-site or
per-bond (consecutive bonds share sites in a chain).

This script just *probes* M's cluster structure for representative 2-body
cases (truly, soft, hard), to see if there's an obvious next pattern to chase.
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np

import framework as fw
from framework.pauli import site_op, _build_bilinear
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]


def build_H(N, terms, J=1.0):
    bonds = chain_bonds(N)
    return _build_bilinear(N, bonds, [(t[0], t[1], J) for t in terms])


def cluster_svs(M, tol=1e-6):
    svs = np.linalg.svd(M, compute_uv=False)
    out = []
    for s in svs:
        placed = False
        for i, (v, c) in enumerate(out):
            if abs(s - v) < tol:
                out[i] = (v, c+1)
                placed = True
                break
        if not placed:
            out.append((s, 1))
    return sorted(out, key=lambda x: -x[0])


def main():
    cases = [
        ("Heisenberg (truly)",   [('X','X'), ('Y','Y'), ('Z','Z')]),
        ("YZ+ZY (soft 2-body)",  [('Y','Z'), ('Z','Y')]),
        ("XX+XY (hard 2-body)",  [('X','X'), ('X','Y')]),
        ("XX (single 2-body)",   [('X','X')]),
        ("YY (single 2-body)",   [('Y','Y')]),
        ("ZZ (single 2-body)",   [('Z','Z')]),
        ("XY+YX (asymmetric)",   [('X','Y'), ('Y','X')]),
        ("YZ (single 2-body)",   [('Y','Z')]),
    ]
    for label, terms in cases:
        print(f"\n=== {label}: terms = {terms} ===")
        for N in [3, 4]:
            chain = fw.ChainSystem(N=N, J=1.0)
            cls = chain.classify_pauli_pair(terms)
            H = build_H(N, terms)
            L = lindbladian_z_dephasing(H, [1.0]*N)
            M = palindrome_residual(L, N*1.0, N)
            clusters = cluster_svs(M)
            norm_sq = float(np.linalg.norm(M)**2)
            print(f"  N={N} ({cls}): ‖M‖²={norm_sq:.0f}, clusters={[(round(v,3), m) for v, m in clusters]}")


if __name__ == '__main__':
    main()
