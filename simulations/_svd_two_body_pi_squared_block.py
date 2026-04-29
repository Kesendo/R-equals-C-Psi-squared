"""Test: do XX+XY hard's two clusters correspond to Π² = ±1 eigenspaces?

Observation: for chain XX+XY at N≥4, M has exactly 2 distinct SVs with equal
multiplicities 4^N/2 each. The Π² eigenspaces V_+ (Π²=+1) and V_- (Π²=-1)
have dimensions 2^(2N-1) each = 4^N/2.

Hypothesis: M restricted to V_+ has uniform SV = σ_+; restricted to V_-
has uniform SV = σ_-. For specific Hamiltonians these σ's may coincide
(N=3 shows uniform 2√2 — accidental collision). At N≥4 they differ.

This script verifies the block-uniform-SV-per-Π²-eigenspace hypothesis.
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np

import framework as fw
from framework.pauli import _build_bilinear, _k_to_indices
from framework.symmetry import pi_squared_eigenvalue
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]


def build_M(N, terms, bonds):
    bilinear = [(t[0], t[1], 1.0) for t in terms]
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [1.0]*N)
    return palindrome_residual(L, N*1.0, N)


def pi_squared_split(N):
    """Returns (idx_plus, idx_minus): row/col indices of M for Π²=+1 and -1."""
    idx_plus = []
    idx_minus = []
    for k in range(4**N):
        eig = pi_squared_eigenvalue(_k_to_indices(k, N))
        if eig == 1:
            idx_plus.append(k)
        else:
            idx_minus.append(k)
    return np.array(idx_plus), np.array(idx_minus)


def block_svs(M, idx):
    block = M[np.ix_(idx, idx)]
    return np.linalg.svd(block, compute_uv=False)


def cluster_svs(svs, tol=1e-6):
    out = []
    for s in svs:
        placed = False
        for i, (v, c) in enumerate(out):
            if abs(s - v) < tol:
                out[i] = (v, c+1); placed = True; break
        if not placed:
            out.append((s, 1))
    return sorted(out, key=lambda x: -x[0])


def test_terms_label(label, terms):
    print(f"  Terms: {label} = {terms}")
    for N in [3, 4, 5]:
        bonds = chain_bonds(N)
        M = build_M(N, terms, bonds)
        idx_p, idx_m = pi_squared_split(N)

        # Check that Π² block is sufficient (M is block-diagonal under Π²?)
        # Off-diagonal: M[idx_p, idx_m] and M[idx_m, idx_p].
        off_pm = float(np.linalg.norm(M[np.ix_(idx_p, idx_m)]))
        off_mp = float(np.linalg.norm(M[np.ix_(idx_m, idx_p)]))

        svs_p = block_svs(M, idx_p)
        svs_m = block_svs(M, idx_m)
        cl_p = cluster_svs(svs_p)
        cl_m = cluster_svs(svs_m)

        full_svs = np.linalg.svd(M, compute_uv=False)
        cl_full = cluster_svs(full_svs)

        print(f"    N={N}: dim V_+={len(idx_p)}, dim V_-={len(idx_m)}")
        print(f"      Off-diagonal blocks: ‖M[+,-]‖={off_pm:.6f}, ‖M[-,+]‖={off_mp:.6f}"
              f" → block-diag={'YES' if off_pm < 1e-9 and off_mp < 1e-9 else 'NO'}")
        print(f"      V_+ block clusters: {[(round(v,4), m) for v, m in cl_p]}")
        print(f"      V_- block clusters: {[(round(v,4), m) for v, m in cl_m]}")
        print(f"      Full M clusters:    {[(round(v,4), m) for v, m in cl_full]}")
        # Combine: block clusters merged
        combined = {}
        for v, m in cl_p + cl_m:
            for k in list(combined.keys()):
                if abs(k - v) < 1e-6:
                    combined[k] += m; break
            else:
                combined[v] = m
        cl_combined = sorted(combined.items(), key=lambda x: -x[0])
        match = (len(cl_combined) == len(cl_full)
                 and all(abs(a[0]-b[0]) < 1e-5 and a[1] == b[1]
                         for a, b in zip(cl_combined, cl_full)))
        print(f"      Block ⊕ block = full M? {'YES' if match else 'NO'}")
    print()


def main():
    print("=" * 78)
    print("Π²-eigenspace block decomposition of M for 2-body bilinears")
    print("=" * 78)
    print()

    cases = [
        ("XX+XY (hard, max-uniform-ish)",  [('X','X'), ('X','Y')]),
        ("YZ (soft, single-term)",         [('Y','Z')]),
        ("YZ+ZY (soft, canonical)",        [('Y','Z'), ('Z','Y')]),
        ("XX (truly)",                     [('X','X')]),
        ("Heisenberg XX+YY+ZZ (truly)",    [('X','X'),('Y','Y'),('Z','Z')]),
    ]
    for label, terms in cases:
        test_terms_label(label, terms)


if __name__ == '__main__':
    main()
