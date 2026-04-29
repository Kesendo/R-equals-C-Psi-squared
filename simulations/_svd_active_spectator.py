"""Test: do M's SVD clusters factor as a sum of single-site M_l operators?

Hypothesis: For bond-summed single-body H = Σ_bonds (I_b·Y_b' + Y_b·I_b') on a
chain (or any topology), H = Σ_l c_l·Y_l where c_l = deg(l). Then:
  L = -i[H, ·] + Σ_l γ_l·(Z_l ρ Z_l - ρ) = Σ_l L_l
where L_l acts only on site l. Likewise Π factors per site, so M factors as
  M = Σ_l M_l ⊗ I_(others)
with M_l = T_Π,l·L_l·T_Π,l⁻¹ + L_l + 2γ_l·I_l, a 4×4 matrix per site.

Eigenvalues of M = sums Σ_l λ_l of single-site eigenvalues.
SVs of M = sums of single-site SVs IF each M_l is normal.

If this prediction matches direct computation → cluster structure is just
the multiplicity-arithmetic of summing single-site spectra.

Run from repo root: python simulations/_svd_active_spectator.py
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np

import framework as fw
from framework.pauli import site_op
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
from itertools import product as iproduct


def chain_bonds(N: int):
    return [(i, i + 1) for i in range(N - 1)]


def star_bonds(N: int):
    """Star: hub at site 0, leaves at sites 1..N-1."""
    return [(0, i) for i in range(1, N)]


def complete_bonds(N: int):
    return [(i, j) for i in range(N) for j in range(i+1, N)]


def topology_degrees(N: int, topology: str) -> list[int]:
    """Site degrees for various topologies."""
    if topology == 'chain':
        if N == 1: return [0]
        return [1] + [2]*(N-2) + [1]
    if topology == 'star':
        return [N-1] + [1]*(N-1)
    if topology == 'complete':
        return [N-1]*N
    raise ValueError(f"unknown topology {topology}")


def build_H_bondsum_iyyi(N: int, bonds=None):
    """H from bond-summed (I,Y) + (Y,I) — i.e., Σ_bonds (I_b·Y_b' + Y_b·I_b').

    Equivalent to Σ_l deg(l)·Y_l (single-body Y on every site, weighted by
    site degree).
    """
    if bonds is None:
        bonds = chain_bonds(N)
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (i, j) in bonds:
        # I_i·Y_j + Y_i·I_j = Y on site j + Y on site i
        H = H + site_op(N, j, 'Y') + site_op(N, i, 'Y')
    return H


def site_degrees_chain(N: int) -> list[int]:
    return topology_degrees(N, 'chain')


def build_M(N: int, H, gamma: float = 1.0) -> np.ndarray:
    L = lindbladian_z_dephasing(H, [gamma] * N)
    return palindrome_residual(L, gamma * N, N)


# ----------------------------------------------------------------------
# Single-site M_l in the per-site Pauli basis (I, X, Y, Z)
# ----------------------------------------------------------------------

def build_single_site_M(c: float, gamma: float = 1.0) -> np.ndarray:
    """M_l for site with Hamiltonian c·Y_l, dephasing γ.

    Per-site Pauli basis = (I, X, Y, Z). Returns the 4×4 matrix
      M_l = T_Π,l · L_l · T_Π,l⁻¹ + L_l + 2γ_l · I_l
    where L_l includes both the site's Hamiltonian commutator and Z-dephasing.
    """
    # L_l in basis (I, X, Y, Z):
    #   -i[c·Y, ·]: X→2c·Z, Z→-2c·X
    #   Z-dephasing: -2γ on X, -2γ on Y, 0 on I, Z
    L_l = np.array([
        [0,    0,    0,    0   ],
        [0,   -2*gamma, 0,   -2*c ],
        [0,    0,   -2*gamma, 0  ],
        [0,    2*c,  0,    0   ],
    ], dtype=complex)

    # Π per site: I↔X (phase 1), Y→iZ, Z→iY
    Pi = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1j],
        [0, 0, 1j, 0],
    ], dtype=complex)
    Pi_inv = np.linalg.inv(Pi)

    return Pi @ L_l @ Pi_inv + L_l + 2 * gamma * np.eye(4)


def is_normal(A: np.ndarray, tol: float = 1e-10) -> bool:
    return np.linalg.norm(A @ A.conj().T - A.conj().T @ A) < tol


def cluster_values(values, tol: float = 1e-6):
    out = []
    for v in values:
        placed = False
        for i, (val, count) in enumerate(out):
            if abs(v - val) < tol:
                out[i] = (val, count + 1)
                placed = True
                break
        if not placed:
            out.append((v, 1))
    return sorted(out, key=lambda x: -abs(x[0]))


def cluster_svs(M: np.ndarray, tol: float = 1e-6):
    return cluster_values(np.linalg.svd(M, compute_uv=False), tol)


def predict_M_eigenvalues(N: int, gamma: float = 1.0, topology: str = 'chain'):
    """Predict eigenvalues of M = Σ_l M_l ⊗ I_others as sums of M_l ev's."""
    degs = topology_degrees(N, topology)
    M_per_site = [build_single_site_M(d, gamma) for d in degs]
    eigs_per_site = [np.linalg.eigvals(M) for M in M_per_site]
    out = []
    for combo in iproduct(*[range(4) for _ in range(N)]):
        ev = sum(eigs_per_site[l][combo[l]] for l in range(N))
        out.append(ev)
    return np.array(out)


def predict_M_svs_via_normality(N: int, gamma: float = 1.0, topology: str = 'chain'):
    return np.abs(predict_M_eigenvalues(N, gamma, topology))


def main():
    print("=" * 78)
    print("Single-site M_l structure (one site, c·Y, dephasing γ=1)")
    print("=" * 78)
    print()

    for c in [1, 2]:
        M_l = build_single_site_M(c, gamma=1.0)
        normal = is_normal(M_l)
        eigs = np.linalg.eigvals(M_l)
        svs = np.linalg.svd(M_l, compute_uv=False)
        print(f"  c={c}: M_l 4×4")
        print(f"    eigenvalues: {sorted([complex(e) for e in eigs], key=lambda z: -abs(z))}")
        print(f"    singular values: {sorted(svs, reverse=True)}")
        print(f"    |eigenvalues|: {sorted([abs(e) for e in eigs], reverse=True)}")
        print(f"    normal? {normal}")
        print()

    print("=" * 78)
    print("Test: M's SVD = |sum of single-site eigenvalues|?")
    print("=" * 78)
    print()

    for topology in ['chain', 'star', 'complete']:
        print(f"--- topology = {topology} ---")
        for N in [3, 4, 5]:
            degs = topology_degrees(N, topology)
            if topology == 'chain':
                bonds = chain_bonds(N)
            elif topology == 'star':
                bonds = star_bonds(N)
            else:
                bonds = complete_bonds(N)
            H = build_H_bondsum_iyyi(N, bonds=bonds)
            M_full = build_M(N, H)
            direct_norm_sq = float(np.linalg.norm(M_full) ** 2)
            direct_clusters = cluster_svs(M_full, tol=1e-6)

            pred_svs = predict_M_svs_via_normality(N, topology=topology)
            pred_clusters = cluster_values(pred_svs, tol=1e-6)
            pred_norm_sq = float(np.sum(pred_svs ** 2))

            match = (len(direct_clusters) == len(pred_clusters)
                     and all(abs(d[0]-p[0]) < 1e-5 and d[1] == p[1]
                             for d, p in zip(direct_clusters, pred_clusters)))
            print(f"  N={N} degs={degs}: ‖M‖²={direct_norm_sq:.0f} (pred={pred_norm_sq:.0f}); "
                  f"clusters={[(round(v,3), m) for v, m in direct_clusters]} "
                  f"match={'YES' if match else 'NO'}")
        print()


if __name__ == '__main__':
    main()
