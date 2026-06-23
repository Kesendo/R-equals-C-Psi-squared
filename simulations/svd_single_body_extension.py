"""Extend the additive single-site decomposition to other single-body Pauli's.

Single-body Pauli letters: X, Y, Z (I trivially gives 0).
Question: For H = Σ_l c_l·P_l (P ∈ {X,Y,Z}), is M = Σ_l M_l with M_l normal?

The N=3,4,5 test for Y was a perfect match. This script checks X and Z.
If both match too, the additive single-site decomposition is universal for
single-body Pauli Hamiltonians (under Z-dephasing).
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np

import framework as fw
from framework.pauli import site_op
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
from itertools import product as iproduct


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]


def site_degrees_chain(N):
    if N == 1: return [0]
    return [1] + [2]*(N-2) + [1]


def build_H_bondsum_single_body(N, P='Y', bonds=None):
    """H = Σ_bonds (I_b·P_b' + P_b·I_b'). For chain → Σ_l deg(l)·P_l."""
    if bonds is None: bonds = chain_bonds(N)
    H = np.zeros((2**N, 2**N), dtype=complex)
    for (i, j) in bonds:
        H = H + site_op(N, j, P) + site_op(N, i, P)
    return H


def build_M(N, H, gamma=1.0):
    L = lindbladian_z_dephasing(H, [gamma]*N)
    return palindrome_residual(L, gamma*N, N)


def build_single_site_L(c, P, gamma=1.0):
    """L_l for site with c·P_l Hamiltonian + γ·Z-dephasing.

    Pauli basis (I, X, Y, Z). Returns 4x4 L_l.
    Hamiltonian commutator -i[c·P, σ] uses Pauli structure constants.
    Dephasing is diagonal with -2γ on X, Y rows; 0 on I, Z rows.
    """
    # Hamiltonian commutator [P, σ]: structure constants
    #   [X, X]=0, [X, Y]=2iZ, [X, Z]=-2iY, [X, I]=0
    #   [Y, X]=-2iZ, [Y, Y]=0, [Y, Z]=2iX, [Y, I]=0
    #   [Z, X]=2iY, [Z, Y]=-2iX, [Z, Z]=0, [Z, I]=0
    # -i[P, σ] gives:
    #   X: I→0, X→0, Y→2Z, Z→-2Y
    #   Y: I→0, X→-2Z, Y→0, Z→2X
    #   Z: I→0, X→2Y, Y→-2X, Z→0
    HamCom = {
        'X': np.array([[0,0,0,0],[0,0,0,0],[0,0,0,-2],[0,0,2,0]], dtype=complex),
        'Y': np.array([[0,0,0,0],[0,0,0,-2],[0,0,0,0],[0,2,0,0]], dtype=complex),
        'Z': np.array([[0,0,0,0],[0,0,-2,0],[0,2,0,0],[0,0,0,0]], dtype=complex),
    }[P]
    # Z-dephasing: diagonal -2γ on X, Y; 0 on I, Z
    Diss = np.diag([0, -2*gamma, -2*gamma, 0]).astype(complex)
    return c * HamCom + Diss


def build_single_site_M(c, P, gamma=1.0):
    L_l = build_single_site_L(c, P, gamma)
    Pi = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,1j,0]], dtype=complex)
    Pi_inv = np.linalg.inv(Pi)
    return Pi @ L_l @ Pi_inv + L_l + 2*gamma*np.eye(4)


def is_normal(A, tol=1e-10):
    return np.linalg.norm(A @ A.conj().T - A.conj().T @ A) < tol


def cluster_values(values, tol=1e-6):
    out = []
    for v in values:
        placed = False
        for i, (val, count) in enumerate(out):
            if abs(v - val) < tol:
                out[i] = (val, count+1)
                placed = True
                break
        if not placed:
            out.append((v, 1))
    return sorted(out, key=lambda x: -abs(x[0]))


def predict_M_eigenvalues(N, P, gamma=1.0):
    degs = site_degrees_chain(N)
    M_per_site = [build_single_site_M(d, P, gamma) for d in degs]
    eigs = [np.linalg.eigvals(M) for M in M_per_site]
    return np.array([sum(eigs[l][combo[l]] for l in range(N))
                     for combo in iproduct(*[range(4) for _ in range(N)])])


def main():
    print("=" * 78)
    print("Single-site M_l for P ∈ {X, Y, Z} (c=1, γ=1)")
    print("=" * 78)
    print()

    for P in ['X', 'Y', 'Z']:
        print(f"  H = c·{P}_l, γ=1:")
        for c in [1, 2]:
            M_l = build_single_site_M(c, P)
            normal = is_normal(M_l)
            eigs = sorted(np.linalg.eigvals(M_l), key=lambda z: (-abs(z), -np.real(z)))
            svs = sorted(np.linalg.svd(M_l, compute_uv=False), reverse=True)
            eigs_str = [f"{complex(e).real:.2f}{complex(e).imag:+.2f}j" for e in eigs]
            print(f"    c={c}: eigs={eigs_str}, "
                  f"svs={[round(s,3) for s in svs]}, normal={normal}")
        print()

    print("=" * 78)
    print("Test: M (full) clusters match additive single-site prediction?")
    print("=" * 78)
    print()

    for P in ['X', 'Y', 'Z']:
        print(f"  H = bond-sum ({P} on every site, weighted by degree):")
        for N in [3, 4]:
            H = build_H_bondsum_single_body(N, P)
            if not np.allclose(H, H.conj().T):
                print(f"    N={N}: H is NOT Hermitian, skipping")
                continue
            M_full = build_M(N, H)
            direct = cluster_values(np.linalg.svd(M_full, compute_uv=False))
            pred = cluster_values(np.abs(predict_M_eigenvalues(N, P)))
            match = (len(direct) == len(pred)
                     and all(abs(a[0]-b[0]) < 1e-5 and a[1] == b[1]
                             for a, b in zip(direct, pred)))
            print(f"    N={N}: direct={[(round(v,3),m) for v,m in direct]}")
            print(f"           pred  ={[(round(v,3),m) for v,m in pred]}")
            print(f"           match: {'YES' if match else 'NO'}")
        print()


if __name__ == '__main__':
    main()
