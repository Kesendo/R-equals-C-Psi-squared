#!/usr/bin/env python3
"""eq014_spectrum_check.py

Verify the C# eigendecomposition data (eq014_eigvals_n7.bin) matches one
of two candidate Hamiltonians at N=7 uniform chain, γ_0 = 0.05:
 (A) Heisenberg  H = Σ_i J (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1})
     (= C# Topology.Chain with PauliTypes=["X","Y","Z"])
 (B) XY (PTF)    H = Σ_i (J/2) (X_i X_{i+1} + Y_i Y_{i+1})
     (= PTF convention, target of the task)

For each candidate we build L_A as a sparse d^2 x d^2 matrix, compute a
small set of slow eigenvalues via ARPACK shift-invert at sigma=0, and
compare to the corresponding slow eigenvalues from the C# run.

If C# matches (A): the engine is wired to Heisenberg, not XY.
If C# matches (B): the engine is wired to XY, and the α_i mismatch is
somewhere else.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

sys.stdout.reconfigure(encoding='utf-8')

RESULTS_DIR = Path(__file__).parent / "results"
N = 7
GAMMA = 0.05
J = 1.0
D = 2 ** N
DD = D * D

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def site_op_sparse(op, site, N):
    op_s = sps.csr_matrix(op, dtype=complex)
    full = sps.eye(1, dtype=complex, format='csr')
    I2_s = sps.eye(2, dtype=complex, format='csr')
    for k in range(N):
        full = sps.kron(full, op_s if k == site else I2_s, format='csr')
    return full


def build_H_heisenberg(J, N):
    """H = J Σ_i (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1})"""
    H = sps.csr_matrix((D, D), dtype=complex)
    for i in range(N - 1):
        for op in (X, Y, Z):
            H = H + J * site_op_sparse(op, i, N) @ site_op_sparse(op, i + 1, N)
    return H.tocsr()


def build_H_xy_ptf(J, N):
    """H = (J/2) Σ_i (X_i X_{i+1} + Y_i Y_{i+1})"""
    H = sps.csr_matrix((D, D), dtype=complex)
    for i in range(N - 1):
        for op in (X, Y):
            H = H + (J / 2) * site_op_sparse(op, i, N) @ site_op_sparse(op, i + 1, N)
    return H.tocsr()


def build_L(H, N, gamma):
    """L = -i (H⊗I - I⊗H^T) + dephasing(γ).
    vec(ρ)[a*d + b] = ρ[a, b]. Dephasing diagonal: -2γ * popcount(a ^ b) in
    big-endian convention (site k at bit N-1-k)."""
    Id = sps.eye(D, dtype=complex, format='csr')
    L_h = -1j * (sps.kron(H, Id, format='csr') - sps.kron(Id, H.T, format='csr'))
    a_idx = np.arange(DD) // D
    b_idx = np.arange(DD) % D
    h_vec = np.array([bin(int(av ^ bv)).count('1') for av, bv in zip(a_idx, b_idx)])
    deph_diag = -2.0 * gamma * h_vec
    L_d = sps.diags(deph_diag.astype(complex), format='csr')
    return (L_h + L_d).tocsr()


def main():
    # Load C# eigenvalues
    csharp = np.fromfile(RESULTS_DIR / "eq014_eigvals_n7.bin", dtype=np.complex128)
    print(f"C# eigenvalues loaded: {len(csharp)}")
    print(f"  Re range: [{csharp.real.min():.4f}, {csharp.real.max():.4f}]")
    print(f"  |Im| max: {np.abs(csharp.imag).max():.4f}")
    print(f"  Stationary (|λ|<1e-10): {np.sum(np.abs(csharp) < 1e-10)}")
    print()

    # Sort C# slow modes (smallest |Re|)
    slow_csharp = np.sort(csharp[np.abs(csharp.real) < 0.2])
    # Sort by (Re, Im) for deterministic comparison
    slow_csharp = slow_csharp[np.lexsort((slow_csharp.imag, slow_csharp.real))]
    print(f"C# slow modes (|Re|<0.2): {len(slow_csharp)}")
    print(f"  First 10: {slow_csharp[:10]}")
    print()

    # Candidate spectra: compare via H's own spectrum (max - min eigenvalue).
    # For L = -i(H⊗I - I⊗H^T): max |Im(λ_L)| = E_max(H) - E_min(H).
    print("=== Candidate A: Heisenberg H = Σ (XX+YY+ZZ) with J=1 ===")
    H_heis = build_H_heisenberg(J, N).toarray()
    E_heis = np.sort(np.linalg.eigvalsh(H_heis))
    print(f"  H_heis spectrum: E_min={E_heis[0]:.4f}, E_max={E_heis[-1]:.4f}")
    print(f"  E_max - E_min = {E_heis[-1] - E_heis[0]:.4f}")
    print(f"  Expected |Im(λ_L)| max = {E_heis[-1] - E_heis[0]:.4f}")

    print("\n=== Candidate B: XY PTF H = Σ (1/2)(XX+YY) with J=1 ===")
    H_xy = build_H_xy_ptf(J, N).toarray()
    E_xy = np.sort(np.linalg.eigvalsh(H_xy))
    print(f"  H_xy spectrum: E_min={E_xy[0]:.4f}, E_max={E_xy[-1]:.4f}")
    print(f"  E_max - E_min = {E_xy[-1] - E_xy[0]:.4f}")
    print(f"  Expected |Im(λ_L)| max = {E_xy[-1] - E_xy[0]:.4f}")

    print("\n=== C# observed |Im(λ)| max = 17.3450 ===")
    diff_heis_match = abs((E_heis[-1] - E_heis[0]) - 17.3450)
    diff_xy_match = abs((E_xy[-1] - E_xy[0]) - 17.3450)
    print(f"  Heisenberg match: |observed - predicted| = {diff_heis_match:.3f}")
    print(f"  XY-PTF match:     |observed - predicted| = {diff_xy_match:.3f}")

    print("\n=== Verdict ===")
    if diff_heis_match < 0.01:
        print("C# data consistent with HEISENBERG (XX+YY+ZZ).")
        print("--> Need to build XY-only path for PTF analysis.")
    elif diff_xy_match < 0.01:
        print("C# data consistent with XY-PTF ((1/2)(XX+YY)).")
        print("--> Hamiltonian is correct; investigate α_i mismatch elsewhere.")
    else:
        print(f"C# data matches NEITHER clean candidate:")
        print(f"  Heisenberg predicts {E_heis[-1] - E_heis[0]:.4f}")
        print(f"  XY-PTF predicts    {E_xy[-1] - E_xy[0]:.4f}")
        print(f"  Observed           17.3450")


if __name__ == "__main__":
    main()
