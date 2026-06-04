#!/usr/bin/env python3
"""The ring dihedral-lock constant c_N: does N->inf give 1/sqrt(2)? No, it gives ln 2.

PROOF_RING_N4_DIHEDRAL_LOCK shows Im_max(L) = DeltaE_max(H), the max eigenvalue GAP of the
isotropic Heisenberg Hamiltonian (the Z-dephasing adds only real decay, never touches Im). So the
dihedral-lock constant

    c_N = Im_max / (J*N) = (E_max - E_min) / (J*N) = 1/4 - E0/(J*N)

is a pure HAMILTONIAN question (2^N), not the full Liouvillian (4^N). E_max = J*N/4 is the
ferromagnet; E0 is the antiferromagnetic Heisenberg-ring ground state. In the limit the per-bond
ground energy is the Hulthen/Bethe value E0/(J*N) -> 1/4 - ln2, so

    c_inf = 1/4 - (1/4 - ln2) = ln 2 = 0.693147...   (NOT 1/sqrt(2) = 0.707107)

The finite-N sequence c_4=3/4, c_6=0.7171, c_8=0.7064 DECREASES through 1/sqrt(2) toward ln2.
1/sqrt(2) is only a value it passes through (the same red-herring lesson as s*=0.709). We confirm
by computing the AFM Heisenberg ring ground state E0(N) up to N=16 (sparse, S_z=0 ground state)
and forming c_N.
"""
from __future__ import annotations

import sys
import math

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def heisenberg_ring_ground_energy(N: int) -> float:
    """E0 of H = (1/4) sum_bonds (X_iX_j + Y_iY_j + Z_iZ_j) on the N-site ring (J=1).
    Built sparsely by bit manipulation over the full 2^N space; ZZ diagonal, XY hops a
    differing pair with matrix element 1/2. eigsh finds the smallest (AFM ground) eigenvalue."""
    dim = 1 << N
    bonds = [(i, (i + 1) % N) for i in range(N)]
    rows, cols, vals = [], [], []
    diag = np.zeros(dim)
    for s in range(dim):
        d = 0.0
        for (i, j) in bonds:
            bi = (s >> i) & 1
            bj = (s >> j) & 1
            d += 0.25 if bi == bj else -0.25            # (1/4) Z_i Z_j
            if bi != bj:                                # (1/4)(X_iX_j+Y_iY_j) hops the differing pair
                t = s ^ (1 << i) ^ (1 << j)
                rows.append(t); cols.append(s); vals.append(0.5)
        diag[s] = d
    H = sp.coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()
    H = H + sp.diags(diag)
    e0 = eigsh(H, k=1, which="SA", return_eigenvectors=False, maxiter=5000)
    return float(e0[0])


def main() -> None:
    LN2 = math.log(2.0)
    INV_SQRT2 = 1.0 / math.sqrt(2.0)
    print("=" * 78)
    print("  Ring dihedral-lock constant c_N = 1/4 - E0/(J*N),  J = 1")
    print(f"  candidates: ln2 = {LN2:.6f}   1/sqrt(2) = {INV_SQRT2:.6f}")
    print("=" * 78)
    print(f"  {'N':>3} {'E0':>11} {'E0/N':>10} {'c_N':>10} {'c_N - ln2':>11} {'c_N - 1/√2':>11}")
    prev = None
    for N in range(4, 17, 2):
        e0 = heisenberg_ring_ground_energy(N)
        cN = 0.25 - e0 / N
        cross = ""
        if prev is not None and (prev - INV_SQRT2) * (cN - INV_SQRT2) < 0:
            cross = "  <- crosses 1/√2 here"
        print(f"  {N:>3} {e0:>11.5f} {e0/N:>10.5f} {cN:>10.5f} {cN-LN2:>+11.5f} {cN-INV_SQRT2:>+11.5f}{cross}")
        prev = cN
    print()
    print("  reading: c_N decreases monotonically and converges to ln2 = 0.6931 (the Hulthen")
    print("  per-bond AFM ground energy 1/4 - ln2), passing THROUGH 1/sqrt(2)=0.7071 on the way.")
    print("  The N->inf limit is ln2, not 1/sqrt(2); 1/sqrt(2) is a value it crosses, like s*.")


if __name__ == "__main__":
    main()
