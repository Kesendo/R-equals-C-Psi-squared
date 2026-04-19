#!/usr/bin/env python3
"""Follow-up: direct computation of stationary-subspace residues in A and B.

If the excitation-sector projectors P_n (n=0..N) are the strict stationary
modes of L (F4 claim), they are J-independent because [H_XY, P_n] = 0
irrespective of the coupling pattern. Then the z->0 residue of G(z) is
J-invariant, and Sum_i ln|<P_i R_stat>| is trivially zero site-by-site.

This isolates whether the resolvent z->0 closure is the TRIVIAL sector
result or a non-trivial statement about eigenvectors beyond the F4 cluster.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, DEFECT_BOND,
    X, Y, Z, I2, site_op,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, bonding_plus_vacuum,
    density_matrix,
)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def compute_stationary_projector(L, tol=1e-10):
    """Return P_stat = projector onto null(L) in operator-vector space, and
    the list of (right, left) eigenvector pairs for strictly stationary modes."""
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    mask = np.abs(eigvals) < tol
    idx = np.where(mask)[0]
    # Build spectral projector sum_s |M_s><W_s|
    P = np.zeros_like(L)
    for s in idx:
        P += np.outer(V_R[:, s], V_Linv[s, :])
    return P, idx, eigvals


def site_pauli_matrix(i, local, N):
    """Return a flat (d^2,) vector for site-i Pauli <local> tensor I elsewhere."""
    P = site_op(local, i, N)
    return P.flatten(order='F')


def main():
    for N in [3, 5]:
        print(f"\n{'='*60}\nN = {N}  (stationary-subspace residue check)\n{'='*60}")
        d = 2**N
        J_A = [J_UNIFORM] * (N - 1)
        L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
        J_B = list(J_A); J_B[DEFECT_BOND[0]] = 1.1
        L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)

        P_stat_A, idx_A, ev_A = compute_stationary_projector(L_A)
        P_stat_B, idx_B, ev_B = compute_stationary_projector(L_B)

        print(f"  strict-stationary count: A = {len(idx_A)}, B = {len(idx_B)}"
              f"  (F4 predicts N+1 = {N+1})")

        # Is P_stat_A == P_stat_B? (J-invariance of stationary subspace)
        diff = np.linalg.norm(P_stat_A - P_stat_B)
        norm_ref = np.linalg.norm(P_stat_A)
        print(f"  || P_stat_A - P_stat_B ||_F = {diff:.3e}"
              f"  (relative {diff/norm_ref:.3e})")

        # Per-site Pauli contractions of P_stat: these capture the site-i
        # footprint of each strict stationary mode.
        print(f"  site-i Pauli | <P_Pauli|P_stat|P_Pauli> |:")
        for i in range(N):
            row_A, row_B = [], []
            for local_name, local in (("X", X), ("Y", Y), ("Z", Z)):
                v = site_pauli_matrix(i, local, N)
                a = float(abs(v.conj() @ P_stat_A @ v))
                b = float(abs(v.conj() @ P_stat_B @ v))
                row_A.append(a); row_B.append(b)
            print(f"    site {i}: A = " + "  ".join(f"{a:.4f}" for a in row_A) +
                  f"     B = " + "  ".join(f"{b:.4f}" for b in row_B))

        # Identity-basis overlap:
        # R_i = sum over Pauli_i of <Pauli_i | P_stat | Pauli_i>
        # If P_stat_A == P_stat_B, these are identical -> trivial closure.
        R_A = np.zeros(N); R_B = np.zeros(N)
        for i in range(N):
            for local in (X, Y, Z):
                v = site_pauli_matrix(i, local, N)
                R_A[i] += float(abs(v.conj() @ P_stat_A @ v))
                R_B[i] += float(abs(v.conj() @ P_stat_B @ v))
        log_ratio = np.log(R_B / R_A)
        print(f"  R_i = sum_pauli |<pauli|P_stat|pauli>|:")
        print(f"    A: " + "  ".join(f"{r:.4f}" for r in R_A))
        print(f"    B: " + "  ".join(f"{r:.4f}" for r in R_B))
        print(f"    Sum_i ln(R_B/R_A) = {np.sum(log_ratio):+.3e}")


if __name__ == "__main__":
    main()
