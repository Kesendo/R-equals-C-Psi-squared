#!/usr/bin/env python3
"""_eq022_b1_step_i_svd_inter_channel.py — heuristic 2-level basis via SVD
of inter-channel coupling.

Step (g): channel-uniform basis trivial V_b. Step (h): slowest-pair-at-Q=0.5
not the EP partners (eigenvalues still in HD=1 cluster, far from coalescence).

This step: define the heuristic 2-level subspace as the SVD of the inter-
channel block of M_H_total. The dominant singular vectors |u_0⟩ ∈ HD=1 and
|v_0⟩ ∈ HD=3 are the modes that maximally couple under H — these ARE the
EP partner pair. The largest singular value σ_0 = g_eff (the heuristic
2-level coupling).

Procedure (c=2 chains N=5..8):
  1. Build full block-L. M_H_total = sum_b M_H_per_bond[b].
  2. Build orthonormal projectors P_HD1 (M×n_HD1), P_HD3 (M×n_HD3) onto
     the two HD subspaces using popcount and bit-flip distance.
  3. Compute V_inter = P_HD1.conj().T @ M_H_total @ P_HD3  (n_HD1 × n_HD3).
  4. SVD: V_inter = U · diag(σ) · V^†. σ_0 = max singular value = g_eff.
  5. EP-partner modes: |u_0⟩ = P_HD1 · U[:, 0]  (in full block, length M).
                      |v_0⟩ = P_HD3 · V[:, 0]
  6. Heuristic 2-level basis: span{|u_0⟩, |v_0⟩}. P_2 = [|u_0⟩, |v_0⟩] (M × 2).
  7. Project M_H_total → 2×2 effective. Should be off-diagonal with phase ±i·σ_0
     (the heuristic same-sign-imaginary form).
  8. Project M_H_per_bond[b] → V_b for each bond. Bond-class averages.
  9. Project Dicke probe → 2-vector. Project S_kernel → 2×2.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def hd_subspace_projector(N, n, hd_value):
    """Orthonormal projector P_HD (M × n_HD) onto the HD-subspace.

    Each column is a basis vector (unit-norm coefficient vector in the flat
    (n, n+1)-block basis, supported only on (p, q) pairs with HD(p, q) = hd_value).
    """
    P_n_states = fw.popcount_states(N, n)
    P_np1_states = fw.popcount_states(N, n + 1)
    Mn = len(P_n_states)
    Mnp1 = len(P_np1_states)
    M = Mn * Mnp1
    p_to_idx = {p: i for i, p in enumerate(P_n_states)}
    q_to_idx = {q: i for i, q in enumerate(P_np1_states)}

    cols = []
    for p in P_n_states:
        for q in P_np1_states:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(M, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
    if not cols:
        return np.zeros((M, 0), dtype=complex)
    return np.column_stack(cols)


def fmt_c(z, prec=4):
    return f"{z.real:+.{prec}f}{z.imag:+.{prec}f}j"


def fmt_M(M, prec=4):
    rows = []
    for i in range(M.shape[0]):
        row = [fmt_c(M[i, j], prec) for j in range(M.shape[1])]
        rows.append("  [" + "  ".join(row) + "]")
    return "\n".join(rows)


def main():
    gamma_0 = 0.05
    cases = [(5, 1), (6, 1), (7, 1), (8, 1)]

    print("# EQ-022 (b1) Step (i): SVD of inter-channel coupling for c=2")
    print("# Goal: find heuristic 2-level basis (EP partner modes) via SVD")
    print()

    for (N, n) in cases:
        print("=" * 72)
        print(f"# c=2, N={N}, n={n}")
        print("=" * 72)

        D_full, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
        M_H_total = sum(M_H_per_bond)

        # HD subspace projectors (orthonormal columns)
        P_HD1 = hd_subspace_projector(N, n, 1)
        P_HD3 = hd_subspace_projector(N, n, 3)
        n_HD1 = P_HD1.shape[1]
        n_HD3 = P_HD3.shape[1]
        Mtot = D_full.shape[0]
        print(f"   Block dim = {Mtot}, HD=1 dim = {n_HD1}, HD=3 dim = {n_HD3}")

        # V_inter = M_H_total restricted between HD=1 and HD=3
        V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
        print(f"   V_inter shape = {V_inter.shape}")

        # SVD
        U, sigma, Vh = np.linalg.svd(V_inter)
        print(f"   Singular values (top 5): {[f'{s:.4f}' for s in sigma[:5]]}")
        print(f"   Number of nonzero (>1e-8) singular values: "
              f"{int((sigma > 1e-8).sum())} of min({n_HD1}, {n_HD3})")

        # Top singular vectors (the EP partner modes)
        u_0 = P_HD1 @ U[:, 0]   # in full block, supported on HD=1
        v_0 = P_HD3 @ Vh.conj().T[:, 0]  # in full block, supported on HD=3
        sigma_0 = sigma[0]
        print(f"   σ_0 = {sigma_0:.6f} (this is the heuristic g_eff)")
        print(f"   |u_0|² = {np.linalg.norm(u_0)**2:.6f} (should be 1)")
        print(f"   |v_0|² = {np.linalg.norm(v_0)**2:.6f}")

        # Build heuristic 2-level basis
        P_2 = np.column_stack([u_0, v_0])  # M × 2
        # M_H_total in this basis
        M_H_2 = P_2.conj().T @ M_H_total @ P_2
        print(f"   M_H_total in heuristic 2-level basis (should be off-diagonal):")
        print(fmt_M(M_H_2))

        # D_full in this basis (should be diagonal with -2γ₀ and -6γ₀)
        D_2 = P_2.conj().T @ D_full @ P_2
        print(f"   D_full in heuristic 2-level basis (should be diag(-{2*gamma_0}, -{6*gamma_0})):")
        print(fmt_M(D_2))

        # Per-bond V_b
        n_bonds = N - 1
        V_b_2 = [P_2.conj().T @ Mb @ P_2 for Mb in M_H_per_bond]

        print(f"   Per-bond V_b in heuristic 2-level basis:")
        for b, V_b in enumerate(V_b_2):
            label = "endpoint" if b in (0, n_bonds - 1) else "interior"
            print(f"   bond {b} ({label}):")
            print(fmt_M(V_b))

        # Bond-class averages
        endpoint_bonds = [0, n_bonds - 1]
        interior_bonds = list(range(1, n_bonds - 1))
        V_int = sum(V_b_2[b] for b in interior_bonds) / max(len(interior_bonds), 1)
        V_end = sum(V_b_2[b] for b in endpoint_bonds) / 2

        print(f"   ⟨V⟩_int (mean of {len(interior_bonds)}):")
        print(fmt_M(V_int))
        print(f"   ⟨V⟩_end (mean of 2):")
        print(fmt_M(V_end))

        # Decompose V_int and V_end into "kinematic" + "off-diagonal"
        # kinematic = diagonal part, coupling = off-diagonal
        print()
        print(f"   V_int diag entries: {fmt_c(V_int[0,0])}, {fmt_c(V_int[1,1])}")
        print(f"   V_int off-diag |V[0,1]| = {abs(V_int[0,1]):.5f}, "
              f"phase = {np.angle(V_int[0,1])/np.pi:.4f}π")
        print(f"   V_end diag entries: {fmt_c(V_end[0,0])}, {fmt_c(V_end[1,1])}")
        print(f"   V_end off-diag |V[0,1]| = {abs(V_end[0,1]):.5f}, "
              f"phase = {np.angle(V_end[0,1])/np.pi:.4f}π")

        # Probe in heuristic basis
        rho0 = fw.dicke_block_probe(N, n)
        probe_2 = P_2.conj().T @ rho0
        print(f"   probe_2 = ({fmt_c(probe_2[0])}, {fmt_c(probe_2[1])})")

        # S_kernel in heuristic basis
        S_full = fw.spatial_sum_coherence_kernel(N, n)
        S_2 = P_2.conj().T @ S_full @ P_2
        print(f"   S_kernel in heuristic 2-level basis:")
        print(fmt_M(S_2))

        # Direct test: build L_eff from this basis and check eigenvalue structure
        # L_eff(Q) = D_2 + Q·γ₀·sum(V_b) ≈ D_2 + Q·γ₀·M_H_2 (uniform J = Qγ₀)
        # At J = 2γ₀/σ_0, EP should occur with eigenvalue Re=-4γ₀
        Q_EP_pred = 2.0 / sigma_0
        J_EP = sigma_0 * gamma_0 * Q_EP_pred  # = 2γ₀ ✓
        L_eff_at_EP = D_2 + J_EP * M_H_2
        evals_at_EP = np.linalg.eigvals(L_eff_at_EP)
        print(f"   At Q_EP_pred = 2/σ_0 = {Q_EP_pred:.4f}, J = {J_EP:.5f}:")
        print(f"     L_eff eigenvalues = {[fmt_c(z) for z in evals_at_EP]}")
        print(f"     (EP would coalesce to -4γ₀ = {-4*gamma_0})")

        print()


if __name__ == "__main__":
    main()
