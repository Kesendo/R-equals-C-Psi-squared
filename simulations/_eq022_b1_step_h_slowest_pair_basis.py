#!/usr/bin/env python3
"""_eq022_b1_step_h_slowest_pair_basis.py — find the basis where the heuristic
same-sign-imaginary 2-level form actually lives.

Step (g) showed: in the channel-uniform basis, V_b for c=2 chains are all
+iα/(N-1)·I (diagonal, identical per bond). No EP, no bond-class distinction.

Conclusion: the heuristic 2-level form is NOT in channel-uniform basis. The
EP must live in the orthogonal complement, and the bond-class distinction
comes from the slowest-pair eigenvectors of the FULL block-L (which mix
channel-uniform with its complement).

Procedure (this step):
  1. Build full block-L at one Q value just below the EP.
  2. Diagonalize → identify the SLOWEST PAIR of eigenvalues (closest to
     Re(λ) = −4γ₀ at the EP).
  3. Get right and left biorthogonal eigenvectors (|v_+⟩, |v_-⟩, ⟨w_+|, ⟨w_-|).
  4. Project full block-L onto this 2D subspace → 2×2 effective matrix.
     Check: does it have the heuristic same-sign-imaginary +iJg_eff form?
  5. Project M_H_per_bond[b] onto the 2D subspace for each bond b.
     Examine: do bond-class averages now differ between Endpoint and Interior?
  6. Project Dicke probe and S_kernel onto the 2D subspace.
  7. Run analytical 2-level prediction with the slowest-pair-basis V_b.
     Compare to full-block-L K_b(Q) curves.

Run on c=2 N=5..8.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

RESULTS_DIR = REPO_ROOT / "simulations" / "results" / "eq022_slowest_pair_basis"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def slowest_pair_basis_at_Q(N, n, gamma_0, Q):
    """Diagonalize full block-L at given Q and return slowest pair structure.

    Slowest pair = the two eigenvalues with largest Re(λ) (smallest decay rate)
    that emerge from the HD=1 / HD=3 channels and approach each other near
    the EP.

    Returns dict with:
      L_full: 50x50 (or larger) full block-L matrix
      M_H_per_bond: list of N-1 per-bond H coefficient matrices
      D_full: J=0 dephasing diagonal
      evals: all eigenvalues sorted by Re(λ) descending
      idx_pair: indices of the two slowest (largest Re) modes
      lam_+, lam_-: the two eigenvalues
      v_+, v_-: right eigenvectors (column vectors)
      w_+, w_-: left eigenvectors (row vectors via biorthogonal of M's right
                eigenvalues)
      P_R: 2×2 column matrix of right eigenvectors
      P_L: 2×2 row matrix of left eigenvectors (biorthonormalized)
      L_eff_2x2: P_L · L_full · P_R (should be diagonal with lam_±)
      V_b_2x2: list of 2×2 V_b in slowest-pair basis
      probe_2x2: probe in slowest-pair basis
      S_kernel_2x2: S_kernel projected onto slowest-pair basis
    """
    D_full, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    J = Q * gamma_0
    L_full = D_full + J * sum(M_H_per_bond)

    evals, R = np.linalg.eig(L_full)
    # Sort by Re(λ) descending (slowest = largest Re, since Re < 0)
    order = np.argsort(-evals.real)
    evals = evals[order]
    R = R[:, order]

    # Slowest pair = indices 0, 1
    idx_pair = (0, 1)
    lam_p = evals[0]
    lam_m = evals[1]
    v_p = R[:, 0]
    v_m = R[:, 1]

    # Left eigenvectors: rows of L_full = inverse(R)
    R_inv = np.linalg.inv(R)
    w_p = R_inv[0, :]  # row vector
    w_m = R_inv[1, :]

    # Biorthonormality check: ⟨w_+|v_+⟩ should be 1, ⟨w_+|v_-⟩ should be 0
    bio_pp = w_p @ v_p
    bio_pm = w_p @ v_m
    bio_mp = w_m @ v_p
    bio_mm = w_m @ v_m

    # Build projection matrices
    P_R = np.column_stack([v_p, v_m])  # right basis (columns)
    P_L = np.vstack([w_p, w_m])         # left basis (rows)

    # 2×2 effective matrices in slowest-pair basis
    L_eff_2x2 = P_L @ L_full @ P_R
    V_b_2x2 = [P_L @ Mb @ P_R for Mb in M_H_per_bond]

    # Probe in slowest-pair basis
    rho0 = fw.dicke_block_probe(N, n)
    probe_2x2 = P_L @ rho0  # 2-vector

    # S_kernel in slowest-pair basis (S is Hermitian PSD on full)
    S_full = fw.spatial_sum_coherence_kernel(N, n)
    # K(t) = 2·Re(⟨ρ(t)| S | ∂ρ/∂J⟩). In slowest-pair basis:
    # ρ(t) = P_R · ρ_2x2(t), so ⟨ρ(t)| S |...⟩ = ρ_2x2(t)^† · (P_R^† · S · P_R) · ...
    # But ρ_2x2 is the LEFT-projected state. To get the standard inner product back,
    # we need to project S consistent with the basis transformations.
    # The cleanest: in slowest-pair basis with right vectors P_R, define
    #   ρ(t) = P_R · a(t)  where a(t) = P_L · ρ_full(t)  (left projection of full state)
    # Then ⟨ρ_full(t)| S |ψ⟩ = a(t)^† · P_R^† · S · ψ (full)
    # If ψ is also expressed via P_R: ψ = P_R · b, so ⟨ρ_full| S |ψ⟩ = a^† · (P_R^† · S · P_R) · b
    # Hence S_kernel_2x2 = P_R^† · S_full · P_R (NOT P_L · S · P_R, since S acts on full ρ).
    S_kernel_2x2 = P_R.conj().T @ S_full @ P_R

    return dict(
        L_full=L_full, M_H_per_bond=M_H_per_bond, D_full=D_full,
        evals=evals, lam_p=lam_p, lam_m=lam_m,
        v_p=v_p, v_m=v_m, w_p=w_p, w_m=w_m,
        bio_pp=bio_pp, bio_pm=bio_pm, bio_mp=bio_mp, bio_mm=bio_mm,
        P_R=P_R, P_L=P_L,
        L_eff_2x2=L_eff_2x2, V_b_2x2=V_b_2x2,
        probe_2x2=probe_2x2, S_kernel_2x2=S_kernel_2x2,
    )


def fmt_complex(z, prec=4):
    return f"{z.real:+.{prec}f}{z.imag:+.{prec}f}j"


def fmt_matrix(M, prec=4):
    rows = []
    for i in range(M.shape[0]):
        row = [fmt_complex(M[i, j], prec) for j in range(M.shape[1])]
        rows.append("  [" + "  ".join(row) + "]")
    return "\n".join(rows)


def main():
    gamma_0 = 0.05
    cases = [(5, 1), (6, 1), (7, 1), (8, 1)]
    print("# EQ-022 (b1) Step (h): slowest-pair basis exploration for c=2")
    print()

    # Use Q = 0.5·Q_EP_estimate as a sample point well below EP.
    # Q_EP is roughly 1.5-2.5 for c=2 chains (close to Q_peak), so Q_sample = 0.5
    # is well below.
    Q_sample = 0.5

    for (N, n) in cases:
        print("=" * 72)
        print(f"# c=2, N={N}, Q_sample = {Q_sample} (below EP for analysis)")
        print("=" * 72)

        d = slowest_pair_basis_at_Q(N, n, gamma_0, Q_sample)

        # Show all eigenvalues (small enough at c=2 to print)
        print(f"   Eigenvalues of L_full (sorted by Re descending, top 6):")
        for k in range(min(6, len(d['evals']))):
            print(f"     λ_{k} = {fmt_complex(d['evals'][k])}")
        print(f"   Slowest pair: λ_+ = {fmt_complex(d['lam_p'])}, "
              f"λ_- = {fmt_complex(d['lam_m'])}")
        print(f"   Re separation |Re(λ_+) - Re(λ_-)| = "
              f"{abs(d['lam_p'].real - d['lam_m'].real):.5f}")
        print(f"   Im separation |Im(λ_+) - Im(λ_-)| = "
              f"{abs(d['lam_p'].imag - d['lam_m'].imag):.5f}")
        print()
        print(f"   Biorthogonality check:")
        print(f"     ⟨w_+|v_+⟩ = {fmt_complex(d['bio_pp'])}, ⟨w_+|v_-⟩ = {fmt_complex(d['bio_pm'])}")
        print(f"     ⟨w_-|v_+⟩ = {fmt_complex(d['bio_mp'])}, ⟨w_-|v_-⟩ = {fmt_complex(d['bio_mm'])}")
        print()
        print(f"   L_eff_2x2 (should be diag with lam_+, lam_-):")
        print(fmt_matrix(d['L_eff_2x2']))
        print()
        print(f"   probe_2x2: ({fmt_complex(d['probe_2x2'][0])}, {fmt_complex(d['probe_2x2'][1])})")
        print(f"   S_kernel_2x2:")
        print(fmt_matrix(d['S_kernel_2x2']))
        print()
        print(f"   Per-bond V_b in slowest-pair basis (sample):")
        n_bonds = N - 1
        for b in range(n_bonds):
            label = "endpoint" if b in (0, n_bonds - 1) else "interior"
            print(f"   bond {b} ({label}):")
            print(fmt_matrix(d['V_b_2x2'][b]))
        print()
        # Bond-class averages
        endpoint_bonds = [0, n_bonds - 1]
        interior_bonds = list(range(1, n_bonds - 1))
        V_int = sum(d['V_b_2x2'][b] for b in interior_bonds) / max(len(interior_bonds), 1)
        V_end = sum(d['V_b_2x2'][b] for b in endpoint_bonds) / 2
        print(f"   ⟨V⟩_int (mean of {len(interior_bonds)} interior bonds):")
        print(fmt_matrix(V_int))
        print(f"   ⟨V⟩_end (mean of 2 endpoint bonds):")
        print(fmt_matrix(V_end))
        # Differences (interior - endpoint)
        delta = V_int - V_end
        print(f"   Δ = ⟨V⟩_int − ⟨V⟩_end:")
        print(fmt_matrix(delta))
        print(f"   |Δ_00|={abs(delta[0,0]):.4e}, |Δ_01|={abs(delta[0,1]):.4e}")
        print(f"   |V_int_00|={abs(V_int[0,0]):.4e}, |V_int_01|={abs(V_int[0,1]):.4e}")
        print(f"   relative |Δ/V_int|: diag={abs(delta[0,0])/max(abs(V_int[0,0]), 1e-12):.3f}, "
              f"offdiag={abs(delta[0,1])/max(abs(V_int[0,1]), 1e-12):.3f}")
        print()


if __name__ == "__main__":
    main()
