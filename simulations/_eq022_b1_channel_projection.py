#!/usr/bin/env python3
"""_eq022_b1_channel_projection.py — EQ-022 (b1) analytical attempt.

Goal: derive Q_peak(c) closed-form from the (n, n+1) block-L structure.

Strategy attempted
------------------
Build the popcount-(n, n+1) coherence block of L = L_D + J·L_H for the
uniform XY chain under Z-dephasing. Project onto the c-dimensional
HD-channel basis (orthonormalized uniform-superposition projectors over
states with HD = 2k+1 for k = 0..c-1). The hope: if the resulting c×c
effective L_eff is N-invariant at fixed c, then Q_peak is a function of
the c×c matrix structure alone and a closed form is reachable.

Result (2026-05-02)
-------------------
The c×c effective M_H_eff (the J-coefficient matrix) turns out to be
**purely diagonal** in the channel-uniform basis: all off-diagonal entries
are exactly zero. The diagonal entries are i·g_k(N, n) with g_k real and
N-dependent (g_k = 0 at the minimal-chromaticity case N = 2c−1, growing
with N).

**Consequence:** in this projection, L_eff(J) = diag(-2γ₀(2k−1) + i·J·g_k),
whose Re(λ) is J-independent. Q_peak as defined via the J-derivative of
the slowest decay rate vanishes identically in this basis. The cross-
channel mixing that drives Q_peak lives in the *orthogonal complement* of
the channel-uniform subspace, NOT in the c-dim channel subspace itself.

**Structural finding worth keeping:** the channel-uniform projectors |c_k⟩
are H-invariant within their own HD-channel (modulo orthogonal complement).
This generalises F73 (which is the c=1 case at the (0, 1) block) to every
chromaticity. The cancellation of off-diagonals between HD-channels is a
sign-symmetry property of single-bond H matrix elements summed over all
HD=k pairs.

**What this leaves open for (b1):** the closed-form Q_peak(c) requires a
richer effective model that includes coupling to the channel-uniform
complement. Candidate next approach: extend the projection to the
two-step-larger sector (basis = channel-uniforms + non-uniform within-
channel modes that H couples them to), giving a > c-dim effective model.
Alternatively: scan full block-L spectrum at growing N, identify which
eigenvectors carry the Dicke-probe weight near Q_peak, and read the
effective coupling from those.

Reference: F74 chromaticity, F86 Q_peak constants (added 2026-05-02),
Q_SCALE_THREE_BANDS, EQ-022 (b1) in review/EMERGING_QUESTIONS.md.

WIP per `_` prefix convention. Exploration tool, not a primitive.

Usage:
    python _eq022_b1_channel_projection.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def build_L_cc_block_split(N, n, gamma_0):
    """Compatibility shim: returns (D, M_H_total, P_n, P_np1) using the
    framework primitive `block_L_split_xy`. Step-b/c scripts import this name.
    """
    D, M_H_per_bond, P_n, P_np1 = fw.block_L_split_xy(N, n, gamma_0)
    return D, sum(M_H_per_bond), P_n, P_np1


def channel_effective_split(N, n, gamma_0):
    """Returns (D_eff, M_H_eff, HDs) such that L_eff(J) = D_eff + J * M_H_eff."""
    D, M_H, _, _ = build_L_cc_block_split(N, n, gamma_0)
    P, HDs = fw.hd_channel_basis(N, n)
    D_eff = P.conj().T @ D @ P
    M_H_eff = P.conj().T @ M_H @ P
    return D_eff, M_H_eff, HDs


# ---------------------------------------------------------------------------
# Q_peak from c×c effective: scan J on the small matrix
# ---------------------------------------------------------------------------
def find_Q_peak_from_split(D_eff, M_H_eff, gamma_0, J_grid):
    """Scan J on c×c matrix L_eff(J) = D_eff + J·M_H_eff. Track slowest
    Re(λ); peak |dRe(λ_slow)/dQ| approximates Q_peak."""
    Qs = J_grid / gamma_0
    slowest_re = np.zeros_like(J_grid)
    for i, J in enumerate(J_grid):
        L = D_eff + J * M_H_eff
        evals = np.linalg.eigvals(L)
        slowest_re[i] = np.max(np.real(evals))
    dQ = Qs[1] - Qs[0]
    deriv = np.gradient(slowest_re, dQ)
    abs_deriv = np.abs(deriv)
    peak_idx = np.argmax(abs_deriv)
    return Qs[peak_idx], abs_deriv[peak_idx], slowest_re


# ---------------------------------------------------------------------------
# Main exploration
# ---------------------------------------------------------------------------
def main():
    gamma_0 = 0.05
    print(f"# EQ-022 (b1): channel projection scan, gamma_0 = {gamma_0}")
    print()

    # --- N-invariance check at fixed c ---
    print("## c×c effective L_eff entries at J = gamma_0 (Q = 1) per (N, n)")
    print()

    test_cases = [
        # c = 2
        (3, 1), (4, 1), (5, 1), (6, 1), (7, 1),
        # c = 3
        (5, 2), (6, 2), (7, 2),
        # c = 4
        (7, 3), (8, 3), (9, 3),
        # c = 5
        (9, 4),
    ]

    by_c = {}
    for (N, n) in test_cases:
        c = fw.chromaticity(N, n)
        D_eff, M_H_eff, HDs = channel_effective_split(N, n, gamma_0)
        by_c.setdefault(c, []).append((N, n, D_eff, M_H_eff, HDs))

    for c, runs in sorted(by_c.items()):
        print(f"### c = {c}")
        for (N, n, D_eff, M_H_eff, HDs) in runs:
            print(f"  (N={N}, n={n}, HDs={HDs}):")
            with np.printoptions(precision=4, suppress=True, linewidth=120):
                print("  D_eff / gamma_0 =")
                print(D_eff / gamma_0)
                print("  M_H_eff (J-coefficient, units 1/J) =")
                print(M_H_eff)
            print()

    # --- Q_peak from channel-effective L ---
    print("\n## Q_peak from channel-effective L (fast J-grid scan)")
    print()
    J_grid = np.linspace(0.0, 0.5, 501)  # Q in [0, 10], dQ = 0.02

    for c, runs in sorted(by_c.items()):
        for (N, n, D_eff, M_H_eff, _) in runs:
            Q_peak, peak_val, _ = find_Q_peak_from_split(D_eff, M_H_eff, gamma_0, J_grid)
            print(f"  c={c}, N={N}, n={n}: Q_peak (effective) = {Q_peak:.3f}, peak |dRe(λ_slow)/dQ| = {peak_val:.3f}")

    # --- Direct comparison to F86 constants ---
    print()
    print("## F86 expected: Q_peak(c=3) = 1.6, Q_peak(c=4) = Q_peak(c=5) = 1.8")
    print("## c=2 wobbles 1.4-1.6 (finite-size sensitive)")

    # --- Sanity: full block-L Q_peak via Dicke probe (no projection) ---
    # The probe S(t) = Σᵢ 2 |(ρᵢ)_{0,1}|² is the F73 spatial-sum coherence.
    # K_CC_pr = ∂S/∂J. We compute Q where |∂S/∂J|_max occurs from full block-L
    # to verify the script reproduces Q_SCALE_THREE_BANDS values.
    print()
    print("## Sanity check: full block-L Q_peak via Dicke probe S(t) at fixed t")
    print("##   (S = ⟨ρ_cc(t), full-block-L eigenvectors⟩ traced into spatial-sum)")
    print()

    from scipy.linalg import expm
    for (N, n) in [(4, 1), (5, 2), (6, 2), (7, 3)]:
        c = fw.chromaticity(N, n)
        D, M_H, _, _ = build_L_cc_block_split(N, n, gamma_0)
        rho0 = fw.dicke_block_probe(N, n)
        S_kernel = fw.spatial_sum_coherence_kernel(N, n)
        t_obs = 1.0  # arbitrary fixed observation time

        Qs = np.linspace(0.05, 5.0, 100)
        J_vals = Qs * gamma_0
        S_vals = np.zeros_like(J_vals)
        for i, J in enumerate(J_vals):
            L = D + J * M_H
            rho_t = expm(L * t_obs) @ rho0
            S_vals[i] = float(np.real(rho_t.conj() @ S_kernel @ rho_t))
        dS_dJ = np.gradient(S_vals, J_vals)
        peak_idx = int(np.argmax(np.abs(dS_dJ)))
        Q_peak = Qs[peak_idx]
        print(f"  full-block S(t={t_obs}): N={N}, n={n}, c={c}: Q_peak ≈ {Q_peak:.3f}, |dS/dJ|_peak = {np.abs(dS_dJ[peak_idx]):.4f}")

    print()
    print("## Note: t_obs = 1.0 is arbitrary; Q_peak in F86 is defined via")
    print("## abs(K_CC_pr)_max which is the J-DERIVATIVE OF S optimised over t.")
    print("## A proper reproduction would scan (J, t) jointly. This is sanity-check only.")


if __name__ == "__main__":
    main()
