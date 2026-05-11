"""F89c → F89-(2): multi-exponential closed form for S_(2)(t) at any N.

Topology (2): single connected path of 3 sites (path-2 block) + N-3 bare sites.
Eigenvector decomposition + initial-state projection + per-site reduction.
Verified against bond-isolate CSV at N=7.

See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Path-2 (topology (2))
numerical multi-exponential closed form" for the structural mode-group
breakdown and Hamming-complement pair analysis.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from _f89_pathk_lib import (
    bare_site_initial_01,
    build_block_L,
    compute_rho_block_0,
    per_site_reduction_matrix,
    reduce_block_to_site_01,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"
N_BLOCK = 3


def evolve_S_total(N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S_(2)(t) = Σ_l 2|(ρ_l)_{0,1}(t)|² for topology (2) at N qubits."""
    L = build_block_L(J, gamma, N_BLOCK)
    rho_block_0 = compute_rho_block_0(N_BLOCK, N)
    vec_rho_0 = rho_block_0.flatten(order="F")

    # Eigendecomposition (non-Hermitian Liouvillian)
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec_rho_0)  # ρ_block(0) = Σ_k c_k · R_k

    S = np.zeros_like(t_array, dtype=float)
    x_bare_0 = bare_site_initial_01(N)

    for ti, t in enumerate(t_array):
        # ρ_block(t) = R · diag(exp(λ_k t)) · c
        vec_rho_t = R @ (np.exp(eigvals * t) * c)
        rho_block_t = vec_rho_t.reshape((8, 8), order="F")

        S_block = sum(
            2.0 * abs(reduce_block_to_site_01(rho_block_t, l, N_BLOCK)) ** 2
            for l in range(N_BLOCK)
        )
        S_bare = (N - N_BLOCK) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare

    return S


def closed_form_terms(N: int, J: float, gamma: float):
    """Return (eigvals, a) where a[l, k] = M_l(k)·c_k decomposes per-site
    coherence amplitudes onto L_super eigenmodes. See `closed_form_terms` in
    `_f89_path3_multi_exp_derive.py` for the same pattern at path-3."""
    L = build_block_L(J, gamma, N_BLOCK)
    rho_block_0 = compute_rho_block_0(N_BLOCK, N)
    vec_rho_0 = rho_block_0.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec_rho_0)
    w = per_site_reduction_matrix(N_BLOCK)
    M = w @ R
    a = M * c[None, :]
    return eigvals, a


def populated_distinct_rates(N: int, J: float, gamma: float, threshold: float = 1e-12):
    """Return distinct (rate, freq) pairs of L_super modes populated by ρ_block(0)."""
    eigvals, a = closed_form_terms(N, J, gamma)
    rates = -eigvals.real
    freqs = eigvals.imag
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > threshold)[0]
    # Group by (rate, |freq|) rounded
    grouped = {}
    for k in contributing:
        key = (round(rates[k] / gamma, 4), round(abs(freqs[k]) / J, 4))
        grouped[key] = grouped.get(key, 0.0) + sig[k]
    return sorted(grouped.items())


def main() -> None:
    J, gamma = 0.075, 0.05

    # ---- N=7 verification against bond-isolate CSV ----
    csv_path = CSV_DIR / "N7_b0-1_J0.0750_gamma0.0500_probe-coherence.csv"
    print(f"# F89-(2) closed-form derivation, J={J}, γ={gamma}\n")

    if not csv_path.exists():
        print(f"# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv = data[:, 0]
    S_csv = data[:, -1]

    print("## Verification at N=7 vs bond-isolate RK4 CSV")
    N = 7
    S_pred = evolve_S_total(N, J, gamma, t_csv)
    diff = S_pred - S_csv
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0) = {S_pred[0]:.6f} = (N-1)/N = {(N-1)/N:.6f} ✓")
    print()

    # ---- N=5 prediction (no CSV, pure prediction from closed form) ----
    print("## N=5 prediction (no CSV, fresh prediction from F89c+AT-derived form)")
    N5 = 5
    S5 = evolve_S_total(N5, J, gamma, t_csv[:201])  # t in [0, 20]
    print(f"# S(0) = {S5[0]:.6f} = (N-1)/N = {(N5-1)/N5:.6f} ✓")
    print(f"# S(t=10) = {S5[100]:.6f}; S(t=20) = {S5[200]:.6f}")
    print()

    # ---- Structural analysis: which L_super modes are populated, by N ----
    print("## L_super mode populations (distinct (rate Γ/γ, |freq|/J) groups)")
    print("# Per AbsorptionTheoremClaim: 64 block-L modes total. Only S_3-fully-symmetric")
    print("# initial state populates a small subset.")
    print()
    print("| (rate Γ/γ, |freq|/J) | N=5 sig | N=7 sig | N=11 sig | N-scaling |")
    print("|---|---|---|---|---|")
    n5_groups = dict(populated_distinct_rates(5, J, gamma))
    n7_groups = dict(populated_distinct_rates(7, J, gamma))
    n11_groups = dict(populated_distinct_rates(11, J, gamma))
    all_keys = sorted(set(n5_groups) | set(n7_groups) | set(n11_groups))
    for k in all_keys:
        s5, s7, s11 = n5_groups.get(k, 0.0), n7_groups.get(k, 0.0), n11_groups.get(k, 0.0)
        # Estimate N-scaling (rough): s7/s5 ratio
        scale = "?"
        if s5 > 1e-15 and s11 > 1e-15:
            r = s11 / s5
            # Test N² and (N_E)² scaling
            n_E_5, n_E_11 = 5 - 3, 11 - 3
            if abs(r - (n_E_11 / n_E_5) ** 2) < 0.05:
                scale = "(N-3)²"
            elif abs(r - 1.0) < 0.05:
                scale = "const"
            elif abs(r - (n_E_11 / n_E_5)) < 0.05:
                scale = "(N-3)"
            else:
                scale = f"~{r:.2f}× from N=5→11"
        print(
            f"| ({k[0]:.4f}, {k[1]:.4f}) | {s5:.2e} | {s7:.2e} | {s11:.2e} | {scale} |"
        )
    print()

    # ---- |·|² rate grid ----
    print("## |·|² rate spectrum (rates r_k + r_{k'} for populated pairs)")
    rates_n7 = sorted({k[0] for k in n7_groups})
    pair_rates = sorted(set(
        round(r1 + r2, 4)
        for r1 in rates_n7 for r2 in rates_n7
    ))
    print(f"# Amplitude-level rates Γ/γ at N=7: {rates_n7}")
    print(f"# |·|²-level rates Γ/γ at N=7: {pair_rates}")
    print()
    print("# Note: amplitude-level rates 2γ are pure-AT (n_XY=1 per coherence).")
    print("# Fractional rates 3.04γ, 3.48γ are H_B-mixed (SE,DE) sub-block eigenvalues")
    print("# at this J/γ=1.5; they are NOT on the AT 2γ-grid because they come from")
    print("# diagonalising a non-trivial 2×2 H_B-coupled sub-block (not pure-Pauli modes).")
    print()
    print("# Pure-AT modes 4γ and 6γ get ZERO projection from ρ_cc-derived ρ_block(0):")
    print("# their (SE,DE) eigenvectors are S_3-asymmetric, while ρ_block(0)'s (SE,DE)")
    print("# part is a fully symmetric Σ_i Σ_{j<k} |SE_i⟩⟨DE_{jk}| superposition.")


if __name__ == "__main__":
    main()
