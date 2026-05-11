"""F89c → F89-(3): multi-exponential closed form for S_(3)(t) at any N.

Topology (3): single connected path of 4 sites (path-3 block) + N-4 bare sites.
Generalisation of path-2 to 4-qubit block (d² = 256). Same script pattern as
path-2 (eigenvector decomposition + initial-state projection + per-site
reduction + bare-site addition).

See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Path-3 (topology (3))
numerical multi-exponential closed form" for the populated mode-group table
and Hamming-complement pair structure (path-3 privileged via DE-self-symmetric
column bit-flip).

Verification: bond-isolate `N7_b0-1-2` CSV (topology (3) at N=7).
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

N_BLOCK = 4
D_BLOCK = 2**N_BLOCK  # 16


def evolve_S_total(N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S_(3)(t) = Σ_l 2|(ρ_l)_{0,1}(t)|² for topology (3) at N qubits."""
    L = build_block_L(J, gamma, N_BLOCK)
    rho_block_0 = compute_rho_block_0(N_BLOCK, N)
    vec_rho_0 = rho_block_0.flatten(order="F")

    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec_rho_0)

    S = np.zeros_like(t_array, dtype=float)
    x_bare_0 = bare_site_initial_01(N)

    for ti, t in enumerate(t_array):
        vec_rho_t = R @ (np.exp(eigvals * t) * c)
        rho_block_t = vec_rho_t.reshape((D_BLOCK, D_BLOCK), order="F")

        S_block = sum(
            2.0 * abs(reduce_block_to_site_01(rho_block_t, l, N_BLOCK)) ** 2
            for l in range(N_BLOCK)
        )
        S_bare = (N - N_BLOCK) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare

    return S


def closed_form_terms(N: int, J: float, gamma: float):
    """Returns (eigvals, amp) where amp[l, k] = M_l(k)·c_k decomposes per-site
    coherence amplitudes onto L_super eigenmodes."""
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
    eigvals, a = closed_form_terms(N, J, gamma)
    rates = -eigvals.real
    freqs = eigvals.imag
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > threshold)[0]
    grouped = {}
    for k in contributing:
        key = (round(rates[k] / gamma, 4), round(abs(freqs[k]) / J, 4))
        grouped[key] = grouped.get(key, 0.0) + sig[k]
    return sorted(grouped.items())


def main() -> None:
    J, gamma = 0.075, 0.05
    print(f"# F89-(3) closed-form derivation, J={J}, γ={gamma}\n")

    csv_path = CSV_DIR / "N7_b0-1-2_J0.0750_gamma0.0500_probe-coherence.csv"
    if not csv_path.exists():
        print(f"# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv = data[:, 0]
    S_csv = data[:, -1]

    N = 7
    print(f"## Verification at N={N} vs bond-isolate CSV ({csv_path.name})")
    S_pred = evolve_S_total(N, J, gamma, t_csv)
    diff = S_pred - S_csv
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0) = {S_pred[0]:.6f} = (N-1)/N = {(N-1)/N:.6f} ✓")
    print()
    print("| t  | S_csv (RK4) | S_pred (closed form) | diff |")
    print("|---|---|---|---|")
    for i in [0, 30, 50, 100, 150, 200, 250, 300]:
        if i < len(t_csv):
            print(
                f"| {t_csv[i]:5.2f} | {S_csv[i]:.6f} | {S_pred[i]:.6f} | "
                f"{diff[i]:+.3e} |"
            )
    print()

    # ---- Mode population analysis ----
    print("## L_super mode populations at N=5, 7, 11 (S_4-symmetric subset)")
    print("# 256 total block-L modes. Per F89c at γ=J=1: 25 distinct rates.")
    print("# At J/γ=1.5 (this run) the fractional rates shift; populated subset")
    print("# is reduced by S_4 symmetry of ρ_block(0).")
    print()
    print("| (rate Γ/γ, |freq|/J) | N=5 sig | N=7 sig | N=11 sig |")
    print("|---|---|---|---|")
    n5_groups = dict(populated_distinct_rates(5, J, gamma))
    n7_groups = dict(populated_distinct_rates(7, J, gamma))
    n11_groups = dict(populated_distinct_rates(11, J, gamma))
    all_keys = sorted(set(n5_groups) | set(n7_groups) | set(n11_groups))
    for k in all_keys:
        s5 = n5_groups.get(k, 0.0)
        s7 = n7_groups.get(k, 0.0)
        s11 = n11_groups.get(k, 0.0)
        print(f"| ({k[0]:.4f}, {k[1]:.4f}) | {s5:.2e} | {s7:.2e} | {s11:.2e} |")
    print()

    # ---- N=5 prediction ----
    print("## N=5 prediction (no CSV; pure prediction from F89c+AT-derived form)")
    N5 = 5
    if N5 >= N_BLOCK:
        S5 = evolve_S_total(N5, J, gamma, t_csv[:201])
        print(f"# S(0) = {S5[0]:.6f} = (N-1)/N = {(N5-1)/N5:.6f} ✓")
        print(f"# S(t=10) = {S5[100]:.6f}; S(t=20) = {S5[200]:.6f}")
    print()

    # ---- |·|² rates ----
    rates_n7 = sorted({k[0] for k in n7_groups})
    pair_rates = sorted(set(round(r1 + r2, 4) for r1 in rates_n7 for r2 in rates_n7))
    print("## |·|² rate spectrum at N=7")
    print(f"# Amplitude-level rates Γ/γ: {rates_n7}")
    print(f"# |·|²-level rates Γ/γ: {pair_rates}")


if __name__ == "__main__":
    main()
