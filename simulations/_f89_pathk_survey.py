"""F89-(k) survey: multi-exponential closed-form derivation across path-k=2..5.

Bond-isolate verification + populated mode-group survey for path-k blocks at
N=7. See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Path-k survey" for
the result table and structural interpretation (path-3 privileged via
DE-self-symmetric Hamming complement at N_block=4).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from _f89_pathk_lib import (
    build_block_L,
    compute_rho_block_0,
    per_site_reduction_matrix,
    reduce_block_to_site_01,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"


def derive_path_k(n_block: int, J: float, gamma: float, N: int):
    """Build L_block, eigendecompose, project ρ_cc-derived ρ_block(0)."""
    D = 2**n_block
    L = build_block_L(J, gamma, n_block)
    rho = compute_rho_block_0(n_block, N)
    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)
    w = per_site_reduction_matrix(n_block)
    M = w @ R
    a = M * c[None, :]
    return L, eigvals, R, c, a, rho


def evolve(n_block: int, J: float, gamma: float, N: int, t_array: np.ndarray) -> np.ndarray:
    L, eigvals, R, c, a, rho = derive_path_k(n_block, J, gamma, N)
    D = 2**n_block
    x_bare_0 = (N - 1) / (2.0 * np.sqrt(N * N * (N - 1) / 2))
    S = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_block = sum(2.0 * abs(reduce_block_to_site_01(rho_t, l, n_block)) ** 2 for l in range(n_block))
        S_bare = (N - n_block) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare
    return S


def main() -> None:
    J, gamma = 0.075, 0.05
    N = 7
    print(f"# F89-(k) path-k closed-form survey, J={J}, γ={gamma}, N={N}\n")

    print("## Bond-isolate verification across path-k")
    print("| Path | bonds | max |diff| | mean |diff| | S(0) match |")
    print("|---|---|---|---|---|")
    for n_block in [3, 4, 5]:
        bonds_str = "-".join(str(b) for b in range(n_block - 1))
        csv = CSV_DIR / f"N{N}_b{bonds_str}_J0.0750_gamma0.0500_probe-coherence.csv"
        if not csv.exists():
            print(f"| path-{n_block-1} | {bonds_str} | (CSV missing) | | |")
            continue
        data = np.loadtxt(csv, delimiter=",", skiprows=1)
        t_csv, S_csv = data[:, 0], data[:, -1]
        S_pred = evolve(n_block, J, gamma, N, t_csv)
        diff = S_pred - S_csv
        match = "✓" if abs(S_pred[0] - (N - 1) / N) < 1e-6 else "✗"
        print(
            f"| path-{n_block-1} | {{{bonds_str}}} | {np.max(np.abs(diff)):.2e} | "
            f"{np.mean(np.abs(diff)):.2e} | {match} |"
        )
    print()

    print("## Populated mode-group survey")
    print("| Path | N_block | d² | Mode-groups | Contributing modes | Pair-sums to 2γ·N_block |")
    print("|---|---|---|---|---|---|")
    for n_block in [3, 4, 5, 6]:
        L, eigvals, R, c, a, rho = derive_path_k(n_block, J, gamma, N)
        sig = np.sum(np.abs(a) ** 2, axis=0)
        contributing = np.where(sig > 1e-12)[0]
        rates = -eigvals.real / gamma
        freqs = np.abs(eigvals.imag) / J
        groups = {}
        for k in contributing:
            key = (round(rates[k], 4), round(freqs[k], 4))
            groups[key] = groups.get(key, 0.0) + sig[k]
        unique_rates = sorted({k[0] for k in groups.keys()})
        target_pair_sum = 2 * n_block
        n_pairs = sum(
            1 for r1 in unique_rates for r2 in unique_rates
            if abs((r1 + r2) - target_pair_sum) < 0.001
        )
        print(
            f"| path-{n_block-1} | {n_block} | {4**n_block} | "
            f"{len(groups)} | {len(contributing)} | "
            f"{n_pairs}{' ✓' if n_pairs > 0 else ''} |"
        )


if __name__ == "__main__":
    main()
