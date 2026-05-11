"""F89 path-6 (single 7-qubit chain at N=7, no bare sites): numerical multi-exp.

Topology (6) at N=7: full-chain block. L_super dim = 4^7 = 16384.
np.linalg.eig at this dim: ~12 GB RAM, ~30-60 min (often longer).

Trivially satisfies the additive identity since m=1 → no subtraction term:
S_(6)(t) = block contribution = full result. See
`experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` for context.
"""

from __future__ import annotations

import sys
import time
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


def main() -> None:
    J, gamma = 0.075, 0.05
    N = 7
    n_block = 7
    D = 2**n_block

    print(f"# F89 path-6 (full chain at N={N}): J={J}, γ={gamma}\n")
    print(f"# L_super dim: {D*D} ({D}×{D} block)")

    t1 = time.time()
    print("# Building L_super (16384×16384)...", flush=True)
    L = build_block_L(J, gamma, n_block)
    print(f"# L_super built in {time.time() - t1:.1f}s; size {L.nbytes / 1024**3:.2f} GB", flush=True)

    print("# Computing ρ_block(0)...", flush=True)
    rho = compute_rho_block_0(n_block, N)
    vec = rho.flatten(order="F")

    t2 = time.time()
    print("# Eigendecomposing L_super (slow: ~30-60 min)...", flush=True)
    eigvals, R = np.linalg.eig(L)
    print(f"# Eigendecomposition done in {time.time() - t2:.1f}s", flush=True)

    t3 = time.time()
    print("# Solving R c = vec(ρ_0)...", flush=True)
    c = np.linalg.solve(R, vec)
    print(f"# Solve done in {time.time() - t3:.1f}s", flush=True)

    csv_path = CSV_DIR / f"N{N}_b0-1-2-3-4-5_J0.0750_gamma0.0500_probe-coherence.csv"
    if not csv_path.exists():
        print(f"\n# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv, S_csv = data[:, 0], data[:, -1]

    t4 = time.time()
    print("# Computing S(t) at all CSV time points...", flush=True)
    S_pred = np.zeros_like(t_csv)
    for ti, t in enumerate(t_csv):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_pred[ti] = sum(2.0 * abs(reduce_block_to_site_01(rho_t, l, n_block)) ** 2 for l in range(n_block))
    print(f"# S(t) computed in {time.time() - t4:.1f}s", flush=True)

    diff = S_pred - S_csv
    print()
    print(f"## N=7 path-6 verification vs bond-isolate CSV")
    print(f"# CSV: {csv_path.name}")
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0): pred={S_pred[0]:.6f}, csv={S_csv[0]:.6f}, expect (N-1)/N={(N-1)/N:.6f}")
    print()
    print("| t  | S_csv (RK4) | S_pred (closed form) | diff |")
    print("|---|---|---|---|")
    for i in [0, 30, 50, 100, 150, 200, 250, 300]:
        if i < len(t_csv):
            print(
                f"| {t_csv[i]:5.2f} | {S_csv[i]:.6f} | {S_pred[i]:.6f} | "
                f"{diff[i]:+.3e} |"
            )

    w = per_site_reduction_matrix(n_block)
    M = w @ R
    a = M * c[None, :]
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > 1e-12)[0]
    rates = -eigvals.real / gamma
    freqs = np.abs(eigvals.imag) / J
    groups = {}
    for k in contributing:
        key = (round(rates[k], 4), round(freqs[k], 4))
        groups[key] = groups.get(key, 0.0) + sig[k]

    print(f"\n## Path-6 mode-group survey")
    print(f"# 16384 total L_super modes, {len(contributing)} populated, {len(groups)} distinct (rate, |freq|) groups")
    print(f"# (Compare path-2..5: 4, 10, 12, 35 mode-groups respectively)")


if __name__ == "__main__":
    main()
