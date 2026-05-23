"""F86e: σ_0(N) for N=4..20 with the 2.86 contract subtracted, the naked change.

Tom 2026-05-20: strip the shared leading part (2.86) from every σ_0(N),
look at the change part alone.
"""
from __future__ import annotations

import sys
import time

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402


def hd_subspace_projector(N: int, n: int, hd_value: int) -> np.ndarray:
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    cols = []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
    return (np.column_stack(cols) if cols else
            np.zeros((len(P_n) * Mnp1, 0), dtype=complex))


def sigma_0(N: int, n: int = 1, gamma_0: float = 0.05) -> float:
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P_HD1 = hd_subspace_projector(N, n, 1)
    P_HD3 = hd_subspace_projector(N, n, 3)
    if P_HD3.shape[1] == 0 or P_HD1.shape[1] == 0:
        return float("nan")
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    return float(np.linalg.svd(V_inter, compute_uv=False)[0])


def main() -> None:
    contract = 2.86
    print("=" * 78)
    print(f"σ_0(N) for N=4..20, contract {contract} subtracted, the naked change")
    print("=" * 78)
    print(f"{'N':>3}  {'σ_0(N)':>16}  {'σ_0 − 2.86':>14}  {'(σ_0−2.86)·10⁴':>16}")
    print("-" * 78)
    naked = {}
    for N in range(4, 21):
        t0 = time.time()
        s = sigma_0(N)
        nk = s - contract
        naked[N] = nk
        elapsed = time.time() - t0
        flag = "" if elapsed < 2 else f"  ({elapsed:.1f}s)"
        print(f"{N:>3}  {s:>16.10f}  {nk:>+14.8f}  {nk*1e4:>+16.4f}{flag}")

    print("-" * 78)
    print()
    print("Naked change ratios (consecutive, step 1):")
    Ns = sorted(naked)
    for i in range(1, len(Ns)):
        N = Ns[i]
        if abs(naked[Ns[i-1]]) > 1e-12:
            r = naked[N] / naked[Ns[i-1]]
            print(f"  N={Ns[i-1]:2d}→{N:2d}: ratio = {r:+.6f}")

    print()
    print("Naked change ratios (step 2):")
    for N in Ns:
        if N + 2 in naked and abs(naked[N]) > 1e-12:
            r = naked[N + 2] / naked[N]
            print(f"  N={N:2d}→{N+2:2d}: ratio = {r:+.6f}")


if __name__ == "__main__":
    main()
