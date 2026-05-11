"""F89 (vac, SE) self-contribution closed form via Parseval orthogonality.

S^(vac,SE)_block(t; k, N) = (k+1)·(N-k-1)²/(N²·(N-1))·exp(-4γ₀ t)

Verification script: compare closed form against numerical (vac, SE)-only
contribution (computed by zeroing the (SE, DE) sector of ρ_block(0)).

See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Path-k (vac, SE)
self-contribution" for derivation and smooth-backbone application.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from _f89_pathk_lib import (
    block_bit_pos,
    build_block_L,
    reduce_block_to_site_01,
    state_idx,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]


def S_vac_se_closed_form(k: int, N: int, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S^(vac,SE)_block(t; k, N) = (k+1)·(N-k-1)²/(N²·(N-1))·exp(-4γ₀ t)."""
    n_block = k + 1
    n_E = N - n_block
    prefactor = n_block * n_E * n_E / (N * N * (N - 1))
    return prefactor * np.exp(-4.0 * gamma * t_array)


def numerical_vac_se_block(k: int, N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """(vac, SE)-only contribution: zero out (SE, DE) terms in ρ_block(0)."""
    n_block = k + 1
    D = 2**n_block
    L = build_block_L(J, gamma, n_block)
    bit_pos = block_bit_pos(n_block)
    N_E = N - n_block
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)

    rho = np.zeros((D, D), dtype=complex)
    for j in range(n_block):
        bits = [0] * n_block
        bits[j] = 1
        rho[0, state_idx(bits, bit_pos)] += pre * N_E
    rho = (rho + rho.conj().T) / 2.0

    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)

    S = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S[ti] = sum(2.0 * abs(reduce_block_to_site_01(rho_t, l, n_block)) ** 2 for l in range(n_block))

    return S


def main() -> None:
    J, gamma = 0.075, 0.05
    t_array = np.linspace(0, 30, 301)

    print(f"# F89 (vac, SE) closed form via Parseval, J={J}, γ={gamma}\n")
    print("# Closed form: S^(vac,SE)_block(t; k, N) = (k+1)(N-k-1)²/(N²(N-1)) · exp(-4γ₀ t)\n")

    print("## Verification: closed form vs numerical (with (SE, DE) zeroed out)")
    print("| k | N | (k+1)(N-k-1)²/(N²(N-1)) | max |closed - numerical| |")
    print("|---|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        for N in [k + 2, k + 4, k + 6]:  # spans different N_E
            if k + 1 > N:
                continue
            S_closed = S_vac_se_closed_form(k, N, gamma, t_array)
            S_numerical = numerical_vac_se_block(k, N, J, gamma, t_array)
            max_diff = np.max(np.abs(S_closed - S_numerical))
            prefactor = (k + 1) * (N - k - 1) ** 2 / (N * N * (N - 1))
            print(f"| {k} | {N} | {prefactor:.6f} | {max_diff:.3e} |")

    print()
    print("## All differences should be at machine precision (~1e-14).")
    print("# Parseval orthogonality is exact, so the closed form is bit-exact.")


if __name__ == "__main__":
    main()
