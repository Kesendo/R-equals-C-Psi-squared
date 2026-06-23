#!/usr/bin/env python3
"""F76 verification: the 0.93 envelope is pure dephasing, not Heisenberg mixing.

Computes MM(t)/MM(0) via single-excitation-sector Lindblad and compares with
the pure-dephasing-only prediction (no Hamiltonian mixing). Agreement within
0.5 percent across N=5..13, k=1..5 at gamma_0 = 0.05, t = 0.1.

The single-excitation sector is only N x N in density-matrix indices (not
2^N x 2^N), making the Lindblad cheap: L is an N^2 x N^2 matrix that can
be diagonalised once per (N, gamma_0).

See docs/ANALYTICAL_FORMULAS.md F76 for the closed-form envelope formula.
"""
import sys
import math
import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def single_excitation_H(N, J=1.0):
    """N x N Heisenberg single-excitation Hamiltonian: hopping + boundary ZZ."""
    H = np.zeros((N, N), dtype=complex)
    for b in range(N - 1):
        H[b, b + 1] += 2.0 * J
        H[b + 1, b] += 2.0 * J
    for j in range(N):
        adj = 1 if j == 0 or j == N - 1 else 2
        H[j, j] = ((N - 1) - 2 * adj) * J
    return H


def single_excitation_L(N, gamma_0, J=1.0):
    """N^2 x N^2 Liouvillian on vectorised rho^(1)."""
    H = single_excitation_H(N, J)
    L = np.zeros((N * N, N * N), dtype=complex)
    for i in range(N):
        for j in range(N):
            idx_ij = i * N + j
            for k in range(N):
                L[idx_ij, k * N + j] += -1j * H[i, k]
                L[idx_ij, i * N + k] -= -1j * H[k, j]
            if i != j:
                L[idx_ij, idx_ij] += -4.0 * gamma_0
    return L


def bonding_k_rho1(N, k):
    """Initial rho^(1) for |psi_k> = sqrt(2/(N+1)) sum_j sin(pi k (j+1)/(N+1)) |1_j>."""
    c = np.array([
        math.sqrt(2.0 / (N + 1)) * math.sin(math.pi * k * (j + 1) / (N + 1))
        for j in range(N)
    ])
    return np.outer(c, c.conj())


def pure_dephasing_rho1(rho1_0, gamma_0, t):
    """rho^(1)(t) under pure Z-dephasing (no H): coherences decay at 4 gamma_0."""
    lam = math.exp(-4 * gamma_0 * t)
    rho = rho1_0.copy()
    for i in range(rho.shape[0]):
        for j in range(rho.shape[0]):
            if i != j:
                rho[i, j] *= lam
    return rho


def h_entropy(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def mirror_pair_mi(rho1, N, ell):
    """MI(site ell, site N-1-ell) from single-excitation rho^(1)."""
    m = N - 1 - ell
    p_ell = np.real(rho1[ell, ell])
    p_m = np.real(rho1[m, m])
    c_em = rho1[ell, m]
    tr = p_ell + p_m
    disc = math.sqrt(max(0.0, ((p_ell - p_m) / 2) ** 2 + abs(c_em) ** 2))
    e1 = max(0.0, 1 - tr)
    e2 = tr / 2 + disc
    e3 = max(0.0, tr / 2 - disc)
    S_ab = 0.0
    for e in (e1, e2, e3):
        if e > 1e-12:
            S_ab -= e * math.log2(e)
    return h_entropy(p_ell) + h_entropy(p_m) - S_ab


def total_mm(rho1, N):
    return sum(mirror_pair_mi(rho1, N, ell) for ell in range(N // 2))


def main() -> None:
    gamma_0 = 0.05
    t = 0.1

    print("F76 verification: MM(t=0.1)/MM(0) vs pure-dephasing prediction")
    print(f"gamma_0 = {gamma_0}, J = 1, t = {t}")
    print()
    print(f"{'N':>3} {'k':>3} | {'MM(0)':>8} {'MM_full/MM(0)':>15} {'MM_deph/MM(0)':>15} {'diff %':>8}")
    print("-" * 60)

    for N in (5, 7, 9, 11, 13):
        for k in (1, 2, 3, 4, 5):
            if k > N:
                continue
            L = single_excitation_L(N, gamma_0)
            rho1_0 = bonding_k_rho1(N, k)
            mm_0 = total_mm(rho1_0, N)

            v_t_full = expm(L * t) @ rho1_0.flatten()
            rho1_t_full = np.real(v_t_full.reshape(N, N))
            mm_full = total_mm(rho1_t_full, N)
            r_full = mm_full / mm_0

            rho1_t_deph = pure_dephasing_rho1(rho1_0, gamma_0, t)
            mm_deph = total_mm(rho1_t_deph, N)
            r_deph = mm_deph / mm_0

            diff_pct = abs(r_full - r_deph) / r_full * 100
            print(f"{N:>3} {k:>3} | {mm_0:>8.4f} {r_full:>15.4f} {r_deph:>15.4f} {diff_pct:>8.3f}")


if __name__ == "__main__":
    main()
