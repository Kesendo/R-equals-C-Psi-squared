#!/usr/bin/env python3
"""EQ-020 triple-painter (k=3) chiral mirror law verification.

Following pair-painter (k=2) verification at N=5, N=7: Σ f_{ij}(ψ_k) =
Σ f_{ij}(ψ_{N+1−k}) machine-precision. Same structural argument should give
chiral mirror at k=3 (triple-painter) and any k:

    Σ f_{i_1 ... i_k}(ψ_k) = Σ f_{i_1 ... i_k}(ψ_{N+1-k})

Mechanism: K_1 ρ K_1† gives ρ_{sites}^K = (∏_s Z_s) ρ_{sites} (∏_s Z_s),
unitary similarity → Tr(ρ_{sites}²) is K_1-invariant. Sine-mode bonding
states map K_1: (vac+ψ_k)/√2 ↔ (vac+ψ_{N+1−k})/√2, so per-triple α and
hence Σ f are equal across the mirror.

Verify at N=5 (10 triples) and N=7 (35 triples).

Also test F70-style visibility: at N=5, the (vac+|S_3>)/√2 state has
|ΔN|=3 coherence, invisible to pair (k=2) but visible to triple (k=3).
Compare site-vs-pair-vs-triple Σ f for this state.
"""
from __future__ import annotations

import math
import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import scipy.sparse as sps

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from _eq014_chiral_mirror_multi_N import (
    build_xy_chain_H, build_L, build_V_L_bond,
    sine_mode_state, rk4_step, fit_alpha,
)
from _eq020_pair_painter_chiral_mirror import (
    evolve_full_traj, per_site_purity_from_traj, per_pair_purity,
)


def k_local_reduced(rho_dxd, sites, N):
    """Trace out all sites except `sites`. Returns 2^k × 2^k matrix."""
    sites = tuple(sorted(sites))
    k = len(sites)
    shape = [2] * (2 * N)
    T = rho_dxd.reshape(shape)
    letters = "abcdefghijklmnopqrstuvwxyz"
    row = list(letters[:N])
    col = list(letters[N:2 * N])
    for s in range(N):
        if s not in sites:
            col[s] = row[s]
    in_spec = "".join(row) + "".join(col)
    out_spec = "".join(row[s] for s in sites) + "".join(col[s] for s in sites)
    out = np.einsum(f"{in_spec}->{out_spec}", T)
    d_sub = 2 ** k
    return out.reshape(d_sub, d_sub)


def per_triple_purity(rho_traj_flat, d, N, sample_indices):
    triples = list(combinations(range(N), 3))
    n_t = len(triples)
    n_s = len(sample_indices)
    out = np.zeros((n_t, n_s))
    for s_idx, t_idx in enumerate(sample_indices):
        rho_dxd = rho_traj_flat[t_idx].reshape(d, d)
        for p_idx, sites in enumerate(triples):
            rho_s = k_local_reduced(rho_dxd, sites, N)
            out[p_idx, s_idx] = float(np.real(np.trace(rho_s @ rho_s)))
    return out, triples


def bonding_mode_density_flat(N, k, d):
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = sine_mode_state(N, k, d)
    phi = (vac + psi) / np.sqrt(2)
    rho = np.outer(phi, phi.conj())
    return rho.flatten()


def run_triple_chiral_mirror(N, J=1.0, gamma=0.05, delta_J=0.01,
                                t_total=None, t_fit=20.0, n_samples=None,
                                bond=(0, 1), dt_small=0.02):
    if t_total is None:
        t_total = 30.0 if N >= 7 else 80.0
    if n_samples is None:
        n_samples = 151 if N >= 7 else 401

    print(f"\n=== N = {N} triple-painter chiral mirror ===")
    d = 2 ** N
    n_triples = comb(N, 3)
    print(f"  {n_triples} triples; d² = {d * d}")

    print(f"  Building L_A, L_B...")
    H_A = build_xy_chain_H(N, J)
    L_A = build_L(H_A, N, gamma)
    V_L = build_V_L_bond(bond, N)
    L_Bp = (L_A + delta_J * V_L).tocsr()
    L_Bm = (L_A - delta_J * V_L).tocsr()
    times = np.linspace(0, t_total, n_samples)
    sample_indices = np.arange(n_samples)

    print(f"\n  {'k':>2} {'k_mirror':>8} {'Σ f_i (site)':>14} "
          f"{'Σ f_{ij} (pair)':>17} {'Σ f_{ijk} (triple)':>20} {'time':>7}")
    print("  " + "-" * 75)

    results = []
    for k in range(1, N + 1):
        rho0 = bonding_mode_density_flat(N, k, d)
        t0 = time.time()
        traj_A = evolve_full_traj(L_A, rho0, times, dt_small=dt_small)
        traj_Bp = evolve_full_traj(L_Bp, rho0, times, dt_small=dt_small)
        traj_Bm = evolve_full_traj(L_Bm, rho0, times, dt_small=dt_small)

        # Site
        P_A_s = per_site_purity_from_traj(traj_A, d, N, sample_indices)
        P_Bp_s = per_site_purity_from_traj(traj_Bp, d, N, sample_indices)
        P_Bm_s = per_site_purity_from_traj(traj_Bm, d, N, sample_indices)
        a_p_s = fit_alpha(P_A_s, P_Bp_s, times, t_max=t_fit)
        a_m_s = fit_alpha(P_A_s, P_Bm_s, times, t_max=t_fit)
        f_site = (np.sum(np.log(a_p_s)) - np.sum(np.log(a_m_s))) / (2 * delta_J)

        # Pair
        P_A_p, _ = per_pair_purity(traj_A, d, N, sample_indices)
        P_Bp_p, _ = per_pair_purity(traj_Bp, d, N, sample_indices)
        P_Bm_p, _ = per_pair_purity(traj_Bm, d, N, sample_indices)
        a_p_p = fit_alpha(P_A_p, P_Bp_p, times, t_max=t_fit)
        a_m_p = fit_alpha(P_A_p, P_Bm_p, times, t_max=t_fit)
        f_pair = (np.sum(np.log(a_p_p)) - np.sum(np.log(a_m_p))) / (2 * delta_J)

        # Triple
        P_A_t, _ = per_triple_purity(traj_A, d, N, sample_indices)
        P_Bp_t, _ = per_triple_purity(traj_Bp, d, N, sample_indices)
        P_Bm_t, _ = per_triple_purity(traj_Bm, d, N, sample_indices)
        a_p_t = fit_alpha(P_A_t, P_Bp_t, times, t_max=t_fit)
        a_m_t = fit_alpha(P_A_t, P_Bm_t, times, t_max=t_fit)
        f_triple = (np.sum(np.log(a_p_t)) - np.sum(np.log(a_m_t))) / (2 * delta_J)

        elapsed = time.time() - t0
        k_mirror = N + 1 - k
        results.append({
            "N": N, "k": k, "k_mirror": k_mirror,
            "f_site": float(f_site), "f_pair": float(f_pair),
            "f_triple": float(f_triple), "elapsed": elapsed,
        })
        print(f"  {k:>2d} {k_mirror:>8d} {f_site:>+14.4f} "
              f"{f_pair:>+17.4f} {f_triple:>+20.4f} {elapsed:>5.1f}s")

    print()
    print(f"  Chiral mirror law verification (TRIPLE painter):")
    for k in range(1, (N + 1) // 2 + 1):
        k_m = N + 1 - k
        if k > k_m:
            continue
        f_k = next(r["f_triple"] for r in results if r["k"] == k)
        f_m = next(r["f_triple"] for r in results if r["k"] == k_m)
        diff = abs(f_k - f_m)
        verdict = "OK" if diff < 1e-3 else f"BREAK ({diff:.2e})"
        if k == k_m:
            print(f"    k={k}: chiral fixed point; Σ f_{{ijk}} = {f_k:+.4f}")
        else:
            print(f"    k={k} ↔ k={k_m}: |Δ Σ f_{{ijk}}| = {diff:.2e}  {verdict}")
    return results


def main():
    print("=" * 80)
    print("EQ-020 triple-painter (k=3): chiral mirror law on Σ f_{ijk}")
    print("=" * 80)

    Ns = [int(x) for x in sys.argv[1:] if x.isdigit()] or [5]
    for N in Ns:
        run_triple_chiral_mirror(N)


if __name__ == "__main__":
    main()
