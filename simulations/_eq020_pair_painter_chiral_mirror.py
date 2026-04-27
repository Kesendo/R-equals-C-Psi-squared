#!/usr/bin/env python3
"""EQ-020: pair-painter PTF extension — chiral mirror law on pair-purity α_{ij}.

Site-local PTF: per-site purity P_i = Tr(ρ_i²), fit α_i so P_B(t) ≈ P_A(α_i t),
closure Σ_i ln α_i. Site-local sees |ΔN|≤1 sector blocks of ρ_0 (F70).

Pair-local: per-pair purity P_{ij} = Tr(ρ_{ij}²), C(N,2) pairs, fit α_{ij},
closure Σ_{(i,j)} ln α_{ij}. Pair-local sees |ΔN|≤2 (F70 generalization).

EQ-014 site-local: Σ f_i obeys the chiral mirror law Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}),
proved structurally from K_1 ψ_k = ψ_{N+1−k} (no sign) and K_1-invariance of
per-site purity.

This script asks: does the chiral mirror law extend to pair-painter Σ f_{ij}?

Mechanism check: K_1 = ⊗_i (-1)^i is a product of single-site Z's. For any
pair-reduced ρ_{ij} = Tr_{¬{i,j}}(ρ), the K_1-rotated state has (K_1)_{i}(K_1)_{j}
applied at the kept sites — i.e. ρ_{ij}^K = (Z_i Z_j) ρ_{ij} (Z_i Z_j). This
is a unitary similarity, so Tr(ρ_{ij}²) = Tr((ρ_{ij}^K)²) — pair-purity is
ALSO K_1-invariant. Hence α_{ij} and Σ f_{ij} should be K_1-mirror-symmetric
under k ↔ N+1-k.

Verification protocol: at N=5 (and N=7 if time allows), for each ψ_k
bonding state |ψ_k_bonding⟩ = (|vac⟩+|ψ_k⟩)/√2:
  - propagate under L_A (uniform), L_B+ (J + δJ at bond (0,1)), L_B-
  - compute α_{ij} for each pair (10 pairs at N=5, 21 at N=7)
  - compute Σ ln α_{ij}, divide by 2δJ → Σ f_{ij}(ψ_k)
  - test mirror law: Σ f_{ij}(k) =? Σ f_{ij}(N+1-k)

Compare against site-local Σ f_i from same trajectories.
"""
from __future__ import annotations

import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

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


def bonding_mode_density_flat(N, k, d):
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = sine_mode_state(N, k, d)
    phi = (vac + psi) / np.sqrt(2)
    rho = np.outer(phi, phi.conj())
    return rho.flatten()


def partial_trace_keep_pair(rho_dxd, i, j, N):
    """Trace out all sites except (i, j); return 4×4 ρ_{ij}.

    Big-endian convention: site i ↔ axis i in shape [2]*2N.
    """
    shape = [2] * (2 * N)
    T = rho_dxd.reshape(shape)
    letters = "abcdefghijklmnop"
    row_labels = list(letters[:N])
    col_labels = list(letters[N:2 * N])
    # Trace any site not in {i, j}: set col_label = row_label
    for s in range(N):
        if s != i and s != j:
            col_labels[s] = row_labels[s]
    in_spec = "".join(row_labels) + "".join(col_labels)
    out_spec = row_labels[i] + row_labels[j] + col_labels[i] + col_labels[j]
    out = np.einsum(f"{in_spec}->{out_spec}", T)
    return out.reshape(4, 4)


def per_pair_purity(rho_traj_flat, d, N, sample_indices):
    """Pair purities over time. Returns (n_pairs, n_samples)."""
    pairs = list(combinations(range(N), 2))
    n_pairs = len(pairs)
    n_samples = len(sample_indices)
    out = np.zeros((n_pairs, n_samples))
    for s_idx, t_idx in enumerate(sample_indices):
        rho_dxd = rho_traj_flat[t_idx].reshape(d, d)
        for p_idx, (i, j) in enumerate(pairs):
            rho_ij = partial_trace_keep_pair(rho_dxd, i, j, N)
            out[p_idx, s_idx] = float(np.real(np.trace(rho_ij @ rho_ij)))
    return out, pairs


def evolve_full_traj(L, rho0_flat, sample_times, dt_small=0.02):
    """RK4 evolve and return full ρ(t) trajectory as (n_samples, d²)."""
    sample_times = np.asarray(sample_times)
    n_samples = len(sample_times)
    d2 = rho0_flat.shape[0]
    out = np.zeros((n_samples, d2), dtype=complex)
    current_rho = rho0_flat.astype(complex).copy()
    t_current = 0.0
    for si, t_target in enumerate(sample_times):
        dt_total = t_target - t_current
        if dt_total > 0:
            n_steps = max(1, int(np.ceil(dt_total / dt_small)))
            dt = dt_total / n_steps
            for _ in range(n_steps):
                current_rho = rk4_step(L, current_rho, dt)
            t_current = t_target
        out[si] = current_rho
    return out


def per_site_purity_from_traj(rho_traj, d, N, sample_indices):
    """Per-site purity along trajectory. Returns (N, n_samples)."""
    n_samples = len(sample_indices)
    out = np.zeros((N, n_samples))
    shape = [2] * (2 * N)
    letters = "abcdefghijklmnop"
    for s_idx, t_idx in enumerate(sample_indices):
        rho_dxd = rho_traj[t_idx].reshape(d, d)
        T = rho_dxd.reshape(shape)
        for i in range(N):
            row = list(letters[:N]); col = list(letters[N:2 * N])
            for s in range(N):
                if s != i:
                    col[s] = row[s]
            in_spec = "".join(row) + "".join(col)
            out_spec = row[i] + col[i]
            rho_i = np.einsum(f"{in_spec}->{out_spec}", T)
            out[i, s_idx] = float(np.real(np.trace(rho_i @ rho_i)))
    return out


def run_N(N, J=1.0, gamma=0.05, delta_J=0.01,
           t_total=None, t_fit=20.0, n_samples=None, bond=(0, 1),
           dt_small=0.02):
    if t_total is None:
        t_total = 30.0 if N >= 7 else 80.0
    if n_samples is None:
        n_samples = 151 if N >= 7 else 401

    print(f"\n=== N = {N} pair-painter chiral mirror ===")
    d = 2 ** N
    print(f"  Building L_A, L_B (d²={d*d})...")
    H_A = build_xy_chain_H(N, J)
    L_A = build_L(H_A, N, gamma)
    V_L = build_V_L_bond(bond, N)
    L_Bp = (L_A + delta_J * V_L).tocsr()
    L_Bm = (L_A - delta_J * V_L).tocsr()
    print(f"  L nnz: {L_A.nnz}")
    times = np.linspace(0, t_total, n_samples)
    sample_indices = np.arange(n_samples)
    n_pairs = N * (N - 1) // 2

    print(f"\n  {'k':>2} {'k_mirror':>8} {'Σ f_i (site)':>14} "
          f"{'Σ f_{ij} (pair)':>17} {'#pairs':>7} {'time':>7}")
    print("  " + "-" * 60)
    results = []
    for k in range(1, N + 1):
        rho0 = bonding_mode_density_flat(N, k, d)
        t0 = time.time()
        traj_A = evolve_full_traj(L_A, rho0, times, dt_small=dt_small)
        traj_Bp = evolve_full_traj(L_Bp, rho0, times, dt_small=dt_small)
        traj_Bm = evolve_full_traj(L_Bm, rho0, times, dt_small=dt_small)

        # Site-local
        P_A_site = per_site_purity_from_traj(traj_A, d, N, sample_indices)
        P_Bp_site = per_site_purity_from_traj(traj_Bp, d, N, sample_indices)
        P_Bm_site = per_site_purity_from_traj(traj_Bm, d, N, sample_indices)
        alpha_p_site = fit_alpha(P_A_site, P_Bp_site, times, t_max=t_fit)
        alpha_m_site = fit_alpha(P_A_site, P_Bm_site, times, t_max=t_fit)
        valid_site = (np.all(np.isfinite(alpha_p_site) & (alpha_p_site > 0))
                      and np.all(np.isfinite(alpha_m_site) & (alpha_m_site > 0)))
        if valid_site:
            sum_ln_p = np.sum(np.log(alpha_p_site))
            sum_ln_m = np.sum(np.log(alpha_m_site))
            f_site = (sum_ln_p - sum_ln_m) / (2 * delta_J)
        else:
            f_site = float('nan')

        # Pair-local
        P_A_pair, pairs_list = per_pair_purity(traj_A, d, N, sample_indices)
        P_Bp_pair, _ = per_pair_purity(traj_Bp, d, N, sample_indices)
        P_Bm_pair, _ = per_pair_purity(traj_Bm, d, N, sample_indices)
        alpha_p_pair = fit_alpha(P_A_pair, P_Bp_pair, times, t_max=t_fit)
        alpha_m_pair = fit_alpha(P_A_pair, P_Bm_pair, times, t_max=t_fit)
        valid_pair = (np.all(np.isfinite(alpha_p_pair) & (alpha_p_pair > 0))
                      and np.all(np.isfinite(alpha_m_pair) & (alpha_m_pair > 0)))
        if valid_pair:
            sum_ln_p_pair = np.sum(np.log(alpha_p_pair))
            sum_ln_m_pair = np.sum(np.log(alpha_m_pair))
            f_pair = (sum_ln_p_pair - sum_ln_m_pair) / (2 * delta_J)
        else:
            f_pair = float('nan')

        elapsed = time.time() - t0
        k_mirror = N + 1 - k
        print(f"  {k:>2d} {k_mirror:>8d} {f_site:>+14.4f} "
              f"{f_pair:>+17.4f} {n_pairs:>7d} {elapsed:>5.1f}s")
        results.append({
            "N": N, "k": k, "k_mirror": k_mirror,
            "f_site": float(f_site), "f_pair": float(f_pair),
            "elapsed": elapsed,
        })

    print()
    print(f"  Chiral mirror law verification (PAIR painter):")
    print(f"  Σ f_{{ij}}(k) =? Σ f_{{ij}}(N+1-k)")
    for k in range(1, (N + 1) // 2 + 1):
        k_m = N + 1 - k
        if k > k_m:
            continue
        f_k = next(r["f_pair"] for r in results if r["k"] == k)
        f_m = next(r["f_pair"] for r in results if r["k"] == k_m)
        diff = abs(f_k - f_m)
        verdict = "OK" if diff < 1e-3 else f"BREAK ({diff:.2e})"
        if k == k_m:
            print(f"    k={k}: chiral fixed point; Σ f_{{ij}} = {f_k:+.4f}")
        else:
            print(f"    k={k} ↔ k={k_m}: |Δ Σ f_{{ij}}| = {diff:.2e}  {verdict}")
    return results


def main():
    print("=" * 80)
    print("EQ-020 pair-painter PTF: chiral mirror law on Σ f_{ij}")
    print("=" * 80)
    Ns = [int(x) for x in sys.argv[1:] if x.isdigit()] or [5]
    for N in Ns:
        run_N(N)


if __name__ == "__main__":
    main()
