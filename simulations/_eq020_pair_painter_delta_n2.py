#!/usr/bin/env python3
"""EQ-020: pair-painter sees ΔN=2 coherences invisible to site-painter.

F70 (kinematic selection rule, proven): k-local observables see only
|ΔN|≤k sector coherences of ρ_0. Site-local sees |ΔN|≤1; pair-local sees
|ΔN|≤2.

Test: for initial states with specific ΔN-coherence content, compare
site-local vs pair-local Σ f. The ΔN=2 boost is the structural signature
of the pair-painter extension.

States tested (N=5):
  - (|vac> + |ψ_1>)/√2:  ΔN=1 (sine bonding, baseline)
  - (|vac> + |S_2>)/√2:  ΔN=2 — site invisible, pair visible
  - (|S_1> + |S_3>)/√2:  ΔN=2 — site invisible, pair visible
  - (|vac> + |S_3>)/√2:  ΔN=3 — invisible to BOTH site and pair (still)
  - (|vac> + |S_4>)/√2:  ΔN=4 — invisible to all k≤4
"""
from __future__ import annotations

import math
import sys
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
    per_pair_purity, evolve_full_traj, per_site_purity_from_traj,
)


N = 5
d = 2 ** N
J, GAMMA, DJ = 1.0, 0.05, 0.01
T_TOTAL, T_FIT = 80.0, 20.0
N_SAMPLES = 401


def dicke_state_ket(n_excitations):
    """|D(N, n)> = (1/sqrt(C(N,n))) Σ_{|x|=n} |x>, big-endian."""
    if n_excitations < 0 or n_excitations > N:
        return np.zeros(d, dtype=complex)
    v = np.zeros(d, dtype=complex)
    norm = 1.0 / math.sqrt(comb(N, n_excitations))
    for idx in range(d):
        if bin(idx).count("1") == n_excitations:
            v[idx] = norm
    return v


def vac_ket():
    v = np.zeros(d, dtype=complex)
    v[0] = 1.0
    return v


def normalize(v):
    return v / np.linalg.norm(v)


def main():
    print("=" * 80)
    print(f"EQ-020 pair-painter Σ f shift on ΔN-coherence states (N = {N})")
    print("=" * 80)
    print()

    print("Building L_A, V_L...")
    H_A = build_xy_chain_H(N, J)
    L_A = build_L(H_A, N, GAMMA)
    V_L = build_V_L_bond((0, 1), N)
    L_Bp = (L_A + DJ * V_L).tocsr()
    L_Bm = (L_A - DJ * V_L).tocsr()

    times = np.linspace(0, T_TOTAL, N_SAMPLES)
    sample_indices = np.arange(N_SAMPLES)

    vac = vac_ket()
    psi_1 = sine_mode_state(N, 1, d)
    S = {n: dicke_state_ket(n) for n in range(N + 1)}

    test_states = [
        ("(|vac>+|ψ_1>)/√2  [ΔN=1, sine bonding]",  normalize(vac + psi_1),  1),
        ("(|vac>+|S_1>)/√2  [ΔN=1, W bonding]",     normalize(vac + S[1]),    1),
        ("(|vac>+|S_2>)/√2  [ΔN=2, pair-only]",     normalize(vac + S[2]),    2),
        ("(|S_1>+|S_3>)/√2  [ΔN=2, pair-only]",     normalize(S[1] + S[3]),   2),
        ("(|vac>+|S_3>)/√2  [ΔN=3, INVISIBLE]",     normalize(vac + S[3]),    3),
        ("(|vac>+|S_4>)/√2  [ΔN=4, INVISIBLE]",     normalize(vac + S[4]),    4),
    ]

    print()
    print(f"  {'state':<45} {'ΔN':>3} {'Σ f_i (site)':>14} "
          f"{'Σ f_{ij} (pair)':>17} {'pair − site':>14}")
    print("  " + "-" * 95)

    for label, ket, dN in test_states:
        rho_0 = np.outer(ket, ket.conj())
        rho_0_flat = rho_0.flatten()

        traj_A = evolve_full_traj(L_A, rho_0_flat, times, dt_small=0.02)
        traj_Bp = evolve_full_traj(L_Bp, rho_0_flat, times, dt_small=0.02)
        traj_Bm = evolve_full_traj(L_Bm, rho_0_flat, times, dt_small=0.02)

        # Site
        P_A_s = per_site_purity_from_traj(traj_A, d, N, sample_indices)
        P_Bp_s = per_site_purity_from_traj(traj_Bp, d, N, sample_indices)
        P_Bm_s = per_site_purity_from_traj(traj_Bm, d, N, sample_indices)
        a_p_s = fit_alpha(P_A_s, P_Bp_s, times, t_max=T_FIT)
        a_m_s = fit_alpha(P_A_s, P_Bm_s, times, t_max=T_FIT)
        f_site = (np.sum(np.log(a_p_s)) - np.sum(np.log(a_m_s))) / (2 * DJ)

        # Pair
        P_A_p, _ = per_pair_purity(traj_A, d, N, sample_indices)
        P_Bp_p, _ = per_pair_purity(traj_Bp, d, N, sample_indices)
        P_Bm_p, _ = per_pair_purity(traj_Bm, d, N, sample_indices)
        a_p_p = fit_alpha(P_A_p, P_Bp_p, times, t_max=T_FIT)
        a_m_p = fit_alpha(P_A_p, P_Bm_p, times, t_max=T_FIT)
        f_pair = (np.sum(np.log(a_p_p)) - np.sum(np.log(a_m_p))) / (2 * DJ)

        diff = f_pair - f_site
        print(f"  {label:<45} {dN:>3d} {f_site:>+14.4f} "
              f"{f_pair:>+17.4f} {diff:>+14.4f}")

    print()
    print("Reading:")
    print("  ΔN=1 states: site & pair both see; values close (small shift).")
    print("  ΔN=2 states: pair-painter sees a NEW coherence block; large shift.")
    print("  ΔN=3, 4 states: invisible to both site and pair-painter (k ≤ 2).")


if __name__ == "__main__":
    main()
