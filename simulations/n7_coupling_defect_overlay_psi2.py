#!/usr/bin/env python3
"""ψ_2 initial-state overlay scan for the N=7 coupling defect.

Runs the same Lindblad RK4 pipeline as n7_coupling_defect_overlay.py but
with the initial state replaced by phi = (|vac> + |psi_2>)/sqrt(2), where
psi_2 is the k=2 single-excitation eigenmode of the uniform XY chain.
psi_2 has a node at site 3 (centre) and anti-nodes at sites 1 and 5.

Defect: bond (0, 1) with J_mod in {0.5, 0.9, 1.1, 1.5}.

This feeds the observer-time falsification test:
  * "intrinsic observer time" predicts the same f_i pattern as psi_1.
  * "slow-mode emergent rescaling" predicts a different f_i pattern,
    in particular the zero-sites (currently 3 and 6 for psi_1) should
    relocate to sites that mirror the psi_2 structure.

Date: 2026-04-18
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from n7_coupling_defect_overlay import (
    N, GAMMA_0, T_FINAL, N_STEPS, DT, J_UNIFORM,
    build_H_XY, build_hamming_matrix,
    rk4_step,
    site_purity_all, pair_cpsi_and_mi_for_refsite, global_l1,
    single_excitation_mode, vacuum,
)

RESULTS_DIR = (Path(__file__).parent / "results"
               / "n7_coupling_defect_overlay_extended" / "psi2_init")

J_MOD_VALUES = [0.5, 0.9, 1.1, 1.5]


def initial_rho_psi2_bonding(N):
    """rho_0 = |phi><phi| with phi = (|vac> + |psi_2>)/sqrt(2)."""
    phi = vacuum(N) + single_excitation_mode(N, 2)
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


def propagate(J_list, init_fn, times, hamming):
    H = build_H_XY(J_list, N)
    rho = init_fn(N)
    T = len(times)
    purity = np.empty((T, N))
    cpsi_0k = np.empty((T, N - 1))
    mi_0k = np.empty((T, N - 1))
    l1 = np.empty(T)
    purity[0] = site_purity_all(rho, N)
    cpsi_0k[0], mi_0k[0] = pair_cpsi_and_mi_for_refsite(rho, 0, N)
    l1[0] = global_l1(rho)
    for step in range(1, T):
        rho = rk4_step(rho, H, hamming, GAMMA_0, DT)
        purity[step] = site_purity_all(rho, N)
        cpsi_0k[step], mi_0k[step] = pair_cpsi_and_mi_for_refsite(rho, 0, N)
        l1[step] = global_l1(rho)
    return purity, cpsi_0k, mi_0k, l1


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    times = np.arange(N_STEPS + 1) * DT
    hamming = build_hamming_matrix(N)
    np.save(RESULTS_DIR / "times.npy", times)

    print(f"ψ_2 initial-state overlay, N = {N}, gamma_0 = {GAMMA_0}")
    print(f"  T = {T_FINAL}, steps = {N_STEPS}, out = {RESULTS_DIR}")

    # Control A
    t0 = time.time()
    purity, cpsi, mi, l1 = propagate([J_UNIFORM] * (N - 1),
                                     initial_rho_psi2_bonding, times, hamming)
    np.savez(RESULTS_DIR / "experiment_A.npz",
             purity=purity, cpsi_0k=cpsi, mi_0k=mi, l1=l1)
    print(f"  A (uniform J): {time.time() - t0:.2f} s, "
          f"purity(t=0) = {purity[0]}")

    for J_mod in J_MOD_VALUES:
        t0 = time.time()
        J_list = [J_mod] + [J_UNIFORM] * (N - 2)
        purity, cpsi, mi, l1 = propagate(J_list, initial_rho_psi2_bonding,
                                         times, hamming)
        label = f"{J_mod:g}"
        np.savez(RESULTS_DIR / f"experiment_B_{label}.npz",
                 purity=purity, cpsi_0k=cpsi, mi_0k=mi, l1=l1)
        print(f"  B (J_mod = {J_mod}): {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
