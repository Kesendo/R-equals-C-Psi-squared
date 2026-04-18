#!/usr/bin/env python3
"""n7_perspectival_extended_states.py

Part C state-independence stress-test for the Perspectival Time Field
closure law Sum_i ln(alpha_i) = 0.

Runs the N=7 XY chain under bond-(0,1) defect with three additional
initial states at J_mod in {0.9, 1.0, 1.1}:
  - psi_3 (k=3, nodes at sites 1 and 5)
  - psi_4 (k=4, node at site 3)
  - |+>^7 (all-excitation-sector uniform superposition)

Writes per-state scans to separate subdirs of
    simulations/results/n7_coupling_defect_overlay_extended/{psi3, psi4, plus7}/
so that observer_time_rescale.py can fit alpha_i and compute
Sum_i ln(alpha_i) per state.

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

ROOT = Path(__file__).parent / "results" / "n7_coupling_defect_overlay_extended"
J_MOD_VALUES = [0.9, 1.0, 1.1]


def initial_rho_psi_k_bonding(N, k):
    phi = vacuum(N) + single_excitation_mode(N, k)
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


def initial_rho_plus_N(N):
    phi = np.ones(2**N, dtype=complex) / np.sqrt(2**N)
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


def run_state_scan(subdir, init_fn, label_state):
    out = ROOT / subdir
    out.mkdir(parents=True, exist_ok=True)
    times = np.arange(N_STEPS + 1) * DT
    hamming = build_hamming_matrix(N)
    np.save(out / "times.npy", times)

    print(f"--- {label_state} init -> {out} ---")
    # Control A
    t0 = time.time()
    purity, c, m, l = propagate([J_UNIFORM] * (N - 1), init_fn, times, hamming)
    np.savez(out / "experiment_A.npz", purity=purity, cpsi_0k=c, mi_0k=m, l1=l)
    print(f"  A: {time.time() - t0:.2f} s, purity(t=0) = {purity[0]}")

    for J_mod in J_MOD_VALUES:
        t0 = time.time()
        J_list = [J_mod] + [J_UNIFORM] * (N - 2)
        purity, c, m, l = propagate(J_list, init_fn, times, hamming)
        np.savez(out / f"experiment_B_{J_mod:g}.npz",
                 purity=purity, cpsi_0k=c, mi_0k=m, l1=l)
        print(f"  B(J_mod={J_mod}): {time.time() - t0:.2f} s")


def main():
    run_state_scan("psi3",
                   lambda N: initial_rho_psi_k_bonding(N, 3),
                   "psi_3 (k=3)")
    run_state_scan("psi4",
                   lambda N: initial_rho_psi_k_bonding(N, 4),
                   "psi_4 (k=4)")
    run_state_scan("plus7",
                   initial_rho_plus_N,
                   "|+>^7")


if __name__ == "__main__":
    main()
