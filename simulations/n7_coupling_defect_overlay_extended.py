#!/usr/bin/env python3
"""N = 7 coupling-defect overlay — extended simulation driver.

Extends `n7_coupling_defect_overlay.py` with three scan families that
feed the observer-time rescale analysis:

  1. Finer J_mod grid at the default defect bond (0, 1):
     {0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0}
  2. Defect-location scan at J_mod in {0.5, 2.0} on bonds (1,2)..(5,6).
  3. Localized initial state |1_3> (centre excitation) at the default
     defect (0, 1) with J_mod in {0.5, 1.5, 2.0}.

Each variant is saved as a .npz in its own subdirectory of
    simulations/results/n7_coupling_defect_overlay_extended/
so that observer_time_rescale.py can load each scan-family independently.

The physics kernel (H_XY, dephasing as Hadamard-damping, RK4) is
imported from the base script.

Date: 2026-04-17 / 2026-04-18 follow-up
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

# Reuse kernels from the base script.
sys.path.insert(0, str(Path(__file__).parent))
from n7_coupling_defect_overlay import (
    N, GAMMA_0, T_FINAL, N_STEPS, DT, J_UNIFORM,
    build_H_XY, build_hamming_matrix,
    initial_rho_psi1_bonding,
    rk4_step,
    site_purity_all, pair_cpsi_and_mi_for_refsite, global_l1,
    single_excitation_mode, vacuum,
)

RESULTS_ROOT = (Path(__file__).parent / "results"
                / "n7_coupling_defect_overlay_extended")
RESULTS_ROOT.mkdir(parents=True, exist_ok=True)


def initial_rho_local_excitation(N, site):
    """rho_0 = |phi><phi| with phi = (|vac> + |1_site>)/sqrt(2).
    Big-endian: |1_site> has state index 2^(N-1-site)."""
    phi = np.zeros(2**N, dtype=complex)
    phi[0] = 1.0
    phi[2**(N - 1 - site)] = 1.0
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


def run_and_save(out_dir, label, J_list, init_fn, hamming, times):
    """Propagate and dump purity (and CPsi, MI, L1) to a .npz file."""
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

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / f"experiment_B_{label}.npz",
             purity=purity, cpsi_0k=cpsi_0k, mi_0k=mi_0k, l1=l1)


def write_control_A(out_dir, init_fn, hamming, times):
    """Write the uniform-J control (Experiment A) into the scan dir."""
    H = build_H_XY([J_UNIFORM] * (N - 1), N)
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
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "experiment_A.npz",
             purity=purity, cpsi_0k=cpsi_0k, mi_0k=mi_0k, l1=l1)
    np.save(out_dir / "times.npy", times)


# ---------------------------------------------------------------------------
# Scan driver
# ---------------------------------------------------------------------------
def run_scan_finer_Jmod(times, hamming):
    """Finer J_mod grid at default defect (0, 1), ψ_1 initial state."""
    out_dir = RESULTS_ROOT / "finer_Jmod"
    write_control_A(out_dir, initial_rho_psi1_bonding, hamming, times)
    J_list_grid = [0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0]
    print(f"[finer_Jmod] out = {out_dir}")
    for J_mod in J_list_grid:
        t0 = time.time()
        J_list = [J_mod] + [J_UNIFORM] * (N - 2)
        # label with explicit digits to avoid ambiguity
        label = f"{J_mod:g}"
        run_and_save(out_dir, label, J_list, initial_rho_psi1_bonding,
                     hamming, times)
        print(f"  J_mod = {J_mod:g}   {time.time() - t0:.2f} s")


def run_scan_defect_location(times, hamming, J_mods=(0.5, 2.0)):
    """Defect at various bonds: (1,2), (2,3), (3,4), (4,5), (5,6)."""
    # The bond index i corresponds to the pair (i, i+1). For N = 7 we have
    # bonds 0..5. We already have bond 0 in the main sweep; here we do 1..5.
    for bond in [1, 2, 3, 4, 5]:
        out_dir = RESULTS_ROOT / f"bond_{bond}_{bond + 1}"
        write_control_A(out_dir, initial_rho_psi1_bonding, hamming, times)
        print(f"[defect bond ({bond},{bond + 1})] out = {out_dir}")
        for J_mod in J_mods:
            t0 = time.time()
            J_list = [J_UNIFORM] * (N - 1)
            J_list[bond] = J_mod
            label = f"{J_mod:g}"
            run_and_save(out_dir, label, J_list, initial_rho_psi1_bonding,
                         hamming, times)
            print(f"  J_mod = {J_mod:g}   {time.time() - t0:.2f} s")


def run_scan_localized_init(times, hamming):
    """Localized initial state |1_3> at default defect (0, 1) with J_mod
    in {0.5, 1.5, 2.0}."""
    out_dir = RESULTS_ROOT / "local_site3"
    init_fn = lambda N: initial_rho_local_excitation(N, 3)
    write_control_A(out_dir, init_fn, hamming, times)
    print(f"[local_site3] out = {out_dir}")
    for J_mod in [0.5, 1.5, 2.0]:
        t0 = time.time()
        J_list = [J_mod] + [J_UNIFORM] * (N - 2)
        label = f"{J_mod:g}"
        run_and_save(out_dir, label, J_list, init_fn, hamming, times)
        print(f"  J_mod = {J_mod:g}   {time.time() - t0:.2f} s")


def main():
    times = np.arange(N_STEPS + 1) * DT
    hamming = build_hamming_matrix(N)
    print("N = 7 coupling-defect overlay, EXTENDED runs")
    print(f"  T = {T_FINAL}, steps = {N_STEPS}, gamma_0 = {GAMMA_0}")
    print(f"  output root: {RESULTS_ROOT}")
    t0 = time.time()
    run_scan_finer_Jmod(times, hamming)
    run_scan_defect_location(times, hamming, J_mods=(0.5, 2.0))
    run_scan_localized_init(times, hamming)
    print(f"\nTotal extended-scan runtime: {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
