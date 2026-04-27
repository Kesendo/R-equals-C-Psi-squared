#!/usr/bin/env python3
"""Snapshot E: framework.py-grounded exact hardware predictions for chiral mirror.

Companion to Snapshot D (Marrakesh hardware finale, April 26):
  Snapshot D: framework predicted ⟨XIZ⟩ = -0.62 for soft case at N=3, hardware
              measured -0.71 to -0.92 across three Heron r2 backends.

Snapshot E (proposed): framework.py Section 16 (added today) gives exact
  predictions for per-site Bloch components and pair-purity trajectories on
  K_1-paired sine-mode bonding states at a target N-qubit chain. The chiral
  mirror identity (today's structural result) gives concrete
  cross-state relations testable on hardware.

Targeted backend parameters:
  ibm_marrakesh path [48, 49, 50, 51] (4-qubit linear sub-path, Heron r2)
  J ≈ 1.84 MHz (native CZ), γ ≈ 0.005 /μs (median dephasing)
  t = {0.5, 1.0, 2.0, 5.0} μs (4 timepoints)
  Initial states: (|vac⟩ + |ψ_2⟩)/√2 and (|vac⟩ + |ψ_3⟩)/√2 (K_1-paired at N=4)
  XX+YY hopping only (preserves K_1 symmetry on OBC chain)

For each timepoint:
  - ⟨X_i⟩(t), ⟨Y_i⟩(t), ⟨Z_i⟩(t), P_i(t)  for i = 0..3
  - K_1-mirror predicted relations:
      ⟨X_i⟩(ψ_3) = +(−1)^i · ⟨X_i⟩(ψ_2)
      ⟨Y_i⟩(ψ_3) = −(−1)^i · ⟨Y_i⟩(ψ_2)
      ⟨Z_i⟩(ψ_3) = +⟨Z_i⟩(ψ_2)
      P_i(ψ_3) = P_i(ψ_2)
  - Pair purity P_{0,3}(t) (endpoint pair)

Output: concrete numeric prediction table, hardware-ready.
"""
from __future__ import annotations

import math
import sys
from itertools import combinations
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

# Use today's framework.py Section 16 primitives directly
import framework as fw
from _eq014_chiral_mirror_multi_N import (
    build_xy_chain_H, build_L, sine_mode_state, rk4_step,
)


# Marrakesh-realistic parameters
N = 4
J = 1.84      # MHz, native CZ angle pi/2 / (2 * gate_time)
GAMMA = 0.005   # /μs, median dephasing across path [48..51]
TIMES_US = [0.5, 1.0, 2.0, 5.0]   # μs
BOND = (0, 1)


def site_op_dense(op, i, N):
    full = np.array([[1.0]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    for j in range(N):
        full = np.kron(full, op if j == i else I2)
    return full


def bonding_density(N, k):
    d = 2 ** N
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = sine_mode_state(N, k, d)
    phi = (vac + psi) / np.sqrt(2)
    return np.outer(phi, phi.conj())


def evolve(L, rho_flat, t_us, dt_us=0.05):
    n_steps = max(1, int(np.ceil(t_us / dt_us)))
    dt = t_us / n_steps
    cur = rho_flat.astype(complex).copy()
    for _ in range(n_steps):
        cur = rk4_step(L, cur, dt)
    return cur


def per_site_bloch(rho_full, N):
    """Per-site (⟨X⟩, ⟨Y⟩, ⟨Z⟩, P_i) using framework's k_local_reduced_density."""
    out = []
    for i in range(N):
        rho_i = fw.k_local_reduced_density(rho_full, [i], N)
        X_op = np.array([[0, 1], [1, 0]], dtype=complex)
        Y_op = np.array([[0, -1j], [1j, 0]], dtype=complex)
        Z_op = np.array([[1, 0], [0, -1]], dtype=complex)
        x = float(np.real(np.trace(rho_i @ X_op)))
        y = float(np.real(np.trace(rho_i @ Y_op)))
        z = float(np.real(np.trace(rho_i @ Z_op)))
        P = float(np.real(np.trace(rho_i @ rho_i)))
        out.append((x, y, z, P))
    return out


def pair_purity(rho_full, sites, N):
    return fw.k_local_purity(rho_full, sites, N)


def main():
    print("=" * 80)
    print(f"Snapshot E framework predictions: chiral mirror at N={N} chain")
    print("=" * 80)
    print(f"  Backend target: ibm_marrakesh path [48, 49, 50, 51] (4-qubit OBC)")
    print(f"  J = {J} MHz native, γ = {GAMMA} /μs median dephasing")
    print(f"  H = J·Σ_<ij> (X_i X_j + Y_i Y_j)/2  (XX+YY hopping, preserves K_1)")
    print(f"  States: (|vac⟩ + |ψ_k⟩)/√2 for k ∈ {{2, 3}} (K_1-paired at N=4)")
    print()

    # Build operators in MHz / μs units (J=1.84 MHz means 1.84 cycles/μs in
    # the sense that this script uses; treated as "per μs" frequency throughout)
    H = build_xy_chain_H(N, J)
    L = build_L(H, N, GAMMA)

    # Initial states
    rho2_0 = bonding_density(N, 2).flatten()
    rho3_0 = bonding_density(N, 3).flatten()

    # Verify K_1 mirror identity at t=0 using framework primitive
    psi_2 = fw.single_excitation_sine_mode(N, 2)
    psi_3 = fw.single_excitation_sine_mode(N, 3)
    K_1_diag = np.array([(-1) ** i for i in range(N)])
    K1_psi2 = K_1_diag * psi_2
    print(f"  Sanity: ‖K_1 ψ_2 − ψ_3‖ = {np.linalg.norm(K1_psi2 - psi_3):.2e}  (must be ~0)")
    print()

    print(f"  {'t (μs)':<8} {'state':<8} ", end="")
    for i in range(N):
        print(f"{'site '+str(i):>10}    ", end="")
    print(f"{'P_{0,N-1}':>10}")
    print(f"  {'':>8} {'':>8} ", end="")
    for i in range(N):
        print(f"{'⟨X⟩,⟨Y⟩,⟨Z⟩':>14}", end="")
    print()
    print("  " + "-" * 80)

    chiral_predictions = []
    for t in TIMES_US:
        rho2_t = evolve(L, rho2_0, t)
        rho3_t = evolve(L, rho3_0, t)
        rho2_full = rho2_t.reshape(2 ** N, 2 ** N)
        rho3_full = rho3_t.reshape(2 ** N, 2 ** N)

        bloch2 = per_site_bloch(rho2_full, N)
        bloch3 = per_site_bloch(rho3_full, N)
        P_pair2 = pair_purity(rho2_full, [0, N - 1], N)
        P_pair3 = pair_purity(rho3_full, [0, N - 1], N)

        # ψ_2 row
        line = f"  {t:<8.2f} {'ψ_2':<8} "
        for (x, y, z, P) in bloch2:
            line += f"{x:>+5.2f},{y:>+5.2f},{z:>+5.2f}"
        line += f"  {P_pair2:>10.4f}"
        print(line)

        # ψ_3 row
        line = f"  {'':<8} {'ψ_3':<8} "
        for (x, y, z, P) in bloch3:
            line += f"{x:>+5.2f},{y:>+5.2f},{z:>+5.2f}"
        line += f"  {P_pair3:>10.4f}"
        print(line)

        # Chiral mirror prediction check
        signs = [(-1) ** i for i in range(N)]
        max_X_diff = max(abs(bloch3[i][0] - signs[i] * bloch2[i][0]) for i in range(N))
        max_Y_diff = max(abs(bloch3[i][1] + signs[i] * bloch2[i][1]) for i in range(N))
        max_Z_diff = max(abs(bloch3[i][2] - bloch2[i][2]) for i in range(N))
        max_P_pair_diff = abs(P_pair3 - P_pair2)
        line = f"  {'':<8} {'mirror':<8} "
        line += f"  X-flip ‖Δ‖={max_X_diff:.1e}, Y-flip ‖Δ‖={max_Y_diff:.1e}, "
        line += f"Z ‖Δ‖={max_Z_diff:.1e}, P_pair ‖Δ‖={max_P_pair_diff:.1e}"
        print(line)
        print()

        chiral_predictions.append({
            "t_us": t, "ψ_2_bloch": bloch2, "ψ_3_bloch": bloch3,
            "P_pair_2": P_pair2, "P_pair_3": P_pair3,
            "max_X_diff": max_X_diff, "max_Y_diff": max_Y_diff,
            "max_Z_diff": max_Z_diff, "max_P_pair_diff": max_P_pair_diff,
        })

    # Concrete numeric prediction table for hardware comparison
    print()
    print("=" * 80)
    print("Hardware-ready prediction table (concrete numbers):")
    print("=" * 80)
    print()
    for cp in chiral_predictions:
        t = cp["t_us"]
        print(f"  t = {t:.2f} μs:")
        print(f"    Initial state (|vac⟩ + |ψ_2⟩)/√2 [N={N}, J={J} MHz, γ={GAMMA} /μs]:")
        for i, (x, y, z, P) in enumerate(cp["ψ_2_bloch"]):
            print(f"      site {i}: ⟨X⟩ = {x:+.5f}, ⟨Y⟩ = {y:+.5f}, "
                  f"⟨Z⟩ = {z:+.5f}, P_i = {P:.5f}")
        print(f"    Endpoint pair purity P_{{0,{N-1}}} = {cp['P_pair_2']:.5f}")
        print(f"    Initial state (|vac⟩ + |ψ_3⟩)/√2:")
        for i, (x, y, z, P) in enumerate(cp["ψ_3_bloch"]):
            print(f"      site {i}: ⟨X⟩ = {x:+.5f}, ⟨Y⟩ = {y:+.5f}, "
                  f"⟨Z⟩ = {z:+.5f}, P_i = {P:.5f}")
        print(f"    Endpoint pair purity P_{{0,{N-1}}} = {cp['P_pair_3']:.5f}")
        print(f"    K_1-mirror identity (framework prediction):")
        print(f"      max |⟨X_i⟩(ψ_3) − (−1)^i ⟨X_i⟩(ψ_2)| = {cp['max_X_diff']:.2e}  (predict 0)")
        print(f"      max |⟨Y_i⟩(ψ_3) + (−1)^i ⟨Y_i⟩(ψ_2)| = {cp['max_Y_diff']:.2e}  (predict 0)")
        print(f"      max |⟨Z_i⟩(ψ_3) − ⟨Z_i⟩(ψ_2)|       = {cp['max_Z_diff']:.2e}  (predict 0)")
        print(f"      max |P_pair(ψ_3) − P_pair(ψ_2)|     = {cp['max_P_pair_diff']:.2e}  (predict 0)")
        print()


if __name__ == "__main__":
    main()
