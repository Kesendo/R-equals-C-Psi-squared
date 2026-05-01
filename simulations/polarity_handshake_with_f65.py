#!/usr/bin/env python3
"""Handshake test: bridge_dynamics' Δt vs F65's analytical eigenmode prediction.

For N=3, the |L, +, L⟩ initial state's middle-site polarity oscillates at
period Δt = π/E_1 where E_1 = 2J·cos(π/(N+1)) is F65's dominant single-
excitation mode energy.

This script tests whether the equivalence holds across N. If it does, the
two tool chains — analytical eigendecomp (F65) and trajectory + sign-change
detection (polarity_crossings) — agree on the same number computed via
different operations. That's the structural form of a handshake: agreement
through independent paths to the same value.

For N=3 the prediction is exact (only k=1 and k=N-1=3 modes contribute,
palindromic pair). For larger N the |L, ..., +, ..., L⟩ state mixes more
modes; if Δt still tracks π/E_1, the dominant single-mode prediction is
robust to multi-mode interference.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def f65_dominant_energy(N, J=1.0):
    """E_1 = 2J·cos(π/(N+1)) — the dominant single-excitation eigenmode."""
    return 2.0 * J * np.cos(np.pi / (N + 1))


def measure_dt_bonding_pair(N, t_max=20.0, n_t=4000, gamma_0=0.05, J=1.0, k=1):
    """Run bloch_trajectory with the bonding-mode pair state (|vac⟩+|ψ_k⟩)/√2.

    The off-diagonal |vac⟩⟨ψ_k|(t) oscillates at frequency E_k (single mode),
    so each site's ⟨X_i⟩ traces a clean cos(E_k t) curve with amplitude ψ_k(i).
    The +→−, −→+ crossing duration should be exactly π/E_k = Δt_F65.
    """
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=gamma_0, J=J)
    psi = fw.bonding_mode_pair_state(N, k)
    t_grid = np.linspace(0.0, t_max, n_t)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    events = fw.polarity_crossings(traj, t_grid)
    sig = fw.bridge_reflection_signature(events)
    if not sig['cycles']:
        return None, len(events)
    # Use the median cycle duration across sites (robust to per-site outliers
    # arising from amplitude-decay effects near the polarity boundary)
    durations = [c['cycle_duration'] for c in sig['cycles']]
    return float(np.median(durations)), len(events)


def main():
    print("Handshake test: F65 analytical Δt_k vs bloch_trajectory measured Δt")
    print("  Initial state: bonding_mode_pair_state(N, k=1) = (|vac⟩ + |ψ_1⟩)/√2")
    print("  Hamiltonian: XY chain, J=1, γ₀=0.05 (framework default)")
    print("  Predicted: Δt_F65 = π/E_k where E_k = 2J·cos(πk/(N+1))")
    print()
    print(f"  {'N':>3s}  {'E_1 (F65)':>10s}  {'Δt_F65':>10s}  "
          f"{'Δt_meas':>10s}  {'rel.err.':>10s}  {'#crossings':>10s}")
    print(f"  {'-' * 3}  {'-' * 10}  {'-' * 10}  "
          f"{'-' * 10}  {'-' * 10}  {'-' * 10}")

    rows = []
    for N in (3, 4, 5, 6):
        E_1 = f65_dominant_energy(N)
        dt_predicted = float(np.pi / E_1)
        print(f"  computing N={N}...", flush=True)
        dt_measured, n_cross = measure_dt_bonding_pair(N, t_max=20.0, n_t=4000, k=1)
        if dt_measured is None:
            print(f"  {N:>3d}  {E_1:>10.4f}  {dt_predicted:>10.4f}  "
                  f"{'(no cycle)':>10s}  {'—':>10s}  {n_cross:>10d}")
            continue
        rel_err = abs(dt_measured - dt_predicted) / dt_predicted
        rows.append((N, E_1, dt_predicted, dt_measured, rel_err))
        print(f"  {N:>3d}  {E_1:>10.4f}  {dt_predicted:>10.4f}  "
              f"{dt_measured:>10.4f}  {rel_err:>10.2e}  {n_cross:>10d}")

    print()
    if rows and all(r[4] < 0.01 for r in rows):
        print("  ✓ Handshake holds across all N: the two tool chains agree to <1%.")
    elif rows:
        max_err = max(r[4] for r in rows)
        print(f"  Max rel.err: {max_err:.2e}")
        if max_err > 0.1:
            print("  ⚠ Larger deviation than expected for single-mode initial state;")
            print("  damping effects or other modes still contributing.")
        else:
            print("  Single-mode F65 prediction tracks measured Δt within a few %.")


if __name__ == '__main__':
    main()
