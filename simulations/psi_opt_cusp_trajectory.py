"""
Phase B: Cusp trajectory for psi_opt, Bell+, W5 under sacrifice Lindblad.
Convention: ibm_april_predictions (np.kron(H, Id) for commutator).

CΨ = Tr(rho²) * L1 / (d-1)  computed on 2-qubit reduced density matrices.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 10, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import numpy as np
from scipy import linalg

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H, build_liouvillian, bell_plus_center, w_full,
    partial_trace_to_pair, wootters_concurrence,
)
from optimal_state_n5_sacrifice import (
    PROFILES, N, gamma_sacrifice, Sg_sac,
)

D = 2**N
SE_IDX = [1 << (N - 1 - k) for k in range(N)]

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "psi_opt_cusp")
os.makedirs(out_dir, exist_ok=True)


# ===================================================================
# CΨ computation
# ===================================================================
def cpsi_pair(rho_2q):
    """CΨ = Tr(rho²) * L1 / (d-1) for a 4x4 density matrix (d=4)."""
    purity = float(np.real(np.trace(rho_2q @ rho_2q)))
    l1 = float(np.sum(np.abs(rho_2q)) - np.sum(np.abs(np.diag(rho_2q))))
    return purity * l1 / 3.0


def cpsi_all_pairs(rho_full, N):
    """CΨ for each adjacent pair. Returns array of N-1 values."""
    vals = []
    for k in range(N - 1):
        rho_2q = partial_trace_to_pair(rho_full, N, k, k+1)
        vals.append(cpsi_pair(rho_2q))
    return np.array(vals)


# ===================================================================
# Time evolution
# ===================================================================
def evolve_trajectory(psi0, profile, times):
    """Returns dict with per-pair CΨ(t), mean, max, also sanity checks."""
    rho0 = np.outer(psi0, psi0.conj())
    c0 = profile['R_inv'] @ rho0.flatten()

    n_pairs = N - 1
    cpsi_pairs = np.zeros((len(times), n_pairs))
    trace_check = np.zeros(len(times))
    min_eig_check = np.zeros(len(times))

    for i, t in enumerate(times):
        rho_vec = profile['R'] @ (c0 * np.exp(profile['eigvals'] * t))
        rho_t = rho_vec.reshape(D, D)
        rho_t = (rho_t + rho_t.conj().T) / 2

        trace_check[i] = float(np.real(np.trace(rho_t)))
        eigs_rho = np.linalg.eigvalsh(rho_t)
        min_eig_check[i] = float(eigs_rho.min())

        cpsi_pairs[i] = cpsi_all_pairs(rho_t, N)

    return {
        'cpsi_pairs': cpsi_pairs,
        'cpsi_mean': cpsi_pairs.mean(axis=1),
        'cpsi_max': cpsi_pairs.max(axis=1),
        'trace': trace_check,
        'min_eig': min_eig_check,
    }


def find_crossings(t, cpsi, level=0.25):
    """Find time points where cpsi crosses level (downward)."""
    crossings = []
    for i in range(len(cpsi) - 1):
        if cpsi[i] >= level and cpsi[i+1] < level:
            # Linear interpolation
            frac = (cpsi[i] - level) / (cpsi[i] - cpsi[i+1])
            t_cross = t[i] + frac * (t[i+1] - t[i])
            # Finite difference for dCPsi/dt at crossing
            dt = t[1] - t[0]
            if i > 0 and i < len(cpsi) - 1:
                dcpsi_dt = (cpsi[i+1] - cpsi[i-1]) / (2 * dt)
            else:
                dcpsi_dt = (cpsi[i+1] - cpsi[i]) / dt
            crossings.append({
                't_cross': float(t_cross),
                'dcpsi_dt': float(dcpsi_dt),
                'dwell_prefactor': float(2.0 / abs(dcpsi_dt)) if abs(dcpsi_dt) > 1e-15 else float('inf'),
            })
    return crossings


# ===================================================================
# State constructors
# ===================================================================
def pure_from_rho(rho):
    w, v = linalg.eigh(rho)
    k = int(np.argmax(w.real))
    return v[:, k] * np.exp(-1j * np.angle(v[0, k] if abs(v[0, k]) > 1e-10 else 1.0))

def make_psi_opt():
    a = np.array([0.099342, 0.238952, 0.427987, 0.571584, 0.650501])
    psi = np.zeros(D, dtype=complex)
    for k in range(N):
        psi[SE_IDX[k]] = a[k]
    return psi / np.linalg.norm(psi)


# ===================================================================
if __name__ == "__main__":
    import time as _time
    t_start = _time.time()

    prof = PROFILES['sacrifice']
    times = np.linspace(0, 20, 2001)  # finer grid for crossing detection

    print("=" * 80)
    print("PHASE B: Cusp trajectories under sacrifice profile")
    print("=" * 80)
    print(f"  CΨ = Tr(rho_pair²) * L1_pair / 3  (pairwise, d=4)")
    print(f"  Pairs: (0,1), (1,2), (2,3), (3,4)")
    print(f"  T_max = {times[-1]}, dt = {times[1]-times[0]:.4f}")

    # Build states
    states = {
        'psi_opt':  make_psi_opt(),
        'Bell+':    pure_from_rho(bell_plus_center(N)),
        'W5_full':  pure_from_rho(w_full(N)),
    }

    # Also add sacrifice_tuned_W5 and W4_1234 as extras
    a_heur = np.sqrt(gamma_sacrifice.min() / gamma_sacrifice)
    psi_heur = np.zeros(D, dtype=complex)
    for k in range(N):
        psi_heur[SE_IDX[k]] = a_heur[k]
    psi_heur /= np.linalg.norm(psi_heur)
    states['sac_tuned_W5'] = psi_heur

    results = {}
    for name, psi in states.items():
        print(f"\n  --- {name} ---")
        traj = evolve_trajectory(psi, prof, times)
        results[name] = traj

        # Sanity checks
        tr_err = np.max(np.abs(traj['trace'] - 1.0))
        min_eig = traj['min_eig'].min()
        print(f"    Trace conservation: max|Tr-1| = {tr_err:.2e}")
        print(f"    Positivity: min eigenvalue = {min_eig:.2e}")

        # Initial CΨ values
        cpsi0 = traj['cpsi_pairs'][0]
        print(f"    CΨ(t=0) per pair: [{', '.join(f'{v:.4f}' for v in cpsi0)}]")
        print(f"    CΨ_mean(t=0) = {traj['cpsi_mean'][0]:.4f}")
        print(f"    CΨ_max(t=0) = {traj['cpsi_max'][0]:.4f}")

        if traj['cpsi_mean'][0] > 1/3:
            print(f"    ** CΨ_mean(0) > 1/3! **")
        if traj['cpsi_max'][0] > 1/3:
            print(f"    ** CΨ_max(0) > 1/3 on at least one pair! **")

        # Crossings of 1/4 for mean CΨ
        crossings_mean = find_crossings(times, traj['cpsi_mean'], 0.25)
        if crossings_mean:
            cr = crossings_mean[0]
            print(f"    CΨ_mean crosses 1/4 at t = {cr['t_cross']:.4f}")
            print(f"    |dCΨ/dt| at crossing = {abs(cr['dcpsi_dt']):.6f}")
            print(f"    Dwell prefactor = {cr['dwell_prefactor']:.4f}")
        else:
            peak = traj['cpsi_mean'].max()
            print(f"    CΨ_mean does NOT cross 1/4 (max = {peak:.4f})")

        # Crossings per pair
        for p in range(N - 1):
            cx = find_crossings(times, traj['cpsi_pairs'][:, p], 0.25)
            if cx:
                print(f"    Pair ({p},{p+1}): crosses 1/4 at t={cx[0]['t_cross']:.4f},"
                      f" |dCΨ/dt|={abs(cx[0]['dcpsi_dt']):.6f}")

    # ==================================================================
    # Comparison table
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("CROSSING COMPARISON")
    print("=" * 80)

    # Bell+ reference speed
    bell_crossings = find_crossings(times, results['Bell+']['cpsi_mean'], 0.25)
    bell_speed = abs(bell_crossings[0]['dcpsi_dt']) if bell_crossings else None

    print(f"\n  {'State':<18}  {'CΨ(0)':>6}  {'crosses?':>8}  {'t_cross':>8}  "
          f"{'|dCΨ/dt|':>10}  {'ratio vs Bell+':>14}")
    print(f"  {'-' * 70}")

    crossing_data = {}
    for name in ['psi_opt', 'Bell+', 'W5_full', 'sac_tuned_W5']:
        traj = results[name]
        cx = find_crossings(times, traj['cpsi_mean'], 0.25)
        cpsi0 = traj['cpsi_mean'][0]
        if cx:
            speed = abs(cx[0]['dcpsi_dt'])
            ratio = speed / bell_speed if bell_speed else float('nan')
            print(f"  {name:<18}  {cpsi0:>6.4f}  {'YES':>8}  {cx[0]['t_cross']:>8.4f}  "
                  f"{speed:>10.6f}  {ratio:>14.4f}")
            crossing_data[name] = {
                't_cross': cx[0]['t_cross'],
                'dcpsi_dt': cx[0]['dcpsi_dt'],
                'dwell_prefactor': cx[0]['dwell_prefactor'],
                'speed_ratio_vs_bell': ratio,
            }
        else:
            peak = traj['cpsi_mean'].max()
            print(f"  {name:<18}  {cpsi0:>6.4f}  {'NO':>8}  {'':>8}  "
                  f"{'':>10}  (max CΨ={peak:.4f})")
            crossing_data[name] = {'crosses': False, 'max_cpsi': float(peak)}

    # ==================================================================
    # Plots
    # ==================================================================
    print(f"\n  Generating plot ...")
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'psi_opt': '#d62728', 'Bell+': '#1f77b4', 'W5_full': '#2ca02c',
              'sac_tuned_W5': '#9467bd'}
    for name in ['psi_opt', 'Bell+', 'W5_full', 'sac_tuned_W5']:
        traj = results[name]
        ax.plot(times, traj['cpsi_mean'], label=name, linewidth=2, color=colors[name])
    ax.axhline(0.25, color='red', linestyle=':', alpha=0.7, label='CΨ = 1/4 (cusp)')
    ax.axhline(1/3, color='orange', linestyle=':', alpha=0.5, label='CΨ = 1/3')
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('CΨ (mean adjacent)')
    ax.set_title('CΨ trajectories under sacrifice profile')
    ax.set_xlim(0, 20)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'trajectories.png'), dpi=150)
    plt.close(fig)
    print("  Plot saved")

    # ==================================================================
    # Save JSON
    # ==================================================================
    traj_json = {}
    for name in states:
        traj = results[name]
        traj_json[name] = {
            't': times.tolist(),
            'cpsi_mean': traj['cpsi_mean'].tolist(),
            'cpsi_max': traj['cpsi_max'].tolist(),
            'cpsi_pairs': traj['cpsi_pairs'].tolist(),
        }
    with open(os.path.join(out_dir, 'trajectories.json'), 'w') as f:
        json.dump(traj_json, f)

    with open(os.path.join(out_dir, 'crossings.json'), 'w') as f:
        json.dump(crossing_data, f, indent=2)

    print(f"\n  Saved trajectories.json and crossings.json")
    print(f"\nTotal runtime: {_time.time() - t_start:.1f}s")
