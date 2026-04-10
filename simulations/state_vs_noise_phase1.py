"""
Phase 1: Time-Evolution and Slow-Band Diagnosis
================================================
Compares 6 initial states under sacrifice and uniform profiles.
Part A: concurrence(t) curves, AUC tables, plots
Part B: W5 slow-band eigenmode diagnosis

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 9, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import time as _time
import numpy as np
from scipy import linalg

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H, build_liouvillian, bell_plus_center, w_full, w_center3_plus,
    partial_trace_to_pair, wootters_concurrence, pauli_basis, n_xy_string,
)
from optimal_state_n5_sacrifice import (
    PROFILES, N, gamma_sacrifice, gamma_uniform, Sg_sac, Sg_uni,
    analyze_state, max_adjacent_concurrence,
)


# ===================================================================
# State constructors
# ===================================================================
def pure_from_rho(rho):
    """Extract pure state from rank-1 density matrix."""
    w, v = linalg.eigh(rho)
    k = int(np.argmax(w.real))
    return v[:, k] * np.exp(-1j * np.angle(v[0, k] if abs(v[0, k]) > 1e-10 else 1.0))


def make_opt_le1():
    """OPT(le1) from cached amplitudes (results file)."""
    psi = np.zeros(2**N, dtype=complex)
    psi[1]  = +0.54356969 + 0.48571456j
    psi[0]  = +0.38343498 + 0.00000000j
    psi[4]  = +0.35124362 - 0.04751746j
    psi[2]  = -0.25122721 - 0.20219736j
    psi[8]  = -0.21216475 + 0.18459728j
    psi[16] = +0.10792562 - 0.03498738j
    psi /= np.linalg.norm(psi)
    return psi


def make_opt_le2():
    """OPT(le2) from cached amplitudes (results file)."""
    psi = np.zeros(2**N, dtype=complex)
    psi[3]  = +0.50697710 - 0.50704621j
    psi[1]  = -0.29531726 - 0.34785300j
    psi[6]  = -0.11072701 + 0.37682604j
    psi[0]  = -0.16649623 + 0.00000000j
    psi[2]  = -0.12074271 - 0.09930390j
    psi[5]  = +0.09193044 + 0.10389792j
    psi[8]  = +0.12757914 - 0.00776383j
    psi[12] = +0.09544184 - 0.03525368j
    psi[10] = -0.05358119 + 0.08543979j
    psi[4]  = -0.02441640 + 0.07760813j
    psi /= np.linalg.norm(psi)
    return psi


def make_plus5():
    """|+>^5 product state."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus.copy()
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    return psi


# ===================================================================
# Time evolution via eigendecomposition
# ===================================================================
def safe_concurrence(rho_ij):
    """Wootters concurrence with NaN -> 0 safety for singular matrices."""
    c = wootters_concurrence(rho_ij)
    return 0.0 if np.isnan(c) else float(c)


def compute_concurrence_curves(psi, profile, times):
    """Concurrence(t) for all adjacent pairs. Returns (C_pairs, C_mean, C_max)."""
    rho0 = np.outer(psi, psi.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    d = 2**N
    n_pairs = N - 1

    C_pairs = np.zeros((len(times), n_pairs))
    for ti, t in enumerate(times):
        rho_vec = profile['R'] @ (c0 * np.exp(profile['eigvals'] * t))
        rho_t = rho_vec.reshape(d, d)
        rho_t = (rho_t + rho_t.conj().T) / 2
        for pi in range(n_pairs):
            rho_ij = partial_trace_to_pair(rho_t, N, pi, pi + 1)
            C_pairs[ti, pi] = safe_concurrence(rho_ij)

    C_mean = C_pairs.mean(axis=1)
    C_max = C_pairs.max(axis=1)
    return C_pairs, C_mean, C_max


# ===================================================================
# Part B: W5 slow-band diagnosis
# ===================================================================
def w5_slow_band_diagnosis(profile, top_k=20):
    """Project W5 onto slow-band eigenmodes, return diagnostic rows."""
    psi_w5 = pure_from_rho(w_full(N))
    rho0 = np.outer(psi_w5, psi_w5.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    probs = np.abs(c0) ** 2
    total = probs.sum()

    strings, P_mats = pauli_basis(N)
    P_flat = np.array([P.flatten() for P in P_mats])

    slow_mask = profile['slow_mask']
    slow_indices = np.where(slow_mask)[0]
    slow_probs = probs[slow_indices]

    order = np.argsort(-slow_probs)[:top_k]

    rows = []
    for o in order:
        k = slow_indices[o]
        ev = profile['eigvals'][k]
        proj = probs[k]

        v = profile['R'][:, k]
        coeffs = P_flat.conj() @ v / (2**N)
        probs_p = np.abs(coeffs)**2
        p_sum = probs_p.sum()
        if p_sum > 1e-30:
            w_norm = probs_p / p_sum
            n_xy = float(sum(wi * n_xy_string(s) for wi, s in zip(w_norm, strings)))
        else:
            n_xy = 0.0

        rows.append({
            'k': int(k),
            're_lambda': float(ev.real),
            'im_lambda': float(ev.imag),
            'n_xy': n_xy,
            'proj': float(proj),
            'norm_pct': float(proj / total * 100),
        })

    return rows, float(total)


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    t_start = _time.time()

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "results", "state_vs_noise_phase1")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("PHASE 1: Time-Evolution and Slow-Band Diagnosis")
    print("=" * 80)
    print(f"  N = {N}")
    print(f"  gamma_sacrifice = {np.round(gamma_sacrifice, 5)}")
    print(f"  Sg = {Sg_sac:.5f}")

    # --- Construct all 6 states ---
    states = [
        ('Bell_plus',  'Bell+|+>',     pure_from_rho(bell_plus_center(N))),
        ('W3_center',  'W3 center|+>', pure_from_rho(w_center3_plus(N))),
        ('W5_full',    'W5 full',      pure_from_rho(w_full(N))),
        ('OPT_le1',   'OPT(le1)',     make_opt_le1()),
        ('OPT_le2',   'OPT(le2)',     make_opt_le2()),
        ('plus5',      '|+>^5',        make_plus5()),
    ]

    # Verify OPT states reproduce known slow_wt
    for skey, sname, psi in states:
        if skey.startswith('OPT'):
            w = analyze_state(psi, PROFILES['sacrifice'])
            print(f"  {sname} slow_wt check: {w['slow_wt']:.2f}%")

    # Verify propagator: R @ c0 should reconstruct rho0.flatten()
    print("\n  Propagator verification:")
    max_dev_all = 0.0
    for plabel, profile in PROFILES.items():
        for skey, sname, psi in states:
            rho0 = np.outer(psi, psi.conj())
            c0 = profile['R_inv'] @ rho0.flatten()
            recon = profile['R'] @ c0
            dev = np.max(np.abs(recon - rho0.flatten()))
            max_dev_all = max(max_dev_all, dev)
            if dev > 1e-8:
                print(f"    WARNING: {sname}/{plabel} max_dev = {dev:.2e}")
    print(f"    All passed: max_dev = {max_dev_all:.2e}")

    # Sanity: |+>^5 concurrence at t=0
    psi_plus5 = states[5][2]
    rho_plus5 = np.outer(psi_plus5, psi_plus5.conj())
    concs_plus5 = []
    for i in range(N - 1):
        rho_ij = partial_trace_to_pair(rho_plus5, N, i, i + 1)
        concs_plus5.append(safe_concurrence(rho_ij))
    print(f"\n  |+>^5 sanity: C(t=0) = {[f'{c:.6f}' for c in concs_plus5]}")
    assert all(c < 1e-6 for c in concs_plus5), f"|+>^5 nonzero concurrence: {concs_plus5}"
    print("  PASSED: all zero")

    # ==================================================================
    # Part A: Time evolution
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("PART A: Time evolution of concurrence")
    print("=" * 80)

    times = np.linspace(0, 30, 301)
    profiles_list = [('sacrifice', PROFILES['sacrifice']),
                     ('uniform',   PROFILES['uniform'])]

    npz_data = {'t': times}
    auc_data = {}

    for skey, sname, psi in states:
        for plabel, profile in profiles_list:
            t0 = _time.time()
            C_pairs, C_mean, C_max = compute_concurrence_curves(psi, profile, times)
            dt = _time.time() - t0

            for pi in range(N - 1):
                npz_data[f'{skey}_{plabel}_C{pi}{pi+1}'] = C_pairs[:, pi]
            npz_data[f'{skey}_{plabel}_Cmean'] = C_mean
            npz_data[f'{skey}_{plabel}_Cmax'] = C_max

            mask2 = times <= 2.0
            mask10 = times <= 10.0
            auc = {
                'T2':  float(np.trapezoid(C_mean[mask2],  times[mask2])),
                'T10': float(np.trapezoid(C_mean[mask10], times[mask10])),
                'T30': float(np.trapezoid(C_mean,         times)),
            }
            auc_data[(skey, plabel)] = auc

            print(f"  {sname:<15} {plabel:<10}  "
                  f"AUC(2)={auc['T2']:.4f}  AUC(10)={auc['T10']:.4f}  "
                  f"AUC(30)={auc['T30']:.4f}  ({dt:.1f}s)")

    npz_path = os.path.join(out_dir, "concurrence_curves.npz")
    np.savez_compressed(npz_path, **npz_data)
    print(f"\n  Saved {npz_path}")

    # --- AUC Tables (formatted) ---
    print(f"\n  {'=' * 65}")
    for plabel in ['sacrifice', 'uniform']:
        print(f"\n  AUC TABLE: {plabel.upper()}")
        print(f"  {'State':<15}  {'AUC(T=2)':>10}  {'AUC(T=10)':>10}  {'AUC(T=30)':>10}")
        print(f"  {'-' * 52}")
        for skey, sname, _ in states:
            d = auc_data[(skey, plabel)]
            print(f"  {sname:<15}  {d['T2']:>10.4f}  {d['T10']:>10.4f}  {d['T30']:>10.4f}")

    # --- Plots ---
    print(f"\n{'=' * 80}")
    print("Generating plots")
    print("=" * 80)

    colors = {
        'Bell_plus': '#1f77b4', 'W3_center': '#ff7f0e', 'W5_full': '#2ca02c',
        'OPT_le1': '#d62728', 'OPT_le2': '#9467bd', 'plus5': '#8c564b',
    }

    # Individual state plots (sacrifice solid, uniform dashed)
    for skey, sname, _ in states:
        fig, ax = plt.subplots(figsize=(8, 5))
        c_sac = npz_data[f'{skey}_sacrifice_Cmean']
        c_uni = npz_data[f'{skey}_uniform_Cmean']
        ax.plot(times, c_sac, label='sacrifice', linewidth=2, color=colors[skey])
        ax.plot(times, c_uni, label='uniform', linewidth=2, linestyle='--',
                color=colors[skey], alpha=0.7)
        ax.set_xlabel('t (dimensionless)')
        ax.set_ylabel('mean adjacent concurrence')
        ax.set_title(sname)
        ax.set_xlim(0, 30)
        ax.set_ylim(bottom=0)
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f'{skey}.png'), dpi=150)
        plt.close(fig)

    # Summary plot: linear scale
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    for ax, (plabel, ptitle) in zip(axes, [('sacrifice', 'Sacrifice profile'),
                                            ('uniform', 'Uniform profile')]):
        for skey, sname, _ in states:
            c = npz_data[f'{skey}_{plabel}_Cmean']
            ax.plot(times, c, label=sname, linewidth=1.5, color=colors[skey])
        ax.set_xlabel('t (dimensionless)')
        ax.set_ylabel('mean adjacent concurrence')
        ax.set_title(ptitle)
        ax.set_xlim(0, 30)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'all_states_Cmean.png'), dpi=150)
    plt.close(fig)

    # Summary plot: log scale
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    for ax, (plabel, ptitle) in zip(axes, [('sacrifice', 'Sacrifice (log)'),
                                            ('uniform', 'Uniform (log)')]):
        for skey, sname, _ in states:
            c = npz_data[f'{skey}_{plabel}_Cmean']
            c_clip = np.maximum(c, 1e-10)
            ax.semilogy(times, c_clip, label=sname, linewidth=1.5, color=colors[skey])
        ax.set_xlabel('t (dimensionless)')
        ax.set_ylabel('mean adjacent concurrence')
        ax.set_title(ptitle)
        ax.set_xlim(0, 30)
        ax.set_ylim(1e-6, 1)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'all_states_Cmean_log.png'), dpi=150)
    plt.close(fig)
    print("  All plots saved")

    # ==================================================================
    # Part B: W5 slow-band diagnosis
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("PART B: W5 slow-band diagnosis")
    print("=" * 80)

    diag_path = os.path.join(out_dir, "w5_slow_band_diagnosis.txt")
    with open(diag_path, 'w', encoding='utf-8') as f:
        f.write("W5 slow-band eigenmode diagnosis\n")
        f.write(f"Computed: {np.datetime64('today')}\n\n")

        for plabel, profile in profiles_list:
            rows, total = w5_slow_band_diagnosis(profile, top_k=20)
            Sg = profile['Sg']
            slow_wt = analyze_state(pure_from_rho(w_full(N)), profile)['slow_wt']

            title = (f"{plabel.upper()} (Sg = {Sg:.5f}, "
                     f"slow band: 0 < rate < {Sg:.3f}, "
                     f"W5 slow_wt = {slow_wt:.2f}%)")
            hdr = (f"  {'#':>3}  {'k':>5}  {'Re(lambda)':>12}  {'Im(lambda)':>12}"
                   f"  {'<n_XY>':>8}  {'|c_k|^2':>12}  {'norm%':>8}")
            sep = f"  {'-' * 70}"

            print(f"\n  {title}")
            print(f"  Total overlap: {total:.6f}")
            print(hdr)
            print(sep)

            f.write(f"{'=' * 75}\n")
            f.write(f"{title}\n")
            f.write(f"{'=' * 75}\n")
            f.write(f"  Total overlap: {total:.6f}\n")
            f.write(hdr + "\n")
            f.write(sep + "\n")

            for i, r in enumerate(rows):
                line = (f"  {i+1:>3}  {r['k']:>5}  {r['re_lambda']:>+12.6f}"
                        f"  {r['im_lambda']:>+12.6f}  {r['n_xy']:>8.3f}"
                        f"  {r['proj']:>12.6f}  {r['norm_pct']:>7.3f}%")
                print(line)
                f.write(line + "\n")
            f.write("\n")

    print(f"\n  Saved {diag_path}")

    t_end = _time.time()
    print(f"\nTotal runtime: {t_end - t_start:.1f}s")
    print("Done.")
