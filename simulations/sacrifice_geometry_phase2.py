"""
Sacrifice Geometry Phase 2: Can We Beat W5?
============================================
Experiment A: systematic scan of ~20 structured candidate states
Experiment B (optional): surrogate-loss optimization on promising subspace

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 9, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import time as _time
import numpy as np
from scipy import linalg
from scipy.optimize import differential_evolution
from itertools import combinations

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

D = 2**N


# ===================================================================
# Utilities
# ===================================================================
def safe_wootters(rho):
    """Wootters concurrence with NaN guard for singular density matrices."""
    c = wootters_concurrence(rho)
    if not np.isfinite(c):
        return 0.0
    return max(0.0, float(c))


def concurrence_trajectory(psi0, profile, times):
    """Returns mean_adjacent_concurrence(t) for a pure initial state."""
    rho0 = np.outer(psi0, psi0.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    conc_traj = np.zeros(len(times))
    for i, t in enumerate(times):
        rho_vec = profile['R'] @ (c0 * np.exp(profile['eigvals'] * t))
        rho_t = rho_vec.reshape(D, D)
        rho_t = (rho_t + rho_t.conj().T) / 2
        concs = []
        for k in range(N - 1):
            rho_kk = partial_trace_to_pair(rho_t, N, k, k + 1)
            concs.append(safe_wootters(rho_kk))
        conc_traj[i] = np.mean(concs)
    return conc_traj


def auc_windows(t, conc_traj, windows=(2.0, 10.0, 30.0)):
    """Returns dict {auc_T: integral} via trapezoidal rule."""
    result = {}
    for T in windows:
        mask = t <= T
        result[f"auc_{int(T)}"] = float(np.trapezoid(conc_traj[mask], t[mask]))
    return result


def dominant_mode_info(psi, profile, overlap_threshold=0.001):
    """Rate, overlap, and <n_XY> of the dominant non-stationary eigenmode."""
    rho0 = np.outer(psi, psi.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    probs = np.abs(c0)**2
    rates = -profile['eigvals'].real
    nonstat = rates > 1e-10

    if not nonstat.any():
        return {'rate': float('nan'), 'overlap': 0.0, 'n_xy': 0.0}

    idxs = np.where(nonstat)[0]
    dom_i = idxs[np.argmax(probs[idxs])]
    return {
        'rate': float(profile['eigvals'][dom_i].real),
        'overlap': float(probs[dom_i]),
    }


def evaluate_state(name, psi, times):
    """Full evaluation of a candidate state under both profiles."""
    rho0 = np.outer(psi, psi.conj())
    c_init_max = float(max_adjacent_concurrence(rho0, N))

    result = {
        'name': name,
        'initial_conc_max': c_init_max,
    }
    for plabel, profile in PROFILES.items():
        traj = concurrence_trajectory(psi, profile, times)
        auc = auc_windows(times, traj)
        dom = dominant_mode_info(psi, profile)
        result[plabel] = {**auc, 'dominant_mode': dom}
        result[f'{plabel}_traj'] = traj  # kept in memory, not serialized

    return result


# ===================================================================
# State constructors
# ===================================================================
def w_on_subset(sites):
    """W-state on given sites, |0> elsewhere."""
    psi = np.zeros(D, dtype=complex)
    for k in sites:
        psi[1 << (N - 1 - k)] = 1.0
    return psi / np.linalg.norm(psi)


def bell_plus_at(i, j):
    """Bell+ pair on adjacent sites (i,j), |+> on all others."""
    assert abs(i - j) == 1 and 0 <= min(i, j) and max(i, j) < N
    i, j = min(i, j), max(i, j)
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    parts = []
    k = 0
    while k < N:
        if k == i:
            parts.append(bell)
            k += 2
        else:
            parts.append(plus)
            k += 1
    psi = parts[0]
    for p in parts[1:]:
        psi = np.kron(psi, p)
    return psi / np.linalg.norm(psi)


def two_excitation_symmetric(pairs=None):
    """Symmetric superposition of two-excitation basis states."""
    psi = np.zeros(D, dtype=complex)
    if pairs is None:
        pairs = list(combinations(range(N), 2))
    for (i, j) in pairs:
        idx = (1 << (N - 1 - i)) | (1 << (N - 1 - j))
        psi[idx] = 1.0
    return psi / np.linalg.norm(psi)


def make_plus5():
    """|+>^5 product state."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus.copy()
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    return psi


def make_w5():
    """W5 as a pure state vector (not density matrix)."""
    return w_on_subset(range(N))


def sacrifice_tuned_w5():
    """W5 with amplitude suppressed on sacrifice qubit (site 0)."""
    weights = np.sqrt(gamma_sacrifice.min() / gamma_sacrifice)
    psi = np.zeros(D, dtype=complex)
    for k in range(N):
        psi[1 << (N - 1 - k)] = weights[k]
    return psi / np.linalg.norm(psi)


# ===================================================================
# Build all candidates
# ===================================================================
def build_candidates():
    """Returns list of (name, psi) tuples."""
    candidates = []

    # Class 1: W5 baseline
    candidates.append(("W5_full", make_w5()))

    # Class 2: W-on-subsets
    for pair in [(0,1), (1,2), (2,3), (3,4)]:
        candidates.append((f"W2_sites_{pair[0]}{pair[1]}", w_on_subset(pair)))
    candidates.append(("W4_sites_0123", w_on_subset([0,1,2,3])))
    candidates.append(("W4_sites_1234", w_on_subset([1,2,3,4])))

    # Class 3: Shifted Bell pairs
    for pair in [(0,1), (1,2), (2,3), (3,4)]:
        candidates.append((f"Bell_at_{pair[0]}{pair[1]}", bell_plus_at(*pair)))

    # Class 4: Symmetric W-superpositions
    w5 = make_w5()
    vacuum = np.zeros(D, dtype=complex)
    vacuum[0] = 1.0
    psi_wv = w5 + vacuum
    psi_wv /= np.linalg.norm(psi_wv)
    candidates.append(("W5_plus_vacuum", psi_wv))

    psi_wp = w5 + make_plus5()
    psi_wp /= np.linalg.norm(psi_wp)
    candidates.append(("W5_plus_plus5", psi_wp))

    # Class 5: Two-excitation symmetric
    candidates.append(("W5_2_all", two_excitation_symmetric()))
    candidates.append(("W5_2_nearest", two_excitation_symmetric(
        [(0,1),(1,2),(2,3),(3,4)])))
    candidates.append(("W5_2_far", two_excitation_symmetric(
        [(0,2),(1,3),(2,4),(0,3),(0,4),(1,4)])))

    # Class 6: Sacrifice-aware
    candidates.append(("sacrifice_tuned_W5", sacrifice_tuned_w5()))

    return candidates


# ===================================================================
# Experiment B: surrogate-loss optimization
# ===================================================================
def effective_rate_surrogate(psi, profile, overlap_threshold=0.05):
    """Negative of effective slow rate with initial-concurrence guard."""
    rho0 = np.outer(psi, psi.conj())
    c0 = profile['R_inv'] @ rho0.flatten()
    probs = np.abs(c0)**2
    total = probs.sum()
    if total < 1e-12:
        return 1e6
    normed = probs / total

    significant = normed > overlap_threshold
    rates = -profile['eigvals'].real
    nonstat = rates > 1e-10

    valid = significant & nonstat
    if not valid.any():
        return 1e6
    effective_rate = float(rates[valid].min())

    conc0 = max_adjacent_concurrence(rho0, N)
    if conc0 < 0.1:
        return effective_rate + 10.0 * (0.1 - conc0)
    return effective_rate


def run_experiment_b(basis_states, profile, seeds=(42, 137, 271),
                     popsize=20, maxiter=500, times=None):
    """Surrogate optimizer over a given basis of pure states."""
    m = len(basis_states)
    n_params = 2 * m - 1

    def state_from_params(params):
        re = params[:m]
        im = np.zeros(m)
        im[1:] = params[m:]
        coeffs = re + 1j * im
        psi = sum(c * v for c, v in zip(coeffs, basis_states))
        nrm = np.linalg.norm(psi)
        if nrm < 1e-15:
            return np.zeros(D, dtype=complex)
        return psi / nrm

    def obj(params):
        psi = state_from_params(params)
        if np.linalg.norm(psi) < 1e-12:
            return 1e6
        return effective_rate_surrogate(psi, profile)

    bounds = [(-2, 2)] * n_params
    results = []

    for seed in seeds:
        print(f"    seed={seed} ...", end='', flush=True)
        t0 = _time.time()
        de_res = differential_evolution(
            obj, bounds, seed=seed, maxiter=maxiter, popsize=popsize,
            tol=1e-10, polish=True, disp=False)
        dt = _time.time() - t0

        best_psi = state_from_params(de_res.x)
        surrogate_val = float(de_res.fun)

        # Real AUC check
        if times is not None:
            traj = concurrence_trajectory(best_psi, profile, times)
            auc = auc_windows(times, traj)
        else:
            auc = {}

        rho0 = np.outer(best_psi, best_psi.conj())
        conc0 = max_adjacent_concurrence(rho0, N)

        results.append({
            'seed': seed,
            'surrogate': surrogate_val,
            'conc0': float(conc0),
            'auc': auc,
            'nfev': int(de_res.nfev),
            'time_s': dt,
            'psi': best_psi,
        })
        auc_str = f"AUC(10)={auc.get('auc_10','?')}" if auc else "no AUC"
        print(f"  surr={surrogate_val:.4f}  C0={conc0:.3f}  {auc_str}  ({dt:.0f}s)")

    return results


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    t_start = _time.time()

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "results", "sacrifice_geometry_phase2")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("SACRIFICE GEOMETRY PHASE 2: Can We Beat W5?")
    print("=" * 80)
    print(f"  N = {N}, chain [80, 8, 79, 53, 85]")
    print(f"  gamma_sacrifice = {np.round(gamma_sacrifice, 5)}")
    print(f"  Sg = {Sg_sac:.5f}")

    times = np.linspace(0, 30, 301)

    # ==================================================================
    # Experiment A: Candidate scan
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("EXPERIMENT A: Structured candidate scan")
    print("=" * 80)

    candidates = build_candidates()
    print(f"  {len(candidates)} candidates constructed\n")

    results = []
    for name, psi in candidates:
        t0 = _time.time()
        r = evaluate_state(name, psi, times)
        dt = _time.time() - t0
        results.append(r)
        sac = r['sacrifice']
        print(f"  {name:<22}  AUC(10)={sac['auc_10']:.4f}  "
              f"AUC(30)={sac['auc_30']:.4f}  "
              f"C0={r['initial_conc_max']:.3f}  "
              f"dom={sac['dominant_mode']['rate']:+.3f}  ({dt:.1f}s)")

    # Sort by sacrifice AUC(T=10), descending
    results.sort(key=lambda r: r['sacrifice']['auc_10'], reverse=True)

    w5_auc10 = None
    for r in results:
        if r['name'] == 'W5_full':
            w5_auc10 = r['sacrifice']['auc_10']
            break

    print(f"\n  W5 baseline AUC(T=10) = {w5_auc10:.4f}")
    print(f"\n  --- SORTED BY AUC(T=10) SACRIFICE ---")
    print(f"  {'#':>2}  {'State':<22}  {'AUC(2)':>8}  {'AUC(10)':>8}  {'AUC(30)':>8}"
          f"  {'C_init':>6}  {'dom_rate':>9}  {'beat W5?':>8}")
    print(f"  {'-' * 85}")
    for i, r in enumerate(results):
        sac = r['sacrifice']
        beat = "YES" if sac['auc_10'] > w5_auc10 + 0.001 else "no"
        print(f"  {i+1:>2}  {r['name']:<22}  {sac['auc_2']:>8.4f}  "
              f"{sac['auc_10']:>8.4f}  {sac['auc_30']:>8.4f}  "
              f"{r['initial_conc_max']:>6.3f}  "
              f"{sac['dominant_mode']['rate']:>+9.4f}  {beat:>8}")

    # --- Persist JSON ---
    json_data = {}
    for r in results:
        entry = {
            'initial_conc_max': r['initial_conc_max'],
            'dominant_mode_sacrifice': r['sacrifice']['dominant_mode'],
        }
        for pl in ['sacrifice', 'uniform']:
            entry[pl] = {k: v for k, v in r[pl].items() if k != 'dominant_mode'}
            entry[pl]['dominant_mode'] = r[pl]['dominant_mode']
        json_data[r['name']] = entry

    json_path = os.path.join(out_dir, "candidate_scan.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(json_data, f, indent=2)
    print(f"\n  Saved {json_path}")

    # --- Persist sorted text ---
    txt_path = os.path.join(out_dir, "candidate_scan_sorted.txt")
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Sacrifice Geometry Phase 2 - Candidate Scan (sorted by sacrifice AUC(T=10))\n")
        f.write(f"Computed: {np.datetime64('today')}\n")
        f.write(f"W5 baseline AUC(T=10) = {w5_auc10:.4f}\n\n")
        f.write(f"{'#':>2}  {'State':<22}  {'AUC(2)':>8}  {'AUC(10)':>8}  {'AUC(30)':>8}"
                f"  {'C_init':>6}  {'dom_rate':>9}  beat_W5\n")
        f.write(f"{'-' * 88}\n")
        for i, r in enumerate(results):
            sac = r['sacrifice']
            beat = "YES" if sac['auc_10'] > w5_auc10 + 0.001 else "no"
            f.write(f"{i+1:>2}  {r['name']:<22}  {sac['auc_2']:>8.4f}  "
                    f"{sac['auc_10']:>8.4f}  {sac['auc_30']:>8.4f}  "
                    f"{r['initial_conc_max']:>6.3f}  "
                    f"{sac['dominant_mode']['rate']:>+9.4f}  {beat}\n")
        # Also uniform
        f.write(f"\n\nUNIFORM PROFILE (same order)\n")
        f.write(f"{'#':>2}  {'State':<22}  {'AUC(2)':>8}  {'AUC(10)':>8}  {'AUC(30)':>8}"
                f"  {'dom_rate':>9}\n")
        f.write(f"{'-' * 70}\n")
        for i, r in enumerate(results):
            uni = r['uniform']
            f.write(f"{i+1:>2}  {r['name']:<22}  {uni['auc_2']:>8.4f}  "
                    f"{uni['auc_10']:>8.4f}  {uni['auc_30']:>8.4f}  "
                    f"{uni['dominant_mode']['rate']:>+9.4f}\n")
    print(f"  Saved {txt_path}")

    # --- Plots ---
    print(f"\n  Generating plots ...")
    colors_top = ['#2ca02c', '#d62728', '#1f77b4', '#ff7f0e', '#9467bd',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # Individual plots for top 5
    for i, r in enumerate(results[:5]):
        fig, ax = plt.subplots(figsize=(8, 5))
        t_sac = r['sacrifice_traj']
        t_uni = r['uniform_traj']
        ax.plot(times, t_sac, label='sacrifice', linewidth=2, color=colors_top[i])
        ax.plot(times, t_uni, label='uniform', linewidth=2, linestyle='--',
                color=colors_top[i], alpha=0.7)
        ax.set_xlabel('t (dimensionless)')
        ax.set_ylabel('mean adjacent concurrence')
        ax.set_title(f"{r['name']}  (AUC(10)={r['sacrifice']['auc_10']:.4f})")
        ax.set_xlim(0, 30)
        ax.set_ylim(bottom=0)
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"concurrence_vs_time_{r['name']}.png"), dpi=150)
        plt.close(fig)

    # Summary: top 5 + W5 reference (sacrifice only)
    fig, ax = plt.subplots(figsize=(10, 6))
    plotted_w5 = False
    for i, r in enumerate(results[:5]):
        t_sac = r['sacrifice_traj']
        lw = 2.5 if r['name'] == 'W5_full' else 1.5
        ls = '-' if r['name'] == 'W5_full' else '-'
        ax.plot(times, t_sac, label=r['name'], linewidth=lw,
                color=colors_top[i], linestyle=ls)
        if r['name'] == 'W5_full':
            plotted_w5 = True
    if not plotted_w5:
        # W5 as dashed reference
        for r in results:
            if r['name'] == 'W5_full':
                ax.plot(times, r['sacrifice_traj'], label='W5_full (ref)',
                        linewidth=2.5, color='#2ca02c', linestyle='--', alpha=0.7)
                break
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('mean adjacent concurrence')
    ax.set_title('Top 5 candidates (sacrifice profile)')
    ax.set_xlim(0, 30)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'top5_comparison.png'), dpi=150)
    plt.close(fig)
    print("  Plots saved")

    # ==================================================================
    # Experiment B decision
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("EXPERIMENT B: Decision")
    print("=" * 80)

    # Check if any candidate beat W5 or came close
    best_non_w5 = None
    for r in results:
        if r['name'] != 'W5_full':
            if best_non_w5 is None or r['sacrifice']['auc_10'] > best_non_w5['sacrifice']['auc_10']:
                best_non_w5 = r

    gap = w5_auc10 - best_non_w5['sacrifice']['auc_10']
    print(f"  Best non-W5: {best_non_w5['name']} with AUC(10) = "
          f"{best_non_w5['sacrifice']['auc_10']:.4f}")
    print(f"  Gap to W5: {gap:.4f}")

    # Determine which class to optimize
    # Run Experiment B if sacrifice-tuned W5 or W4_sites_1234 is competitive
    # or if any weighted-W variant is within 20% of W5
    run_b = False
    opt_basis_label = None
    opt_basis_states = None

    # Check if weighted-W variants are worth refining
    sacrifice_tuned = None
    w4_1234 = None
    for r in results:
        if r['name'] == 'sacrifice_tuned_W5':
            sacrifice_tuned = r
        if r['name'] == 'W4_sites_1234':
            w4_1234 = r

    if sacrifice_tuned and sacrifice_tuned['sacrifice']['auc_10'] > w5_auc10 * 0.8:
        run_b = True
        opt_basis_label = "weighted single-excitation (5 sites + vacuum)"
        # Basis: |00000>, |10000>, |01000>, |00100>, |00010>, |00001>
        opt_basis_states = [np.zeros(D, dtype=complex) for _ in range(N + 1)]
        opt_basis_states[0][0] = 1.0  # vacuum
        for k in range(N):
            opt_basis_states[k + 1][1 << (N - 1 - k)] = 1.0
        print(f"  Running Experiment B: {opt_basis_label}")
    else:
        # Also check two-excitation states
        best_2exc = max((r for r in results if r['name'].startswith('W5_2')),
                        key=lambda r: r['sacrifice']['auc_10'], default=None)
        if best_2exc and best_2exc['sacrifice']['auc_10'] > w5_auc10 * 0.5:
            run_b = True
            opt_basis_label = "weight-1 + weight-2 subspace"
            basis = []
            for i in range(D):
                hw = bin(i).count('1')
                if hw <= 2:
                    v = np.zeros(D, dtype=complex)
                    v[i] = 1.0
                    basis.append(v)
            opt_basis_states = basis
            print(f"  Running Experiment B: {opt_basis_label} ({len(basis)} states)")

    if not run_b:
        print("  No candidate close enough to W5. Skipping Experiment B.")

    expb_results = None
    if run_b:
        print(f"\n{'=' * 80}")
        print(f"EXPERIMENT B: Surrogate optimization ({opt_basis_label})")
        print("=" * 80)
        expb_results = run_experiment_b(
            opt_basis_states, PROFILES['sacrifice'],
            seeds=(42, 137, 271), popsize=20, maxiter=500, times=times)

        # Summary
        print(f"\n  --- Experiment B results ---")
        print(f"  {'seed':>6}  {'surrogate':>10}  {'C0':>6}  {'AUC(10)':>10}  {'beat W5?':>8}")
        for eb in expb_results:
            beat = "YES" if eb['auc'].get('auc_10', 0) > w5_auc10 + 0.001 else "no"
            print(f"  {eb['seed']:>6}  {eb['surrogate']:>10.4f}  {eb['conc0']:>6.3f}  "
                  f"{eb['auc'].get('auc_10', 0):>10.4f}  {beat:>8}")

        # Save
        eb_json = []
        for eb in expb_results:
            eb_json.append({
                'seed': eb['seed'], 'surrogate': eb['surrogate'],
                'conc0': eb['conc0'], 'auc': eb['auc'],
                'nfev': eb['nfev'], 'time_s': eb['time_s'],
            })
        eb_path = os.path.join(out_dir, "optimizer_refinement.json")
        with open(eb_path, 'w', encoding='utf-8') as f:
            json.dump(eb_json, f, indent=2)
        print(f"  Saved {eb_path}")

    # ==================================================================
    # Final summary
    # ==================================================================
    print(f"\n{'=' * 80}")
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"  W5 AUC(T=10) = {w5_auc10:.4f}  (sacrifice profile)")
    print(f"  Best non-W5:   {best_non_w5['name']} = "
          f"{best_non_w5['sacrifice']['auc_10']:.4f}")
    if expb_results:
        best_eb = max(expb_results, key=lambda x: x['auc'].get('auc_10', 0))
        print(f"  Best optimizer: seed={best_eb['seed']}  AUC(10)="
              f"{best_eb['auc'].get('auc_10', 0):.4f}  surr={best_eb['surrogate']:.4f}")

    beaten = any(r['sacrifice']['auc_10'] > w5_auc10 + 0.001
                 for r in results if r['name'] != 'W5_full')
    if beaten:
        winners = [r['name'] for r in results
                   if r['name'] != 'W5_full' and r['sacrifice']['auc_10'] > w5_auc10 + 0.001]
        print(f"\n  STATES THAT BEAT W5: {winners}")
    else:
        print(f"\n  NO STATE BEATS W5 on AUC(T=10) under sacrifice.")

    t_end = _time.time()
    print(f"\nTotal runtime: {t_end - t_start:.1f}s")
