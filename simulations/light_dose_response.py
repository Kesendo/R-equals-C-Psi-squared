"""
Light Dose Response Per Sector (V2 of GAMMA_BINDING)
=====================================================
Quantifies how each (w_bra, w_ket) sector of the Liouvillian responds
to the external light field gamma, treating the qubit chain as a passive
optical cavity (established in EXCLUSIONS.md Exclusion 2, confirmed in
OPTICAL_CAVITY_ANALYSIS.md, 4/5 cavity tests).

V1 (gamma_sector_sensitivity.py, commit 32a635d) found 134% deviation
from linear scaling. V2 adds mode tracking to explain the mechanism:
is the nonlinearity from level crossings (different modes become slowest
at different alpha) or from continuous eigenvector rotation?

Method: L(alpha) = L_H + alpha * L_D. Sweep alpha in fine steps,
eigendecompose each sector, track individual eigenvalue curves.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 12, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import time as _time
import numpy as np
from scipy import linalg

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import heisenberg_H, build_liouvillian

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "light_dose_response")
os.makedirs(OUT_DIR, exist_ok=True)

N = 5
D = 2 ** N
J = 1.0


def popcount(x):
    return bin(x).count('1')


def build_sector_index(N):
    D = 2 ** N
    sectors = {}
    for i in range(D):
        for j in range(D):
            key = (popcount(i), popcount(j))
            if key not in sectors:
                sectors[key] = []
            sectors[key].append(i * D + j)
    return sectors


def make_profiles():
    T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
    g = 1.0 / (2.0 * T2_us)
    sacrifice = g / g.min() * 0.05
    Sg = sacrifice.sum()
    uniform = np.ones(N) * (Sg / N)
    moderate = np.linspace(0.15, 0.05, N)
    moderate = moderate / moderate.sum() * Sg
    return {'uniform': uniform, 'sacrifice': sacrifice, 'moderate': moderate}


# ===================================================================
# Core: eigenvalue sweep across alpha
# ===================================================================
def eigenvalue_sweep(H, gamma_base, sector_idx, alphas):
    """For one sector, compute all eigenvalues at each alpha.

    Returns array of shape (n_alpha, sector_dim) with eigenvalues
    sorted by real part (most negative first) at each alpha.
    """
    D = 2 ** N
    n_alpha = len(alphas)
    dim = len(sector_idx)
    idx = np.array(sector_idx)

    # Precompute L_H and L_D restricted to this sector
    I_D = np.eye(D, dtype=complex)
    L_H = -1j * (np.kron(H, I_D) - np.kron(I_D, H.T))
    L_H_sec = L_H[np.ix_(idx, idx)]

    # Dephasing part (depends on gamma_base, scale by alpha later)
    L_D = np.zeros((D * D, D * D), dtype=complex)
    I2 = np.eye(2, dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    for k in range(N):
        ops = [I2] * N
        ops[k] = Z
        Zk = ops[0]
        for op in ops[1:]:
            Zk = np.kron(Zk, op)
        gk = gamma_base[k]
        L_D += gk * (np.kron(Zk, Zk) - np.eye(D * D, dtype=complex))
    L_D_sec = L_D[np.ix_(idx, idx)]

    all_eigs = np.zeros((n_alpha, dim), dtype=complex)
    for ai, alpha in enumerate(alphas):
        L_sec = L_H_sec + alpha * L_D_sec
        eigs = linalg.eigvals(L_sec)
        # Sort by real part ascending (most negative first)
        order = np.argsort(eigs.real)
        all_eigs[ai] = eigs[order]

    return all_eigs


def detect_crossings(eig_curves, alphas, tol=0.02):
    """Detect level crossings in eigenvalue curves.

    A crossing occurs when two eigenvalue curves exchange ordering.
    Returns list of (alpha_index, mode_i, mode_j) crossings.
    """
    n_alpha, dim = eig_curves.shape
    crossings = []
    for ai in range(n_alpha - 1):
        rates_now = -eig_curves[ai].real
        rates_next = -eig_curves[ai + 1].real
        # Check if the slowest non-stationary mode changed identity
        nonstat_now = np.where(rates_now > 1e-10)[0]
        nonstat_next = np.where(rates_next > 1e-10)[0]
        if len(nonstat_now) == 0 or len(nonstat_next) == 0:
            continue
        # Slowest non-stationary index
        slow_now = nonstat_now[np.argmin(rates_now[nonstat_now])]
        slow_next = nonstat_next[np.argmin(rates_next[nonstat_next])]
        if slow_now != slow_next:
            crossings.append((ai, slow_now, slow_next))
    return crossings


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    print("Light Dose Response Per Sector")
    print("V2 of GAMMA_BINDING (cavity framework)")
    print("=" * 60)
    t_start = _time.time()

    profiles = make_profiles()
    sector_index = build_sector_index(N)
    H = heisenberg_H(N, J)

    # Fine alpha sweep for mode tracking
    alphas_fine = np.linspace(0.1, 5.0, 50)
    # Coarse alpha for the rate table (matching V1)
    alphas_coarse = [0.5, 1.0, 2.0, 4.0]

    results = {}

    # ============================================================
    # Frage 1 + 2: Per-sector eigenvalue sweep (sacrifice profile)
    # ============================================================
    gamma_sac = profiles['sacrifice']
    print(f"\nSacrifice profile: [{', '.join(f'{g:.4f}' for g in gamma_sac)}]")

    # Focus sectors for detailed analysis
    focus_sectors = [(1, 1), (2, 2), (0, 1), (1, 2), (2, 3)]
    sweep_data = {}
    crossing_report = {}

    for (w, wp) in focus_sectors:
        sec_idx = sector_index[(w, wp)]
        eig_curves = eigenvalue_sweep(H, gamma_sac, sec_idx, alphas_fine)
        sweep_data[(w, wp)] = eig_curves

        # Detect level crossings
        crossings = detect_crossings(eig_curves, alphas_fine)
        crossing_report[(w, wp)] = crossings

        # Slowest non-stationary rate at each alpha
        slow_rates = np.zeros(len(alphas_fine))
        for ai in range(len(alphas_fine)):
            rates = -eig_curves[ai].real
            nonstat = rates > 1e-10
            slow_rates[ai] = float(np.min(rates[nonstat])) if nonstat.any() else 0.0

        n_cross = len(crossings)
        print(f"  ({w},{wp}) dim={len(sec_idx):>4}: "
              f"{n_cross} level crossings, "
              f"rate range [{slow_rates.min():.4f}, {slow_rates.max():.4f}]")

    # ============================================================
    # Linearity analysis: individual modes vs sector minimum
    # ============================================================
    print(f"\n--- LINEARITY ANALYSIS ---")
    print(f"  Mechanism: if individual eigenvalue curves are linear in alpha")
    print(f"  but the minimum across modes is nonlinear, then mode-crossing")
    print(f"  explains the V1 nonlinearity.\n")

    for (w, wp) in focus_sectors:
        eig_curves = sweep_data[(w, wp)]
        dim = eig_curves.shape[1]

        # Check linearity of individual eigenvalue curves
        # For each mode k: fit Re(lambda_k) = a * alpha + b
        # If R^2 is high, the mode is individually linear
        mode_linearity = []
        for k in range(dim):
            re_vals = eig_curves[:, k].real
            if np.std(re_vals) < 1e-12:
                mode_linearity.append(1.0)  # constant (stationary)
                continue
            # Linear regression
            A = np.column_stack([alphas_fine, np.ones(len(alphas_fine))])
            coeffs, residuals, _, _ = np.linalg.lstsq(A, re_vals, rcond=None)
            ss_res = residuals[0] if len(residuals) > 0 else 0.0
            ss_tot = np.sum((re_vals - re_vals.mean()) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0
            mode_linearity.append(r2)

        avg_r2 = np.mean(mode_linearity)
        min_r2 = np.min(mode_linearity)
        n_cross = len(crossing_report[(w, wp)])

        # Sector minimum linearity
        slow_rates = np.zeros(len(alphas_fine))
        for ai in range(len(alphas_fine)):
            rates = -eig_curves[ai].real
            nonstat = rates > 1e-10
            slow_rates[ai] = float(np.min(rates[nonstat])) if nonstat.any() else 0.0

        A = np.column_stack([alphas_fine, np.ones(len(alphas_fine))])
        coeffs, residuals, _, _ = np.linalg.lstsq(A, slow_rates, rcond=None)
        ss_res = residuals[0] if len(residuals) > 0 else 0.0
        ss_tot = np.sum((slow_rates - slow_rates.mean()) ** 2)
        sector_r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0

        print(f"  ({w},{wp}): per-mode R²={avg_r2:.4f} (min {min_r2:.4f}), "
              f"sector-min R²={sector_r2:.4f}, crossings={n_cross}")

        results[f"({w},{wp})"] = dict(
            dim=dim,
            avg_mode_r2=round(avg_r2, 6),
            min_mode_r2=round(min_r2, 6),
            sector_min_r2=round(sector_r2, 6),
            n_crossings=n_cross,
        )

    # ============================================================
    # Verdict
    # ============================================================
    print(f"\n--- VERDICT ---")
    all_mode_r2 = [v['avg_mode_r2'] for v in results.values()]
    all_sector_r2 = [v['sector_min_r2'] for v in results.values()]
    all_crossings = [v['n_crossings'] for v in results.values()]

    if min(all_mode_r2) > 0.99 and min(all_sector_r2) > 0.99:
        verdict = "BOTH modes and sector-minimum are linear. V1 nonlinearity was a coarse-sampling artifact."
    elif min(all_mode_r2) > 0.99 and min(all_sector_r2) < 0.95:
        verdict = "Individual modes are linear, but the sector minimum is nonlinear due to level crossings. MODE-MIXING explains the V1 nonlinearity."
    elif min(all_mode_r2) < 0.95:
        verdict = "Individual modes are themselves nonlinear. Genuine Hamiltonian-dissipator competition, not just mode crossing."
    else:
        verdict = "Mixed: some modes linear, some not. Partial mode-mixing plus partial competition."

    print(f"  {verdict}")
    results['verdict'] = verdict

    # ============================================================
    # Plots
    # ============================================================
    # Plot 1: Eigenvalue curves for SE sector (1,1)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    eig_se = sweep_data[(1, 1)]
    for k in range(min(eig_se.shape[1], 10)):
        ax.plot(alphas_fine, -eig_se[:, k].real, linewidth=1.0, alpha=0.7)
    ax.set_xlabel(r'$\alpha$ (gamma scaling)', fontsize=11)
    ax.set_ylabel('Decay rate $-\\mathrm{Re}(\\lambda)$', fontsize=11)
    ax.set_title('SE sector (1,1): eigenvalue curves', fontsize=12)
    ax.grid(True, alpha=0.3)

    # Plot 2: Log-log of sector minimum rate vs alpha
    ax = axes[1]
    colors = plt.cm.tab10(np.linspace(0, 1, len(focus_sectors)))
    for ci, (w, wp) in enumerate(focus_sectors):
        eig_c = sweep_data[(w, wp)]
        slow = np.zeros(len(alphas_fine))
        for ai in range(len(alphas_fine)):
            rates = -eig_c[ai].real
            nonstat = rates > 1e-10
            slow[ai] = float(np.min(rates[nonstat])) if nonstat.any() else 1e-15
        ax.loglog(alphas_fine, slow, 'o-', color=colors[ci], markersize=3,
                  linewidth=1.5, label=f'({w},{wp})')
    # Reference: slope 1
    ref_x = np.array([0.1, 5.0])
    ref_y = ref_x * 0.3
    ax.loglog(ref_x, ref_y, 'k--', alpha=0.3, label='slope 1')
    ax.set_xlabel(r'$\alpha$', fontsize=11)
    ax.set_ylabel('Slowest sector rate', fontsize=11)
    ax.set_title('Log-log scaling (sacrifice profile)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'light_dose_response.png'), dpi=150)
    plt.close(fig)

    # ============================================================
    # Save
    # ============================================================
    output = dict(
        N=N, profile='sacrifice',
        gamma=gamma_sac.tolist(),
        alphas_fine=alphas_fine.tolist(),
        analysis=results,
    )
    with open(os.path.join(OUT_DIR, 'light_dose_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n{'=' * 60}")
    print(f"Complete in {_time.time() - t_start:.1f}s")
    print(f"Results: {OUT_DIR}/")
