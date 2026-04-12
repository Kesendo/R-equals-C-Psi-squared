"""
Gamma as Binding Parameter: per-sector rate sensitivity (EQ-008)
================================================================
For N=5 Heisenberg chain under Z-dephasing, compute the slowest
decay rate in each (w_bra, w_ket) sector as a function of:
  1. Three gamma profiles (uniform, sacrifice, moderate gradient)
  2. Four global scaling factors alpha in {0.5, 1.0, 2.0, 4.0}

Tests whether gamma is a simple linear scale parameter or whether
different sectors respond with differential sensitivity.

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
                       "results", "gamma_binding")
os.makedirs(OUT_DIR, exist_ok=True)

N = 5
D = 2 ** N
J = 1.0


def popcount(x):
    return bin(x).count('1')


# Precompute sector membership for each basis element |i><j|
def build_sector_index(N):
    """Map each (w_bra, w_ket) to list of flat indices i*D+j."""
    D = 2 ** N
    sectors = {}
    for i in range(D):
        for j in range(D):
            w = popcount(i)
            wp = popcount(j)
            key = (w, wp)
            if key not in sectors:
                sectors[key] = []
            sectors[key].append(i * D + j)
    return sectors


def slowest_rate_per_sector(L_full, sector_index):
    """Extract slowest non-stationary rate for each (w, w') sector."""
    results = {}
    for (w, wp), flat_idx in sorted(sector_index.items()):
        idx = np.array(flat_idx)
        L_sec = L_full[np.ix_(idx, idx)]
        eigvals = linalg.eigvals(L_sec)
        rates = -eigvals.real
        nonstat = rates > 1e-10
        if nonstat.any():
            slowest = float(np.min(rates[nonstat]))
        else:
            slowest = 0.0  # stationary sector (diagonal, w=w')
        results[(w, wp)] = slowest
    return results


# Gamma profiles
def make_profiles():
    """Three gamma profiles for N=5, all with same sum for comparability."""
    # IBM Torino T2-derived sacrifice
    T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
    g = 1.0 / (2.0 * T2_us)
    sacrifice = g / g.min() * 0.05

    # Uniform (matched sum)
    Sg = sacrifice.sum()
    uniform = np.ones(N) * (Sg / N)

    # Moderate gradient
    moderate = np.linspace(0.15, 0.05, N)
    moderate = moderate / moderate.sum() * Sg  # match sum

    return {
        'uniform': uniform,
        'sacrifice': sacrifice,
        'moderate': moderate,
    }


if __name__ == "__main__":
    print("Gamma as Binding Parameter: Per-Sector Sensitivity")
    print("=" * 60)
    t_start = _time.time()

    profiles = make_profiles()
    sector_index = build_sector_index(N)
    alphas = [0.5, 1.0, 2.0, 4.0]
    H = heisenberg_H(N, J)

    all_results = {}

    # ============================================================
    # FRAGE 1: Per-sector rates for three profiles (alpha=1.0)
    # ============================================================
    print("\n--- FRAGE 1: Per-sector rates (alpha=1.0) ---\n")

    for pname, gamma in profiles.items():
        L = build_liouvillian(H, gamma)
        rates = slowest_rate_per_sector(L, sector_index)
        all_results[pname] = {'alpha_1.0': rates}
        print(f"  Profile: {pname}")
        print(f"  gamma = [{', '.join(f'{g:.4f}' for g in gamma)}]")
        print(f"  Sg = {gamma.sum():.4f}")

    # Print comparison table for diagonal sectors (the "exit" sectors)
    print(f"\n  Diagonal sector rates (w = w', these contain the exits):")
    print(f"  {'w':>3}  {'uniform':>10}  {'sacrifice':>10}  {'moderate':>10}")
    print(f"  {'-' * 37}")
    for w in range(N + 1):
        vals = [all_results[p]['alpha_1.0'][(w, w)] for p in profiles]
        print(f"  {w:>3}  {vals[0]:>10.6f}  {vals[1]:>10.6f}  {vals[2]:>10.6f}")

    # Off-diagonal sectors (coherence sectors)
    print(f"\n  Off-diagonal sector rates (w != w', cross-sector coherences):")
    print(f"  {"(w,w')":>6}  {'uniform':>10}  {'sacrifice':>10}  {'moderate':>10}")
    print(f"  {'-' * 40}")
    for w in range(N + 1):
        for wp in range(N + 1):
            if w == wp:
                continue
            vals = [all_results[p]['alpha_1.0'][(w, wp)] for p in profiles]
            # Only print if at least one is nonzero
            if max(vals) > 1e-10:
                print(f"  ({w},{wp})  {vals[0]:>10.6f}  "
                      f"{vals[1]:>10.6f}  {vals[2]:>10.6f}")

    # ============================================================
    # FRAGE 2: Scaling with alpha
    # ============================================================
    print(f"\n--- FRAGE 2: Scaling with alpha ---\n")

    scaling_data = {}
    for pname, gamma_base in profiles.items():
        scaling_data[pname] = {}
        for alpha in alphas:
            gamma = alpha * gamma_base
            L = build_liouvillian(H, gamma)
            rates = slowest_rate_per_sector(L, sector_index)
            scaling_data[pname][alpha] = rates

    # Check linearity: rate(alpha) / rate(1.0) should equal alpha
    print(f"  Linearity check: rate(alpha) / (alpha * rate(1.0))")
    print(f"  If ratio = 1.0 everywhere, scaling is perfectly linear.\n")

    deviations = {}
    for pname in profiles:
        base_rates = scaling_data[pname][1.0]
        max_dev = 0.0
        for alpha in alphas:
            if alpha == 1.0:
                continue
            for (w, wp), rate in scaling_data[pname][alpha].items():
                base = base_rates[(w, wp)]
                if base > 1e-10 and rate > 1e-10:
                    ratio = rate / (alpha * base)
                    dev = abs(ratio - 1.0)
                    if dev > max_dev:
                        max_dev = dev
        deviations[pname] = max_dev
        print(f"  {pname}: max deviation from linearity = {max_dev:.6f}")

    # Detailed linearity table for sacrifice profile
    print(f"\n  Sacrifice profile: rate(alpha) / (alpha * rate(1))")
    print(f"  {"(w,w')":>6}  {'a=0.5':>8}  {'a=2.0':>8}  {'a=4.0':>8}")
    print(f"  {'-' * 35}")
    base = scaling_data['sacrifice'][1.0]
    for w in range(N + 1):
        for wp in range(w, N + 1):
            b = base[(w, wp)]
            if b < 1e-10:
                continue
            ratios = []
            for alpha in [0.5, 2.0, 4.0]:
                r = scaling_data['sacrifice'][alpha][(w, wp)]
                ratios.append(r / (alpha * b) if alpha * b > 1e-10 else 0.0)
            if any(abs(r - 1.0) > 0.001 for r in ratios):
                marker = " <-- nonlinear"
            else:
                marker = ""
            print(f"  ({w},{wp})  {ratios[0]:>8.4f}  {ratios[1]:>8.4f}  "
                  f"{ratios[2]:>8.4f}{marker}")

    # ============================================================
    # Plots
    # ============================================================
    # Plot 1: Heatmap of per-sector rates (sacrifice, alpha=1.0)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, pname in zip(axes, profiles):
        rates = all_results[pname]['alpha_1.0']
        grid = np.zeros((N + 1, N + 1))
        for (w, wp), val in rates.items():
            grid[w, wp] = val
        im = ax.imshow(grid, origin='lower', cmap='viridis',
                       extent=[-0.5, N + 0.5, -0.5, N + 0.5])
        ax.set_title(f'{pname}', fontsize=11)
        ax.set_xlabel('w_ket')
        ax.set_ylabel('w_bra')
        plt.colorbar(im, ax=ax, label='Slowest rate')
    fig.suptitle('Per-sector slowest decay rate (N=5)', fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'sector_rate_heatmaps.png'), dpi=150)
    plt.close(fig)

    # Plot 2: Log-log scaling for selected sectors
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    ci = 0
    for (w, wp) in [(1, 1), (2, 2), (0, 1), (1, 2), (2, 3)]:
        rates_at_alpha = []
        for alpha in alphas:
            rates_at_alpha.append(scaling_data['sacrifice'][alpha][(w, wp)])
        ax.loglog(alphas, rates_at_alpha, 'o-', color=colors[ci],
                  label=f'({w},{wp})', linewidth=1.5, markersize=6)
        ci += 1
    # Reference line: perfect linearity
    ref = [alphas[1] * scaling_data['sacrifice'][1.0][(1, 1)] * a
           for a in [0.3, 5]]
    ax.loglog([0.3, 5], ref, 'k--', alpha=0.3, label='slope=1 (linear)')
    ax.set_xlabel('Scaling factor alpha', fontsize=11)
    ax.set_ylabel('Slowest sector rate', fontsize=11)
    ax.set_title('Rate scaling with gamma (sacrifice profile)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'scaling_loglog.png'), dpi=150)
    plt.close(fig)

    # ============================================================
    # Save results
    # ============================================================
    def key_to_str(d):
        """Convert tuple keys to strings for JSON."""
        if isinstance(d, dict):
            return {str(k): key_to_str(v) for k, v in d.items()}
        return d

    output = dict(
        N=N, profiles={p: g.tolist() for p, g in profiles.items()},
        alphas=alphas,
        rates_alpha1={p: key_to_str(all_results[p]['alpha_1.0'])
                      for p in profiles},
        scaling=key_to_str(scaling_data),
        linearity_deviations=deviations,
    )
    with open(os.path.join(OUT_DIR, 'gamma_binding_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n{'=' * 60}")
    print(f"Complete in {_time.time() - t_start:.1f}s")
    print(f"Results saved to {OUT_DIR}/")
