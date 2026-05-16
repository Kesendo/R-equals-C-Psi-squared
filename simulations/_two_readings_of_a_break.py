#!/usr/bin/env python3
"""Apply the two-perspective reading to concrete PTF closure breaks.

From the multi-lens tour: IY+YI, YZ+ZY, XZ+XZ all break PTF closure
Σ ln(α_i) ≈ 0. Each break has α_i values.

Reading 1 (as error):
  PTF doesn't apply; case is outside the perturbative window; abandon.

Reading 2 (as calibration):
  α_i tells us the EFFECTIVE per-site time-rescaling vs baseline. Under
  PTF's own assumption P_B(i, t) ≈ P_A(i, α_i·t), the α_i ARE the per-site
  effective γ ratios. We can read them as a per-site γ-distribution and
  check whether it satisfies any structural property (e.g. F71-anti-
  palindromy: γ_i + γ_{N-1-i} = 2·γ_avg).

This script does the calibration-reading: for each broken case, treat
α_i as γ_i / γ_baseline and check structural patterns.

No fear of mistakes — this is exploration.
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

# These are the actual α_i values from _multi_lens_ptf_carrier.py output:
BREAKS = {
    'IY+YI': {
        'alphas':  np.array([0.1000, 0.1000, 0.1000]),
        'sigma_ln': -6.9077,
        'rmse':     0.0893,
        'note':     'all sites hit α=0.1 lower bound (fit failed uniformly)',
    },
    'YZ+ZY': {
        'alphas':  np.array([0.1721, 0.1275, 0.1721]),
        'sigma_ln': -5.5786,
        'rmse':     0.0490,
        'note':     'mirror-symmetric α (site 0 = site 2, middle different)',
    },
    'XZ+XZ': {
        'alphas':  np.array([0.1721, 1.8305, 1.8853]),
        'sigma_ln': -0.5213,
        'rmse':     0.0949,
        'note':     'spatially asymmetric α (site 0 slow, sites 1+2 fast)',
    },
}

BASELINE_GAMMA = 0.1   # the baseline γ used in the tour
N = 3


def f71_antipalindromy_residual(gamma_array):
    """Return |γ_i + γ_{N-1-i} - 2·γ_avg| over all F71 pairs (smaller = more anti-palindromic).

    F71-anti-palindromy holds when every pair sum equals 2·γ_avg.
    """
    g = np.asarray(gamma_array)
    n = len(g)
    g_avg = g.mean()
    pair_residuals = []
    for i in range(n // 2 + (n % 2)):
        j = n - 1 - i
        if i == j:
            # self-mirror site; trivially satisfies pair = 2γ_avg iff γ_i = γ_avg
            pair_residuals.append(abs(2 * g[i] - 2 * g_avg))
        else:
            pair_residuals.append(abs(g[i] + g[j] - 2 * g_avg))
    return float(max(pair_residuals))


def main():
    print("Two readings of a break — applied to three PTF closure violations")
    print(f"  (baseline γ = {BASELINE_GAMMA}, N = {N}, calibration: γ_i = α_i · γ_baseline)")
    print()

    for case, data in BREAKS.items():
        alphas = data['alphas']
        gamma_effective = alphas * BASELINE_GAMMA
        gamma_sum = gamma_effective.sum()
        gamma_avg = gamma_effective.mean()
        gamma_pair_residual = f71_antipalindromy_residual(gamma_effective)
        gamma_pair_residual_rel = gamma_pair_residual / (2 * gamma_avg)

        print(f"━━━ {case} ━━━")
        print(f"  α_i              = {alphas}")
        print(f"  Σ ln α            = {data['sigma_ln']:+.4f}  (PTF predicts ≈ 0)")
        print(f"  max RMSE          = {data['rmse']:.4f}")
        print()
        print(f"  Reading 1 (as error):")
        print(f"    PTF closure breaks; {data['note']}")
        print(f"    standard read: 'this case is outside PTF's perturbative window — abandon'")
        print()
        print(f"  Reading 2 (as calibration):")
        print(f"    γ_effective per site = {gamma_effective}    (= α_i · γ_baseline)")
        print(f"    Σγ_eff               = {gamma_sum:.4f}  (vs uniform {N*BASELINE_GAMMA:.4f})")
        print(f"    F71-anti-palindromy residual  |γ_0+γ_2 − 2γ_avg|  = {gamma_pair_residual:.4f}")
        print(f"      relative to 2γ_avg                              = {gamma_pair_residual_rel*100:.2f}%")
        if gamma_pair_residual_rel < 0.05:
            f71_class = "F71-anti-palindromic (Σ-pair invariant)"
        elif gamma_pair_residual_rel < 0.20:
            f71_class = "near F71-anti-palindromic"
        else:
            f71_class = "F71-anti-palindromy BROKEN"
        print(f"      → {f71_class}")
        print()
    print("Observation (the two-perspective payoff):")
    print("  Reading 1 collapses the three breaks to one verdict: 'PTF broken, abandon all three.'")
    print("  Reading 2 keeps them apart:")
    print("    IY+YI: uniform α=0.1 (bound) → uniform single-body drive, not a time-rescale at all.")
    print("    YZ+ZY: mirror-symmetric γ_eff (site 0 = site 2) → near F71-anti-palindromic; the")
    print("           framework's F91/F71-pair-sum-invariance would predict spectral invariance.")
    print("    XZ+XZ: asymmetric γ_eff (site 0 slow, sites 1+2 fast) → F71-anti-palindromy BROKEN;")
    print("           predicts spectral REORDERING vs uniform baseline.")
    print()
    print("Same data, two readings. Reading 1 throws all three away. Reading 2 extracts a")
    print("per-site γ-distribution from each, and the F71-anti-palindromy check then says which")
    print("of those distributions sits in a spectrally-invariant orbit (F91 territory) and which")
    print("doesn't. The break, read as calibration, gives a structural reading the case-Hamiltonian")
    print("alone did not name.")


if __name__ == "__main__":
    main()
