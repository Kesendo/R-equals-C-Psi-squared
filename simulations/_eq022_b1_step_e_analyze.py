#!/usr/bin/env python3
"""_eq022_b1_step_e_analyze.py — analyze the resonance-shape scan for
universality.

Reads `simulations/results/eq022_resonance_shape/curves.json` (produced by
`_eq022_b1_step_e_resonance_shape.py`) and tests several normalisation
schemes for the abs(K_CC_pr)(Q) curves to see which (if any) collapse the
curves across different (c, N) onto a universal shape.

Normalisations tested:
  (A)  x = Q − Q_peak,           y = K(Q) / |K|max
  (B)  x = (Q − Q_peak)/Q_peak,  y = K(Q) / |K|max
  (C)  x = (Q − Q_peak)/HWHM_avg, y = K(Q) / |K|max  (HWHM_avg = mean of
       left & right half-widths)

For each normalisation, compute pairwise overlap of the curves on a
common x-grid via interpolation; report mean residual and worst case.
The lowest residual is the candidate universal scheme.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "simulations" / "results" / "eq022_resonance_shape"


def find_peak_and_hwhm(Q_grid, K_curve):
    Q_grid = np.array(Q_grid)
    K = np.array(K_curve)
    i_max = int(np.argmax(K))
    K_max = float(K[i_max])
    # Parabolic interp
    if 0 < i_max < len(Q_grid) - 1:
        x = Q_grid[i_max - 1: i_max + 2]
        y = K[i_max - 1: i_max + 2]
        coefs = np.polyfit(x, y, 2)
        if coefs[0] < 0:
            Q_star = -coefs[1] / (2 * coefs[0])
            K_max_interp = coefs[2] - coefs[1] ** 2 / (4 * coefs[0])
            if abs(Q_star - Q_grid[i_max]) <= (Q_grid[i_max + 1] - Q_grid[i_max - 1]):
                K_max = float(K_max_interp)
            else:
                Q_star = float(Q_grid[i_max])
        else:
            Q_star = float(Q_grid[i_max])
    else:
        Q_star = float(Q_grid[i_max])

    half = K_max / 2.0
    hwhm_l = None
    hwhm_r = None
    for i in range(i_max, -1, -1):
        if K[i] < half:
            x0, x1 = Q_grid[i], Q_grid[i + 1]
            y0, y1 = K[i], K[i + 1]
            if y1 == y0:
                x_half = x0
            else:
                x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_l = float(Q_star - x_half)
            break
    for i in range(i_max, len(Q_grid)):
        if K[i] < half:
            x0, x1 = Q_grid[i - 1], Q_grid[i]
            y0, y1 = K[i - 1], K[i]
            if y1 == y0:
                x_half = x1
            else:
                x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_r = float(x_half - Q_star)
            break
    return Q_star, K_max, hwhm_l, hwhm_r


def normalize_curve(Q_grid, K_curve, scheme, params):
    """Return (x_normed, y_normed) for given normalisation scheme."""
    Q_star = params['Q_star']
    K_max = params['K_max']
    hwhm_avg = (params['hwhm_l'] + params['hwhm_r']) / 2.0 if params['hwhm_l'] and params['hwhm_r'] else None

    Q = np.array(Q_grid)
    K = np.array(K_curve)
    y = K / K_max

    if scheme == 'A':
        x = Q - Q_star
    elif scheme == 'B':
        x = (Q - Q_star) / Q_star
    elif scheme == 'C':
        if hwhm_avg is None:
            return None, None
        x = (Q - Q_star) / hwhm_avg
    else:
        raise ValueError(scheme)
    return x, y


def pairwise_residual(curves, scheme):
    """Compute pairwise mean-square residual on a common x-grid (linear
    interp). Returns (mean_residual, max_residual, common_x_grid, normed_dict)."""
    normed = {}
    for case_key, c in curves.items():
        Q_g = c['Q_grid']
        K_c = c['K_interior']
        Q_star, K_max, hl, hr = find_peak_and_hwhm(Q_g, K_c)
        params = {'Q_star': Q_star, 'K_max': K_max, 'hwhm_l': hl, 'hwhm_r': hr}
        x, y = normalize_curve(Q_g, K_c, scheme, params)
        if x is None:
            continue
        normed[case_key] = (x, y, params)

    if not normed:
        return None, None, None, normed

    # Common x range = intersection of all x ranges
    x_min = max(arr[0][0] for arr in normed.values())
    x_max = min(arr[0][-1] for arr in normed.values())
    x_common = np.linspace(x_min, x_max, 100)

    interp_curves = {}
    for k, (x, y, params) in normed.items():
        interp_curves[k] = np.interp(x_common, x, y)

    keys = list(interp_curves.keys())
    n = len(keys)
    if n < 2:
        return 0, 0, x_common, normed

    residuals = []
    for i in range(n):
        for j in range(i + 1, n):
            res = float(np.mean((interp_curves[keys[i]] - interp_curves[keys[j]]) ** 2))
            residuals.append(res)

    return float(np.mean(residuals)), float(np.max(residuals)), x_common, normed


def main():
    curves_path = RESULTS_DIR / "curves.json"
    if not curves_path.exists():
        print(f"# {curves_path} does not exist yet — scan still running?")
        return
    with open(curves_path) as f:
        curves = json.load(f)

    print(f"# Loaded {len(curves)} cases from {curves_path}")
    for case_key, c in curves.items():
        print(f"#   {case_key}: c={c['c']}, N={c['N']}, n={c['n']}, dim={c['block_dim']}")
    print()

    # Per-case summary
    print(f"# Per-case Interior-bond peak summary")
    print(f"{'case':<14} {'Q*':>8} {'|K|max':>8} {'HWHM-':>8} {'HWHM+':>8} {'asym':>6}")
    for case_key, c in curves.items():
        Q_star, K_max, hl, hr = find_peak_and_hwhm(c['Q_grid'], c['K_interior'])
        if hl and hr:
            print(f"{case_key:<14} {Q_star:>8.4f} {K_max:>8.5f} {hl:>8.4f} {hr:>8.4f} {hr/hl:>6.3f}")
        else:
            print(f"{case_key:<14} {Q_star:>8.4f} {K_max:>8.5f} {'n/a':>8} {'n/a':>8}")
    print()

    # Universality test for each scheme
    for scheme, name in [('A', 'absolute (Q − Q*)'),
                         ('B', 'relative ((Q − Q*)/Q*)'),
                         ('C', 'HWHM-scaled ((Q − Q*)/HWHM_avg)')]:
        mean_res, max_res, x_common, normed = pairwise_residual(curves, scheme)
        if mean_res is None:
            print(f"# Scheme {scheme} ({name}): not computable")
            continue
        print(f"# Scheme {scheme} ({name}):")
        print(f"#   mean pairwise residual: {mean_res:.6f}")
        print(f"#   max pairwise residual:  {max_res:.6f}")
        print()


if __name__ == "__main__":
    main()
