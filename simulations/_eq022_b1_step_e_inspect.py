#!/usr/bin/env python3
"""_eq022_b1_step_e_inspect.py — inspect the resonance-shape curves.

Pulls the curves saved by `_eq022_b1_step_e_resonance_shape.py`, prints
K at key Q points (peak ± fractions of Q_peak), and tests the relative-Q
universal-shape hypothesis with explicit point-by-point comparisons.
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


def find_peak_parabolic(Q_grid, K_curve):
    Q = np.array(Q_grid)
    K = np.array(K_curve)
    i_max = int(np.argmax(K))
    K_max = float(K[i_max])
    if 0 < i_max < len(Q) - 1:
        x = Q[i_max - 1: i_max + 2]
        y = K[i_max - 1: i_max + 2]
        coefs = np.polyfit(x, y, 2)
        if coefs[0] < 0:
            Q_star = -coefs[1] / (2 * coefs[0])
            K_max = coefs[2] - coefs[1] ** 2 / (4 * coefs[0])
            return float(Q_star), float(K_max)
    return float(Q[i_max]), float(K_max)


def normalise_relative(Q_grid, K_curve):
    Q_star, K_max = find_peak_parabolic(Q_grid, K_curve)
    Q = np.array(Q_grid)
    K = np.array(K_curve)
    x = (Q - Q_star) / Q_star
    y = K / K_max
    return x, y, Q_star, K_max


def main():
    curves_path = RESULTS_DIR / "curves.json"
    with open(curves_path) as f:
        curves = json.load(f)

    # Sample relative-Q points to compare curves
    sample_x = np.array([-0.6, -0.4, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2,
                          0.4, 0.6, 0.8, 1.0, 1.5])

    # Two bond classes
    for bond_kind, label in [('K_interior', 'INTERIOR'), ('K_endpoint', 'ENDPOINT')]:
        print(f"\n## {label} relative-Q shape comparison")
        print(f"  Each row = K/|K|max evaluated at (Q-Q*)/Q* = sample_x point.")
        print(f"  If shape universal, columns should be near-constant.")
        print()
        header = "  case          Q*   |K|max  " + "  ".join(f"{x:+.2f}" for x in sample_x)
        print(header)

        for case_key, c in curves.items():
            Q_grid = c['Q_grid']
            K_curve = c[bond_kind]
            x, y, Q_star, K_max = normalise_relative(Q_grid, K_curve)
            sampled = np.interp(sample_x, x, y)
            row = f"  {case_key:<10} {Q_star:6.3f} {K_max:7.4f}  "
            row += "  ".join(f"{v:+.3f}" for v in sampled)
            print(row)

    print()
    print("\n## Interior shape: HWHM and tail behaviour from interpolated curves")
    print("  Looking for x where y = 0.5 (HWHM in relative-Q units).")
    for case_key, c in curves.items():
        x, y, Q_star, K_max = normalise_relative(c['Q_grid'], c['K_interior'])
        # Find x where y = 0.5 on each side
        i_max = int(np.argmax(y))
        # Left side
        x_half_left = None
        for i in range(i_max, -1, -1):
            if y[i] < 0.5:
                # Linear interp
                x0, x1 = x[i], x[i + 1]
                y0, y1 = y[i], y[i + 1]
                x_half_left = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
                break
        # Right side
        x_half_right = None
        for i in range(i_max, len(x)):
            if y[i] < 0.5:
                x0, x1 = x[i - 1], x[i]
                y0, y1 = y[i - 1], y[i]
                x_half_right = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
                break
        # Tail at x=1.0 and x=2.0
        y_at_1 = float(np.interp(1.0, x, y))
        y_at_2 = float(np.interp(2.0, x, y)) if x[-1] >= 2.0 else None

        l_str = f"{x_half_left:+.3f}" if x_half_left is not None else "  n/a "
        r_str = f"{x_half_right:+.3f}" if x_half_right is not None else "  n/a "
        y2_str = f"{y_at_2:.3f}" if y_at_2 is not None else " n/a"
        print(f"  {case_key:<10}  HWHM_left = {l_str}  HWHM_right = {r_str}  y(x=1) = {y_at_1:.3f}  y(x=2) = {y2_str}")


if __name__ == "__main__":
    main()
