#!/usr/bin/env python3
"""F97 lens on the Kingston angle-steering data (2026-05-16).

F97 says the Mandelbrot main cardioid is the locus in complex-c where the
period-1 fixed point z* of z² + c has magnitude exactly b = 1/2. The Kingston
trajectory c(t) = CΨ_com(t) traces spirals in the complex-c plane. Two F97
questions on the data:

  (Q1) Does c(t) cross the cardioid boundary at any t? At what t and arg?
       The Kingston spirals can be checked against the cardioid curve, not
       just against the cusp at c = 1/4.

  (Q2) For each c(t) along the trajectory, what is the corresponding F97
       fixed point z*(c) = (1 - sqrt(1 - 4c)) / 2 and how does |z*(t)| evolve?
       On the cardioid |z*| = 1/2 invariant; off the cardioid it varies.
       Does the trajectory drive |z*| toward 1/2 (cardioid attractor) or
       away (cardioid repeller)?

  (Q3) Two crossings: the cardioid's |c| ranges from 1/4 at cusp (phi=0) to
       3/4 at tail (phi=pi). At phi = arccos(1/4) ≈ 75.5°, |c| = 1/2 on
       the cardioid. So the trajectory may cross |c| = 1/4 (cusp threshold,
       F95's b²) AND |c| = 1/2 (Half threshold, F97's b) at distinct times.
       Two F97-anchor-typed crossings per trajectory.

Reads `data/ibm_f95_angle_steering_may2026/*.json`. Reuses the canonical
CΨ_complex and trajectory-extraction helpers from `hardware_cpsi_cplane.py`.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np

from hardware_cpsi_cplane import extract_trajectory

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def f97_fixed_point(c, b=0.5):
    """Period-1 attracting fixed point of z² + c: z* = (1 - sqrt(1 - 4c)) / 2.

    Returns the inner branch (smaller |z*| for c inside the cardioid;
    z* on the cardioid for c on the boundary).
    """
    disc = b ** 2 - c  # 1/4 - c at b = 1/2
    return b - np.sqrt(disc)


def linear_interp(t_a, t_b, y_a, y_b, target):
    """Find t* such that linear interp between (t_a, y_a) and (t_b, y_b) hits target."""
    if y_b == y_a:
        return None
    return t_a + (target - y_a) / (y_b - y_a) * (t_b - t_a)


def analyze_pair(pair_data, pair_label, omega):
    """Compute F97-lens analysis on one pair's trajectory."""
    tr = extract_trajectory(pair_data)
    times = tr["t_us"]
    c_values_list = tr["cpsi_complex"]
    pair_info = tr["pair"] or pair_data.get("pair", {})

    print(f"--- {pair_label} (Ω = {omega} rad/μs) ---")
    print(f"  qubits: {pair_info.get('qubits', '?')}")
    if "gamma_per_us" in pair_info:
        print(f"  γ = {pair_info['gamma_per_us']:.3e} μs⁻¹")
    print()

    print(f"  {'t (μs)':>8} {'|CΨ_com|':>10} {'arg(c) (°)':>12} "
          f"{'|z*(c)|':>10} {'gap to b=½':>14} {'on cardioid?':>14}")
    print("  " + "-" * 78)

    c_values = []
    for t, c in zip(times, c_values_list):
        z_star = f97_fixed_point(c)
        c_values.append((t, c, z_star))
        gap = abs(z_star) - 0.5  # On the cardioid <=> |z*| ≈ 1/2 (b)
        on_card = "yes" if abs(gap) < 0.01 else f"({gap:+.3f})"
        print(f"  {t:>8.3f} {abs(c):>10.5f} {np.degrees(np.angle(c)):>+12.2f} "
              f"{abs(z_star):>10.5f} {gap:>+14.6f} {on_card:>14}")
    print()

    # Q3: find crossings of |c| = 1/4 (cusp) and |c| = 1/2 (Half on cardioid)
    print("  F97-anchor crossings (linear interp between adjacent delays):")
    targets = [("|c| = 1/4 (Quarter cusp)", 0.25), ("|c| = 1/2 (Half cardioid)", 0.5)]
    for label, target in targets:
        found = False
        for i in range(len(c_values) - 1):
            t_a, c_a, _ = c_values[i]
            t_b, c_b, _ = c_values[i + 1]
            ma, mb = abs(c_a), abs(c_b)
            if (ma - target) * (mb - target) < 0:
                t_cross = linear_interp(t_a, t_b, ma, mb, target)
                arg_a, arg_b = np.angle(c_a), np.angle(c_b)
                # Unwrap across the ±π discontinuity
                if abs(arg_a - arg_b) > np.pi:
                    if arg_b < arg_a:
                        arg_b += 2 * np.pi
                    else:
                        arg_a += 2 * np.pi
                arg_cross = arg_a + (target - ma) / (mb - ma) * (arg_b - arg_a)
                arg_cross_deg = np.degrees(arg_cross)
                while arg_cross_deg > 180:
                    arg_cross_deg -= 360
                while arg_cross_deg <= -180:
                    arg_cross_deg += 360
                print(f"    {label}: t = {t_cross:.4f} μs, "
                      f"arg = {arg_cross_deg:+.2f}°")
                found = True
                break
        if not found:
            mags = [abs(c) for _, c, _ in c_values]
            print(f"    {label}: no crossing in delay window "
                  f"(|c| ranges {min(mags):.4f} to {max(mags):.4f})")
    print()


def main():
    repo_data = Path("data/ibm_f95_angle_steering_may2026")
    files = sorted(repo_data.glob("cusp_complex_phase_hardware_ibm_kingston_omega*.json"))
    if not files:
        print(f"No data files found in {repo_data}")
        return

    print("=" * 80)
    print("F97 Lens on Kingston Angle-Steering Data (2026-05-16)")
    print("=" * 80)
    print()

    for f in files:
        with open(f) as fh:
            data = json.load(fh)
        omega = data["omega_per_us"]
        print(f"### File: {f.name}")
        print(f"### Ω = {omega} rad/μs, phi_0 = {data['phi_0_deg']}°")
        print()
        for pair_label, pair_data in data["pair_runs"].items():
            analyze_pair(pair_data, pair_label, omega)
        print("=" * 80)
        print()


if __name__ == "__main__":
    main()
