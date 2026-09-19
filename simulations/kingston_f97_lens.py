#!/usr/bin/env python3
"""Offline radial/phase reader for the Kingston 2026-05-16 archive.

The stable filename is historical. The reconstructed `Cpsi_complex` is a
density-matrix radial diagnostic, not the F97 parameter `c=z-z²` and not an
F95 root. The reader reports finite magnitude/argument rows and interpolation
times at selected radial levels. Only recurrence `c=+1/4` is the period-one
cardioid cusp; a complex diagnostic with magnitude 1/4 is a different object.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np

from hardware_cpsi_cplane import extract_trajectory

def linear_interp(t_a, t_b, y_a, y_b, target):
    """Find t* such that linear interp between (t_a, y_a) and (t_b, y_b) hits target."""
    if y_b == y_a:
        return None
    return t_a + (target - y_a) / (y_b - y_a) * (t_b - t_a)


def analyze_pair(pair_data, pair_label, omega):
    """Print the saved trajectory's radial diagnostic and phase."""
    tr = extract_trajectory(pair_data)
    times = tr["t_us"]
    diagnostic_values = tr["cpsi_complex"]
    pair_info = tr["pair"] or pair_data.get("pair", {})

    print(f"--- {pair_label} (Ω = {omega} rad/μs) ---")
    print(f"  qubits: {pair_info.get('qubits', '?')}")
    if "gamma_per_us" in pair_info:
        print(f"  γ = {pair_info['gamma_per_us']:.3e} μs⁻¹")
    print()

    print(f"  {'t (μs)':>8} {'|CΨ_com|':>10} {'arg(CΨ_com) (°)':>18}")
    print("  " + "-" * 42)

    rows = []
    for t, value in zip(times, diagnostic_values):
        rows.append((t, value))
        print(f"  {t:>8.3f} {abs(value):>10.5f} "
              f"{np.degrees(np.angle(value)):>+18.2f}")
    print()

    print("  Radial diagnostic interpolations between adjacent delays:")
    targets = [("|CΨ_com| = 1/4", 0.25), ("|CΨ_com| = 1/2", 0.5)]
    for label, target in targets:
        found = False
        for i in range(len(rows) - 1):
            t_a, value_a = rows[i]
            t_b, value_b = rows[i + 1]
            ma, mb = abs(value_a), abs(value_b)
            if (ma - target) * (mb - target) < 0:
                t_cross = linear_interp(t_a, t_b, ma, mb, target)
                arg_a, arg_b = np.angle(value_a), np.angle(value_b)
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
            mags = [abs(value) for _, value in rows]
            print(f"    {label}: no crossing in delay window "
                  f"(|CΨ_com| ranges {min(mags):.4f} to {max(mags):.4f})")
    print()


def main():
    if sys.platform == "win32":
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass

    repo_data = (Path(__file__).resolve().parents[1] / "data"
                 / "ibm_f95_angle_steering_may2026")
    files = sorted(repo_data.glob("cusp_complex_phase_hardware_ibm_kingston_omega*.json"))
    if not files:
        print(f"No data files found in {repo_data}")
        return

    print("=" * 80)
    print("Kingston radial/phase archive reader (2026-05-16)")
    print("=" * 80)
    print()

    for f in files:
        with open(f, encoding="utf-8") as fh:
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
