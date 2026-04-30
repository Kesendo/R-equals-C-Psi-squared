#!/usr/bin/env python3
"""Demo of the Cockpit predictive panel: cockpit_panel().

Adds two predictive structures to the existing 3-observable Decoherence
Cockpit (Purity, Concurrence, L1-coherence):
  1. Lebensader (skeleton + trace + rating)
  2. Cusp-pattern classifier

The existing Cockpit MEASURES; these two PREDICT, from the framework's
operator structure, what the trajectory will do before running.

Tested cases (4 bond-flipped softs at N=3, |+−+⟩, plus truly XX+YY):
"""
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.lebensader import cockpit_panel as lebensader_cockpit_panel


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),  # bond-flipped Z-free, drop=1
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),  # bond-flipped Z-free, drop=0
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),  # bond-flipped Z-content, drop=28
        ('XZ+ZX',       [('X', 'Z', J), ('Z', 'X', J)]),  # bond-flipped Z-content, drop=29
        ('XZ+XZ',       [('X', 'Z', J), ('X', 'Z', J)]),  # repeat, drop=40
    ]

    print(f"Cockpit predictive panel test")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}")
    print()

    print(f"{'case':<14s}  {'drop':>4s}  {'n_cross':>7s}  "
          f"{'tail (0<θ<5°)':>14s}  {'α':>7s}  {'pattern':<11s}  "
          f"{'mode type':<40s}  {'Lebensader rating':<35s}")
    print('-' * 145)

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        panel = lebensader_cockpit_panel(
            H, [GAMMA_DEPH] * N, rho_0, N,
            gamma_t1_l=[GAMMA_T1] * N, t_max=8.0, dt=0.005,
        )

        skel = panel['lebensader']['skeleton']
        tr = panel['lebensader']['trace']
        cusp = panel['cusp']
        rating = panel['lebensader']['rating']

        alpha_str = f"{tr['alpha_descent']:>+7.3f}" if tr['alpha_descent'] is not None else "    —  "
        print(f"{label:<14s}  {skel['drop']:>4d}  {tr['n_crossings']:>7d}  "
              f"{tr['tail_duration_sub5deg']:>13.4f}  {alpha_str}  "
              f"{cusp['pattern']:<11s}  "
              f"{cusp['mode_type']:<40s}  {rating:<35s}")

    print()
    print("Reading guide:")
    print("  Lebensader rating:")
    print("    'intact'    — bond-flipped Z-free pairs (XY+YX, IY+YI)")
    print("    'partial'   — skeleton holds OR trace persists, not both")
    print("    'collapsed' — Z-containing fragile cases under T1")
    print()
    print("  Cusp axes (orthogonal):")
    print("    pattern    'monotonic' (1 crossing) | 'heartbeat' (>1 crossing)")
    print("    mode type  'real-decay' | 'oscillatory'")
    print("               | 'steady-state + real-decay sub-mode' (settling)")
    print("               | 'steady-state + oscillatory sub-mode' (factorising)")


if __name__ == "__main__":
    main()
