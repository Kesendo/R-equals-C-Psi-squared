#!/usr/bin/env python3
"""Full Cockpit validation: Lebensader + Cusp + Chiral + Y-parity panels
on the same 6+2 cases at N=3, with two initial states.

Initial 1: |+−+⟩ — the existing Snapshot-D state.
  ρ_0 has Y-parity 0 (no Y components), but is NOT a K-eigenstate
  (K|+−+⟩ = |+++⟩). Y-parity protection active for "1 Y per term" H's;
  K protection inactive.

Initial 2: (|+−+⟩+|+++⟩)/√2 — K-symmetric superposition.
  ρ_0 IS a K-eigenstate (w_+ = 1, w_- = 0). Both Y-parity (some) AND
  K (full) protection can act, depending on H's K-status.
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

    psi_xneel = np.kron(plus, np.kron(minus, plus))
    psi_K_sym = (psi_xneel + np.kron(plus, np.kron(plus, plus))) / math.sqrt(2)
    psi_K_sym = psi_K_sym / np.linalg.norm(psi_K_sym)

    states = {
        '|+−+⟩': np.outer(psi_xneel, psi_xneel.conj()),
        '(|+−+⟩+|+++⟩)/√2': np.outer(psi_K_sym, psi_K_sym.conj()),
    }

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XY+YZ',       [('X', 'Y', J), ('Y', 'Z', J)]),
        ('XX+YZ',       [('X', 'X', J), ('Y', 'Z', J)]),  # hard, fragile (Y-mixed)
    ]

    for state_label, rho_0 in states.items():
        print(f"=" * 100)
        print(f"Initial state: {state_label}")
        print(f"=" * 100)
        print()

        print(f"  {'case':<14s}  {'Leb':<10s}  "
              f"{'Π drop':>6s}  "
              f"{'Cusp pat':<11s}  "
              f"{'K_stat':<8s}  {'K-prot':>7s}  "
              f"{'Y-pres?':>7s}  {'Y-prot':>7s}")
        print('-' * 100)

        for label, terms in cases:
            H = fw._build_bilinear(N, bonds, terms)
            panel = lebensader_cockpit_panel(
                H, [GAMMA_DEPH] * N, rho_0, N,
                gamma_t1_l=[GAMMA_T1] * N, t_max=8.0, dt=0.005,
            )

            leb = panel['lebensader']
            rating_short = leb['rating'].split(' ')[0]  # 'intact' / 'partial' / 'collapsed'
            cusp_pat = panel['cusp']['pattern']
            chiral = panel['chiral']
            yp = panel['y_parity']

            n_K_dyn = len(chiral['chiral_protected_dynamic'])
            y_preserved = "yes" if yp['L_preserves_Y_parity'] else "no"
            n_Y_prot = len(yp['Y_parity_protected'])
            print(f"  {label:<14s}  {rating_short:<10s}  "
                  f"{leb['skeleton']['drop']:>6d}  "
                  f"{cusp_pat:<11s}  "
                  f"{chiral['K_status']:<8s}  {n_K_dyn:>7d}  "
                  f"{y_preserved:>7s}  {n_Y_prot:>7d}")
        print()

    print("=" * 100)
    print("Reading the four panels:")
    print("=" * 100)
    print()
    print("Lebensader (Π skeleton + θ trace):")
    print("  intact = drop ≤ 1 AND tail > 0.05 — the framework's primary")
    print("  T1-robustness signal.")
    print()
    print("Cusp pattern:")
    print("  monotonic (1 crossing) | heartbeat (>1 crossing) | never crosses")
    print()
    print("Chiral (K_full = ⊗_{odd i} Z_i):")
    print("  K_stat = K-even / K-odd / K-mixed (commutator status with H)")
    print("  K-prot = # K-protected Pauli observables (requires K-symmetric L")
    print("           AND ρ_0 in K-eigenstate)")
    print()
    print("Y-parity (bit_a · bit_b grading):")
    print("  Y-pres? = does L block-diagonalise in Y-parity?")
    print("  Y-prot = # Y-parity-protected Pauli observables (requires both)")


if __name__ == "__main__":
    main()
