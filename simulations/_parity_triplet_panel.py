#!/usr/bin/env python3
"""3-bit parity signature panel: bit_a, bit_b, bit_a XOR bit_b parities.

Together with Y-parity (= bit_a · bit_b parity at single-site, NOT
multiplicative), these give a complete classification of (H, ρ_0) by
which Z₂ symmetries L preserves.

For each test case, report:
  bit_a parity:  preserved? + protected count
  bit_b parity:  preserved? + protected count
  bit_a XOR b:   preserved? + protected count
  Y-parity:      preserved? + protected count
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


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XY+YZ',       [('X', 'Y', J), ('Y', 'Z', J)]),
        ('XX+YZ',       [('X', 'X', J), ('Y', 'Z', J)]),
        ('XX+ZZ',       [('X', 'X', J), ('Z', 'Z', J)]),  # 0-Y case
        ('IY+IZ',       [('I', 'Y', J), ('I', 'Z', J)]),  # mixed Y/Z
    ]

    print(f"3-bit parity signature panel, N={N}, |+−+⟩, +T1")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}")
    print()

    print(f"  {'case':<14s}  ", end='')
    for kind in ['bit_a', 'bit_b', 'bit_aXb', 'Y-parity']:
        print(f"{kind:>9s} {'#prot':>5s}  ", end='')
    print(f"  {'signature':>12s}")
    print('-' * 100)

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)

        results = {}
        for kind in ['bit_a', 'bit_b', 'bit_a_xor_b']:
            r = fw.parity_panel(H, [GAMMA] * N, rho_0, N, kind,
                                  gamma_t1_l=[GAMMA_T1] * N)
            results[kind] = r
        # Y-parity (existing)
        yp = fw.y_parity_panel(H, [GAMMA] * N, rho_0, N,
                                 gamma_t1_l=[GAMMA_T1] * N)
        results['y_parity'] = yp

        # Build 4-bit signature: 1 if preserved, 0 if not
        sig = ""
        for kind in ['bit_a', 'bit_b', 'bit_a_xor_b']:
            sig += "1" if results[kind]['L_preserves_parity'] else "0"
        sig += "1" if results['y_parity']['L_preserves_Y_parity'] else "0"

        print(f"  {label:<14s}  ", end='')
        for kind, k_label in [('bit_a', 'bit_a'),
                                ('bit_b', 'bit_b'),
                                ('bit_a_xor_b', 'bit_aXb')]:
            r = results[kind]
            preserved = "yes" if r['L_preserves_parity'] else "no"
            n_prot = len(r['protected'])
            print(f"{preserved:>9s} {n_prot:>5d}  ", end='')
        # Y-parity
        yp_pres = "yes" if results['y_parity']['L_preserves_Y_parity'] else "no"
        n_yp = len(results['y_parity']['Y_parity_protected'])
        print(f"{yp_pres:>9s} {n_yp:>5d}  ", end='')
        print(f"  {sig:>12s}")

    print()
    print("Reading guide:")
    print("  4-bit signature: bit_a, bit_b, bit_a XOR b, Y-parity")
    print("  preservation order. '1' means L preserves that parity.")
    print()
    print("  Z-dephasing alone preserves all 3 multiplicative parities.")
    print("  T1 preserves bit_a always (σ⁻ has pure bit_a=1).")
    print("  T1 mixes bit_b and bit_a XOR b (σ⁻ has mixed bit_b content).")
    print("  Combined L preserves a parity iff [H,·] also preserves it,")
    print("  i.e., iff H is parity-even in its Pauli decomposition.")


if __name__ == "__main__":
    main()
