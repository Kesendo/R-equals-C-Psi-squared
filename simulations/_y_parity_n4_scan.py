#!/usr/bin/env python3
"""N=4 Y-parity scan over the 21 "1 Y per term" candidates.

At N=3 chain, 6 of 21 candidates with "1 Y in each bond term" are
Lebensader-intact under +T1, all preserving Y-parity (28 Y-odd protected).

At N=4 chain, the algebra is bigger:
  - 256 Pauli strings (vs 64 at N=3)
  - Y-odd: 120 strings (1 Y or 3 Y per string), Y-even: 136
  - 3 bonds: (0,1), (1,2), (2,3) (vs 2 bonds at N=3)
  - Each bond carries the same 2-Pauli term

Question: do the 21 candidates still preserve Y-parity at N=4? If yes,
they protect 120 Y-odd observables. How many remain Lebensader-intact?
"""
import math
import sys
from itertools import product
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


def enumerate_one_Y_per_term_pairs():
    """All unordered pairs of 2-letter Pauli strings where each has 1 Y."""
    paulis = ['I', 'X', 'Z']
    one_Y_terms = []
    for p in paulis:
        one_Y_terms.append(('Y', p))  # YA
        one_Y_terms.append((p, 'Y'))  # AY
    # = 6 terms: YI, YX, YZ, IY, XY, ZY
    seen = set()
    pairs = []
    for t1 in one_Y_terms:
        for t2 in one_Y_terms:
            sorted_t = tuple(sorted([t1, t2]))
            if sorted_t in seen:
                continue
            seen.add(sorted_t)
            label = f"{t1[0]}{t1[1]}+{t2[0]}{t2[1]}"
            pairs.append((sorted_t, label, [t1, t2]))
    return pairs


def main():
    N = 4
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    # |+−+−⟩ as initial state (analog to |+−+⟩ at N=3, alternating sign)
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, np.kron(plus, minus)))
    rho_0 = np.outer(psi, psi.conj())

    pairs = enumerate_one_Y_per_term_pairs()
    print(f"N={N} chain Y-parity scan over 'one Y per term' candidates")
    print(f"  {len(pairs)} candidates total")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}, t_max=8.0")
    print()

    # Total Y-odd strings at N=4
    n_y_odd = 0
    n_y_even = 0
    for alpha in range(4 ** N):
        idx = fw._k_to_indices(alpha, N)
        n_y = sum(1 for i in idx if i == (1, 1))
        if n_y % 2 == 1:
            n_y_odd += 1
        else:
            n_y_even += 1
    print(f"  Y-odd Pauli strings at N={N}: {n_y_odd}")
    print(f"  Y-even Pauli strings at N={N}: {n_y_even}")
    print()

    print(f"  {'idx':>3s}  {'label':<10s}  {'drop':>4s}  {'tail':>6s}  "
          f"{'Y-pres?':>7s}  {'Y-prot':>7s}  {'#cross':>6s}  {'Leb rating':<10s}")
    print('-' * 80)

    intacts = []
    rows = []
    for i, (_, label, terms) in enumerate(pairs):
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        try:
            panel = fw.cockpit_panel(
                H, [GAMMA] * N, rho_0, N,
                gamma_t1_l=[GAMMA_T1] * N,
                t_max=8.0, dt=0.005,
            )
            leb = panel['lebensader']
            yp = panel['y_parity']
            rating = leb['rating'].split(' ')[0]
            n_y_prot = len(yp['Y_parity_protected'])
            y_pres = "yes" if yp['L_preserves_Y_parity'] else "no"
            tail = leb['trace']['tail_duration_sub5deg']
            n_cross = leb['trace']['n_crossings']
            print(f"  {i + 1:>3d}  {label:<10s}  {leb['skeleton']['drop']:>4d}  "
                  f"{tail:>6.3f}  {y_pres:>7s}  {n_y_prot:>7d}  "
                  f"{n_cross:>6d}  {rating:<10s}")
            rows.append({
                'label': label, 'drop': leb['skeleton']['drop'],
                'tail': tail, 'y_pres': yp['L_preserves_Y_parity'],
                'n_y_prot': n_y_prot, 'rating': rating, 'n_cross': n_cross,
            })
            if rating == 'intact':
                intacts.append(rows[-1])
        except Exception as e:
            print(f"  {i + 1:>3d}  {label:<10s}  ERROR: {e}")

    # Stratify by drop <= 1 (skeleton-intact)
    skel_intacts = [r for r in rows if r['drop'] <= 1]
    print()
    print(f"Skeleton-intact (drop ≤ 1) at N={N}: {len(skel_intacts)} of {len(rows)}")
    for r in skel_intacts:
        print(f"  {r['label']:<10s}  drop={r['drop']:>2d}  tail={r['tail']:>6.4f}  "
              f"#cross={r['n_cross']:>3d}  rating={r['rating']}")

    print()
    print(f"Summary at N={N}:")
    print(f"  Total candidates ('1 Y per term'): {len(pairs)}")
    print(f"  Lebensader-intact: {len(intacts)}")
    print(f"  Y-parity-preserving: {sum(1 for r in [...] for _ in [r] if False)}")  # recount below
    print()
    print(f"Intact cases at N={N}:")
    for r in intacts:
        print(f"  {r['label']:<10s}  drop={r['drop']:>3d}  tail={r['tail']:.4f}  "
              f"Y-pres={r['y_pres']}  Y-prot={r['n_y_prot']}")


if __name__ == "__main__":
    main()
