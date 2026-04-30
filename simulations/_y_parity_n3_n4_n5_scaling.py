#!/usr/bin/env python3
"""Y-parity scan over N = 3, 4, 5 chain — scaling of skeleton-intact /
trace-dwell with system size, plus threshold-design data.

For each N ∈ {3, 4, 5}:
  - Initial state |+−+−...⟩ (alternating)
  - 21 'one Y per term' candidate Hamiltonians on the chain
  - Cockpit panel: drop, tail, n_crossings, Y-parity preservation

Tabulate: skeleton-intact count (drop ≤ 1), tail distribution, max tail,
median tail. Use the data to design an N-scaled tail threshold for the
Lebensader 'intact' rating.
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
from framework.lebensader import cockpit_panel as lebensader_cockpit_panel


def enumerate_one_Y_per_term_pairs():
    paulis = ['I', 'X', 'Z']
    one_Y_terms = []
    for p in paulis:
        one_Y_terms.append(('Y', p))
        one_Y_terms.append((p, 'Y'))
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


def alternating_xneel(N):
    """|+−+−...⟩ with N qubits, alternating + and -."""
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = plus.copy()
    for k in range(1, N):
        psi = np.kron(psi, minus if k % 2 == 1 else plus)
    return psi


def n_y_odd_at_N(N):
    """Number of Y-parity-odd Pauli strings on N qubits."""
    n = 0
    for alpha in range(4 ** N):
        idx = fw._k_to_indices(alpha, N)
        n_y = sum(1 for i in idx if i == (1, 1))
        if n_y % 2 == 1:
            n += 1
    return n


def main():
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0

    pairs = enumerate_one_Y_per_term_pairs()
    print(f"Y-parity scaling scan over N = 3, 4, 5 chain")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}, t_max=8.0")
    print(f"  21 'one Y per term' candidates")
    print()

    summary_per_N = {}

    for N in [3, 4, 5]:
        bonds = [(i, i + 1) for i in range(N - 1)]
        psi = alternating_xneel(N)
        rho_0 = np.outer(psi, psi.conj())

        n_y_odd = n_y_odd_at_N(N)
        print(f"--- N = {N} ---")
        print(f"  Y-odd Pauli strings: {n_y_odd}")
        print()

        rows = []
        for i, (_, label, terms) in enumerate(pairs):
            bilinear = [(t[0], t[1], J) for t in terms]
            H = fw._build_bilinear(N, bonds, bilinear)
            try:
                panel = lebensader_cockpit_panel(
                    H, [GAMMA] * N, rho_0, N,
                    gamma_t1_l=[GAMMA_T1] * N,
                    t_max=8.0, dt=0.005,
                )
                leb = panel['lebensader']
                yp = panel['y_parity']
                tail = leb['trace']['tail_duration_sub5deg']
                drop = leb['skeleton']['drop']
                rows.append({
                    'label': label, 'drop': drop, 'tail': tail,
                    'n_y_prot': len(yp['Y_parity_protected']),
                    'y_pres': yp['L_preserves_Y_parity'],
                })
            except Exception as e:
                print(f"  ERROR for {label}: {e}")
                continue

        # Summarize
        skel_intact = [r for r in rows if r['drop'] <= 1]
        y_pres_all = [r for r in rows if r['y_pres']]
        max_tail = max(r['tail'] for r in rows) if rows else 0.0
        median_tail = float(np.median([r['tail'] for r in rows])) if rows else 0.0
        max_tail_skel_intact = max((r['tail'] for r in skel_intact), default=0.0)
        median_tail_skel_intact = float(np.median([r['tail'] for r in skel_intact])) if skel_intact else 0.0

        summary_per_N[N] = {
            'rows': rows,
            'n_y_odd': n_y_odd,
            'n_skel_intact': len(skel_intact),
            'n_y_pres': len(y_pres_all),
            'max_tail': max_tail,
            'median_tail': median_tail,
            'max_tail_skel_intact': max_tail_skel_intact,
            'median_tail_skel_intact': median_tail_skel_intact,
        }
        print(f"  Skeleton-intact (drop ≤ 1):  {len(skel_intact)} of {len(rows)}")
        print(f"  Y-parity-preserving:         {len(y_pres_all)} of {len(rows)}")
        print(f"  Max tail (all):              {max_tail:.4f}")
        print(f"  Median tail (all):           {median_tail:.4f}")
        print(f"  Max tail (skel-intact):      {max_tail_skel_intact:.4f}")
        print(f"  Median tail (skel-intact):   {median_tail_skel_intact:.4f}")
        print()

    print("=" * 80)
    print("Scaling table:")
    print("=" * 80)
    print()
    print(f"  {'N':>2s}  {'#Y-odd':>6s}  {'#skel-int':>10s}  {'#Y-pres':>8s}  "
          f"{'max τ':>8s}  {'med τ':>8s}  {'max τ_skel':>11s}  {'med τ_skel':>11s}")
    for N in [3, 4, 5]:
        s = summary_per_N[N]
        print(f"  {N:>2d}  {s['n_y_odd']:>6d}  {s['n_skel_intact']:>10d}  "
              f"{s['n_y_pres']:>8d}  {s['max_tail']:>8.4f}  {s['median_tail']:>8.4f}  "
              f"{s['max_tail_skel_intact']:>11.4f}  "
              f"{s['median_tail_skel_intact']:>11.4f}")

    print()
    print("Threshold design:")
    print("  Current 'intact' rule: drop ≤ 1 AND tail > 0.05")
    print()
    print("  Looking at max τ_skel (the largest tail among skeleton-intact cases)")
    print("  per N as a sensible scaling target:")
    for N in [3, 4, 5]:
        s = summary_per_N[N]
        print(f"    N={N}: max τ_skel = {s['max_tail_skel_intact']:.4f}")


if __name__ == "__main__":
    main()
