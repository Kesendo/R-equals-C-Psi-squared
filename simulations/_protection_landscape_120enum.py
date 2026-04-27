#!/usr/bin/env python3
"""Protection landscape over all 120 Pauli-pair Hamiltonians at N=3.

For each H ∈ 120-enum:
  - Classify (truly / soft / hard) via spectrum + palindrome residual
  - Run recommender across the full catalog of initial states
  - Record max n_protected and which ρ_0 achieves it
  - Compare to theoretical max (4^N − 2^N = 56 at N=3)

Cross-tabulate: how does V-Effect's category structure correlate with
the protection landscape? Do truly cases always saturate? Are some hard
cases unexpectedly close to saturation?
"""
import math
import sys
from collections import Counter
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


def classify(H, N, gamma):
    """V-Effect category: truly (M=0), soft (spectrum-paired, M≠0), hard."""
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    if float(np.linalg.norm(M)) < 1e-10:
        return 'truly'
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    sigma_g = N * gamma
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_g
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        used[i] = True
        if best_j != i:
            used[best_j] = True
        max_err = max(max_err, float(dists[best_j]))
    return 'soft' if max_err < 1e-6 else 'hard'


def enumerate_pairs():
    paulis = ['I', 'X', 'Y', 'Z']
    seen = set()
    pairs = []
    for term1 in product(paulis, repeat=2):
        for term2 in product(paulis, repeat=2):
            if term1 == ('I', 'I') or term2 == ('I', 'I'):
                continue
            sorted_terms = tuple(sorted([term1, term2]))
            if sorted_terms in seen:
                continue
            seen.add(sorted_terms)
            label = f"{term1[0]}{term1[1]}+{term2[0]}{term2[1]}"
            pairs.append((sorted_terms, label, [term1, term2]))
    return pairs


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]
    max_pure = fw.theoretical_max_protected_pure(N)

    pairs = enumerate_pairs()
    print(f"Protection landscape over 120-enum at N={N}")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}")
    print(f"  Theoretical max (pure ρ_0): {max_pure}")
    print(f"  Catalog size: {2 ** N + 6} = 14 candidates per H")
    print(f"  Running... (estimated 4-5 min)")
    print()

    rows = []
    for i, (_, label, terms) in enumerate(pairs):
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        cat = classify(H, N, GAMMA)
        try:
            rec = fw.recommend_initial_state(
                H, [GAMMA] * N, N, gamma_t1_l=[GAMMA_T1] * N, top_k=3,
            )
            best_n = rec['pure_protection_best']['n_protected']
            best_label = rec['pure_protection_best']['label']
            quantum_best = rec.get('quantum_best')
            quantum_n = quantum_best['n_protected'] if quantum_best else 0
            quantum_label = quantum_best['label'] if quantum_best else '—'
            saturation = best_n / max_pure
        except Exception as e:
            print(f"  ERROR for {label}: {e}")
            continue
        rows.append({
            'label': label, 'category': cat,
            'best_n': best_n, 'best_label': best_label,
            'quantum_n': quantum_n, 'quantum_label': quantum_label,
            'saturation': saturation,
        })
        if (i + 1) % 30 == 0:
            print(f"  ... {i + 1}/{len(pairs)} done")
    print()

    # ===== Cross-tabulation =====
    print(f"=" * 100)
    print(f"Saturation distribution by V-Effect category")
    print(f"=" * 100)
    print()

    print(f"  {'cat':<6s}  {'n':>3s}  {'sat 100%':>9s}  {'sat 75-100':>11s}  "
          f"{'sat 50-75':>10s}  {'sat 25-50':>10s}  {'sat <25':>8s}  "
          f"{'mean sat':>9s}")
    print('-' * 100)
    for cat in ['truly', 'soft', 'hard']:
        cat_rows = [r for r in rows if r['category'] == cat]
        if not cat_rows:
            continue
        n_full = sum(1 for r in cat_rows if r['saturation'] >= 0.99)
        n_75 = sum(1 for r in cat_rows if 0.75 <= r['saturation'] < 0.99)
        n_50 = sum(1 for r in cat_rows if 0.50 <= r['saturation'] < 0.75)
        n_25 = sum(1 for r in cat_rows if 0.25 <= r['saturation'] < 0.50)
        n_less = sum(1 for r in cat_rows if r['saturation'] < 0.25)
        mean_sat = float(np.mean([r['saturation'] for r in cat_rows]))
        print(f"  {cat:<6s}  {len(cat_rows):>3d}  {n_full:>9d}  {n_75:>11d}  "
              f"{n_50:>10d}  {n_25:>10d}  {n_less:>8d}  {mean_sat:>8.0%}")
    print()

    # ===== Outliers =====
    print(f"=" * 100)
    print(f"Outliers — hards/softs that saturate (overlap with truly's regime)")
    print(f"=" * 100)
    print()

    saturating_non_truly = [r for r in rows if r['category'] != 'truly'
                              and r['saturation'] >= 0.99]
    print(f"  Non-truly cases reaching ≥99% saturation: {len(saturating_non_truly)}")
    for r in sorted(saturating_non_truly, key=lambda x: x['label']):
        print(f"    {r['label']:<10s}  cat={r['category']:<6s}  "
              f"n_prot={r['best_n']}  via {r['best_label']}")
    print()

    # ===== Trulys that DON'T saturate =====
    truly_not_full = [r for r in rows if r['category'] == 'truly'
                       and r['saturation'] < 0.99]
    print(f"  Truly cases NOT reaching 99% saturation: {len(truly_not_full)}")
    for r in sorted(truly_not_full, key=lambda x: -x['saturation'])[:10]:
        print(f"    {r['label']:<10s}  saturation={r['saturation']:.0%}  "
              f"n_prot={r['best_n']}  via {r['best_label']}")
    print()

    # ===== Summary =====
    print(f"=" * 100)
    print(f"Quantum-best (excluding classical states) saturation by category")
    print(f"=" * 100)
    print()
    for cat in ['truly', 'soft', 'hard']:
        cat_rows = [r for r in rows if r['category'] == cat]
        if not cat_rows:
            continue
        quantum_sats = [r['quantum_n'] / max_pure for r in cat_rows
                         if r['quantum_n'] > 0]
        if quantum_sats:
            mean_q = float(np.mean(quantum_sats))
            max_q = max(quantum_sats)
            print(f"  {cat:<6s}: mean quantum-sat = {mean_q:.0%}, max = {max_q:.0%}")


if __name__ == "__main__":
    main()
