#!/usr/bin/env python3
"""Cockpit panel run over the full 120-enum at N=3.

For each of the 120 unordered Pauli-pair Hamiltonians at N=3, |+−+⟩,
γ_deph=0.1, γ_T1=0.01:
  - Compute the cockpit_panel
  - Categorise (truly / soft / hard) via spectrum-pairing + palindrome_residual
  - Cross-tabulate: (category, Lebensader rating, cusp pattern, mode type)

This validates the panel at scale and reveals the structural taxonomy
across all 2-Pauli Hamiltonians.
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
from framework.lebensader import cockpit_panel as lebensader_cockpit_panel


def spectrum_pair_max_err(L, sigma_gamma):
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_gamma
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        used[i] = True
        if best_j != i:
            used[best_j] = True
        max_err = max(max_err, float(dists[best_j]))
    return max_err


def classify(H, N, gamma):
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    if float(np.linalg.norm(M)) < 1e-10:
        return 'truly'
    if spectrum_pair_max_err(L, N * gamma) < 1e-6:
        return 'soft'
    return 'hard'


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

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    pairs = enumerate_pairs()
    print(f"Cockpit panel over the 120-enum at N={N}")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}, t_max=8.0, dt=0.005")
    print(f"  Running... (~6-8 min)")
    print()

    rows = []
    for i, (_, label, terms) in enumerate(pairs):
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        cat = classify(H, N, GAMMA)
        try:
            panel = lebensader_cockpit_panel(
                H, [GAMMA] * N, rho_0, N, gamma_t1_l=[GAMMA_T1] * N,
                t_max=8.0, dt=0.005,
            )
            skel = panel['lebensader']['skeleton']
            tr = panel['lebensader']['trace']
            cusp = panel['cusp']
            row = {
                'label': label, 'category': cat,
                'drop': skel['drop'],
                'n_protected_pure': skel['n_protected_pure'],
                'n_protected_t1': skel['n_protected_t1'],
                'n_crossings': tr['n_crossings'],
                'tail': tr['tail_duration_sub5deg'],
                'alpha': tr['alpha_descent'],
                'pattern': cusp['pattern'],
                'mode_type': cusp['mode_type'],
                'rating': panel['lebensader']['rating'],
            }
            rows.append(row)
        except Exception as e:
            rows.append({'label': label, 'category': cat, 'error': str(e)})
        if (i + 1) % 20 == 0:
            print(f"  ... {i + 1} / {len(pairs)} done")

    print()
    print(f"=== Cross-tabulations ===")
    print()

    # Category × Lebensader rating
    cat_rating = Counter((r['category'], r.get('rating', 'error')) for r in rows)
    print(f"Category × Lebensader rating:")
    print(f"  {'cat':<6s}  {'intact':>7s}  {'partial':>10s}  {'collapsed':>10s}")
    for cat in ['truly', 'soft', 'hard']:
        n_intact = sum(1 for r in rows if r['category'] == cat and 'intact' in r.get('rating', ''))
        n_partial = sum(1 for r in rows if r['category'] == cat and 'partial' in r.get('rating', ''))
        n_collapsed = sum(1 for r in rows if r['category'] == cat and 'collapsed' in r.get('rating', ''))
        print(f"  {cat:<6s}  {n_intact:>7d}  {n_partial:>10d}  {n_collapsed:>10d}")
    print()

    # Cusp pattern × mode type, all 120
    print(f"Cusp pattern × mode type (all 120):")
    pat_mode = Counter((r.get('pattern', 'error'), r.get('mode_type', 'error'))
                       for r in rows if 'error' not in r)
    print(f"  {'pattern':<11s}  {'mode_type':<42s}  {'count':>5s}")
    for (p, m), c in sorted(pat_mode.items(), key=lambda x: (-x[1], x[0])):
        print(f"  {p:<11s}  {m:<42s}  {c:>5d}")
    print()

    # Soft only — what's the breakdown?
    print(f"Soft category breakdown (46 cases):")
    soft_rows = [r for r in rows if r['category'] == 'soft']
    print(f"  By Lebensader rating:")
    soft_ratings = Counter(r.get('rating', 'error') for r in soft_rows)
    for rating, c in sorted(soft_ratings.items(), key=lambda x: -x[1]):
        print(f"    {rating:<50s}  {c:>3d}")
    print(f"  By cusp pattern × mode_type:")
    soft_cusps = Counter((r.get('pattern', 'err'), r.get('mode_type', 'err'))
                          for r in soft_rows)
    for (p, m), c in sorted(soft_cusps.items(), key=lambda x: -x[1]):
        print(f"    {p:<11s} | {m:<40s}  {c:>3d}")
    print()

    # Intact softs — list them
    intact = [r for r in rows if 'intact' in r.get('rating', '')]
    print(f"All 'intact' cases ({len(intact)}):")
    for r in intact:
        print(f"  {r['label']:<10s}  cat={r['category']:<6s}  drop={r['drop']:>3d}  "
              f"n_cross={r['n_crossings']:>3d}  tail={r['tail']:.4f}  "
              f"pattern={r['pattern']}, mode_type={r['mode_type']}")


if __name__ == "__main__":
    main()
