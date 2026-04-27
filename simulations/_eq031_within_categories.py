#!/usr/bin/env python3
"""EQ-031 follow-up: within-category fine structure across N1 → N2.

EQ-031 showed all 120 Pauli-pair Hamiltonians keep their truly/soft/hard
category across N=3 → N=4. This script asks: within each category, does
the *quantitative* fine structure shift, and is the operator-norm
amplification factor universal across N or N-pair-specific?

For each Hamiltonian, extract (op_norm, spec_err, n_protected) at N1 and
N2 on the |+−+−...⟩ initial state. Tabulate per-category Spearman rank
correlation and op_norm ratio statistics.

Usage:
    python _eq031_within_categories.py            # default N1=3, N2=4
    python _eq031_within_categories.py 4 5        # custom pair
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


GAMMA = 0.1
J = 1.0


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
        if best_j != i:
            used[i] = True
            used[best_j] = True
        else:
            used[i] = True
        max_err = max(max_err, float(dists[best_j]))
    return max_err


def alternating_state(N):
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = plus.copy()
    for k in range(1, N):
        psi = np.kron(psi, minus if k % 2 == 1 else plus)
    return np.outer(psi, psi.conj())


def measure_at_N(N, terms):
    bonds = [(i, i + 1) for i in range(N - 1)]
    bilinear = [(t[0], t[1], J) for t in terms]
    H = fw._build_bilinear(N, bonds, bilinear)
    L = fw.lindbladian_z_dephasing(H, [GAMMA] * N)
    Sigma_gamma = N * GAMMA

    M = fw.palindrome_residual(L, Sigma_gamma, N)
    op_norm = float(np.linalg.norm(M))
    spec_err = spectrum_pair_max_err(L, Sigma_gamma)

    rho_0 = alternating_state(N)
    result = fw.pi_protected_observables(H, [GAMMA] * N, rho_0, N)
    n_prot = len(result['protected'])
    n_act = len(result['active'])

    op_ok = op_norm < 1e-10
    spec_ok = spec_err < 1e-6
    if op_ok:
        verdict = "truly"
    elif spec_ok:
        verdict = "soft"
    else:
        verdict = "hard"

    return {
        'op_norm': op_norm,
        'spec_err': spec_err,
        'n_protected': n_prot,
        'n_active': n_act,
        'verdict': verdict,
    }


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


def spearman(xs, ys):
    """Spearman rank correlation; returns NaN if degenerate."""
    n = len(xs)
    if n < 2:
        return float('nan')
    rx = np.argsort(np.argsort(xs))
    ry = np.argsort(np.argsort(ys))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt((rx * rx).sum() * (ry * ry).sum())
    if denom == 0.0:
        return float('nan')
    return float((rx * ry).sum() / denom)


def main():
    N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 4

    pairs = enumerate_pairs()
    print(f"EQ-031 within-category analysis ({len(pairs)} Pauli-pair Hamiltonians)")
    print(f"  N={N1} → N={N2}")
    print(f"  init |+−+...⟩, γ={GAMMA}, J={J}")
    print()

    rows = []
    for sorted_terms, label, terms in pairs:
        m3 = measure_at_N(N1, terms)
        m4 = measure_at_N(N2, terms)
        rows.append({
            'label': label,
            'verdict': m3['verdict'],
            'op_norm_3': m3['op_norm'], 'op_norm_4': m4['op_norm'],
            'spec_err_3': m3['spec_err'], 'spec_err_4': m4['spec_err'],
            'n_prot_3': m3['n_protected'], 'n_prot_4': m4['n_protected'],
            'n_act_3': m3['n_active'], 'n_act_4': m4['n_active'],
        })
        # Sanity check: category invariance per EQ-031 / 079c7ce
        assert m3['verdict'] == m4['verdict'], (
            f"category drift by {label}: N={N1}→{m3['verdict']}, N={N2}→{m4['verdict']}"
        )

    # Per-category analysis
    for cat in ['truly', 'soft', 'hard']:
        sub = [r for r in rows if r['verdict'] == cat]
        if not sub:
            continue
        print(f"=== {cat.upper()} ({len(sub)} Hamiltonians) ===")

        op3 = np.array([r['op_norm_3'] for r in sub])
        op4 = np.array([r['op_norm_4'] for r in sub])
        sp3 = np.array([r['spec_err_3'] for r in sub])
        sp4 = np.array([r['spec_err_4'] for r in sub])
        np3 = np.array([r['n_prot_3'] for r in sub])
        np4 = np.array([r['n_prot_4'] for r in sub])
        na3 = np.array([r['n_act_3'] for r in sub])
        na4 = np.array([r['n_act_4'] for r in sub])

        def stats(arr, label, fmt='.3e'):
            print(f"  {label:>12s}:  min={arr.min():{fmt}}  "
                  f"max={arr.max():{fmt}}  "
                  f"mean={arr.mean():{fmt}}  std={arr.std():{fmt}}")

        if cat in ('soft', 'hard'):
            stats(op3, f'op_norm N={N1}')
            stats(op4, f'op_norm N={N2}')
            r = spearman(op3, op4)
            ratio = op4 / np.where(op3 > 0, op3, np.nan)
            print(f"  op_norm Spearman(N={N1}, N={N2}): {r:+.4f}")
            if np.isfinite(ratio).any():
                ratio_sq = ratio ** 2
                print(f"  op_norm ratio N={N2}/N={N1}: "
                      f"min={np.nanmin(ratio):.3f}, max={np.nanmax(ratio):.3f}, "
                      f"mean={np.nanmean(ratio):.3f}, std={np.nanstd(ratio):.3f}")
                # Compare against framework closed-form predictions
                try:
                    pred_main = fw.palindrome_residual_norm_ratio_squared(N1, N2, 'main')
                    pred_sb = fw.palindrome_residual_norm_ratio_squared(N1, N2, 'single_body')
                    print(f"  ratio² range: min={np.nanmin(ratio_sq):.4f}, "
                          f"max={np.nanmax(ratio_sq):.4f}")
                    print(f"  framework prediction: main={pred_main:.4f}, "
                          f"single_body={pred_sb:.4f}")
                except (NotImplementedError, ValueError):
                    pass

        if cat == 'hard':
            stats(sp3, f'spec_err N={N1}')
            stats(sp4, f'spec_err N={N2}')
            r = spearman(sp3, sp4)
            print(f"  spec_err Spearman(N={N1}, N={N2}): {r:+.4f}")

        stats(np3.astype(float), f'n_prot N={N1}', fmt='.1f')
        stats(np4.astype(float), f'n_prot N={N2}', fmt='.1f')
        stats(na3.astype(float), f'n_act  N={N1}', fmt='.1f')
        stats(na4.astype(float), f'n_act  N={N2}', fmt='.1f')
        r = spearman(np3, np4)
        print(f"  n_prot Spearman(N={N1}, N={N2}): {r:+.4f}")
        # Pauli-basis fraction: 4^N - 1 = total non-identity Paulis
        frac3 = np3 / (4**N1 - 1)
        frac4 = np4 / (4**N2 - 1)
        print(f"  n_prot / (4^N − 1):  N={N1} mean={frac3.mean():.3f}, "
              f"N={N2} mean={frac4.mean():.3f}  "
              f"(diff={frac4.mean() - frac3.mean():+.3f})")
        print()

    # Spotlight: any Hamiltonians where op_norm increases sharply or n_prot collapses
    print("=== Spotlight: largest fine-structure changes ===")
    print()
    print("Top-5 op_norm growth (within soft+hard):")
    sub = [r for r in rows if r['verdict'] in ('soft', 'hard') and r['op_norm_3'] > 0]
    for r in sorted(sub, key=lambda x: x['op_norm_4'] / x['op_norm_3'], reverse=True)[:5]:
        ratio = r['op_norm_4'] / r['op_norm_3']
        print(f"  {r['label']:>10s} ({r['verdict']:<5s})  "
              f"op_norm: {r['op_norm_3']:.3e} → {r['op_norm_4']:.3e}  "
              f"(ratio {ratio:.2f})")
    print()
    print("Top-5 op_norm shrink:")
    for r in sorted(sub, key=lambda x: x['op_norm_4'] / x['op_norm_3'])[:5]:
        ratio = r['op_norm_4'] / r['op_norm_3']
        print(f"  {r['label']:>10s} ({r['verdict']:<5s})  "
              f"op_norm: {r['op_norm_3']:.3e} → {r['op_norm_4']:.3e}  "
              f"(ratio {ratio:.2f})")


if __name__ == "__main__":
    main()
