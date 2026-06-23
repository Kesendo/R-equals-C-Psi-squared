#!/usr/bin/env python3
"""EQ-030 follow-up: is soft's T1-robustness universal across all 46 soft cases,
or specific to XY+YX?

The XY+YX result was striking: 30 protected under pure Z-dephasing, 29 under
γ_T1 = 0.1·γ_deph. Only one observable leaked. That looked nothing like
truly's collapse (32 → 1) or hard's collapse (32 → 0).

This script enumerates ALL 46 soft Hamiltonians at N=3, runs the T1 sweep
on each, and reports the distribution of "drops" (how many protected
observables leak when T1 is turned on at 10% of γ_dephasing).

If most soft cases drop by 0-1, soft's T1-robustness is a universal
category property. If the distribution is wide, the robustness is
specific to certain soft sub-orbits.
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
    """Return 'truly' | 'soft' | 'hard'."""
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    op_norm = float(np.linalg.norm(M))
    if op_norm < 1e-10:
        return 'truly'
    spec_err = spectrum_pair_max_err(L, N * gamma)
    return 'soft' if spec_err < 1e-6 else 'hard'


def count_protected(H, gamma_l, gamma_t1_l, rho_0, N,
                    threshold=1e-9, cluster_tol=1e-8):
    """Count Π-protected Pauli observables under L = L_H + L_Z + L_T1.

    Uses pure Z-dephasing if all gamma_t1 = 0, otherwise lindbladian_z_plus_t1.
    """
    if any(g != 0 for g in gamma_t1_l):
        L = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = fw.lindbladian_z_dephasing(H, gamma_l)

    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)
    rho_pauli = fw.pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli

    n_eig = len(evals)
    used = np.zeros(n_eig, dtype=bool)
    clusters = []
    for i in range(n_eig):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n_eig):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    n_protected = 0
    for alpha in range(1, 4 ** N):
        max_S = 0.0
        for cl in clusters:
            S = sum(V[alpha, k] * c[k] for k in cl)
            max_S = max(max_S, abs(S))
        if max_S < threshold:
            n_protected += 1
    return n_protected


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
    GAMMA_DEPH = 0.1
    T1_RATIO = 0.1
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    print(f"Soft T1-robustness universality at N={N}, |+−+⟩, "
          f"γ_deph={GAMMA_DEPH}, γ_T1/γ_deph={T1_RATIO}")
    print(f"Enumerating all 120 unordered Pauli-pair Hamiltonians,")
    print(f"filtering to soft, then measuring protected-count drop under T1.")
    print()

    pairs = enumerate_pairs()
    soft_results = []
    truly_results = []
    hard_results = []

    for sorted_terms, label, terms in pairs:
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        cat = classify(H, N, GAMMA_DEPH)
        n_pure = count_protected(
            H, [GAMMA_DEPH] * N, [0.0] * N, rho_0, N
        )
        n_t1 = count_protected(
            H, [GAMMA_DEPH] * N, [T1_RATIO * GAMMA_DEPH] * N, rho_0, N
        )
        record = {
            'label': label,
            'terms': terms,
            'n_pure': n_pure,
            'n_t1': n_t1,
            'drop': n_pure - n_t1,
        }
        if cat == 'soft':
            soft_results.append(record)
        elif cat == 'truly':
            truly_results.append(record)
        else:
            hard_results.append(record)

    print(f"Counts: truly={len(truly_results)}, soft={len(soft_results)}, "
          f"hard={len(hard_results)} (total {len(pairs)})")
    print()

    print("=" * 78)
    print("SOFT category — full enumeration")
    print("=" * 78)
    print(f"  {'label':<12s}  {'pure-Z':>8s}  {'+ T1':>8s}  {'drop':>6s}")
    print(f"  {'-' * 12}  {'-' * 8}  {'-' * 8}  {'-' * 6}")
    soft_results.sort(key=lambda r: (r['drop'], -r['n_pure']))
    for r in soft_results:
        print(f"  {r['label']:<12s}  {r['n_pure']:>8d}  {r['n_t1']:>8d}  {r['drop']:>6d}")

    print()
    drop_hist = {}
    for r in soft_results:
        drop_hist[r['drop']] = drop_hist.get(r['drop'], 0) + 1
    print("Drop distribution (soft only):")
    for drop in sorted(drop_hist.keys()):
        bar = "█" * drop_hist[drop]
        print(f"  drop={drop:>3d}: {drop_hist[drop]:>3d} cases  {bar}")

    n_robust = sum(1 for r in soft_results if r['drop'] <= 1)
    n_total = len(soft_results)
    print()
    print(f"Soft cases with drop ≤ 1 (T1-robust): {n_robust}/{n_total} "
          f"({100 * n_robust / n_total:.1f}%)")

    print()
    print("=" * 78)
    print("Reference: TRULY and HARD drop distributions")
    print("=" * 78)
    for cat_name, results in [('truly', truly_results), ('hard', hard_results)]:
        print(f"\n{cat_name} ({len(results)} cases):")
        drop_hist_ref = {}
        for r in results:
            drop_hist_ref[r['drop']] = drop_hist_ref.get(r['drop'], 0) + 1
        for drop in sorted(drop_hist_ref.keys()):
            bar = "█" * drop_hist_ref[drop]
            print(f"  drop={drop:>3d}: {drop_hist_ref[drop]:>3d} cases  {bar}")

    print()
    print("Verdict:")
    if n_robust == n_total:
        print(f"  Soft's T1-robustness is UNIVERSAL: all {n_total} soft cases "
              f"keep ≥ {min(r['n_t1'] for r in soft_results)} protected under T1.")
    elif n_robust >= 0.8 * n_total:
        print(f"  Soft's T1-robustness is MOSTLY universal "
              f"({n_robust}/{n_total}); a minority leaks more than 1.")
    else:
        print(f"  Soft's T1-robustness is SUB-CLUSTER-SPECIFIC "
              f"(only {n_robust}/{n_total} robust).")


if __name__ == "__main__":
    main()
