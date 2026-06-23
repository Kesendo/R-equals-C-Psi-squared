#!/usr/bin/env python3
"""EQ-026: Soft sub-cluster symmetry at N=5 — applying today's V-Effect layers.

The 46 soft cases at N=5 distribute into discrete protected-count clusters:
  992 (1): YZ+ZY
  862 (2): IY+IY, YI+YI
  781 (1): IY+YI
  772 (6): XY+XY, XZ+XZ, YX+YX, YZ+YZ, ZX+ZX, ZY+ZY
  512 (18): ...

Question: which structural feature orders these sub-clusters?

This script checks each of the 46 soft cases against today's combinatorial
features:
  - bit_a-parity preservation by L (Section 14)
  - bit_b-parity preservation
  - Y-parity preservation (Section 13)
  - bond-flip
  - Z-align
  - 1-BPE / 0-BPE / both-single classification
  - protected count via pi_protected_observables
  - saturation vs theoretical max (4^N - 2^N = 992 at N=5)

Cross-tabulate against the empirical sub-cluster grouping.
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


GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
BPE = {('I', 'I'), ('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}


def n_I(T):
    return T.count('I')


def active_letters(T):
    return [c for c in T if c != 'I']


def classify_combinatorial_pair(T1, T2):
    n1, n2 = n_I(T1), n_I(T2)
    if n1 == 1 and n2 == 1:
        letters = active_letters(T1) + active_letters(T2)
        unique_letters = set(letters)
        if 'Z' in unique_letters: return 'hard'
        if unique_letters == {'X'}: return 'truly'
        return 'soft'
    if (n1 == 1) != (n2 == 1):
        single = T1 if n1 == 1 else T2
        double = T2 if n1 == 1 else T1
        single_letter = active_letters(single)[0]
        if single_letter == 'Z': return 'hard'
        if single_letter == 'X':
            if double in BPE: return 'truly'
            if double in {('Y','Z'),('Z','Y')}: return 'soft'
            return 'hard'
        if single_letter == 'Y':
            if double in BPE or double in {('X','Z'),('Z','X')}: return 'soft'
            return 'hard'
    if n1 == 0 and n2 == 0:
        in_BPE_1 = T1 in BPE; in_BPE_2 = T2 in BPE
        if in_BPE_1 and in_BPE_2: return 'truly'
        if T1 == T2:
            if T1 in BPE: return 'truly'
            return 'soft'
        if in_BPE_1 ^ in_BPE_2:
            bpe = T1 if in_BPE_1 else T2
            non_bpe = T2 if in_BPE_1 else T1
            a = bpe[0]
            if a == 'Z': return 'soft'
            partner = 'Y' if a == 'X' else 'X'
            if (non_bpe[0] == a and non_bpe[1] == partner) or (non_bpe[1] == a and non_bpe[0] == partner):
                return 'hard'
            return 'soft'
        if T1[::-1] == T2: return 'soft'
        if (T1[0] == 'Z' and T2[0] == 'Z') or (T1[1] == 'Z' and T2[1] == 'Z'):
            return 'soft'
        return 'hard'
    return 'soft'


def alternating_xneel(N):
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = plus.copy()
    for k in range(1, N):
        psi = np.kron(psi, minus if k % 2 == 1 else plus)
    return psi


def compute_n_protected(H, gamma_l, gamma_t1_l, rho_0, N):
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
        if used[i]: continue
        cl = [i]; used[i] = True
        for j in range(i+1, n_eig):
            if not used[j] and abs(evals[j] - evals[i]) < 1e-8:
                cl.append(j); used[j] = True
        clusters.append(cl)
    n_protected = 0
    for alpha in range(1, 4 ** N):
        max_S = max(abs(sum(V[alpha, k] * c[k] for k in cl)) for cl in clusters)
        if max_S < 1e-9:
            n_protected += 1
    return n_protected


def main():
    N = 5
    bonds = [(i, i+1) for i in range(N-1)]
    psi = alternating_xneel(N)
    rho_0 = np.outer(psi, psi.conj())

    paulis = ['I','X','Y','Z']
    seen = set()
    pairs = []
    for t1 in product(paulis, repeat=2):
        for t2 in product(paulis, repeat=2):
            if t1 == ('I','I') or t2 == ('I','I'): continue
            sorted_t = tuple(sorted([t1, t2]))
            if sorted_t in seen: continue
            seen.add(sorted_t)
            pairs.append(sorted_t)

    # Filter to 46 soft cases via combinatorial classifier
    soft_pairs = [p for p in pairs if classify_combinatorial_pair(p[0], p[1]) == 'soft']
    print(f"N={N}, |+−+−+⟩, computing protected counts for {len(soft_pairs)} soft cases (this may take a few min)...")
    print()

    max_pure = 4**N - 2**N
    print(f"  Theoretical max (pure ρ_0): {max_pure}")
    print()

    rows = []
    for i, pair in enumerate(soft_pairs):
        t1, t2 = pair
        terms = [(t1[0], t1[1], J), (t2[0], t2[1], J)]
        H = fw._build_bilinear(N, bonds, terms)

        # n_protected under +T1
        n_t1 = compute_n_protected(H, [GAMMA]*N, [GAMMA_T1]*N, rho_0, N)
        # n_protected under pure-Z
        n_pure = compute_n_protected(H, [GAMMA]*N, [0.0]*N, rho_0, N)

        # Today's structural features
        n1, n2 = n_I(t1), n_I(t2)
        is_bond_flipped = (t1[::-1] == t2)
        is_z_aligned = ((t1[0] == 'Z' and t2[0] == 'Z') or
                         (t1[1] == 'Z' and t2[1] == 'Z'))
        in_bpe_count = int(t1 in BPE) + int(t2 in BPE)
        # Layer signature
        if n1 == 0 and n2 == 0:
            cls = '0-BPE' if in_bpe_count == 0 else '1-BPE' if in_bpe_count == 1 else '2-BPE'
        elif n1 == 1 and n2 == 1:
            cls = 'both-single'
        else:
            cls = 'mixed'

        label = f"{t1[0]}{t1[1]}+{t2[0]}{t2[1]}"
        rows.append({
            'label': label, 't1': t1, 't2': t2,
            'n_pure': n_pure, 'n_t1': n_t1,
            'drop': n_pure - n_t1,
            'class': cls, 'bf': is_bond_flipped, 'z_align': is_z_aligned,
        })
        if (i+1) % 10 == 0:
            print(f"  ... {i+1}/{len(soft_pairs)} done")

    # Group by n_pure (sub-cluster identity)
    print()
    print(f"=" * 100)
    print(f"Sub-clusters by n_protected (pure-Z) at N={N}:")
    print(f"=" * 100)
    print()
    by_npure = {}
    for r in rows:
        by_npure.setdefault(r['n_pure'], []).append(r)
    for n_p in sorted(by_npure.keys(), reverse=True):
        rs = by_npure[n_p]
        # Structural features summary
        classes = sorted(set(r['class'] for r in rs))
        n_bf = sum(1 for r in rs if r['bf'])
        n_za = sum(1 for r in rs if r['z_align'])
        members = [r['label'] for r in rs]
        sat = n_p / max_pure * 100
        print(f"  n_protected = {n_p:>4d} ({sat:>5.1f}% of max), {len(rs)} cases, "
              f"classes={classes}, bond-flipped={n_bf}, Z-aligned={n_za}")
        for label in members[:5]:
            print(f"    {label}")
        if len(members) > 5:
            print(f"    ... ({len(members) - 5} more)")
        print()


if __name__ == "__main__":
    main()
