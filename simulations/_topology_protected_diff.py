#!/usr/bin/env python3
"""Diff the protected Pauli sets across topologies for XY+YX and XY+YZ.

From the prior topology comparison:
  XY+YX  drop:  chain=1  ring=4  star=4   ← 3 extra leaks vs chain
  XY+YZ  drop:  chain=0  ring=0  star=4   ← star-specific leak

This script:
  1. For each (H, topology, channel), compute the full protected set
  2. Show the Y-parity composition (Y-odd vs Y-even) of each set
  3. Diff: what does chain protect that ring/star doesn't (and vice versa)?
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


def get_protected_set(H, gamma_l, gamma_t1_l, rho_0, N,
                      threshold=1e-9, cluster_tol=1e-8):
    """Return set of Pauli labels (e.g., 'XYZ') protected under L."""
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

    protected = set()
    for alpha in range(1, 4 ** N):
        max_S = 0.0
        for cl in clusters:
            S = sum(V[alpha, k] * c[k] for k in cl)
            max_S = max(max_S, abs(S))
        if max_S < threshold:
            label = ''.join(fw.PAULI_LABELS[idx]
                             for idx in fw._k_to_indices(alpha, N))
            protected.add(label)
    return protected


def y_parity(label):
    return label.count('Y') % 2


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    topologies = {
        'chain': [(0, 1), (1, 2)],
        'ring':  [(0, 1), (1, 2), (2, 0)],
        'star':  [(0, 1), (0, 2)],
    }

    cases = [
        ('XY+YX', [('X', 'Y', J), ('Y', 'X', J)]),
        ('XY+YZ', [('X', 'Y', J), ('Y', 'Z', J)]),
    ]

    for label, terms in cases:
        print(f"=" * 90)
        print(f"Hamiltonian: {label}")
        print(f"=" * 90)
        print()

        sets_pure = {}
        sets_t1 = {}
        for top, bonds in topologies.items():
            H = fw._build_bilinear(N, bonds, terms)
            sets_pure[top] = get_protected_set(
                H, [GAMMA] * N, [0.0] * N, rho_0, N
            )
            sets_t1[top] = get_protected_set(
                H, [GAMMA] * N, [GAMMA_T1] * N, rho_0, N
            )

        print(f"  {'Topology':<8s}  {'pure-Z |total':>14s}  {'pure-Z |Y-even':>14s}  "
              f"{'pure-Z |Y-odd':>13s}  {'+T1 |total':>10s}  "
              f"{'+T1 |Y-even':>11s}  {'+T1 |Y-odd':>10s}  {'drop':>4s}")
        for top in topologies:
            sp = sets_pure[top]
            st1 = sets_t1[top]
            sp_even = sum(1 for x in sp if y_parity(x) == 0)
            sp_odd = sum(1 for x in sp if y_parity(x) == 1)
            st_even = sum(1 for x in st1 if y_parity(x) == 0)
            st_odd = sum(1 for x in st1 if y_parity(x) == 1)
            drop = len(sp) - len(st1)
            print(f"  {top:<8s}  {len(sp):>14d}  {sp_even:>14d}  {sp_odd:>13d}  "
                  f"{len(st1):>10d}  {st_even:>11d}  {st_odd:>10d}  {drop:>4d}")
        print()

        # Diff under +T1
        print(f"  +T1 protected sets diff:")
        chain_t1 = sets_t1['chain']
        ring_t1 = sets_t1['ring']
        star_t1 = sets_t1['star']

        # in chain but not in ring
        c_minus_r = chain_t1 - ring_t1
        r_minus_c = ring_t1 - chain_t1
        print(f"    chain \\ ring   (in chain, NOT in ring):")
        print(f"      {sorted(c_minus_r)}  ({len(c_minus_r)} entries)")
        print(f"    ring \\ chain   (in ring, NOT in chain):")
        print(f"      {sorted(r_minus_c)}  ({len(r_minus_c)} entries)")
        print()

        c_minus_s = chain_t1 - star_t1
        s_minus_c = star_t1 - chain_t1
        print(f"    chain \\ star   (in chain, NOT in star):")
        print(f"      {sorted(c_minus_s)}  ({len(c_minus_s)} entries)")
        print(f"    star \\ chain   (in star, NOT in chain):")
        print(f"      {sorted(s_minus_c)}  ({len(s_minus_c)} entries)")
        print()

        # Common across all 3
        common = chain_t1 & ring_t1 & star_t1
        print(f"  Common to all 3 topologies under +T1: {len(common)} entries")
        # Y-classification of common
        c_even = sum(1 for x in common if y_parity(x) == 0)
        c_odd = sum(1 for x in common if y_parity(x) == 1)
        print(f"    Y-odd: {c_odd}, Y-even: {c_even}")
        print()


if __name__ == "__main__":
    main()
