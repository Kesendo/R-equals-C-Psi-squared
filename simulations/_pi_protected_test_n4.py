#!/usr/bin/env python3
"""EQ-026 attack: extend the 14/19/3 trichotomy to N=4 using pi_protected_observables.

For each of the 36 ordered two-term Pauli-pair Hamiltonians at N=4
(σ_a σ_b + σ_c σ_d on each of the three bonds (0,1), (1,2), (2,3)),
classify by:
  1. Operator residual ‖M‖ via framework.palindrome_residual
  2. Eigenvalue spectrum pairing under λ ↔ −λ − 2Σγ
  3. Number of Π-protected Pauli-string observables on |+−+−⟩

If the trichotomy structure (truly / soft / hard) holds at N=4, the
counts should cluster into three groups: high-protected (truly),
intermediate (soft), low (hard).

The init state |+−+−⟩ generalises Snapshot D's N=3 |+−+⟩.
"""
import math
import sys
from itertools import combinations_with_replacement, product
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
    """V-Effect-style: max |λ_i + λ_j + 2Σγ| over best pairing of L's eigenvalues."""
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


def main():
    N = 4
    GAMMA = 0.1
    J = 1.0
    Σγ = N * GAMMA

    # |+−+−⟩ initial state
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = plus.copy()
    for sign in [-1, +1, -1]:  # signs[1..N-1]
        psi = np.kron(psi, plus if sign > 0 else minus)
    rho_0 = np.outer(psi, psi.conj())

    bonds = [(i, i + 1) for i in range(N - 1)]
    paulis = ['I', 'X', 'Y', 'Z']

    # All 36 two-term Pauli pairs (unordered terms)
    all_pairs = []
    for term1 in product(paulis, repeat=2):
        for term2 in product(paulis, repeat=2):
            if term1 == ('I', 'I') or term2 == ('I', 'I'):
                continue  # exclude identity-only terms
            # Take unordered pairs (term1, term2)
            sorted_terms = tuple(sorted([term1, term2]))
            if sorted_terms in [t for (t1, t2), _ in all_pairs for t in [(t1, t2)]]:
                # Already seen
                continue
            all_pairs.append((sorted_terms, [term1, term2]))
    # Deduplicate by sorted-pair-tuple
    seen = set()
    pairs = []
    for sorted_terms, terms in all_pairs:
        if sorted_terms in seen:
            continue
        seen.add(sorted_terms)
        pairs.append((sorted_terms, terms))

    print(f"Testing pi_protected_observables at N={N}, γ={GAMMA}, |+−+−⟩")
    print(f"{len(pairs)} unordered two-term Pauli-pair Hamiltonians.")
    print()
    print(f"  {'H = J(...)':>12s}  {'‖M‖_op':>10s}  {'spec_err':>10s}  "
          f"{'protected':>10s}  {'active':>8s}  {'verdict':>10s}")
    print(f"  {'-' * 12}  {'-' * 10}  {'-' * 10}  {'-' * 10}  {'-' * 8}  {'-' * 10}")

    results = []
    for sorted_terms, terms in pairs:
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        L = fw.lindbladian_z_dephasing(H, [GAMMA] * N)

        # Operator residual
        M = fw.palindrome_residual(L, Σγ, N)
        op_norm = float(np.linalg.norm(M))

        # Spectrum pairing
        spec_err = spectrum_pair_max_err(L, Σγ)

        # Π-protected observables
        result = fw.pi_protected_observables(H, [GAMMA] * N, rho_0, N)
        n_prot = len(result['protected'])
        n_act = len(result['active'])

        # Classify
        spec_ok = spec_err < 1e-6
        op_ok = op_norm < 1e-10
        if op_ok:
            verdict = "truly"
        elif spec_ok:
            verdict = "soft"
        else:
            verdict = "hard"

        label = f"{terms[0][0]}{terms[0][1]}+{terms[1][0]}{terms[1][1]}"
        results.append({
            'pair': sorted_terms, 'label': label,
            'op_norm': op_norm, 'spec_err': spec_err,
            'n_protected': n_prot, 'n_active': n_act,
            'verdict': verdict,
        })

    # Sort by verdict, then by n_protected
    order = {'truly': 0, 'soft': 1, 'hard': 2}
    results.sort(key=lambda r: (order[r['verdict']], -r['n_protected']))
    for r in results:
        print(f"  {r['label']:>12s}  {r['op_norm']:>10.2e}  {r['spec_err']:>10.2e}  "
              f"{r['n_protected']:>10d}  {r['n_active']:>8d}  {r['verdict']:>10s}")

    # Summary
    n_truly = sum(1 for r in results if r['verdict'] == 'truly')
    n_soft = sum(1 for r in results if r['verdict'] == 'soft')
    n_hard = sum(1 for r in results if r['verdict'] == 'hard')

    print()
    print(f"Summary at N={N}:")
    print(f"  truly:  {n_truly:>3d}  (operator equation exact, spectrum paired)")
    print(f"  soft:   {n_soft:>3d}  (operator equation broken, spectrum still paired)")
    print(f"  hard:   {n_hard:>3d}  (both broken)")
    print(f"  total:  {n_truly + n_soft + n_hard:>3d}")
    print()
    print(f"Protected-observable count distribution:")
    n_prot_truly = [r['n_protected'] for r in results if r['verdict'] == 'truly']
    n_prot_soft = [r['n_protected'] for r in results if r['verdict'] == 'soft']
    n_prot_hard = [r['n_protected'] for r in results if r['verdict'] == 'hard']
    if n_prot_truly:
        print(f"  truly:  {min(n_prot_truly):>3d} - {max(n_prot_truly):>3d}  (range)")
    if n_prot_soft:
        print(f"  soft:   {min(n_prot_soft):>3d} - {max(n_prot_soft):>3d}  (range)")
    if n_prot_hard:
        print(f"  hard:   {min(n_prot_hard):>3d} - {max(n_prot_hard):>3d}  (range)")

    print()
    print("Comparison to N=3 (V_EFFECT_FINE_STRUCTURE.md):")
    print(f"  N=3:  3 truly,  19 soft,  14 hard")
    print(f"  N=4:  {n_truly} truly, {n_soft} soft, {n_hard} hard")


if __name__ == "__main__":
    main()
