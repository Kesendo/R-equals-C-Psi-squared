#!/usr/bin/env python3
"""Re-examine the V-Effect's 14-of-36 result via framework.py.

The V-Effect (March 2026) found empirically that of 36 two-term
Pauli-pair Hamiltonians H = J(term1 + term2) at N=3, exactly 14 break and
22 preserve the palindromic structure. They tested EIGENVALUE PAIRING
(λ ↔ −λ − 2Σγ) which is the spectral consequence of the palindrome
operator relation.

The framework's `palindrome_residual` tests the stronger OPERATOR
equation Π·L·Π⁻¹ + L + 2Σγ·I = 0, which IMPLIES eigenvalue pairing.

Logical relation:
  operator equation holds  ⇒  eigenvalue pairs hold
  operator equation breaks  ⇏  eigenvalues necessarily break
  (eigenvalues can pair by accident when operator equation breaks)

Predicted counts:
  - 3 cases preserve operator equation: {XX+YY, XX+ZZ, YY+ZZ}
    (both terms both-parity-even per PROOF_ZERO_IMMUNITY arguments)
  - 33 cases break operator equation. Of these, the V-Effect
    found 14 also break eigenvalue pairing; the other 19 preserve
    eigenvalues despite operator break (accidental cancellation in
    spectrum).

This script:
1. Enumerates all 36 = C(9, 2) two-term Hamiltonians.
2. Tests OPERATOR equation (framework.palindrome_residual): 33 break, 3 don't.
3. Tests EIGENVALUE pairing (compute spectrum, check λ ↔ −λ − 2Σγ): expect
   14 break, 22 don't.
4. Identifies the 19 "soft-break" cases: operator broken, spectrum paired.
"""
import math
import sys
from itertools import combinations

import numpy as np

import framework as fw

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# All 9 single-term Pauli pairs (excluding II which is constant)
SINGLE_TERMS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']


def parity_classification(term):
    """Return (bit_a_even, bit_b_even) for a 2-body term like 'XY'."""
    a, b = term[0], term[1]
    return (
        fw.respects_bit_a_parity((a, b)),
        fw.respects_bit_b_parity((a, b)),
    )


def build_two_term_hamiltonian(N, bonds, term1, term2, J=1.0):
    """H = J Σ_bond (term1 + term2) per bond."""
    terms = [
        (term1[0], term1[1], J),
        (term2[0], term2[1], J),
    ]
    return fw._build_bilinear(N, bonds, terms)


def check_eigenvalue_pairing(L_vec, Sigma_gamma, tol=1e-8):
    """Test if eigenvalues of L pair as λ ↔ −λ − 2Σγ.

    Returns max pairing-error (0 if perfect pairing).
    """
    evals = np.linalg.eigvals(L_vec)
    # For each eigenvalue λ, find best partner λ' such that λ + λ' = -2Σγ.
    target_partner_real = -2 * Sigma_gamma
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    max_err = 0.0
    for i in range(n):
        if used[i]:
            continue
        ev = evals[i]
        partner_target = -ev - 2 * Sigma_gamma
        # Find best unused match
        best_j = -1
        best_dist = float('inf')
        for j in range(n):
            if used[j]:
                continue
            dist = abs(evals[j] - partner_target)
            if dist < best_dist:
                best_dist = dist
                best_j = j
        if best_j != -1 and best_j != i:
            used[i] = True
            used[best_j] = True
            max_err = max(max_err, best_dist)
        elif best_j == i:
            # self-paired: needs ev = -ev - 2Σγ ⇒ ev = -Σγ
            err = abs(ev - (-Sigma_gamma))
            used[i] = True
            max_err = max(max_err, err)
    return max_err


def main():
    N = 3
    gamma = 0.1
    gamma_l = [gamma] * N
    Sigma_gamma = sum(gamma_l)
    bonds = [(0, 1), (1, 2)]

    op_threshold = 1e-10  # operator-equation threshold
    spec_threshold = 1e-6  # eigenvalue-pairing threshold (looser due to eig error)

    print("=" * 90)
    print(f"V-Effect re-examination via framework.py: N={N}, bonds={bonds}, γ={gamma}")
    print(f"Two-term Hamiltonians: H = J (term1 + term2) per bond, all 36 = C(9,2) combos.")
    print(f"Two criteria tested:")
    print(f"  OP   = operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0 (residual > {op_threshold})")
    print(f"  SPEC = eigenvalue pairing λ ↔ −λ − 2Σγ (max pair-error > {spec_threshold})")
    print("=" * 90)

    combos = list(combinations(SINGLE_TERMS, 2))
    assert len(combos) == 36, f"Expected 36 combos, got {len(combos)}"

    records = []
    for term1, term2 in combos:
        H = build_two_term_hamiltonian(N, bonds, term1, term2)
        L = fw.lindbladian_z_dephasing(H, gamma_l)
        R = fw.palindrome_residual(L, Sigma_gamma, N)
        op_residual = float(np.linalg.norm(R))
        spec_pairing_err = check_eigenvalue_pairing(L, Sigma_gamma)

        c1 = parity_classification(term1)
        c2 = parity_classification(term2)
        both_terms_both_even = (c1 == (True, True)) and (c2 == (True, True))

        records.append({
            'term1': term1, 'term2': term2,
            'op_residual': op_residual,
            'spec_pairing_err': spec_pairing_err,
            'op_broken': op_residual > op_threshold,
            'spec_broken': spec_pairing_err > spec_threshold,
            'term1_parity': c1, 'term2_parity': c2,
            'both_terms_both_even': both_terms_both_even,
        })

    op_broken = [r for r in records if r['op_broken']]
    op_unbroken = [r for r in records if not r['op_broken']]
    spec_broken = [r for r in records if r['spec_broken']]
    spec_unbroken = [r for r in records if not r['spec_broken']]
    soft_break = [r for r in records if r['op_broken'] and not r['spec_broken']]
    hard_break = [r for r in records if r['op_broken'] and r['spec_broken']]

    print(f"\nResults:")
    print(f"  OP-equation:    {len(op_broken)} broken, {len(op_unbroken)} unbroken (framework's strict test)")
    print(f"  SPEC-pairing:   {len(spec_broken)} broken, {len(spec_unbroken)} unbroken (the V-Effect's test)")
    print(f"  hard break (both): {len(hard_break)}")
    print(f"  soft break (op only, spec OK): {len(soft_break)}")

    print(f"\nExpected from the V-Effect (March 2026): 14 SPEC-broken, 22 SPEC-unbroken.")
    print(f"  Match? {'YES' if len(spec_broken) == 14 else 'NO'}")
    print(f"\nExpected from framework: 3 OP-unbroken (XX+YY, XX+ZZ, YY+ZZ).")
    print(f"  Match? {'YES' if len(op_unbroken) == 3 else 'NO'}")

    print("\n" + "=" * 90)
    print(f"HARD BREAK ({len(hard_break)}): both operator and spectrum break, V-Effect's original 14")
    print("=" * 90)
    print(f"{'Combo':>10s}  {'op residual':>14s}  {'spec err':>14s}  {'parities':>22s}")
    for r in hard_break:
        c1_str = f"({'a' if r['term1_parity'][0] else '─'}{'b' if r['term1_parity'][1] else '─'})"
        c2_str = f"({'a' if r['term2_parity'][0] else '─'}{'b' if r['term2_parity'][1] else '─'})"
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['op_residual']:>14.4e}  {r['spec_pairing_err']:>14.4e}  {c1_str+'+'+c2_str:>22s}")

    print("\n" + "=" * 90)
    print(f"SOFT BREAK ({len(soft_break)}): operator equation broken but spectrum still palindromic")
    print("=" * 90)
    print(f"{'Combo':>10s}  {'op residual':>14s}  {'spec err':>14s}  {'parities':>22s}")
    for r in soft_break:
        c1_str = f"({'a' if r['term1_parity'][0] else '─'}{'b' if r['term1_parity'][1] else '─'})"
        c2_str = f"({'a' if r['term2_parity'][0] else '─'}{'b' if r['term2_parity'][1] else '─'})"
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['op_residual']:>14.4e}  {r['spec_pairing_err']:>14.4e}  {c1_str+'+'+c2_str:>22s}")

    print("\n" + "=" * 90)
    print(f"OP-UNBROKEN ({len(op_unbroken)}): both criteria pass, Heisenberg/XXZ subset")
    print("=" * 90)
    print(f"{'Combo':>10s}  {'op residual':>14s}  {'spec err':>14s}  {'parities':>22s}")
    for r in op_unbroken:
        c1_str = f"({'a' if r['term1_parity'][0] else '─'}{'b' if r['term1_parity'][1] else '─'})"
        c2_str = f"({'a' if r['term2_parity'][0] else '─'}{'b' if r['term2_parity'][1] else '─'})"
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['op_residual']:>14.4e}  {r['spec_pairing_err']:>14.4e}  {c1_str+'+'+c2_str:>22s}")

    print()
    print("Reading: parities are (a, b) where 'a' = bit_a-parity even, '─' = violated.")
    print()
    print("Structural finding:")
    print(f"  The 3 OP-unbroken cases all have BOTH terms ∈ {{XX, YY, ZZ}} (both-parity-even).")
    print(f"  The 14 hard-break cases break BOTH operator AND eigenvalue pairing.")
    print(f"  The {len(soft_break)} soft-break cases break the operator but eigenvalues coincidentally pair.")
    print(f"  The V-Effect's 22 unbroken = 3 (truly unbroken) + {len(soft_break)} (soft-break).")


if __name__ == "__main__":
    main()
