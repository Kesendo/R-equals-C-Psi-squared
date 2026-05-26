"""Polarity probe #9: search for k_max(N), the Pauli-rank boundary of balance.

Post-probes-7+8 conjecture (2026-05-26):
    There exists k_max(N) such that if c is supported on <= k_max Pauli strings
    (with arbitrary complex coefficients), then asymmetry = 0; if c supported
    on > k_max Pauli strings, asymmetry != 0.

Empirical anchors so far:
    k_pauli = 1 (probe D): BALANCE
    k_pauli = 2 (probe F, X + iZ): BALANCE
    k_pauli = 4^N (probe 7 random): BROKEN

This probe sweeps k_pauli geometrically and finds the boundary.

For each N in {2, 3}:
    For k_pauli in {1, 2, 4, 8, ..., 4^N}:
        For trial in range(5):
            Pick k_pauli random Pauli strings from the 4^N total
            Random complex coefficients a_alpha
            c = sum a_alpha * P_alpha
            Test asymmetry

We also do a parallel sweep over a SINGLE site (sticking to Pauli strings on
site 0 only with identity elsewhere) to test "site spread" vs "Pauli rank"
separately.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import product

from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform, site_op

PAULI_2X2 = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string(letters):
    """Return the 2^N x 2^N matrix for the Pauli string given as a list of letters."""
    op = PAULI_2X2[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI_2X2[l])
    return op


def enumerate_pauli_strings(N):
    """All 4^N Pauli strings as letter tuples."""
    return list(product(['I', 'X', 'Y', 'Z'], repeat=N))


def build_L_with_c(H, c, gamma, sigma_target=None):
    """Standard Lindblad with a single jump operator c."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    L = L + gamma * np.kron(c, c.conj())
    return L


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def test_c_at_pauli_rank(N, k_pauli, all_paulis, rng, fixed_site_mask=None):
    """Build c = sum_{alpha in S} a_alpha * P_alpha for random subset S of size k_pauli.

    If fixed_site_mask is given (a tuple like ('Z', 'I', 'I') means: only consider
    Pauli strings that are supported on site 0 only -- letters at sites where mask
    has 'I' must be 'I'), restrict the sampling to that support pattern.
    """
    d = 2 ** N

    if fixed_site_mask is not None:
        # Restrict to Pauli strings consistent with fixed identity sites
        eligible = [p for p in all_paulis
                    if all(p[i] == 'I' for i, m in enumerate(fixed_site_mask)
                           if m == 'I')]
    else:
        eligible = all_paulis

    if k_pauli > len(eligible):
        return None  # not enough Paulis

    # Pick k_pauli random Paulis (without replacement)
    indices = rng.choice(len(eligible), size=k_pauli, replace=False)
    selected = [eligible[i] for i in indices]

    # Random complex coefficients
    coeffs = rng.normal(size=k_pauli) + 1j * rng.normal(size=k_pauli)

    c = np.zeros((d, d), dtype=complex)
    for coef, p in zip(coeffs, selected):
        c = c + coef * pauli_string(p)

    # Use a fixed simple Hermitian H (Heisenberg-like) so the only varying part is c
    # For N=2: XX + YY + ZZ on the single bond
    # For N=3: chain XX+YY+ZZ on each of 2 bonds
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for letter in ['X', 'Y', 'Z']:
            P_l = np.eye(1, dtype=complex)
            for site in range(N):
                if site == b or site == b + 1:
                    P_l = np.kron(P_l, PAULI_2X2[letter])
                else:
                    P_l = np.kron(P_l, PAULI_2X2['I'])
            H = H + P_l

    sigma = 0.1
    L_vec = build_L_with_c(H, c, gamma=0.1)
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = result['norm_sq']['M']
    asym = result['asymmetry']
    rel_asym = abs(asym) / max(ns_M, 1e-15)
    return ns_M, asym, rel_asym


def main():
    rng_master = np.random.default_rng(seed=2026)

    for N in [2, 3]:
        d = 2 ** N
        total_paulis = 4 ** N
        all_paulis = enumerate_pauli_strings(N)
        print(f"=" * 90)
        print(f"N = {N}, d = {d}, total Pauli strings = {total_paulis}")
        print(f"=" * 90)

        # Sweep 1: full Pauli basis support, vary k_pauli
        ks = [1, 2, 4, 8, 16, 32, 64, 128, 256]
        ks = [k for k in ks if k <= total_paulis]
        if total_paulis not in ks:
            ks.append(total_paulis)

        print(f"\n--- Sweep 1: full-spread Pauli support, k_pauli sweep ---\n")
        print(f"{'k_pauli':>8}  {'trial':>5}  {'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 70)
        for k_pauli in ks:
            for trial in range(5):
                rng = np.random.default_rng(seed=int(rng_master.integers(0, 10**9)))
                ns_M, asym, rel_asym = test_c_at_pauli_rank(N, k_pauli, all_paulis, rng)
                marker = "BALANCE" if rel_asym < 1e-10 else "BROKEN "
                print(f"{k_pauli:>8}  {trial:>5}  {ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
            print()

        # Sweep 2: restricted to single-site c (mask = X on site 0, I elsewhere) -- max 4 Paulis
        # Wait, we'd need to pick Paulis where only site 0 is non-identity.
        # That gives 4 Paulis: III..., XII..., YII..., ZII... (4 total for any N)
        print(f"--- Sweep 2: single-site c (Pauli at site 0 only), k_pauli sweep ---\n")
        single_site_mask = tuple('I' if i > 0 else '*' for i in range(N))
        # Eligible: Paulis where sites 1..N-1 are all 'I' -> 4 such Paulis
        eligible_single = [p for p in all_paulis if all(p[i] == 'I' for i in range(1, N))]
        print(f"Eligible single-site Paulis (sites 1..{N-1} all I): {len(eligible_single)}")
        for p in eligible_single:
            print(f"  {p}")
        print()

        ks_single = [1, 2, 3, 4]
        print(f"{'k_pauli':>8}  {'trial':>5}  {'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 70)
        for k_pauli in ks_single:
            for trial in range(5):
                rng = np.random.default_rng(seed=int(rng_master.integers(0, 10**9)))
                result = test_c_at_pauli_rank(N, k_pauli, eligible_single, rng)
                if result is None:
                    print(f"{k_pauli:>8}  {trial:>5}  insufficient eligible Paulis")
                    continue
                ns_M, asym, rel_asym = result
                marker = "BALANCE" if rel_asym < 1e-10 else "BROKEN "
                print(f"{k_pauli:>8}  {trial:>5}  {ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
            print()
        print()


if __name__ == '__main__':
    main()
