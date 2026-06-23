"""Polarity probe #12: verify the Z2-cubed sub-cell pattern at N >= 3.

Probe 11 (N=2) found two classes of Pauli-string pairs:
  Class A (same Z_2^3 sub-cell):     BALANCE for ALL coefficient choices
  Class B (cross Z_2^3 sub-cell):    BALANCE only for phase-matched coefficients

Z_2^3 cell = (Klein_a, Klein_b, y_par) = ((#X+#Y) mod 2, (#Y+#Z) mod 2, #Y mod 2)

This probe extends the test to N=3 and N=4 to verify the pattern scales.

Sweeps:
  Sweep A: within-cell c (k_pauli = 2, 3, 4 Pauli strings, all in same cell)
           with random complex coefficients. Conjecture: BALANCE always.
  Sweep B: cross-cell c (k_pauli >= 2 Pauli strings spanning >= 2 cells)
           with random complex coefficients. Conjecture: BROKEN mostly.
  Sweep C: same as A but with random REAL coefficients (sanity check that
           ANY coefs work for same-cell).

If pattern scales: F112 candidate (StandardLindbladZ2CubedPiBalance) is real.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import product
from collections import defaultdict

from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform

PAULI_2X2 = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string(letters):
    op = PAULI_2X2[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI_2X2[l])
    return op


def z2cubed_cell(p):
    """(Klein_a, Klein_b, y_par) for Pauli string p."""
    nx = p.count('X')
    ny = p.count('Y')
    nz = p.count('Z')
    return ((nx + ny) % 2, (ny + nz) % 2, ny % 2)


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def heisenberg_H(N):
    """H = sum_b (XX + YY + ZZ) on bonds."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for letter in ['X', 'Y', 'Z']:
            letters = ['I'] * N
            letters[b] = letter
            letters[b + 1] = letter
            H = H + pauli_string(letters)
    return H


def test_c(N, c, H, sigma):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T)) + 0.1 * np.kron(c, c.conj())
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = result['norm_sq']['M']
    rel_asym = abs(result['asymmetry']) / max(ns_M, 1e-15)
    return ns_M, rel_asym


def main():
    for N in [2, 3, 4]:
        d = 2 ** N
        sigma = 0.1
        H = heisenberg_H(N)

        # Group all 4^N Pauli strings by Z_2^3 cell
        all_paulis = list(product(['I', 'X', 'Y', 'Z'], repeat=N))
        cells = defaultdict(list)
        for p in all_paulis:
            cells[z2cubed_cell(p)].append(p)

        print("=" * 100)
        print(f"N = {N}, d = {d}, total Pauli strings = {len(all_paulis)}, Z_2^3 cells = {len(cells)}")
        print("=" * 100)
        for cell, members in sorted(cells.items()):
            sample = members[:3] + (['...'] if len(members) > 3 else [])
            print(f"  cell {cell}: {len(members):>3} members, e.g. "
                  f"{['' .join(p) if isinstance(p, tuple) else p for p in sample]}")
        print()

        rng = np.random.default_rng(seed=2026 + N * 1000)

        # ---- Sweep A: within-cell c, k_pauli in {2, 3, 4} ----
        print(f"--- Sweep A (within-cell c, random complex coefs, 3 trials) ---")
        print(f"{'cell':<12}  {'k':<3}  {'trial':<5}  {'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 75)
        for cell, members in sorted(cells.items()):
            if len(members) < 2:
                continue
            for k_pauli in [2, 3, 4]:
                if k_pauli > len(members):
                    continue
                for trial in range(3):
                    indices = rng.choice(len(members), size=k_pauli, replace=False)
                    selected = [members[i] for i in indices]
                    coeffs = rng.normal(size=k_pauli) + 1j * rng.normal(size=k_pauli)
                    c = sum(coef * pauli_string(p) for coef, p in zip(coeffs, selected))
                    ns_M, rel_asym = test_c(N, c, H, sigma)
                    marker = "BAL" if rel_asym < 1e-10 else "BREAK"
                    print(f"{str(cell):<12}  {k_pauli:<3}  {trial:<5}  {ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
        print()

        # ---- Sweep B: cross-cell c (pick 1 from each of 2 cells), random coefs ----
        print(f"--- Sweep B (cross-cell c, 1 Pauli from each of 2 cells, random complex coefs, 5 trials) ---")
        print(f"{'cell_a':<10}  {'cell_b':<10}  {'p_a':<8}  {'p_b':<8}  {'trial':<5}  "
              f"{'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 90)
        cell_list = list(cells.keys())
        for i, cell_a in enumerate(cell_list):
            for cell_b in cell_list[i+1:]:
                if not cells[cell_a] or not cells[cell_b]:
                    continue
                # Pick one representative each
                p_a = cells[cell_a][0]
                p_b = cells[cell_b][0]
                for trial in range(3):
                    a, b = (complex(rng.normal(), rng.normal()),
                            complex(rng.normal(), rng.normal()))
                    c = a * pauli_string(p_a) + b * pauli_string(p_b)
                    ns_M, rel_asym = test_c(N, c, H, sigma)
                    marker = "BAL" if rel_asym < 1e-10 else "BREAK"
                    print(f"{str(cell_a):<10}  {str(cell_b):<10}  "
                          f"{''.join(p_a):<8}  {''.join(p_b):<8}  {trial:<5}  "
                          f"{ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
        print()

        # ---- Sweep C: within-cell with random REAL coefs (extra sanity) ----
        print(f"--- Sweep C (within-cell c, random REAL coefs, 3 trials per cell, k=4) ---")
        print(f"{'cell':<12}  {'trial':<5}  {'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 60)
        for cell, members in sorted(cells.items()):
            if len(members) < 4:
                continue
            for trial in range(3):
                indices = rng.choice(len(members), size=4, replace=False)
                selected = [members[i] for i in indices]
                coeffs = rng.normal(size=4)  # REAL coefficients
                c = sum(coef * pauli_string(p) for coef, p in zip(coeffs, selected))
                ns_M, rel_asym = test_c(N, c, H, sigma)
                marker = "BAL" if rel_asym < 1e-10 else "BREAK"
                print(f"{str(cell):<12}  {trial:<5}  {ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
        print()
        print()


if __name__ == '__main__':
    main()
