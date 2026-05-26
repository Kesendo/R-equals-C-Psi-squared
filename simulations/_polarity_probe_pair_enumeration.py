"""Polarity probe #10: enumerate all Pauli string PAIRS at N=2, find pattern.

Post-probe-9 finding (2026-05-26):
    k_pauli = 1 always preserves balance.
    k_pauli >= 2 is selection-dependent: some Pauli string pairs preserve
    balance, others break it.

This probe enumerates ALL C(16, 2) = 120 pairs at N=2 (+ 16 self-pairs).
For each pair (P_alpha, P_beta), test c = P_alpha + i*P_beta and classify
preserve/broken. Then look for the algebraic pattern in preserved pairs:
  - Klein index relationship
  - Pi^2 parity relationship
  - Pi-orbit relationship
  - Same site vs different sites
  - Number-of-Y relationship

Hopefully a clean pattern emerges.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import combinations_with_replacement, product

from framework.symmetry import build_pi_full
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


def klein_index(p):
    """(bit_a, bit_b) = (#X + #Y mod 2, #Y + #Z mod 2)"""
    nx = p.count('X')
    ny = p.count('Y')
    nz = p.count('Z')
    return ((nx + ny) % 2, (ny + nz) % 2)


def n_y(p):
    return p.count('Y')


def n_nonI(p):
    return sum(1 for l in p if l != 'I')


def build_L_with_c(H, c, gamma):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    L = L + gamma * np.kron(c, c.conj())
    return L


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def test_pair(N, p_alpha, p_beta, H, sigma):
    """Test c = P_alpha + i * P_beta. Returns rel_asym."""
    d = 2 ** N
    c = pauli_string(p_alpha) + 1j * pauli_string(p_beta)
    L_vec = build_L_with_c(H, c, gamma=0.1)
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = result['norm_sq']['M']
    asym = result['asymmetry']
    rel_asym = abs(asym) / max(ns_M, 1e-15)
    return ns_M, rel_asym


def main():
    N = 2
    d = 2 ** N
    sigma = 0.1

    # Fixed H = Heisenberg XX + YY + ZZ
    H = np.zeros((d, d), dtype=complex)
    for letter in ['X', 'Y', 'Z']:
        H = H + pauli_string([letter, letter])

    all_paulis = list(product(['I', 'X', 'Y', 'Z'], repeat=N))
    print(f"N = {N}, total Pauli strings = {len(all_paulis)}")
    print()
    print(f"Testing c = P_alpha + i * P_beta for all (alpha, beta) including alpha = beta")
    print(f"H = XX + YY + ZZ (Heisenberg), gamma = 0.1, sigma = {sigma}")
    print()

    preserve = []
    broken = []

    for p_alpha in all_paulis:
        for p_beta in all_paulis:
            if p_alpha > p_beta:
                continue  # unordered pairs (with self)
            ns_M, rel_asym = test_pair(N, p_alpha, p_beta, H, sigma)
            tag = (p_alpha, p_beta, ns_M, rel_asym)
            if rel_asym < 1e-10:
                preserve.append(tag)
            else:
                broken.append(tag)

    print(f"Total ordered pairs tested: {len(preserve) + len(broken)}")
    print(f"  Preserved: {len(preserve)}")
    print(f"  Broken:    {len(broken)}")
    print()

    # Look for pattern in preserved pairs
    print("=" * 80)
    print(f"Preserved pairs (sorted by P_alpha then P_beta):")
    print("=" * 80)
    print(f"{'alpha':<8} {'beta':<8} {'klein_a':<8} {'klein_b':<8} {'#Y_a':<5} {'#Y_b':<5}")
    for p_a, p_b, ns_M, rel in preserve:
        k_a = klein_index(p_a)
        k_b = klein_index(p_b)
        print(f"{''.join(p_a):<8} {''.join(p_b):<8} {str(k_a):<8} {str(k_b):<8} {n_y(p_a):<5} {n_y(p_b):<5}")
    print()

    # Pattern analysis: in preserved pairs, what's the relationship between alpha and beta?
    print("=" * 80)
    print("Pattern analysis (preserved pairs):")
    print("=" * 80)

    # 1. Same Klein?
    same_klein = sum(1 for p_a, p_b, _, _ in preserve if klein_index(p_a) == klein_index(p_b))
    print(f"  Same Klein index:  {same_klein} / {len(preserve)}")

    # 2. Both have same #Y parity?
    same_y_par = sum(1 for p_a, p_b, _, _ in preserve if (n_y(p_a) % 2) == (n_y(p_b) % 2))
    print(f"  Same y_par (n_Y mod 2): {same_y_par} / {len(preserve)}")

    # 3. alpha = beta (self-pair)?
    self_pair = sum(1 for p_a, p_b, _, _ in preserve if p_a == p_b)
    print(f"  Self pair (alpha = beta): {self_pair} / {len(preserve)}")

    # 4. One of them is I*N?
    one_is_identity = sum(1 for p_a, p_b, _, _ in preserve
                          if all(l == 'I' for l in p_a) or all(l == 'I' for l in p_b))
    print(f"  One is I*N (identity Pauli): {one_is_identity} / {len(preserve)}")

    # 5. Same n_nonI?
    same_nonI = sum(1 for p_a, p_b, _, _ in preserve if n_nonI(p_a) == n_nonI(p_b))
    print(f"  Same n_nonI: {same_nonI} / {len(preserve)}")

    # 6. Self-pair (alpha = beta) only?
    print(f"  Just self-pairs would be at most: 16; observed self-pairs preserved: {self_pair}")

    # 7. By relative property: alpha = beta OR alpha and beta are related how?
    print()
    print("Non-self preserved pairs:")
    for p_a, p_b, _, _ in preserve:
        if p_a != p_b:
            print(f"  {''.join(p_a)} + i*{''.join(p_b)}  klein={klein_index(p_a)},{klein_index(p_b)}  y_par={n_y(p_a)%2},{n_y(p_b)%2}")


if __name__ == '__main__':
    main()
