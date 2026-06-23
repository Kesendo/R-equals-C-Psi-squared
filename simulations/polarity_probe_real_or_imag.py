"""Polarity probe #8: refined hypothesis after probe 7.

Probe 7 broke the StandardLindbladPiBalance conjecture: random complex c, H
give nonzero asymmetry. But probes 1-5 (Pauli-structured) preserve balance.

The structural difference: Pauli matrices satisfy P^* = +/- P (X, Z, I are
real; Y is pure imaginary). Their tensor products have entries that are
either real OR pure imaginary, never fully complex. Random complex matrices
have fully complex entries.

Refined conjecture: if L is entrywise "real OR pure imaginary" (i.e.,
L_ij is real or pure imaginary for each entry), then asymmetry = 0.

This probe tests:
  A) c real, H real Hermitian: expect asymmetry = 0
  B) c pure imaginary, H real Hermitian: expect asymmetry = 0
  C) c complex (=A+iB random), H real Hermitian: expect asymmetry != 0
  D) c Pauli string (real or pure imaginary), H Pauli sum (real coefs): expect asymmetry = 0

Each tested at N=2, N=3 with multiple random seeds.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np

from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform, site_op


def build_L_standard(H, jump_ops, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(jump_ops, gammas):
        L = L + g * np.kron(c, c.conj())
    return L


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def test_config(label, H, jump_ops, gammas, sigma, N):
    L_vec = build_L_standard(H, jump_ops, gammas)
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = result['norm_sq']['M']
    asym = result['asymmetry']
    rel_asym = abs(asym) / max(ns_M, 1e-15)
    marker = "BALANCE" if rel_asym < 1e-10 else "BROKEN "
    print(f"  {label:<60}  ||M||^2={ns_M:>12.4f}  rel_asym={rel_asym:>10.3e}  [{marker}]")


def main():
    rng = np.random.default_rng(seed=2026)

    for N in [2, 3]:
        d = 2 ** N
        sigma = 0.05 * N
        print(f"=== N = {N}, d = {d}, sigma = {sigma} ===\n")

        for trial in range(3):
            seed = 2026 + trial * 1000 + N * 10000
            sub_rng = np.random.default_rng(seed)
            print(f"Trial {trial} (seed={seed}):")

            # A) c real, H real Hermitian
            H_real = sub_rng.normal(size=(d, d))
            H_real = (H_real + H_real.T) / 2  # real Hermitian
            c_real = sub_rng.normal(size=(d, d))  # real, any structure
            test_config("A: c real, H real Hermitian",
                        H_real.astype(complex), [c_real.astype(complex)],
                        [0.1], sigma, N)

            # B) c pure imaginary, H real Hermitian
            c_imag = 1j * sub_rng.normal(size=(d, d))  # pure imaginary entries
            test_config("B: c pure imaginary, H real Hermitian",
                        H_real.astype(complex), [c_imag],
                        [0.1], sigma, N)

            # C) c complex (real + i*real), H real Hermitian
            c_complex = sub_rng.normal(size=(d, d)) + 1j * sub_rng.normal(size=(d, d))
            test_config("C: c fully complex, H real Hermitian",
                        H_real.astype(complex), [c_complex],
                        [0.1], sigma, N)

            # D) c Pauli string (single site or product), H Pauli sum
            # Use site_op which returns Pauli string tensor products
            if N >= 2:
                c_pauli_X0 = site_op(N, 0, 'X')  # X on site 0
                c_pauli_Y1 = site_op(N, 1, 'Y')  # Y on site 1
                H_pauli = site_op(N, 0, 'X') @ site_op(N, 1, 'X')  # XX (real-entries)
                if N >= 3:
                    H_pauli = H_pauli + site_op(N, 0, 'Y') @ site_op(N, 1, 'Y')  # + YY (real)
                test_config("D: c Pauli strings, H Pauli sum",
                            H_pauli.astype(complex),
                            [c_pauli_X0.astype(complex), c_pauli_Y1.astype(complex)],
                            [0.1, 0.05], sigma, N)

            # E) c is a Pauli string with imaginary coefficient (still entrywise real-or-imag)
            c_pauli_iY = 1j * site_op(N, 0, 'Z')  # i * Z, pure imaginary entries
            test_config("E: c = i*Z (pure imaginary entries), H real-entry Pauli",
                        H_pauli.astype(complex), [c_pauli_iY],
                        [0.1], sigma, N)

            # F) c with mixed real and pure-imag Pauli sum (entries can be complex)
            c_mixed = site_op(N, 0, 'X') + 1j * site_op(N, 0, 'Z')  # X + i*Z
            test_config("F: c = X + i*Z (mixed complex entries), H real-entry Pauli",
                        H_pauli.astype(complex), [c_mixed.astype(complex)],
                        [0.1], sigma, N)

            print()


if __name__ == '__main__':
    main()
