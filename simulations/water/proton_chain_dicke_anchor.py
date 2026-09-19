"""Abstract qubit-chain calculation for a Dicke-state parity readout.

The model is an XXX Hamiltonian with uniform all-site Z dephasing on abstract
two-level systems. It is not a proton or hydrogen-bond model. Water supplies a
selected qubit-model analogy only; this source contains no molecular geometry,
electronic surface, bath calibration, or chemistry confirmation.

Two results must remain separate:

1. The initial Dicke-superposition letter-parity fraction follows from the
   X^N overlap and is checked for named finite N=3..6 states.
2. The long-time fraction (N+2)/(4(N+1)) follows, for the stated stationary
   mixture, from conserved excitation weight and sector-uniform stationary
   diagonal together with an exact Krawtchouk sum. The finite N=4..16
   combinatorial checks do not depend on the numerical dynamics.

The limiting scalar 1/4 is a limit of this selected construction. It does not
establish a Mandelbrot or C-Psi boundary, a material mechanism, or a universal
open-system transition. This heavy historical diagnostic is syntax/source-copy
only in Task 6 and must not be imported or executed by its verification gate.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import comb

import numpy as np


I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = (I2, SX, SY, SZ)
LETTER_PARITY = (0, 0, 1, 1)


def kron_n(matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = np.kron(result, matrix)
    return result


def site_op(operator, site, N):
    return kron_n([operator if index == site else I2 for index in range(N)])


def heisenberg_chain(N, J=1.0):
    """Return H=(J/4) sum(XX+YY+ZZ), the selected open XXX chain."""
    hamiltonian = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for bond in range(N - 1):
        for operator in (SX, SY, SZ):
            hamiltonian += (
                (J / 4.0)
                * site_op(operator, bond, N)
                @ site_op(operator, bond + 1, N)
            )
    return hamiltonian


def z_dephasing_lindblad(hamiltonian, gamma, N):
    """Return the column-vectorized generator with local Z dephasing."""
    dimension = 2 ** N
    identity = np.eye(dimension, dtype=complex)
    generator = -1j * (
        np.kron(identity, hamiltonian)
        - np.kron(hamiltonian.T, identity)
    )
    super_identity = np.kron(identity, identity)
    for site in range(N):
        z_site = site_op(SZ, site, N)
        generator += gamma * (
            np.kron(z_site.conj(), z_site) - super_identity
        )
    return generator


def dicke_state(N, weight):
    """Return the normalized computational-weight Dicke state."""
    state = np.zeros(2 ** N, dtype=complex)
    amplitude = 1.0 / np.sqrt(comb(N, weight))
    for basis_index in range(2 ** N):
        if basis_index.bit_count() == weight:
            state[basis_index] = amplitude
    return state


def dicke_superposition(N, weight):
    """Return (|D_weight>+|D_(weight+1)>)/sqrt(2)."""
    if not 0 <= weight < N:
        raise ValueError("requires 0 <= weight < N")
    return (
        dicke_state(N, weight) + dicke_state(N, weight + 1)
    ) / np.sqrt(2.0)


def density_matrix(state):
    return np.outer(state, state.conj())


def letter_parity_split(operator, N):
    """Split an operator by the chosen Y/Z Pauli-letter parity."""
    dimension = 2 ** N
    even = np.zeros_like(operator)
    odd = np.zeros_like(operator)
    for encoded in range(4 ** N):
        value = encoded
        indices = []
        for _ in range(N):
            indices.append(value & 3)
            value >>= 2
        basis_operator = kron_n([PAULIS[index] for index in indices])
        coefficient = np.trace(basis_operator @ operator) / dimension
        if abs(coefficient) < 1e-14:
            continue
        if sum(LETTER_PARITY[index] for index in indices) % 2:
            odd += coefficient * basis_operator
        else:
            even += coefficient * basis_operator
    total_norm = np.linalg.norm(operator, "fro") ** 2
    odd_norm = np.linalg.norm(odd, "fro") ** 2
    return even, odd, total_norm, odd_norm


def x_global_overlap(state, N):
    x_global = kron_n([SX] * N)
    return float(np.real(state.conj() @ x_global @ state))


def classify_dicke_anchor(N, weight):
    if 2 * weight + 1 == N:
        return "Mirror"
    if N % 2 == 0 and weight in (N // 2 - 1, N // 2):
        return "KIntermediate"
    return "Generic"


def alpha_total_predicted(anchor):
    return {"Mirror": 0.0, "KIntermediate": 3.0 / 8.0, "Generic": 0.5}[anchor]


def gamma_predicted(anchor):
    return {"Mirror": 1.0, "KIntermediate": 0.5, "Generic": 0.0}[anchor]


def alpha_infinity_closed_form(N):
    """Exact selected-mixture limit for even N."""
    if N < 2 or N % 2:
        raise ValueError("requires positive even N")
    return Fraction(N + 2, 4 * (N + 1))


def krawtchouk_subset_sum(N, weight, subset_size):
    """Sum (-1)^intersection over weight-bit strings for a fixed subset."""
    lower = max(0, weight - (N - subset_size))
    upper = min(weight, subset_size)
    return sum(
        (-1) ** intersection
        * comb(subset_size, intersection)
        * comb(N - subset_size, weight - intersection)
        for intersection in range(lower, upper + 1)
    )


def odd_krawtchouk_norm(N, weight):
    """Exact odd-letter projector norm from the Krawtchouk census."""
    numerator = sum(
        comb(N, subset_size)
        * krawtchouk_subset_sum(N, weight, subset_size) ** 2
        for subset_size in range(1, N + 1, 2)
    )
    return Fraction(numerator, 2 ** N)


def evolve_from_diagonalization(rho0, generator, time):
    values, vectors = np.linalg.eig(generator)
    inverse = np.linalg.inv(vectors)
    initial = rho0.flatten(order="F")
    evolved = vectors @ (np.exp(values * time) * (inverse @ initial))
    matrix = evolved.reshape(rho0.shape, order="F")
    return (matrix + matrix.conj().T) / 2


def analyse_anchor(N, weight, J=1.0, gamma=0.05, times=None):
    """Print one finite numerical trajectory in the stipulated N<=6 model."""
    if times is None:
        times = (0.0, 1.0, 5.0, 10.0, 25.0)
    anchor = classify_dicke_anchor(N, weight)
    state = dicke_superposition(N, weight)
    overlap = x_global_overlap(state, N)
    rho0 = density_matrix(state)
    _, _, total0, odd0 = letter_parity_split(rho0, N)
    alpha0 = odd0 / total0
    assert abs(overlap - gamma_predicted(anchor)) < 1e-12
    assert abs(alpha0 - alpha_total_predicted(anchor)) < 1e-12

    generator = z_dephasing_lindblad(heisenberg_chain(N, J), gamma, N)
    print(f"N={N}, weight={weight}, class={anchor}, alpha(0)={alpha0:.12f}")
    for time in times:
        rho = evolve_from_diagonalization(rho0, generator, time)
        _, _, total, odd = letter_parity_split(rho, N)
        print(f"  t={time:6.1f}: selected letter-parity fraction={odd / total:.12f}")


def run():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    print("=== ABSTRACT QUBIT-CHAIN CALCULATION ===")
    print("XXX Hamiltonian plus all-site Z dephasing.")
    print("This is not a proton or hydrogen-bond model.")
    print("Water is a selected qubit-model analogy; there is no chemistry confirmation.")

    for N, weight in ((3, 1), (4, 0), (4, 1), (4, 2), (5, 2), (6, 2)):
        analyse_anchor(N, weight)

    print("\nExact F98-side Krawtchouk combinatorics, separate from dynamics:")
    for N in range(4, 17, 2):
        weight = N // 2 - 1
        observed_norm = odd_krawtchouk_norm(N, weight)
        expected_norm = Fraction(comb(N, weight), 2)
        assert observed_norm == expected_norm
        f_pred = alpha_infinity_closed_form(N)
        alpha_predicted = float(f_pred)
        print(
            f"  N={N}: odd norm={observed_norm}; "
            f"selected-mixture alpha(infinity)={f_pred} ({alpha_predicted:.12f})"
        )

    print("The limit is 1/4 in this selected qubit construction.")
    print("It is not a Mandelbrot or C-Psi boundary and gives no material claim.")


if __name__ == "__main__":
    run()
