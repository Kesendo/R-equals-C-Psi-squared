"""Finite tolerance-based centre-eigenvalue census for the N=4 XY chain.

The bare uniform-Z dephasing basis layer n_XY=N/2 has the exact combinatorial
dimension C(4,2) 2^4 = 96. After the XY Hamiltonian is included, this script
finds 94 eigenvalues on Re(lambda)=-sum(gamma_l) at tolerance 1e-9. Thus:

* 94 observed at tolerance 1e-9; basis-layer reference is 96; arithmetic gap 2
  is not a leakage count. A continuation of individual modes would be needed
  before assigning such a mechanism.
* The 18 eigenvalues with additionally |Im(lambda)|<1e-9 and the R-parity
  counts are finite numerical census results.
* Site-reflection R parity does not establish Pi-fixed vectors. The composite
  conjugating fold fixes the real centre line, while the linear spectral fold
  fixes only lambda=-sum(gamma_l).
* Raw eigenvector weights are basis-dependent in a degenerate subspace.
  Majorana language is analogy only; no individual vector theorem follows.

The basis-layer reference dimension is 96. This bounded diagnostic evaluates
N=4 only when run directly; importing it performs no eigendecomposition.

Run: python simulations/axis_modes_n4.py
"""

import itertools
import sys
from math import comb

import numpy as np


J = 1.0
GAMMA = 0.05

PAULI = np.array([
    [[1, 0], [0, 1]],
    [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]],
    [[1, 0], [0, -1]],
], dtype=complex)


def pauli_string(alpha):
    result = PAULI[alpha[0]]
    for a in alpha[1:]:
        result = np.kron(result, PAULI[a])
    return result


def n_xy(alpha):
    return sum(1 for a in alpha if a in (1, 2))


def build_liouvillian(N, gamma, J):
    """Return the XY-chain Liouvillian with uniform local Z dephasing."""
    dim = 2 ** N
    hamiltonian = np.zeros((dim, dim), dtype=complex)
    for site in range(N - 1):
        for axis_index in (1, 2):
            alpha = [0] * N
            alpha[site] = axis_index
            alpha[site + 1] = axis_index
            hamiltonian += (J / 2) * pauli_string(alpha)

    identity = np.eye(dim, dtype=complex)
    hamiltonian_part = -1j * (
        np.kron(identity, hamiltonian)
        - np.kron(hamiltonian.T, identity)
    )
    dissipator = np.zeros((dim * dim, dim * dim), dtype=complex)
    super_identity = np.eye(dim * dim, dtype=complex)
    for site in range(N):
        alpha = [0] * N
        alpha[site] = 3
        z_site = pauli_string(alpha)
        dissipator += gamma * (np.kron(z_site, z_site) - super_identity)
    return hamiltonian_part + dissipator


def axis_census(eigvals, N, gamma, tolerance):
    """Return observed centre count, bare-layer reference, and arithmetic gap."""
    if N != 4:
        raise ValueError("axis_census is restricted to the N=4 diagnostic")
    centre = -N * gamma
    observed = int(np.sum(np.abs(eigvals.real - centre) < tolerance))
    reference = comb(N, N // 2) * (2 ** N)
    return observed, reference, reference - observed


def reverse_bits(value, count):
    result = 0
    for bit in range(count):
        if (value >> bit) & 1:
            result |= 1 << (count - 1 - bit)
    return result


def reflection_parity_counts(vectors, reflection, span_tolerance=1e-10,
                             parity_tolerance=1e-6):
    """Count +/- site-reflection eigenvalues on the spanned numerical subspace."""
    left_vectors, singular_values, _ = np.linalg.svd(
        vectors, full_matrices=False,
    )
    rank = int(np.sum(singular_values > span_tolerance))
    if rank == 0:
        return 0, 0, 0
    orthonormal = left_vectors[:, :rank]
    restricted = orthonormal.conj().T @ reflection @ orthonormal
    restricted = (restricted + restricted.conj().T) / 2
    values = np.linalg.eigvalsh(restricted)
    even = int(np.sum(np.abs(values - 1) < parity_tolerance))
    odd = int(np.sum(np.abs(values + 1) < parity_tolerance))
    return even, odd, len(values) - even - odd


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    N = 4
    dim = 2 ** N
    tolerance = 1e-9
    liouvillian = build_liouvillian(N, GAMMA, J)
    eigvals, eigvecs = np.linalg.eig(liouvillian)
    centre = -N * GAMMA
    mask = np.abs(eigvals.real - centre) < tolerance
    axis_indices = np.where(mask)[0]
    observed_axis, layer_reference, reference_minus_observed = axis_census(eigvals, N, GAMMA, 1e-9)
    assert (observed_axis, layer_reference, reference_minus_observed) == (94, 96, 2)

    print("=== N=4 FINITE TOLERANCE-SELECTED CENTRE EIGENSPACE ===")
    print(f"XY chain, J={J}, uniform gamma={GAMMA}, tolerance={tolerance:g}")
    print(f"basis-layer reference dimension: {layer_reference}")
    print(f"{observed_axis} observed at tolerance 1e-9; basis-layer reference is {layer_reference}; arithmetic gap {reference_minus_observed} is not a leakage count")
    print("Site-reflection R counts below are subspace counts, not Pi-fixed vectors.")
    print("The composite conjugating fold fixes the real centre line; the "
          "linear spectral fold fixes only lambda=-sum(gamma_l).")
    print("Raw eigenvector weights are basis-dependent in a degenerate subspace.")

    rounded_imaginary = np.round(eigvals.imag[mask], 6)
    unique_imaginary, imaginary_counts = np.unique(
        rounded_imaginary, return_counts=True,
    )
    print(f"centre-line Im clusters: {len(unique_imaginary)}")
    for value, count in zip(unique_imaginary, imaginary_counts):
        print(f"  Im(lambda)={value:+.6f}: {count}")

    im_zero_indices = np.where(mask & (np.abs(eigvals.imag) < tolerance))[0]
    print(f"finite |Im(lambda)|<1e-9 count: {len(im_zero_indices)}")
    assert len(im_zero_indices) == 18

    all_alphas = list(itertools.product(range(4), repeat=N))
    pauli_basis = np.zeros((dim * dim, dim * dim), dtype=complex)
    for index, alpha in enumerate(all_alphas):
        pauli_basis[:, index] = pauli_string(alpha).flatten("F")
    nxy_per_alpha = np.array([n_xy(alpha) for alpha in all_alphas])

    mean_support = np.zeros(N + 1)
    for index in axis_indices:
        coefficients = pauli_basis.conj().T @ eigvecs[:, index] / (2 ** N)
        weights = np.abs(coefficients) ** 2
        weights /= weights.sum()
        for layer in range(N + 1):
            mean_support[layer] += weights[nxy_per_alpha == layer].sum()
    mean_support /= len(axis_indices)
    print("mean light-content support of this numerical eigenbasis:")
    for layer, weight in enumerate(mean_support):
        print(f"  n_XY={layer}: {weight:.6f}")

    reflection_hilbert = np.zeros((dim, dim), dtype=complex)
    for basis_index in range(dim):
        reflection_hilbert[reverse_bits(basis_index, N), basis_index] = 1.0
    reflection_operator = np.kron(reflection_hilbert, reflection_hilbert)

    centre_even, centre_odd, centre_ambiguous = reflection_parity_counts(
        eigvecs[:, axis_indices], reflection_operator,
    )
    print(
        "site-reflection R parity on the tolerance-selected centre eigenspace: "
        f"even={centre_even}, odd={centre_odd}, ambiguous={centre_ambiguous}"
    )
    assert (centre_even, centre_odd, centre_ambiguous) == (58, 36, 0)

    zero_even, zero_odd, zero_ambiguous = reflection_parity_counts(
        eigvecs[:, im_zero_indices], reflection_operator,
    )
    print(
        "site-reflection R parity on the finite real-centre subset: "
        f"even={zero_even}, odd={zero_odd}, ambiguous={zero_ambiguous}"
    )
    print("These finite subspace counts do not select individual eigenvectors.")


if __name__ == "__main__":
    main()
