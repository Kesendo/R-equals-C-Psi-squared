"""Finite tolerance-selected centre-eigenvalue census, restricted to N=4 only.

This historical parameterized script once treated the bare-dephasing
n_XY=N/2 basis-layer dimension as a prediction for the interacting spectrum.
They are different objects. The basis-layer dimension reference is exactly 96
at N=4, while the interacting XY-chain calculation selects 94 eigenvalues at
tolerance 1e-9. The arithmetic difference is not a leakage count.

Site-reflection R parity is a subspace classification and does not imply
Pi-fixed vectors. Raw eigenvector weights are basis-dependent in a degenerate
subspace. Majorana terminology is an analogy, not an individual-vector result.

This compatibility entry point is import-silent and deliberately N=4 only.
Task 6 validates its source and import seam but does not execute it; the bounded
reference run lives in ``axis_modes_n4.py``.
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
    for value in alpha[1:]:
        result = np.kron(result, PAULI[value])
    return result


def build_liouvillian(N, gamma, coupling):
    """Construct the finite XY-chain Liouvillian used by this diagnostic."""
    dimension = 2 ** N
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for site in range(N - 1):
        for axis_index in (1, 2):
            alpha = [0] * N
            alpha[site] = axis_index
            alpha[site + 1] = axis_index
            hamiltonian += (coupling / 2) * pauli_string(alpha)

    identity = np.eye(dimension, dtype=complex)
    generator = -1j * (
        np.kron(identity, hamiltonian)
        - np.kron(hamiltonian.T, identity)
    )
    super_identity = np.eye(dimension * dimension, dtype=complex)
    for site in range(N):
        alpha = [0] * N
        alpha[site] = 3
        z_site = pauli_string(alpha)
        generator += gamma * (np.kron(z_site, z_site) - super_identity)
    return generator


def centre_census(eigenvalues, N, gamma, tolerance):
    centre = -N * gamma
    observed = int(np.sum(np.abs(eigenvalues.real - centre) < tolerance))
    reference = comb(N, N // 2) * (2 ** N)
    return observed, reference


def main(argv=None):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    arguments = sys.argv[1:] if argv is None else list(argv)
    N = int(arguments[0]) if arguments else 4
    if N != 4:
        raise ValueError("This repaired compatibility diagnostic is N=4 only")

    eigenvalues = np.linalg.eigvals(build_liouvillian(N, GAMMA, J))
    observed, reference = centre_census(eigenvalues, N, GAMMA, 1e-9)
    print("=== N=4 TOLERANCE-SELECTED CENTRE EIGENSPACE ===")
    print(f"observed at tolerance 1e-9: {observed}")
    print(f"basis-layer dimension reference: {reference}")
    print("Site-reflection R parity is a subspace label, not Pi-fixed vectors.")
    print("Raw eigenvector weights in degenerate spaces are basis-dependent.")
    print("No individual Majorana or vector identity is inferred.")


if __name__ == "__main__":
    main()
