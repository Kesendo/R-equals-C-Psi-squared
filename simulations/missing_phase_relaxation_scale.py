"""Invariant reduction for the N=7 centre-watched missing-phase A block."""

import numpy as np
import sympy


N = 7
SEAT = 3


def _exact_scalar(value):
    value = sympy.sympify(value)
    if value.has(sympy.Float):
        raise TypeError("exact reduction requires exact scalar inputs")
    return value


def hopping_float(epsilon: float, both_ends: bool = False) -> np.ndarray:
    """Return the N=7 path Hamiltonian with the requested end detuning."""
    r = 1.0 + epsilon
    bonds = [2.0 * r, 2.0, 2.0, 2.0, 2.0, 2.0 * r if both_ends else 2.0]
    h = np.zeros((N, N), dtype=complex)
    for site, bond in enumerate(bonds):
        h[site, site + 1] = h[site + 1, site] = bond
    return h


def blind_vector_float(epsilon: float) -> np.ndarray:
    """Return the normalized reflection-broken blind vector."""
    r = 1.0 + epsilon
    v = np.array([1.0, 0.0, -r, 0.0, r, 0.0, -r], dtype=complex)
    return v / np.linalg.norm(v)


def watched_involution_float(seat: int = SEAT) -> np.ndarray:
    """Return the diagonal Z involution that watches one seat."""
    z = np.eye(N, dtype=complex)
    z[seat, seat] = -1.0
    return z


def k_reduced_float(epsilon: float, gamma: float, *, seat: int = SEAT) -> np.ndarray:
    """Return K = -i h - 2 gamma |seat><seat| for the invariant right factor."""
    projector = np.zeros((N, N), dtype=complex)
    projector[seat, seat] = 1.0
    return -1j * hopping_float(epsilon) - 2.0 * gamma * projector


def a_generator_float(epsilon: float, gamma: float, seat: int = SEAT) -> np.ndarray:
    """Assemble the 49-dimensional A generator for row-major stacking."""
    h = hopping_float(epsilon)
    identity = np.eye(N, dtype=complex)
    z = watched_involution_float(seat)
    return (
        -1j * (np.kron(h, identity) - np.kron(identity, h.T))
        + gamma * (np.kron(z, z.T) - np.eye(N * N, dtype=complex))
    )


def row_stack_right_embedding(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Row-stack the right invariant embedding u v-dagger."""
    return np.kron(u, np.conj(v))


def row_stack_adjoint_embedding(v: np.ndarray, u_plus: np.ndarray) -> np.ndarray:
    """Row-stack the adjoint invariant embedding v u_plus-dagger."""
    return np.kron(v, np.conj(u_plus))


def exact_reduction_objects(epsilon, gamma):
    """Return h, z, the blind vector, and K over exact SymPy scalars."""
    epsilon = _exact_scalar(epsilon)
    gamma = _exact_scalar(gamma)
    r = 1 + epsilon
    h = sympy.zeros(N, N)
    for site, bond in enumerate((2 * r, 2, 2, 2, 2, 2)):
        h[site, site + 1] = h[site + 1, site] = bond

    z = sympy.eye(N)
    z[SEAT, SEAT] = -1
    v = sympy.Matrix([1, 0, -r, 0, r, 0, -r])
    projector = sympy.zeros(N, N)
    projector[SEAT, SEAT] = 1
    k = -sympy.I * h - 2 * gamma * projector
    return h, z, v, k


def exact_a_action(h, z, gamma, x):
    """Apply the centre-watched A generator in exact arithmetic."""
    action = -sympy.I * (h * x - x * h) + gamma * (z * x * z - x)
    return action.applyfunc(sympy.expand)


def reduction_mutation_residuals(epsilon, gamma):
    """Return exact squared residuals for four reductions-breaking mutations."""
    h, z, v, k = exact_reduction_objects(epsilon, gamma)
    u = sympy.Matrix([1, sympy.I, 2, 0, -1, 3 * sympy.I, 4])
    x = u * sympy.conjugate(v.T)
    target = k * u * sympy.conjugate(v.T)

    watched_seat = sympy.eye(N)
    watched_seat[2, 2] = -1
    wrong_side = v * sympy.conjugate(u.T)
    mutants = {
        "commutator_sign": (
            sympy.I * (h * x - x * h) + gamma * (z * x * z - x)
        ),
        "operator_side": exact_a_action(h, z, gamma, wrong_side),
        "watched_seat": exact_a_action(h, watched_seat, gamma, x),
        "missing_minus_identity": -sympy.I * (h * x - x * h) + gamma * z * x * z,
    }

    def squared_frobenius(delta):
        return sympy.simplify(sum(sympy.conjugate(entry) * entry for entry in delta))

    return {
        name: squared_frobenius(mutant - target)
        for name, mutant in mutants.items()
    }


def characteristic_polynomial(lam, r, gamma):
    """Return det(lambda I - K) from the exact reduced 7x7 construction."""
    lam = _exact_scalar(lam)
    r = _exact_scalar(r)
    gamma = _exact_scalar(gamma)
    _, _, _, k = exact_reduction_objects(r - 1, gamma)
    return sympy.expand((lam * sympy.eye(N) - k).det())


def minus_branch_series(epsilon, gamma):
    """Return the cubic germ of the branch from lambda = -2 sqrt(2) i."""
    epsilon = _exact_scalar(epsilon)
    gamma = _exact_scalar(gamma)
    return (
        -2 * sympy.sqrt(2) * sympy.I
        - sympy.I * epsilon / sympy.sqrt(2)
        - (gamma / 2 + 3 * sympy.sqrt(2) * sympy.I / 16) * epsilon**2
        - (gamma / 2 - 19 * sympy.sqrt(2) * sympy.I / 64) * epsilon**3
    )


def gap_series(epsilon, gamma):
    """Return the cubic local relaxation-gap germ."""
    epsilon = _exact_scalar(epsilon)
    gamma = _exact_scalar(gamma)
    return gamma * epsilon**2 / 2 + gamma * epsilon**3 / 2
