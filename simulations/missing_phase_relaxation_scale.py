"""Invariant reduction for the N=7 centre-watched missing-phase A block."""

import numpy as np
import sympy
from sympy.polys.matrices import DomainMatrix


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


def _row_stack_exact(matrix):
    return sympy.Matrix(
        [matrix[row, column] for row in range(N) for column in range(N)]
    )


def _exact_generator_from_h_z(h, z, gamma):
    columns = []
    for coordinate in range(N * N):
        basis = sympy.zeros(N, N)
        basis[coordinate // N, coordinate % N] = 1
        columns.append(_row_stack_exact(exact_a_action(h, z, gamma, basis)))
    return sympy.Matrix.hstack(*columns)


def exact_a_generator(epsilon, gamma):
    """Return the exact 49 by 49 row-stacked A generator."""
    h, z, _, _ = exact_reduction_objects(epsilon, gamma)
    return _exact_generator_from_h_z(h, z, gamma)


def _exact_rank(matrix):
    """Use SymPy's exact domain backend rather than heuristic simplification."""
    return DomainMatrix.from_Matrix(matrix).rank()


def _normalized(vector):
    norm_squared = sympy.simplify((sympy.conjugate(vector.T) * vector)[0])
    return vector.applyfunc(lambda entry: sympy.simplify(entry / sympy.sqrt(norm_squared)))


def _uniform_blind_energy_basis():
    """Derive the ordered (d_minus, d_zero, d_plus) basis from ker K^T."""
    h, _, _, _ = exact_reduction_objects(sympy.Integer(0), sympy.Integer(1))
    watched = sympy.eye(N).col(SEAT)
    krylov = sympy.Matrix.hstack(*[(h**power) * watched for power in range(N)])
    blind_columns = krylov.T.nullspace()
    if len(blind_columns) != 3:
        raise AssertionError("uniform centre-seat blind space must have dimension three")

    blind = sympy.Matrix.hstack(*blind_columns)
    gram = sympy.conjugate(blind.T) * blind
    restriction = gram.inv() * sympy.conjugate(blind.T) * h * blind
    eigenspaces = []
    for energy, multiplicity, coefficient_vectors in restriction.eigenvects():
        if multiplicity != 1 or len(coefficient_vectors) != 1:
            raise AssertionError("blind Hamiltonian must have three simple energies")
        vector = _normalized(blind * coefficient_vectors[0])
        eigenspaces.append((sympy.simplify(energy), vector))
    eigenspaces.sort(key=lambda pair: float(sympy.N(pair[0])))
    return h, eigenspaces, blind, krylov


def _peripheral_bases():
    h, eigenspaces, _, _ = _uniform_blind_energy_basis()
    blind_projector = sympy.zeros(N, N)
    for _, vector in eigenspaces:
        blind_projector += vector * sympy.conjugate(vector.T)
    outer_identity = sympy.eye(N) - blind_projector

    frequencies = {
        "-4sqrt2i": -4 * sympy.sqrt(2) * sympy.I,
        "-2sqrt2i": -2 * sympy.sqrt(2) * sympy.I,
        "0": sympy.Integer(0),
        "+2sqrt2i": 2 * sympy.sqrt(2) * sympy.I,
        "+4sqrt2i": 4 * sympy.sqrt(2) * sympy.I,
    }
    operators = {key: [] for key in frequencies}
    for left_energy, left in eigenspaces:
        for right_energy, right in eigenspaces:
            frequency = sympy.simplify(-sympy.I * (left_energy - right_energy))
            key = next(
                name
                for name, expected in frequencies.items()
                if sympy.simplify(frequency - expected) == 0
            )
            operators[key].append(left * sympy.conjugate(right.T))
    operators["0"].append(outer_identity / 2)
    return h, eigenspaces, outer_identity, frequencies, operators


def exact_uniform_peripheral_census(gamma):
    """Count the complete uniform peripheral algebra in exact arithmetic."""
    gamma = _exact_scalar(gamma)
    h, eigenspaces, outer_identity, frequencies, operators = _peripheral_bases()
    _, z, _, _ = exact_reduction_objects(sympy.Integer(0), gamma)

    counted = {}
    for key, basis in operators.items():
        frequency = frequencies[key]
        for operator in basis:
            if exact_a_action(h, z, gamma, operator) != frequency * operator:
                raise AssertionError("peripheral operator failed its exact frequency gate")
        counted[key] = len(basis)

    generator = exact_a_generator(sympy.Integer(0), gamma)
    kernel_dimension = N * N - _exact_rank(generator)
    if _exact_rank(outer_identity) != 4:
        raise AssertionError("outer Hilbert complement must have dimension four")
    return {
        "blind_dimension": len(eigenspaces),
        "end_blind_algebra_dimension": len(eigenspaces) ** 2,
        "identity_outer_dimension": 1,
        "kernel_dimension": kernel_dimension,
        "peripheral_dimension": sum(counted.values()),
        "frequency_multiplicities": counted,
    }


def _domain_solve(matrix, rhs):
    """Solve a small exact block over an explicit algebraic function field."""
    symbols = sorted(
        set().union(
            *(entry.free_symbols for entry in matrix),
            *(entry.free_symbols for entry in rhs),
        ),
        key=str,
    )
    domain = sympy.QQ.algebraic_field(sympy.I, sympy.sqrt(2))
    if symbols:
        domain = domain.frac_field(*symbols)
    domain_matrix = DomainMatrix.from_Matrix(matrix).convert_to(domain)
    domain_rhs = DomainMatrix.from_Matrix(rhs).convert_to(domain)
    numerator, denominator = domain_matrix.solve_den(domain_rhs, method="rref")
    return numerator.to_Matrix() / domain.to_sympy(denominator)


def _algebraic_function_rank(matrix):
    symbols = sorted(
        set().union(*(entry.free_symbols for entry in matrix)), key=str
    )
    domain = sympy.QQ.algebraic_field(sympy.I, sympy.sqrt(2))
    if symbols:
        domain = domain.frac_field(*symbols)
    return DomainMatrix.from_Matrix(matrix).convert_to(domain).rank()


def _resolvent_blocks(l0, frequency, projector):
    """Describe the full 49-dimensional reduced resolvent by invariant blocks."""
    shifted = l0 - frequency * sympy.eye(N * N) + projector
    blind = range(3)
    outer = range(3, 7)
    index_blocks = []
    index_blocks.extend([[left * N + right for right in outer] for left in blind])
    index_blocks.extend([[left * N + right for left in outer] for right in blind])
    index_blocks.append([left * N + right for left in blind for right in blind])
    index_blocks.append([left * N + right for left in outer for right in outer])
    blocks = [
        {
            "indices": tuple(indices),
            "shifted": shifted.extract(indices, indices),
            "unshifted": (l0 - frequency * sympy.eye(N * N)).extract(indices, indices),
            "projector": projector.extract(indices, indices),
        }
        for indices in index_blocks
    ]
    for block in blocks:
        block["shifted_is_invertible"] = (
            _algebraic_function_rank(block["shifted"])
            == len(block["indices"])
        )
    return blocks


def _apply_reduced_resolvent(blocks, projector, source):
    """Apply G to source via exact local solves; skip identically zero blocks."""
    projected_source = (sympy.eye(N * N) - projector) * source
    image = sympy.zeros(N * N, source.cols)
    for block in blocks:
        indices = list(block["indices"])
        local_rhs = projected_source.extract(indices, range(source.cols))
        if local_rhs == sympy.zeros(len(indices), source.cols):
            continue
        local_image = _domain_solve(block["shifted"], local_rhs)
        for local_row, global_row in enumerate(indices):
            for column in range(source.cols):
                image[global_row, column] = local_image[local_row, column]
    return image


def exact_effective_operators(gamma):
    """Derive first/second Kato operators and their exact certificates."""
    gamma = _exact_scalar(gamma)
    h, eigenspaces, _, frequencies, _ = _peripheral_bases()

    # A sparse reflection-even Krylov basis completes the three derived blind
    # energy vectors without introducing nested trigonometric radicals.
    outer_eigenvectors = [
        sympy.eye(N).col(SEAT),
        (sympy.eye(N).col(2) + sympy.eye(N).col(4)) / sympy.sqrt(2),
        (sympy.eye(N).col(1) + sympy.eye(N).col(5)) / sympy.sqrt(2),
        (sympy.eye(N).col(0) + sympy.eye(N).col(6)) / sympy.sqrt(2),
    ]
    hilbert_basis = sympy.Matrix.hstack(
        *(vector for _, vector in eigenspaces), *outer_eigenvectors
    )
    if sympy.simplify(sympy.conjugate(hilbert_basis.T) * hilbert_basis) != sympy.eye(N):
        raise AssertionError("Hilbert basis must be exactly orthonormal")

    _, site_z, _, _ = exact_reduction_objects(sympy.Integer(0), gamma)
    energy_h = (
        sympy.conjugate(hilbert_basis.T) * h * hilbert_basis
    ).applyfunc(sympy.simplify)
    energy_z = (
        sympy.conjugate(hilbert_basis.T) * site_z * hilbert_basis
    ).applyfunc(sympy.simplify)
    l0 = _exact_generator_from_h_z(energy_h, energy_z, gamma)

    epsilon = sympy.Symbol("epsilon")
    symbolic_h, _, _, _ = exact_reduction_objects(epsilon, gamma)
    h1 = symbolic_h.diff(epsilon).subs(epsilon, 0)
    energy_h1 = (
        sympy.conjugate(hilbert_basis.T) * h1 * hilbert_basis
    ).applyfunc(sympy.simplify)
    l1 = _exact_generator_from_h_z(energy_h1, sympy.eye(N), sympy.Integer(0))

    first = {}
    second = {}
    polynomials = {"first": {}, "second": {}}
    certificates = {}
    for key, frequency in frequencies.items():
        vectors = []
        for left in range(3):
            for right_index in range(3):
                if sympy.simplify(
                    -sympy.I * (energy_h[left, left] - energy_h[right_index, right_index])
                    - frequency
                ) == 0:
                    vector = sympy.zeros(N * N, 1)
                    vector[left * N + right_index] = 1
                    vectors.append(vector)
        if key == "0":
            outer_identity = sympy.zeros(N * N, 1)
            for index in range(3, 7):
                outer_identity[index * N + index] = sympy.Rational(1, 2)
            vectors.append(outer_identity)
        right = sympy.Matrix.hstack(*vectors)
        left = right
        projector = right * sympy.conjugate(left.T)
        blocks = _resolvent_blocks(l0, frequency, projector)
        f1 = (sympy.conjugate(left.T) * l1 * right).applyfunc(sympy.simplify)
        resolvent_image = _apply_reduced_resolvent(
            blocks, projector, l1 * right
        )
        f2 = (
            -sympy.conjugate(left.T) * l1 * resolvent_image
        ).applyfunc(sympy.simplify)
        first[key] = f1
        second[key] = f2
        variable = sympy.Dummy("lambda")
        polynomials["first"][key] = f1.charpoly(variable)
        polynomials["second"][key] = f2.charpoly(variable)
        certificates[key] = {
            "generator": l0,
            "frequency": frequency,
            "right_basis": right,
            "left_basis": left,
            "projector": projector,
            "resolvent_blocks": blocks,
            "resolvent_image_of_l1_right": resolvent_image,
            "basis_dimension": sum(len(block["indices"]) for block in blocks),
        }
    return {
        "first": first,
        "second": second,
        "characteristic_polynomials": polynomials,
        "certificates": certificates,
    }


def _nonzero_row_basis(matrix):
    reduced, _ = matrix.rref()
    rows = [reduced.row(index) for index in range(reduced.rows) if any(reduced.row(index))]
    return sympy.Matrix.vstack(*rows) if rows else sympy.zeros(0, matrix.cols)


def _largest_invariant_subspace(operator, forbidden_rows):
    """Return intersection_k ker(forbidden_rows * operator**k), exactly."""
    constraints = _nonzero_row_basis(forbidden_rows)
    for _ in range(N * N):
        previous_rank = constraints.rows
        constraints = _nonzero_row_basis(
            sympy.Matrix.vstack(constraints, constraints * operator)
        )
        if constraints.rows == previous_rank:
            nullspace = constraints.nullspace()
            return (
                sympy.Matrix.hstack(*nullspace)
                if nullspace
                else sympy.zeros(N * N, 0)
            )
    raise AssertionError("invariant-subspace iteration did not stabilize")


def exact_punctured_kernel_certificate(epsilon, gamma):
    """Certify the rational punctured cases without a floating eigensolver."""
    epsilon = _exact_scalar(epsilon)
    gamma = _exact_scalar(gamma)
    if epsilon == 0:
        raise ValueError("punctured certificate requires nonzero epsilon")
    h, z, blind_vector, _ = exact_reduction_objects(epsilon, gamma)
    generator = exact_a_generator(epsilon, gamma)

    blind_vector = _normalized(blind_vector)
    stationary = [
        blind_vector * sympy.conjugate(blind_vector.T),
        sympy.eye(N) - blind_vector * sympy.conjugate(blind_vector.T),
    ]
    stationary_basis = sympy.Matrix.hstack(
        *[_row_stack_exact(operator) for operator in stationary]
    )
    if generator * stationary_basis != sympy.zeros(N * N, 2):
        raise AssertionError("displayed stationary operators must lie in ker A")

    commutator_columns = []
    dissipator_columns = []
    for coordinate in range(N * N):
        basis_matrix = sympy.zeros(N, N)
        basis_matrix[coordinate // N, coordinate % N] = 1
        commutator_columns.append(
            _row_stack_exact(-sympy.I * (h * basis_matrix - basis_matrix * h))
        )
        dissipator_columns.append(
            _row_stack_exact(z * basis_matrix * z - basis_matrix)
        )
    commutator = sympy.Matrix.hstack(*commutator_columns)
    dissipator = sympy.Matrix.hstack(*dissipator_columns)
    charged_coordinates = [
        index for index in range(N * N) if dissipator[index, index] != 0
    ]
    forbidden_rows = sympy.zeros(len(charged_coordinates), N * N)
    for row, coordinate in enumerate(charged_coordinates):
        forbidden_rows[row, coordinate] = 1
    invariant_basis = _largest_invariant_subspace(commutator, forbidden_rows)

    gram = sympy.conjugate(invariant_basis.T) * invariant_basis
    restricted = gram.inv() * sympy.conjugate(invariant_basis.T) * commutator * invariant_basis
    variable = sympy.Dummy("lambda")
    restricted_polynomial = restricted.charpoly(variable)
    kernel_dimension = N * N - _exact_rank(generator)
    stationary_rank = _exact_rank(stationary_basis)
    same_span = (
        _exact_rank(sympy.Matrix.hstack(invariant_basis, stationary_basis))
        == invariant_basis.cols
        == stationary_rank
    )
    if not same_span:
        raise AssertionError("zero-cost invariant space must equal the stationary span")
    nonzero_imaginary_axis_dimension = invariant_basis.cols - restricted.nullspace().__len__()
    return {
        "kernel_dimension": kernel_dimension,
        "stationary_rank": stationary_rank,
        "zero_cost_invariant_dimension": invariant_basis.cols,
        "nonzero_imaginary_axis_dimension": nonzero_imaginary_axis_dimension,
        "restricted_characteristic_polynomial": restricted_polynomial,
        "stationary_basis": stationary_basis,
        "zero_cost_invariant_basis": invariant_basis,
    }
