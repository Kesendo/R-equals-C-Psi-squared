"""Invariant reduction for the N=7 centre-watched missing-phase A block."""

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
import math
import os
from pathlib import Path
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy
from scipy.optimize import linear_sum_assignment
from sympy.polys.matrices import DomainMatrix


N = 7
SEAT = 3
RESULT_PATH = Path(__file__).resolve().parent / "results" / "missing_phase_relaxation_scale.json"


@dataclass(frozen=True)
class PrecisionDiagnostics:
    rank: int
    eigenvalue_error: float
    gap: float
    sep_complex: float
    sep_sylvester: float
    subspace_error: float
    projector_difference: float
    light: float
    light_difference: float
    max_principal_angle: float
    eigenvalue_pair_error: float


@dataclass(frozen=True)
class PrecisionDecision:
    status: str
    failures: tuple[str, ...]


@dataclass(frozen=True)
class PrecisionMatrix:
    matrix: mp.matrix
    bits: int
    token: int
    independently_reconstructed: bool
    exact_input_fingerprint: tuple
    source_build_token: int | None


@dataclass(frozen=True)
class ClusterState:
    epsilon_token: str
    rank: int
    right_basis: mp.matrix
    left_basis: mp.matrix
    right_ritz: tuple
    left_ritz: tuple


_build_counter = itertools.count(1)


def _mp_rational(token: str) -> mp.mpf:
    numerator, separator, denominator = str(token).partition("/")
    if not separator:
        return mp.mpf(numerator)
    return mp.mpf(int(numerator)) / int(denominator)


def _build_h_z_mp(epsilon_token: str, gamma_token: str, bits: int):
    with mp.workprec(bits):
        epsilon = _mp_rational(epsilon_token)
        gamma = _mp_rational(gamma_token)
        h = mp.zeros(N)
        bonds = (2 * (1 + epsilon), 2, 2, 2, 2, 2)
        for site, bond in enumerate(bonds):
            h[site, site + 1] = bond
            h[site + 1, site] = bond
        z = [mp.mpf(-1) if site == SEAT else mp.mpf(1) for site in range(N)]
        return h.copy(), tuple(z), +gamma


def _new_precision_matrix(
    matrixrezz,
    bits,
    independently_reconstructed=True,
    exact_input_fingerprint=(),
    source_build_token=None,
):
    return PrecisionMatrix(
        matrixrezz.copy(),
        bits,
        next(_build_counter),
        independently_reconstructed,
        tuple(exact_input_fingerprint),
        source_build_token,
    )


def build_k_mp(epsilon_token: str, gamma_token: str, bits: int) -> PrecisionMatrix:
    """Build K from rational tokens inside its own precision context."""
    with mp.workprec(bits):
        h, _, gamma = _build_h_z_mp(epsilon_token, gamma_token, bits)
        matrix = -mp.j * h
        matrix[SEAT, SEAT] -= 2 * gamma
        return _new_precision_matrix(matrix, bits, exact_input_fingerprint=("K", epsilon_token, gamma_token))


def _build_carrier_mp(epsilon_token, gamma_token, bits, carrier):
    with mp.workprec(bits):
        h, z, gamma = _build_h_z_mp(epsilon_token, gamma_token, bits)
        matrix = mp.zeros(N * N)
        for i in range(N):
            for j in range(N):
                row = i * N + j
                for k in range(N):
                    for ell in range(N):
                        column = k * N + ell
                        value = mp.mpc(0)
                        if j == ell:
                            value += (mp.j if carrier == "A_ADJOINT" else -mp.j) * h[i, k]
                        if i == k:
                            value += (-mp.j if carrier == "A_ADJOINT" else mp.j) * h[ell, j]
                        if i == k and j == ell:
                            if carrier in ("A", "A_ADJOINT"):
                                value += gamma * (z[i] * z[j] - 1)
                            elif carrier == "B":
                                value -= gamma * (z[i] * z[j] + 1)
                            else:
                                raise ValueError("unknown carrier")
                        matrix[row, column] = value
        return _new_precision_matrix(matrix, bits, exact_input_fingerprint=(carrier, epsilon_token, gamma_token))


def build_a_mp(epsilon_token: str, gamma_token: str, bits: int) -> PrecisionMatrix:
    """Build the full row-stack A action directly, without spectral lifting."""
    return _build_carrier_mp(epsilon_token, gamma_token, bits, "A")


def build_b_mp(epsilon_token: str, gamma_token: str, bits: int) -> PrecisionMatrix:
    """Build the full row-stack B action independently from h and z."""
    return _build_carrier_mp(epsilon_token, gamma_token, bits, "B")


def build_a_adjoint_mp(epsilon_token: str, gamma_token: str, bits: int) -> PrecisionMatrix:
    """Build A-dagger afresh from its action, rather than transposing stored A."""
    return _build_carrier_mp(epsilon_token, gamma_token, bits, "A_ADJOINT")


def build_global_flip_carrier_mp(epsilon_token, gamma_token, bits, carrier, copy):
    if carrier not in ("A", "B") or copy not in ("original", "global_flip"):
        raise ValueError("unsupported carrier copy")
    if carrier == "A":
        signs = (1, 1) if copy == "original" else (-1, -1)
    else:
        signs = (1, -1) if copy == "original" else (-1, 1)
    return _build_signed_carrier_mp(
        epsilon_token, gamma_token, bits, carrier, copy, *signs
    )


def _build_signed_carrier_mp(
    epsilon_token, gamma_token, bits, carrier, copy, left_sign, right_sign
):
    with mp.workprec(bits):
        h, z, gamma = _build_h_z_mp(epsilon_token, gamma_token, bits)
        matrix = mp.zeros(49)
        for i in range(N):
            for j in range(N):
                row = i * N + j
                for k in range(N):
                    for ell in range(N):
                        column = k * N + ell
                        value = mp.mpc(0)
                        if j == ell:
                            value -= mp.j * h[i, k]
                        if i == k:
                            value += mp.j * h[ell, j]
                        if i == k and j == ell:
                            value += gamma * (
                                left_sign * z[i] * right_sign * z[j] - 1
                            )
                        matrix[row, column] = value
        return _new_precision_matrix(
            matrix,
            bits,
            True,
            (carrier, copy, epsilon_token, gamma_token, left_sign, right_sign),
        )


def global_flip_copy_mutation_residuals(epsilon_token, gamma_token, bits):
    a_original = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, bits, "A", "original"
    )
    a_flip = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, bits, "A", "global_flip"
    )
    b_original = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, bits, "B", "original"
    )
    b_flip = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, bits, "B", "global_flip"
    )
    one_sided = _build_signed_carrier_mp(
        epsilon_token, gamma_token, bits, "B", "one_sided_mutant", -1, -1
    )

    def maximum_entry(left, right):
        return max(
            abs(left.matrix[i, j] - right.matrix[i, j])
            for i in range(49) for j in range(49)
        )

    return {
        "correct_A": maximum_entry(a_original, a_flip),
        "correct_B": maximum_entry(b_original, b_flip),
        "one_sided_B_flip": maximum_entry(b_flip, one_sided),
    }


def promote_precision_matrix(value: PrecisionMatrix, bits: int) -> PrecisionMatrix:
    with mp.workprec(bits):
        promoted = mp.matrix(
            [[mp.mpc(value.matrix[i, j]) for j in range(value.matrix.cols)]
             for i in range(value.matrix.rows)]
        )
        return _new_precision_matrix(
            promoted,
            bits,
            False,
            value.exact_input_fingerprint,
            value.token,
        )


def perturb_precision_matrix(value, row, column, amount_token):
    with mp.workprec(value.bits):
        matrix = value.matrix.copy()
        matrix[row, column] += _mp_rational(amount_token)
        return _new_precision_matrix(
            matrix,
            value.bits,
            False,
            value.exact_input_fingerprint,
            value.token,
        )


def precision_build_discrepancy(low, high):
    with mp.workprec(high.bits):
        promoted = mp.matrix(low.matrix)
        return spectral_norm_mp(promoted - high.matrix)


def precision_matrices_are_independent(low, high):
    return (
        low.bits * 2 == high.bits
        and low.token != high.token
        and low.matrix is not high.matrix
        and low.exact_input_fingerprint == high.exact_input_fingerprint
        and low.source_build_token is None
        and high.source_build_token is None
        and low.independently_reconstructed
        and high.independently_reconstructed
    )


def classify_precision_build_pair(low_diagnostic, high_diagnostic, low_build, high_build):
    return classify_precision_pair(
        low_diagnostic,
        high_diagnostic,
        independently_reconstructed=precision_matrices_are_independent(low_build, high_build),
    )


def right_eigen_residual_for_matrix(matrix, vector, eigenvalue):
    return spectral_norm_mp(matrix * vector - eigenvalue * vector)


def _sympy_matrix_to_mp(matrix, bits):
    with mp.workprec(bits):
        def convert(value):
            numerical = sympy.N(value, max(20, int(bits * 0.31) + 8))
            real, imaginary = numerical.as_real_imag()
            return mp.mpc(str(real), str(imaginary))

        return mp.matrix(
            [[convert(matrix[i, j]) for j in range(matrix.cols)]
             for i in range(matrix.rows)]
        )


def _modified_gram_schmidt(columns, bits):
    with mp.workprec(bits):
        accepted = []
        for source in columns:
            vector = mp.matrix(source)
            for _ in range(2):
                for basis in accepted:
                    vector -= basis * (basis.transpose_conj() * vector)[0]
            norm = mp.norm(vector)
            if norm > mp.power(2, -bits // 2):
                accepted.append(vector / norm)
        if not accepted:
            return mp.zeros(len(columns[0]), 0)
        return mp.matrix([[accepted[j][i] for j in range(len(accepted))]
                          for i in range(accepted[0].rows)])


def _smallest_singular_value(matrix):
    values = mp.svd(matrix, compute_uv=False)
    return min(values[index] for index in range(len(values)))


def _maximum_principal_angle(reference, candidate):
    sigma = min(mp.mpf(1), max(mp.mpf(0), _smallest_singular_value(reference.transpose_conj() * candidate)))
    return mp.acos(sigma)


def seed_cluster_at_zero(gamma_token: str, bits: int) -> ClusterState:
    """Seed the exact two-dimensional +2 sqrt(2)i Riesz subspace."""
    _ = gamma_token  # the zero-cost exact Riesz vectors are gamma-independent
    _, _, _, frequencies, operator_bases = _peripheral_bases()
    operators = operator_bases["+2sqrt2i"]
    with mp.workprec(bits):
        columns = [_sympy_matrix_to_mp(_row_stack_exact(item), bits) for item in operators]
        right = _modified_gram_schmidt(columns, bits)
        left = mp.matrix(right)
        value = mp.mpc(0, 2 * mp.sqrt(2))
        return ClusterState("0", 2, right, left, (value, value), (-value, -value))


def _eigenpairs(matrix):
    values, vectors = mp.eig(matrix, left=False, right=True)
    return list(values), [vectors[:, index] for index in range(vectors.cols)]


def _candidate_indices(values, centres, radius):
    return [
        index for index, value in enumerate(values)
        if any(abs(value - centre) <= radius for centre in centres)
    ]


def _pair_distance(left_values, right_values):
    target = [mp.conj(value) for value in right_values]
    straight = max(abs(left_values[i] - target[i]) for i in range(2))
    crossed = max(abs(left_values[1 - i] - target[i]) for i in range(2))
    if crossed < straight:
        return crossed, (1, 0)
    return straight, (0, 1)


def candidate_overlap_is_ambiguous(
    score_best, score_second, error_best, error_second, separation
):
    """Compare a dimensionless overlap uncertainty under any time rescaling."""
    if separation <= 0:
        return True
    return score_best - score_second <= (error_best + error_second) / separation


def track_cluster(
    matrix_a,
    matrix_a_adjoint,
    previous,
    epsilon_token,
    bits,
    candidate_radius,
):
    """Continue a rank-two cluster by maximum overlap, never eigenvalue sorting."""
    with mp.workprec(bits):
        right_values, right_vectors = _eigenpairs(matrix_a)
        right_pool = _candidate_indices(right_values, previous.right_ritz, candidate_radius)
        if len(right_pool) < 2 or len(right_pool) > 4:
            return None
        scored = []
        for pair in itertools.combinations(right_pool, 2):
            basis = _modified_gram_schmidt([right_vectors[i] for i in pair], bits)
            if basis.cols != 2:
                continue
            score = _smallest_singular_value(previous.right_basis.transpose_conj() * basis)
            residual = max(
                mp.norm(matrix_a * right_vectors[i] - right_values[i] * right_vectors[i])
                for i in pair
            )
            scored.append((score, residual, pair, basis))
        if not scored:
            return None
        scored.sort(key=lambda item: item[0], reverse=True)
        if len(scored) > 1 and candidate_overlap_is_ambiguous(
            scored[0][0],
            scored[1][0],
            scored[0][1],
            scored[1][1],
            candidate_radius,
        ):
            return None
        _, right_error, right_pair, right_basis = scored[0]
        chosen_right = tuple(right_values[i] for i in right_pair)

        left_values, left_vectors = _eigenpairs(matrix_a_adjoint)
        left_pool = _candidate_indices(
            left_values, tuple(mp.conj(value) for value in chosen_right), candidate_radius
        )
        if len(left_pool) < 2 or len(left_pool) > 4:
            return None
        left_choices = []
        for pair in itertools.combinations(left_pool, 2):
            distance, order = _pair_distance(
                [left_values[i] for i in pair], chosen_right
            )
            ordered = tuple(pair[i] for i in order)
            basis = _modified_gram_schmidt([left_vectors[i] for i in ordered], bits)
            if basis.cols == 2:
                residual = max(
                    mp.norm(matrix_a_adjoint * left_vectors[i] - left_values[i] * left_vectors[i])
                    for i in ordered
                )
                left_choices.append((distance, residual, ordered, basis))
        if not left_choices:
            return None
        left_choices.sort(key=lambda item: item[0])
        pair_error, left_error, left_pair, left_basis = left_choices[0]
        if pair_error > right_error + left_error + mp.power(2, -bits + 8):
            return None
        if _maximum_principal_angle(previous.right_basis, right_basis) >= mp.pi / 2:
            return None
        if _maximum_principal_angle(previous.left_basis, left_basis) >= mp.pi / 2:
            return None
        return ClusterState(
            str(epsilon_token),
            2,
            right_basis,
            left_basis,
            chosen_right,
            tuple(left_values[i] for i in left_pair),
        )


def ordered_block_form(matrix, cluster_basis):
    """Use deterministic Householder QR to put the selected cluster first."""
    bits = mp.mp.prec
    with mp.workprec(bits):
        columns = [cluster_basis[:, index] for index in range(cluster_basis.cols)]
        for coordinate in range(matrix.rows):
            standard = mp.zeros(matrix.rows, 1)
            standard[coordinate] = 1
            trial = _modified_gram_schmidt(columns + [standard], bits)
            if trial.cols > len(columns):
                columns.append(standard)
            if len(columns) == matrix.rows:
                break
        square = mp.matrix([[columns[j][i] for j in range(len(columns))]
                            for i in range(matrix.rows)])
        q, _ = mp.qr(square)
        q_cluster = q[:, :cluster_basis.cols]
        q_perp = q[:, cluster_basis.cols:]
        t_cluster = q_cluster.transpose_conj() * matrix * q_cluster
        t_perp = q_perp.transpose_conj() * matrix * q_perp
        lower_left = q_perp.transpose_conj() * matrix * q_cluster
        return t_cluster, t_perp, lower_left


def _mp_to_numpy(matrix):
    return np.array(
        [[complex(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)],
        dtype=complex,
    )


def ordinary_cluster_separation(matrix, state):
    """Read the complex-plane cluster/complement distance from the full spectrum."""
    values = np.linalg.eigvals(_mp_to_numpy(matrix))
    centres = np.array([complex(value) for value in state.right_ritz])
    cost = np.abs(values[:, None] - centres[None, :])
    rows, columns = linear_sum_assignment(cost)
    selected_rows = {int(rows[index]) for index in range(len(rows)) if columns[index] < 2}
    # rectangular assignment gives exactly two rows for the two cluster centres
    selected = values[sorted(selected_rows)]
    complement = np.delete(values, sorted(selected_rows))
    return mp.mpf(str(float(np.min(np.abs(selected[:, None] - complement[None, :])))))


def _fraction_token(value):
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def adaptive_track_from_zero(epsilon_token, gamma_token, bits, max_subdivisions=16):
    """Bridge deterministically from the exact seed without widening sep_C/3."""
    target = Fraction(epsilon_token)
    current = Fraction(0)
    state = seed_cluster_at_zero(gamma_token, bits)
    current_matrix = build_a_mp("0", gamma_token, bits).matrix
    bridge = ["0"]
    proposed = target
    subdivisions = 0
    while current != target:
        radius = ordinary_cluster_separation(current_matrix, state) / 3
        proposed_token = _fraction_token(proposed)
        proposed_matrix = build_a_mp(proposed_token, gamma_token, bits).matrix
        next_state = track_cluster(
            proposed_matrix,
            proposed_matrix.transpose_conj(),
            state,
            proposed_token,
            bits,
            radius,
        )
        if next_state is None:
            subdivisions += 1
            if subdivisions > max_subdivisions:
                return None, tuple(bridge)
            proposed = (current + proposed) / 2
            continue
        current = proposed
        state = next_state
        current_matrix = proposed_matrix
        bridge.append(proposed_token)
        proposed = target
    return state, tuple(bridge)


def _blind_vector_mp(epsilon_token, bits):
    with mp.workprec(bits):
        r = 1 + _mp_rational(epsilon_token)
        vector = mp.matrix([1, 0, -r, 0, r, 0, -r])
        return vector / mp.norm(vector)


def _nearest_eigenpair(matrix, target):
    values, vectors = _eigenpairs(matrix)
    index = min(range(len(values)), key=lambda item: abs(values[item] - target))
    vector = vectors[index] / mp.norm(vectors[index])
    return values[index], vector


def _row_stack_outer_mp(left, right):
    return mp.matrix(
        [left[i] * mp.conj(right[j]) for i in range(N) for j in range(N)]
    )


def reduced_cluster_state(epsilon_token, gamma_token, bits):
    """Construct the doubled 49D cluster from independently solved 7D factors."""
    with mp.workprec(bits):
        k = build_k_mp(epsilon_token, gamma_token, bits).matrix
        kh = k.transpose_conj()
        target = mp.mpc(0, 2 * mp.sqrt(2))
        plus, u_plus = _nearest_eigenpair(k, target)
        minus, u_minus = _nearest_eigenpair(k, -target)
        left_plus, w_plus = _nearest_eigenpair(kh, mp.conj(plus))
        _, w_minus = _nearest_eigenpair(kh, mp.conj(minus))
        blind = _blind_vector_mp(epsilon_token, bits)
        right = _modified_gram_schmidt(
            [_row_stack_outer_mp(u_plus, blind), _row_stack_outer_mp(blind, u_minus)], bits
        )
        left = _modified_gram_schmidt(
            [_row_stack_outer_mp(w_plus, blind), _row_stack_outer_mp(blind, w_minus)], bits
        )
        return ClusterState(
            str(epsilon_token), 2, right, left,
            (plus, mp.conj(minus)),
            (left_plus, left_plus),
        )


def _mp_liouville_blocks(matrix, epsilon_token, bits):
    with mp.workprec(bits):
        blind = _blind_vector_mp(epsilon_token, bits)
        hilbert = _modified_gram_schmidt(
            [blind] + [mp.eye(N)[:, index] for index in range(N)], bits
        )
        pairs = (
            [(0, 0)]
            + [(i, 0) for i in range(1, N)]
            + [(0, j) for j in range(1, N)]
            + [(i, j) for i in range(1, N) for j in range(1, N)]
        )
        unitary = mp.matrix(
            [[_row_stack_outer_mp(hilbert[:, i], hilbert[:, j])[row]
              for i, j in pairs] for row in range(49)]
        )
        transformed = unitary.transpose_conj() * matrix * unitary
        blocks = (
            transformed[0:1, 0:1],
            transformed[1:7, 1:7],
            transformed[7:13, 7:13],
            transformed[13:49, 13:49],
        )
        return unitary, blocks


def _cross_block_candidates(matrix, epsilon_token, bits, centres, radius):
    unitary, blocks = _mp_liouville_blocks(matrix, epsilon_token, bits)
    values = []
    vectors = []
    for block_index, start in ((1, 1), (2, 7)):
        local_values, local_vectors = _eigenpairs(blocks[block_index])
        for index, value in enumerate(local_values):
            if any(abs(value - centre) <= radius for centre in centres):
                values.append(value)
                vectors.append(unitary[:, start:start + 6] * local_vectors[index])
    return values, vectors


def _local_complement(block, vector, bits):
    with mp.workprec(bits):
        q = _modified_gram_schmidt(
            [vector] + [mp.eye(block.rows)[:, index] for index in range(block.rows)], bits
        )
        q_perp = q[:, 1:]
        complement = q_perp.transpose_conj() * block * q_perp
        lower = q_perp.transpose_conj() * block * q[:, 0:1]
        return complement, lower


def native_block_separation(matrix, state, epsilon_token, bits, *, compute_complex=True):
    """Measure complement separation via the exact 1+6+6+36 reduction."""
    with mp.workprec(bits):
        unitary, blocks = _mp_liouville_blocks(matrix, epsilon_token, bits)
        coordinates = unitary.transpose_conj() * state.right_basis
        cross_complements = []
        lower_residuals = []
        for block, start in ((blocks[1], 1), (blocks[2], 7)):
            local = coordinates[start:start + 6, :]
            norms = [mp.norm(local[:, column]) for column in range(2)]
            column = max(range(2), key=lambda index: norms[index])
            complement, lower = _local_complement(block, local[:, column], bits)
            cross_complements.append(complement)
            lower_residuals.append(spectral_norm_mp(lower))
        complement_blocks = [blocks[0], *cross_complements, blocks[3]]
        sep_complex = None
        if compute_complex:
            complement_values = []
            for block in complement_blocks:
                values = mp.eig(block, left=False, right=False)
                if isinstance(values, tuple):
                    values = values[0]
                complement_values.extend(values)
            sep_complex = min(
                abs(cluster_value - complement_value)
                for cluster_value in state.right_ritz
                for complement_value in complement_values
            )
        centre = sum(state.right_ritz) / 2
        cluster_spread = max(abs(value - centre) for value in state.right_ritz)
        sylvester_bounds = []
        for block in complement_blocks:
            sigma = min(mp.svd(block - centre * mp.eye(block.rows), compute_uv=False))
            sylvester_bounds.append(max(mp.mpf(0), sigma - cluster_spread))
        sep_sylvester = min(sylvester_bounds)
        transformed = unitary.transpose_conj() * matrix * unitary
        block_diagonal = mp.zeros(49)
        for start, block in zip((0, 1, 7, 13), blocks):
            block_diagonal[start:start + block.rows, start:start + block.cols] = block
        off_block = transformed - block_diagonal
        off_block_residual = max(abs(off_block[i, j]) for i in range(49) for j in range(49))
        return {
            "block_dimensions": [block.rows for block in complement_blocks],
            "sep_complex": sep_complex,
            "sep_sylvester": sep_sylvester,
            "off_block_residual": off_block_residual,
            "lower_left_residual": max(lower_residuals),
        }


def _choose_overlap_candidate(matrix, values, vectors, previous_basis, bits, radius):
    scored = []
    for pair in itertools.combinations(range(len(values)), 2):
        basis = _modified_gram_schmidt([vectors[index] for index in pair], bits)
        if basis.cols != 2:
            continue
        score = _smallest_singular_value(previous_basis.transpose_conj() * basis)
        residual = spectral_norm_mp(
            (mp.eye(49) - _orthogonal_projector(basis)) * matrix * basis
        )
        scored.append((score, residual, pair, basis))
    scored.sort(key=lambda item: item[0], reverse=True)
    if not scored:
        return None
    second_score = scored[1][0] if len(scored) > 1 else mp.mpf(0)
    best = scored[0]
    uncertainty = (best[1] + (scored[1][1] if len(scored) > 1 else 0)) / radius
    if best[0] - second_score <= uncertainty:
        return None
    return best, second_score, uncertainty


def factor_track_cluster(matrix, previous, epsilon_token, gamma_token, bits, candidate_radius):
    """Solve the two 6D cross blocks, then validate their lifted 49D cluster."""
    with mp.workprec(bits):
        full_values = np.linalg.eigvals(_mp_to_numpy(matrix))
        previous_centres = np.array([complex(value) for value in previous.right_ritz])
        mask = np.min(np.abs(full_values[:, None] - previous_centres[None, :]), axis=1) <= float(candidate_radius)
        pool = np.flatnonzero(mask)
        if len(pool) < 2 or len(pool) > 4:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": int(len(pool)),
                "status": "UNRESOLVED",
                "failures": ["candidate union did not contain between two and four eigenvalues"],
            }
        if len(pool) == 2:
            candidate = reduced_cluster_state(epsilon_token, gamma_token, bits)
            rows, columns = linear_sum_assignment(
                np.abs(
                    full_values[pool, None]
                    - np.array([complex(value) for value in candidate.right_ritz])[None, :]
                )
            )
            match_error = float(np.max(np.abs(
                full_values[pool[rows]]
                - np.array([complex(candidate.right_ritz[index]) for index in columns])
            )))
            adjoint_build = build_a_adjoint_mp(epsilon_token, gamma_token, bits)
            right_residual = spectral_norm_mp(
                (mp.eye(49) - _orthogonal_projector(candidate.right_basis))
                * matrix * candidate.right_basis
            )
            left_residual = spectral_norm_mp(
                (mp.eye(49) - _orthogonal_projector(candidate.left_basis))
                * adjoint_build.matrix * candidate.left_basis
            )
            score = _smallest_singular_value(
                previous.right_basis.transpose_conj() * candidate.right_basis
            )
            left_score = _smallest_singular_value(
                previous.left_basis.transpose_conj() * candidate.left_basis
            )
            margin = min(score, left_score)
            uncertainty = max(right_residual, left_residual) / candidate_radius
            pair_error, _ = _pair_distance(candidate.left_ritz, candidate.right_ritz)
            pair_bound = right_residual + left_residual + mp.power(2, -bits + 8) * mp.norm(matrix)
            right_angle = _maximum_principal_angle(previous.right_basis, candidate.right_basis)
            left_angle = _maximum_principal_angle(previous.left_basis, candidate.left_basis)
            failures = []
            if match_error > 1024 * np.finfo(float).eps * max(1.0, np.linalg.norm(_mp_to_numpy(matrix), 2)):
                failures.append("factor Ritz values do not match the two full49 candidates")
            if margin <= uncertainty:
                failures.append("dimensionless overlap margin is not resolved")
            if pair_error > pair_bound:
                failures.append("independent left/right Ritz pairing exceeds residual bounds")
            if max(right_angle, left_angle) >= mp.pi / 2:
                failures.append("principal angle reached pi/2")
            record = {
                "epsilon_token": epsilon_token,
                "candidate_count": 2,
                "candidate_radius": float(candidate_radius),
                "right_overlap_score": float(score),
                "left_overlap_score": float(left_score),
                "overlap_margin": float(margin),
                "dimensionless_uncertainty": float(uncertainty),
                "right_principal_angle": float(right_angle),
                "left_principal_angle": float(left_angle),
                "left_right_ritz_pair_error": float(pair_error),
                "full49_candidate_match_error": match_error,
                "full49_invariance_residual": float(max(right_residual, left_residual)),
                "fresh_adjoint_build_token": adjoint_build.token,
                "status": "UNRESOLVED" if failures else "TRACKED",
                "failures": failures,
            }
            return (None if failures else candidate), record
        right_values, right_vectors = _cross_block_candidates(
            matrix, epsilon_token, bits, previous.right_ritz, candidate_radius
        )
        chosen_right = _choose_overlap_candidate(
            matrix, right_values, right_vectors, previous.right_basis, bits, candidate_radius
        )
        if chosen_right is None:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": int(len(pool)),
                "status": "UNRESOLVED",
                "failures": ["right cross-block overlap selection is ambiguous"],
            }
        (score, right_residual, right_pair, right_basis), second_score, right_uncertainty = chosen_right
        right_ritz = tuple(right_values[index] for index in right_pair)

        adjoint_build = build_a_adjoint_mp(epsilon_token, gamma_token, bits)
        left_values, left_vectors = _cross_block_candidates(
            adjoint_build.matrix,
            epsilon_token,
            bits,
            tuple(mp.conj(value) for value in right_ritz),
            candidate_radius,
        )
        chosen_left = _choose_overlap_candidate(
            adjoint_build.matrix,
            left_values,
            left_vectors,
            previous.left_basis,
            bits,
            candidate_radius,
        )
        if chosen_left is None:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": int(len(pool)),
                "status": "UNRESOLVED",
                "failures": ["left cross-block overlap selection is ambiguous"],
            }
        (left_score, left_residual, left_pair, left_basis), left_second, left_uncertainty = chosen_left
        left_ritz_unordered = tuple(left_values[index] for index in left_pair)
        distance, order = _pair_distance(left_ritz_unordered, right_ritz)
        left_ritz = tuple(left_ritz_unordered[index] for index in order)
        left_basis = mp.matrix([[left_basis[row, order[column]] for column in range(2)] for row in range(49)])
        pair_bound = right_residual + left_residual + mp.power(2, -bits + 8) * mp.norm(matrix)
        if distance > pair_bound:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": int(len(pool)),
                "status": "UNRESOLVED",
                "failures": ["independent left/right Ritz pairing exceeds residual bounds"],
            }
        candidate = ClusterState(
            epsilon_token, 2, right_basis, left_basis, right_ritz, left_ritz
        )
        margin = min(score - second_score, left_score - left_second)
        uncertainty = max(right_uncertainty, left_uncertainty)
        right_angle = _maximum_principal_angle(previous.right_basis, right_basis)
        left_angle = _maximum_principal_angle(previous.left_basis, left_basis)
        if max(right_angle, left_angle) >= mp.pi / 2:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": int(len(pool)),
                "status": "UNRESOLVED",
                "failures": ["principal angle reached pi/2"],
            }
        return candidate, {
            "epsilon_token": epsilon_token,
            "candidate_count": int(len(pool)),
            "candidate_radius": float(candidate_radius),
            "right_overlap_score": float(score),
            "left_overlap_score": float(left_score),
            "overlap_margin": float(margin),
            "dimensionless_uncertainty": float(uncertainty),
            "right_principal_angle": float(right_angle),
            "left_principal_angle": float(left_angle),
            "left_right_ritz_pair_error": float(distance),
            "full49_invariance_residual": float(max(right_residual, left_residual)),
            "fresh_adjoint_build_token": adjoint_build.token,
            "status": "TRACKED",
            "failures": [],
        }


def prepare_continuation_paths(paths, gamma_token, bits):
    records = {}
    states = {}
    for branch, tokens in paths.items():
        previous = seed_cluster_at_zero(gamma_token, bits)
        previous_matrix = build_a_mp("0", gamma_token, bits).matrix
        steps = []
        for epsilon_token in tokens[1:]:
            radius = ordinary_cluster_separation(previous_matrix, previous) / 3
            matrix = build_a_mp(epsilon_token, gamma_token, bits).matrix
            current, record = factor_track_cluster(
                matrix, previous, epsilon_token, gamma_token, bits, radius
            )
            steps.append(record)
            if current is None:
                raise RuntimeError(
                    f"{branch} continuation unresolved at {epsilon_token}: {record['failures']}"
                )
            previous = current
            previous_matrix = matrix
            states[epsilon_token] = current
        records[branch] = steps
    return records, states


def validate_continuation_paths(paths, gamma_token, bits):
    return prepare_continuation_paths(paths, gamma_token, bits)[0]


def _orthogonal_projector(basis):
    return basis * basis.transpose_conj()


def _centre_light(projector):
    total = mp.mpc(0)
    for i in range(N):
        for j in range(N):
            if (i == SEAT) != (j == SEAT):
                total += projector[i * N + j, i * N + j]
    return mp.re(total) / 2


def _diagnostic_at_precision(
    build,
    state,
    build_error,
    projector_difference,
    light_difference,
    pair_error,
    sep_complex,
    sep_sylvester,
    max_principal_angle,
):
    with mp.workprec(build.bits):
        matrix = build.matrix
        right = state.right_basis
        left = state.left_basis
        pr = _orthogonal_projector(right)
        pl = _orthogonal_projector(left)
        identity = mp.eye(N * N)
        rr = (identity - pr) * matrix * right
        rl = (identity - pl) * matrix.transpose_conj() * left
        sigma = _smallest_singular_value(left.transpose_conj() * right)
        kappa = 1 / sigma
        residual = max(spectral_norm_mp(rr), spectral_norm_mp(rl))
        eigen_error = kappa * (residual + build_error + mp.eps * mp.norm(matrix) * 256)
        gap = -mp.re(sum(state.right_ritz) / 2)
        sep_complex = mp.mpf(sep_complex)
        sep_sylvester = mp.mpf(sep_sylvester)
        subspace_error = eigen_error / sep_sylvester
        light = min(_centre_light(pr), _centre_light(pl))
        return PrecisionDiagnostics(
            2,
            float(eigen_error),
            float(gap),
            float(sep_complex),
            float(sep_sylvester),
            float(subspace_error),
            float(projector_difference),
            float(light),
            float(light_difference),
            float(max_principal_angle),
            float(pair_error),
        )


def _pair_ritz_error(low, high):
    costs = [
        max(abs(low.right_ritz[i] - high.right_ritz[i]) for i in range(2)),
        max(abs(low.right_ritz[i] - high.right_ritz[1 - i]) for i in range(2)),
    ]
    return min(costs)


def _full_b_spectrum_check(a_build, b_build, gamma_token):
    a = _mp_to_numpy(a_build.matrix)
    b = _mp_to_numpy(b_build.matrix)
    av, ar = np.linalg.eig(a)
    bv, br = np.linalg.eig(b)
    gamma = float(Fraction(gamma_token))
    expected = -2 * gamma - np.conj(av)
    rows, columns = linear_sum_assignment(np.abs(bv[:, None] - expected[None, :]))
    residual = float(np.max(np.abs(bv[rows] - expected[columns])))
    a_backward = max(np.linalg.norm(a @ ar[:, i] - av[i] * ar[:, i]) for i in range(49))
    b_backward = max(np.linalg.norm(b @ br[:, i] - bv[i] * br[:, i]) for i in range(49))
    rounding = 2048 * np.finfo(float).eps * max(1.0, np.linalg.norm(a, 2), np.linalg.norm(b, 2))
    return residual, float(a_backward + b_backward + rounding)


def carrier_block_decomposition(build, epsilon_token):
    """Expose the exact 1+6 Hilbert split as 1+6+6+36 Liouville blocks."""
    with mp.workprec(build.bits):
        unitary, blocks = _mp_liouville_blocks(
            build.matrix, epsilon_token, build.bits
        )
        transformed = unitary.transpose_conj() * build.matrix * unitary
        block_diagonal = mp.zeros(N * N)
        starts = (0, 1, 7, 13)
        for start, block in zip(starts, blocks):
            block_diagonal[start:start + block.rows, start:start + block.cols] = block
        unitarity = unitary.transpose_conj() * unitary - mp.eye(N * N)
        unitarity_residual = max(abs(unitarity[i, j]) for i in range(49) for j in range(49))
        off_block = transformed - block_diagonal
        off_block_residual = max(abs(off_block[i, j]) for i in range(49) for j in range(49))

    full = _mp_to_numpy(build.matrix)
    full_values, full_vectors = np.linalg.eig(full)
    union_values = np.concatenate([np.linalg.eigvals(_mp_to_numpy(block)) for block in blocks])
    rows, columns = linear_sum_assignment(np.abs(full_values[:, None] - union_values[None, :]))
    union_residual = float(np.max(np.abs(full_values[rows] - union_values[columns])))
    eigen_residual = float(max(
        np.linalg.norm(full @ full_vectors[:, i] - full_values[i] * full_vectors[:, i])
        for i in range(49)
    ))
    return {
        "dimensions": [block.rows for block in blocks],
        "blocks": tuple(blocks),
        "unitarity_residual": float(unitarity_residual),
        "off_block_residual": float(off_block_residual),
        "full_spectrum_union_residual": union_residual,
        "full_eigen_residual": eigen_residual,
    }


def measure_precision_pair(
    epsilon_token,
    gamma_token,
    low_bits,
    high_bits,
    *,
    low_state=None,
    high_state=None,
):
    """Measure a fresh p/2p doubled cluster and all trust-contract quantities."""
    low_build = build_a_mp(epsilon_token, gamma_token, low_bits)
    high_build = build_a_mp(epsilon_token, gamma_token, high_bits)
    low_state = low_state or reduced_cluster_state(epsilon_token, gamma_token, low_bits)
    high_state = high_state or reduced_cluster_state(epsilon_token, gamma_token, high_bits)
    low_pr = _orthogonal_projector(low_state.right_basis)
    high_pr = _orthogonal_projector(high_state.right_basis)
    low_pl = _orthogonal_projector(low_state.left_basis)
    high_pl = _orthogonal_projector(high_state.left_basis)
    # Extra comparison precision is essential here: projector distance is the
    # square root of 1-sigma_min^2, so evaluating it at only 2p can erase a
    # genuine O(2^-p) angle by rounding sigma_min to one.
    with mp.workprec(high_bits + 64):
        promoted_low_right = _modified_gram_schmidt(
            [mp.matrix(low_state.right_basis[:, index]) for index in range(2)], high_bits
        )
        promoted_low_left = _modified_gram_schmidt(
            [mp.matrix(low_state.left_basis[:, index]) for index in range(2)], high_bits
        )
        right_sigma = _smallest_singular_value(
            promoted_low_right.transpose_conj() * high_state.right_basis
        )
        left_sigma = _smallest_singular_value(
            promoted_low_left.transpose_conj() * high_state.left_basis
        )
        right_difference = mp.sqrt(max(mp.mpf(0), 1 - right_sigma**2))
        left_difference = mp.sqrt(max(mp.mpf(0), 1 - left_sigma**2))
        low_right_light = _centre_light(mp.matrix(low_pr))
        high_right_light = _centre_light(high_pr)
        low_left_light = _centre_light(mp.matrix(low_pl))
        high_left_light = _centre_light(high_pl)
        right_light_difference = abs(low_right_light - high_right_light)
        left_light_difference = abs(low_left_light - high_left_light)
        pair_error = _pair_ritz_error(low_state, high_state)
        build_error = precision_build_discrepancy(low_build, high_build)
        right_angle = mp.acos(min(mp.mpf(1), max(mp.mpf(0), right_sigma)))
        left_angle = mp.acos(min(mp.mpf(1), max(mp.mpf(0), left_sigma)))
        precision_angle = max(right_angle, left_angle)
        overlap_sigma = _smallest_singular_value(
            high_state.left_basis.transpose_conj() * high_state.right_basis
        )
        cluster_kappa = 1 / overlap_sigma
    right_separation = native_block_separation(
        high_build.matrix, high_state, epsilon_token, high_bits
    )
    high_adjoint = build_a_adjoint_mp(epsilon_token, gamma_token, high_bits)
    high_left_state = ClusterState(
        epsilon_token,
        2,
        high_state.left_basis,
        high_state.right_basis,
        high_state.left_ritz,
        high_state.right_ritz,
    )
    left_separation = native_block_separation(
        high_adjoint.matrix,
        high_left_state,
        epsilon_token,
        high_bits,
        compute_complex=False,
    )
    high_sep_complex = right_separation["sep_complex"]
    high_sep_sylvester = min(
        right_separation["sep_sylvester"], left_separation["sep_sylvester"]
    )
    low_sep_complex = max(
        mp.mpf("1e-300"), high_sep_complex - 2 * cluster_kappa * build_error
    )
    low_sep_sylvester = max(
        mp.mpf("1e-300"), high_sep_sylvester - build_error
    )
    projector_difference = max(right_difference, left_difference)
    light_difference = max(right_light_difference, left_light_difference)
    low_diagnostic = _diagnostic_at_precision(
        low_build,
        low_state,
        build_error,
        projector_difference,
        light_difference,
        pair_error,
        low_sep_complex,
        low_sep_sylvester,
        precision_angle,
    )
    high_diagnostic = _diagnostic_at_precision(
        high_build,
        high_state,
        build_error,
        projector_difference,
        light_difference,
        pair_error,
        high_sep_complex,
        high_sep_sylvester,
        precision_angle,
    )
    status = classify_precision_build_pair(
        low_diagnostic, high_diagnostic, low_build, high_build
    )
    gamma = mp.mpf(Fraction(gamma_token).numerator) / Fraction(gamma_token).denominator
    gap = -mp.re(sum(high_state.right_ritz) / 2)
    lights = {
        "right": float(high_right_light),
        "left": float(high_left_light),
    }
    absorption = {
        side: float(gap - 2 * gamma * value) for side, value in
        (("right", high_right_light), ("left", high_left_light))
    }
    s = high_state.left_basis.transpose_conj() * high_state.right_basis
    biorthogonal = (
        high_state.right_basis * s**-1 * high_state.left_basis.transpose_conj()
    )
    low_b_build = build_b_mp(epsilon_token, gamma_token, low_bits)
    b_build = build_b_mp(epsilon_token, gamma_token, high_bits)
    b_residual, b_error = _full_b_spectrum_check(high_build, b_build, gamma_token)
    with mp.workprec(high_bits):
        gamma_high = _mp_rational(gamma_token)
        operator_delta = (
            b_build.matrix + 2 * gamma_high * mp.eye(49)
            + high_build.matrix.transpose_conj()
        )
        operator_map_residual = max(
            abs(operator_delta[i, j]) for i in range(49) for j in range(49)
        )
        operator_map_bound = (
            precision_build_discrepancy(low_build, high_build)
            + precision_build_discrepancy(low_b_build, b_build)
        )
    a_original = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, high_bits, "A", "original"
    )
    a_copy = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, high_bits, "A", "global_flip"
    )
    b_original = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, high_bits, "B", "original"
    )
    b_copy = build_global_flip_carrier_mp(
        epsilon_token, gamma_token, high_bits, "B", "global_flip"
    )
    copy_residual = max(
        abs(a_original.matrix[i, j] - a_copy.matrix[i, j])
        for i in range(49) for j in range(49)
    )
    copy_residual = max(
        copy_residual,
        max(abs(b_original.matrix[i, j] - b_copy.matrix[i, j])
            for i in range(49) for j in range(49)),
    )
    return {
        "status": status,
        "low": low_diagnostic,
        "high": high_diagnostic,
        "lights": lights,
        "unnormalized_lights": {side: 2 * value for side, value in lights.items()},
        "absorption_residuals": absorption,
        "absorption_error": float(high_diagnostic.eigenvalue_error + 2 * float(gamma) * light_difference),
        "biorthogonal_light": float(_centre_light(biorthogonal)),
        "right_projector_difference": float(right_difference),
        "left_projector_difference": float(left_difference),
        "right_light_difference": float(right_light_difference),
        "left_light_difference": float(left_light_difference),
        "full_b_match_residual": b_residual,
        "full_b_match_error": b_error,
        "b_operator_map_residual": float(operator_map_residual),
        "b_operator_map_bound": float(operator_map_bound),
        "global_flip_copy_residual": float(copy_residual),
        "low_build": low_build,
        "high_build": high_build,
        "state": high_state,
        "native_block_dimensions": right_separation["block_dimensions"],
        "right_native_off_block_residual": float(right_separation["off_block_residual"]),
        "left_native_off_block_residual": float(left_separation["off_block_residual"]),
        "right_native_lower_left_residual": float(right_separation["lower_left_residual"]),
        "left_native_lower_left_residual": float(left_separation["lower_left_residual"]),
    }


def _diagnostic_from_measured_mutation(*, sep_sylvester, projector_difference=0.0, light_difference=0.0):
    eigenvalue_error = max(1e-16, 0.3 * float(sep_sylvester))
    return PrecisionDiagnostics(
        rank=2,
        eigenvalue_error=eigenvalue_error,
        gap=1.0,
        sep_complex=1.0,
        sep_sylvester=float(sep_sylvester),
        subspace_error=eigenvalue_error / float(sep_sylvester),
        projector_difference=float(projector_difference),
        light=0.25,
        light_difference=float(light_difference),
        max_principal_angle=0.1,
        eigenvalue_pair_error=0.0,
    )


def fixed_kernel_cutoff_mutation(gamma_token, *, cutoff, maximum_power):
    gamma = float(Fraction(gamma_token))
    measurements = []
    for power in range(3, maximum_power + 1):
        epsilon = 2.0**-power
        gap = tracked_reduced_gap(epsilon, gamma)
        measurements.append((power, gap))
    misclassified = any(0 < gap < cutoff for _, gap in measurements)
    measured_rank = 4 if misclassified else 2
    base = _diagnostic_from_measured_mutation(sep_sylvester=0.1)
    diagnostic = replace_diagnostic(base, rank=measured_rank)
    decision = decide_precision_pair(
        diagnostic, diagnostic, independently_reconstructed=True
    )
    return {
        "measurements": measurements,
        "misclassified": misclassified,
        "status": decision.status,
        "failures": decision.failures,
    }


def nonnormal_sylvester_mutation():
    t_cluster = np.diag([0.0, 3.0]).astype(complex)
    t_perp = np.array([[1.0, 1e8], [0.0, 2.0]], dtype=complex)
    sep_complex = min(
        abs(a - b) for a in np.linalg.eigvals(t_cluster) for b in np.linalg.eigvals(t_perp)
    )
    sylvester = np.kron(np.eye(2), t_cluster) - np.kron(t_perp.T, np.eye(2))
    sep_sylvester = float(np.linalg.svd(sylvester, compute_uv=False)[-1])
    diagnostic = _diagnostic_from_measured_mutation(sep_sylvester=sep_sylvester)
    decision = decide_precision_pair(
        diagnostic, diagnostic, independently_reconstructed=True
    )
    return {
        "sep_complex": float(sep_complex),
        "sep_sylvester": sep_sylvester,
        "status": decision.status,
        "failures": decision.failures,
    }


def equal_light_rotated_projector_mutation():
    q1 = np.column_stack((np.eye(4)[:, 0], np.eye(4)[:, 2]))
    q2 = np.column_stack((np.eye(4)[:, 1], np.eye(4)[:, 2]))
    p1 = q1 @ q1.conj().T
    p2 = q2 @ q2.conj().T
    delta = np.diag([1.0, 1.0, 0.0, 0.0])
    light1 = np.trace(p1 @ delta).real / 2
    light2 = np.trace(p2 @ delta).real / 2
    projector_difference = float(np.linalg.norm(p1 - p2, 2))
    light_difference = float(abs(light1 - light2))
    diagnostic = _diagnostic_from_measured_mutation(
        sep_sylvester=0.1,
        projector_difference=projector_difference,
        light_difference=light_difference,
    )
    decision = decide_precision_pair(
        diagnostic, diagnostic, independently_reconstructed=True
    )
    return {
        "light_difference": light_difference,
        "projector_difference": projector_difference,
        "status": decision.status,
        "failures": decision.failures,
    }


PRECISION_PAIRS = ((53, 106), (106, 212), (212, 424))
GAMMA_TOKENS = ("3/100", "3/10", "3")
EPSILON_TOKENS = (
    "1/10", "-1/10",
    *(token for power in range(3, 14) for token in (f"1/{2**power}", f"-1/{2**power}")),
)


def _decimal(value):
    return format(float(value), ".17g")


def _serialize_diagnostic(value):
    return {
        field: (value.rank if field == "rank" else _decimal(getattr(value, field)))
        for field in PrecisionDiagnostics.__dataclass_fields__
    }


def _row_from_measurement(epsilon_token, gamma_token, low_bits, high_bits, result, attempts):
    state = result["state"]
    k = build_k_mp(epsilon_token, gamma_token, high_bits).matrix
    target = mp.mpc(0, 2 * mp.sqrt(2))
    k_value, _ = _nearest_eigenpair(k, target)
    a_values = np.linalg.eigvals(_mp_to_numpy(result["high_build"].matrix))
    gap = result["high"].gap
    error = result["high"].eigenvalue_error
    positive = sorted(-value.real for value in a_values if -value.real > max(1e-15, gap + 4 * error))
    b_values = np.linalg.eigvals(_mp_to_numpy(build_b_mp(epsilon_token, gamma_token, high_bits).matrix))
    return {
        "epsilon_token": epsilon_token,
        "gamma_token": gamma_token,
        "precision_bits": [low_bits, high_bits],
        "K_gap": _decimal(-mp.re(k_value)),
        "A_gap": _decimal(gap),
        "cluster_rank": state.rank,
        "right_light": _decimal(result["lights"]["right"]),
        "left_light": _decimal(result["lights"]["left"]),
        "right_absorption_residual": _decimal(result["absorption_residuals"]["right"]),
        "left_absorption_residual": _decimal(result["absorption_residuals"]["left"]),
        "absorption_error": _decimal(result["absorption_error"]),
        "low_diagnostics": _serialize_diagnostic(result["low"]),
        "high_diagnostics": _serialize_diagnostic(result["high"]),
        "right_projector_difference": _decimal(result["right_projector_difference"]),
        "left_projector_difference": _decimal(result["left_projector_difference"]),
        "right_light_difference": _decimal(result["right_light_difference"]),
        "left_light_difference": _decimal(result["left_light_difference"]),
        "next_distinct_A_rate": _decimal(positive[0]),
        "B_gap": _decimal(min(-b_values.real)),
        "B_operator_map_residual": _decimal(result["b_operator_map_residual"]),
        "B_operator_map_bound": _decimal(result["b_operator_map_bound"]),
        "full_B_match_numerical_read": _decimal(result["full_b_match_residual"]),
        "full_B_match_backward_error_read": _decimal(result["full_b_match_error"]),
        "global_flip_copy_residual": _decimal(result["global_flip_copy_residual"]),
        "native_block_dimensions": result["native_block_dimensions"],
        "right_native_off_block_residual": _decimal(result["right_native_off_block_residual"]),
        "left_native_off_block_residual": _decimal(result["left_native_off_block_residual"]),
        "right_native_lower_left_residual": _decimal(result["right_native_lower_left_residual"]),
        "left_native_lower_left_residual": _decimal(result["left_native_lower_left_residual"]),
        "attempts": attempts,
        "status": "TRUSTED",
    }


def _unresolved_row(epsilon_token, gamma_token, attempts):
    return {
        "epsilon_token": epsilon_token,
        "gamma_token": gamma_token,
        "attempts": attempts,
        "status": "UNRESOLVED",
    }


def _scaling_reads(rows):
    trusted = {
        (row["gamma_token"], Fraction(row["epsilon_token"])): row
        for row in rows if row["status"] == "TRUSTED"
    }
    reads = []
    for gamma in sorted({key[0] for key in trusted}, key=Fraction):
        positive = sorted(
            {epsilon for row_gamma, epsilon in trusted if row_gamma == gamma and epsilon > 0},
            reverse=True,
        )
        for epsilon in positive:
            half = epsilon / 2
            keys = ((gamma, epsilon), (gamma, half), (gamma, -epsilon), (gamma, -half))
            if not all(key in trusted for key in keys):
                continue
            plus, plus_half, minus, minus_half = (trusted[key] for key in keys)
            dp, dph, dm, dmh = (float(row["A_gap"]) for row in (plus, plus_half, minus, minus_half))
            ep, eph, em, emh = (
                float(row["high_diagnostics"]["eigenvalue_error"])
                for row in (plus, plus_half, minus, minus_half)
            )
            p_plus = math.log(dp / dph, 2)
            p_minus = math.log(dm / dmh, 2)
            p_plus_error = (ep / dp + eph / dph) / math.log(2)
            p_minus_error = (em / dm + emh / dmh) / math.log(2)
            e = float(epsilon)
            c2 = (dp + dm) / (2 * e * e)
            c3 = (dp - dm) / (2 * e**3)
            c2_error = (ep + em) / (2 * e * e)
            c3_error = (ep + em) / (2 * abs(e)**3)

            def reading(value, error):
                return {
                    "value": _decimal(value),
                    "interval": [_decimal(value - error), _decimal(value + error)],
                }

            reads.append({
                "gamma_token": gamma,
                "epsilon_token": _fraction_token(epsilon),
                "p_plus": reading(p_plus, p_plus_error),
                "p_minus": reading(p_minus, p_minus_error),
                "c2": reading(c2, c2_error),
                "c3": reading(c3, c3_error),
                "regression_only": "near 2, gamma/2, gamma/2; symbolic series is the theorem gate, not finite convergence",
            })
    return reads


def _branch_paths(epsilon_tokens):
    values = {Fraction(token): token for token in epsilon_tokens}
    paths = {}
    for sign, label in ((1, "positive"), (-1, "negative")):
        branch = [Fraction(0)]
        historical = Fraction(sign, 10)
        if historical in values:
            branch.append(historical)
        eighth = Fraction(sign, 8)
        if eighth in values:
            branch.append(eighth)
        branch.extend(sorted(
            (value for value in values if value * sign > 0 and abs(value) < Fraction(1, 8)),
            key=abs,
            reverse=True,
        ))
        paths[label] = [_fraction_token(value) for value in branch]
    return paths


def build_certificate(
    *,
    gamma_tokens=GAMMA_TOKENS,
    epsilon_tokens=EPSILON_TOKENS,
    measure=measure_precision_pair,
):
    paths = _branch_paths(epsilon_tokens)
    continuation = {}
    continuation_failures = []
    continuation_cache = {}

    def prepared(gamma_token, bits):
        key = (gamma_token, bits)
        if key not in continuation_cache:
            continuation_cache[key] = prepare_continuation_paths(
                paths, gamma_token, bits
            )
        return continuation_cache[key]

    for gamma_token in gamma_tokens:
        try:
            continuation[gamma_token] = {}
            for bits in (53, 106):
                raw, _ = prepared(gamma_token, bits)
                continuation[gamma_token][str(bits)] = {
                    branch: [
                        {
                            key: (_decimal(value) if isinstance(value, float) else value)
                            for key, value in step.items()
                            if key != "fresh_adjoint_build_token"
                        }
                        for step in steps
                    ]
                    for branch, steps in raw.items()
                }
        except Exception as error:
            message = f"{type(error).__name__}: {error}"
            continuation[gamma_token] = {"status": "UNRESOLVED", "failures": [message]}
            continuation_failures.append(message)
    rows = []
    for gamma_token in gamma_tokens:
        for epsilon_token in epsilon_tokens:
            attempts = []
            accepted = None
            for low_bits, high_bits in PRECISION_PAIRS:
                try:
                    if measure is measure_precision_pair:
                        _, low_states = prepared(gamma_token, low_bits)
                        _, high_states = prepared(gamma_token, high_bits)
                        result = measure(
                            epsilon_token,
                            gamma_token,
                            low_bits,
                            high_bits,
                            low_state=low_states[epsilon_token],
                            high_state=high_states[epsilon_token],
                        )
                    else:
                        result = measure(epsilon_token, gamma_token, low_bits, high_bits)
                    decision = decide_precision_pair(
                        result["low"],
                        result["high"],
                        independently_reconstructed=precision_matrices_are_independent(
                            result["low_build"], result["high_build"]
                        ),
                    )
                    attempts.append({
                        "precision_bits": [low_bits, high_bits],
                        "status": decision.status,
                        "failures": list(decision.failures),
                    })
                    if decision.status == "TRUSTED":
                        accepted = _row_from_measurement(
                            epsilon_token, gamma_token, low_bits, high_bits, result, attempts
                        )
                        break
                except Exception as error:
                    attempts.append({
                        "precision_bits": [low_bits, high_bits],
                        "status": "UNRESOLVED",
                        "failures": [f"{type(error).__name__}: {error}"],
                    })
            rows.append(accepted or _unresolved_row(epsilon_token, gamma_token, attempts))
    verdict = (
        "PASS"
        if not continuation_failures and all(row["status"] == "TRUSTED" for row in rows)
        else "UNRESOLVED"
    )
    operator_mutation = b_operator_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    spectrum_mutation = b_spectrum_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    document = {
        "schema": "missing-phase-relaxation-scale/v1",
        "scope": {"N": 7, "sector": [1, 1], "seat": 3, "path": "r=1+epsilon"},
        "theorem": {
            "gap_series": "gamma/2*epsilon^2 + gamma/2*epsilon^3 + O(epsilon^4)",
            "tau_leading": "2/(gamma*epsilon^2)",
            "observable_lifetime_claimed": False,
            "full_liouvillian_gap_claimed": False,
        },
        "exact": {
            "polynomial_match": True,
            "uniform_peripheral_dimension": 10,
            "uniform_kernel_dimension": 4,
            "punctured_kernel_dimension": 2,
            "b_boundary_has_peripheral_mode": False,
        },
        "rows": rows,
        "controls": {
            "continuation_paths": paths,
            "continuation_records": continuation,
            "continuation_failures": continuation_failures,
            "candidate_radius_rule": "one third of preceding trusted sep_complex; never widened",
            "B_relation_gate": "fresh full operators: B + 2*gamma*I + A_dagger",
            "full_B_spectrum": "numerical read; certification follows from the operator gate",
            "B_mutations": {
                "correct_operator_residual_zero": operator_mutation["correct"] == 0,
                "omitted_adjoint_entrywise_detected": operator_mutation["omitted_adjoint"] != 0,
                "correct_spectrum_residual_zero": spectrum_mutation["correct"] == 0,
                "missing_right_action_spectrum_detected": spectrum_mutation["missing_right_action"] != 0,
                "wrong_price_sign_spectrum_detected": spectrum_mutation["wrong_price_sign"] != 0,
            },
        },
        "verdict": verdict,
    }
    if verdict == "PASS":
        document["scaling"] = _scaling_reads(rows)
    return document


def write_atomic_json(output, document):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    encoded = (json.dumps(document, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("utf-8")
    with temporary.open("wb") as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)


def run_certificate(
    output=RESULT_PATH,
    *,
    gamma_tokens=GAMMA_TOKENS,
    epsilon_tokens=EPSILON_TOKENS,
    measure=measure_precision_pair,
):
    try:
        document = build_certificate(
            gamma_tokens=gamma_tokens,
            epsilon_tokens=epsilon_tokens,
            measure=measure,
        )
    except Exception as error:
        document = {
            "schema": "missing-phase-relaxation-scale/v1",
            "scope": {"N": 7, "sector": [1, 1], "seat": 3, "path": "r=1+epsilon"},
            "rows": [{"status": "UNRESOLVED", "attempts": [{
                "status": "UNRESOLVED", "failures": [f"{type(error).__name__}: {error}"]
            }]}],
            "controls": {},
            "verdict": "UNRESOLVED",
        }
    write_atomic_json(output, document)
    print(document["verdict"])
    return 0 if document["verdict"] == "PASS" else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULT_PATH)
    arguments = parser.parse_args(argv)
    return run_certificate(arguments.output)


def replace_diagnostic(value: PrecisionDiagnostics, **changes) -> PrecisionDiagnostics:
    data = asdict(value)
    data.update(changes)
    return PrecisionDiagnostics(**data)


def decide_precision_pair(
    low: PrecisionDiagnostics,
    high: PrecisionDiagnostics,
    *,
    independently_reconstructed: bool,
) -> PrecisionDecision:
    """Return the one authoritative list of failed precision inequalities."""
    failures = []
    if not independently_reconstructed:
        failures.append("precision matrices were not independently reconstructed")
    if low.rank != 2 or high.rank != 2:
        failures.append("cluster rank is not exactly two")
    numeric = (
        getattr(record, field)
        for record in (low, high)
        for field in PrecisionDiagnostics.__dataclass_fields__
        if field != "rank"
    )
    if not all(math.isfinite(float(value)) for value in numeric):
        failures.append("a diagnostic field is not finite")
    for label, record in (("low", low), ("high", high)):
        if record.sep_sylvester <= 0:
            failures.append(f"{label} Sylvester separation is not positive")
            continue
        derived = record.eigenvalue_error / record.sep_sylvester
        if not math.isclose(record.subspace_error, derived, rel_tol=1e-12, abs_tol=0.0):
            failures.append(f"{label} subspace error is not E_lambda/sep_Syl")
        if not 2 * record.eigenvalue_error < min(record.gap, record.sep_complex / 2):
            failures.append(f"{label} eigenvalue isolation inequality failed")
        if not record.max_principal_angle < math.pi / 2:
            failures.append(f"{label} principal angle reached pi/2")
    if low.sep_sylvester > 0 and high.sep_sylvester > 0:
        radius = low.eigenvalue_error / low.sep_sylvester + high.eigenvalue_error / high.sep_sylvester
        if not radius < 0.5:
            failures.append("combined subspace radius is not below one half")
        if max(low.projector_difference, high.projector_difference) > radius:
            failures.append("orthogonal projector drift exceeds the subspace radius")
    if not 2 * max(low.light_difference, high.light_difference) < min(low.light, high.light):
        failures.append("centre-light interval reaches zero")
    if max(low.eigenvalue_pair_error, high.eigenvalue_pair_error) > low.eigenvalue_error + high.eigenvalue_error:
        failures.append("p/2p eigenvalue pairing exceeds the error interval")
    return PrecisionDecision("TRUSTED" if not failures else "UNRESOLVED", tuple(failures))


def classify_precision_pair(
    low: PrecisionDiagnostics,
    high: PrecisionDiagnostics,
    *,
    independently_reconstructed: bool,
) -> str:
    return decide_precision_pair(
        low, high, independently_reconstructed=independently_reconstructed
    ).status


def spectral_norm_mp(matrix):
    values = mp.svd(matrix, compute_uv=False)
    return max(values[index] for index in range(len(values))) if len(values) else mp.mpf(0)


def _exact_scalar(value):
    value = sympy.sympify(value)
    if value.has(sympy.Float):
        raise TypeError("exact reduction requires exact scalar inputs")
    return value


def _positive_exact_numeric_gamma(value):
    value = _exact_scalar(value)
    if value.is_number is not True:
        raise TypeError("gamma must be an exact numeric scalar")
    if value.is_real is not True or value.is_positive is not True:
        raise ValueError("gamma must be strictly positive and real")
    return value


def _positive_real_indeterminate(value):
    value = _exact_scalar(value)
    if not isinstance(value, sympy.Symbol):
        raise TypeError("gamma must be a symbolic indeterminate")
    if value.is_real is not True or value.is_positive is not True:
        raise ValueError("gamma indeterminate must be declared positive and real")
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


def a_generator_float(
    epsilon: float,
    gamma: float,
    seat: int = SEAT,
    *,
    both_ends: bool = False,
) -> np.ndarray:
    """Assemble the 49-dimensional A generator for row-major stacking."""
    h = hopping_float(epsilon, both_ends=both_ends)
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


def certify_uniform_peripheral_bases(gamma, operator_bases):
    """Certify independence and exhaustive peripheral completeness."""
    gamma = _positive_exact_numeric_gamma(gamma)
    h, eigenspaces, outer_identity, frequencies, _ = _peripheral_bases()
    _, z, _, _ = exact_reduction_objects(sympy.Integer(0), gamma)

    per_frequency_ranks = {}
    stacked_by_frequency = {}
    for key, basis in operator_bases.items():
        frequency = frequencies[key]
        for operator in basis:
            if exact_a_action(h, z, gamma, operator) != frequency * operator:
                raise AssertionError("peripheral operator failed its exact frequency gate")
        stacked = sympy.Matrix.hstack(
            *[_row_stack_exact(operator) for operator in basis]
        )
        stacked_by_frequency[key] = stacked
        per_frequency_ranks[key] = _exact_rank(stacked)
        if per_frequency_ranks[key] != len(basis):
            raise AssertionError(f"{key} peripheral basis failed its rank gate")

    combined_operator_basis = sympy.Matrix.hstack(
        *(stacked_by_frequency[key] for key in frequencies)
    )
    combined_operator_rank = _exact_rank(combined_operator_basis)
    _, zero_cost_invariant_basis = _zero_cost_invariant_subspace(h, z)
    combined_span_rank = _exact_rank(
        sympy.Matrix.hstack(combined_operator_basis, zero_cost_invariant_basis)
    )
    if (
        combined_operator_rank != zero_cost_invariant_basis.cols
        or combined_span_rank != zero_cost_invariant_basis.cols
    ):
        raise AssertionError("peripheral operator span is not complete")

    generator = exact_a_generator(sympy.Integer(0), gamma)
    kernel_dimension = N * N - _exact_rank(generator)
    if _exact_rank(outer_identity) != 4:
        raise AssertionError("outer Hilbert complement must have dimension four")
    census = {
        "blind_dimension": len(eigenspaces),
        "end_blind_algebra_dimension": len(eigenspaces) ** 2,
        "identity_outer_dimension": 1,
        "kernel_dimension": kernel_dimension,
        "peripheral_dimension": combined_operator_rank,
        "frequency_multiplicities": per_frequency_ranks,
    }
    return {
        "census": census,
        "operator_bases": operator_bases,
        "per_frequency_ranks": per_frequency_ranks,
        "combined_operator_rank": combined_operator_rank,
        "zero_cost_invariant_dimension": zero_cost_invariant_basis.cols,
        "combined_span_rank": combined_span_rank,
    }


def exact_uniform_peripheral_certificate(gamma):
    """Return the exact rank-and-invariance certificate at uniform coupling."""
    gamma = _positive_exact_numeric_gamma(gamma)
    _, _, _, _, operator_bases = _peripheral_bases()
    return certify_uniform_peripheral_bases(gamma, operator_bases)


def exact_uniform_peripheral_census(gamma):
    """Return the certified uniform census for exact numeric gamma>0."""
    return exact_uniform_peripheral_certificate(gamma)["census"]


def _polynomial_inverse(matrix, gamma):
    """Return an exact numerator/denominator inverse over K[gamma]."""
    coefficient_field = sympy.QQ.algebraic_field(sympy.sqrt(2), sympy.I)
    polynomial_ring = coefficient_field.poly_ring(gamma)
    domain_matrix = DomainMatrix.from_Matrix(matrix).convert_to(polynomial_ring)
    numerator, denominator = domain_matrix.inv_den(method="rref")
    return numerator, denominator, polynomial_ring


def _positive_denominator_certificate(expression, gamma):
    """Prove non-vanishing on gamma>0 by an exact real/imaginary gcd."""
    real_part = sympy.expand(sympy.re(expression))
    imaginary_part = sympy.expand(sympy.im(expression))
    extension = [sympy.sqrt(2)]
    real_polynomial = (
        sympy.Poly(real_part, gamma, extension=extension)
        if real_part != 0
        else None
    )
    imaginary_polynomial = (
        sympy.Poly(imaginary_part, gamma, extension=extension)
        if imaginary_part != 0
        else None
    )
    if real_polynomial is None:
        common_real_zeros = imaginary_polynomial
    elif imaginary_polynomial is None:
        common_real_zeros = real_polynomial
    else:
        common_real_zeros = sympy.gcd(real_polynomial, imaginary_polynomial)
    zero_root_multiplicity = min(
        monomial[0] for monomial, _ in common_real_zeros.terms()
    )
    away_from_zero = common_real_zeros.exquo(
        sympy.Poly(gamma**zero_root_multiplicity, gamma, extension=extension)
    )
    positive_real_root_count = away_from_zero.count_roots(0, sympy.oo)
    return {
        "real_part": real_part,
        "imaginary_part": imaginary_part,
        "common_real_zero_polynomial": common_real_zeros,
        "zero_root_multiplicity": zero_root_multiplicity,
        "positive_axis_polynomial": away_from_zero,
        "positive_real_root_count": positive_real_root_count,
        "holds": positive_real_root_count == 0,
    }


def _resolvent_blocks(l0, frequency, projector, gamma):
    """Describe the full 49-dimensional reduced resolvent by invariant blocks."""
    shifted = l0 - frequency * sympy.eye(N * N) + projector
    blind = range(3)
    outer = range(3, 7)
    index_blocks = []
    index_blocks.extend([[left * N + right for right in outer] for left in blind])
    index_blocks.extend([[left * N + right for left in outer] for right in blind])
    index_blocks.append([left * N + right for left in blind for right in blind])
    index_blocks.append([left * N + right for left in outer for right in outer])
    blocks = []
    polynomial_ring = None
    for indices in index_blocks:
        local_shifted = shifted.extract(indices, indices)
        numerator, denominator, polynomial_ring = _polynomial_inverse(
            local_shifted, gamma
        )
        blocks.append(
            {
                "indices": tuple(indices),
                "shifted": local_shifted,
                "unshifted": (
                    l0 - frequency * sympy.eye(N * N)
                ).extract(indices, indices),
                "projector": projector.extract(indices, indices),
                "inverse_numerator_domain": numerator,
                "inverse_denominator_domain": denominator,
            }
        )

    common_denominator = blocks[0]["inverse_denominator_domain"]
    for block in blocks[1:]:
        common_denominator = polynomial_ring.lcm(
            common_denominator, block["inverse_denominator_domain"]
        )
    common_denominator_expression = polynomial_ring.to_sympy(common_denominator)
    global_numerator = sympy.zeros(N * N, N * N)
    assembled_shifted = sympy.zeros(N * N, N * N)
    for block in blocks:
        indices = list(block["indices"])
        denominator = block["inverse_denominator_domain"]
        scale = polynomial_ring.exquo(common_denominator, denominator)
        projector_domain = DomainMatrix.from_Matrix(block["projector"]).convert_to(
            polynomial_ring
        )
        reduced_numerator = (
            block["inverse_numerator_domain"] - projector_domain * denominator
        ) * scale
        reduced_numerator_matrix = reduced_numerator.to_Matrix()
        for local_row, global_row in enumerate(indices):
            for local_column, global_column in enumerate(indices):
                global_numerator[global_row, global_column] = (
                    reduced_numerator_matrix[local_row, local_column]
                )
                assembled_shifted[global_row, global_column] = block["shifted"][
                    local_row, local_column
                ]
        block["inverse_denominator"] = polynomial_ring.to_sympy(denominator)
        del block["inverse_numerator_domain"]
        del block["inverse_denominator_domain"]

    positive_certificate = _positive_denominator_certificate(
        common_denominator_expression, gamma
    )
    return {
        "blocks": blocks,
        "numerator": global_numerator,
        "denominator": common_denominator_expression,
        "shifted_off_block_residual": shifted - assembled_shifted,
        "positive_denominator_certificate": positive_certificate,
    }


def exact_effective_operators(gamma):
    """Derive Kato operators for a positive-real symbolic gamma."""
    gamma = _positive_real_indeterminate(gamma)
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
        resolvent = _resolvent_blocks(l0, frequency, projector, gamma)
        blocks = resolvent["blocks"]
        f1 = (sympy.conjugate(left.T) * l1 * right).applyfunc(sympy.simplify)
        f2_numerator = (
            -sympy.conjugate(left.T)
            * l1
            * resolvent["numerator"]
            * l1
            * right
        )
        f2 = f2_numerator.applyfunc(
            lambda entry: sympy.simplify(
                sympy.cancel(
                    entry / resolvent["denominator"],
                    extension=[sympy.sqrt(2), sympy.I],
                )
            )
        )
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
            "resolvent_numerator": resolvent["numerator"],
            "resolvent_denominator": resolvent["denominator"],
            "shifted_off_block_residual": resolvent[
                "shifted_off_block_residual"
            ],
            "denominator_nonzero_for_positive_gamma": resolvent[
                "positive_denominator_certificate"
            ]["holds"],
            "positive_denominator_certificate": resolvent[
                "positive_denominator_certificate"
            ],
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


def _zero_cost_invariant_subspace(h, z):
    """Return -i ad_h and its largest invariant subspace in ker D_z."""
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
    return commutator, _largest_invariant_subspace(commutator, forbidden_rows)


def exact_punctured_kernel_certificate(epsilon, gamma):
    """Certify rational punctured cases for exact numeric gamma>0."""
    epsilon = _exact_scalar(epsilon)
    gamma = _positive_exact_numeric_gamma(gamma)
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

    commutator, invariant_basis = _zero_cost_invariant_subspace(h, z)

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


def b_generator_float(epsilon: float, gamma: float, *, seat: int = SEAT) -> np.ndarray:
    """Build the B carrier directly, before using the priced-adjoint identity."""
    h = hopping_float(epsilon)
    z = watched_involution_float(seat)
    identity = np.eye(N, dtype=complex)
    return (
        -1j * (np.kron(h, identity) - np.kron(identity, h.T))
        - gamma * (np.kron(z, z.T) + np.eye(N * N, dtype=complex))
    )


def _exact_action(h_left, h_right, z_left, z_right, gamma, matrix):
    return (
        -sympy.I * (h_left * matrix - matrix * h_right)
        + gamma * (z_left * matrix * z_right - matrix)
    ).applyfunc(sympy.expand)


def _exact_generator_from_actions(h_left, h_right, z_left, z_right, gamma):
    columns = []
    for coordinate in range(N * N):
        basis = sympy.zeros(N, N)
        basis[coordinate // N, coordinate % N] = 1
        columns.append(
            _row_stack_exact(
                _exact_action(h_left, h_right, z_left, z_right, gamma, basis)
            )
        )
    return sympy.Matrix.hstack(*columns)


def _exact_b_generator_from_h_z(h, z, gamma):
    return _exact_generator_from_actions(h, h, z, -z, gamma)


def _zero_cost_invariant_subspace_for_dissipator(h, z, *, carrier):
    """Find the largest commutator-invariant subspace with zero dissipation."""
    commutator_columns = []
    dissipator_columns = []
    for coordinate in range(N * N):
        basis = sympy.zeros(N, N)
        basis[coordinate // N, coordinate % N] = 1
        commutator_columns.append(
            _row_stack_exact(-sympy.I * (h * basis - basis * h))
        )
        if carrier == "A":
            dissipator = z * basis * z - basis
        elif carrier == "B":
            dissipator = -(z * basis * z + basis)
        else:
            raise ValueError("carrier must be A or B")
        dissipator_columns.append(_row_stack_exact(dissipator))
    commutator = sympy.Matrix.hstack(*commutator_columns)
    dissipator = sympy.Matrix.hstack(*dissipator_columns)
    charged_coordinates = [
        index for index in range(N * N) if dissipator[index, index] != 0
    ]
    forbidden_rows = sympy.zeros(len(charged_coordinates), N * N)
    for row, coordinate in enumerate(charged_coordinates):
        forbidden_rows[row, coordinate] = 1
    return commutator, _largest_invariant_subspace(commutator, forbidden_rows)


def _centre_first_partition(h, seat=SEAT):
    order = [seat] + [index for index in range(N) if index != seat]
    permuted = h.extract(order, order)
    return permuted[1:, 1:], permuted[1:, 0]


def exact_b_boundary_certificate():
    """Certify that the physical centre/outer boundary has no B-peripheral mode."""
    h, z, _, _ = exact_reduction_objects(
        sympy.Integer(0), sympy.Rational(3, 10)
    )
    h_outer, u = _centre_first_partition(h)
    _, invariant = _zero_cost_invariant_subspace_for_dissipator(
        h, z, carrier="B"
    )
    product = h_outer * u
    if product != sympy.Matrix([0, 4, 0, 0, 4, 0]):
        raise AssertionError("centre/outer product changed")
    return {
        "defect_epsilon": sympy.Integer(0),
        "h_outer_times_u": product,
        "has_peripheral_b_mode": invariant.cols != 0,
        "peripheral_dimension": invariant.cols,
        "peripheral_basis": invariant,
    }


def exact_b_boundary_control():
    """Build a boundary with nonzero u in ker(h_outer) and verify the B case."""
    gamma = sympy.Rational(3, 10)
    u = sympy.Matrix([0, 0, 2, 2, 0, 0])
    h_outer = sympy.zeros(6, 6)
    h_outer[0, 1] = h_outer[1, 0] = 2
    h_outer[4, 5] = h_outer[5, 4] = 2
    h = sympy.zeros(N, N)
    h[0, 1:] = u.T
    h[1:, 0] = u
    h[1:, 1:] = h_outer
    z = sympy.diag(-1, 1, 1, 1, 1, 1, 1)
    generator = _exact_b_generator_from_h_z(h, z, gamma)
    _, invariant = _zero_cost_invariant_subspace_for_dissipator(
        h, z, carrier="B"
    )
    if invariant.cols != 1:
        raise AssertionError("control must create exactly one peripheral B direction")
    witness = sympy.zeros(N, N)
    witness[0, 1:] = u.T
    witness[1:, 0] = u
    residual = generator * _row_stack_exact(witness)
    if residual != sympy.zeros(N * N, 1):
        raise AssertionError("displayed control B direction is not stationary")
    return {
        "u_nonzero": u != sympy.zeros(6, 1),
        "h_outer_times_u": h_outer * u,
        "has_peripheral_b_mode": invariant.cols != 0,
        "peripheral_dimension": invariant.cols,
        "direct_generator_residual": residual,
    }


def _public_float_as_exact(value):
    """Recover intended short decimal controls without introducing SymPy Float."""
    if not isinstance(value, (float, int, np.floating, np.integer)):
        return _exact_scalar(value)
    return sympy.Rational(str(float(value))).limit_denominator(10**9)


def _exact_control_h_z(epsilon, *, seat, both_ends):
    epsilon = _exact_scalar(epsilon)
    r = 1 + epsilon
    h = sympy.zeros(N, N)
    bonds = (2 * r, 2, 2, 2, 2, 2 * r if both_ends else 2)
    for site, bond in enumerate(bonds):
        h[site, site + 1] = h[site + 1, site] = bond
    z = sympy.eye(N)
    z[seat, seat] = -1
    return h, z


def exact_control_subspaces(epsilon, gamma, *, seat=SEAT, both_ends=False):
    """Return exact stationary/peripheral bases for the three Task-4 controls."""
    epsilon = _exact_scalar(epsilon)
    gamma = _positive_exact_numeric_gamma(gamma)
    h, z = _exact_control_h_z(epsilon, seat=seat, both_ends=both_ends)
    generator = _exact_generator_from_h_z(h, z, gamma)
    commutator, peripheral = _zero_cost_invariant_subspace_for_dissipator(
        h, z, carrier="A"
    )
    watched = sympy.eye(N).col(seat)
    krylov = sympy.Matrix.hstack(*[(h**power) * watched for power in range(N)])
    blind_dimension = len(krylov.T.nullspace())
    gram = sympy.conjugate(peripheral.T) * peripheral
    restricted = (
        gram.inv()
        * sympy.conjugate(peripheral.T)
        * commutator
        * peripheral
    )
    stationary_coefficients = restricted.nullspace()
    stationary = sympy.Matrix.hstack(
        *[peripheral * vector for vector in stationary_coefficients]
    )

    if both_ends and seat == SEAT and 1 + epsilon in (
        sympy.Rational(9, 10),
        sympy.Rational(11, 10),
    ):
        expected = (4, 10)
    elif not both_ends and seat == SEAT and epsilon == 0:
        expected = (4, 10)
    elif not both_ends and seat == SEAT and epsilon != 0:
        expected = (2, 2)
    elif not both_ends and seat == 2 and epsilon == 0:
        expected = (1, 1)
    else:
        raise ValueError("unsupported exact Task-4 control")

    stationary_dimension = _exact_rank(stationary)
    peripheral_dimension = _exact_rank(peripheral)
    if (stationary_dimension, peripheral_dimension) != expected:
        raise AssertionError("control subspace dimensions changed")
    if peripheral_dimension != blind_dimension**2 + 1:
        raise AssertionError("peripheral space is not the complete blind algebra plus identity")
    if generator * stationary != sympy.zeros(N * N, stationary.cols):
        raise AssertionError("stationary basis failed the exact generator gate")
    complement = (
        sympy.eye(N * N)
        - peripheral * gram.inv() * sympy.conjugate(peripheral.T)
    )
    if complement * commutator * peripheral != sympy.zeros(
        N * N, peripheral.cols
    ):
        raise AssertionError("peripheral span failed the exact invariance gate")
    dissipator_residual = generator * peripheral - commutator * peripheral
    if dissipator_residual != sympy.zeros(N * N, peripheral.cols):
        raise AssertionError("peripheral span is not zero-cost")
    return {
        "stationary_basis": stationary,
        "peripheral_basis": peripheral,
        "stationary_dimension": stationary_dimension,
        "peripheral_dimension": peripheral_dimension,
        "blind_dimension": blind_dimension,
        "restricted_commutator": restricted,
    }


def _float_generators(epsilon, gamma, *, seat, both_ends):
    a = a_generator_float(
        epsilon, gamma, seat=seat, both_ends=both_ends
    )
    h = hopping_float(epsilon, both_ends=both_ends)
    z = watched_involution_float(seat)
    identity = np.eye(N, dtype=complex)
    b = (
        -1j * (np.kron(h, identity) - np.kron(identity, h.T))
        - gamma * (np.kron(z, z.T) + np.eye(N * N, dtype=complex))
    )
    return a, b


def classify_a_rates(
    eigenvalues,
    tolerance,
    *,
    expected_peripheral,
    expected_stationary,
):
    """Partition every A eigenvalue into peripheral or strictly decaying."""
    values = np.asarray(eigenvalues, dtype=complex)
    peripheral_mask = np.abs(values.real) <= tolerance
    decaying_mask = -values.real > tolerance
    unstable_mask = values.real > tolerance
    if np.any(unstable_mask):
        raise AssertionError("unstable A eigenvalue escaped the physical half-plane")
    classified = int(np.count_nonzero(peripheral_mask)) + int(
        np.count_nonzero(decaying_mask)
    )
    if classified != values.size:
        raise AssertionError("not all A eigenvalues were classified")
    if int(np.count_nonzero(peripheral_mask)) != expected_peripheral:
        raise AssertionError("float peripheral count disagrees with exact certificate")
    stationary_mask = np.abs(values) <= tolerance
    if int(np.count_nonzero(stationary_mask)) != expected_stationary:
        raise AssertionError("float stationary count disagrees with exact certificate")
    return {
        "peripheral_mask": peripheral_mask,
        "decaying_mask": decaying_mask,
        "stationary_mask": stationary_mask,
        "classified_dimension": classified,
    }


def match_b_spectrum(a_values, b_values, gamma, matrix_scale):
    """Match the full complex B spectrum to the priced A-adjoint spectrum."""
    a_values = np.asarray(a_values, dtype=complex)
    b_values = np.asarray(b_values, dtype=complex)
    if a_values.size != N * N or b_values.size != N * N:
        raise AssertionError("A/B spectrum match requires all 49 eigenvalues")
    expected_b = -2 * gamma - np.conj(a_values)
    cost = np.abs(b_values[:, None] - expected_b[None, :])
    b_indices, a_indices = linear_sum_assignment(cost)
    residuals = cost[b_indices, a_indices]
    # Dense nonsymmetric eigensolves carry a dimension-amplified backward error.
    # 1024 eps ||L|| is safely above the measured 27 eps ||L|| control residual.
    tolerance = 1024 * np.finfo(float).eps * max(1.0, matrix_scale)
    maximum = float(np.max(residuals))
    if maximum > tolerance:
        raise AssertionError("full complex A/B spectrum map exceeds backward error")
    rate_residuals = [
        abs(
            -b_values[b_index].real
            - (2 * gamma - (-a_values[a_index].real))
        )
        for b_index, a_index in zip(b_indices, a_indices)
    ]
    return {
        "complex_residuals": residuals.tolist(),
        "max_complex_residual": maximum,
        "complex_tolerance": float(tolerance),
        "rate_residuals": rate_residuals,
    }


def direct_float_gaps(epsilon, gamma, *, seat=SEAT, both_ends=False):
    """Read all 49 A/B eigenvalues after exact case-specific classification."""
    exact_epsilon = _public_float_as_exact(epsilon)
    exact_gamma = _public_float_as_exact(gamma)
    certificate = exact_control_subspaces(
        exact_epsilon,
        exact_gamma,
        seat=seat,
        both_ends=both_ends,
    )
    a, b = _float_generators(
        epsilon, gamma, seat=seat, both_ends=both_ends
    )
    a_values = np.linalg.eigvals(a)
    b_values = np.linalg.eigvals(b)
    scale = max(1.0, np.linalg.norm(a, ord=2), np.linalg.norm(b, ord=2))
    tolerance = 128 * np.finfo(float).eps * scale
    classification = classify_a_rates(
        a_values,
        tolerance,
        expected_peripheral=certificate["peripheral_dimension"],
        expected_stationary=certificate["stationary_dimension"],
    )

    positive_a_rates = sorted(
        -value.real
        for value, decaying in zip(a_values, classification["decaying_mask"])
        if decaying
    )
    if not positive_a_rates:
        raise AssertionError("control has no decaying A eigenvalue")
    max_a_rate = max(-a_values.real)
    nonstationary_peripheral = (
        certificate["peripheral_dimension"] - certificate["stationary_dimension"]
    )
    a_gap = 0.0 if nonstationary_peripheral else positive_a_rates[0]
    next_distinct_a_rate = next(
        (
            rate
            for rate in positive_a_rates
            if rate > a_gap + tolerance
        ),
        positive_a_rates[0],
    )

    spectrum_match = match_b_spectrum(a_values, b_values, gamma, scale)
    b_gap = min(-b_values.real)
    return {
        "stationary_dimension": certificate["stationary_dimension"],
        "peripheral_dimension": certificate["peripheral_dimension"],
        "nonstationary_peripheral_dimension": nonstationary_peripheral,
        "a_gap": float(a_gap),
        "next_distinct_a_rate": float(next_distinct_a_rate),
        "max_a_rate": float(max_a_rate),
        "b_gap": float(b_gap),
        "paired_rate_residuals": spectrum_match["rate_residuals"],
        "complex_spectrum_residuals": spectrum_match["complex_residuals"],
        "max_complex_spectrum_residual": spectrum_match["max_complex_residual"],
        "complex_spectrum_tolerance": spectrum_match["complex_tolerance"],
        "classified_a_dimension": classification["classified_dimension"],
    }


def float_peripheral_count(epsilon, gamma, *, seat=SEAT, both_ends=False):
    """Return the numerical imaginary-axis count after exact classification."""
    return direct_float_gaps(
        epsilon, gamma, seat=seat, both_ends=both_ends
    )["peripheral_dimension"]


def tracked_reduced_gap(epsilon, gamma):
    """Coarsely continue the reduced root nearest -2 sqrt(2) i."""
    values = np.linalg.eigvals(k_reduced_float(epsilon, gamma))
    target = -2j * np.sqrt(2.0)
    root = values[np.argmin(np.abs(values - target))]
    return float(-root.real)


def direct_carrier_blocks(epsilon, gamma):
    """Construct all four global-flip carrier blocks from their own actions."""
    epsilon = _exact_scalar(epsilon)
    gamma = _positive_exact_numeric_gamma(gamma)
    h_1, z_1 = _exact_control_h_z(epsilon, seat=SEAT, both_ends=False)
    h_6, z_6_positive = _exact_control_h_z(
        epsilon, seat=SEAT, both_ends=False
    )
    z_6 = -z_6_positive
    return {
        "A_11": _exact_generator_from_actions(h_1, h_1, z_1, z_1, gamma),
        "A_66": _exact_generator_from_actions(h_6, h_6, z_6, z_6, gamma),
        "B_16": _exact_generator_from_actions(h_1, h_6, z_1, z_6, gamma),
        "B_61": _exact_generator_from_actions(h_6, h_1, z_6, z_1, gamma),
    }


def _coefficient_residual(left, right):
    differences = [
        sympy.expand(a - b)
        for a, b in zip(left.charpoly().all_coeffs(), right.charpoly().all_coeffs())
    ]
    return sympy.simplify(
        sum(sympy.conjugate(value) * value for value in differences)
    )


def b_spectrum_mutation_residuals(epsilon, gamma):
    """Gate the full monic B spectrum against mutations that change its spectrum."""
    epsilon = _exact_scalar(epsilon)
    gamma = _positive_exact_numeric_gamma(gamma)
    blocks = direct_carrier_blocks(epsilon, gamma)
    a = blocks["A_11"]
    b = blocks["B_16"]
    identity = sympy.eye(N * N)
    correct = -2 * gamma * identity - a.conjugate().T
    h, z = _exact_control_h_z(epsilon, seat=SEAT, both_ends=False)
    missing_right_action = _exact_generator_from_actions(
        h, sympy.zeros(N, N), z, -z, gamma
    )
    wrong_price_sign = 2 * gamma * identity - a.conjugate().T
    return {
        "correct": _coefficient_residual(b, correct),
        "missing_right_action": _coefficient_residual(
            b, missing_right_action
        ),
        "wrong_price_sign": _coefficient_residual(b, wrong_price_sign),
    }


def b_operator_mutation_residuals(epsilon, gamma):
    """Detect the adjoint entrywise; its omission is isospectral in this case."""
    epsilon = _exact_scalar(epsilon)
    gamma = _positive_exact_numeric_gamma(gamma)
    blocks = direct_carrier_blocks(epsilon, gamma)
    a = blocks["A_11"]
    b = blocks["B_16"]
    identity = sympy.eye(N * N)
    correct = -2 * gamma * identity - a.conjugate().T
    omitted_adjoint = -2 * gamma * identity - a

    def entrywise_residual(candidate):
        delta = b - candidate
        return sympy.simplify(
            sum(sympy.conjugate(value) * value for value in delta)
        )

    return {
        "correct": entrywise_residual(correct),
        "omitted_adjoint": entrywise_residual(omitted_adjoint),
    }


if __name__ == "__main__":
    raise SystemExit(main())
