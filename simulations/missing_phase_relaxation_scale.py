"""Invariant reduction for the N=7 centre-watched missing-phase A block."""

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math
import os
import platform
from pathlib import Path
from fractions import Fraction

import mpmath as mp
import numpy as np
import scipy
import sympy
from scipy.optimize import linear_sum_assignment
from sympy.polys.matrices import DomainMatrix


N = 7
SEAT = 3
RESULT_PATH = Path(__file__).resolve().parent / "results" / "missing_phase_relaxation_scale.json"

REQUIRED_B_MUTATION_GATES = frozenset({
    "correct_operator_residual_zero",
    "omitted_adjoint_entrywise_detected",
    "correct_spectrum_residual_zero",
    "missing_right_action_spectrum_detected",
    "wrong_price_sign_spectrum_detected",
})
REQUIRED_EXACT_GATES = frozenset({
    "polynomial_match",
    "gap_series_match",
    "effective_operator_match",
    "uniform_peripheral_dimension_match",
    "uniform_kernel_dimension_match",
    "punctured_kernel_dimension_match",
    "b_boundary_match",
})


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
class _BuildRecord:
    build_id: int
    matrix_id: int
    bits: int
    exact_input_fingerprint: tuple
    source_build_token: int | None


@dataclass(frozen=True)
class _MeasurementRecord:
    result_id: int
    fingerprint: str


@dataclass(frozen=True)
class ClusterState:
    epsilon_token: str
    rank: int
    right_basis: mp.matrix
    left_basis: mp.matrix
    right_ritz: tuple
    left_ritz: tuple


_build_counter = itertools.count(1)
_build_registry = {}
_measurement_counter = itertools.count(1)
_measurement_registry = {}


def _measurement_fingerprint_value(value):
    """Freeze a measured attempt without trusting caller-visible metadata."""
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return ("float", repr(value))
    if isinstance(value, Fraction):
        return ("fraction", value.numerator, value.denominator)
    if isinstance(value, mp.mpf):
        return ("mpf", value._mpf_)
    if isinstance(value, mp.mpc):
        return ("mpc", value._mpc_)
    if isinstance(value, mp.matrix):
        return (
            "matrix",
            value.rows,
            value.cols,
            tuple(
                _measurement_fingerprint_value(value[row, column])
                for row in range(value.rows)
                for column in range(value.cols)
            ),
        )
    if isinstance(value, PrecisionMatrix):
        return (
            "PrecisionMatrix",
            id(value),
            id(value.matrix),
            value.bits,
            value.token,
            value.independently_reconstructed,
            _measurement_fingerprint_value(value.exact_input_fingerprint),
            value.source_build_token,
            _measurement_fingerprint_value(value.matrix),
        )
    if isinstance(value, ClusterState):
        return (
            "ClusterState",
            value.epsilon_token,
            value.rank,
            id(value.right_basis),
            _measurement_fingerprint_value(value.right_basis),
            id(value.left_basis),
            _measurement_fingerprint_value(value.left_basis),
            _measurement_fingerprint_value(value.right_ritz),
            _measurement_fingerprint_value(value.left_ritz),
        )
    if isinstance(value, PrecisionDiagnostics):
        return ("PrecisionDiagnostics",) + tuple(
            _measurement_fingerprint_value(field_value)
            for field_value in asdict(value).values()
        )
    if isinstance(value, PrecisionDecision):
        return (
            "PrecisionDecision",
            value.status,
            _measurement_fingerprint_value(value.failures),
        )
    if isinstance(value, dict):
        return tuple(
            (key, _measurement_fingerprint_value(item))
            for key, item in sorted(value.items())
            if key != "_measurement_token"
        )
    if isinstance(value, (list, tuple)):
        return tuple(_measurement_fingerprint_value(item) for item in value)
    raise TypeError(f"unsupported measurement evidence type: {type(value).__name__}")


def _measurement_fingerprint(result):
    frozen = _measurement_fingerprint_value(result)
    return hashlib.sha256(repr(frozen).encode("utf-8")).hexdigest()


def _register_measurement_result(result):
    token = next(_measurement_counter)
    result["_measurement_token"] = token
    _measurement_registry[token] = _MeasurementRecord(
        result_id=id(result),
        fingerprint=_measurement_fingerprint(result),
    )


def _seal_measurement_result(result):
    token = result.get("_measurement_token")
    record = _measurement_registry.get(token)
    if record is None or record.result_id != id(result):
        raise AssertionError("measurement result cannot be sealed without provenance")
    _measurement_registry[token] = _MeasurementRecord(
        result_id=id(result),
        fingerprint=_measurement_fingerprint(result),
    )


def _measurement_result_is_authentic(result):
    if not isinstance(result, dict):
        return False
    record = _measurement_registry.get(result.get("_measurement_token"))
    if record is None or record.result_id != id(result):
        return False
    try:
        return record.fingerprint == _measurement_fingerprint(result)
    except (AttributeError, TypeError, ValueError):
        return False


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
    matrix_copy = matrixrezz.copy()
    token = next(_build_counter)
    result = PrecisionMatrix(
        matrix_copy,
        bits,
        token,
        source_build_token is None,
        tuple(exact_input_fingerprint),
        source_build_token,
    )
    _build_registry[token] = _BuildRecord(
        id(result),
        id(matrix_copy),
        bits,
        tuple(exact_input_fingerprint),
        source_build_token,
    )
    return result


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
    low_record = _build_registry.get(low.token)
    high_record = _build_registry.get(high.token)
    if low_record is None or high_record is None:
        return False
    return (
        low_record.build_id == id(low)
        and high_record.build_id == id(high)
        and low_record.matrix_id == id(low.matrix)
        and high_record.matrix_id == id(high.matrix)
        and low_record.bits == low.bits
        and high_record.bits == high.bits
        and low_record.bits * 2 == high_record.bits
        and low.token != high.token
        and low.matrix is not high.matrix
        and low_record.exact_input_fingerprint == low.exact_input_fingerprint
        and high_record.exact_input_fingerprint == high.exact_input_fingerprint
        and low_record.exact_input_fingerprint == high_record.exact_input_fingerprint
        and low_record.source_build_token is None
        and high_record.source_build_token is None
    )


def classify_precision_build_pair(low_diagnostic, high_diagnostic, low_build, high_build):
    return classify_precision_pair(
        low_diagnostic,
        high_diagnostic,
        independently_reconstructed=precision_matrices_are_independent(low_build, high_build),
    )


def decide_certificate_attempt(result):
    if not _measurement_result_is_authentic(result):
        return PrecisionDecision(
            "UNRESOLVED", ("measurement result provenance is invalid",)
        )
    precision = decide_precision_pair(
        result["low"],
        result["high"],
        independently_reconstructed=precision_matrices_are_independent(
            result["low_build"], result["high_build"]
        ),
    )
    failures = list(precision.failures)
    a_build_error = precision_build_discrepancy(
        result["low_build"], result["high_build"]
    )
    a_component_fields = (
        "right_cluster_residual",
        "left_cluster_residual",
        "biorthogonal_sigma_min",
        "condition_number",
        "construction_discrepancy",
        "rounding_term",
        "eigenvalue_error",
    )
    for label in ("low", "high"):
        state = result.get(f"{label}_state")
        stored = result.get(f"{label}_diagnostic_components", {})
        if state is None:
            failures.append(f"{label} A diagnostic has no bound cluster state")
            continue
        actual = _a_residual_components(
            result[f"{label}_build"], state, a_build_error
        )
        if (
            any(
                field not in stored
                or not math.isfinite(float(stored[field]))
                or not math.isclose(
                    float(stored[field]),
                    float(actual[field]),
                    rel_tol=1e-12,
                    abs_tol=0.0,
                )
                for field in a_component_fields
            )
            or not math.isclose(
                result[label].eigenvalue_error,
                actual["eigenvalue_error"],
                rel_tol=1e-12,
                abs_tol=0.0,
            )
        ):
            failures.append(f"{label} A diagnostic components are not bound to its build")
    k_builds = result.get("k_builds", ())
    k_builds_valid = (
        len(k_builds) == 2 and precision_matrices_are_independent(*k_builds)
    )
    if not k_builds_valid:
        failures.append("fresh K p/2p construction provenance is invalid")
    k_pair_error = result.get("k_eigenvalue_pair_error", float("nan"))
    k_pair_bound = result.get("k_eigenvalue_pair_bound", float("nan"))
    k_build_discrepancy = result.get("k_build_discrepancy", float("nan"))
    recomputed_k = {}
    if k_builds_valid:
        actual_k_discrepancy = precision_build_discrepancy(*k_builds)
        for label, build in zip(("low", "high"), k_builds):
            with mp.workprec(build.bits):
                actual_value, actual_components = _k_eigenvalue_diagnostic(
                    build,
                    mp.mpc(0, 2 * mp.sqrt(2)),
                    actual_k_discrepancy,
                )
            recomputed_k[label] = (actual_value, actual_components)
        if not math.isclose(
            float(k_build_discrepancy),
            float(actual_k_discrepancy),
            rel_tol=1e-12,
            abs_tol=0.0,
        ):
            failures.append("K p/2p build discrepancy is not bound to registered builds")
    if (
        not math.isfinite(float(k_pair_error))
        or not math.isfinite(float(k_pair_bound))
        or not math.isfinite(float(k_build_discrepancy))
        or k_pair_bound < 0
        or k_build_discrepancy < 0
        or k_build_discrepancy > k_pair_bound
        or k_pair_error > k_pair_bound
    ):
        failures.append("K p/2p discrepancy or eigenvalue pair exceeds its certified bound")
    for label in ("low", "high"):
        components = result.get(f"{label}_k_diagnostic_components", {})
        agreement = result.get(f"{label}_a_k_eigenvalue_agreement", float("nan"))
        agreement_bound = result.get(
            f"{label}_a_k_eigenvalue_agreement_bound", float("nan")
        )
        required_components = (
            "right_residual",
            "left_residual",
            "left_right_pair_residual",
            "biorthogonal_sigma_min",
            "condition_number",
            "construction_discrepancy",
            "rounding_term",
            "eigenvalue_error",
            "gap",
            "sep_complex",
        )
        finite_components = all(
            name in components and math.isfinite(float(components[name]))
            for name in required_components
        )
        derived_k_error = (
            components.get("condition_number", float("nan"))
            * (
                max(
                    components.get("right_residual", float("nan")),
                    components.get("left_residual", float("nan")),
                    components.get("left_right_pair_residual", float("nan")),
                )
                + components.get("construction_discrepancy", float("nan"))
                + components.get("rounding_term", float("nan"))
            )
        )
        expected_agreement_bound = (
            getattr(result[label], "eigenvalue_error", float("nan"))
            + components.get("eigenvalue_error", float("nan"))
        )
        actual_value, actual_components = recomputed_k.get(
            label, (mp.mpc("nan"), {})
        )
        stored_value = result.get(f"{label}_k_eigenvalue", mp.mpc("nan"))
        components_match_build = all(
            name in actual_components
            and name in components
            and math.isclose(
                float(components[name]),
                float(actual_components[name]),
                rel_tol=1e-12,
                abs_tol=0.0,
            )
            for name in required_components
        )
        state = result.get(f"{label}_state")
        if state is not None and k_builds_valid:
            build = k_builds[0 if label == "low" else 1]
            with mp.workprec(build.bits + 32):
                actual_agreement = max(
                    abs(value - actual_value) for value in state.right_ritz
                )
        else:
            actual_agreement = mp.inf
        actual_agreement_bound = (
            getattr(result[label], "eigenvalue_error", float("nan"))
            + actual_components.get("eigenvalue_error", float("nan"))
        )
        k_checks = (
            (finite_components, "K diagnostic components are finite"),
            (components_match_build, "K diagnostic components match the registered build"),
            (
                abs(mp.mpc(stored_value) - actual_value)
                <= k_pair_bound,
                "K Ritz value matches the registered build",
            ),
            (components.get("biorthogonal_sigma_min", 0) > 0, "K sigma_min is positive"),
            (components.get("condition_number", 0) >= 1, "K condition number is at least one"),
            (components.get("eigenvalue_error", -1) >= 0, "K eigenvalue error is nonnegative"),
            (
                math.isclose(
                    components.get("construction_discrepancy", float("nan")),
                    k_build_discrepancy,
                    rel_tol=1e-12,
                    abs_tol=0.0,
                ),
                "K construction discrepancy matches p/2p builds",
            ),
            (
                math.isclose(
                    components.get("eigenvalue_error", float("nan")),
                    derived_k_error,
                    rel_tol=1e-12,
                    abs_tol=0.0,
                ),
                "K eigenvalue error is derived from K residuals",
            ),
            (
                2 * components.get("eigenvalue_error", float("inf")) < min(
                    components.get("gap", float("nan")),
                    components.get("sep_complex", float("nan")) / 2,
                ),
                "K isolation inequality holds",
            ),
            (math.isfinite(float(agreement)), "A/K agreement is finite"),
            (
                math.isfinite(float(agreement_bound)) and agreement_bound >= 0,
                "A/K agreement bound is finite and nonnegative",
            ),
            (
                math.isclose(
                    agreement_bound,
                    expected_agreement_bound,
                    rel_tol=1e-12,
                    abs_tol=0.0,
                ),
                "A/K agreement bound is derived from stored diagnostics",
            ),
            (
                abs(float(agreement) - float(actual_agreement))
                <= float(actual_agreement_bound + k_pair_bound),
                "A/K agreement matches registered build and state",
            ),
            (
                actual_agreement <= actual_agreement_bound + k_pair_bound,
                "registered A/K Ritz values agree within p/2p bounds",
            ),
            (
                math.isclose(
                    float(agreement_bound),
                    float(actual_agreement_bound),
                    rel_tol=1e-12,
                    abs_tol=0.0,
                ),
                "A/K agreement bound matches registered diagnostics",
            ),
            (agreement <= agreement_bound, "A/K agreement lies within its bound"),
        )
        failures.extend(
            f"{label} {description} failed"
            for passed, description in k_checks if not passed
        )
    expected_k_pair_bound = sum(
        recomputed_k.get(label, (None, {}))[1].get(
            "eigenvalue_error", float("nan")
        ) for label in ("low", "high")
    )
    expected_k_pair_error = (
        abs(recomputed_k["low"][0] - recomputed_k["high"][0])
        if len(recomputed_k) == 2 else mp.nan
    )
    if not math.isclose(
        k_pair_bound, expected_k_pair_bound, rel_tol=1e-12, abs_tol=0.0
    ) or not math.isclose(
        float(k_pair_error), float(expected_k_pair_error), rel_tol=1e-12, abs_tol=0.0
    ):
        failures.append("K p/2p eigenvalue pair bound is not derived from K diagnostics")
    b_builds = result.get("b_builds", ())
    a_fingerprint = result["low_build"].exact_input_fingerprint
    expected_b_fingerprint = ("B", *a_fingerprint[1:])
    b_builds_valid = (
        len(b_builds) == 2
        and precision_matrices_are_independent(*b_builds)
        and all(build.matrix.rows == build.matrix.cols == 49 for build in b_builds)
        and all(
            build.exact_input_fingerprint == expected_b_fingerprint
            for build in b_builds
        )
    )
    if not b_builds_valid:
        failures.append("fresh B p/2p construction provenance is invalid")
    flip_pairs = result.get("global_flip_build_pairs", ())
    epsilon_token, gamma_token = a_fingerprint[1:]
    expected_flip_fingerprints = (
        ("A", "original", epsilon_token, gamma_token, 1, 1),
        ("A", "global_flip", epsilon_token, gamma_token, -1, -1),
        ("B", "original", epsilon_token, gamma_token, 1, -1),
        ("B", "global_flip", epsilon_token, gamma_token, -1, 1),
    )
    flip_pairs_valid = not (
        len(flip_pairs) != 4
        or any(
            len(pair) != 2 or not precision_matrices_are_independent(*pair)
            for pair in flip_pairs
        )
        or any(
            pair[0].exact_input_fingerprint != fingerprint
            or pair[1].exact_input_fingerprint != fingerprint
            or pair[0].matrix.rows != 49
            or pair[0].matrix.cols != 49
            or pair[1].matrix.rows != 49
            or pair[1].matrix.cols != 49
            for pair, fingerprint in zip(flip_pairs, expected_flip_fingerprints)
        )
    )
    if not flip_pairs_valid:
        failures.append("global-flip carrier construction provenance is invalid")
    absorption_error = result["absorption_error"]
    if not math.isfinite(float(absorption_error)) or absorption_error < 0:
        failures.append("absorption interval is not finite and nonnegative")
    for side, residual in result["absorption_residuals"].items():
        if not math.isfinite(float(residual)) or abs(residual) > absorption_error:
            failures.append(f"{side} absorption identity exceeds its interval")
    if b_builds_valid:
        b_build_error = precision_build_discrepancy(*b_builds)
        with mp.workprec(b_builds[1].bits):
            gamma = _mp_rational(gamma_token)
            operator_delta = (
                b_builds[1].matrix
                + 2 * gamma * mp.eye(49)
                + result["high_build"].matrix.transpose_conj()
            )
            actual_b_residual = (
                mp.mpf(0) if not any(operator_delta)
                else spectral_norm_mp(operator_delta)
            )
            actual_b_bound = a_build_error + b_build_error
        if (
            not math.isclose(
                float(result["b_operator_map_residual"]),
                float(actual_b_residual),
                rel_tol=1e-12,
                abs_tol=0.0,
            )
            or not math.isclose(
                float(result["b_operator_map_bound"]),
                float(actual_b_bound),
                rel_tol=1e-12,
                abs_tol=0.0,
            )
        ):
            failures.append("fresh B operator evidence is not bound to registered builds")
    if (
        not math.isfinite(float(result["b_operator_map_residual"]))
        or not math.isfinite(float(result["b_operator_map_bound"]))
        or result["b_operator_map_bound"] < 0
        or result["b_operator_map_residual"] > result["b_operator_map_bound"]
    ):
        failures.append("fresh B operator identity exceeds the construction bound")
    if flip_pairs_valid and b_builds_valid:
        (_, a_original), (_, a_copy), (_, b_original), (_, b_copy) = flip_pairs
        with mp.workprec(result["high_build"].bits):
            a_delta = a_original.matrix - a_copy.matrix
            b_delta = b_original.matrix - b_copy.matrix
            actual_copy_residual = max(
                mp.mpf(0) if not any(a_delta) else spectral_norm_mp(a_delta),
                mp.mpf(0) if not any(b_delta) else spectral_norm_mp(b_delta),
            )
            pair_drifts = [
                mp.norm(mp.matrix(low.matrix) - high.matrix)
                for low, high in flip_pairs
            ]
            actual_copy_bound = max(
                pair_drifts[0] + pair_drifts[1],
                pair_drifts[2] + pair_drifts[3],
            )
        if (
            not math.isclose(
                float(result["global_flip_copy_residual"]),
                float(actual_copy_residual),
                rel_tol=1e-12,
                abs_tol=0.0,
            )
            or not math.isclose(
                float(result["global_flip_copy_bound"]),
                float(actual_copy_bound),
                rel_tol=1e-12,
                abs_tol=0.0,
            )
        ):
            failures.append("global-flip evidence is not bound to registered carriers")
    if (
        not math.isfinite(float(result["global_flip_copy_residual"]))
        or not math.isfinite(float(result["global_flip_copy_bound"]))
        or result["global_flip_copy_bound"] < 0
        or result["global_flip_copy_residual"] > result["global_flip_copy_bound"]
    ):
        failures.append("global-flip carrier copies exceed the construction bound")
    decomposition_gates = result.get("decomposition_gates", ())
    expected_bits = {result["low_build"].bits, result["high_build"].bits}
    if (
        len(decomposition_gates) != 2
        or {gate.get("precision_bits") for gate in decomposition_gates} != expected_bits
    ):
        failures.append("decomposition gates do not cover both p and 2p builds")
    for gate in decomposition_gates:
        for field in (
            "unitarity_residual",
            "off_block_residual",
            "full_operator_reconstruction_residual",
            "full_eigen_residual",
        ):
            if (
                not math.isfinite(float(gate[field]))
                or not math.isfinite(float(gate["bound"]))
                or gate["bound"] < 0
                or gate[field] > gate["bound"]
            ):
                failures.append(
                    f"{gate['precision_bits']}-bit decomposition {field} exceeds its bound"
                )
    kernel_bounds = {
        gate.get("precision_bits"): gate.get("bound", float("nan"))
        for gate in decomposition_gates
    }
    for label, bits in (
        ("low", result["low_build"].bits),
        ("high", result["high_build"].bits),
    ):
        kernel_bound = kernel_bounds.get(bits, float("nan"))
        residual = result.get(f"{label}_kernel_residual", float("nan"))
        lower_left = result.get(f"{label}_kernel_lower_left_residual", float("nan"))
        if (
            result.get(f"{label}_exact_stationary_dimension") != 2
            or not math.isfinite(float(kernel_bound))
            or not math.isfinite(float(residual))
            or not math.isfinite(float(lower_left))
            or residual > kernel_bound
            or lower_left > kernel_bound
        ):
            failures.append(f"{label} exact stationary-space projection gate failed")
    next_rate = result.get("next_distinct_a_rate", float("nan"))
    rate_tolerance = result.get("next_distinct_a_rate_tolerance", float("nan"))
    if (
        not math.isfinite(float(next_rate))
        or not math.isfinite(float(rate_tolerance))
        or rate_tolerance <= 0
        or next_rate <= result["high"].gap + rate_tolerance
    ):
        failures.append("next distinct A rate is not separated from the tracked rate class")
    continuation = result.get("continuation_records", ())
    if len(continuation) != 2 or any(
        record is None or record.get("status") != "TRACKED"
        for record in continuation
    ):
        failures.append("p/2p continuation record is not TRACKED")
    return PrecisionDecision(
        "TRUSTED" if not failures else "UNRESOLVED", tuple(failures)
    )


def certificate_verdict(rows, continuation_failures, b_mutations, exact_gates):
    failures = list(continuation_failures)
    if not rows:
        failures.append("certificate contains no measured rows")
    elif any(row.get("status") != "TRUSTED" for row in rows):
        failures.append("one or more certificate rows are UNRESOLVED")
    elif any(not row.get("measurement_registry_bound", False) for row in rows):
        failures.append("one or more certificate rows lack authoritative measurement provenance")
    missing_b = REQUIRED_B_MUTATION_GATES.difference(b_mutations)
    if missing_b:
        failures.append(f"missing B mutation gates: {', '.join(sorted(missing_b))}")
    for name, passed in b_mutations.items():
        if not passed:
            failures.append(f"B mutation gate failed: {name}")
    missing_exact = REQUIRED_EXACT_GATES.difference(exact_gates)
    if missing_exact:
        failures.append(f"missing exact theorem gates: {', '.join(sorted(missing_exact))}")
    for name, passed in exact_gates.items():
        if not passed:
            failures.append(f"exact theorem gate failed: {name}")
    return PrecisionDecision("PASS" if not failures else "UNRESOLVED", tuple(failures))


def next_distinct_a_rate(eigenvalues, cluster_ritz, *, rate_tolerance):
    """Remove the continued rate class from a stationary-projected spectrum."""
    tolerance = mp.mpf(rate_tolerance)
    if not mp.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("rate_tolerance must be finite and positive")
    values = [mp.mpc(value) for value in eigenvalues]
    targets = [mp.mpc(value) for value in cluster_ritz]
    targets.extend(mp.conj(value) for value in cluster_ritz)
    rows = _minimum_cost_target_assignment_mp(values, targets)
    if len(rows) != 4 or max(
        abs(values[row] - targets[column])
        for column, row in enumerate(rows)
    ) > tolerance:
        raise AssertionError("continued conjugate rate class was not resolved")
    remaining = [
        value for index, value in enumerate(values) if index not in set(rows)
    ]
    cluster_rate = -mp.re(sum(targets[:2]) / 2)
    rates = [-mp.re(value) for value in remaining]
    if any(rate < -tolerance for rate in rates):
        raise AssertionError("unstable A mode remains after exact stationary projection")
    if any(rate <= tolerance for rate in rates):
        raise AssertionError("stationary A mode remains after exact stationary projection")
    if any(abs(rate - cluster_rate) <= tolerance for rate in rates):
        raise AssertionError("continued rate class was not completely removed")
    positive = sorted(rates)
    if not positive:
        raise AssertionError("no positive complementary A rate remains")
    return float(positive[0])


def _minimum_cost_target_assignment_mp(values, targets):
    """Assign a small target list without demoting costs to float64."""
    full_mask = (1 << len(targets)) - 1
    states = {0: (mp.mpf(0), ())}
    for value_index, value in enumerate(values):
        updated = dict(states)
        for mask, (cost, chosen) in states.items():
            for target_index, target in enumerate(targets):
                bit = 1 << target_index
                if mask & bit:
                    continue
                candidate = (
                    cost + abs(value - target),
                    chosen + ((target_index, value_index),),
                )
                prior = updated.get(mask | bit)
                if prior is None or candidate[0] < prior[0]:
                    updated[mask | bit] = candidate
        states = updated
    if full_mask not in states:
        raise AssertionError("not enough eigenvalues for target assignment")
    assignment = [None] * len(targets)
    for target_index, value_index in states[full_mask][1]:
        assignment[target_index] = value_index
    return tuple(assignment)


def structured_spectrum_mp(build, epsilon_token):
    """Read the complete spectrum from mp-native 1+6+6+36 blocks."""
    with mp.workprec(build.bits):
        _, blocks = _mp_liouville_blocks(build.matrix, epsilon_token, build.bits)
        values = []
        for block in blocks:
            block_values = mp.eig(block, left=False, right=False)
            if isinstance(block_values, tuple):
                block_values = block_values[0]
            values.extend(block_values)
        if len(values) != 49:
            raise AssertionError("structured spectrum did not contain 49 modes")
        return tuple(values)


def _stationary_complement_blocks_mp(build, epsilon_token):
    """Remove the exact two-dimensional stationary space before eigensolving."""
    with mp.workprec(build.bits):
        unitary, blocks = _mp_liouville_blocks(
            build.matrix, epsilon_token, build.bits
        )
        outer_identity = mp.zeros(36, 1)
        for index in range(6):
            outer_identity[index * 6 + index] = 1
        outer_identity /= mp.sqrt(6)
        outer_complement, lower_left = _local_complement(
            blocks[3], outer_identity, build.bits
        )
        stationary_basis = mp.matrix(49, 2)
        stationary_basis[:, 0] = unitary[:, 0]
        stationary_basis[:, 1] = unitary[:, 13:49] * outer_identity
        kernel_residual = spectral_norm_mp(build.matrix * stationary_basis)
        evidence = {
            "exact_stationary_dimension": 2,
            "kernel_residual": float(kernel_residual),
            "lower_left_residual": float(spectral_norm_mp(lower_left)),
        }
        return (blocks[1], blocks[2], outer_complement), evidence


def stationary_complement_spectrum_mp(build, epsilon_token):
    """Return all 47 nonstationary eigenvalues after exact-space projection."""
    with mp.workprec(build.bits):
        blocks, evidence = _stationary_complement_blocks_mp(build, epsilon_token)
        values = []
        for block in blocks:
            block_values = mp.eig(block, left=False, right=False)
            if isinstance(block_values, tuple):
                block_values = block_values[0]
            values.extend(block_values)
        if len(values) != 47:
            raise AssertionError("stationary complement spectrum did not contain 47 modes")
        return tuple(values), evidence


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
    bits = mp.mp.prec
    with mp.workprec(bits):
        _, blocks = _mp_liouville_blocks(matrix, state.epsilon_token, bits)
        values = []
        for block in blocks:
            block_values = mp.eig(block, left=False, right=False)
            if isinstance(block_values, tuple):
                block_values = block_values[0]
            values.extend(block_values)
        best = min(
            (
                max(
                    abs(values[first] - state.right_ritz[0]),
                    abs(values[second] - state.right_ritz[1]),
                ),
                first,
                second,
            )
            for first in range(len(values))
            for second in range(len(values))
            if first != second
        )
        selected_indices = {best[1], best[2]}
        selected = [values[index] for index in selected_indices]
        complement = [
            value for index, value in enumerate(values)
            if index not in selected_indices
        ]
        return min(abs(left - right) for left in selected for right in complement)


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


def _k_eigenvalue_diagnostic(build, target, build_discrepancy):
    """Bound one K eigenvalue using its own residuals and p/2p drift."""
    with mp.workprec(build.bits):
        values, vectors = _eigenpairs(build.matrix)
        index = min(range(len(values)), key=lambda item: abs(values[item] - target))
        value = values[index]
        right = vectors[index] / mp.norm(vectors[index])
        sep_complex = min(
            abs(value - candidate)
            for candidate_index, candidate in enumerate(values)
            if candidate_index != index
        )
        left_value, left = _nearest_eigenpair(
            build.matrix.transpose_conj(), mp.conj(value)
        )
        sigma = abs((left.transpose_conj() * right)[0])
        if sigma == 0:
            raise AssertionError("K left/right eigenvectors are biorthogonally singular")
        kappa = 1 / sigma
        right_residual = mp.norm(build.matrix * right - value * right)
        left_residual = mp.norm(
            build.matrix.transpose_conj() * left - left_value * left
        )
        left_right_pair_residual = abs(left_value - mp.conj(value))
        rounding = mp.eps * mp.norm(build.matrix) * 128
        error = kappa * (
            max(right_residual, left_residual, left_right_pair_residual)
            + build_discrepancy
            + rounding
        )
        return value, {
            "right_residual": float(right_residual),
            "left_residual": float(left_residual),
            "left_right_pair_residual": float(left_right_pair_residual),
            "biorthogonal_sigma_min": float(sigma),
            "condition_number": float(kappa),
            "construction_discrepancy": float(build_discrepancy),
            "rounding_term": float(rounding),
            "eigenvalue_error": float(error),
            "gap": float(-mp.re(value)),
            "sep_complex": float(sep_complex),
        }


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
        t_cluster = state.right_basis.transpose_conj() * matrix * state.right_basis
        sylvester_values = []
        for block in complement_blocks:
            sylvester_values.append(sylvester_separation_mp(t_cluster, block))
        sep_sylvester = min(sylvester_values)
        transformed = unitary.transpose_conj() * matrix * unitary
        block_diagonal = mp.zeros(49)
        for start, block in zip((0, 1, 7, 13), blocks):
            block_diagonal[start:start + block.rows, start:start + block.cols] = block
        off_block = transformed - block_diagonal
        off_block_residual = mp.norm(off_block)
        return {
            "block_dimensions": [block.rows for block in complement_blocks],
            "sep_complex": sep_complex,
            "sep_sylvester": sep_sylvester,
            "t_cluster": t_cluster,
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
        right_values, right_vectors = _cross_block_candidates(
            matrix, epsilon_token, bits, previous.right_ritz, candidate_radius
        )
        candidate_count = len(right_values)
        if candidate_count < 2 or candidate_count > 4:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": candidate_count,
                "status": "UNRESOLVED",
                "failures": ["candidate union did not contain between two and four eigenvalues"],
            }
        chosen_right = _choose_overlap_candidate(
            matrix, right_values, right_vectors, previous.right_basis, bits, candidate_radius
        )
        if chosen_right is None:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": candidate_count,
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
        if len(left_values) < 2 or len(left_values) > 4:
            return None, {
                "epsilon_token": epsilon_token,
                "candidate_count": candidate_count,
                "left_candidate_count": len(left_values),
                "status": "UNRESOLVED",
                "failures": ["left candidate union did not contain between two and four eigenvalues"],
            }
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
                "candidate_count": candidate_count,
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
                "candidate_count": candidate_count,
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
                "candidate_count": candidate_count,
                "status": "UNRESOLVED",
                "failures": ["principal angle reached pi/2"],
            }
        return candidate, {
            "epsilon_token": epsilon_token,
            "candidate_count": candidate_count,
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


def _a_residual_components(build, state, build_error):
    with mp.workprec(build.bits):
        matrix = build.matrix
        right = state.right_basis
        left = state.left_basis
        identity = mp.eye(N * N)
        right_projector = _orthogonal_projector(right)
        left_projector = _orthogonal_projector(left)
        right_residual = spectral_norm_mp(
            (identity - right_projector) * matrix * right
        )
        left_residual = spectral_norm_mp(
            (identity - left_projector) * matrix.transpose_conj() * left
        )
        sigma = _smallest_singular_value(left.transpose_conj() * right)
        kappa = 1 / sigma
        rounding_term = mp.eps * mp.norm(matrix) * 256
        eigenvalue_error = kappa * (
            max(right_residual, left_residual) + build_error + rounding_term
        )
        return {
            "right_cluster_residual": float(right_residual),
            "left_cluster_residual": float(left_residual),
            "biorthogonal_sigma_min": float(sigma),
            "condition_number": float(kappa),
            "construction_discrepancy": float(build_error),
            "rounding_term": float(rounding_term),
            "eigenvalue_error": float(eigenvalue_error),
        }


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
        right = state.right_basis
        left = state.left_basis
        pr = _orthogonal_projector(right)
        pl = _orthogonal_projector(left)
        components = _a_residual_components(build, state, build_error)
        eigen_error = components["eigenvalue_error"]
        gap = -mp.re(sum(state.right_ritz) / 2)
        sep_complex = mp.mpf(sep_complex)
        sep_sylvester = mp.mpf(sep_sylvester)
        subspace_error = eigen_error / sep_sylvester
        light = min(_centre_light(pr), _centre_light(pl))
        diagnostic = PrecisionDiagnostics(
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
        return diagnostic, components


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
        unitarity_residual = mp.norm(unitarity)
        off_block = transformed - block_diagonal
        off_block_residual = mp.norm(off_block)
        # This is an operator reconstruction residual.  The separately stored
        # lifted-eigenpair residual is the only spectral read made here; no
        # forward nonnormal spectral perturbation bound is claimed.
        union_residual = mp.norm(off_block)
        eigen_residuals = []
        for block, start in zip(blocks, starts):
            values, vectors = _eigenpairs(block)
            carrier = unitary[:, start:start + block.rows]
            for value, vector in zip(values, vectors):
                lifted = carrier * vector
                lifted /= mp.norm(lifted)
                eigen_residuals.append(
                    spectral_norm_mp(build.matrix * lifted - value * lifted)
                )
        return {
            "dimensions": [block.rows for block in blocks],
            "blocks": tuple(blocks),
            "unitarity_residual": float(unitarity_residual),
            "off_block_residual": float(off_block_residual),
            "full_operator_reconstruction_residual": float(union_residual),
            "full_eigen_residual": float(max(eigen_residuals)),
        }


def measure_precision_pair(
    epsilon_token,
    gamma_token,
    low_bits,
    high_bits,
    *,
    low_state=None,
    high_state=None,
    low_track_record=None,
    high_track_record=None,
):
    """Measure a fresh p/2p doubled cluster and all trust-contract quantities."""
    target_fraction = Fraction(epsilon_token)
    sign = 1 if target_fraction > 0 else -1
    support_tokens = tuple(dict.fromkeys((
        _fraction_token(Fraction(sign, 10)),
        _fraction_token(Fraction(sign, 8)),
        epsilon_token,
    )))
    local_paths = _branch_paths(support_tokens)

    def tracked_at(bits):
        records, states = prepare_continuation_paths(local_paths, gamma_token, bits)
        record = next(
            step for steps in records.values() for step in steps
            if step["epsilon_token"] == epsilon_token
        )
        return states[epsilon_token], record

    if low_state is None or low_track_record is None:
        low_state, low_track_record = tracked_at(low_bits)
    if high_state is None or high_track_record is None:
        high_state, high_track_record = tracked_at(high_bits)
    low_build = build_a_mp(epsilon_token, gamma_token, low_bits)
    high_build = build_a_mp(epsilon_token, gamma_token, high_bits)
    low_k_build = build_k_mp(epsilon_token, gamma_token, low_bits)
    high_k_build = build_k_mp(epsilon_token, gamma_token, high_bits)
    with mp.workprec(high_bits + 32):
        k_build_discrepancy = precision_build_discrepancy(
            low_k_build, high_k_build
        )
    with mp.workprec(low_bits):
        low_k_value, low_k_components = _k_eigenvalue_diagnostic(
            low_k_build,
            mp.mpc(0, 2 * mp.sqrt(2)),
            k_build_discrepancy,
        )
    with mp.workprec(high_bits):
        high_k_value, high_k_components = _k_eigenvalue_diagnostic(
            high_k_build,
            mp.mpc(0, 2 * mp.sqrt(2)),
            k_build_discrepancy,
        )
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
    low_adjoint = build_a_adjoint_mp(epsilon_token, gamma_token, low_bits)
    high_adjoint = build_a_adjoint_mp(epsilon_token, gamma_token, high_bits)
    low_left_state = ClusterState(
        epsilon_token, 2, low_state.left_basis, low_state.right_basis,
        low_state.left_ritz, low_state.right_ritz,
    )
    high_left_state = ClusterState(
        epsilon_token, 2, high_state.left_basis, high_state.right_basis,
        high_state.left_ritz, high_state.right_ritz,
    )
    low_right_separation = native_block_separation(
        low_build.matrix, low_state, epsilon_token, low_bits
    )
    low_left_separation = native_block_separation(
        low_adjoint.matrix, low_left_state, epsilon_token, low_bits
    )
    high_right_separation = native_block_separation(
        high_build.matrix, high_state, epsilon_token, high_bits
    )
    high_left_separation = native_block_separation(
        high_adjoint.matrix,
        high_left_state,
        epsilon_token,
        high_bits,
    )
    low_sep_complex = min(
        low_right_separation["sep_complex"], low_left_separation["sep_complex"]
    )
    low_sep_sylvester = min(
        low_right_separation["sep_sylvester"], low_left_separation["sep_sylvester"]
    )
    high_sep_complex = min(
        high_right_separation["sep_complex"], high_left_separation["sep_complex"]
    )
    high_sep_sylvester = min(
        high_right_separation["sep_sylvester"], high_left_separation["sep_sylvester"]
    )
    projector_difference = max(right_difference, left_difference)
    light_difference = max(right_light_difference, left_light_difference)
    low_diagnostic, low_components = _diagnostic_at_precision(
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
    high_diagnostic, high_components = _diagnostic_at_precision(
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
        operator_map_residual = (
            mp.mpf(0) if not any(operator_delta)
            else spectral_norm_mp(operator_delta)
        )
        operator_map_bound = (
            precision_build_discrepancy(low_build, high_build)
            + precision_build_discrepancy(low_b_build, b_build)
        )
    flip_pairs = tuple(
        (
            build_global_flip_carrier_mp(
                epsilon_token, gamma_token, low_bits, carrier, copy
            ),
            build_global_flip_carrier_mp(
                epsilon_token, gamma_token, high_bits, carrier, copy
            ),
        )
        for carrier, copy in (
            ("A", "original"),
            ("A", "global_flip"),
            ("B", "original"),
            ("B", "global_flip"),
        )
    )
    (_, a_original), (_, a_copy), (_, b_original), (_, b_copy) = flip_pairs
    with mp.workprec(high_bits):
        a_copy_delta = a_original.matrix - a_copy.matrix
        b_copy_delta = b_original.matrix - b_copy.matrix
        copy_residual = max(
            mp.mpf(0) if not any(a_copy_delta) else spectral_norm_mp(a_copy_delta),
            mp.mpf(0) if not any(b_copy_delta) else spectral_norm_mp(b_copy_delta),
        )
        copy_bounds = []
        for low_flip, high_flip in flip_pairs:
            copy_bounds.append(mp.norm(mp.matrix(low_flip.matrix) - high_flip.matrix))
        copy_bound = max(
            copy_bounds[0] + copy_bounds[1],
            copy_bounds[2] + copy_bounds[3],
        )
    decomposition_gates = []
    for current_build in (low_build, high_build):
        decomposition = carrier_block_decomposition(current_build, epsilon_token)
        with mp.workprec(current_build.bits):
            bound = (
                8192 * mp.power(2, -current_build.bits)
                * max(mp.mpf(1), mp.norm(current_build.matrix))
            )
        decomposition_gates.append({
            "precision_bits": current_build.bits,
            "unitarity_residual": decomposition["unitarity_residual"],
            "off_block_residual": decomposition["off_block_residual"],
            "full_operator_reconstruction_residual": decomposition["full_operator_reconstruction_residual"],
            "full_eigen_residual": decomposition["full_eigen_residual"],
            "bound": float(bound),
        })
    with mp.workprec(high_bits + 32):
        k_pair_error = abs(mp.mpc(low_k_value) - mp.mpc(high_k_value))
    k_pair_bound = (
        low_k_components["eigenvalue_error"]
        + high_k_components["eigenvalue_error"]
    )
    with mp.workprec(low_bits + 32):
        low_a_k_agreement = max(
            abs(value - low_k_value) for value in low_state.right_ritz
        )
    with mp.workprec(high_bits + 32):
        high_a_k_agreement = max(
            abs(value - high_k_value) for value in high_state.right_ritz
        )
    low_a_k_bound = (
        low_diagnostic.eigenvalue_error
        + low_k_components["eigenvalue_error"]
    )
    high_a_k_bound = (
        high_diagnostic.eigenvalue_error
        + high_k_components["eigenvalue_error"]
    )
    _, low_kernel_evidence = _stationary_complement_blocks_mp(
        low_build, epsilon_token
    )
    a_values, high_kernel_evidence = stationary_complement_spectrum_mp(
        high_build, epsilon_token
    )
    with mp.workprec(high_bits):
        spectrum_scale = max(mp.mpf(1), *(abs(value) for value in a_values))
        rate_tolerance = max(
            mp.mpf(8) * high_diagnostic.eigenvalue_error,
            mp.mpf(8) * high_diagnostic.eigenvalue_pair_error,
            mp.power(2, -high_bits + 8) * spectrum_scale,
        )
        next_rate = next_distinct_a_rate(
            a_values, high_state.right_ritz, rate_tolerance=rate_tolerance
        )
    result = {
        "low": low_diagnostic,
        "high": high_diagnostic,
        "low_diagnostic_components": low_components,
        "high_diagnostic_components": high_components,
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
        "global_flip_copy_bound": float(copy_bound),
        "b_builds": (low_b_build, b_build),
        "k_builds": (low_k_build, high_k_build),
        "low_k_eigenvalue": low_k_value,
        "high_k_eigenvalue": high_k_value,
        "k_eigenvalue_pair_error": float(k_pair_error),
        "k_eigenvalue_pair_bound": float(k_pair_bound),
        "k_build_discrepancy": float(k_build_discrepancy),
        "low_k_diagnostic_components": low_k_components,
        "high_k_diagnostic_components": high_k_components,
        "low_a_k_eigenvalue_agreement": float(low_a_k_agreement),
        "high_a_k_eigenvalue_agreement": float(high_a_k_agreement),
        "low_a_k_eigenvalue_agreement_bound": float(low_a_k_bound),
        "high_a_k_eigenvalue_agreement_bound": float(high_a_k_bound),
        "next_distinct_a_rate": float(next_rate),
        "next_distinct_a_rate_tolerance": float(rate_tolerance),
        "low_exact_stationary_dimension": low_kernel_evidence["exact_stationary_dimension"],
        "high_exact_stationary_dimension": high_kernel_evidence["exact_stationary_dimension"],
        "low_kernel_residual": low_kernel_evidence["kernel_residual"],
        "high_kernel_residual": high_kernel_evidence["kernel_residual"],
        "low_kernel_lower_left_residual": low_kernel_evidence["lower_left_residual"],
        "high_kernel_lower_left_residual": high_kernel_evidence["lower_left_residual"],
        "global_flip_build_pairs": flip_pairs,
        "low_build": low_build,
        "high_build": high_build,
        "low_state": low_state,
        "high_state": high_state,
        "state": high_state,
        "native_block_dimensions": high_right_separation["block_dimensions"],
        "right_native_off_block_residual": float(high_right_separation["off_block_residual"]),
        "left_native_off_block_residual": float(high_left_separation["off_block_residual"]),
        "right_native_lower_left_residual": float(high_right_separation["lower_left_residual"]),
        "left_native_lower_left_residual": float(high_left_separation["lower_left_residual"]),
        "separation_calls": [
            [low_bits, "right"], [low_bits, "left"],
            [high_bits, "right"], [high_bits, "left"],
        ],
        "decomposition_gates": decomposition_gates,
        "continuation_records": [low_track_record, high_track_record],
    }
    _register_measurement_result(result)
    result["decision"] = decide_certificate_attempt(result)
    result["status"] = result["decision"].status
    _seal_measurement_result(result)
    return result


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
    if (
        not _measurement_result_is_authentic(result)
        or result.get("status") != "TRUSTED"
        or not isinstance(result.get("decision"), PrecisionDecision)
        or result["decision"].status != "TRUSTED"
    ):
        raise AssertionError("measurement provenance is not sealed as TRUSTED")
    state = result["state"]
    gap = result["high"].gap
    next_rate = result["next_distinct_a_rate"]
    rate_tolerance = result["next_distinct_a_rate_tolerance"]
    b_values = np.linalg.eigvals(_mp_to_numpy(build_b_mp(epsilon_token, gamma_token, high_bits).matrix))
    return {
        "epsilon_token": epsilon_token,
        "gamma_token": gamma_token,
        "precision_bits": [low_bits, high_bits],
        "K_gap": _decimal(-mp.re(result["high_k_eigenvalue"])),
        "low_K_gap": _decimal(-mp.re(result["low_k_eigenvalue"])),
        "K_eigenvalue_pair_error": _decimal(result["k_eigenvalue_pair_error"]),
        "K_eigenvalue_pair_bound": _decimal(result["k_eigenvalue_pair_bound"]),
        "K_build_discrepancy": _decimal(result["k_build_discrepancy"]),
        "low_K_diagnostic_components": {
            key: _decimal(value)
            for key, value in result["low_k_diagnostic_components"].items()
        },
        "high_K_diagnostic_components": {
            key: _decimal(value)
            for key, value in result["high_k_diagnostic_components"].items()
        },
        "low_A_K_eigenvalue_agreement": _decimal(
            result["low_a_k_eigenvalue_agreement"]
        ),
        "high_A_K_eigenvalue_agreement": _decimal(
            result["high_a_k_eigenvalue_agreement"]
        ),
        "low_A_K_eigenvalue_agreement_bound": _decimal(
            result["low_a_k_eigenvalue_agreement_bound"]
        ),
        "high_A_K_eigenvalue_agreement_bound": _decimal(
            result["high_a_k_eigenvalue_agreement_bound"]
        ),
        "build_registry_gates": {
            "A_p2p_independent": precision_matrices_are_independent(
                result["low_build"], result["high_build"]
            ),
            "K_p2p_independent": precision_matrices_are_independent(
                *result["k_builds"]
            ),
            "B_p2p_independent": precision_matrices_are_independent(
                *result["b_builds"]
            ),
            "global_flip_p2p_independent": all(
                precision_matrices_are_independent(*pair)
                for pair in result["global_flip_build_pairs"]
            ),
        },
        "measurement_registry_bound": _measurement_result_is_authentic(result),
        "A_gap": _decimal(gap),
        "cluster_rank": state.rank,
        "right_light": _decimal(result["lights"]["right"]),
        "left_light": _decimal(result["lights"]["left"]),
        "right_absorption_residual": _decimal(result["absorption_residuals"]["right"]),
        "left_absorption_residual": _decimal(result["absorption_residuals"]["left"]),
        "absorption_error": _decimal(result["absorption_error"]),
        "low_diagnostics": _serialize_diagnostic(result["low"]),
        "high_diagnostics": _serialize_diagnostic(result["high"]),
        "low_diagnostic_components": {
            key: _decimal(value)
            for key, value in result["low_diagnostic_components"].items()
        },
        "high_diagnostic_components": {
            key: _decimal(value)
            for key, value in result["high_diagnostic_components"].items()
        },
        "right_projector_difference": _decimal(result["right_projector_difference"]),
        "left_projector_difference": _decimal(result["left_projector_difference"]),
        "right_light_difference": _decimal(result["right_light_difference"]),
        "left_light_difference": _decimal(result["left_light_difference"]),
        "next_distinct_A_rate": _decimal(next_rate),
        "next_distinct_A_rate_tolerance": _decimal(rate_tolerance),
        "low_exact_stationary_dimension": result["low_exact_stationary_dimension"],
        "high_exact_stationary_dimension": result["high_exact_stationary_dimension"],
        "low_kernel_residual": _decimal(result["low_kernel_residual"]),
        "high_kernel_residual": _decimal(result["high_kernel_residual"]),
        "low_kernel_lower_left_residual": _decimal(
            result["low_kernel_lower_left_residual"]
        ),
        "high_kernel_lower_left_residual": _decimal(
            result["high_kernel_lower_left_residual"]
        ),
        "B_gap": _decimal(min(-b_values.real)),
        "B_operator_map_residual": _decimal(result["b_operator_map_residual"]),
        "B_operator_map_bound": _decimal(result["b_operator_map_bound"]),
        "full_B_match_numerical_read": _decimal(result["full_b_match_residual"]),
        "full_B_match_backward_error_read": _decimal(result["full_b_match_error"]),
        "global_flip_copy_residual": _decimal(result["global_flip_copy_residual"]),
        "global_flip_copy_bound": _decimal(result["global_flip_copy_bound"]),
        "decomposition_gates": [
            {
                key: (value if key == "precision_bits" else _decimal(value))
                for key, value in gate.items()
            }
            for gate in result["decomposition_gates"]
        ],
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
            (
                value
                for value in values
                if value * sign > 0
                and abs(value) < Fraction(1, 8)
                and value != historical
            ),
            key=abs,
            reverse=True,
        ))
        paths[label] = [_fraction_token(value) for value in branch]
    return paths


def _expected_effective_operators_for_certificate(gamma):
    root2 = sympy.sqrt(2)
    return {
        "first": {
            "-4sqrt2i": sympy.Matrix([[-sympy.I * root2]]),
            "-2sqrt2i": -sympy.I / root2 * sympy.eye(2),
            "0": sympy.zeros(4),
            "+2sqrt2i": sympy.I / root2 * sympy.eye(2),
            "+4sqrt2i": sympy.Matrix([[sympy.I * root2]]),
        },
        "second": {
            "-4sqrt2i": sympy.Matrix([[-gamma - 3 * root2 * sympy.I / 8]]),
            "-2sqrt2i": (-gamma / 2 - 3 * root2 * sympy.I / 16) * sympy.eye(2),
            "0": sympy.Matrix([
                [-gamma, 0, 0, gamma / 2],
                [0, 0, 0, 0],
                [0, 0, -gamma, gamma / 2],
                [gamma / 2, 0, gamma / 2, -gamma / 2],
            ]),
            "+2sqrt2i": (-gamma / 2 + 3 * root2 * sympy.I / 16) * sympy.eye(2),
            "+4sqrt2i": sympy.Matrix([[-gamma + 3 * root2 * sympy.I / 8]]),
        },
    }


def _serialize_sympy_matrix(matrix):
    return [[str(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)]


def exact_top_level_gates():
    """Execute every exact claim serialized by the certificate header."""
    lam, r, gamma = sympy.symbols("lambda r gamma")
    determinant = sympy.expand(characteristic_polynomial(lam, r, gamma))
    displayed_q = (
        lam**6 + 2 * gamma * lam**5 + (4 * r**2 + 20) * lam**4
        + (8 * gamma * r**2 + 24 * gamma) * lam**3
        + (64 * r**2 + 96) * lam**2
        + (64 * gamma * r**2 + 64 * gamma) * lam
        + 192 * r**2 + 64
    )
    polynomial_match = sympy.expand(determinant - lam * displayed_q) == 0

    epsilon = sympy.symbols("epsilon", real=True)
    gamma_real = sympy.symbols("gamma_positive", positive=True, real=True)
    branch = minus_branch_series(epsilon, gamma_real)
    substituted = sympy.expand(
        characteristic_polynomial(lam, r, gamma_real).subs(
            {lam: branch, r: 1 + epsilon}
        )
    )
    series_gate = all(
        sympy.simplify(substituted.coeff(epsilon, order)) == 0
        for order in range(4)
    ) and sympy.simplify(gap_series(epsilon, gamma_real) + sympy.re(branch)) == 0

    effective = exact_effective_operators(gamma_real)
    expected_effective = _expected_effective_operators_for_certificate(gamma_real)
    effective_gate = (
        all(
            effective[order][frequency] == expected_matrix
            and effective["characteristic_polynomials"][order][frequency]
            == expected_matrix.charpoly(
                effective["characteristic_polynomials"][order][frequency].gen
            )
            for order in ("first", "second")
            for frequency, expected_matrix in expected_effective[order].items()
        )
        and all(
            certificate["denominator_nonzero_for_positive_gamma"]
            and certificate["shifted_off_block_residual"] == sympy.zeros(49)
            for certificate in effective["certificates"].values()
        )
    )

    uniform_certificate = exact_uniform_peripheral_certificate(
        sympy.Rational(3, 10)
    )
    uniform = uniform_certificate["census"]
    punctured = exact_punctured_kernel_certificate(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    punctured_polynomial = punctured["restricted_characteristic_polynomial"]
    boundary = exact_b_boundary_certificate()
    denominator_gates = {
        key: {
            "denominator": str(certificate["resolvent_denominator"]),
            "zero_root_multiplicity": int(
                certificate["positive_denominator_certificate"][
                    "zero_root_multiplicity"
                ]
            ),
            "positive_axis_polynomial": str(
                certificate["positive_denominator_certificate"][
                    "positive_axis_polynomial"
                ].as_expr()
            ),
            "positive_real_root_count": int(
                certificate["positive_denominator_certificate"][
                    "positive_real_root_count"
                ]
            ),
        }
        for key, certificate in effective["certificates"].items()
    }
    values = {
        "polynomial_match": bool(polynomial_match),
        "characteristic_polynomial": str(determinant),
        "uniform_peripheral_dimension": int(uniform["peripheral_dimension"]),
        "uniform_kernel_dimension": int(uniform["kernel_dimension"]),
        "frequency_multiplicities": {
            key: int(value)
            for key, value in uniform["frequency_multiplicities"].items()
        },
        "uniform_combined_operator_rank": int(
            uniform_certificate["combined_operator_rank"]
        ),
        "uniform_zero_cost_invariant_dimension": int(
            uniform_certificate["zero_cost_invariant_dimension"]
        ),
        "effective_denominator_gates": denominator_gates,
        "effective_operators": {
            order: {
                frequency: _serialize_sympy_matrix(matrix)
                for frequency, matrix in effective[order].items()
            }
            for order in ("first", "second")
        },
        "effective_characteristic_polynomials": {
            order: {
                frequency: str(polynomial.as_expr())
                for frequency, polynomial in
                effective["characteristic_polynomials"][order].items()
            }
            for order in ("first", "second")
        },
        "punctured_kernel_dimension": int(punctured["kernel_dimension"]),
        "punctured_stationary_rank": int(punctured["stationary_rank"]),
        "punctured_zero_cost_invariant_dimension": int(
            punctured["zero_cost_invariant_dimension"]
        ),
        "punctured_nonzero_imaginary_axis_dimension": int(
            punctured["nonzero_imaginary_axis_dimension"]
        ),
        "punctured_restricted_characteristic_polynomial": str(
            punctured["restricted_characteristic_polynomial"].as_expr()
        ),
        "b_boundary_has_peripheral_mode": bool(boundary["has_peripheral_b_mode"]),
        "B_boundary_peripheral_dimension": int(boundary["peripheral_dimension"]),
        "B_boundary_h_outer_times_u": [
            str(value) for value in boundary["h_outer_times_u"]
        ],
    }
    checks = {
        "polynomial_match": values["polynomial_match"],
        "gap_series_match": bool(series_gate),
        "effective_operator_match": bool(effective_gate),
        "uniform_peripheral_dimension_match": values["uniform_peripheral_dimension"] == 10,
        "uniform_kernel_dimension_match": values["uniform_kernel_dimension"] == 4,
        "punctured_kernel_dimension_match": (
            values["punctured_kernel_dimension"] == 2
            and values["punctured_stationary_rank"] == 2
            and values["punctured_zero_cost_invariant_dimension"] == 2
            and values["punctured_nonzero_imaginary_axis_dimension"] == 0
            and punctured_polynomial.as_expr() == punctured_polynomial.gen**2
        ),
        "b_boundary_match": (
            not values["b_boundary_has_peripheral_mode"]
            and values["B_boundary_peripheral_dimension"] == 0
        ),
    }
    return checks, values


def producer_provenance():
    source = Path(__file__).read_text(encoding="utf-8").replace("\r\n", "\n")
    return {
        "producer_revision": "source-sha256-normalized-lf",
        "source_sha256_normalized_lf": hashlib.sha256(
            source.encode("utf-8")
        ).hexdigest(),
        "dependencies": {
            "python": platform.python_version(),
            "mpmath": mp.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sympy.__version__,
        },
        "exact_producers": [
            "characteristic_polynomial",
            "minus_branch_series",
            "exact_effective_operators",
            "exact_uniform_peripheral_certificate",
            "exact_punctured_kernel_certificate",
            "exact_b_boundary_certificate",
        ],
    }


def build_certificate(
    *,
    gamma_tokens=GAMMA_TOKENS,
    epsilon_tokens=EPSILON_TOKENS,
    measure=measure_precision_pair,
):
    exact_gates, exact_values = exact_top_level_gates()
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

    def track_record(records, epsilon_token):
        return next(
            step
            for steps in records.values()
            for step in steps
            if step["epsilon_token"] == epsilon_token
        )

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
                        low_records, low_states = prepared(gamma_token, low_bits)
                        high_records, high_states = prepared(gamma_token, high_bits)
                        result = measure(
                            epsilon_token,
                            gamma_token,
                            low_bits,
                            high_bits,
                            low_state=low_states[epsilon_token],
                            high_state=high_states[epsilon_token],
                            low_track_record=track_record(low_records, epsilon_token),
                            high_track_record=track_record(high_records, epsilon_token),
                        )
                    else:
                        result = measure(epsilon_token, gamma_token, low_bits, high_bits)
                    decision = decide_certificate_attempt(result)
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
    operator_mutation = b_operator_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    spectrum_mutation = b_spectrum_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    b_mutations = {
        "correct_operator_residual_zero": operator_mutation["correct"] == 0,
        "omitted_adjoint_entrywise_detected": operator_mutation["omitted_adjoint"] != 0,
        "correct_spectrum_residual_zero": spectrum_mutation["correct"] == 0,
        "missing_right_action_spectrum_detected": spectrum_mutation["missing_right_action"] != 0,
        "wrong_price_sign_spectrum_detected": spectrum_mutation["wrong_price_sign"] != 0,
    }
    top_decision = certificate_verdict(
        rows, continuation_failures, b_mutations, exact_gates
    )
    verdict = top_decision.status
    document = {
        "schema": "missing-phase-relaxation-scale/v1",
        "scope": {"N": 7, "sector": [1, 1], "seat": 3, "path": "r=1+epsilon"},
        "theorem": {
            "gap_series": "gamma/2*epsilon^2 + gamma/2*epsilon^3 + O(epsilon^4)",
            "tau_leading": "2/(gamma*epsilon^2)",
            "observable_lifetime_claimed": False,
            "full_liouvillian_gap_claimed": False,
        },
        "exact": exact_values,
        "provenance": producer_provenance(),
        "rows": rows,
        "controls": {
            "continuation_paths": paths,
            "continuation_records": continuation,
            "continuation_failures": continuation_failures,
            "verdict_failures": list(top_decision.failures),
            "candidate_radius_rule": "one third of preceding trusted sep_complex; never widened",
            "B_relation_gate": "fresh full operators: B + 2*gamma*I + A_dagger",
            "full_B_spectrum": "numerical read; certification follows from the operator gate",
            "B_mutations": b_mutations,
            "exact_gates": exact_gates,
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


def kronecker_mp(left, right):
    result = mp.zeros(left.rows * right.rows, left.cols * right.cols)
    for i in range(left.rows):
        for j in range(left.cols):
            for k in range(right.rows):
                for ell in range(right.cols):
                    result[i * right.rows + k, j * right.cols + ell] = left[i, j] * right[k, ell]
    return result


def sylvester_separation_mp(cluster, complement):
    operator = (
        kronecker_mp(mp.eye(complement.rows), cluster)
        - kronecker_mp(complement.transpose(), mp.eye(cluster.rows))
    )
    return min(mp.svd(operator, compute_uv=False))


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
