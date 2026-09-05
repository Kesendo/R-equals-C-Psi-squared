import pathlib
import subprocess
import sys

import numpy as np
import pytest

NEURAL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NEURAL))

from neural_palindrome import (
    exact_ensemble_census,
    fitted_offdiagonal_residual,
    is_involution,
    make_exact_network,
    partner_subspace_error,
    scalar_center_residual,
    spectral_pairing_error,
)


def test_exact_complex_pair_has_zero_scalar_residual_and_pairing_error():
    J = np.array([[-0.5, -0.25], [0.25, -0.25]])
    perm = np.array([1, 0])
    assert is_involution(perm)
    assert scalar_center_residual(J, perm, s=0.375) == 0.0
    assert spectral_pairing_error(np.linalg.eigvals(J), s=0.375) < 1e-13
    assert np.max(np.abs(np.linalg.eigvals(J).imag)) > 0.2


def test_fixed_seat_is_rejected_by_scalar_residual_but_not_old_fit():
    J = np.diag([-0.5, -0.25, -0.5])
    perm = np.array([1, 0, 2])
    assert fitted_offdiagonal_residual(J, perm) == 0.0
    assert scalar_center_residual(J, perm, s=0.375) > 0.2


def test_pairing_error_preserves_multiplicity():
    values = np.array([-0.5 + 0j, -0.25 + 0j, -0.25 + 0j])
    assert spectral_pairing_error(values, s=0.375) > 0.2


def test_pairing_error_keeps_imaginary_parts():
    values = np.array([-0.5 + 1j, -0.25 - 2j])
    # Both assigned differences are exactly -1j; real parts alone would pair.
    assert spectral_pairing_error(values, s=0.375) == 1.0


def test_exact_ensemble_census_reproduces_all_200_seed_rows():
    assert exact_ensemble_census() == [
        (0.5, 200, 24, 0),
        (1.5, 200, 110, 1),
        (3.0, 200, 149, 15),
        (5.0, 200, 159, 24),
        (10.0, 200, 167, 45),
    ]


def test_exact_network_has_scalar_identity_and_partner_subspaces():
    J, perm, s = make_exact_network(
        n=10, tau_e=5, tau_i=10, alpha=0.5, seed=42
    )
    assert scalar_center_residual(J, perm, s) < 1e-13
    assert partner_subspace_error(J, perm, s) < 1e-8


def test_exact_network_positional_seed_matches_keyword_seed():
    positional_J, positional_perm, positional_s = make_exact_network(10, 5, 10, 0.5, 0)
    keyword_J, keyword_perm, keyword_s = make_exact_network(
        n=10, tau_e=5, tau_i=10, alpha=0.5, seed=0
    )
    assert np.array_equal(positional_J, keyword_J)
    assert np.array_equal(positional_perm, keyword_perm)
    assert positional_s == keyword_s


@pytest.mark.parametrize("field,value", [
    ("n", 0), ("n", -2), ("n", 3), ("n", 10.0), ("n", True),
    ("tau_e", np.nan), ("tau_e", np.inf), ("tau_e", 0),
    ("tau_i", np.nan), ("tau_i", -np.inf), ("tau_i", -1),
    ("alpha", np.nan), ("alpha", np.inf), ("alpha", -np.inf),
    ("density", np.nan), ("density", -1), ("density", 2),
    ("seed", 0.5), ("seed", None), ("seed", True),
])
def test_exact_network_rejects_invalid_inputs(field, value):
    with pytest.raises(ValueError, match=field):
        make_exact_network(**{field: value})


def test_complex_transport_distinguishes_reflection_from_conjugated_target():
    J = np.array([[-0.5, -0.25], [0.25, -0.25]])
    assert partner_subspace_error(J, [1, 0], s=0.375) < 1e-12
    # Here the two roots are conjugates about -s. The incorrect target
    # -conj(mu)-2*s is mu itself, so it would falsely accept identity Q.
    assert partner_subspace_error(J, [0, 1], s=0.375) > 0.8


def test_degenerate_transport_compares_subspaces_not_individual_vectors():
    J = -0.375 * np.eye(4)
    perm = np.array([1, 0, 3, 2])
    assert partner_subspace_error(J, perm, s=0.375) < 1e-12


def test_transport_rejects_wrong_permutation_even_when_spectrum_pairs():
    J = np.diag([-0.5, -0.5, -0.25, -0.25])
    assert spectral_pairing_error(np.linalg.eigvals(J), s=0.375) == 0.0
    assert partner_subspace_error(J, [2, 3, 0, 1], s=0.375) < 1e-12
    assert partner_subspace_error(J, [1, 0, 3, 2], s=0.375) > 0.9


def test_transport_uses_algebraic_subspaces_for_defective_clusters():
    block = np.array([[-0.5, 0.125], [0.0, -0.5]])
    J = np.zeros((4, 4))
    J[:2, :2] = block
    J[2:, 2:] = -block - 0.75 * np.eye(2)
    assert partner_subspace_error(J, [2, 3, 0, 1], s=0.375) < 1e-12


def test_transport_rejects_missing_algebraic_partner_multiplicity():
    with pytest.raises(ValueError, match="algebraic cluster multiplicity"):
        partner_subspace_error(np.diag([-0.5, -0.25, -0.25]), [1, 0, 2], 0.375)


def test_transport_empty_matrix_has_no_angles():
    assert partner_subspace_error(np.empty((0, 0)), np.array([], dtype=int), 0.375) == 0.0


def test_translation_cli_runs_all_named_gates():
    result = subprocess.run(
        [sys.executable, str(NEURAL / "neural_translation_gate.py")],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stderr == ""
    assert "(10.0, 200, 167, 45)" in result.stdout
    assert result.stdout.splitlines()[-1] == "ALL NEURAL TRANSLATION GATES PASS"


def test_optimized_translation_cli_rejects_bad_scalar_reading():
    result = subprocess.run(
        [sys.executable, "-O", "-c",
         "import neural_translation_gate as gate; "
         "gate.scalar_center_residual = lambda *args: 1.0; gate.main()"],
        cwd=NEURAL, capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, result.stdout
    assert "PASS" not in result.stdout


def test_drive_fixed_point_reports_equation_residual():
    from neural_palindrome import make_balanced_dale_network, solve_fixed_point

    W, signs = make_balanced_dale_network(n=50, n_exc=25, density=0.3, seed=42)
    x, residual = solve_fixed_point(W, signs, alpha=0.3, drive=4,
                                    tol=1e-12, max_iter=5000)
    slopes = np.where(signs > 0, 1.3, 2.0)
    thresholds = np.where(signs > 0, 4.0, 3.7)
    rhs = 1 / (1 + np.exp(-slopes * (0.3 * W @ x + 4 - thresholds)))
    assert residual < 1e-10
    assert np.max(np.abs(x - rhs)) < 1e-10


def test_drive_fixed_point_rejects_unconverged_iterate():
    from neural_palindrome import make_balanced_dale_network, solve_fixed_point

    W, signs = make_balanced_dale_network(n=50, n_exc=25, seed=42)
    with pytest.raises(RuntimeError, match=r"last residual[=: ]+[0-9]"):
        solve_fixed_point(W, signs, alpha=0.3, drive=4,
                          tol=1e-12, max_iter=1)


@pytest.mark.parametrize("module", [
    "veffect_exact", "veffect_and_heat", "cpsi_two_perspectives",
    "celegans_trichotomy", "validation_checks", "find_quarter",
])
def test_producer_import_is_silent(module):
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"], cwd=NEURAL,
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""


@pytest.mark.parametrize("field,value", [
    ("n", 0), ("n", 3), ("n", True), ("n_exc", 24),
    ("density", np.nan), ("density", -1), ("seed", None),
])
def test_balanced_dale_builder_validates_inputs(field, value):
    from neural_palindrome import make_balanced_dale_network

    with pytest.raises(ValueError, match=field):
        make_balanced_dale_network(**{field: value})


@pytest.mark.parametrize("field,value", [
    ("W", np.ones((2, 3))), ("W", np.full((2, 2), np.nan)),
    ("signs", [1]), ("signs", [1, 0]),
    ("alpha", np.inf), ("drive", np.nan), ("tol", 0),
    ("max_iter", 0), ("max_iter", 1.5), ("max_iter", True),
])
def test_fixed_point_validates_inputs(field, value):
    from neural_palindrome import solve_fixed_point

    args = dict(W=np.zeros((2, 2)), signs=[1, -1], alpha=0.3,
                drive=4, tol=1e-12, max_iter=5000)
    args[field] = value
    with pytest.raises(ValueError, match=field):
        solve_fixed_point(**args)


@pytest.mark.parametrize("module", ["veffect_exact", "veffect_and_heat"])
def test_frequency_counts_share_explicit_resolution(module):
    import importlib

    producer = importlib.import_module(module)
    values = np.array([0.104j, -0.104j, 0.106j, -0.106j])
    activity, correlation = producer.frequency_counts(values, frequency_tolerance=0.01)
    assert activity == 2
    assert correlation == 1


def test_exact_producer_legacy_helpers_keep_balanced_linear_callers():
    from veffect_exact import (
        build_exact_palindromic_network, build_linear_jacobian,
        count_freqs, count_corr_freqs, palindrome_residual,
    )

    W, signs, _ = build_exact_palindromic_network(10, 5, 5, 10, seed=42)
    J = build_linear_jacobian(W, signs, 5, 10, 0.5)
    assert palindrome_residual(J, signs) == 0.0
    assert count_freqs(np.array([0.2j, -0.2j])) == 1
    assert count_corr_freqs(np.array([0.2j, -0.2j])) == 1


def test_legacy_residual_wrapper_rejects_odd_seat():
    from veffect_exact import palindrome_residual

    with pytest.raises(ValueError, match="scalar_center_residual"):
        palindrome_residual(np.diag([-0.5, -0.25, -0.5]), np.array([1, -1, 1]))


def test_legacy_residual_wrapper_rejects_nonuniform_type_diagonals():
    from veffect_exact import palindrome_residual

    with pytest.raises(ValueError, match="scalar_center_residual"):
        palindrome_residual(np.diag([-0.5, -0.25, -0.4, -0.25]),
                            np.array([1, -1, 1, -1]))


def test_fixed_point_rechecks_equation_after_small_expanding_update():
    from neural_palindrome import solve_fixed_point

    W = np.array([[0.0, 100.0], [100.0, 0.0]])
    drive = 4 + np.log(0.3001 / (1 - 0.3001)) / 1.3 - 30
    x, residual = solve_fixed_point(W, [1, 1], alpha=1, drive=drive,
                                    tol=0.001, max_iter=100)
    rhs = 1 / (1 + np.exp(-1.3 * (W @ x + drive - 4)))
    assert residual <= 0.001
    assert np.max(np.abs(x - rhs)) <= 0.001


def test_fixed_point_exhaustion_reports_fresh_equation_residual():
    from neural_palindrome import solve_fixed_point

    W = np.array([[0.0, 100.0], [100.0, 0.0]])
    drive = 4 + np.log(0.3001 / (1 - 0.3001)) / 1.3 - 30
    candidate = 1 / (1 + np.exp(-1.3 * (W @ np.full(2, 0.3) + drive - 4)))
    rhs = 1 / (1 + np.exp(-1.3 * (W @ candidate + drive - 4)))
    expected = np.max(np.abs(candidate - rhs))
    assert expected > 0.001
    with pytest.raises(RuntimeError, match="last residual=") as exc:
        solve_fixed_point(W, [1, 1], alpha=1, drive=drive,
                          tol=0.001, max_iter=1)
    reported = float(str(exc.value).split("last residual=")[1])
    assert abs(reported - expected) <= 4 * np.finfo(float).eps * expected


def test_crossing_value_uses_the_time_interpolation_fraction():
    from cpsi_two_perspectives import interpolated_value

    # An asymmetric crossing: the midpoint would be 0.375, not 0.5.
    assert interpolated_value(0.125, 0.625, fraction=0.75) == 0.5
    # Same path away from the crossing guards against returning constant 0.5.
    assert interpolated_value(0.125, 0.625, fraction=0.25) == 0.25
