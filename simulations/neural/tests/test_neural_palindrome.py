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
