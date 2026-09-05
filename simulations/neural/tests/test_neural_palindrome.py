import pathlib
import sys

import numpy as np

NEURAL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NEURAL))

from neural_palindrome import (
    fitted_offdiagonal_residual,
    is_involution,
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
