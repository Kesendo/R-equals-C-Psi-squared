#!/usr/bin/env python3
"""Run the conditional neural scalar-pairing and mode-transport gates."""

import numpy as np

if __package__:
    from .neural_palindrome import (
        exact_ensemble_census, fitted_offdiagonal_residual, is_involution,
        make_exact_network, partner_subspace_error, scalar_center_residual,
        spectral_pairing_error,
    )
else:
    from neural_palindrome import (
        exact_ensemble_census, fitted_offdiagonal_residual, is_involution,
        make_exact_network, partner_subspace_error, scalar_center_residual,
        spectral_pairing_error,
    )


def main():
    J = np.array([[-0.5, -0.25], [0.25, -0.25]])
    perm = np.array([1, 0])
    values = np.linalg.eigvals(J)
    assert is_involution(perm)
    assert scalar_center_residual(J, perm, 0.375) == 0.0
    assert spectral_pairing_error(values, 0.375) < 1e-13
    assert np.max(np.abs(values.imag)) > 0.2
    print("PASS exact complex dyadic two-seat gate")

    fixed = np.diag([-0.5, -0.25, -0.5])
    assert fitted_offdiagonal_residual(fixed, [1, 0, 2]) == 0.0
    assert scalar_center_residual(fixed, [1, 0, 2], 0.375) > 0.2
    print("PASS fixed-seat scalar rejection; fitted off-diagonal residual = 0")

    assert spectral_pairing_error([-0.5, -0.25, -0.25], 0.375) > 0.2
    print("PASS multiplicity rejection")
    assert spectral_pairing_error([-0.5 + 1j, -0.25 - 2j], 0.375) == 1.0
    print("PASS full-complex negative control")

    expected = [
        (0.5, 200, 24, 0), (1.5, 200, 110, 1), (3.0, 200, 149, 15),
        (5.0, 200, 159, 24), (10.0, 200, 167, 45),
    ]
    rows = exact_ensemble_census()
    print("Ensemble (alpha, scalar-pass, oscillatory, unstable), 200 seeds per alpha:")
    for row, target in zip(rows, expected):
        print(row)
        assert row == target, (row, target)
    assert rows == expected

    J, perm, s = make_exact_network()
    scalar = scalar_center_residual(J, perm, s)
    transport = partner_subspace_error(J, perm, s)
    assert scalar < 1e-13
    assert transport < 1e-8
    print(f"PASS exact network: scalar residual={scalar:.6e}, transport sine={transport:.6e}")

    degenerate = partner_subspace_error(-0.375 * np.eye(4), [1, 0, 3, 2], 0.375)
    assert degenerate < 1e-12
    print(f"PASS degenerate invariant-subspace transport: sine={degenerate:.6e}")
    wrong = partner_subspace_error(np.diag([-0.5, -0.5, -0.25, -0.25]), [1, 0, 3, 2], 0.375)
    assert wrong > 0.9
    print(f"PASS wrong-permutation transport rejection: sine={wrong:.6e}")
    print("ALL NEURAL TRANSLATION GATES PASS")


if __name__ == "__main__":
    main()
