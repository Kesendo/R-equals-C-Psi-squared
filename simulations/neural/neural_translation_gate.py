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


def _require(condition, message):
    """Keep verdict checks active even when Python runs with optimization."""
    if not condition:
        raise RuntimeError(message)


def main():
    J = np.array([[-0.5, -0.25], [0.25, -0.25]])
    perm = np.array([1, 0])
    values = np.linalg.eigvals(J)
    _require(is_involution(perm), "two-seat permutation must be an involution")
    _require(scalar_center_residual(J, perm, 0.375) == 0.0,
             "dyadic scalar-center residual must be zero")
    _require(spectral_pairing_error(values, 0.375) < 1e-13,
             "dyadic complex spectrum must pair")
    _require(np.max(np.abs(values.imag)) > 0.2, "dyadic pair must be nonreal")
    print("PASS exact complex dyadic two-seat gate")

    fixed = np.diag([-0.5, -0.25, -0.5])
    _require(fitted_offdiagonal_residual(fixed, [1, 0, 2]) == 0.0,
             "fixed-seat fitted off-diagonal residual must be zero")
    _require(scalar_center_residual(fixed, [1, 0, 2], 0.375) > 0.2,
             "scalar-center residual must reject the fixed seat")
    print("PASS fixed-seat scalar rejection; fitted off-diagonal residual = 0")

    _require(spectral_pairing_error([-0.5, -0.25, -0.25], 0.375) > 0.2,
             "spectral pairing must reject mismatched multiplicity")
    print("PASS multiplicity rejection")
    _require(spectral_pairing_error([-0.5 + 1j, -0.25 - 2j], 0.375) == 1.0,
             "spectral pairing must retain imaginary parts")
    print("PASS full-complex negative control")

    expected = [
        (0.5, 200, 24, 0), (1.5, 200, 110, 1), (3.0, 200, 149, 15),
        (5.0, 200, 159, 24), (10.0, 200, 167, 45),
    ]
    rows = exact_ensemble_census()
    print("Ensemble (alpha, scalar-pass, oscillatory, unstable), 200 seeds per alpha:")
    for row, target in zip(rows, expected):
        print(row)
        _require(row == target, f"ensemble row {row} differs from {target}")
    _require(rows == expected, "ensemble census differs from the five expected rows")

    J, perm, s = make_exact_network()
    scalar = scalar_center_residual(J, perm, s)
    transport = partner_subspace_error(J, perm, s)
    _require(scalar < 1e-13, f"exact-network scalar residual too large: {scalar}")
    _require(transport < 1e-8, f"exact-network transport error too large: {transport}")
    print(f"PASS exact network: scalar residual={scalar:.6e}, transport sine={transport:.6e}")

    degenerate = partner_subspace_error(-0.375 * np.eye(4), [1, 0, 3, 2], 0.375)
    _require(degenerate < 1e-12, f"degenerate transport error too large: {degenerate}")
    print(f"PASS degenerate invariant-subspace transport: sine={degenerate:.6e}")
    wrong = partner_subspace_error(np.diag([-0.5, -0.5, -0.25, -0.25]), [1, 0, 3, 2], 0.375)
    _require(wrong > 0.9, f"wrong permutation was not rejected: {wrong}")
    print(f"PASS wrong-permutation transport rejection: sine={wrong:.6e}")
    print("ALL NEURAL TRANSLATION GATES PASS")


if __name__ == "__main__":
    main()
