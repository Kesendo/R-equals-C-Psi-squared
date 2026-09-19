#!/usr/bin/env python3
"""Numerical checks for the F97 period-one cardioid parametrization.

For the selected fixed point z*=exp(i phi)/2 of z -> z^2+c,

    c(phi) = exp(i phi)/2 - exp(2 i phi)/4.

Only the selected root has multiplier modulus one along this curve. The other
root is 1-z* and is not generally marginal: at phi=pi/2 its multiplier has
modulus sqrt(5). Also, c=0 has distinct roots 0 and 1. The checks below use a
numerical tolerance; they are not bit-exact floating-point identities.

F95 is a sibling real-parameter region, not an ancestor of F97. At b=1/2 they
share the single real coordinate c=+1/4, and only c=+1/4 is the cardioid cusp.
The radial value |c| is not constant around the cardioid.
"""

from __future__ import annotations

import sys

import numpy as np


def period_one_cardioid_point(phi: float) -> tuple[complex, complex]:
    """Return the selected marginal fixed point and its period-one c value."""
    z_star = np.exp(1j * phi) / 2
    c = z_star - z_star * z_star
    return z_star, c


def fixed_points_at(c: complex) -> tuple[complex, complex]:
    """Return both roots of z^2-z+c=0 in minus/plus square-root order."""
    root = np.sqrt(1 - 4 * c + 0j)
    return ((1 - root) / 2, (1 + root) / 2)


def fixed_point_at(c: complex, target: complex) -> complex:
    """Select the algebraic root nearest a supplied continuation target."""
    return min(fixed_points_at(c), key=lambda root: abs(root - target))


def wrapped_angle(z: complex) -> float:
    """Return arg(z) in [0,2*pi)."""
    angle = np.angle(z)
    return (angle + 2 * np.pi) % (2 * np.pi)


def dense_cardioid_control(sample_count: int = 1000):
    """Reconstruct the selected fixed-point identities on a dense phi grid."""
    if sample_count < 2:
        raise ValueError("sample_count must be at least 2")

    phis = np.linspace(0.0, 2 * np.pi, sample_count)
    identity_residuals = []
    fixed_point_magnitudes = []
    for phi in phis:
        z_star, c = period_one_cardioid_point(phi)
        identity_residuals.append(abs(c - z_star * (1 - z_star)))
        fixed_point_magnitudes.append(abs(z_star))

    evaluated_count = len(identity_residuals)
    max_identity_residual = max(identity_residuals)
    magnitude_spread = max(fixed_point_magnitudes) - min(fixed_point_magnitudes)
    magnitude_std = float(np.std(fixed_point_magnitudes))
    return (evaluated_count, max_identity_residual, magnitude_spread, magnitude_std)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    tolerance = 1e-12
    print("=== F97 PERIOD-ONE CARDIOID: SELECTED AND OTHER ROOTS ===")
    print("The selected fixed point has unit multiplier modulus.")
    print("The other-root multiplier is generally different; at pi/2 it is sqrt(5).")
    print("F95 sibling contact occurs only at the real cusp c=+1/4.")

    phis = (0.0, np.pi / 6, np.pi / 2, np.pi, 2 * np.pi - 1e-12)
    magnitudes = []
    for phi in phis:
        z_star, c = period_one_cardioid_point(phi)
        selected = fixed_point_at(c, z_star)
        other = 1 - z_star
        residual = selected * selected - selected + c
        assert abs(residual) < tolerance
        assert abs(abs(2 * selected) - 1) < tolerance
        magnitudes.append(abs(c))
        print(
            f"phi={phi:.12f}: c={c.real:+.12f}{c.imag:+.12f}i; "
            f"selected |2z|={abs(2 * selected):.12f}; "
            f"other |2z|={abs(2 * other):.12f}"
        )

    z_star, c = period_one_cardioid_point(np.pi / 2)
    other = 1 - z_star
    assert abs(abs(2 * other) - np.sqrt(5)) < tolerance
    print(f"pi/2 other-root multiplier: {abs(2 * other):.12f} = sqrt(5)")

    zero_roots = fixed_points_at(0j)
    assert zero_roots == (0j, 1 + 0j)
    print("c=0 has distinct roots 0 and 1.")

    cusp_z, cusp_c = period_one_cardioid_point(0.0)
    assert abs(cusp_c - 0.25) < tolerance
    assert abs(cusp_z - 0.5) < tolerance
    assert abs((max(magnitudes) - min(magnitudes)) - 0.5) < tolerance
    print("Only c=+1/4 is the cardioid cusp; |c| varies along the curve.")

    dense_sample_count = 1000
    (dense_count, max_identity_residual,
     magnitude_spread, magnitude_std) = dense_cardioid_control(
         dense_sample_count
     )
    assert dense_count == dense_sample_count == 1000
    assert max_identity_residual < tolerance
    assert magnitude_spread < tolerance
    assert magnitude_std < tolerance
    print(
        f"1000-point dense control: count={dense_count}; "
        f"max identity residual={max_identity_residual:.2e}; "
        f"|z*| spread={magnitude_spread:.2e}; std={magnitude_std:.1e}"
    )
    print("All displayed identities were checked with numerical tolerance 1e-12.")


if __name__ == "__main__":
    main()
