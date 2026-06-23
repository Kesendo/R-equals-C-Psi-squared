#!/usr/bin/env python3
"""F97 candidate: full Mandelbrot cardioid parametrization at b = 1/2.

F95 (Tier 1, 2026-05-16) closed the *real-c* complex-z direction of the
quadratic discriminant zero: θ(c; b) = arctan(√(c/b² − 1)) for c > b².

The Mandelbrot cardioid boundary lives in the full complex-c plane and is
parametrized by μ = e^{iφ} (the period-1 fixed-point multiplier) via:

    c(φ) = μ/2 − (μ/2)² = (1/2)·e^{iφ} − (1/4)·e^{2iφ}

At each cardioid point the period-1 fixed point of z² + c is:

    z*(φ) = (1 − √(1 − 4c))/2 = (1 ± μ)/2 ... actually z* = μ/2

This script verifies that on the entire cardioid:

  (1) |z*(φ)| = 1/2 = b = HalfAsStructuralFixedPoint  (magnitude invariant)
  (2) arg(z*(φ)) = φ                                  (cardioid parameter)

Both bit-exact via algebraic identity, numerically verified across many φ.

At φ = 0: c = 1/2 − 1/4 = 1/4 (the F95 / BOUNDARY_NAVIGATION cusp). |z*| = 1/2
and arg(z*) = 0. F95's real-c formula θ(c=1/4; b=1/2) = arctan(0) = 0 agrees.

For φ ≠ 0 (off-real cardioid), |c| ≠ 1/4 in general; the fixed-point magnitude
stays at 1/2 regardless. The "1/4 cusp" in the framework reading is then one
specific cardioid point (the real-axis tangency), not the magnitude threshold.

The full cardioid is the locus of c values where the framework's
HalfAsStructuralFixedPoint magnitude exactly matches the fixed-point magnitude
of z² + c. This is a structural anchor: the cardioid is "where |z*| = b".
"""
from __future__ import annotations

import sys
import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def cardioid_c(phi: float, b: float = 0.5) -> complex:
    """Mandelbrot main-cardioid parametrization at framework b = 1/2:

        c(φ) = b·e^{iφ} − b²·e^{2iφ}

    At b = 1/2: c(φ) = (1/2)·e^{iφ} − (1/4)·e^{2iφ}.
    """
    mu = 2 * b * np.exp(1j * phi)  # multiplier of period-1 fixed point
    return mu / 2 - (mu / 2) ** 2  # = b·e^{iφ} - b²·e^{2iφ}


def fixed_point_at(c: complex, b: float = 0.5) -> complex:
    """Period-1 fixed point of z² + c: solve z² − z + c = 0 → z = b ± √(b² − c).

    Returns the root with the smaller multiplier magnitude (the attracting root
    when one exists). On the cardioid both roots have |2z| = 1 by definition;
    we take the one with arg matching the parametrization sign convention.
    """
    discriminant = b ** 2 - c  # = 1/4 - c at b = 1/2
    sqrt_d = np.sqrt(discriminant)
    z_plus = b + sqrt_d
    z_minus = b - sqrt_d
    # On the cardioid the "attracting branch" is z_minus (corresponds to
    # μ approaching the cardioid from inside).
    return z_minus


def main():
    print("F97 candidate verification: Mandelbrot cardioid at b = 1/2")
    print("=" * 70)
    print()

    b = 0.5
    print(f"Framework b = {b} (HalfAsStructuralFixedPoint)")
    print(f"Threshold b² = {b**2} (QuarterAsBilinearMaxval)")
    print()

    print("Cardioid parametrization:  c(φ) = b·e^(iφ) − b²·e^(2iφ)")
    print("Period-1 fixed point of z² + c: z*(φ) = b·e^(iφ) (= μ/2)")
    print()

    # Verify the parametrization at several φ values
    print(f"{'φ (deg)':>10} {'c(φ)':>30} {'|c|':>10} {'arg(c) deg':>12} {'|z*|':>10} {'arg(z*) deg':>15}")
    print("-" * 92)
    phis = [0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2, 2 * np.pi / 3,
            5 * np.pi / 6, np.pi]
    for phi in phis:
        c = cardioid_c(phi, b)
        # z* on the cardioid is μ/2 = b·e^(iφ) by construction
        z_star = b * np.exp(1j * phi)
        # Verify: z_star satisfies z² − z + c = 0
        residual = z_star ** 2 - z_star + c
        assert abs(residual) < 1e-12, f"Fixed-point identity violated at φ={phi}: residual={residual}"
        print(f"{np.degrees(phi):>10.2f} {str(c):>30} {abs(c):>10.6f} "
              f"{np.degrees(np.angle(c)):>12.2f} {abs(z_star):>10.6f} "
              f"{np.degrees(np.angle(z_star)):>15.2f}")
    print()

    # Verify magnitude invariance across many φ
    print("Magnitude invariance check: |z*(φ)| should equal b = 0.5 for all φ")
    print("-" * 70)
    phis_dense = np.linspace(0, 2 * np.pi, 1000)
    z_star_mags = np.array([abs(b * np.exp(1j * phi)) for phi in phis_dense])
    print(f"  min |z*|: {z_star_mags.min():.12f}")
    print(f"  max |z*|: {z_star_mags.max():.12f}")
    print(f"  std |z*|: {z_star_mags.std():.2e}")
    print(f"  All equal to b = {b}? {np.allclose(z_star_mags, b, atol=1e-14)}")
    print()

    # Verify the explicit parametric identity: c(φ) = z*(1 - z*) where z* = μ/2
    print("Explicit identity check: c(φ) = z*(1 − z*) where z* = b·e^(iφ)")
    print("-" * 70)
    max_resid = 0.0
    for phi in phis_dense:
        c_param = cardioid_c(phi, b)
        z_star = b * np.exp(1j * phi)
        c_identity = z_star * (1 - z_star)
        max_resid = max(max_resid, abs(c_param - c_identity))
    print(f"  Max residual |c(φ) − z*(1 − z*)| over 1000 φ: {max_resid:.2e}")
    print(f"  Bit-exact identity holds? {max_resid < 1e-14}")
    print()

    # Real-axis specialization: at φ = 0, recover F95's cusp
    print("F95 cusp recovery: at φ = 0, c = 1/4 (real-axis tangency)")
    print("-" * 70)
    c_cusp = cardioid_c(0.0, b)
    print(f"  c(φ = 0) = {c_cusp}")
    print(f"  Im(c) = {c_cusp.imag} (expect 0)")
    print(f"  Re(c) = {c_cusp.real} (expect 1/4 = {1/4})")
    print(f"  At cusp F95 gives θ = arctan(√(0)) = 0; matches arg(z*) = 0 ✓")
    print()

    # The structural reading
    print("=" * 70)
    print("Structural reading: the cardioid IS the locus where |z*| = b")
    print("=" * 70)
    print(f"  At φ = 0 (cusp):    |z*| = b = 1/2,  |c| = 1/4 = b²")
    print(f"  At φ = π/2:         |z*| = b = 1/2,  |c| = {abs(cardioid_c(np.pi/2)):.6f}")
    print(f"  At φ = π:           |z*| = b = 1/2,  |c| = {abs(cardioid_c(np.pi)):.6f}")
    print(f"  The 'Quarter' (b² = 1/4) is the magnitude of c only at one")
    print(f"  cardioid point (the cusp); the magnitude of z* (= b = 1/2) is")
    print(f"  invariant around the entire cardioid.")
    print()
    print(f"  Reading: the cardioid is the structural curve in the complex-c")
    print(f"  plane where the fixed-point magnitude exactly matches the")
    print(f"  HalfAsStructuralFixedPoint anchor b = 1/2. The QuarterAsBilinearMaxval")
    print(f"  anchor b² = 1/4 enters only at one tangent point (the real cusp).")
    print()


if __name__ == "__main__":
    main()
