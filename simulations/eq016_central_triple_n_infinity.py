#!/usr/bin/env python3
"""EQ-016: analytical N→∞ limit of central-Dicke-triple max pair-CPsi.

In the N→∞ limit at half-filling k = N/2, all matrix elements C(N-2, n-s)/C(N,n)
of the |D_n>'s pair-AB diagonal converge to 1/4. The off-diagonal between
adjacent Dicke states converges to specific simple forms.

Result: at large N with |psi> = a|D_{k-1}> + b|D_k> + a|D_{k+1}>, the limit
ρ_AB matrix is:

    ρ_∞ = [1/4,    ab/2,   ab/2,   a²/4]
          [ab/2,   1/4,    1/4,    ab/2]
          [ab/2,   1/4,    1/4,    ab/2]
          [a²/4,   ab/2,   ab/2,   1/4]

(symmetric structure from k = N/2 — particle-hole D_{k-1} ↔ D_{k+1}).

Tr(ρ_∞²) = 4·(1/4)² + 4·(ab/2)² + 2·(a²/4)² + 4·(ab/2)² + 2·(1/4)²
        = 3/8 + 2a²b² + a⁴/8

L1_off_∞ = 4·|ab/2| + 2·|a²/4| + 2·|1/4| + 4·|ab/2|
        = 4·ab + a²/2 + 1/2     (positive a, b)

cpsi_∞(a, b) = Tr(ρ_∞²) · L1_off_∞ / 3.

Constraint: 2a² + b² = 1 (norm). Maximize over (a, b).

This script:
  1. Numerically maximize cpsi_∞ to high precision.
  2. Check candidate algebraic forms for the optimum (a, b, cpsi_∞).
  3. Compare to numerical asymptote from finite-N data.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, brentq

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def cpsi_inf(x):
    """Asymptotic cpsi at x = a², b² = 1 - 2x. Param x ∈ [0, 1/2]."""
    if x < 0 or x > 0.5:
        return -1.0
    b2 = 1 - 2 * x
    if b2 < 0:
        return -1.0
    s = math.sqrt(x * b2)  # = ab if a, b > 0
    a4 = x * x
    a2b2 = x * b2
    Tr2 = 3.0 / 8 + 2 * a2b2 + a4 / 8
    L1 = 1.0 / 2 + 4 * s + x / 2
    return Tr2 * L1 / 3.0


def main():
    print("=" * 80)
    print("Central-Dicke-triple: N→∞ asymptotic optimum")
    print("=" * 80)
    print()
    print("Limit ρ_∞: all diagonals 1/4, e_01=e_12=ab/2, e_02=a²/4, d_1_off=1/4.")
    print("Tr(ρ_∞²) = 3/8 + 2a²b² + a⁴/8")
    print("L1_off_∞ = 1/2 + 4ab + a²/2")
    print("cpsi_∞(a,b) = Tr(ρ_∞²) · L1_off_∞ / 3   subject to 2a²+b²=1.")
    print()

    # Parametrise by x = a^2 ∈ [0, 1/2]; b² = 1 - 2x.
    # Brute force: very fine grid + scipy refine.
    xs = np.linspace(0.001, 0.499, 100000)
    cps = np.array([cpsi_inf(float(x)) for x in xs])
    best_idx = int(np.argmax(cps))
    x_opt = float(xs[best_idx])
    cps_opt = float(cps[best_idx])

    # Refine
    res = minimize(lambda v: -cpsi_inf(float(v[0])), [x_opt],
                    method="Nelder-Mead",
                    options={"xatol": 1e-15, "fatol": 1e-18})
    x_opt = float(res.x[0])
    cps_opt = -res.fun

    a_opt = math.sqrt(x_opt)
    b_opt = math.sqrt(1 - 2 * x_opt)

    print(f"Optimum:")
    print(f"  x = a² = {x_opt:.12f}")
    print(f"  a = {a_opt:.12f}")
    print(f"  b = {b_opt:.12f}")
    print(f"  a/b = {a_opt/b_opt:.12f}")
    print(f"  a²/b² = {x_opt/(1 - 2*x_opt):.12f}")
    print(f"  cpsi_∞ = {cps_opt:.12f}")
    print()

    # Compare to candidate algebraic forms
    print("Candidate algebraic forms for x = a²:")
    candidates = [
        ("9/34",   9/34),
        ("4/15",   4/15),
        ("5/19",   5/19),
        ("(3-√3)/4", (3 - math.sqrt(3))/4),
        ("(7-√29)/2", (7 - math.sqrt(29))/2),
    ]
    for name, val in candidates:
        diff = abs(val - x_opt)
        print(f"  {name} = {val:.12f}, diff to x_opt = {diff:.2e}")
    print()

    print("Candidate algebraic forms for cpsi_∞:")
    candidates_cp = [
        ("3/7",      3/7),
        ("0.4313",   0.4313),
        ("271189/628864", 271189/628864),
        ("17/(2·19) = 17/38", 17/38),
        ("(3+√5)/12", (3 + math.sqrt(5))/12),
    ]
    for name, val in candidates_cp:
        diff = abs(val - cps_opt)
        print(f"  {name} = {val:.12f}, diff = {diff:.2e}")
    print()

    # Compare to numerical finite-N data
    print(f"Compare: numerical cpsi at N=10000 was 0.4313396050 (diff = "
          f"{abs(cps_opt - 0.4313396050):.4e})")
    print()

    # The exact optimum x satisfies dcpsi_∞/dx = 0.
    # Set up: A(x) = 3/8 + 2x - 31x²/8, B(x) = 1/2 + 4√(x(1-2x)) + x/2.
    # cpsi = A·B/3. Critical point: A·dB/dx + B·dA/dx = 0.
    # Squaring to eliminate √ gives a polynomial equation.
    print("Setting up exact critical-point equation (squared, polynomial in x):")
    print("  Let P(x) = 3 + 16x - 31x²,  R(x) = 1 + x,  S(x) = 8 - 31x.")
    print("  Critical point: (x - 2x²)·(19 - 30x - 93x²)² = (12 + 144x - 1132x² + 1488x³)²")
    print()

    def lhs_minus_rhs(x):
        p1 = (x - 2 * x ** 2) * (19 - 30 * x - 93 * x ** 2) ** 2
        p2 = (12 + 144 * x - 1132 * x ** 2 + 1488 * x ** 3) ** 2
        return p1 - p2

    # Find roots
    xs_check = np.linspace(0.001, 0.499, 1000)
    poly_vals = np.array([lhs_minus_rhs(x) for x in xs_check])
    sign_changes = []
    for i in range(len(poly_vals) - 1):
        if poly_vals[i] * poly_vals[i + 1] < 0:
            try:
                root = brentq(lhs_minus_rhs, xs_check[i], xs_check[i + 1])
                sign_changes.append(root)
            except Exception:
                pass
    print(f"  Polynomial roots in (0, 0.5): {sign_changes}")
    for r in sign_changes:
        cpsi_at_r = cpsi_inf(r)
        print(f"    x = {r:.12f}, cpsi_∞ = {cpsi_at_r:.12f}")


if __name__ == "__main__":
    main()
