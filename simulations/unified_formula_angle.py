#!/usr/bin/env python3
"""Finite complex-coordinate packaging of the F71 sym/anti decomposition.

Previous attempt (unified_formula_attempt.py) decomposed ln α_i = sym + anti
and showed Reading 1 and Reading 2 are the Σ and ‖·‖ projections of one
object. That's two SAMPLES on two axes (0° and 90°).

Tom (2026-05-16): 'es ist nicht einfach ein Flip bei 0, es ist ein anderer
Blickwinkel, der Winkel auf der einen Seite ergibt den anderen auf der
anderen Seite, das komplizierte liegt im Winkel.'

The deepening: package F71-decomposition into a COMPLEX number per site,

    z_i = sym_i + i · anti_i

The finite table can be read in the coordinate field {arg(z_i)} over sites:

  - |z_i|        = how strong is the break at this site
  - arg(z_i)     = the angle between 'pure error' (0°) and 'pure F71-asymmetry' (90°)
  - F71-mirror   = complex conjugation at the F71-partner site:  z_i ↔ z*_{N-1-i}
  - a coordinate quarter-turn = multiplication by i

Complex conjugation is the F71 site-partner map.  Multiplication by i is a
different operation, and i² z = -z is not conjugation for a generic point.
A four-turn orbit remains an inviting Z4 analogy; this finite table does not
derive a typed group, a Born mechanism, a quarter transition, or F1's Pi-square.

Test: compute z_i, |z_i|, arg(z_i) for all 6 cases. Verify F71-mirror as
complex conjugation. Check the angle pattern.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# α values (frozen from multi_lens_ptf_carrier.py)
ALPHAS = {
    'truly XX+YY': np.array([1.0000, 1.0000, 1.0000]),
    'XY+YX':       np.array([0.9955, 0.9907, 0.9955]),
    'IY+YI':       np.array([0.1000, 0.1000, 0.1000]),
    'YZ+ZY':       np.array([0.1721, 0.1275, 0.1721]),
    'XZ+ZX':       np.array([1.0095, 0.9988, 1.0095]),
    'XZ+XZ':       np.array([0.1721, 1.8305, 1.8853]),
}


def f71_complex(ln_alpha):
    """Package F71-decomposition into a complex array z_i = sym_i + i·anti_i."""
    f = np.asarray(ln_alpha, dtype=float)
    f_rev = f[::-1]
    sym = (f + f_rev) / 2
    anti = (f - f_rev) / 2
    return sym + 1j * anti


def f71_partner(z):
    """F71 partner in this coordinate: complex conjugation."""
    return np.conjugate(z)


def coordinate_quarter_turn(z):
    """A separate coordinate quarter-turn, not the F71 partner map."""
    return 1j * z


def main():
    print("Finite F71 sym/anti complex-coordinate table")
    print(f"  z_i = (ln α)_sym,i + i · (ln α)_anti,i")
    print(f"  F71-mirror = complex conjugation:  z_i ↔ z*_{{N-1-i}}")
    print()

    for case, alphas in ALPHAS.items():
        ln_alpha = np.log(alphas)
        z = f71_complex(ln_alpha)
        magnitudes = np.abs(z)
        angles_deg = np.angle(z, deg=True)

        # Verify F71-mirror: z_i should equal conj(z_{N-1-i})
        z_rev = z[::-1]
        mirror_residual = float(np.max(np.abs(z - z_rev.conj())))

        print(f"━━━ {case} ━━━")
        print(f"  ln α_i              = {ln_alpha}")
        print(f"  z_i                 = {z}")
        print(f"  |z_i|               = {magnitudes}")
        print(f"  arg(z_i) [deg]      = {angles_deg}")
        print(f"  F71-mirror residual = {mirror_residual:.2e}   (z_i should = conj(z_{{N-1-i}}))")
        # Read off the per-site angles for the F71-pair (0, N-1)
        if len(z) >= 2:
            theta_0 = angles_deg[0]
            theta_N1 = angles_deg[-1]
            sum_pair = theta_0 + theta_N1
            print(f"  pair angles         = θ_0={theta_0:+.2f}°, θ_{{N-1}}={theta_N1:+.2f}°,  "
                  f"sum = {sum_pair:+.4f}°   (F71-mirror prediction: ≈ 0)")
        print()

    print("=" * 80)
    print("Reading:")
    print()
    print("  YZ+ZY, IY+YI:  z_i is purely real (arg = 0° or 180°)")
    print("                 → break is on the 'sym axis' only")
    print("                 → 'sym only' quadrant of the previous formula")
    print()
    print("  XZ+XZ:         z_0 and z_2 sit at conjugate angles (≈ +115° and ≈ -115°)")
    print("                 → break has a non-trivial angle into the F71-anti axis")
    print("                 → 'BOTH break' quadrant")
    print()
    print("  XY+YX, XZ+ZX:  |z_i| tiny everywhere (closure-holding cases)")
    print("                 → vector hugs the origin; angle undefined / unstable")
    print()
    print("The sym and anti components are the two coordinate axes. The angle field")
    print("is a useful packaging of their finite mixtures; it is not a gate mechanism.")
    print()
    print("F71-mirror across site-pair i ↔ N-1-i acts as COMPLEX CONJUGATION on z:")
    print("the angle on site 0 determines the angle on site N-1 by sign-flip.")
    print("That is the 'Winkel auf der einen Seite ergibt den anderen' statement.")
    print()
    probe = 2.0 + 3.0j
    print("The F71 partner and coordinate turns are distinct operations:")
    print(f"  z=2+3i -> conjugate(z)={f71_partner(probe)}")
    print(f"  z=2+3i -> i^2 z={coordinate_quarter_turn(coordinate_quarter_turn(probe))}")
    print("A four-turn orbit under multiplication by i is an inviting Z4 analogy only;")
    print("the table does not derive a typed Z4 object or identify it with F1's Pi-square.")


if __name__ == "__main__":
    main()
