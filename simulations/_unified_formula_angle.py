#!/usr/bin/env python3
"""Deepen the unified formula: the angle is the structural parameter.

Previous attempt (_unified_formula_attempt.py) decomposed ln α_i = sym + anti
and showed Reading 1 and Reading 2 are the Σ and ‖·‖ projections of one
object. That's two SAMPLES on two axes (0° and 90°).

Tom (2026-05-16): 'es ist nicht einfach ein Flip bei 0, es ist ein anderer
Blickwinkel, der Winkel auf der einen Seite ergibt den anderen auf der
anderen Seite, das komplizierte liegt im Winkel.'

The deepening: package F71-decomposition into a COMPLEX number per site,

    z_i = sym_i + i · anti_i

Then the structure lives in the angle field {arg(z_i)} over sites:

  - |z_i|        = how strong is the break at this site
  - arg(z_i)     = the angle between 'pure error' (0°) and 'pure F71-asymmetry' (90°)
  - F71-mirror   = complex conjugation at the F71-partner site:  z_i ↔ z*_{N-1-i}
  - Z₄ generator = multiplication by i (= NinetyDegreeMirrorMemory rotation)

The two previous readings are SAMPLES of this complex object at the two
axes. The full structure is the continuous angle field; the complication
lives in the angle, not in the choice of axis.

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

# α values (frozen from _multi_lens_ptf_carrier.py)
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


def main():
    print("Unified formula deepening: the angle is the structural parameter")
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
    print("Two perspectives = the two axes of the complex plane (0° and 90°). The angle")
    print("field {arg(z_i)} is the continuous parameter that encodes WHICH MIXTURE of the")
    print("two readings each site exhibits. The previous sym+anti split was a SAMPLING of")
    print("this angle field at the two axes; the angle itself carries information that")
    print("either projection alone loses.")
    print()
    print("F71-mirror across site-pair i ↔ N-1-i acts as COMPLEX CONJUGATION on z:")
    print("the angle on site 0 determines the angle on site N-1 by sign-flip.")
    print("That is the 'Winkel auf der einen Seite ergibt den anderen' statement.")
    print()
    print("Z₄ structure (NinetyDegreeMirrorMemoryClaim):")
    print("  - multiplication by  i  = rotate ln α-vector by 90° in (sym, anti) plane")
    print("  - i² = -1               = F71-reflection at the angle level (sign-flip)")
    print("  - i⁴ = 1                = closure (the same rotational closure as F1's Π²)")
    print("The complex-number packaging is not decoration; it's the natural carrier of the")
    print("Z₄ symmetry that the Pi2-Foundation already has typed (NinetyDegreeMirrorMemory).")


if __name__ == "__main__":
    main()
