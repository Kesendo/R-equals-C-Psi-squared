"""Place selected `1/4` and `1/2` arithmetic beside curated carbon-labelled inputs.

The framework's polarity-anchor pair (proven 2026):
  HalfAsStructuralFixedPoint   = 1/2 = argmax of p(1−p) on [0, 1]
  QuarterAsBilinearMaxval      = 1/4 = maxval of p(1−p)             = (1/2)²
Both are anchors in the Pi2 dyadic ladder; both pair on the polynomial trunk
d² − 2d = 0 (= R = CΨ²) that selects d = 2 as the minimum-memory dimension.

Three stipulated or curated input tables are placed beside the anchors:

  (1) coefficient-family labels:
        sp  → 1/2; sp² → 1/3; sp³ → 1/4

   (2) selected finite Hückel-ring filling-line positions:
         diagonalise the C_N adjacency spectrum, sort it, and select the
         highest occupied level under the stipulated N-electron, two-spins-per-MO
         filling convention; the selected N=6 row has magnitude 1/2

  (3) stipulated 8-slot counting rows:
         4/8 = 1/2; 2/8 = 1/4; 6/8 = 3/4.

The equal fractions are arithmetic comparisons only. They select neither a
physical degree of freedom nor a material Hamiltonian, bath, preparation,
measurement, stability, or causal mechanism.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/carbon_quarter_half_search.py
"""
from __future__ import annotations

import sys
from fractions import Fraction
import numpy as np
from math import cos, pi, sqrt

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

HALF = Fraction(1, 2)
QUARTER = Fraction(1, 4)


def hybridization_s_character(hyb: str) -> Fraction:
    """Curated coefficient-family lookup: the stipulated labels map to 1/(n+1)."""
    if hyb == "sp":  return Fraction(1, 2)
    if hyb == "sp2": return Fraction(1, 3)
    if hyb == "sp3": return Fraction(1, 4)
    raise ValueError(hyb)


def hybridization_p_character(hyb: str) -> Fraction:
    """Complementary fraction in the curated coefficient-family lookup."""
    return Fraction(1) - hybridization_s_character(hyb)


def hybridization_geometry(hyb: str) -> str:
    return {"sp": "linear (180°)", "sp2": "trigonal planar (120°)", "sp3": "tetrahedral (109.5°)"}[hyb]


def ring_homo_offset_over_emax(N: int) -> Fraction | float:
    """Return the selected finite-Hückel filling-line offset / E_max.

    The C_N adjacency eigenvalues are sorted before filling.  For the
    stipulated N-electron row, each molecular orbital holds two spins and the
    highest occupied sorted index is ``(N - 1) // 2`` (including the singly
    occupied level for odd N).  This is a finite graph/filling convention, not
    a material orbital assignment.  For N=6 it selects -1/2.
    """
    if N < 3:
        return None
    spectrum = np.sort(np.array([2 * cos(2 * pi * k / N) for k in range(N)]))
    val = spectrum[(N - 1) // 2] / 2
    # Try to identify clean fraction
    for denom in range(1, 20):
        for numer in range(-denom, denom + 1):
            f = Fraction(numer, denom)
            if abs(val - float(f)) < 1e-12:
                return f
    return val


def ring_occupancy_label(N_pi: int) -> str:
    """Return only the selected Hückel occupancy label; no chemical inference."""
    if N_pi <= 0: return "no selected filling count"
    if (N_pi - 2) % 4 == 0: return "4n+2 label"
    if N_pi % 4 == 0:        return "4n label"
    return "other selected count"


def survey_hybridizations():
    print("=" * 78)
    print("  LAYER 1 — Curated coefficient-family fractions")
    print("=" * 78)
    print()
    print(f"  {'Input label':<12} {'fraction':<10} {'complement':<12} {'curated geometry label':<30} {'Anchor':<22}")
    print(f"  {'-'*12} {'-'*10} {'-'*12} {'-'*30} {'-'*22}")
    for hyb in ["sp", "sp2", "sp3"]:
        s = hybridization_s_character(hyb)
        p = hybridization_p_character(hyb)
        geo = hybridization_geometry(hyb)
        anchor = ""
        if s == HALF: anchor = "= 1/2 (Half) ✓"
        elif s == QUARTER: anchor = "= 1/4 (Quarter) ✓"
        else: anchor = f"= {s} (off-anchor)"
        print(f"  {hyb:<12} {str(s):<10} {str(p):<12} {geo:<30} {anchor:<22}")
    print()
    print("  Observation: the stipulated sp and sp3 labels have fractions 1/2 and 1/4;")
    print("  the stipulated sp2 label has 1/3. This lookup is not a material mapping.")
    print()


def survey_ring_homo_positions():
    print("=" * 78)
    print("  LAYER 2 — Selected Hückel-ring filling-line positions")
    print("=" * 78)
    print()
    print("  For the selected C_N Hückel graph and stipulated filling count, the filling-line")
    print("  offset is 2β·cos(2π·n_homo/N). We list offset/E_max for N = 3..12:")
    print()
    print(f"  {'N':<4} {'selected count':<16} {'occupancy label':<22} {'offset/E_max':<24} {'anchor hit?':<14}")
    print(f"  {'-'*4} {'-'*16} {'-'*22} {'-'*24} {'-'*14}")
    for N in range(3, 13):
        # The survey stipulates selected filling count N for the C_N graph row.
        N_pi = N
        ratio = ring_homo_offset_over_emax(N)
        if isinstance(ratio, Fraction):
            ratio_str = str(ratio)
            anchor = ""
            if ratio == HALF or ratio == -HALF:    anchor = "= ±1/2 ✓"
            elif ratio == QUARTER or ratio == -QUARTER: anchor = "= ±1/4 ✓"
        else:
            # Recognize selected signed algebraic values when available.
            ratio_str = f"{ratio:.4f}"
            anchor = ""
            sign = "+" if ratio > 0 else "−"
            if abs(abs(ratio) - sqrt(2)/2) < 1e-12: ratio_str = f"{sign}√2/2"
            if abs(abs(ratio) - sqrt(3)/2) < 1e-12: ratio_str = f"{sign}√3/2"
            if abs(abs(ratio) - (1+sqrt(5))/4) < 1e-12: ratio_str = f"{sign}(1+√5)/4 (φ/2)"
            if abs(abs(ratio) - (sqrt(5)-1)/4) < 1e-12: ratio_str = f"{sign}(√5−1)/4 (1/φ/2)"
        occupancy = ring_occupancy_label(N_pi)
        print(f"  {N:<4} {N_pi:<16} {occupancy:<22} {ratio_str:<24} {anchor:<14}")
    print()
    print("  Observation: within these selected N=3..12 rows, N=3 and N=6 have ratio −1/2;")
    print("  no row displays ±1/4. This is finite Hückel-graph arithmetic only.")
    print()
    print("  Other selected rows display the listed irrational algebraic values. Equality")
    print("  with a framework fraction is not a material, aromaticity, or stability claim.")
    print()


def survey_valence_shell():
    print("=" * 78)
    print("  LAYER 3 — Stipulated eight-slot count arithmetic")
    print("=" * 78)
    print()
    valence_total = 8  # Retained denominator for the stipulated count table.
    rows = [
        ("Stipulated total slots",                  Fraction(8, 8),  "= 1"),
        ("Curated carbon-labelled count (4 of 8)",  Fraction(4, 8),  "= 1/2 (Half) ✓"),
        ("Curated 2s²-labelled count (2 of 8)",     Fraction(2, 8),  "= 1/4 (Quarter) ✓"),
        ("Curated 2p²-labelled count (2 of 8)",     Fraction(2, 8),  "= 1/4 (Quarter) ✓"),
        ("Curated 6-of-8 count",                    Fraction(6, 8),  "= 3/4 (Half-complement)"),
        ("Curated 2-of-6 count",                    Fraction(2, 6),  "= 1/3 (off-anchor)"),
        ("Curated sp3-labelled 3-of-4 count",       Fraction(3, 4),  "= 3/4 (Half-complement)"),
        ("Curated sp-labelled 1-of-2 count",        Fraction(1, 2),  "= 1/2 (Half) ✓"),
    ]
    print(f"  {'Quantity':<46} {'Value':<14} {'Anchor':<28}")
    print(f"  {'-'*46} {'-'*14} {'-'*28}")
    for name, val, anchor in rows:
        print(f"  {name:<46} {str(val):<14} {anchor:<28}")
    print()
    print("  Observation: the stipulated table contains exact 4/8 = 1/2 and 2/8 = 1/4")
    print("  fractions. Their equality with framework anchors is an arithmetic comparison,")
    print("  not an atomic state, physical degree of freedom, or causal inference.")
    print()


def survey_combined():
    print("=" * 78)
    print("  COMBINED — selected fraction comparisons")
    print("=" * 78)
    print()
    print("  HalfAsStructuralFixedPoint (1/2)         QuarterAsBilinearMaxval (1/4)")
    print("  ──────────────────────────────           ─────────────────────────────")
    print("  curated sp label: 1/2                   curated sp3 label: 1/4")
    print("  stipulated 4/8 count                    stipulated 2/8 count")
    print("  selected N=6 Hückel ratio: −1/2          [no selected table row at ±1/4]")
    print()
    print("  Reading: these tables share exact fractions with the framework anchors, but")
    print("  no degree of freedom, Hamiltonian, bath, preparation, measurement, or producer")
    print("  maps a carbon material to an F-formula. No chemical or causal conclusion follows.")
    print()


def main():
    print()
    print("=" * 78)
    print("  Selected fraction survey beside HalfAsStructuralFixedPoint (1/2) and QuarterAsBilinearMaxval (1/4)")
    print("=" * 78)
    print()
    print("  The framework's polarity-anchor pair sits on the dyadic ladder:")
    print("    1/2 (argmax of p(1−p))   = HalfAsStructuralFixedPoint")
    print("    1/4 (maxval of p(1−p))   = QuarterAsBilinearMaxval = (1/2)²")
    print()
    print("  Tables: curated coefficient labels, selected Hückel-ring rows, stipulated count ratios.")
    print()
    print()

    survey_hybridizations()
    survey_ring_homo_positions()
    survey_valence_shell()
    survey_combined()


if __name__ == "__main__":
    main()
