"""Constructible-angle α values and selected-model V-gain arithmetic.

The script evaluates the F86b expression `α = sin²(θ)/2` for a finite list
of constructible angles and labels the resulting values as rational or
algebraic-irrational. It also evaluates the selected-ring-model quantity
`V(N) = 2·cos²(π/(2N))` used in `THERMAL_BREAKING.md`, and an even-ring cosine
comparison.

The shared trigonometric constants are arithmetic observations. This script
does not identify a carbon degree of freedom, a material ring, a thermal
channel, aromaticity, Jahn–Teller behavior, or a chemical mechanism. Such a
bridge remains conditional on a specified model and producer.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/off_niven_angles_as_wave_breaking.py
"""
from __future__ import annotations

import sys
from math import cos, sin, pi, sqrt
from fractions import Fraction

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------- Algebraic identification helper ----------------------

def identify_algebraic(value: float, tol: float = 1e-10) -> str:
    """Try to recognise common algebraic constants in the value."""
    # Rationals up to denominator 32
    f = Fraction(value).limit_denominator(32)
    if abs(float(f) - value) < tol:
        return f"{f} (rational)"

    # Common algebraic constants
    candidates = [
        (sqrt(2) / 2,           "√2/2"),
        (sqrt(3) / 2,           "√3/2"),
        ((sqrt(5) - 1) / 4,     "(√5−1)/4  (=1/(2φ), golden inverse)"),
        ((1 + sqrt(5)) / 4,     "(1+√5)/4  (=φ/2, golden ratio half)"),
        ((sqrt(6) - sqrt(2)) / 4, "(√6−√2)/4 (=sin(15°))"),
        ((sqrt(6) + sqrt(2)) / 4, "(√6+√2)/4 (=cos(15°))"),
        (sqrt((2 - sqrt(2)) / 4), "√(2−√2)/2 (=sin(22.5°))"),
        (sqrt((2 + sqrt(2)) / 4), "√(2+√2)/2 (=cos(22.5°))"),
        ((2 - sqrt(3)) / 4,     "(2−√3)/4"),
        ((3 - sqrt(5)) / 8,     "(3−√5)/8"),
        ((5 - sqrt(5)) / 8,     "(5−√5)/8"),
        ((5 + sqrt(5)) / 8,     "(5+√5)/8"),
        (1 + sqrt(2),           "1 + √2 (silver ratio)"),
        (1 + sqrt(3),           "1 + √3"),
        (2 + sqrt(3),           "2 + √3"),
        (1 + sqrt(5),           "1 + √5"),
    ]
    for cval, label in candidates:
        if abs(cval - value) < tol:
            return f"{value:.10f} = {label}"
    return f"{value:.10f}"


def survey_niven_vs_off_niven():
    """For each canonical angle θ in [0°, 90°] at 7.5° resolution, compute
    sin²(θ) and α = sin²(θ)/2, mark as rational F99-canonical vs
    off-Niven algebraic."""
    print("=" * 86)
    print("  Constructible angles in [0°, 90°]: F99 canonical vs off-Niven algebraic")
    print("=" * 86)
    print()
    print(f"  {'θ':>5} {'sin²(θ)':<38} {'α = sin²(θ)/2':<28} {'Class'}")
    print(f"  {'-'*5} {'-'*38} {'-'*28} {'-'*20}")

    angles_deg = [0, 7.5, 15, 18, 22.5, 30, 36, 45, 54, 60, 67.5, 72, 75, 82.5, 90]
    for deg in angles_deg:
        theta = deg * pi / 180
        sin_sq = sin(theta) ** 2
        alpha = sin_sq / 2

        f_sinsq = Fraction(sin_sq).limit_denominator(16)
        is_niven = abs(float(f_sinsq) - sin_sq) < 1e-12

        if is_niven:
            sin_str = str(f_sinsq)
            alpha_str = str(Fraction(alpha).limit_denominator(16))
            cls = "F99 canonical (rational)"
        else:
            sin_str = identify_algebraic(sin_sq)
            alpha_str = identify_algebraic(alpha)
            cls = "off-Niven (algebraic)"
        print(f"  {deg:>5}° {sin_str:<38} {alpha_str:<28} {cls}")
    print()


def survey_v_effect_angles():
    """The selected-model V-gain V(N) = 2·cos²(π/(2N)) at N=2..8.
    π/(2N) for N=2..8 = 45°, 30°, 22.5°, 18°, 15°, ≈12.86°, ≈11.25°."""
    print("=" * 86)
    print("  Selected-model V-gain angles π/(2N)")
    print("=" * 86)
    print()
    print("  Per THERMAL_BREAKING.md: V(N) = ω_max(N) / ω_max(N=2), on the (0,1) block")
    print("                          = 2·cos²(π/(2N))")
    print()
    print(f"  {'N':>3} {'π/(2N)':>10} {'V(N) = 2cos²(π/(2N))':<32} {'Niven?'}")
    print(f"  {'-'*3} {'-'*10} {'-'*32} {'-'*32}")

    for N in [2, 3, 4, 5, 6, 7, 8]:
        theta_deg = 90 / N
        theta = pi / (2 * N)
        V = 2 * cos(theta) ** 2
        V_str = identify_algebraic(V, tol=1e-10)
        f_V = Fraction(V).limit_denominator(32)
        is_rational = abs(float(f_V) - V) < 1e-12
        niven_label = f"★ Niven (V = {f_V})" if is_rational else "off-Niven (irrational algebraic)"
        print(f"  {N:>3} {theta_deg:>9.4f}° {V_str:<32} {niven_label}")
    print()
    print("  Observation: the displayed N = 2,3 values are rational; the displayed")
    print("  N ≥ 4 values are algebraic-irrational. This is selected-model arithmetic,")
    print("  not a material or thermal-mechanism assignment.")
    print()


def survey_even_ring_cosines():
    """Evaluate |cos(π/N)| for the selected even-ring cosine comparison."""
    print("=" * 86)
    print("  Even-ring cosine values across N")
    print("=" * 86)
    print()
    print("  Selected free-ring cosine comparison: |cos(π/N)|.")
    print()
    print(f"  {'N':>3} {'π/N':>10} {'|cos(π/N)|':<32} {'Class'}")
    print(f"  {'-'*3} {'-'*10} {'-'*32} {'-'*30}")

    for N in [4, 6, 8, 10, 12]:
        theta_deg = 180 / N
        theta = pi / N
        c = abs(cos(theta))
        c_str = identify_algebraic(c, tol=1e-10)
        f_c = Fraction(c).limit_denominator(32)
        is_rational = abs(float(f_c) - c) < 1e-12
        niven_label = "★ Niven (rational)" if is_rational else "off-Niven (algebraic)"
        print(f"  {N:>3} {theta_deg:>9.4f}° {c_str:<32} {niven_label}")
    print()
    print("  These are trigonometric values only; no material-ring interpretation is assigned.")
    print()


def show_summary():
    print("=" * 86)
    print("  SUMMARY — constructible-angle arithmetic")
    print("=" * 86)
    print()
    print("  The finite survey displays two value classes on the α-axis:")
    print()
    print("  (1) F99 canonical angles {0°, 30°, 45°, 60°, 90°} — rational α")
    print("      → named anchors {0, 1/8, 1/4, 3/8, 1/2}")
    print()
    print("  (2) Additional surveyed constructible angles — algebraic-irrational α")
    print("      → √2, √3, and √5 algebraic families in the displayed values")
    print("      → selected-model V(N) values at π/(2N) for the stated model")
    print()
    print("  Both classes arise from the same F86b expression α = sin²(θ)/2 and")
    print("  differ here only by the rationality of the displayed values.")
    print()
    print("  Any bridge from these shared constants to wave-breaking, heat, chemistry,")
    print("  or a material system remains conditional on a separately specified model.")
    print()


def main():
    print()
    print("=" * 86)
    print("  Off-Niven constructible-angle survey")
    print("=" * 86)
    print()
    print("  Exact trigonometric evaluations plus selected-model V-gain arithmetic.")
    print()
    print()
    survey_niven_vs_off_niven()
    survey_v_effect_angles()
    survey_even_ring_cosines()
    show_summary()


if __name__ == "__main__":
    main()
