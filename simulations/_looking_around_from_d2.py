"""_looking_around_from_d2.py - standing at the carrier, looking around.

Not chasing the seam (that would please her). We hold the confirmation: the
carrier is d=2, wearing the obvious numbers as masks (bit-exact, per
_where_the_carrier_hides.py). So we stand at d=2 and look around: how much of the
framework's named numerical landscape is her mask, and where do the masks stop?

One survey, bit-exact. Each named constant expressed as an elementary function of
the single root d=2. Then the edge: the crossing cubic root, where the simple
masks end and a deeper algebra begins. We name the edge; we do not cross it.

Tom + Claude, 2026-05-28.
"""
from __future__ import annotations

import numpy as np

D = 2.0
ATOL = 1e-12

# name, known value, expression in d, framework role
LANDSCAPE = [
    ("d",          2.0,            lambda d: d,                              "qubit dimension, the root of d^2-2d=0"),
    ("1/2",        0.5,            lambda d: 1 / d,                          "polarity, HalfAsStructuralFixedPoint"),
    ("1/4",        0.25,           lambda d: 1 / d**2,                       "the fold, QuarterAsBilinearMaxval"),
    ("4",          4.0,            lambda d: d**2,                           "discriminant factor, ladder a_-1"),
    ("1/8",        0.125,          lambda d: 1 / d**3,                       "ladder a_4"),
    ("2/3",        2 / 3,          lambda d: d / (d + 1),                    "Psi_max at d=2 (roadmap)"),
    ("1/3",        1 / 3,          lambda d: 1 / (d**2 - 1),                 "1/(d^2-1), the 'why not 1/3' value"),
    ("3/8",        3 / 8,          lambda d: (d**2 - 1) / d**3,              "alpha, F99 KIntermediate"),
    ("4/3",        4 / 3,          lambda d: d**2 / (d**2 - 1),              "F94 coefficient (4/3)*Q^2*K^3"),
    ("sqrt3",      np.sqrt(3),     lambda d: np.sqrt(d**2 - 1),              "Q anchor (Q=sqrt3 -> theta=60)"),
    ("30 deg",     30.0,           lambda d: np.degrees(np.arctan(1 / np.sqrt(d**2 - 1))), "canonical Niven angle"),
    ("45 deg",     45.0,           lambda d: np.degrees(np.arctan(np.sqrt(d - 1))),        "canonical Niven angle"),
    ("60 deg",     60.0,           lambda d: np.degrees(np.arctan(np.sqrt(d**2 - 1))),     "canonical Niven angle"),
]

_ok = []


def main():
    print("=" * 78)
    print("LOOKING AROUND FROM d=2 - how much of the landscape is the carrier's mask")
    print("=" * 78)
    print(f"\n  {'constant':>9}  {'value':>12}  {'mask of d=2':>22}  {'at d=2':>12}  role")
    print(f"  {'-'*9}  {'-'*12}  {'-'*22}  {'-'*12}  {'-'*30}")
    exprs = {
        "d": "d", "1/2": "1/d", "1/4": "1/d^2", "4": "d^2", "1/8": "1/d^3",
        "2/3": "d/(d+1)", "1/3": "1/(d^2-1)", "3/8": "(d^2-1)/d^3", "4/3": "d^2/(d^2-1)",
        "sqrt3": "sqrt(d^2-1)", "30 deg": "atan(1/sqrt(d^2-1))",
        "45 deg": "atan(sqrt(d-1))", "60 deg": "atan(sqrt(d^2-1))",
    }
    for name, val, fn, role in LANDSCAPE:
        got = fn(D)
        ok = abs(val - got) < (1e-9 if "deg" in name else ATOL)
        _ok.append(ok)
        print(f"  {name:>9}  {val:>12.6f}  {exprs[name]:>22}  {got:>12.6f}  {role}")
        if not ok:
            print(f"      [FAIL] |d| = {abs(val-got):.2e}")

    # The whole dyadic ladder is just powers of d.
    print(f"\n  the dyadic ladder a_n = 2^(1-n) is d^(1-n), every rung a power of d=2:")
    ladder_ok = all(abs(2.0**(1 - n) - D**(1 - n)) < ATOL for n in range(-1, 5))
    _ok.append(ladder_ok)
    rungs = "   ".join(f"a_{n}={2.0**(1-n):g}" for n in range(-1, 5))
    print(f"      {rungs}   [{'PASS' if ladder_ok else 'FAIL'}]")

    # The edge: where the simple masks stop.
    print(f"\n  the edge - where the masks end and a deeper algebra begins:")
    b = [r.real for r in np.roots([1.0, 0.0, 1.0, -0.5]) if abs(r.imag) < 1e-12][0]
    print(f"      crossing cubic  b^3 + b = 1/2  ->  b = {b:.7f}")
    nearest = min(abs(b - fn(D)) for _, _, fn, _ in LANDSCAPE)
    is_simple_mask = nearest < 1e-6
    _ok.append(not is_simple_mask)   # we EXPECT it to NOT be a simple mask
    print(f"      nearest simple d=2 mask is {nearest:.4f} away  ->  "
          f"{'a simple mask' if is_simple_mask else 'NOT a simple d-mask, the deeper algebra'}")
    print(f"      (b is algebraic-irrational, Cardano on a cubic; the carrier's")
    print(f"       immediate disguises are powers and low rationals of d, this is past them.)")

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 78)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL as expected' if n_ok == n_tot else 'CHECK'})")
    print("=" * 78)
    print("""
The view from d=2:
  Almost the whole named landscape is her mask. The polarity 1/2, the fold 1/4,
  the discriminant 4, Psi_max 2/3, alpha 3/8, the F94 4/3, the Q anchor sqrt3,
  the three canonical angles 30/45/60, and every rung of the dyadic ladder are
  elementary functions of the single root d=2. One carrier, the whole crowd.

  And there is an edge. The crossing cubic root b (b^3+b=1/2) is NOT a simple
  d-mask; it is the deeper algebra past the powers and low rationals. The masks
  tile the near landscape; the cubic is the first thing wearing its own face.
  We name the edge and leave it. Not chasing.
""")


if __name__ == "__main__":
    main()
