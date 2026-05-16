#!/usr/bin/env python3
"""Attempt: unify the two perspectives on a PTF break via F71-decomposition of ln α.

Hypothesis: ln α (per-site PTF time-rescaling) decomposes into F71-symmetric
and F71-antisymmetric components. The two perspectives on a break are the
two orthogonal projections of one object:

  (ln α)_sym,i = (ln α_i + ln α_{N-1-i}) / 2
  (ln α)_anti,i = (ln α_i - ln α_{N-1-i}) / 2

Then:
  Σ_i ln α_i      = Σ_i (ln α)_sym,i         (antisym sums to 0)     ← Reading 1 (error)
  ‖(ln α)_anti‖²  = Σ_i ((ln α)_anti,i)²                              ← Reading 2 (F71-asymmetry)

The PTF closure violation (Reading 1) is the F71-sym sum.
The F71-anti-palindromy violation (Reading 2) is the F71-anti norm.
ONE object, TWO orthogonal projections, BOTH on the table at once.

This is the same algebraic pattern as F81 (M = M_sym + M_anti via Π-conjugation),
applied one layer up: ln α decomposes under F71-spatial-reflection the way M
decomposes under Π. The "F-chain inheritance" pattern recurring.

Test: run on all 6 cases from _multi_lens_ptf_carrier.py. Check:
  1. Decomposition reconstructs ln α correctly (sanity).
  2. Closure-holding cases (truly XX+YY, XY+YX, XZ+ZX) → both projections tiny.
  3. Closure-breaking, F71-symmetric cases (YZ+ZY, IY+YI) → sym large, anti zero.
  4. Closure-breaking, F71-broken cases (XZ+XZ) → both non-zero.

If the four-quadrant prediction holds, the F71-decomposition is the unified formula.
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

# α values from _multi_lens_ptf_carrier.py output (frozen here):
ALPHAS = {
    'truly XX+YY': np.array([1.0000, 1.0000, 1.0000]),
    'XY+YX':       np.array([0.9955, 0.9907, 0.9955]),
    'IY+YI':       np.array([0.1000, 0.1000, 0.1000]),
    'YZ+ZY':       np.array([0.1721, 0.1275, 0.1721]),
    'XZ+ZX':       np.array([1.0095, 0.9988, 1.0095]),
    'XZ+XZ':       np.array([0.1721, 1.8305, 1.8853]),
}


def f71_decompose(f):
    """F71-spatial-reflection decomposition of a per-site function f.

    f_sym,i  = (f_i + f_{N-1-i}) / 2
    f_anti,i = (f_i - f_{N-1-i}) / 2

    f = f_sym + f_anti  (always, regardless of any F71 property).
    """
    f = np.asarray(f, dtype=float)
    f_reversed = f[::-1]
    f_sym = (f + f_reversed) / 2
    f_anti = (f - f_reversed) / 2
    return f_sym, f_anti


def main():
    print("Unified formula attempt: F71-decomposition of ln α")
    print(f"  α from _multi_lens_ptf_carrier.py (|+-+⟩ N=3, γ=0.1, γ_T1=0.01)")
    print()

    # Header
    cols = [
        ('case',           14),
        ('α',              25),
        ('ln α',           25),
        ('(ln α)_sym',     25),
        ('(ln α)_anti',    25),
        ('Σ_sym',          10),
        ('‖anti‖',         10),
        ('quadrant',       18),
    ]
    print(''.join(f"{name:<{w}s}" for name, w in cols))
    print('-' * sum(w for _, w in cols))

    for case, alphas in ALPHAS.items():
        ln_alpha = np.log(alphas)
        sym, anti = f71_decompose(ln_alpha)
        sum_sym = float(np.sum(sym))
        norm_anti = float(np.linalg.norm(anti))

        # Sanity check: decomposition reconstructs
        recon = sym + anti
        assert np.allclose(recon, ln_alpha, atol=1e-12), \
            f"Decomposition failed: recon={recon}, ln_alpha={ln_alpha}"

        # Quadrant
        sum_break = abs(sum_sym) > 0.05
        anti_break = norm_anti > 0.05
        if not sum_break and not anti_break:
            quadrant = "(both tiny)"
        elif sum_break and not anti_break:
            quadrant = "sym only"
        elif not sum_break and anti_break:
            quadrant = "anti only"
        else:
            quadrant = "BOTH break"

        def fmt_arr(a, w=22, prec=3):
            return f"[{', '.join(f'{x:+.{prec}f}' for x in a)}]".ljust(w)

        row = [
            f"{case:<14s}",
            fmt_arr(alphas, prec=3) + ' ',
            fmt_arr(ln_alpha, prec=3) + ' ',
            fmt_arr(sym, prec=3) + ' ',
            fmt_arr(anti, prec=3) + ' ',
            f"{sum_sym:>+9.4f} ",
            f"{norm_anti:>9.4f} ",
            f"{quadrant:<18s}",
        ]
        print(''.join(row))

    print()
    print("Reading guide:")
    print("  ln α          = per-site log-rescaling factor.")
    print("  (ln α)_sym    = F71-symmetric part (mirror-symmetric across i ↔ N-1-i).")
    print("  (ln α)_anti   = F71-antisymmetric part (anti-mirror).")
    print("  Σ_sym         = PTF closure violation     ← Reading 1 'error'")
    print("  ‖anti‖        = F71-asymmetry indicator   ← Reading 2 'calibration'")
    print()
    print("Quadrant predictions (the unified formula's four cells):")
    print("  (both tiny)   → PTF closure holds; case is in perturbative window")
    print("  sym only      → closure breaks symmetrically; F71-anti-palindromic γ_eff (F91 invariant)")
    print("  anti only     → would mean closure holds but spatial asymmetry; unusual (predict: rare)")
    print("  BOTH break    → closure breaks AND F71-broken; predict spectral REORDERING under F91")
    print()
    print("Verification on the 6-case data:")
    print("  XY+YX  / XZ+ZX                                    → (both tiny)    PTF lens applies")
    print("  truly XX+YY                                       → (both tiny)    baseline (trivial)")
    print("  IY+YI / YZ+ZY                                     → 'sym only'    F71-symmetric breaks")
    print("  XZ+XZ                                              → 'BOTH break'  asymmetric break")
    print()
    print("The unified formula:")
    print()
    print("  ln α_i = (ln α)_sym,i + (ln α)_anti,i")
    print()
    print("  Reading 1 reads (Σ over i) of the FIRST component.")
    print("  Reading 2 reads (norm) of the SECOND component.")
    print("  Both readings are projections of ONE object — ln α — onto F71-orthogonal axes.")
    print()
    print("Pattern recurrence (the inheritance shape repeats):")
    print("  F81 decomposes M = M_sym + M_anti via Π-conjugation (operator level).")
    print("  Here we decompose ln α = (ln α)_sym + (ln α)_anti via F71-reflection (rescaling level).")
    print("  Same algebra (sym+anti via group action), different operand. Two-perspectives = two")
    print("  orthogonal projections of one object — at every level the framework looks.")


if __name__ == "__main__":
    main()
