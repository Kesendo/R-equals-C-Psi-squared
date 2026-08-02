#!/usr/bin/env python3
"""How wide is F115's valuation criterion? The census behind the arc `f115_valuation_gate_width`.

ANALYTICAL_FORMULAS F115 statement 1 says a Z-dephasing diagonal-cell Mixed PAIR is hard iff the
two X/Y masks have different (1+x)-adic valuations, with no y-parity condition. The shipped
certifier (PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation) additionally requires the
two terms to agree in Y-count parity, and returns no verdict otherwise. This script counts both
populations so the gap is a number rather than an impression:

    hard MASK-pairs                     = A203241 = (4^(k-1) - 3*2^(k-1) + 2)/3
    hard STRING-pairs by the criterion  = twice the certified family
    of those, y-parity-homogeneous      = 2^(2k-3) * A203241, the registry's dressed hard count

so the registry's count is the certified HALF. Whether the other half is genuinely hard (making
the certifier conservative) or not (making statement 1 too wide) is the open question in the arc;
this script does not answer it, it sizes it.

TRAP, worth the line: the enumeration must run over the FULL alphabet including the identity
letter. Over XYZ alone it gives 6 / 20 / 60 strings and matches nothing, silently.
"""
from __future__ import annotations
import itertools
import sys

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

ALPHABET = "IXYZ"


def is_diagonal_cell_mixed(s: str) -> bool:
    """The CellTerms / TryDiagonalCellMixedMask gate: X/Y count even and >= 2, #(Y/Z) odd."""
    na = sum(1 for c in s if c in "XY")
    nb = sum(1 for c in s if c in "YZ")
    return na % 2 == 0 and na >= 2 and nb % 2 == 1


def xy_mask(s: str) -> tuple[int, ...]:
    return tuple(1 if c in "XY" else 0 for c in s)


def y_parity(s: str) -> int:
    return sum(1 for c in s if c == "Y") % 2


def valuation_at_one_plus_x(mask) -> int:
    """How many times (1+x) divides the GF(2) polynomial with these coefficients."""
    p, v = list(mask), 0
    while True:
        suffix, acc = [], 0
        for c in reversed(p):
            acc ^= c
            suffix.append(acc)
        if acc != 0:                      # p(1) != 0, so (1+x) does not divide
            return v
        p = list(reversed(suffix))[1:]
        v += 1
        if not p:
            return v


def a203241(k: int) -> int:
    return (4 ** (k - 1) - 3 * 2 ** (k - 1) + 2) // 3


def census(k: int):
    strings = [w for w in map("".join, itertools.product(ALPHABET, repeat=k))
               if is_diagonal_cell_mixed(w)]
    masks = sorted({xy_mask(w) for w in strings})
    hard_masks = sum(1 for a, b in itertools.combinations(masks, 2)
                     if valuation_at_one_plus_x(a) != valuation_at_one_plus_x(b))
    hard_pairs = [(a, b) for a, b in itertools.combinations(strings, 2)
                  if valuation_at_one_plus_x(xy_mask(a)) != valuation_at_one_plus_x(xy_mask(b))]
    homogeneous = [(a, b) for a, b in hard_pairs if y_parity(a) == y_parity(b)]
    return strings, masks, hard_masks, hard_pairs, homogeneous


if __name__ == "__main__":
    print(f"{'k':>2} {'strings':>8} {'masks':>6} {'hard masks':>11} {'A203241':>8} "
          f"{'hard pairs':>11} {'y-par homog':>12} {'2^(2k-3)*A':>11}")
    ok = True
    for k in (3, 4, 5):
        strings, masks, hard_masks, hard_pairs, homog = census(k)
        dressed = 2 ** (2 * k - 3) * hard_masks
        ok &= hard_masks == a203241(k)
        ok &= len(homog) == dressed
        ok &= len(hard_pairs) == 2 * len(homog)
        print(f"{k:>2} {len(strings):>8} {len(masks):>6} {hard_masks:>11} {a203241(k):>8} "
              f"{len(hard_pairs):>11} {len(homog):>12} {dressed:>11}")
    print("\nhard masks = A203241, y-parity-homogeneous pairs = 2^(2k-3)*A203241, "
          "and the criterion marks exactly twice the certified family"
          if ok else "\nMISMATCH")
    sys.exit(0 if ok else 1)
