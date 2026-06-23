#!/usr/bin/env python3
"""F89 Path-D denominator: the bonus-free clean form (Blickwinkel, 2026-06-04).

The proof PROOF_F89_PATH_D_CLOSED_FORM.md writes

    D_k = odd(k)^2 * 2^E(k),   E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)

and calls the last term the "deep-2-power bonus", whose standalone derivation it
flags as open ("2-adic over-divisibility chain / doubling-formula tree", Angle D,
f89_pathk_symbolic_derivation.py lines 513-516). This script shows the bonus is a
PARAMETERISATION ARTIFACT, not a mechanism. Writing D_k with odd(k)^2 as the base
forces the odd/2 split, which over-removes 2^(2 v2) from k^2; the bonus 2^(v2-2) is
exactly the term that adds back to reach k^2. Collapsing it:

    v2(k) + max(0, v2(k)-2) == 2*v2(k) - min(v2(k), 2)        (identity, all v2)

so   odd(k)^2 * 2^(v2 + max(0,v2-2)) = k^2 / 2^min(v2,2)      and the clean form is

    D_k = 2^max(0, floor((k-5)/2)) * k^2 / 2^min(v2(k), 2).

The k^2 is the eigenvector 1/sqrt(k) normalisation squared (in p_n = |S_c|^2 ||Mv||^2 / 2);
2^min(v2,2) is a cap-at-2 cancellation; 2^floor((k-5)/2) is the untouched Chebyshev
degree-growth term. No third "deep" term. Verified bit-exact against the full proof
table k=3..300, and the equivalence proved algebraically for k=3..600.
"""
from __future__ import annotations

import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def v2(n: int) -> int:
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def odd_part(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def poly_deg_pow(k: int) -> int:
    """The Chebyshev degree-growth 2-power exponent, untouched by this reframe."""
    return max(0, (k - 5) // 2)


def D_threeterm(k: int) -> int:
    """The proof's form: odd(k)^2 * 2^E(k) with the three-term E(k)."""
    vk = v2(k)
    E = poly_deg_pow(k) + vk + max(0, vk - 2)
    return odd_part(k) ** 2 * (2 ** E)


def D_clean(k: int) -> int:
    """The bonus-free form: 2^polydeg * k^2 / 2^min(v2,2)."""
    cap = min(v2(k), 2)
    # k^2 is always divisible by 2^cap (since v2(k^2) = 2 v2(k) >= cap), so this is exact.
    return (2 ** poly_deg_pow(k)) * (k * k) // (2 ** cap)


# Full denominator table transcribed from PROOF_F89_PATH_D_CLOSED_FORM.md
# (k=3..24 the 22-point verification table; k=25..300 the pipeline-extended table).
DOC_TABLE = {
    3: 9, 4: 4, 5: 25, 6: 18, 7: 98, 8: 32, 9: 324, 10: 200, 11: 968,
    12: 288, 13: 2704, 14: 1568, 15: 7200, 16: 2048, 17: 18496, 18: 10368,
    19: 46208, 20: 12800, 21: 112896, 22: 61952, 23: 270848, 24: 73728,
    25: 640000, 30: 1843200, 40: 52428800, 46: 1109393408, 50: 5242880000,
    60: 120795955200, 75: 193273528320000, 100: 351843720888320000,
    150: 53126622932283508654080000,
    200: 1584563250285286751870879006720000,
    300: 4014134135735512165476429289076705071076474880000,
}


def main() -> None:
    print("=" * 92)
    print("F89 Path-D denominator: the deep-2-power bonus is a parameterisation artifact")
    print("=" * 92)
    print()

    # 1. Algebraic equivalence: the two forms agree for every k.
    bad = [k for k in range(3, 601) if D_threeterm(k) != D_clean(k)]
    assert not bad, f"three-term vs clean disagree at k={bad[:10]}"
    print(f"  [OK] D_threeterm(k) == D_clean(k) bit-exact for all k = 3..600 ({600 - 3 + 1} values)")

    # 2. The collapsing identity that powers it, term by term in v2.
    bad_id = [vk for vk in range(0, 12) if vk + max(0, vk - 2) != 2 * vk - min(vk, 2)]
    assert not bad_id, f"identity fails at v2={bad_id}"
    print("  [OK] identity  v2 + max(0, v2-2) == 2*v2 - min(v2, 2)  holds for v2 = 0..11")
    print()

    # 3. Both forms vs the proof's own table.
    print("  cross-check against PROOF_F89_PATH_D_CLOSED_FORM.md table:")
    print(f"  {'k':>4} {'v2':>3} {'doc D_k':>34} {'clean form':>34} {'match':>6}")
    for k in sorted(DOC_TABLE):
        dc = D_clean(k)
        ok = (dc == DOC_TABLE[k] == D_threeterm(k))
        assert ok, f"mismatch at k={k}: doc={DOC_TABLE[k]} clean={dc} three={D_threeterm(k)}"
        print(f"  {k:>4} {v2(k):>3} {DOC_TABLE[k]:>34} {dc:>34} {'OK' if ok else 'FAIL':>6}")
    print()
    print(f"  [OK] all {len(DOC_TABLE)} tabulated denominators reproduced by 2^polydeg * k^2 / 2^min(v2,2)")
    print()

    # 4. Show how the spurious "bonus" term arises only in the odd(k)^2 base.
    print("  where the 'bonus' comes from (the odd(k)^2 base over-removes, the bonus adds back):")
    print(f"  {'v2':>3} {'odd-base 2-power':>18} {'= v2 + bonus':>14} {'clean: 2v2 - min(v2,2)':>24}")
    for vk in range(0, 7):
        base = vk + max(0, vk - 2)
        bonus = max(0, vk - 2)
        clean = 2 * vk - min(vk, 2)
        print(f"  {vk:>3} {base:>18} {f'{vk} + {bonus}':>14} {f'2*{vk} - {min(vk,2)} = {clean}':>24}")
    print()
    print("  reading: the three-term E(k) splits k^2 into odd(k)^2 * 2^(2 v2), then the")
    print("  k-self term (v2) and the deep-bonus (max(0,v2-2)) together restore 2v2 - min(v2,2)")
    print("  powers of 2 -- i.e. they rebuild k^2 / 2^min(v2,2). In the clean form k^2 is never")
    print("  split, so neither the k-self term nor the 'deep 2-power bonus' appears at all.")
    print()
    print("  D_k  =  2^max(0, floor((k-5)/2))  *  k^2  /  2^min(v2(k), 2)")
    print()
    print("  Remaining genuine question (sharpened): why exactly 2^min(v2,2)? The k^2 is the")
    print("  eigenvector 1/sqrt(k) normalisation squared; the cancellation of 2-powers from k^2")
    print("  caps at 4. One concrete cap, not a v2-growing 'over-ramification'.")


if __name__ == "__main__":
    main()
