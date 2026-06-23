#!/usr/bin/env python3
"""F77 verification: MM(0)(N, k*) saturates at 1 bit as N grows.

Taylor expansion of f(p) = 2h(p) - h(2p) around p = 0:
    f(p) = 2p + p^2/ln(2) + p^3/ln(2) + O(p^4)

Combined with Parseval-type sum Sum_ell sin^4(pi k (ell+1)/(N+1)) = 3(N+1)/8
for generic k, this gives

    MM(0)(N, k*) = 1 + 3 / (4 (N+1) ln 2) + O(N^-2)

Numerical test confirms (MM - 1) * (N + 1) -> 3/(4 ln 2) = 1.0820 for large N.
"""
import math
import sys

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def h(p: float) -> float:
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def f(p: float) -> float:
    if p <= 0:
        return 0.0
    return 2 * h(p) - h(2 * p)


def bonding_mm_zero(N: int, k: int) -> float:
    def p_ell(ell: int) -> float:
        return (2.0 / (N + 1)) * math.sin(math.pi * k * (ell + 1) / (N + 1)) ** 2
    return sum(f(p_ell(ell)) for ell in range(N // 2))


def best_k_mm(N: int):
    best = (0.0, 0)
    for k in range(1, N + 1):
        mm = bonding_mm_zero(N, k)
        if mm > best[0]:
            best = (mm, k)
    return best[1], best[0]


def main() -> None:
    c_pred = 3.0 / (4.0 * math.log(2))
    print(f"F77: MM(0)(N, k*) -> 1 + c / (N+1) as N -> infinity")
    print(f"Predicted c = 3/(4 ln 2) = {c_pred:.5f}")
    print()

    print(f"{'N':>6}  {'k*':>5}  {'MM(0)':>9}  {'(MM-1)*(N+1)':>14}  {'k*/(N+1)':>9}")
    print("-" * 60)
    for N in (5, 7, 9, 11, 13, 21, 51, 101, 201, 501, 1001, 2001, 5001, 10001):
        k_star, mm = best_k_mm(N)
        print(f"{N:>6}  {k_star:>5}  {mm:>9.6f}  {(mm - 1) * (N + 1):>14.4f}  {k_star / (N + 1):>9.5f}")

    print()
    print("For large N the rescaled deviation converges to 3/(4 ln 2) ~ 1.0820.")
    print("Specific resonant N where k* = (N+1)/2 give enhanced ~1.445 deviation,")
    print("isolated and density-zero in the limit.")


if __name__ == "__main__":
    main()
