"""F89 path-3: probe whether the degree-8 octic factor splits over Q[√5] or via
substitution (e.g. λ = μ - 4 to centre on rate 4γ).
"""

from __future__ import annotations

import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def main() -> None:
    q = sp.symbols("q", positive=True, real=True)
    lam = sp.symbols("lambda")
    mu = sp.symbols("mu")
    sqrt5 = sp.sqrt(5)

    # Octic from path-3 sympy run
    P = (
        lam**8 + 32 * lam**7 + lam**6 * (72 * q**2 + 432)
        + lam**5 * (-64 * sp.I * q**3 + 1728 * q**2 + 3200)
        + lam**4 * (1200 * q**4 - 1280 * sp.I * q**3 + 16608 * q**2 + 14176)
        + lam**3 * (-1024 * sp.I * q**5 + 19200 * q**4 - 9728 * sp.I * q**3 + 81408 * q**2 + 38400)
        + lam**2 * (4480 * q**6 - 12288 * sp.I * q**5 + 110720 * q**4 - 34816 * sp.I * q**3 + 213888 * q**2 + 62208)
        + lam * (-7168 * sp.I * q**7 + 35840 * q**6 - 43008 * sp.I * q**5 + 271360 * q**4 - 58368 * sp.I * q**3 + 285696 * q**2 + 55296)
        - 2816 * q**8 - 28672 * sp.I * q**7 + 73216 * q**6 - 40960 * sp.I * q**5 + 238336 * q**4 - 36864 * sp.I * q**3 + 152064 * q**2 + 20736
    )

    # Sum of roots = -(coef of lambda^7) = -32, so center = -32/8 = -4
    print("# Octic centered: trace = -32, so 8 roots avg λ = -4 (rate 4γ)")
    print("# Try substitution λ = μ - 4 and see if depressed form has structure.\n")

    P_centered = sp.expand(P.subs(lam, mu - 4))
    print("## Depressed octic in μ = λ + 4:")
    P_centered_collected = sp.collect(P_centered, mu)
    print(P_centered_collected)
    print()

    # Try factor in Q
    print("## Factor over Q (default):")
    f1 = sp.factor(P)
    print(f1)
    print()

    # Try factor over Q[√5]
    print("## Factor over Q[√5]:")
    f2 = sp.factor(P, extension=[sqrt5])
    print(f2)
    print()

    # Try factor over Q[i]
    print("## Factor over Q[i]:")
    f3 = sp.factor(P, extension=[sp.I])
    print(f3)
    print()

    # Try factor over Q[√5, i]
    print("## Factor over Q[i, √5]:")
    f4 = sp.factor(P, extension=[sp.I, sqrt5])
    print(f4)


if __name__ == "__main__":
    main()
