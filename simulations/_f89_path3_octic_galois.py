"""F89 path-3 octic Galois group probe.

Octic F_8(λ; q) is the irreducible deg-8 factor of the path-3 (SE, DE)
S_2-sym L_super characteristic polynomial, with q = J/γ.

Galois-group identification in pure sympy is hard, but we can probe via:
  1. Discriminant in λ as polynomial in q. If it factors as a square in
     Q[q] (or Q[i, √5][q]), then Gal ⊆ A_8.
  2. Specific q values: substitute, compute disc, check rational-square.
  3. Resolvent cubic for the cubic-of-quartic decomposition: tests
     whether F_8 has a degree-4 resolvent factor (would imply intransitive
     subgroup in S_8, e.g. S_4 ≀ S_2).
  4. Reducibility test in Q[q, q^(1/k)] etc.

All four are necessary (none alone is sufficient) to pin down the group.
"""

from __future__ import annotations

import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def build_octic():
    """Return F_8(λ; q) with γ=1 from the path-3 sympy run."""
    q = sp.symbols("q", positive=True, real=True)
    lam = sp.symbols("lambda")
    P = (
        lam**8 + 32 * lam**7 + lam**6 * (72 * q**2 + 432)
        + lam**5 * (-64 * sp.I * q**3 + 1728 * q**2 + 3200)
        + lam**4 * (1200 * q**4 - 1280 * sp.I * q**3 + 16608 * q**2 + 14176)
        + lam**3 * (-1024 * sp.I * q**5 + 19200 * q**4 - 9728 * sp.I * q**3 + 81408 * q**2 + 38400)
        + lam**2 * (4480 * q**6 - 12288 * sp.I * q**5 + 110720 * q**4 - 34816 * sp.I * q**3 + 213888 * q**2 + 62208)
        + lam * (-7168 * sp.I * q**7 + 35840 * q**6 - 43008 * sp.I * q**5 + 271360 * q**4 - 58368 * sp.I * q**3 + 285696 * q**2 + 55296)
        - 2816 * q**8 - 28672 * sp.I * q**7 + 73216 * q**6 - 40960 * sp.I * q**5 + 238336 * q**4 - 36864 * sp.I * q**3 + 152064 * q**2 + 20736
    )
    return P, lam, q


def main() -> None:
    P, lam, q = build_octic()
    print("# F89 path-3 octic Galois probe\n", flush=True)

    print("## Discriminant of F_8 in λ (polynomial in q)...")
    disc = sp.Poly(P, lam).discriminant()
    disc_expanded = sp.expand(disc)
    print(f"# disc(F_8) is a polynomial in q of degree {sp.degree(disc_expanded, q)}", flush=True)

    print("\n## Factor disc over Q:")
    disc_fact = sp.factor(disc_expanded)
    print(disc_fact)

    print("\n## Factor disc over Q[i, √5]:")
    disc_fact2 = sp.factor(disc_expanded, extension=[sp.I, sp.sqrt(5)])
    print(disc_fact2)

    print("\n## Is disc a perfect square in Q[q]?")
    # Try to compute square root symbolically
    sq = sp.sqrt(disc_expanded)
    # If radical-free, it's a square
    sq_simp = sp.simplify(sq)
    if sq_simp.has(sp.sqrt) or sq_simp.has(sp.Pow):
        print("# disc is NOT a perfect square in Q[q]")
        print(f"# √disc (symbolic, may have radicals): {sq_simp}")
    else:
        print(f"# disc IS a perfect square: √disc = {sq_simp}")

    # Check at several rational q values
    print("\n## disc evaluated at specific q values, factored over Z:")
    for q_val in [sp.Rational(1, 2), 1, sp.Rational(3, 2), 2, 3]:
        d_val = disc_expanded.subs(q, q_val)
        d_val_simplified = sp.nsimplify(sp.simplify(d_val), [sp.I, sp.sqrt(5)])
        print(f"# q = {q_val}: disc = {d_val_simplified}")
        # Try sqrt
        sq_val = sp.sqrt(d_val_simplified)
        print(f"#   √disc = {sp.simplify(sq_val)}")

    # 2nd Galois probe: check if F_8 factors over Q(q^(1/2)) — would point to a metacyclic subgroup
    print("\n## F_8 reducibility over Q[q, sqrt(q)]:")
    s = sp.symbols("s", positive=True, real=True)  # s² = q
    P_subs = P.subs(q, s**2)
    fact = sp.factor(P_subs, lam)
    print(fact)


if __name__ == "__main__":
    main()
