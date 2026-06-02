#!/usr/bin/env python3
"""Test the closed form for the star compactness upper bound before baking it in.

Claim:  Σ M_λ(star, N) = 4 · [x^{N-1}] (1-x)^{-4} (1-x^2)^{-6}
        (the Littlewood principal specialisation Σ_λ dim S^λ(C^d) x^|λ| =
         (1-x)^{-d} (1-x^2)^{-C(d,2)}, here d=4, times the hub factor 4).

Tested four INDEPENDENT ways (not just against the same combinatorics):
  A  hook-content sum   : Σ_{λ⊢(N-1)} 4·dim S^λ(C^4) via repn_prediction (the
                          rep-theory side, hook-content formula for dim S^λ).
  B  binomial GF sum    : 4·Σ_{m'} C(m'+5,5)·C((N-1-2m')+3,3) (the closed form,
                          coefficient extraction by hand).
  C  sympy series       : 4·[x^{N-1}] of the symbolic rational function (catches
                          any error in the hand convolution B); skipped if no sympy.
  D  character theory   : for S_3 (N=4), the multiplicity m_λ of each irrep in
                          (C^4)^⊗3 via m_λ = (1/|G|)Σ_g χ_λ(g)·4^{cycles(g)},
                          confirming mult = dim S^λ(C^4) (the premise itself) and
                          Σ 4·m_λ = the closed form.

If A == B == C across a wide N range and D confirms the premise, the closed form
(and the 4400 at N=8) is solid: a structured generating-function coefficient, not
a code limit or an artifact.
"""
from __future__ import annotations

import sys
from math import comb

sys.path.insert(0, 'simulations')
from star_degeneracy_repn import repn_prediction, dim_schur  # noqa: E402


def gf_binomial(m, d=4):
    """[x^m] (1-x)^-d (1-x^2)^-C(d,2)  by hand convolution (the closed form)."""
    c2 = comb(d, 2)
    # (1-x)^-d = Σ_k C(k+d-1, d-1) x^k ; (1-x^2)^-c2 = Σ_p C(p+c2-1, c2-1) x^{2p}
    return sum(comb(p + c2 - 1, c2 - 1) * comb((m - 2 * p) + d - 1, d - 1)
               for p in range(m // 2 + 1))


def gf_sympy(m, d=4):
    """[x^m] of the symbolic rational function (independent of gf_binomial)."""
    try:
        from sympy import symbols, series, O
    except Exception:
        return None
    x = symbols('x')
    c2 = comb(d, 2)
    expr = (1 - x) ** (-d) * (1 - x ** 2) ** (-c2)
    s = expr.series(x, 0, m + 1).removeO()
    return int(s.coeff(x, m))


# --- D: S_3 character table (classes: identity, transposition, 3-cycle) ---
# cycle types and class sizes for S_3
_S3_CLASSES = [
    # (class_size, n_cycles, partition-of-the-cycle-type)
    (1, 3),   # identity 1^3 -> 3 cycles
    (3, 2),   # transposition 2.1 -> 2 cycles
    (2, 1),   # 3-cycle -> 1 cycle
]
# irrep characters on (id, transp, 3-cycle), keyed by partition of 3
_S3_CHAR = {
    (3,):       [1, 1, 1],     # trivial
    (2, 1):     [2, 0, -1],    # standard
    (1, 1, 1):  [1, -1, 1],    # sign
}


def s3_multiplicity(lam, d=4):
    """Multiplicity of S_3-irrep lam in (C^d)^⊗3, via m = (1/6) Σ_g χ(g) d^{cyc(g)}."""
    chars = _S3_CHAR[lam]
    total = sum(size * chars[i] * d ** ncyc
                for i, (size, ncyc) in enumerate(_S3_CLASSES))
    assert total % 6 == 0
    return total // 6


def main():
    print("Testing  Σ M_λ(star, N) = 4·[x^(N-1)] (1-x)^-4 (1-x^2)^-6")
    print()
    have_sympy = gf_sympy(2) is not None
    print(f"  sympy available: {have_sympy}")
    print()
    print(f"  {'N':>2s} {'A hook-content':>15s} {'B binomial GF':>14s} "
          f"{'C sympy series':>15s} {'all agree':>10s}")
    print(f"  {'-'*2} {'-'*15} {'-'*14} {'-'*15} {'-'*10}")

    for N in range(3, 17):
        m = N - 1
        A = repn_prediction(N)[0]
        B = 4 * gf_binomial(m)
        C = gf_sympy(m)
        C4 = 4 * C if C is not None else None
        agree = (A == B) and (C4 is None or A == C4)
        assert A == B, f"A != B at N={N}: {A} vs {B}"
        if C4 is not None:
            assert A == C4, f"A != C(sympy) at N={N}: {A} vs {C4}"
        cstr = str(C4) if C4 is not None else "(skipped)"
        print(f"  {N:>2d} {A:>15d} {B:>14d} {cstr:>15s} {'YES' if agree else 'NO':>10s}")

    print()
    print("  D: S_3 character-theory check (N=4), independent of hook-content:")
    print(f"     {'λ':>10s} {'m_λ (characters)':>17s} {'dim S^λ(C^4)':>13s} {'match':>6s}")
    tot = 0
    for lam in [(3,), (2, 1), (1, 1, 1)]:
        m_char = s3_multiplicity(lam)
        m_hook = dim_schur(lam, 4)
        tot += 4 * m_char
        match = 'OK' if m_char == m_hook else 'FAIL'
        assert m_char == m_hook, f"character != hook-content for {lam}: {m_char} vs {m_hook}"
        print(f"     {str(lam):>10s} {m_char:>17d} {m_hook:>13d} {match:>6s}")
    print(f"     Σ 4·m_λ = {tot}  (closed form at N=4 = {4*gf_binomial(3)})")
    assert tot == 4 * gf_binomial(3) == repn_prediction(4)[0] == 176

    print()
    print("  ✓ ALL CHECKS PASS. The closed form holds across N=3..16 by three")
    print("  independent routes (hook-content, binomial GF, sympy series), and the")
    print("  S_3 character computation confirms the multiplicity premise itself.")
    print("  4400 at N=8 = 4·1100 is a structured GF coefficient, not a code limit.")


if __name__ == "__main__":
    main()
