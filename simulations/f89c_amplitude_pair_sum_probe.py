"""F89c amplitude-pair-sum probe: structural anchor at n ↔ k+2−n.

For even k, the S_2-anti Bloch orbit {2, 4, ..., k} closes under the chiral
involution n ↔ k+2−n with y_{k+2−n} = −y_n (since cos(π−θ) = −cos(θ)). The
pair-sum σ_n + σ_{k+2−n} therefore reduces to the EVEN-degree-only part of
P_path(y) evaluated at y_n² — half the polynomial complexity.

Tests at even k where the closed form is tabulated (k=4, 6, 8):
- Pair-sum identity σ_n + σ_{k+2−n} = 2·P_even(y_n²) / [D·N²(N−1)]
- Whether pair-sum value is RATIONAL (Galois invariant under chiral pair)
- Whether denominator structure reveals D_k pattern directly
- Whether the FIXED POINT at n = (k+2)/2 (when present) carries special structure

If σ_n + σ_{k+2−n} has a closed form whose denominator IS D_k (or scales
cleanly with D_k), that would be the F89c-amplitude-layer anchor sought
in PROOF_F89_PATH_D_CLOSED_FORM.md Open Questions.
"""
import sys

import sympy as sp
from sympy import cos, pi, Rational, Symbol, simplify, factor, factorint

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PATH_TABLE = {
    3: ([47, 14], 9),
    4: ([25, 10], 4),
    5: ([129, 82, 13], 25),
    6: ([80, 72, 17], 18),
    7: ([382, 292, 130, 21], 98),
    8: ([110, 68, 54, 13], 32),
    9: ([1476, 440, 288, 190, 31], 324),
}

y = Symbol('y')


def v2(n: int) -> int:
    n = abs(int(n))
    v = 0
    while n > 0 and n % 2 == 0:
        n //= 2
        v += 1
    return v


def odd_part(n: int) -> int:
    n = abs(int(n))
    while n > 0 and n % 2 == 0:
        n //= 2
    return max(n, 1)


def D_k(k: int) -> int:
    vk = v2(k)
    E = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
    return odd_part(k) ** 2 * (2 ** E)


def closed_form_sigma_pair_sum(k: int) -> dict:
    """Symbolically compute σ_n + σ_{k+2−n} for every chiral pair at even k."""
    P_coefs, D = PATH_TABLE[k]
    P = sum(c * y**i for i, c in enumerate(P_coefs))
    chainN = k + 1
    denom = D * chainN ** 2 * (chainN - 1)

    N_block = chainN
    orbit = [2 * i for i in range(1, N_block // 2 + 1)]

    pairs = []
    fixed = None
    seen = set()
    for n in orbit:
        if n in seen:
            continue
        m = (k + 2) - n
        if m in orbit:
            if n < m:
                pairs.append((n, m))
                seen.add(n)
                seen.add(m)
            elif n == m:
                fixed = n
                seen.add(n)

    out = {"k": k, "orbit": orbit, "pairs": [], "fixed": None, "D_k": D, "D_pred": D_k(k)}

    for (n, m) in pairs:
        y_n = 4 * cos(pi * n / (k + 2))
        y_m = 4 * cos(pi * m / (k + 2))
        sigma_n = P.subs(y, y_n) / denom
        sigma_m = P.subs(y, y_m) / denom
        pair_sum = simplify(sigma_n + sigma_m)
        # As a rational: numerator and denominator after simplify
        pair_sum_rational = sp.nsimplify(pair_sum, rational=True)
        out["pairs"].append({
            "n": n,
            "m": m,
            "y_n": sp.nsimplify(sp.simplify(y_n), rational=False),
            "y_m": sp.nsimplify(sp.simplify(y_m), rational=False),
            "pair_sum": pair_sum_rational,
            "numerator": pair_sum_rational.p if hasattr(pair_sum_rational, 'p') else None,
            "denominator": pair_sum_rational.q if hasattr(pair_sum_rational, 'q') else None,
        })

    if fixed is not None:
        y_f = 4 * cos(pi * fixed / (k + 2))
        sigma_f = P.subs(y, y_f) / denom
        sigma_f_rational = sp.nsimplify(simplify(sigma_f), rational=True)
        out["fixed"] = {
            "n": fixed,
            "y": sp.simplify(y_f),
            "sigma": sigma_f_rational,
            "numerator": sigma_f_rational.p if hasattr(sigma_f_rational, 'p') else None,
            "denominator": sigma_f_rational.q if hasattr(sigma_f_rational, 'q') else None,
        }

    return out


def factor_str(n: int) -> str:
    if n is None:
        return "n/a"
    n = abs(int(n))
    if n == 1:
        return "1"
    facs = sorted(factorint(n).items())
    return "·".join(f"{p}^{e}" if e > 1 else str(p) for p, e in facs)


print("# F89c amplitude pair-sum probe: σ_n + σ_{k+2−n} for even k")
print()

for k in [4, 6, 8]:
    result = closed_form_sigma_pair_sum(k)
    print(f"# k = {k}, orbit = {result['orbit']}, D_k = {result['D_k']}")
    for pair in result["pairs"]:
        ps = pair["pair_sum"]
        num = pair["numerator"]
        den = pair["denominator"]
        print(f"#   pair (n={pair['n']}, m={pair['m']}): σ_n + σ_m = {ps}")
        print(f"#     numerator   = {num} = {factor_str(num)}")
        print(f"#     denominator = {den} = {factor_str(den)}")
        print(f"#     denom / D_k = {Rational(den, result['D_k']) if den else 'n/a'}")
    if result["fixed"]:
        f = result["fixed"]
        print(f"#   FIXED n={f['n']}, y=0: σ_fix = {f['sigma']}")
        print(f"#     numerator   = {f['numerator']} = {factor_str(f['numerator'])}")
        print(f"#     denominator = {f['denominator']} = {factor_str(f['denominator'])}")
        print(f"#     denom / D_k = {Rational(f['denominator'], result['D_k']) if f['denominator'] else 'n/a'}")
    print()

# Cross-table: pair-sum denominators vs D_k
print("# Cross-table: denominators of σ pair-sums vs D_k 2-adic structure")
print(f"# {'k':>3} {'D_k':>10} {'D_k factored':<18} {'pair-sum denoms':<40}")
for k in [4, 6, 8]:
    result = closed_form_sigma_pair_sum(k)
    denoms = [p["denominator"] for p in result["pairs"]]
    if result["fixed"]:
        denoms.append(result["fixed"]["denominator"])
    denom_factors = [factor_str(d) for d in denoms]
    print(f"  {k:>3} {result['D_k']:>10} {factor_str(result['D_k']):<18} {' | '.join(denom_factors):<40}")
