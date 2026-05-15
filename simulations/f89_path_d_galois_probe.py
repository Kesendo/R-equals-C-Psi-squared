"""F89 path-D structural probe: Galois-theoretic angle on D_k denominator.

Tests Attack Path 2 from PROOF_F89_PATH_D_CLOSED_FORM.md "Candidate Attack Paths":
"Cyclotomic Galois ring-of-integers: compute the index [O_K : Z[2·cos(π/(k+2))]]
for k=3..24 and check whether its 2-adic content matches the deep-2-power
bonus pattern."

Tests:
1. Index [O_K : Z[α]] where α = 2·cos(π/(k+2)). By Washington's Theorem
   (Introduction to Cyclotomic Fields, Z[ζ + ζ⁻¹] = O_{K_+} for ζ primitive
   n-th root), this index is identically 1 for the maximal real cyclotomic
   subfield. We confirm explicitly via disc(p_α) computation.
2. 2-adic content v₂(disc(p_α)) vs v₂(D_k) = E(k). Look for any pattern
   beyond the well-known disc >> D_k from earlier Angle B.
3. Ratio disc(p_α) / D_k^2 — would be the index^2 if D_k were the index;
   tests whether D_k is even in the right ballpark for a Galois-theoretic
   invariant.

Empirical D_k from F89UnifiedFaClosedFormClaim:
  D_k = odd(k)² · 2^E(k), E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k)-2)
"""
import sys

import sympy as sp
from sympy import Symbol, cos, pi, minimal_polynomial, discriminant, factorint

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def v2(n: int) -> int:
    n = abs(n)
    v = 0
    while n > 0 and n % 2 == 0:
        n //= 2
        v += 1
    return v


def odd_part(n: int) -> int:
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    return n


def D_k(k: int) -> int:
    vk = v2(k)
    E = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
    return odd_part(k) ** 2 * (2 ** E)


def E_k(k: int) -> int:
    vk = v2(k)
    return max(0, (k - 5) // 2) + vk + max(0, vk - 2)


x = Symbol('x')
print("# k=3..14 probe (heavier k cost sympy time; min poly degree grows ~φ(2(k+2))/2)")
print(f"# {'k':>3} {'m=k+2':>5} {'deg':>4} {'D_k':>10} {'E(k)':>5} "
      f"{'v2(D_k)':>7} {'disc(p_α)':>32} {'v2(disc)':>9} {'v2_diff':>8}")

for k in range(3, 15):
    m = k + 2
    alpha = 2 * cos(pi / m)
    try:
        p_alpha = minimal_polynomial(alpha, x)
        deg = int(sp.degree(p_alpha))
        disc = sp.discriminant(p_alpha, x)
        disc_int = int(abs(disc))

        D = int(D_k(k))
        E = int(E_k(k))
        v_D = v2(D)
        v_disc = v2(disc_int)

        print(f"  {k:>3} {m:>5} {deg:>4} {D:>10} {E:>5} {v_D:>7} "
              f"{disc_int:>32} {v_disc:>9} {v_disc - v_D:>+8}")
    except Exception as e:
        print(f"  k={k} m={m}: error {type(e).__name__}: {e}")

print()
print("# Interpretation:")
print("# - 'v2_diff' = v2(disc(p_α)) − v2(D_k)")
print("# - If v2_diff is constant or follows a clean pattern, the 2-adic content")
print("#   of the discriminant is structurally related to E(k).")
print("# - Washington's theorem: index [O_K : Z[α]] = 1 for α = 2·cos(2π/n) on")
print("#   K_+ = max real subfield. So disc(Z[α]) = disc(K), and any relation")
print("#   between disc and D_k is a relation between the field discriminant")
print("#   and the amplitude-layer denominator directly.")

# ---- Factorisation of disc/D_k² ratio (the index² candidate) ----
print()
print("# Ratio disc(p_α) / D_k^2 factored (the candidate 'index²' if D_k were the index):")
print(f"# {'k':>3} {'D_k^2':>14} {'ratio':>32} {'ratio factored':<40}")
for k in range(3, 15):
    m = k + 2
    alpha = 2 * cos(pi / m)
    try:
        p_alpha = minimal_polynomial(alpha, x)
        disc = int(abs(sp.discriminant(p_alpha, x)))
        D = D_k(k)
        D_sq = D * D
        ratio = sp.Rational(disc, D_sq)
        ratio_int = ratio if ratio.is_integer else None
        if ratio_int is not None:
            factored = sorted(factorint(int(ratio)).items())
            factored_str = "·".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factored) if factored else "1"
            print(f"  {k:>3} {D_sq:>14} {str(ratio):>32} {factored_str:<40}")
        else:
            print(f"  {k:>3} {D_sq:>14} {str(ratio):>32} (non-integer)")
    except Exception as e:
        print(f"  k={k}: error {type(e).__name__}: {e}")
