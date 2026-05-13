"""g_eff(N, b) symbolic closed form via F89 AT-locked F_a sigmas (Phase 2, 2026-05-13).

For each path k in {3, 4, 5, 6} (N_block = k+1, N = k+2):
1. sigs[F_a:n](N) = P_path(y_n) / [D_path * N^2 * (N-1)]
2. y_n = 4*cos(pi*n/(N_block+1)), n in S_2-anti orbit {2, 4, ..., 2*floor(N_block/2)}
3. Sum over n: total_F_a(k) = sum_n sigs[F_a:n]

PathPolynomial table from F89UnifiedFaClosedFormClaim.cs:43-62:
  path-3: P=14y+47, D=9
  path-4: P=10y+25, D=4
  path-5: P=13y^2+82y+129, D=25
  path-6: P=17y^2+72y+80, D=18

Output: symbolic + numeric total_F_a per path; compare against Phase 1 empirical
g_eff per sub-class.
"""
import sympy as sp

PATH_POLY_TABLE = {
    3: ([47, 14], 9),       # 14y + 47
    4: ([25, 10], 4),       # 10y + 25
    5: ([129, 82, 13], 25), # 13y^2 + 82y + 129
    6: ([80, 72, 17], 18),  # 17y^2 + 72y + 80
}

# Phase 1 anchors for cross-check
PHASE1_GEFF_MEAN = {
    # Mean per-bond g_eff at N from Phase 1 output (non-escape bonds only).
    # Source: simulations/_f86_hwhm_subclass_stratification.py output, 2026-05-13.
    3: 1.1270,  # N=5: all 4 bonds (no escape sub-class at N=5)
    4: 1.1258,  # N=6: all 5 bonds (no escape sub-class at N=6)
    5: 1.1056,  # N=7: 4 non-escape bonds (skip Orbit2Escape at b=1, b=4)
    6: 1.1109,  # N=8: 4 non-escape bonds (skip Orbit2Escape at b=1, b=5; CentralEscapeOrbit3 at b=3)
}


def s_2_anti_orbit(n_block: int) -> list:
    return [n for n in range(2, n_block + 1, 2)]


def y_n_sym(n_block: int, n: int):
    return 4 * sp.cos(sp.pi * n / (n_block + 1))


def sigma_F_a_sym(k: int, n: int):
    coefs, denom = PATH_POLY_TABLE[k]
    n_block = k + 1
    N = k + 2
    y = sp.Symbol("y")
    poly = sum(c * y**i for i, c in enumerate(coefs))
    return poly.subs(y, y_n_sym(n_block, n)) / (denom * N**2 * (N - 1))


def total_F_a_sym(k: int):
    return sum(sigma_F_a_sym(k, n) for n in s_2_anti_orbit(k + 1))


def main():
    print("=" * 100)
    print("F89 AT-locked F_a sigmas: symbolic + numeric (Phase 2)")
    print("=" * 100)
    print(f"{'k':>3} {'N':>3} {'orbit':<14} {'symbolic total_F_a':<50} {'numeric':>11} {'phase1 g_eff':>13}")
    print("-" * 100)
    for k in [3, 4, 5, 6]:
        N = k + 2
        n_block = k + 1
        orbit = s_2_anti_orbit(n_block)
        total = total_F_a_sym(k)
        total_simplified = sp.nsimplify(sp.simplify(total), rational=False)
        numeric = float(total)
        phase1 = PHASE1_GEFF_MEAN.get(k, float("nan"))
        print(f"{k:>3} {N:>3} {str(orbit):<14} {str(total_simplified):<50} "
              f"{numeric:>11.6f} {phase1:>13.6f}")


if __name__ == "__main__":
    main()
