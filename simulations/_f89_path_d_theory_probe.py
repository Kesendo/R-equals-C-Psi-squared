"""F89 path-D theory probe: structural derivation of D_k denominator (2026-05-13).

Investigates three theoretical angles for the empirical closed form:
  D_k = (odd(k))^2 * 2^E(k)
  E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)

Verified bit-exact at k=3..24 (22 data points) by prior probe scripts.

Theoretical angles probed here:

ANGLE A: Free-fermion structural argument (main result).
  p_n = sigma_n * N^2*(N-1) = |S_c(n)|^2 * ||Mv(n)||^2 / 2
  where:
    S_c(n) = sum of all entries of the normalized F_a eigenvector v_n
    Mv(n) = w @ v_n, w = per-site reduction matrix
             (w[l, b] = 1 if basis element b = (i, jk) with i in jk and other DE site = l)
  The F_a eigenvector at Bloch mode n lives in the overlap subspace.
  S_c(n) and ||Mv(n)||^2 each carry factors from Bloch amplitude sums over (k+2) sites.
  D_k = odd(k)^2 from the fact that |S_c(n)|^2 * ||Mv(n)||^2 has denominator k^2 in Z[y_n]
  after the Bloch-mode normalization.
  The 2-power E(k) comes from the 2J hopping coefficient and degree growth.

ANGLE B: Cyclotomic / ring-of-integers approach.
  Check: does D_k relate to disc(minimal poly of 4*cos(pi/m)) or index [O_K:Z[y]]?
  Conclusion: no clean divisibility; the cyclotomic discriminant is much larger than D_k.

ANGLE C: Vandermonde determinant squared.
  |V_det|^2 for the S_2-anti orbit grows much faster than D_k.
  Algebraic cancellations (p_n has same algebraic structure as y_n) reduce |V_det|^2
  to D_k in the Cramer denominator.

ANGLE D: 2-adic bonus max(0, v2(k)-2).
  This kicks in at k=8 (v2=3, bonus=1), k=16 (v2=4, bonus=2), k=24 (v2=3, bonus=1).
  Structural interpretation: when k=2^a*odd, the Bloch momentum grid has 2^a-fold
  over-divisibility by 2 in the hopping matrix coefficients, creating extra 2-power factors
  in the sigma numerators that grow as 2^{max(0,v2(k)-2)} above the base 2^{v2(k)} level.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def v2(n: int) -> int:
    """2-adic valuation of n."""
    if n <= 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def odd_part(n: int) -> int:
    """Odd part of n."""
    if n <= 0:
        return abs(n)
    while n % 2 == 0:
        n //= 2
    return n


def factorize(n: int) -> dict[int, int]:
    """Trial division prime factorization."""
    factors: dict[int, int] = {}
    d, tmp = 2, abs(n)
    while d * d <= tmp:
        while tmp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            tmp //= d
        d += 1
    if tmp > 1:
        factors[tmp] = factors.get(tmp, 0) + 1
    return factors


def factor_str(n: int) -> str:
    """Compact factorization: '2^3*5^2'."""
    if n == 0:
        return "0"
    n = abs(n)
    if n == 1:
        return "1"
    facs = factorize(n)
    return "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(facs.items()))


def predicted_D(k: int) -> tuple[int, int]:
    """Empirical closed form: D_k = odd(k)^2 * 2^E(k).

    E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2).
    Returns (D_k, E).
    """
    vk = v2(k)
    E = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
    return odd_part(k) ** 2 * (2 ** E), E


# ---------------------------------------------------------------------------
# Ground-truth D values from prior probe scripts (k=3..24, verified bit-exact)
# ---------------------------------------------------------------------------

VERIFIED_D: dict[int, int] = {
    3: 9,
    4: 4,
    5: 25,
    6: 18,
    7: 98,
    8: 32,
    9: 324,
    10: 200,
    11: 968,
    12: 288,
    13: 2704,
    14: 1568,
    15: 7200,
    16: 2048,
    17: 18496,
    18: 10368,
    19: 46208,
    20: 12800,
    21: 112896,
    22: 61952,
    23: 270848,
    24: 73728,
}


# ---------------------------------------------------------------------------
# MAIN TABLE: Full 22-point data set with formula verification
# ---------------------------------------------------------------------------

def print_full_table() -> None:
    print("=" * 100)
    print("FULL 22-POINT VERIFICATION TABLE: D_k = odd(k)^2 * 2^E(k)")
    print("E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)")
    print("=" * 100)
    print()
    hdr = (
        f"{'k':>4}  {'v2(k)':>5}  {'odd(k)':>6}  {'FA':>4}  {'deg':>4}  "
        f"{'E(k)':>5}  {'D_pred':>10}  {'D_verify':>10}  {'match':>6}  factored"
    )
    print(hdr)
    print("-" * len(hdr))
    all_match = True
    for k in range(3, 25):
        n_block = k + 1
        fa = n_block // 2
        deg = fa - 1
        vk = v2(k)
        ok = odd_part(k)
        E = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
        D_pred = ok ** 2 * (2 ** E)
        D_ver = VERIFIED_D.get(k, None)
        match = (D_pred == D_ver) if D_ver is not None else None
        match_str = "OK" if match else ("??" if match is None else "FAIL")
        if match is False:
            all_match = False
        print(
            f"{k:>4}  {vk:>5}  {ok:>6}  {fa:>4}  {deg:>4}  {E:>5}  "
            f"{D_pred:>10}  {D_ver or 0:>10}  {match_str:>6}  {factor_str(D_pred)}"
        )
    print()
    print(f"All formula predictions match verified D values (k=3..24): {all_match}")
    print()


# ---------------------------------------------------------------------------
# ANGLE A: Structural derivation -- why odd(k)^2 is the odd part of D
# ---------------------------------------------------------------------------

def _build_se_de_L_correct(n_block: int, J_val: float, gamma_val: float):
    """Build (SE,DE) sub-block Liouvillian with 2J hopping (matching structure probe)."""
    import numpy as np
    de_pairs = [(j, k) for j in range(n_block) for k in range(j + 1, n_block)]
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    n_basis = len(basis)
    M_SE = np.zeros((n_block, n_block))
    for a in range(n_block):
        for b in range(n_block):
            if abs(a - b) == 1:
                M_SE[a, b] = 2.0 * J_val
    n_de = len(de_pairs)
    M_DE = np.zeros((n_de, n_de))
    for idx, (j, k) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < n_block and new_j != k:
                new_pair = tuple(sorted([new_j, k]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2.0 * J_val
        for new_k in [k - 1, k + 1]:
            if 0 <= new_k < n_block and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - k) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2.0 * J_val
    L = np.zeros((n_basis, n_basis), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE[i2, i] != 0.0:
                L[basis.index((i2, jk)), idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(n_de):
            if M_DE[jk_idx, jk2_idx] != 0.0:
                L[basis.index((i, de_pairs[jk2_idx])), idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2.0 * gamma_val if i in jk else -6.0 * gamma_val
    return L, basis, de_pairs


def angle_a_structural_argument() -> None:
    print("=" * 100)
    print("ANGLE A: Free-fermion structural argument for D_k = odd(k)^2 * 2^E(k)")
    print("=" * 100)
    print()

    print("SETUP:")
    print("  For path-k (k-bond OBC chain with k+1 sites, N_block = k+1, m = k+2):")
    print("  - S_2-anti Bloch orbit: n in {2, 4, ..., 2*floor((k+1)/2)}, FA = floor((k+1)/2) elements")
    print("  - Bloch energies: y_n = 4*J*cos(pi*n/m) = 4*J*cos(pi*n/(k+2))")
    print("  - F_a eigenvalue: lambda_n = -2*gamma + i*y_n (AT-locked at rate 2*gamma)")
    print()
    print("KEY STRUCTURAL RESULT (numerically verified for k=3..6):")
    print("  p_n = sigma_n * N^2*(N-1) = |S_c(n)|^2 * ||Mv(n)||^2 / 2")
    print("  where:")
    print("    S_c(n) = vdot(v_n, ones) = sum of all entries of normalized F_a eigenvector v_n")
    print("    Mv(n) = w @ v_n, with per-site reduction matrix w:")
    print("      w[l, b] = 1 if basis element b = (i, jk) with i in jk and other DE site = l")
    print()
    print("DERIVATION:")
    print("  1. rho_flat in the (SE,DE) block = pre/2 on all dim entries (uniform).")
    print("     pre = sqrt(2 / (N^2*(N-1))) so that rho_flat is normalized.")
    print()
    print("  2. Inner product: c_n = <v_n | rho_flat> = (pre/2) * sum_j conj(v_n[j]).")
    print("     F_a eigenvectors live in the overlap subspace (diagonal -2*gamma).")
    print("     The basis has dim = n_block * C(n_block, 2) elements (i, jk);")
    print("     the F_a eigenvectors span the overlap subspace where i in {j,k}.")
    print("     Summing all entries: c_n = (pre/2) * S_c(n)* (complex conjugate).")
    print("     |c_n|^2 = (pre^2/4) * |S_c(n)|^2.")
    print()
    print("  3. Per-site reduction: (Mv)[l] = sum_{b: other-DE-site=l} v_n[b].")
    print("     F_a modes spread over the full overlap subspace (not concentrated on one site).")
    print("     ||Mv||^2 = sum_l |(Mv)[l]|^2 is the squared norm of the reduced vector.")
    print()
    print("  4. sigma_n = |c_n|^2 * ||Mv||^2 = (pre^2/4) * |S_c(n)|^2 * ||Mv(n)||^2.")
    print()
    print("  5. p_n = sigma_n * N^2*(N-1).")
    print("     pre^2 = 2/(N^2*(N-1)), so p_n = (2/(4*N^2*(N-1))) * |S_c(n)|^2 * ||Mv(n)||^2 * N^2*(N-1)")
    print("          = |S_c(n)|^2 * ||Mv(n)||^2 / 2.")
    print()
    print("  6. For path-3 (k=3), n=2 (y_2 = sqrt(5)-1 = 1.2361):")
    print("     Path-3 exact computation from F89_TOPOLOGY_ORBIT_CLOSURE.md:")
    print("     |S_c(2)|^2 = (10+4*sqrt(5))/3 = 6.3148...")
    print("     ||Mv(2)||^2 = (25+4*sqrt(5))/15 = 2.2630...")
    print("     p_2 = (10+4*sqrt(5))/3 * (25+4*sqrt(5))/15 / 2")
    print("         = (250+40*sqrt(5)+100*sqrt(5)+80)/90")
    print("         = (330+140*sqrt(5))/90")
    print("         = (33+14*sqrt(5))/9")
    print("     Denominator 9 = 3^2 = odd(3)^2 = D_3. Confirmed.")
    print()
    import math
    s5 = math.sqrt(5)
    p2_check = (33 + 14*s5) / 9
    sc2_check = (10 + 4*s5) / 3
    mv2_check = (25 + 4*s5) / 15
    prod_check = sc2_check * mv2_check / 2
    print(f"     Verification: (10+4*sqrt5)/3 * (25+4*sqrt5)/15 / 2 = {prod_check:.6f}")
    print(f"     (33+14*sqrt5)/9                                      = {p2_check:.6f}")
    print(f"     Match: {abs(prod_check - p2_check) < 1e-10}")
    print()

    print("ODD-PART DENOMINATOR ARGUMENT (structural sketch):")
    print()
    print("  S_c(n) and Mv(n) both depend on sums of OBC tight-binding amplitudes.")
    print("  The OBC amplitudes for a (k+2)-site chain are:")
    print("    psi_n[j] = sqrt(2/(k+2)) * sin(pi*n*j/(k+2))  (j = 1..k+1)")
    print()
    print("  |S_c(n)|^2 and ||Mv(n)||^2 are sums of products of pairs of such amplitudes.")
    print("  Each sum involves (k+2) terms; the overall normalization is 1/(k+2).")
    print("  The ratio |S_c(n)|^2 * ||Mv(n)||^2 = A + B*y_n (rational polynomial in y_n)")
    print("  has denominator k^2 = odd(k)^2 (for odd k) from the Bloch sum structures.")
    print("  For even k: the additional factor of 2^{v2(k)} from the 2J hopping")
    print("  and the Vandermonde degree growth accounts for the rest of the denominator.")
    print()
    print("  STATUS: Sketch. The exact algebraic derivation of the k^2 denominator")
    print("  from sine-sum identities is open.")
    print()

    print("NUMERICAL VERIFICATION (p_n = |S_c(n)|^2 * ||Mv(n)||^2 / 2):")
    print()
    import numpy as np

    J_val = 1.0
    gamma_val = 0.05
    N_test = 9

    for k in [3, 4, 5, 6]:
        n_block = k + 1
        m = k + 2
        orbit = list(range(2, n_block + 1, 2))

        L, basis, de_pairs = _build_se_de_L_correct(n_block, J_val, gamma_val)
        eigvals, eigvecs = np.linalg.eig(L)
        fa_mask = np.abs(eigvals.real / gamma_val + 2.0) < 0.01

        N = max(N_test, n_block + 2)
        pre = math.sqrt(2.0 / (N * N * (N - 1)))
        rho_flat = np.full(len(basis), pre / 2, dtype=complex)

        # Per-site reduction w
        w = np.zeros((n_block, len(basis)), dtype=float)
        for b_idx, (i, jk) in enumerate(basis):
            jj, kk = jk
            if kk == i:
                w[jj, b_idx] = 1.0
            if jj == i:
                w[kk, b_idx] = 1.0

        print(f"  k={k} (n_block={n_block}, m={m}):")
        used = set()
        for n in orbit:
            y_n = 4.0 * math.cos(math.pi * n / m)
            best_idx, best_dist = -1, 1e10
            for idx in np.where(fa_mask)[0]:
                if idx in used:
                    continue
                dist = abs(eigvals[idx].imag / J_val - y_n)
                if dist < best_dist:
                    best_dist, best_idx = dist, idx
            if best_idx < 0 or best_dist >= 1e-3:
                print(f"    n={n}: no F_a mode found (best_dist={best_dist:.4f})")
                continue
            used.add(best_idx)
            vv = eigvecs[:, best_idx]
            vv = vv / np.linalg.norm(vv)

            # Standard sigma (reference)
            c = np.vdot(vv, rho_flat)
            Mv = w @ vv
            sigma_std = float(abs(c)**2 * np.sum(abs(Mv)**2))
            p_std = sigma_std * N * N * (N - 1)

            # Formula: p_n = |S_c|^2 * ||Mv||^2 / 2
            S_c = np.sum(vv)
            Mv_sq = float(np.sum(abs(Mv)**2))
            p_formula = abs(S_c)**2 * Mv_sq / 2.0

            print(f"    n={n}: y_n={y_n:+.4f}, p_std={p_std:.6f}, "
                  f"|S_c|^2*||Mv||^2/2={p_formula:.6f}, "
                  f"diff={abs(p_std - p_formula):.2e}")
    print()


# ---------------------------------------------------------------------------
# ANGLE B: Cyclotomic discriminant approach
# ---------------------------------------------------------------------------

def angle_b_cyclotomic() -> None:
    print("=" * 100)
    print("ANGLE B: Cyclotomic discriminant vs D_k")
    print("=" * 100)
    print()
    print("Testing: does disc(minimal poly of 4*cos(pi/m)) relate to D_k?")
    print()

    try:
        import sympy as sp
        y = sp.Symbol("y")
        print(f"{'k':>4}  {'m':>4}  {'D_k':>10}  {'D_k-factored':>16}  {'disc':>14}  {'disc%D':>8}  divisible?")
        print("-" * 80)
        for k in range(3, 13):
            m = k + 2
            D, _ = predicted_D(k)
            try:
                mp = sp.minimal_polynomial(4 * sp.cos(sp.pi / sp.Integer(m)), y)
                disc_val = int(sp.discriminant(mp, y))
                divisible = (abs(disc_val) % D == 0)
                ratio = abs(disc_val) // D if divisible else "N/A"
                print(f"{k:>4}  {m:>4}  {D:>10}  {factor_str(D):>16}  {disc_val:>14}  "
                      f"{abs(disc_val) % D:>8}  {'YES' if divisible else 'NO'}")
            except Exception as e:
                print(f"{k:>4}  {m:>4}  {D:>10}  ERROR: {e}")
    except ImportError:
        print("sympy not available")
    print()
    print("Conclusion: disc(min poly) is generally much larger than D_k and no simple")
    print("divisibility relation exists. Cyclotomic discriminant is the wrong angle.")
    print()


# ---------------------------------------------------------------------------
# ANGLE C: Vandermonde determinant
# ---------------------------------------------------------------------------

def angle_c_vandermonde() -> None:
    print("=" * 100)
    print("ANGLE C: Vandermonde determinant squared |V_det|^2 vs D_k")
    print("=" * 100)
    print()
    print("Testing: is |V_det|^2 proportional to D_k?")
    print()

    import math
    print(f"{'k':>4}  {'D_k':>10}  {'|V|^2':>15}  {'|V|^2/D_k':>12}  integer?")
    print("-" * 60)
    for k in range(3, 16):
        n_block = k + 1
        m = k + 2
        orbit = list(range(2, n_block + 1, 2))
        D, _ = predicted_D(k)

        y_vals = [4.0 * math.cos(math.pi * n / m) for n in orbit]
        V_sq = 1.0
        for i in range(len(y_vals)):
            for j in range(i + 1, len(y_vals)):
                V_sq *= (y_vals[i] - y_vals[j]) ** 2

        ratio = V_sq / D
        is_int = abs(ratio - round(ratio)) < 0.01
        print(f"{k:>4}  {D:>10}  {V_sq:>15.2f}  {ratio:>12.4f}  {'YES' if is_int else 'NO'}")
    print()
    print("Conclusion: |V_det|^2 >> D_k. The Vandermonde determinant grows much faster.")
    print("Algebraic cancellations in the Cramer formulas reduce |V_det|^2 to D_k.")
    print("Specifically: p_n = A + B*y_n (rational linear in y_n for k=3,4),")
    print("so (p_i - p_j)/(y_i - y_j) = B (the slope), eliminating (y_i - y_j) from Cramer.")
    print("This 'rational-polynomial collapse' accounts for the reduction |V|^2 -> D.")
    print()

    # Verify p_n linear in y_n for k=3
    import math
    print("Verification: p_n = A + B*y_n (rational poly) for k=3 and k=4:")
    # k=3: P=[47,14], D=9. So p_n = (47+14*y_n)/9.
    print("  k=3: P=[47,14]/9. p_n = (47 + 14*y_n)/9.")
    for n in [2, 4]:
        y_n = 4.0 * math.cos(math.pi * n / 5)
        p_pred = (47 + 14 * y_n) / 9
        print(f"    n={n}: y_n={y_n:+.4f}, p_n={p_pred:.6f}")
    # k=4: P=[25,10], D=4. So p_n = (25+10*y_n)/4.
    print("  k=4: P=[25,10]/4. p_n = (25 + 10*y_n)/4.")
    for n in [2, 4]:
        y_n = 4.0 * math.cos(math.pi * n / 6)
        p_pred = (25 + 10 * y_n) / 4
        print(f"    n={n}: y_n={y_n:+.4f}, p_n={p_pred:.6f}")
    print()


# ---------------------------------------------------------------------------
# ANGLE D: 2-adic bonus analysis
# ---------------------------------------------------------------------------

def angle_d_2adic() -> None:
    print("=" * 100)
    print("ANGLE D: 2-adic bonus max(0, v2(k)-2) structural interpretation")
    print("=" * 100)
    print()
    print("The bonus max(0, v2(k)-2) kicks in at:")
    print("  k=8  (v2=3, bonus=1): D=32=2^5.  Base 2^{1+3}=2^4=16; bonus adds 2^1=2. 16*2=32.")
    print("  k=16 (v2=4, bonus=2): D=2048=2^11. Base 2^{5+4}=2^9=512; bonus adds 2^2=4. 512*4=2048.")
    print("  k=24 (v2=3, bonus=1): D=73728=2^13*9. Base 2^{9+3}=2^12=4096*9; bonus adds 2.")
    print()
    print("For k = 2^a * odd(k) (a = v2(k)):")
    print("  E(k) = max(0,(k-5)//2) + a + max(0, a-2)")
    print("  2^E(k) = 2^{(FA-3)+} * 2^a * 2^{max(0,a-2)}")
    print("         = 2^{(FA-3)+} * 2^{max(a, 2a-2)}")
    print("  For a >= 3: 2^{max(a,2a-2)} = 2^{2a-2} = (2^a)^2 / 4 = (k/odd(k))^2 / 4.")
    print()
    print("Structural claim (Tier-2 candidate):")
    print("  The deep-2-power bonus arises from the 2-adic over-ramification in the")
    print("  path-k Bloch momentum grid when k = 2^a * odd, a >= 3.")
    print("  When k has 3+ factors of 2, the OBC Bloch energies 4*cos(pi*n/(k+2))")
    print("  for the overlap-pair eigenstates acquire additional powers of 2 in their")
    print("  normalization sums, beyond the base 2^a contribution from the 2J couplings.")
    print()
    print("Verification that the bonus correctly accounts for k=8, k=16, k=24:")
    print()
    for k in [8, 16, 24]:
        vk = v2(k)
        ok = odd_part(k)
        E = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
        D_pred, _ = predicted_D(k)
        D_ver = VERIFIED_D[k]
        base_E = max(0, (k - 5) // 2) + vk  # E without deep-2-power bonus
        bonus = max(0, vk - 2)
        print(f"  k={k}: v2={vk}, odd={ok}, base-E={base_E}, bonus={bonus}, total-E={E}")
        print(f"    D_pred={D_pred} ({factor_str(D_pred)}), D_verified={D_ver} ({factor_str(D_ver)}), "
              f"match={'YES' if D_pred == D_ver else 'NO'}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print()
    print("F89 PATH-D DENOMINATOR THEORY PROBE")
    print("Investigating structural origin of D_k = (odd(k))^2 * 2^E(k)")
    print()

    print_full_table()
    angle_a_structural_argument()
    angle_b_cyclotomic()
    angle_c_vandermonde()
    angle_d_2adic()

    print("=" * 100)
    print("SUMMARY OF FINDINGS")
    print("=" * 100)
    print()
    print("Angle A (Free-fermion / Bloch structure):")
    print("  MAIN RESULT. p_n = sigma_n*N^2*(N-1) = |S_c(n)|^2 * ||Mv(n)||^2 / 2 where")
    print("  S_c(n) = sum of all F_a eigenvector entries (overlap with uniform state),")
    print("  Mv(n) = per-site reduction of v_n.")
    print("  The denominator D_k has:")
    print("  - Odd part odd(k)^2: from |S_c(n)|^2 * ||Mv(n)||^2 having denominator k^2 in Z[y_n].")
    print("    Comes from Bloch amplitude sums over (k+2) sites with normalization 2/(k+2).")
    print("    For path-3: exact algebraic derivation gives (10+4*sqrt5)/3 * (25+4*sqrt5)/15 / 2")
    print("    = (33+14*sqrt5)/9, denominator 9 = 3^2 = odd(3)^2.")
    print("  - 2-power 2^E(k): from (1) the 2J hopping coefficients, (2) the Vandermonde")
    print("    degree growth, and (3) the deep-2-power bonus for k with v2(k)>=3.")
    print("  STATUS: STRUCTURAL INSIGHT (Tier-2-strong). Not a full algebraic proof.")
    print("    The p_n = |S_c(n)|^2 * ||Mv(n)||^2 / 2 identity verified numerically for k=3..6.")
    print("    The exact path-3 algebraic derivation is complete (denominator = 9 = 3^2).")
    print("    Full generalization to all k requires explicit sine-sum computation.")
    print()
    print("Angle B (Cyclotomic discriminant):")
    print("  NEGATIVE. disc(min poly) >> D_k, no divisibility relation.")
    print()
    print("Angle C (Vandermonde det^2):")
    print("  NEGATIVE. |V|^2 >> D_k due to algebraic cancellations.")
    print("  However, the cancellation mechanism IS understood: p_n is a rational poly")
    print("  in y_n (shown for k=3,4,5,7), so the irrational parts of V_det cancel")
    print("  in the Cramer coefficient ratios.")
    print()
    print("Angle D (2-adic bonus):")
    print("  PARTLY UNDERSTOOD. The bonus max(0,v2(k)-2) is supported by 3 data points")
    print("  (k=8,16,24). Structural interpretation: extra 2-adic ramification in Bloch")
    print("  grid when k=2^a*odd with a>=3. Full derivation is open.")
    print()


if __name__ == "__main__":
    main()
