"""F89 path-D structural probe: cyclotomic/Galois correlates for the D_path denominator.

Investigates the structural origin of D_path in:
  sigs[F_a:n](N) = P_path(y_n) / [D_path * N^2 * (N-1)]
  y_n = 4*cos(pi*n/(N_block+1))   for n in SE-anti Bloch orbit

Known data from F89UnifiedFaClosedFormClaim (path-3..7):
  k=3: (P=[47,14],     D=9  = 3^2)       nBlock=4
  k=4: (P=[25,10],     D=4  = 2^2)       nBlock=5
  k=5: (P=[129,82,13], D=25 = 5^2)       nBlock=6
  k=6: (P=[80,72,17],  D=18 = 2*3^2)     nBlock=7
  k=7: (P=[382,292,130,21], D=98=2*7^2)  nBlock=8

This probe extracts (P, D) numerically for paths 3..15 using the (SE,DE) sub-block
Liouvillian eigendecomposition, then pairs results with cyclotomic/Galois invariants
for m = nBlock+1 via sympy.

Method (validated against path-3..6 typed claims in F89UnifiedFaClosedFormClaim):
  1. Build (SE,DE) sub-block L for path-k (nBlock = k+1 sites, SE-anti orbit).
  2. Eigendecompose. Filter F_a modes by Re(eig) near -2*gamma.
  3. Find mode per Bloch index n by matching Im(eig)/J to 4*cos(pi*n/(nBlock+1)).
  4. Compute sigma_n = |c_n|^2 * sum_l |Mv_n[l]|^2, where:
       v_n = normalized F_a eigenvector for orbit index n
       c_n = vdot(v_n, rho_flat)  with rho_flat = pre/2 * ones, pre = sqrt(2/(N^2*(N-1)))
       Mv = w @ v, w = per-site reduction matrix for (SE,DE) basis
  5. sigma_n * N^2*(N-1) is N-invariant; call it p_n.
  6. The Vandermonde system P(y_n) = D * p_n is solved for integer coefs by finding
     smallest D in [1, D_MAX] making D*p_n near-integer for all n simultaneously.
  7. Verify Vandermonde fit residual < 1e-6.
  8. Pair with cyclotomic invariants from sympy.

Runtime:
  nBlock=4..9 (SE,DE) sub-block dimension = nBlock * C(nBlock,2) ranges 12..360.
  nBlock=10..16 extends to dimension 1980 (nBlock=16). All eigendecomps in dense
  complex128 via np.linalg.eig; estimated ~5-60 s per path. Total ~5 min.

Extended range: paths 3..15 (nBlock 4..16).

Hypotheses under test (2-adic structure of D):
  Odd-k candidate closed form: D_odd_k = 2^max(0,(k-5)/2) * k^2  (integer division)
    Predicts: k=11 -> 968, k=13 -> 2704, k=15 -> 7200
  Even-k hypotheses (fits k=4,6; differ at k=8,10,12,14):
    (a) v2(D) = v2(k)
    (b) v2(D) = v2(k) + max(0, FA-3)
    (c) ad-hoc special bonus for pure 2-powers
  Odd-part rule (all k): odd(D) = (odd(k))^2  -- verified k=3..9.
"""
from __future__ import annotations

import sys
import math
from pathlib import Path
from itertools import combinations

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

# Ensure ASCII-safe output (avoid encoding errors on Windows console)
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
J = 1.0
GAMMA = 0.05
N_FIXED = 9  # fixed N for the c_n inner-product (N-invariant after scaling)

# F_a filter: |Re(eig)/gamma + 2| < FA_TOL_REL
FA_TOL_REL = 0.01

# Bloch-index matching tolerance: |Im(eig)/J - y_n| < BLOCH_TOL
BLOCH_TOL = 1e-3

# D search range
D_MAX = 200000

# Integer rationalization tolerance
RATIONAL_TOL = 1e-3

# Typed reference (path-3..9) from F89UnifiedFaClosedFormClaim.cs + prior run
TYPED = {
    3: (9,   [47, 14]),
    4: (4,   [25, 10]),
    5: (25,  [129, 82, 13]),
    6: (18,  [80,  72, 17]),
    7: (98,  [382, 292, 130, 21]),
    8: (32,  None),   # D=32=2^5 from prior run; coefs not typed
    9: (324, None),   # D=324=2^2*3^4 from prior run; coefs not typed
}

# Path range: 3..K_MAX (inclusive)
K_MAX = 15

# ---------------------------------------------------------------------------
# 2-adic / odd-part helpers
# ---------------------------------------------------------------------------

def v2(n: int) -> int:
    """2-adic valuation of n (v2(0) defined as 0)."""
    if n <= 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def odd_part(n: int) -> int:
    """Odd part of n: n / 2^v2(n)."""
    if n <= 0:
        return n
    while n % 2 == 0:
        n //= 2
    return n


def predicted_D_odd_k(k: int) -> int:
    """Candidate odd-k closed form: 2^max(0,(k-5)//2) * k^2."""
    exp = max(0, (k - 5) // 2)
    return (2 ** exp) * k * k


# ---------------------------------------------------------------------------
# SE-anti Bloch orbit
# ---------------------------------------------------------------------------
def se_anti_orbit(n_block: int) -> list[int]:
    """S_2-anti Bloch orbit: even n in [2, n_block], step 2."""
    return list(range(2, n_block + 1, 2))


def bloch_y(n: int, n_block: int) -> float:
    """y_n = 4*cos(pi*n/(n_block+1))."""
    return 4.0 * math.cos(math.pi * n / (n_block + 1))


# ---------------------------------------------------------------------------
# (SE,DE) sub-block Liouvillian
# ---------------------------------------------------------------------------
def build_se_de_L(n_block: int, J_val: float, gamma_val: float):
    """Build (SE,DE) sub-block Liouvillian for path-(n_block-1).

    Basis: all (SE-site i, DE-pair (j,k)) with j < k.
    Diagonal: -2*gamma if i in {j,k} (HD=1), else -6*gamma (HD=3).
    Off-diagonal sign convention (matches path-3..6 AT-locked amplitude scripts):
      SE-site hop i->i2: coefficient -i*M_SE[i2,i] = -i*2J
      DE-pair hop jk->jk2: coefficient +i*M_DE[jk2,jk] = +i*2J

    Returns (L, basis, de_pairs).
    """
    de_pairs = [(j, k) for j in range(n_block) for k in range(j + 1, n_block)]
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    n_basis = len(basis)

    # SE hopping matrix (coupling amplitude 2J between adjacent sites)
    M_SE = np.zeros((n_block, n_block))
    for a in range(n_block):
        for b in range(n_block):
            if abs(a - b) == 1:
                M_SE[a, b] = 2.0 * J_val

    # DE hopping matrix (coupling amplitude 2J between DE pairs sharing one site)
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
        # SE-site hopping: off-diagonal -i*M_SE[i2, i]
        for i2 in range(n_block):
            if M_SE[i2, i] != 0.0:
                idx2 = basis.index((i2, jk))
                L[idx2, idx] += -1j * M_SE[i2, i]
        # DE-pair hopping: off-diagonal +i*M_DE[jk2_idx, jk_idx]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(n_de):
            if M_DE[jk_idx, jk2_idx] != 0.0:
                jk2 = de_pairs[jk2_idx]
                idx2 = basis.index((i, jk2))
                L[idx2, idx] += 1j * M_DE[jk_idx, jk2_idx]
        # Diagonal: decay rate
        L[idx, idx] += -2.0 * gamma_val if i in jk else -6.0 * gamma_val

    return L, basis, de_pairs


def per_site_reduction_se_de(basis: list, de_pairs: list, n_block: int) -> np.ndarray:
    """Per-site reduction matrix w for (SE,DE) basis.

    w[l, b_idx] = 1 if basis[b_idx] = (i, jk) with:
      l in {j, k} = jk AND the other element of {j,k} equals i.
    This selects pairs where the SE excitation overlaps the DE pair at site l.
    """
    w = np.zeros((n_block, len(basis)), dtype=float)
    for b_idx, (i, jk) in enumerate(basis):
        j, k = jk
        # Contribute to w[j] if k == i (SE site is k, DE site is j)
        if k == i:
            w[j, b_idx] = 1.0
        # Contribute to w[k] if j == i (SE site is j, DE site is k)
        if j == i:
            w[k, b_idx] = 1.0
    return w


# ---------------------------------------------------------------------------
# Sigma extraction per Bloch index via (SE,DE) eigendecomposition
# ---------------------------------------------------------------------------
def extract_sigma_scaled(n_block: int) -> dict[int, float]:
    """Extract p_n = sigma_n * N^2*(N-1) for each SE-anti Bloch index n.

    Uses the (SE,DE) sub-block eigendecomposition with N = N_FIXED.
    The product p_n is N-invariant (verified for path-3..6 in prior scripts).

    Returns {n: p_n}.
    """
    N = max(N_FIXED, n_block + 2)  # need N > n_block for at least one bare site
    orbit = se_anti_orbit(n_block)
    y_targets = {n: bloch_y(n, n_block) for n in orbit}

    L, basis, de_pairs = build_se_de_L(n_block, J, GAMMA)
    eigvals, eigvecs = np.linalg.eig(L)

    # Filter F_a: Re(eig) near -2*gamma
    fa_mask = np.abs(eigvals.real / GAMMA + 2.0) < FA_TOL_REL

    # Uniform flat rho vector: all (SE,DE) entries = pre/2
    pre = math.sqrt(2.0 / (N * N * (N - 1)))
    rho_flat = np.full(len(basis), pre / 2, dtype=complex)

    w = per_site_reduction_se_de(basis, de_pairs, n_block)

    # For each Bloch index n, find the matching F_a eigenmode and accumulate sigma.
    # Matching is DIRECT: Im(eig)/J must equal y_n (which can be positive or negative).
    # y_n = 4*cos(pi*n/(nBlock+1)) for n in SE-anti orbit (even n, 2..nBlock).
    # The modes at Im/J = y_n (not conjugate -y_n) carry the dominant sigma.
    sigma_scaled: dict[int, float] = {}
    used_idx: set[int] = set()

    for n, y_n in sorted(y_targets.items()):
        best_idx = -1
        best_dist = 1e10
        for idx in np.where(fa_mask)[0]:
            if idx in used_idx:
                continue
            # Direct match only: Im(eig)/J vs y_n
            dist = abs(eigvals[idx].imag / J - y_n)
            if dist < best_dist:
                best_dist = dist
                best_idx = idx

        if best_idx < 0 or best_dist >= BLOCH_TOL:
            sigma_scaled[n] = float("nan")
            continue

        v = eigvecs[:, best_idx]
        v = v / np.linalg.norm(v)
        c = np.vdot(v, rho_flat)
        Mv = w @ v
        sig = abs(c) ** 2 * np.sum(np.abs(Mv) ** 2)
        sigma_scaled[n] = float(sig.real) * N * N * (N - 1)
        used_idx.add(best_idx)

    return sigma_scaled


# ---------------------------------------------------------------------------
# Rationalization: find D via Vandermonde coefficient rationalization
# ---------------------------------------------------------------------------
def rationalize_polynomial_D(
    p_vals: list[float],
    y_vals: list[float],
    d_max: int = D_MAX,
) -> tuple[int, list[int]]:
    """Find smallest D making all Vandermonde polynomial coefficients near-integer.

    Given sigma-scaled values p_n and Bloch y_n, fit P(y) with D=1 first (so P(y_n) = p_n),
    then find smallest D in [1, d_max] such that all coefficients of D*P(y) are near-integer.

    Returns (D, coefs_low_to_high) or (-1, []) if not found.
    """
    deg = len(y_vals) - 1
    V = np.vander(y_vals, N=deg + 1, increasing=True)
    coefs_d1, _, _, _ = np.linalg.lstsq(V, np.array(p_vals, dtype=float), rcond=None)

    for D in range(1, d_max + 1):
        scaled = D * coefs_d1
        if all(abs(s - round(s)) < RATIONAL_TOL for s in scaled):
            int_coefs = [round(float(s)) for s in scaled]
            # Verify the fit: check P(y_n) = D * p_n for all n
            residuals = [
                abs(sum(int_coefs[i] * (y ** i) for i in range(len(int_coefs))) - D * p)
                for y, p in zip(y_vals, p_vals)
            ]
            if max(residuals) < 1e-1:
                return D, int_coefs
    return -1, []


# ---------------------------------------------------------------------------
# Cyclotomic/Galois invariants via sympy
# ---------------------------------------------------------------------------
def cyclotomic_invariants(m: int) -> dict:
    """Cyclotomic/Galois invariants for Q(zeta_m).

    m = n_block + 1.
    Returns dict with phi, real_sub_deg, min_poly_4cos_y1, disc_formula.
    """
    try:
        import sympy as sp
        y = sp.Symbol("y")
        phi = int(sp.totient(m))
        real_sub_deg = phi // 2

        # Minimal polynomial of 4*cos(pi/m) = y_1 (first SE-anti Bloch frequency)
        try:
            min_poly_sym = sp.minimal_polynomial(4 * sp.cos(sp.pi / sp.Integer(m)), y)
            min_poly_str = str(min_poly_sym)
        except Exception:
            min_poly_str = "?"

        # Discriminant formula for Q(zeta_m):
        # disc(Q(zeta_n)) = (-1)^(phi(n)/2) * n^phi(n) / prod_{p|n} p^(phi(n)/(p-1))
        m_factors = sp.factorint(m)
        disc_sign = (-1) ** (phi // 2)
        disc_num = m ** phi
        disc_den = 1
        for p_factor in m_factors:
            disc_den *= p_factor ** (phi // (p_factor - 1))
        disc_val = disc_sign * disc_num // disc_den
        disc_formula = f"(-1)^({phi//2}) * {m}^{phi} / {disc_den} = {disc_val}"

        return {
            "phi": phi,
            "real_sub_deg": real_sub_deg,
            "min_poly_4cos_y1": min_poly_str,
            "disc_formula": disc_formula,
            "disc_val": disc_val,
        }
    except ImportError:
        return {
            "phi": -1, "real_sub_deg": -1,
            "min_poly_4cos_y1": "(sympy not available)",
            "disc_formula": "?", "disc_val": 0,
        }


# ---------------------------------------------------------------------------
# Prime factorization helpers
# ---------------------------------------------------------------------------
def factorize(n: int) -> dict[int, int]:
    """Trial division prime factorization."""
    factors: dict[int, int] = {}
    d = 2
    tmp = n
    while d * d <= tmp:
        while tmp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            tmp //= d
        d += 1
    if tmp > 1:
        factors[tmp] = factors.get(tmp, 0) + 1
    return factors


def factor_str(n: int) -> str:
    """Compact factorization string: '2*3^2'."""
    if n <= 0:
        return str(n)
    if n == 1:
        return "1"
    facs = factorize(n)
    parts = []
    for p in sorted(facs):
        e = facs[p]
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts)


def largest_prime_factor(n: int) -> int:
    """Largest prime factor of n."""
    facs = factorize(n)
    return max(facs.keys()) if facs else 1


def is_prime(n: int) -> bool:
    """Simple primality test."""
    if n < 2:
        return False
    facs = factorize(n)
    return len(facs) == 1 and list(facs.values())[0] == 1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 110)
    print("F89 path-D structure probe: cyclotomic/Galois correlates for D_path (extended paths 3..15)")
    print("=" * 110)
    print(f"Parameters: J={J}, gamma={GAMMA}, N_fixed={N_FIXED}")
    print(f"D search range: [1, {D_MAX}], rationalization tolerance: {RATIONAL_TOL}")
    print()

    results: dict[int, dict] = {}

    print(f"Extracting sigma_n for paths 3..{K_MAX} ...")
    print()

    for k in range(3, K_MAX + 1):
        n_block = k + 1
        orbit = se_anti_orbit(n_block)
        fa_count = len(orbit)
        m = n_block + 1
        dim = n_block * (n_block * (n_block - 1) // 2)
        print(
            f"  path-{k}: nBlock={n_block}, m={m}, orbit={orbit} "
            f"(FA count={fa_count}, (SE,DE) dim={dim})",
            end=" ... ", flush=True,
        )

        sigma_scaled = extract_sigma_scaled(n_block)
        inv = cyclotomic_invariants(m)

        # Check for NaN
        nan_ns = [n for n, vv in sigma_scaled.items() if math.isnan(vv)]
        if nan_ns:
            print(f"WARNING: NaN sigma at n={nan_ns}")
            results[k] = {"error": f"NaN at n={nan_ns}", "n_block": n_block, "m": m,
                          "fa_count": fa_count, "inv": inv, "sigma_scaled": sigma_scaled}
            continue

        p_vals = [sigma_scaled[n] for n in orbit]
        y_vals = [bloch_y(n, n_block) for n in orbit]

        D_found, coefs = rationalize_polynomial_D(p_vals, y_vals, D_MAX)

        if D_found < 0:
            print(f"D NOT FOUND in [1,{D_MAX}]  p_vals={[f'{vv:.6f}' for vv in p_vals]}")
            results[k] = {"error": "D not found", "n_block": n_block, "m": m,
                          "fa_count": fa_count, "inv": inv, "sigma_scaled": sigma_scaled,
                          "D": -1, "coefs": [], "y_vals": y_vals}
            continue

        # Verify against typed data (D only for k=8,9 where coefs are not typed)
        if k in TYPED:
            typed_D, typed_coefs = TYPED[k]
            D_ok = (D_found == typed_D)
            if typed_coefs is not None and coefs:
                coefs_ok = (coefs == typed_coefs)
                status = "OK" if (D_ok and coefs_ok) else f"MISMATCH(D:{D_found}vs{typed_D},c:{coefs}vs{typed_coefs})"
            else:
                status = "D-OK" if D_ok else f"D-MISMATCH({D_found}vs{typed_D})"
            print(f"D={D_found} [{factor_str(D_found)}]  coefs={coefs}  verify={status}")
        else:
            print(f"D={D_found} [{factor_str(D_found)}]  coefs={coefs}  (NEW DATA)")

        p = largest_prime_factor(k)
        results[k] = {
            "n_block": n_block, "m": m, "fa_count": fa_count,
            "orbit": orbit, "sigma_scaled": sigma_scaled,
            "D": D_found, "coefs": coefs, "y_vals": y_vals,
            "inv": inv,
            "p_lpf_k": p,
            "d_pred_p2": p * p,
            "d_pred_2p2": 2 * p * p,
        }

    print()

    # ---------------------------------------------------------------------------
    # Print the main table (extended with hypothesis-test columns)
    # ---------------------------------------------------------------------------
    print("=" * 110)
    print("FULL TABLE: paths 3..15 with 2-adic hypothesis columns")
    print("=" * 110)

    # Columns: k, nBlock, FA-cnt, D, D-factored, v2(D), v2(k), odd(D), (odd(k))^2,
    #          odd-part-rule-ok, pred_D_odd_k (odd k only), odd-k-formula-match
    col_w = [4, 6, 6, 8, 14, 6, 5, 8, 9, 14, 14, 18]
    headers = [
        "k", "nBlock", "FA", "D", "D-factored",
        "v2(D)", "v2(k)", "odd(D)", "(odd(k))^2",
        "odd-part-rule?", "pred_D_odd_k", "odd-k-formula?"
    ]
    sep = "  "
    header_line = sep.join(h.ljust(w) for h, w in zip(headers, col_w))
    print(header_line)
    print("-" * len(header_line))

    for k in range(3, K_MAX + 1):
        if "error" in results.get(k, {}):
            print(f"  k={k}: ERROR - {results[k]['error']}")
            continue
        r = results.get(k)
        if r is None:
            print(f"  k={k}: no result")
            continue
        D = r["D"]
        fa = r["fa_count"]

        v2_D = v2(D)
        v2_k = v2(k)
        odd_D = odd_part(D)
        ok = k % 2  # 1 if k is odd
        odd_k = odd_part(k)
        odd_k_sq = odd_k * odd_k
        odd_part_rule_ok = "YES" if odd_D == odd_k_sq else f"NO({odd_D}!={odd_k_sq})"

        if k % 2 == 1:  # odd k
            pred = predicted_D_odd_k(k)
            match = "YES" if D == pred else f"NO(pred={pred})"
            pred_str = str(pred)
        else:
            pred_str = "n/a"
            match = "skip"

        row = [
            str(k), str(r["n_block"]), str(fa), str(D), factor_str(D),
            str(v2_D), str(v2_k), str(odd_D), str(odd_k_sq),
            odd_part_rule_ok, pred_str, match,
        ]
        print(sep.join(v.ljust(w) for v, w in zip(row, col_w)))

    print()
    print("Polynomial coefs P(y) low->high for each path:")
    for k in range(3, K_MAX + 1):
        r = results.get(k, {})
        if "error" in r:
            print(f"  k={k}: {r.get('error', 'unknown error')}")
        elif r:
            typed_ref = ""
            if k in TYPED:
                tD, tC = TYPED[k]
                typed_ref = f"  [typed D={tD}" + (f", coefs={tC}]" if tC else "]")
            else:
                typed_ref = "  [new]"
            print(f"  k={k}: D={r['D']}, P(y)={r['coefs']}{typed_ref}")

    print()
    print("Cyclotomic minimal polynomial of 4*cos(pi/m) over Q (= y_1, first orbit element):")
    for k in range(3, K_MAX + 1):
        r = results.get(k, {})
        inv = r.get("inv", {})
        m = r.get("m", k + 2)
        print(f"  k={k}, m={m}: {inv.get('min_poly_4cos_y1', '?')}")

    print()
    print("Discriminant of Q(zeta_m):")
    for k in range(3, K_MAX + 1):
        r = results.get(k, {})
        inv = r.get("inv", {})
        m = r.get("m", k + 2)
        print(f"  k={k}, m={m}: {inv.get('disc_formula', '?')}")

    print()

    # ---------------------------------------------------------------------------
    # Even-k hypothesis disambiguation table
    # ---------------------------------------------------------------------------
    print("=" * 110)
    print("EVEN-k HYPOTHESIS DISAMBIGUATION: v2(D) vs hypotheses (a), (b), (c)")
    print("  (a) v2(D) = v2(k)")
    print("  (b) v2(D) = v2(k) + max(0, FA-3)")
    print("  (c) ad-hoc: bonus when k is a pure 2-power and FA >= 4")
    print("=" * 110)
    even_ks = [k for k in range(3, K_MAX + 1) if k % 2 == 0]
    col_e = [5, 5, 4, 6, 13, 25, 20, 8]
    hdr_e = ["k", "v2(k)", "FA", "v2(D)", "hyp(a)=v2(k)", "hyp(b)=v2(k)+max(0,FA-3)", "hyp(c)-special?", "status"]
    print(sep.join(h.ljust(w) for h, w in zip(hdr_e, col_e)))
    print("-" * sum(col_e + [2 * (len(col_e) - 1)]))

    hyp_results = {}  # k -> {a_ok, b_ok, c_ok}
    for k in even_ks:
        r = results.get(k, {})
        if "error" in r or not r:
            row_e = [str(k), str(v2(k)), "?", "?", "?", "?", "?", r.get("error", "no data")]
            print(sep.join(v.ljust(w) for v, w in zip(row_e, col_e)))
            continue
        D = r["D"]
        fa = r["fa_count"]
        v2k = v2(k)
        v2D = v2(D)
        # Hypothesis predictions
        pred_a = v2k
        pred_b = v2k + max(0, fa - 3)
        # (c): pure-2-power special rule: if k = 2^e for some e and FA >= 4, add extra
        k_facs = factorize(k)
        is_pure2power = (len(k_facs) == 1 and 2 in k_facs)
        pred_c_str = f"special({v2D})" if is_pure2power and fa >= 4 else f"same-as-(a)={pred_a}"

        a_ok = (v2D == pred_a)
        b_ok = (v2D == pred_b)
        c_ok = v2D == (v2D if (is_pure2power and fa >= 4) else pred_a)  # tautological for special

        a_str = f"{pred_a} {'YES' if a_ok else 'NO'}"
        b_str = f"{pred_b} {'YES' if b_ok else 'NO'}"

        hyp_results[k] = {"a": a_ok, "b": b_ok, "c_special": is_pure2power and fa >= 4,
                          "v2k": v2k, "v2D": v2D, "fa": fa, "D": D}

        status = []
        if a_ok:
            status.append("(a)")
        if b_ok:
            status.append("(b)")
        if is_pure2power and fa >= 4:
            status.append("(c)-case")
        status_str = "+".join(status) if status else "none"

        row_e = [str(k), str(v2k), str(fa), str(v2D), a_str, b_str, pred_c_str, status_str]
        print(sep.join(v.ljust(w) for v, w in zip(row_e, col_e)))

    print()
    # Summary
    print("Even-k hypothesis summary:")
    all_even_valid = [k for k in even_ks if k in hyp_results]
    hyp_a_all = all(hyp_results[k]["a"] for k in all_even_valid)
    hyp_b_all = all(hyp_results[k]["b"] for k in all_even_valid)
    print(f"  Hypothesis (a) [v2(D)=v2(k)] holds for ALL even k: {hyp_a_all}")
    print(f"  Hypothesis (b) [v2(D)=v2(k)+max(0,FA-3)] holds for ALL even k: {hyp_b_all}")
    if not hyp_a_all:
        fails_a = [k for k in all_even_valid if not hyp_results[k]["a"]]
        print(f"    (a) fails at k={fails_a}: "
              + ", ".join(f"k={k} v2(D)={hyp_results[k]['v2D']} predicted={hyp_results[k]['v2k']}"
                          for k in fails_a))
    if not hyp_b_all:
        fails_b = [k for k in all_even_valid if not hyp_results[k]["b"]]
        print(f"    (b) fails at k={fails_b}: "
              + ", ".join(
                  f"k={k} v2(D)={hyp_results[k]['v2D']} predicted={hyp_results[k]['v2k']+max(0,hyp_results[k]['fa']-3)}"
                  for k in fails_b))

    print()

    # ---------------------------------------------------------------------------
    # Odd-k closed-form formula verification summary
    # ---------------------------------------------------------------------------
    print("=" * 110)
    print("ODD-k CLOSED-FORM FORMULA CHECK: D = 2^max(0,(k-5)//2) * k^2")
    print("=" * 110)
    odd_ks = [k for k in range(3, K_MAX + 1) if k % 2 == 1]
    all_odd_match = True
    for k in odd_ks:
        r = results.get(k, {})
        if "error" in r or not r:
            print(f"  k={k}: ERROR - {r.get('error', 'no data')}")
            all_odd_match = False
            continue
        D = r["D"]
        pred = predicted_D_odd_k(k)
        exp = max(0, (k - 5) // 2)
        match = (D == pred)
        if not match:
            all_odd_match = False
        print(f"  k={k}: predicted 2^{exp}*{k}^2 = {pred}, actual D={D} [{factor_str(D)}]  {'MATCH' if match else 'MISMATCH'}")
    print(f"  All odd k (3..{K_MAX}) match closed form: {all_odd_match}")

    print()

    # ---------------------------------------------------------------------------
    # Pattern observations
    # ---------------------------------------------------------------------------
    print("=" * 110)
    print("PATTERN OBSERVATIONS (observational, not claims)")
    print("=" * 110)

    valid_ks = [k for k in range(3, K_MAX + 1) if k in results and "error" not in results[k]]

    print()
    print("1. Odd-part rule: odd(D) = (odd(k))^2?")
    all_odd_part_ok = True
    for k in valid_ks:
        r = results[k]
        D = r["D"]
        odd_D = odd_part(D)
        odd_k_sq = odd_part(k) ** 2
        ok_flag = (odd_D == odd_k_sq)
        if not ok_flag:
            all_odd_part_ok = False
        print(f"   k={k}: D={D}={factor_str(D)}, odd(D)={odd_D}, (odd({k}))^2={odd_k_sq}, rule: {'YES' if ok_flag else 'NO'}")
    print(f"   Odd-part rule holds for ALL valid k: {all_odd_part_ok}")

    print()
    print("2. FA count equals real_subfield_degree phi(m)/2?")
    for k in valid_ks:
        r = results[k]
        inv = r["inv"]
        fa = r["fa_count"]
        rsd = inv["real_sub_deg"]
        print(f"   k={k}: FA_count={fa}, phi({r['m']})/2={rsd}, equal: {fa == rsd}")

    print()
    print("3. Polynomial degree = FA_count - 1?")
    for k in valid_ks:
        r = results[k]
        coefs = r["coefs"]
        deg = len(coefs) - 1 if coefs else -1
        fa = r["fa_count"]
        print(f"   k={k}: poly_deg={deg}, FA_count-1={fa-1}, consistent: {deg == fa - 1}")

    print()
    print("4. Factor-2 structure in D and k structure:")
    for k in valid_ks:
        r = results[k]
        D = r["D"]
        has2 = (D % 2 == 0) if D > 0 else None
        k_is_prime = is_prime(k)
        k_is_prime_power = (len(factorize(k)) == 1)
        print(
            f"   k={k}={factor_str(k)}: D={D}={factor_str(D)}, "
            f"v2(D)={v2(D)}, v2(k)={v2(k)}, k-prime: {k_is_prime}, k-prime-power: {k_is_prime_power}"
        )

    print()
    print("5. Discriminant of Q(zeta_m) vs D:")
    for k in valid_ks:
        r = results[k]
        inv = r["inv"]
        D = r["D"]
        disc = abs(inv.get("disc_val", 0))
        if disc > 0 and D > 0:
            divides = disc % D == 0
            ratio = disc // D if divides else "N/A"
        else:
            divides = "?"
            ratio = "?"
        print(f"   k={k}: D={D}, |disc(Q(zeta_{r['m']}))| = {disc}, D|disc: {divides}, ratio={ratio}")

    print()
    print("6. P(y_n) values: D * p_n evaluated via fitted polynomial:")
    for k in valid_ks:
        r = results[k]
        orbit = r["orbit"]
        p_vals_map = {n: r["sigma_scaled"][n] for n in orbit}
        n_block = r["n_block"]
        D = r["D"]
        coefs = r["coefs"]
        print(f"   k={k}: p_n = {[f'{p_vals_map[n]:.6f}' for n in orbit]}")
        if D > 0 and coefs:
            P_vals = []
            for n in orbit:
                y_n = bloch_y(n, n_block)
                P_n = sum(coefs[i] * (y_n ** i) for i in range(len(coefs)))
                P_vals.append(round(P_n))
            print(f"          P(y_n) = D*[fit(p_n)] = {P_vals}")

    print()


if __name__ == "__main__":
    main()
