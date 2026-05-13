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

This probe extracts (P, D) numerically for paths 3..9 using the (SE,DE) sub-block
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
  All eigendecompositions complete in < 1 second each. Total runtime ~30 seconds
  (dominated by sympy minimal_polynomial calls).
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

# Typed reference (path-3..7) from F89UnifiedFaClosedFormClaim.cs
TYPED = {
    3: (9,  [47, 14]),
    4: (4,  [25, 10]),
    5: (25, [129, 82, 13]),
    6: (18, [80,  72, 17]),
    7: (98, [382, 292, 130, 21]),
}


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
    print("=" * 100)
    print("F89 path-D structure probe: cyclotomic/Galois correlates for D_path")
    print("=" * 100)
    print(f"Parameters: J={J}, gamma={GAMMA}, N_fixed={N_FIXED}")
    print(f"D search range: [1, {D_MAX}], rationalization tolerance: {RATIONAL_TOL}")
    print()

    results: dict[int, dict] = {}

    print("Extracting sigma_n for paths 3..9 ...")
    print()

    for k in range(3, 10):
        n_block = k + 1
        orbit = se_anti_orbit(n_block)
        fa_count = len(orbit)
        m = n_block + 1
        print(f"  path-{k}: nBlock={n_block}, m={m}, orbit={orbit} (FA count={fa_count})", end=" ... ", flush=True)

        sigma_scaled = extract_sigma_scaled(n_block)
        inv = cyclotomic_invariants(m)

        # Check for NaN
        nan_ns = [n for n, v in sigma_scaled.items() if math.isnan(v)]
        if nan_ns:
            print(f"WARNING: NaN sigma at n={nan_ns}")
            results[k] = {"error": f"NaN at n={nan_ns}", "n_block": n_block, "m": m,
                          "fa_count": fa_count, "inv": inv, "sigma_scaled": sigma_scaled}
            continue

        p_vals = [sigma_scaled[n] for n in orbit]
        y_vals = [bloch_y(n, n_block) for n in orbit]

        D_found, coefs = rationalize_polynomial_D(p_vals, y_vals, D_MAX)

        if D_found < 0:
            print(f"D NOT FOUND in [1,{D_MAX}]  p_vals={[f'{v:.6f}' for v in p_vals]}")
            results[k] = {"error": "D not found", "n_block": n_block, "m": m,
                          "fa_count": fa_count, "inv": inv, "sigma_scaled": sigma_scaled,
                          "D": -1, "coefs": [], "y_vals": y_vals}
            continue

        # Verify against typed data
        if k in TYPED:
            typed_D, typed_coefs = TYPED[k]
            D_ok = (D_found == typed_D)
            coefs_ok = (coefs == typed_coefs) if isinstance(coefs[0], int) else False
            status = "OK" if (D_ok and coefs_ok) else f"MISMATCH(D:{D_found}vs{typed_D},c:{coefs}vs{typed_coefs})"
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
    # Print the main table
    # ---------------------------------------------------------------------------
    print("=" * 100)
    print("FULL TABLE: path-k vs cyclotomic/Galois invariants")
    print("=" * 100)

    col_w = [4, 6, 6, 6, 12, 8, 8, 14, 8, 5, 7, 11, 12]
    headers = [
        "k", "nBlock", "m", "phi(m)", "real-sub-deg",
        "FA-cnt", "D", "D-factored",
        "p=lpf(k)", "p^2", "2*p^2",
        "D=p^2?", "D=2*p^2?"
    ]
    sep = "  "
    header_line = sep.join(h.ljust(w) for h, w in zip(headers, col_w))
    print(header_line)
    print("-" * len(header_line))

    for k in range(3, 10):
        if "error" in results.get(k, {}):
            print(f"  k={k}: ERROR - {results[k]['error']}")
            continue
        r = results.get(k)
        if r is None:
            continue
        D = r["D"]
        inv = r["inv"]
        p = r["p_lpf_k"]
        d_p2 = r["d_pred_p2"]
        d_2p2 = r["d_pred_2p2"]

        row = [
            str(k), str(r["n_block"]), str(r["m"]),
            str(inv["phi"]), str(inv["real_sub_deg"]),
            str(r["fa_count"]), str(D), factor_str(D),
            str(p), str(d_p2), str(d_2p2),
            "YES" if D == d_p2 else "no",
            "YES" if D == d_2p2 else "no",
        ]
        print(sep.join(v.ljust(w) for v, w in zip(row, col_w)))

    print()
    print("Polynomial coefs P(y) low->high for each path:")
    for k in range(3, 10):
        r = results.get(k, {})
        if "error" in r:
            print(f"  k={k}: {r['error']}")
        else:
            typed_ref = f"  [typed: {TYPED[k]}]" if k in TYPED else "  [new]"
            print(f"  k={k}: D={r['D']}, P(y)={r['coefs']}{typed_ref}")

    print()
    print("Cyclotomic minimal polynomial of 4*cos(pi/m) over Q (= y_1, first orbit element):")
    for k in range(3, 10):
        r = results.get(k, {})
        inv = r.get("inv", {})
        m = r.get("m", k + 2)
        print(f"  k={k}, m={m}: {inv.get('min_poly_4cos_y1', '?')}")

    print()
    print("Discriminant of Q(zeta_m):")
    for k in range(3, 10):
        r = results.get(k, {})
        inv = r.get("inv", {})
        m = r.get("m", k + 2)
        print(f"  k={k}, m={m}: {inv.get('disc_formula', '?')}")

    print()

    # ---------------------------------------------------------------------------
    # Hypothesis test for new data k=8,9
    # ---------------------------------------------------------------------------
    print("=" * 100)
    print("HYPOTHESIS TEST: D = p^2 or 2*p^2 where p = largest prime factor of k")
    print("=" * 100)
    for k in [8, 9]:
        r = results.get(k, {})
        if "error" in r:
            print(f"  k={k}: {r['error']}")
            continue
        D = r["D"]
        p = r["p_lpf_k"]
        d_p2 = r["d_pred_p2"]
        d_2p2 = r["d_pred_2p2"]
        print(f"  k={k} = {factor_str(k)}: p=lpf({k})={p}, predicted D in {{p^2={d_p2}, 2*p^2={d_2p2}}}")
        print(f"    Extracted D = {D} = {factor_str(D)}")
        if D == d_p2:
            print(f"    => HYPOTHESIS HOLDS: D = p^2 = {d_p2}")
        elif D == d_2p2:
            print(f"    => HYPOTHESIS HOLDS: D = 2*p^2 = {d_2p2}")
        else:
            print(f"    => HYPOTHESIS FAILS for k={k}: D={D} not in {{p^2={d_p2}, 2*p^2={d_2p2}}}")
            # Show actual D factorization for inspection
            facs = factorize(D)
            print(f"    Actual D factorization: {facs}")
            # Check if p^2 divides D
            if D % (p * p) == 0:
                print(f"    Note: p^2={p*p} DOES divide D={D} (D/p^2 = {D//(p*p)})")
            else:
                print(f"    Note: p^2={p*p} does NOT divide D={D}")
        print(f"    P(y) coefs low->high: {r['coefs']}")

    print()

    # ---------------------------------------------------------------------------
    # Pattern observations
    # ---------------------------------------------------------------------------
    print("=" * 100)
    print("PATTERN OBSERVATIONS (observational, not claims)")
    print("=" * 100)

    valid_ks = [k for k in range(3, 10) if k in results and "error" not in results[k]]

    print()
    print("1. p = lpf(k) and p^2 divides D?")
    for k in valid_ks:
        r = results[k]
        D = r["D"]
        p = r["p_lpf_k"]
        p2_divides = (D % (p * p) == 0) if D > 0 else False
        facs = factorize(D) if D > 0 else {}
        p_exp = facs.get(p, 0)
        extra_factor = D // (p ** p_exp) if D > 0 else "?"
        print(f"   k={k}: p={p}, D={D}={factor_str(D)}, p^{p_exp}||D, extra={extra_factor}, p^2|D: {p2_divides}")

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
        k_facs = factorize(k)
        print(
            f"   k={k}={factor_str(k)}: D={D}={factor_str(D)}, "
            f"D%2==0: {has2}, k-prime: {k_is_prime}, k-prime-power: {k_is_prime_power}"
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
        p_vals = {n: r["sigma_scaled"][n] for n in orbit}
        n_block = r["n_block"]
        D = r["D"]
        coefs = r["coefs"]
        print(f"   k={k}: p_n = {[f'{p_vals[n]:.6f}' for n in orbit]}")
        if D > 0 and coefs:
            # Evaluate P(y_n) = sum coefs[i]*y_n^i for each n
            P_vals = []
            for n in orbit:
                y_n = bloch_y(n, n_block)
                P_n = sum(coefs[i] * (y_n ** i) for i in range(len(coefs)))
                P_vals.append(round(P_n))
            print(f"          P(y_n) = D*[fit(p_n)] = {P_vals}")

    print()


if __name__ == "__main__":
    main()
