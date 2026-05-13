"""F89 path-D closed-form verification at k=16 and k=17 (critical test, 2026-05-13).

Tests the candidate formula extracted from k=3..15 data:
  D_k = (odd(k))^2 * 2^E(k)
  E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)

Predictions:
  k=16: v2=4, odd=1, E=5+4+2=11, D = 1^2 * 2^11 = 2048
         (critical test for max(0, v2(k)-2) bonus at v2(k)=4)
  k=17: v2=0, odd=17, E=6+0+0=6, D = 17^2 * 2^6 = 18496

Block dimensions for (SE,DE) sub-block:
  k=16: nBlock=17, dim = 17 * C(17,2) = 17*136 = 2312
  k=17: nBlock=18, dim = 18 * C(18,2) = 18*153 = 2754

Method (identical to _f89_path_d_structure_probe.py which validated k=3..15):
  1. Build (SE,DE) sub-block L via build_se_de_L.
  2. Eigendecompose, filter F_a modes by Re(eig) near -2*gamma.
  3. Match modes to SE-anti Bloch orbit (even n, n=2..nBlock step 2).
  4. Compute sigma_n = |c_n|^2 * sum_l |Mv_n[l]|^2, scaled by N^2*(N-1).
  5. Vandermonde-fit polynomial P(y); find smallest D making all coefs integral.
  6. Compare actual D vs predicted D from candidate formula.

k=16 runs FIRST (early-exit on failure to guide refinement if formula is broken).
"""
from __future__ import annotations

import sys
import math
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

# Ensure UTF-8 output on Windows (checkmarks, etc.)
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Re-use the build_se_de_L + per_site_reduction_se_de from the probe script
# (inline them here to keep this script fully self-contained and avoid import ambiguity)

# ---------------------------------------------------------------------------
# 2-adic helpers
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
    d = 2
    tmp = abs(n)
    while d * d <= tmp:
        while tmp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            tmp //= d
        d += 1
    if tmp > 1:
        factors[tmp] = factors.get(tmp, 0) + 1
    return factors


def factor_str(n: int) -> str:
    """Compact factorization: '2^5*3^2'."""
    if n == 0:
        return "0"
    if n == 1:
        return "1"
    facs = factorize(n)
    parts = []
    for p in sorted(facs):
        e = facs[p]
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts)


def predicted_d(k: int) -> tuple[int, int]:
    """Candidate formula: D_k = (odd(k))^2 * 2^E(k).

    E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)

    Returns (D_predicted, E).
    """
    vk = v2(k)
    e = max(0, (k - 5) // 2) + vk + max(0, vk - 2)
    return (odd_part(k) ** 2) * (2 ** e), e


# ---------------------------------------------------------------------------
# (SE,DE) sub-block Liouvillian (copied verbatim from _f89_path_d_structure_probe.py)
# ---------------------------------------------------------------------------

def build_se_de_L(n_block: int, J_val: float, gamma_val: float):
    """Build (SE,DE) sub-block Liouvillian for path-(n_block-1).

    Returns (L, basis, de_pairs).
    """
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
    for idx, (j, kk) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < n_block and new_j != kk:
                new_pair = tuple(sorted([new_j, kk]))
                if abs(new_j - j) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2.0 * J_val
        for new_k in [kk - 1, kk + 1]:
            if 0 <= new_k < n_block and new_k != j:
                new_pair = tuple(sorted([j, new_k]))
                if abs(new_k - kk) == 1 and new_pair in de_pairs:
                    M_DE[de_pairs.index(new_pair), idx] += 2.0 * J_val

    L = np.zeros((n_basis, n_basis), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE[i2, i] != 0.0:
                idx2 = basis.index((i2, jk))
                L[idx2, idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(n_de):
            if M_DE[jk_idx, jk2_idx] != 0.0:
                jk2 = de_pairs[jk2_idx]
                idx2 = basis.index((i, jk2))
                L[idx2, idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2.0 * gamma_val if i in jk else -6.0 * gamma_val

    return L, basis, de_pairs


def per_site_reduction_se_de(basis: list, de_pairs: list, n_block: int) -> np.ndarray:
    """Per-site reduction matrix w for (SE,DE) basis."""
    w = np.zeros((n_block, len(basis)), dtype=float)
    for b_idx, (i, jk) in enumerate(basis):
        j, kk = jk
        if kk == i:
            w[j, b_idx] = 1.0
        if j == i:
            w[kk, b_idx] = 1.0
    return w


# ---------------------------------------------------------------------------
# Parameters (match _f89_path_d_structure_probe.py exactly)
# ---------------------------------------------------------------------------
J = 1.0
GAMMA = 0.05
N_FIXED = 9          # base N; raised to n_block + 2 if n_block + 2 > N_FIXED

FA_TOL_REL = 0.01    # |Re(eig)/gamma + 2| < FA_TOL_REL
BLOCH_TOL = 1e-3     # |Im(eig)/J - y_n| < BLOCH_TOL
D_MAX = 5_000_000
RATIONAL_TOL = 1e-3


# ---------------------------------------------------------------------------
# Sigma extraction (SE-anti orbit, same method as probe)
# ---------------------------------------------------------------------------

def se_anti_orbit(n_block: int) -> list[int]:
    return list(range(2, n_block + 1, 2))


def bloch_y(n: int, n_block: int) -> float:
    return 4.0 * math.cos(math.pi * n / (n_block + 1))


def extract_sigma_scaled(n_block: int) -> dict[int, float]:
    """Extract p_n = sigma_n * N^2*(N-1) for SE-anti orbit."""
    N = max(N_FIXED, n_block + 2)
    orbit = se_anti_orbit(n_block)
    y_targets = {n: bloch_y(n, n_block) for n in orbit}

    print(f"    building (SE,DE) sub-block (dim={n_block * n_block * (n_block - 1) // 2})...",
          flush=True)
    L, basis, de_pairs = build_se_de_L(n_block, J, GAMMA)

    print(f"    eigendecomposing {L.shape[0]}x{L.shape[0]} matrix...", flush=True)
    t0 = time.perf_counter()
    eigvals, eigvecs = np.linalg.eig(L)
    dt = time.perf_counter() - t0
    print(f"    eigendecomp done in {dt:.1f}s.", flush=True)

    fa_mask = np.abs(eigvals.real / GAMMA + 2.0) < FA_TOL_REL
    n_fa = int(fa_mask.sum())
    print(f"    F_a modes found: {n_fa} (expected {len(orbit)})", flush=True)

    pre = math.sqrt(2.0 / (N * N * (N - 1)))
    rho_flat = np.full(len(basis), pre / 2, dtype=complex)
    w = per_site_reduction_se_de(basis, de_pairs, n_block)

    sigma_scaled: dict[int, float] = {}
    used_idx: set[int] = set()

    for n, y_n in sorted(y_targets.items()):
        best_idx = -1
        best_dist = 1e10
        for idx in np.where(fa_mask)[0]:
            if idx in used_idx:
                continue
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
# Rationalization: Vandermonde + smallest integer D
# ---------------------------------------------------------------------------

def rationalize_polynomial_D(
    p_vals: list[float],
    y_vals: list[float],
    d_max: int = D_MAX,
) -> tuple[int, list[int]]:
    """Find smallest D in [1, d_max] making all Vandermonde coefs near-integer."""
    deg = len(y_vals) - 1
    V = np.vander(y_vals, N=deg + 1, increasing=True)
    coefs_d1, _, _, _ = np.linalg.lstsq(V, np.array(p_vals, dtype=float), rcond=None)

    for D in range(1, d_max + 1):
        scaled = D * coefs_d1
        if all(abs(s - round(s)) < RATIONAL_TOL for s in scaled):
            int_coefs = [round(float(s)) for s in scaled]
            residuals = [
                abs(sum(int_coefs[i] * (y ** i) for i in range(len(int_coefs))) - D * p)
                for y, p in zip(y_vals, p_vals)
            ]
            if max(residuals) < 1e-1:
                return D, int_coefs
    return -1, []


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_path(k: int) -> tuple[str, int | None, int, list[int]]:
    """Run extraction for a single path k.

    Returns (verdict_symbol, actual_d, pred_d, coefs).
    """
    n_block = k + 1
    orbit = se_anti_orbit(n_block)
    pred_d_val, e_val = predicted_d(k)

    print(f"\n  --- path k={k} (nBlock={n_block}, orbit={orbit}) ---", flush=True)
    vk = v2(k)
    print(f"  Predicted: E={e_val} = max(0,{(k-5)//2}) + {vk} + max(0,{vk}-2) = "
          f"{max(0,(k-5)//2)} + {vk} + {max(0,vk-2)}", flush=True)
    print(f"  Predicted D = (odd({k}))^2 * 2^{e_val} = {odd_part(k)}^2 * 2^{e_val} = {pred_d_val}",
          flush=True)

    t_start = time.perf_counter()
    sigma_scaled = extract_sigma_scaled(n_block)
    nan_ns = [n for n, vv in sigma_scaled.items() if math.isnan(vv)]
    if nan_ns:
        print(f"  ERROR: NaN sigma at orbit indices {nan_ns}", flush=True)
        return "ERROR", None, pred_d_val, []

    p_vals = [sigma_scaled[n] for n in orbit]
    y_vals = [bloch_y(n, n_block) for n in orbit]
    print(f"  p_n (scaled sigmas): {[f'{vv:.6f}' for vv in p_vals]}", flush=True)

    print(f"  Rationalizing polynomial (D_max={D_MAX})...", flush=True)
    actual_d, coefs = rationalize_polynomial_D(p_vals, y_vals, D_MAX)
    t_elapsed = time.perf_counter() - t_start

    if actual_d < 0:
        print(f"  ERROR: No D <= {D_MAX} found.", flush=True)
        return "NOFOUND", None, pred_d_val, []

    match = (actual_d == pred_d_val)
    verdict = "CONFIRMED" if match else "BROKEN"
    symbol = "+" if match else "X"

    print(f"  Extracted D = {actual_d} = {factor_str(actual_d)}", flush=True)
    print(f"  Polynomial coefs (low->high): {coefs}", flush=True)
    print(f"  Wall time: {t_elapsed:.1f}s", flush=True)
    if match:
        print(f"  Result: {symbol} FORMULA {verdict} (predicted {pred_d_val}, actual {actual_d})",
              flush=True)
    else:
        print(f"  Result: {symbol} FORMULA {verdict} -- predicted {pred_d_val} [{factor_str(pred_d_val)}],"
              f" actual {actual_d} [{factor_str(actual_d)}]", flush=True)
        print(f"  Deviation: actual/predicted = {actual_d}/{pred_d_val}", flush=True)
        print(f"  v2(actual)={v2(actual_d)}, v2(predicted)={v2(pred_d_val)}", flush=True)
        print(f"  odd(actual)={odd_part(actual_d)}, odd(predicted)={odd_part(pred_d_val)}", flush=True)

    return symbol, actual_d, pred_d_val, coefs


def main() -> None:
    print("=" * 72)
    print("F89 path-D closed-form verification: k=16 and k=17")
    print("=" * 72)
    print()
    print("Candidate formula (extracted from k=3..15 data, 13 points):")
    print("  D_k = (odd(k))^2 * 2^E(k)")
    print("  E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)")
    print()
    print("Critical test: k=16 has v2(k)=4, bonus=max(0,4-2)=2.")
    print("  This bonus is currently supported by only ONE data point (k=8, v2=3, bonus=1).")
    print("  k=16 is the first v2=4 case -- if actual v2(D) != 11, the bonus formula is wrong.")
    print()
    print("Parameters:")
    print(f"  J={J}, gamma={GAMMA}, N_fixed={N_FIXED} (raised if n_block+2 > N_fixed)")
    print(f"  D_max={D_MAX}, rationalization tol={RATIONAL_TOL}, Bloch tol={BLOCH_TOL}")
    print()

    # Pre-print predictions
    print("Predictions:")
    for k in [16, 17]:
        d_pred, e = predicted_d(k)
        vk = v2(k)
        print(f"  k={k}: E = max(0,{(k-5)//2}) + {vk} + max(0,{vk}-2)"
              f" = {max(0,(k-5)//2)} + {vk} + {max(0,vk-2)} = {e}")
        print(f"         D = {odd_part(k)}^2 * 2^{e} = {odd_part(k)**2} * {2**e} = {d_pred}")
    print()

    results: dict[int, tuple[str, int | None, int, list[int]]] = {}

    # Run k=16 FIRST (critical test)
    print("=" * 72)
    print("Running k=16 (CRITICAL TEST for v2(k)=4 bonus) ...")
    print("=" * 72)
    t_wall_k16 = time.perf_counter()
    results[16] = run_path(16)
    t_wall_k16 = time.perf_counter() - t_wall_k16

    symbol16, actual16, pred16, coefs16 = results[16]

    print()
    print(f"k=16 DONE ({t_wall_k16:.1f}s). Verdict: {symbol16}")
    if symbol16 not in ("CONFIRMED", "+"):
        print("  k=16 formula BROKEN. Proceeding to k=17 for additional data.")

    # Run k=17
    print()
    print("=" * 72)
    print("Running k=17 ...")
    print("=" * 72)
    t_wall_k17 = time.perf_counter()
    results[17] = run_path(17)
    t_wall_k17 = time.perf_counter() - t_wall_k17

    symbol17, actual17, pred17, coefs17 = results[17]

    print()
    print(f"k=17 DONE ({t_wall_k17:.1f}s). Verdict: {symbol17}")

    # ---------------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------------
    print()
    print("=" * 72)
    print("SUMMARY: F89 path-D closed-form verification k=16, k=17")
    print("=" * 72)
    print()
    print(f"{'k':>4}  {'nBlock':>6}  {'pred D':>10}  {'pred-factored':>18}  {'actual D':>10}  {'actual-factored':>18}  verdict")
    print("-" * 95)
    for k, (sym, actual_d, pred_d_val, _) in results.items():
        n_block = k + 1
        pred_f = factor_str(pred_d_val)
        actual_f = factor_str(actual_d) if actual_d is not None else "N/A"
        actual_s = str(actual_d) if actual_d is not None else "N/A"
        v = f"{sym} formula confirmed" if sym == "+" else (
            f"X formula BROKEN" if sym == "X" else f"  {sym}")
        print(f"{k:>4}  {n_block:>6}  {pred_d_val:>10}  {pred_f:>18}  {actual_s:>10}  {actual_f:>18}  {v}")

    # Overall verdict
    both_confirmed = (symbol16 == "+") and (symbol17 == "+")
    any_broken = (symbol16 == "X") or (symbol17 == "X")

    print()
    if both_confirmed:
        print("OVERALL: + FORMULA CONFIRMED at k=16 and k=17")
        print("  Candidate D_k = (odd(k))^2 * 2^E(k) now validated on 15 data points (k=3..17).")
        print("  The v2(k)=4 bonus (max(0,v2(k)-2)=2) is confirmed by k=16.")
    elif any_broken:
        print("OVERALL: X FORMULA BROKEN at one or both test points.")
        print("  See per-path output above for deviation analysis.")
        # Suggest a refined formula if odd-part rule holds but v2 is off
        if actual16 is not None and odd_part(actual16) == odd_part(pred16):
            actual_v2_16 = v2(actual16)
            pred_v2_16 = v2(pred16)
            print(f"  k=16: odd-part rule holds (odd(D)={odd_part(actual16)}=1^2=1), "
                  f"but v2(D)={actual_v2_16} vs predicted {pred_v2_16}.")
            print(f"  Refined bonus for v2(k)=4: max(0,v2(k)-2) should be {actual_v2_16 - (max(0,(16-5)//2) + v2(16))}?")
        if actual17 is not None and odd_part(actual17) == odd_part(pred17):
            actual_v2_17 = v2(actual17)
            pred_v2_17 = v2(pred17)
            print(f"  k=17: odd-part rule holds (odd(D)={odd_part(actual17)}=17^2={17**2}), "
                  f"but v2(D)={actual_v2_17} vs predicted {pred_v2_17}.")
    else:
        print("OVERALL: INCONCLUSIVE (extraction errors).")

    print()
    print(f"Wall times: k=16: {t_wall_k16:.1f}s, k=17: {t_wall_k17:.1f}s")


if __name__ == "__main__":
    main()
