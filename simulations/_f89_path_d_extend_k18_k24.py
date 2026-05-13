"""F89 path-D closed-form belt-and-suspenders: k=18..24 (2026-05-13).

Tests the candidate formula on six additional sequential data points (k=18..23)
and a CRITICAL second v2=3 test at k=24.

Candidate formula (verified k=3..17, 15 data points):
  D_k = (odd(k))^2 * 2^E(k)
  E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)

Equivalently: E(k) = max(0, floor((k-5)/2)) + max(v2(k), 2*v2(k) - 2).

Predictions:
  k=18: v2=1, odd=9,  E = 6+1+0 = 7,  D = 81  * 128   = 10368
  k=19: v2=0, odd=19, E = 7+0+0 = 7,  D = 361 * 128   = 46208
  k=20: v2=2, odd=5,  E = 7+2+0 = 9,  D = 25  * 512   = 12800
  k=21: v2=0, odd=21, E = 8+0+0 = 8,  D = 441 * 256   = 112896
  k=22: v2=1, odd=11, E = 8+1+0 = 9,  D = 121 * 512   = 61952
  k=23: v2=0, odd=23, E = 9+0+0 = 9,  D = 529 * 512   = 270848
  k=24: v2=3, odd=3,  E = 9+3+1 = 13, D = 9   * 8192  = 73728   (CRITICAL: v2=3 bonus again)

k=24 is the second v2=3 data point (first: k=8). If actual D_24=73728, the
max(0, v2(k)-2) bonus is confirmed at TWO independent v2=3 points.

Block dimensions for (SE,DE) sub-block:
  k=18: nBlock=19, dim = 19 * C(19,2) = 19*171 = 3249
  k=19: nBlock=20, dim = 20 * C(20,2) = 20*190 = 3800
  k=20: nBlock=21, dim = 21 * C(21,2) = 21*210 = 4410
  k=21: nBlock=22, dim = 22 * C(22,2) = 22*231 = 5082
  k=22: nBlock=23, dim = 23 * C(23,2) = 23*253 = 5819
  k=23: nBlock=24, dim = 24 * C(24,2) = 24*276 = 6624
  k=24: nBlock=25, dim = 25 * C(25,2) = 25*300 = 7500

Method (identical to _f89_path_d_verify_k16_k17.py which validated k=16,17):
  1. Build (SE,DE) sub-block L via build_se_de_L.
  2. Eigendecompose, filter F_a modes by Re(eig) near -2*gamma.
  3. Match modes to SE-anti Bloch orbit (even n, n=2..nBlock step 2).
  4. Compute sigma_n = |c_n|^2 * sum_l |Mv_n[l]|^2, scaled by N^2*(N-1).
  5. Vandermonde-fit polynomial P(y); find smallest D making all coefs integral.
  6. Compare actual D vs predicted D from candidate formula.

Output saved to simulations/results/_f89_path_d_extend_k18_k24.txt.
"""
from __future__ import annotations

import sys
import math
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# 2-adic helpers
# ---------------------------------------------------------------------------

def v2(n: int) -> int:
    """2-adic valuation of n."""
    if n <= 0:
        return 0
    count = 0
    while n % 2 == 0:
        n //= 2
        count += 1
    return count


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
    """Compact factorization: '2^7*3^4'."""
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
# (SE,DE) sub-block Liouvillian (verbatim from _f89_path_d_verify_k16_k17.py)
# ---------------------------------------------------------------------------

def build_se_de_L(n_block: int, J_val: float, gamma_val: float):
    """Build (SE,DE) sub-block Liouvillian for path-(n_block-1).

    Returns (L, basis, de_pairs).
    """
    de_pairs = [(j, kk) for j in range(n_block) for kk in range(j + 1, n_block)]
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
# Parameters (match prior scripts exactly)
# ---------------------------------------------------------------------------
J = 1.0
GAMMA = 0.05
N_FIXED = 9       # base N; raised to n_block + 2 if n_block + 2 > N_FIXED

FA_TOL_REL = 0.01
BLOCH_TOL = 1e-3
D_MAX = 5_000_000
RATIONAL_TOL = 1e-3


# ---------------------------------------------------------------------------
# Sigma extraction (SE-anti orbit)
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

    dim = n_block * (n_block * (n_block - 1) // 2)
    print(f"    building (SE,DE) sub-block (dim={dim})...", flush=True)
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
# Per-path run function
# ---------------------------------------------------------------------------

def run_path(k: int) -> tuple[str, int | None, int, int, list[int], float]:
    """Run extraction for a single path k.

    Returns (verdict_symbol, actual_d, pred_d, E, coefs, wall_time_s).
    verdict_symbol: '+' = confirmed, 'X' = broken, 'ERROR' / 'NOFOUND' = failure.
    """
    n_block = k + 1
    orbit = se_anti_orbit(n_block)
    pred_d_val, e_val = predicted_d(k)
    vk = v2(k)

    print(f"\n  --- path k={k} (nBlock={n_block}, orbit size={len(orbit)}) ---", flush=True)
    print(
        f"  Predicted: E = max(0,{(k-5)//2}) + {vk} + max(0,{vk}-2)"
        f" = {max(0,(k-5)//2)} + {vk} + {max(0,vk-2)} = {e_val}",
        flush=True,
    )
    print(
        f"  Predicted D = (odd({k}))^2 * 2^{e_val}"
        f" = {odd_part(k)}^2 * {2**e_val} = {pred_d_val} = {factor_str(pred_d_val)}",
        flush=True,
    )

    t_start = time.perf_counter()
    sigma_scaled = extract_sigma_scaled(n_block)
    nan_ns = [n for n, vv in sigma_scaled.items() if math.isnan(vv)]
    if nan_ns:
        print(f"  ERROR: NaN sigma at orbit indices {nan_ns}", flush=True)
        return "ERROR", None, pred_d_val, e_val, [], time.perf_counter() - t_start

    p_vals = [sigma_scaled[n] for n in orbit]
    y_vals = [bloch_y(n, n_block) for n in orbit]
    print(f"  p_n (scaled sigmas): {[f'{vv:.6f}' for vv in p_vals]}", flush=True)

    print(f"  Rationalizing polynomial (D_max={D_MAX})...", flush=True)
    actual_d, coefs = rationalize_polynomial_D(p_vals, y_vals, D_MAX)
    wall_time = time.perf_counter() - t_start

    if actual_d < 0:
        print(f"  ERROR: No D <= {D_MAX} found.", flush=True)
        return "NOFOUND", None, pred_d_val, e_val, [], wall_time

    match = (actual_d == pred_d_val)
    symbol = "+" if match else "X"

    print(f"  Extracted D = {actual_d} = {factor_str(actual_d)}", flush=True)
    print(f"  Polynomial coefs (low->high): {coefs}", flush=True)
    print(f"  Wall time: {wall_time:.1f}s", flush=True)
    if match:
        print(
            f"  Result: + FORMULA CONFIRMED (predicted {pred_d_val} [{factor_str(pred_d_val)}],"
            f" actual {actual_d})",
            flush=True,
        )
    else:
        print(
            f"  Result: X FORMULA BROKEN"
            f" -- predicted {pred_d_val} [{factor_str(pred_d_val)}],"
            f" actual {actual_d} [{factor_str(actual_d)}]",
            flush=True,
        )
        print(f"  Deviation: actual/predicted = {actual_d}/{pred_d_val}", flush=True)
        print(f"  v2(actual)={v2(actual_d)}, v2(predicted)={v2(pred_d_val)}", flush=True)
        print(
            f"  odd(actual)={odd_part(actual_d)}, odd(predicted)={odd_part(pred_d_val)}",
            flush=True,
        )

    return symbol, actual_d, pred_d_val, e_val, coefs, wall_time


# ---------------------------------------------------------------------------
# Target k values and pre-computed predictions table
# ---------------------------------------------------------------------------

TARGETS = [18, 19, 20, 21, 22, 23, 24]

# Pre-computed expected values (for documentation; verified by predicted_d() at runtime)
_EXPECTED = {
    18: (7,  10368,  "v2=1, max(0,1-2)=0"),
    19: (7,  46208,  "v2=0, bonus=0"),
    20: (9,  12800,  "v2=2, max(0,2-2)=0"),
    21: (8,  112896, "v2=0, bonus=0"),
    22: (9,  61952,  "v2=1, bonus=0"),
    23: (9,  270848, "v2=0, bonus=0"),
    24: (13, 73728,  "v2=3, max(0,3-2)=1 CRITICAL"),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    header = "F89 path-D closed-form belt-and-suspenders: k=18..24"
    print("=" * 80)
    print(header)
    print("=" * 80)
    print()
    print("Candidate formula (verified k=3..17, 15 data points):")
    print("  D_k = (odd(k))^2 * 2^E(k)")
    print("  E(k) = max(0, floor((k-5)/2)) + v2(k) + max(0, v2(k)-2)")
    print()
    print("Predictions:")
    print(f"  {'k':>4}  {'v2(k)':>5}  {'odd(k)':>6}  {'E(k)':>5}  {'Predicted D':>12}  note")
    print("  " + "-" * 65)
    for k in TARGETS:
        d_pred, e = predicted_d(k)
        vk = v2(k)
        note = _EXPECTED[k][2]
        print(f"  {k:>4}  {vk:>5}  {odd_part(k):>6}  {e:>5}  {d_pred:>12}  {note}")
    print()
    print("k=24 is the CRITICAL test: second v2=3 data point (first was k=8).")
    print("  If actual D_24 = 73728 = 2^13*9, the max(0,v2-2) bonus is confirmed")
    print("  at TWO independent v2=3 data points.")
    print()
    print(f"Parameters: J={J}, gamma={GAMMA}, N_fixed={N_FIXED} (raised if needed)")
    print(f"  D_max={D_MAX}, rationalization tol={RATIONAL_TOL}, Bloch tol={BLOCH_TOL}")
    print()

    results: dict[int, tuple[str, int | None, int, int, list[int], float]] = {}

    for k in TARGETS:
        print("=" * 80)
        label = "(CRITICAL: v2=3 bonus second data point)" if k == 24 else ""
        print(f"Running k={k} {label}")
        print("=" * 80)
        results[k] = run_path(k)
        sym = results[k][0]
        wt = results[k][5]
        print(f"\nk={k} DONE ({wt:.1f}s). Verdict: {sym}", flush=True)
        print()

    # ---------------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------------
    print()
    print("=" * 80)
    print("SUMMARY: F89 path-D closed-form verification k=18..24")
    print("=" * 80)
    print()
    col_fmt = f"  {'k':>4}  {'nBlock':>6}  {'E':>4}  {'Predicted D':>12}  {'Pred factored':>16}  {'Actual D':>12}  {'Actual factored':>16}  {'Wall(s)':>7}  verdict"
    print(col_fmt)
    print("  " + "-" * 105)

    all_confirmed = True
    broken_at: list[tuple[int, int, int]] = []

    for k in TARGETS:
        sym, actual_d, pred_d_val, e_val, _coefs, wt = results[k]
        n_block = k + 1
        pred_f = factor_str(pred_d_val)
        actual_s = str(actual_d) if actual_d is not None else "N/A"
        actual_f = factor_str(actual_d) if actual_d is not None else "N/A"
        if sym == "+":
            verdict = "+ confirmed"
        elif sym == "X":
            verdict = "X BROKEN"
            all_confirmed = False
            broken_at.append((k, pred_d_val, actual_d or 0))
        else:
            verdict = f"  {sym}"
            all_confirmed = False
        print(
            f"  {k:>4}  {n_block:>6}  {e_val:>4}  {pred_d_val:>12}  {pred_f:>16}"
            f"  {actual_s:>12}  {actual_f:>16}  {wt:>7.1f}  {verdict}"
        )

    print()

    # Special highlight for k=24 v2=3 bonus
    sym24, actual24, pred24, _, _, _ = results[24]
    v2_bonus_confirmed = (sym24 == "+")
    print("Critical k=24 (v2=3 bonus) result:")
    if v2_bonus_confirmed:
        print(f"  + CONFIRMED: actual D={actual24} = {factor_str(actual24 or 0)}")
        print("  The max(0, v2(k)-2) bonus is now confirmed at TWO v2=3 data points:")
        print("    k=8  (v2=3, bonus=1, D=32=2^5)")
        print("    k=24 (v2=3, bonus=1, D=73728=2^13*9)")
    else:
        print(f"  X NOT CONFIRMED: predicted={pred24}, actual={actual24}")
        if actual24 is not None:
            print(f"  v2(actual)={v2(actual24)}, odd(actual)={odd_part(actual24)}")
            if odd_part(actual24 or 0) == odd_part(pred24):
                delta_v2 = v2(actual24 or 0) - (v2(24) + max(0, (24-5)//2))
                print(f"  Odd-part rule holds; v2(D)={v2(actual24)}, need bonus={delta_v2} at v2(k)=3")

    print()

    # Overall verdict line
    n_data_points = 17 + 7  # k=3..17 prior + k=18..24 new
    if all_confirmed:
        print(
            f"Candidate formula confirmed at all {len(TARGETS)} additional data points"
            f" (k=18..24)."
        )
        print(
            f"Formula now validated on {n_data_points} data points total (k=3..24)."
        )
    else:
        fail_strs = [f"k={k} (predicted {p}, actual {a})" for k, p, a in broken_at]
        print(f"Formula broken at: {'; '.join(fail_strs)}")

    # Wall-time summary
    total_wt = sum(results[k][5] for k in TARGETS)
    print()
    print("Wall-time per path:")
    for k in TARGETS:
        wt = results[k][5]
        print(f"  k={k}: {wt:.1f}s")
    print(f"  Total: {total_wt:.1f}s")


if __name__ == "__main__":
    # Tee output to file
    out_dir = Path(__file__).parent / "results"
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / "_f89_path_d_extend_k18_k24.txt"

    class Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, data: str) -> int:
            for s in self._streams:
                s.write(data)
            return len(data)

        def flush(self) -> None:
            for s in self._streams:
                s.flush()

        @property
        def encoding(self):
            return self._streams[0].encoding

        @property
        def errors(self):
            return self._streams[0].errors

    with open(out_file, "w", encoding="utf-8") as fh:
        tee = Tee(sys.stdout, fh)
        sys.stdout = tee  # type: ignore[assignment]
        try:
            main()
        finally:
            sys.stdout = tee._streams[0]  # type: ignore[assignment]

    print(f"\nOutput written to: {out_file}", flush=True)
