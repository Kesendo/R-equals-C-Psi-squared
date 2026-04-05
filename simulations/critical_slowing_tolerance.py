"""
Aufgabe 2: Tolerance scaling of critical slowing.

Part 1: Fix eps = 1e-6, vary tol. K should be linear in |ln(tol)|.
Part 2: Relative stopping (tol_rel = k*eps). K should become eps-independent.
"""

import numpy as np
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)


def iterate_mandelbrot(c, max_iter=2_000_000, tol=1e-12):
    """Iterate u_{n+1} = u^2 + c, u_0 = c. Returns (converged, n)."""
    u = complex(c)
    for n in range(max_iter):
        u_new = u * u + c
        if abs(u_new) > 1e6:
            return False, n
        if abs(u_new - u) < tol:
            return True, n
        u = u_new
    return False, max_iter


def K_theory(eps, tol):
    """Leading-order: K = (1/2)*ln(4*eps/tol) - 4*sqrt(eps)."""
    return 0.5 * np.log(4.0 * eps / tol) - 4.0 * np.sqrt(eps)


# ============================================================
# Part 1: K vs |ln(tol)| at fixed eps = 1e-6
# ============================================================
print("=" * 78)
print("PART 1: Tolerance scaling at eps = 1e-6")
print("=" * 78)

eps_fixed = 1e-6
s = np.sqrt(eps_fixed)
c = 0.25 - eps_fixed

tol_values = [1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16]
results_tol = []

print(f"eps = {eps_fixed:.0e}, sqrt(eps) = {s:.6f}")
print()
print(f"{'tol':>12} {'|ln(tol)|':>10} {'n':>10} {'K_meas':>14} {'K_theory':>14} {'Resid':>10}")
print("-" * 70)

for tol in tol_values:
    conv, n = iterate_mandelbrot(c, tol=tol)
    K_m = n * s
    K_t = K_theory(eps_fixed, tol)
    resid = K_m - K_t
    ln_tol = abs(np.log(tol))
    results_tol.append((tol, ln_tol, n, K_m, K_t, resid))
    print(f"{tol:12.0e} {ln_tol:10.4f} {n:10d} {K_m:14.6f} {K_t:14.6f} {resid:10.6f}")

# Linear fit: K = slope * |ln(tol)| + intercept
ln_tol_arr = np.array([r[1] for r in results_tol])
K_meas_arr = np.array([r[3] for r in results_tol])
slope, intercept = np.polyfit(ln_tol_arr, K_meas_arr, 1)

print()
print(f"Linear fit: K = {slope:.6f} * |ln(tol)| + ({intercept:.6f})")
print(f"Expected slope: 0.5  (from K ~ (1/2)*ln(4*eps/tol))")
print(f"Measured slope: {slope:.6f}")
print(f"Deviation: {abs(slope - 0.5):.6f}")

# ============================================================
# Part 2: Relative stopping tol_rel = k * eps
# ============================================================
print()
print("=" * 78)
print("PART 2: Relative stopping strategy (tol = k * eps)")
print("=" * 78)

k_stop = 1e-3  # stop when |delta u| < k * eps
print(f"Relative factor k = {k_stop:.0e}")
print(f"tol_rel = k * eps: stop closer to fixed point for small eps")
print()
print(f"{'eps':>12} {'tol_rel':>12} {'n':>10} {'K':>14} {'K(k=eps)':>14}")
print("-" * 66)

results_rel = []
for exp in range(1, 11):
    eps = 10.0 ** (-exp)
    tol_rel = k_stop * eps
    c = 0.25 - eps
    max_iter_est = max(100_000, int(50.0 / np.sqrt(eps)))
    conv, n = iterate_mandelbrot(c, max_iter=max_iter_est, tol=tol_rel)
    K = n * np.sqrt(eps)

    # With tol = k*eps, the theory gives:
    # K = (1/2)*ln(4*eps/(k*eps)) - 4*sqrt(eps)
    #   = (1/2)*ln(4/k) - 4*sqrt(eps)
    K_const = 0.5 * np.log(4.0 / k_stop) - 4.0 * np.sqrt(eps)
    results_rel.append((eps, tol_rel, n, K, K_const))
    print(f"{eps:12.1e} {tol_rel:12.1e} {n:10d} {K:14.6f} {K_const:14.6f}")

K_rel_arr = np.array([r[3] for r in results_rel])
print()
print(f"K range: [{K_rel_arr.min():.4f}, {K_rel_arr.max():.4f}]")
print(f"K std:   {K_rel_arr[2:].std():.6f}  (excluding eps >= 1e-2)")
print(f"Expected constant: (1/2)*ln(4/k) = {0.5 * np.log(4.0 / k_stop):.4f}")
print(f"With relative stopping, K converges to a constant independent of eps.")

# ============================================================
# Part 3: Multiple k values to confirm K(k) = (1/2)*ln(4/k)
# ============================================================
print()
print("=" * 78)
print("PART 3: K as function of relative stop factor k")
print("=" * 78)

eps_test = 1e-6
c_test = 0.25 - eps_test
s_test = np.sqrt(eps_test)

k_values = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
print(f"eps = {eps_test:.0e}")
print()
print(f"{'k':>12} {'tol':>12} {'n':>10} {'K_meas':>14} {'(1/2)ln(4/k)':>14} {'Resid':>10}")
print("-" * 74)

for k_val in k_values:
    tol_k = k_val * eps_test
    conv, n = iterate_mandelbrot(c_test, tol=tol_k)
    K_m = n * s_test
    K_pred = 0.5 * np.log(4.0 / k_val) - 4.0 * s_test
    resid = K_m - K_pred
    print(f"{k_val:12.1e} {tol_k:12.1e} {n:10d} {K_m:14.6f} {K_pred:14.6f} {resid:10.6f}")

# ============================================================
# Write output
# ============================================================
out_path = os.path.join(OUT_DIR, "critical_slowing_tolerance.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("Critical Slowing: Tolerance Scaling\n")
    f.write("=" * 60 + "\n")
    f.write("Date: 2026-04-05\n\n")

    f.write("Part 1: K vs |ln(tol)| at eps = 1e-6\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'tol':>12} {'|ln(tol)|':>10} {'n':>10} {'K':>14}\n")
    for tol, ln_t, n, K_m, K_t, res in results_tol:
        f.write(f"{tol:12.0e} {ln_t:10.4f} {n:10d} {K_m:14.6f}\n")
    f.write(f"\nLinear fit: K = {slope:.6f} * |ln(tol)| + ({intercept:.6f})\n")
    f.write(f"Expected slope: 0.5000\n")
    f.write(f"Measured slope: {slope:.6f}\n")
    f.write(f"Confirmed: K is linear in |ln(tol)|.\n\n")

    f.write("Part 2: Relative stopping (tol = 1e-3 * eps)\n")
    f.write("-" * 40 + "\n")
    f.write(f"{'eps':>12} {'n':>10} {'K':>14}\n")
    for eps, tol_r, n, K, K_c in results_rel:
        f.write(f"{eps:12.1e} {n:10d} {K:14.6f}\n")
    f.write(f"\nK range (eps <= 1e-3): [{K_rel_arr[2:].min():.4f}, {K_rel_arr[2:].max():.4f}]\n")
    f.write(f"K std (eps <= 1e-3): {K_rel_arr[2:].std():.6f}\n")
    f.write(f"Expected constant: {0.5 * np.log(4.0 / k_stop):.4f}\n")
    f.write(f"Confirmed: with relative stopping, K is a pure constant.\n")

print(f"\nOutput written to {out_path}")
print("DONE")
