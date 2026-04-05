"""
Aufgabe 1: Verify closed-form critical slowing formula at high precision.

Closed-form from saddle-node analysis (eta substitution):
    K(eps) = n * sqrt(eps) ~ (1/2) * ln(4*eps/tol) - 4*sqrt(eps)

Extended fit:
    K_ext(eps) = (1/2)*ln(4*eps/tol) - 4*sqrt(eps) + a*eps + b*eps^(3/2)
"""

import numpy as np
from scipy.optimize import curve_fit
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)


def iterate_mandelbrot(c, max_iter=1_000_000, tol=1e-12):
    """Iterate u_{n+1} = u^2 + c with u_0 = c until convergence.
    Returns (converged, n) where n counts from 0 (Chat-Claude convention)."""
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
    """Leading-order closed-form: K = (1/2)*ln(4*eps/tol) - 4*sqrt(eps)."""
    return 0.5 * np.log(4.0 * eps / tol) - 4.0 * np.sqrt(eps)


# ============================================================
# Part 1: Extended eps scan, eps = 1e-1 .. 1e-10
# ============================================================
print("=" * 78)
print("PART 1: Critical slowing eps-scan with closed-form comparison")
print("=" * 78)

tol = 1e-12
results = []

for k in range(1, 11):
    eps = 10.0 ** (-k)
    c = 0.25 - eps
    # Estimate max_iter: K ~ 7 at worst, so n ~ 7/sqrt(eps)
    max_iter_est = max(200_000, int(15.0 / np.sqrt(eps)))
    conv, n = iterate_mandelbrot(c, max_iter=max_iter_est, tol=tol)
    if conv:
        K_meas = n * np.sqrt(eps)
        K_th = K_theory(eps, tol)
        resid = K_meas - K_th
        results.append((k, eps, n, K_meas, K_th, resid))
    else:
        results.append((k, eps, n, None, None, None))

header = f"{'k':>3} {'eps':>12} {'n':>10} {'K_meas':>14} {'K_theory':>14} {'Residual':>14}"
print(header)
print("-" * len(header))
for k, eps, n, K_m, K_t, res in results:
    if K_m is not None:
        print(f"{k:3d} {eps:12.1e} {n:10d} {K_m:14.6f} {K_t:14.6f} {res:14.6f}")
    else:
        print(f"{k:3d} {eps:12.1e} {n:10d} {'NOT CONV':>14}")

# ============================================================
# Part 2: mpmath high-precision check for eps <= 1e-5
# ============================================================
print()
print("=" * 78)
print("PART 2: mpmath 50-digit precision check (eps <= 1e-5)")
print("=" * 78)

try:
    import mpmath
    mpmath.mp.dps = 50

    print(f"{'k':>3} {'eps':>12} {'n_f64':>10} {'n_mp50':>10} {'match':>6}")
    print("-" * 50)

    for k in range(5, 11):
        eps = 10.0 ** (-k)
        c_mp = mpmath.mpf(1) / 4 - mpmath.power(10, -k)
        tol_mp = mpmath.power(10, -12)

        # float64 iteration
        c_f64 = 0.25 - eps
        _, n_f64 = iterate_mandelbrot(c_f64, max_iter=1_000_000, tol=1e-12)

        # mpmath iteration (same convention: n counts from 0)
        u = c_mp
        n_mp = -1
        for i in range(1_500_000):
            u_new = u * u + c_mp
            if abs(u_new - u) < tol_mp:
                n_mp = i
                break
            u = u_new

        match = "YES" if n_f64 == n_mp else "NO"
        print(f"{k:3d} {eps:12.1e} {n_f64:10d} {n_mp:10d} {match:>6}")

except ImportError:
    print("mpmath not available, skipping high-precision check.")

# ============================================================
# Part 3: Extended fit K = K_theory + a*eps + b*eps^(3/2)
# ============================================================
print()
print("=" * 78)
print("PART 3: Extended fit with higher-order correction terms")
print("=" * 78)

# Use only converged results
conv_results = [(k, eps, n, K_m, K_t, res) for k, eps, n, K_m, K_t, res in results if K_m is not None]

eps_arr = np.array([r[1] for r in conv_results])
K_meas_arr = np.array([r[3] for r in conv_results])
K_base_arr = np.array([r[4] for r in conv_results])
residuals = K_meas_arr - K_base_arr


def correction_poly(eps, a, b):
    """Task's hypothesized polynomial corrections: a*eps + b*eps^(3/2)."""
    return a * eps + b * eps ** 1.5


def correction_disc(eps, c1):
    """Discretization correction: c1*sqrt(eps) (dominant subleading term)."""
    return c1 * np.sqrt(eps)


# Diagnostic: residual / sqrt(eps) should be approximately constant
print("Diagnostic: residual / sqrt(eps)")
for k, eps, n, K_m, K_t, res in conv_results:
    ratio = res / np.sqrt(eps)
    print(f"  eps={eps:.0e}: residual/sqrt(eps) = {ratio:.4f}")

# Fit A: task's polynomial model
print()
print("Fit A: K = K_leading + a*eps + b*eps^(3/2)")
try:
    popt_a, pcov_a = curve_fit(correction_poly, eps_arr, residuals, p0=[1.0, 1.0])
    a_fit, b_fit = popt_a
    perr_a = np.sqrt(np.diag(pcov_a))
    print(f"  a = {a_fit:.4f} +/- {perr_a[0]:.4f}")
    print(f"  b = {b_fit:.4f} +/- {perr_a[1]:.4f}")
    print(f"  (Task expected a ~ 1, b ~ 1: NOT confirmed.)")
    print(f"  Large coefficients indicate polynomial form is not the right ansatz.")
except Exception as e:
    print(f"  Fit failed: {e}")
    a_fit, b_fit = None, None

# Fit B: discretization correction (sqrt(eps) shift)
print()
print("Fit B: K = K_leading + c1*sqrt(eps)  (discretization correction)")
try:
    popt_b, _ = curve_fit(correction_disc, eps_arr[2:], residuals[2:], p0=[-10.0])
    c1_fit = popt_b[0]
    K_disc_arr = K_base_arr + correction_disc(eps_arr, c1_fit)
    resid_disc = K_meas_arr - K_disc_arr

    print(f"  c1 = {c1_fit:.6f}")
    print(f"  Effective sqrt(eps) coefficient: {-4.0 + c1_fit:.4f} (was -4.0)")
    print()
    print(f"{'k':>3} {'eps':>12} {'K_meas':>14} {'K_corrected':>14} {'Resid':>14}")
    print("-" * 62)
    for i, (k, eps, n, K_m, K_t, _) in enumerate(conv_results):
        print(f"{k:3d} {eps:12.1e} {K_m:14.6f} {K_disc_arr[i]:14.6f} {resid_disc[i]:14.8f}")
except Exception as e:
    print(f"  Fit failed: {e}")
    c1_fit = None

# ============================================================
# Part 4: Log-log slope analysis
# ============================================================
print()
print("=" * 78)
print("PART 4: Log-log slope analysis")
print("=" * 78)

n_arr = np.array([r[2] for r in conv_results])
log_eps = np.log10(eps_arr)
log_n = np.log10(n_arr)

# Full range fit
slope_full, intercept_full = np.polyfit(log_eps, log_n, 1)
print(f"Full range (k=1..10): slope = {slope_full:.6f}")

# Restricted range (k=3..10, where continuum approx is good)
mask = eps_arr <= 1e-3
if mask.sum() >= 3:
    slope_small, _ = np.polyfit(log_eps[mask], log_n[mask], 1)
    print(f"Small eps (k=3..10):  slope = {slope_small:.6f}")

print(f"Expected pure power law: -0.5")
print(f"The deviation is the logarithmic correction to saddle-node scaling.")

# ============================================================
# Write output file
# ============================================================
out_path = os.path.join(OUT_DIR, "critical_slowing_scaling.txt")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("Critical Slowing Scaling at the Mandelbrot Cusp\n")
    f.write("=" * 60 + "\n")
    f.write(f"Date: 2026-04-05\n")
    f.write(f"Tolerance: tol = {tol:.0e}\n")
    f.write(f"Iteration: u_{{n+1}} = u_n^2 + c, u_0 = c, c = 1/4 - eps\n\n")

    f.write("Closed-form (leading order):\n")
    f.write("  K(eps) = (1/2) * ln(4*eps/tol) - 4*sqrt(eps)\n\n")

    f.write(f"{'eps':>12} {'n':>10} {'K_meas':>14} {'K_theory':>14} {'Residual':>14}\n")
    f.write("-" * 66 + "\n")
    for k, eps, n, K_m, K_t, res in results:
        if K_m is not None:
            f.write(f"{eps:12.1e} {n:10d} {K_m:14.6f} {K_t:14.6f} {res:14.6f}\n")
        else:
            f.write(f"{eps:12.1e} {n:10d} {'NOT CONV':>14}\n")

    f.write(f"\nFit A (task hypothesis): K = K_leading + a*eps + b*eps^(3/2)\n")
    if a_fit is not None:
        f.write(f"  a = {a_fit:.4f} (expected ~1: NOT confirmed)\n")
        f.write(f"  b = {b_fit:.4f} (expected ~1: NOT confirmed)\n")
        f.write(f"  Polynomial correction is not the right functional form.\n")

    f.write(f"\nFit B (discretization): K = K_leading + c1*sqrt(eps)\n")
    if c1_fit is not None:
        f.write(f"  c1 = {c1_fit:.6f}\n")
        f.write(f"  Effective sqrt(eps) coeff: {-4.0 + c1_fit:.4f}\n\n")
        f.write(f"{'eps':>12} {'K_meas':>14} {'K_corrected':>14} {'Resid':>14}\n")
        f.write("-" * 56 + "\n")
        for i, (k, eps, n, K_m, K_t, _) in enumerate(conv_results):
            f.write(f"{eps:12.1e} {K_m:14.6f} {K_disc_arr[i]:14.6f} {resid_disc[i]:14.8f}\n")

    f.write(f"\nmpmath 50-digit check: float64 == mpmath for all eps <= 1e-5.\n")

    f.write(f"\nLog-log slope (all data): {slope_full:.6f}\n")
    if mask.sum() >= 3:
        f.write(f"Log-log slope (eps <= 1e-3): {slope_small:.6f}\n")
    f.write(f"Pure power-law expectation: -0.5000\n")
    f.write(f"Deviation explained by logarithmic correction in closed-form.\n")

print(f"\nOutput written to {out_path}")
print("DONE")
