"""Critical slowing scan near c = 1/4."""
import numpy as np

def iterate_mandelbrot(c, max_iter=50000, bailout=1e6, tol=1e-12):
    u = complex(c)
    for n in range(max_iter):
        u_new = u*u + c
        if abs(u_new) > bailout:
            return False, u_new, n
        if abs(u_new - u) < tol:
            return True, u_new, n
        u = u_new
    return False, u, max_iter

print()
print("=" * 70)
print("Critical slowing: n(eps) and K = n*sqrt(eps), eps = 1/4 - c")
print("=" * 70)
print(f"{'k':>3} {'eps':>12} {'n':>8} {'sqrt(eps)':>12} {'K':>14}")
print("-" * 60)

results = []
for k in range(1, 8):
    eps = 10.0**(-k)
    c = 0.25 - eps
    conv, u_inf, n = iterate_mandelbrot(c, max_iter=100000, tol=1e-12)
    if conv:
        K = n * np.sqrt(eps)
        results.append((k, eps, n, K))
        print(f"{k:3d} {eps:12.1e} {n:8d} {np.sqrt(eps):12.6e} {K:14.6f}")
    else:
        print(f"{k:3d} {eps:12.1e}  NOT CONVERGED (n={n})")

# Log-log fit: log(n) vs log(eps)
if len(results) >= 3:
    eps_arr = np.array([r[1] for r in results])
    n_arr = np.array([r[2] for r in results])
    slope, intercept = np.polyfit(np.log10(eps_arr), np.log10(n_arr), 1)
    print()
    print(f"Log-log fit: log(n) = {slope:.6f} * log(eps) + {intercept:.6f}")
    print(f"  -> n ~ {10**intercept:.4f} * eps^{slope:.4f}")
    print(f"  Expected saddle-node: slope = -0.5")
    print(f"  Deviation from -1/2: {slope - (-0.5):.6f}")

# Theoretische Vorhersage fuer K
# delta_{n+1} = delta_n + delta_n^2 - eps (near saddle node)
# passage time n ~ (pi/(2*sqrt(eps))) for full traversal from -infty to +infty
# For our initial u_0 = c = 1/4 - eps, delta_0 = -1/4 - eps, stopping at tol 1e-12
# Standard result: n * sqrt(eps) -> pi/2 ~ 1.5708 for full traversal
# But our n counts only until |u_{n+1} - u_n| < tol, which is earlier (near stable FP)
# So K should be smaller than pi and depend on tol logarithmically
#
# More precisely: near stable FP delta* = -sqrt(eps), convergence is linear
# with rate (1 + 2*delta*) = (1 - 2*sqrt(eps)). For tol convergence,
# we need ~ -ln(tol) / (2*sqrt(eps)) extra iterations. That is the dominant term.
# Predicted K = n*sqrt(eps) ~ -ln(tol)/2 = -ln(1e-12)/2 = 27.63/2 = 13.82
print()
print(f"Theoretical K from linear-convergence tail with tol=1e-12:")
print(f"  K_theory = -ln(1e-12)/2 = {-np.log(1e-12)/2:.4f}")
print()
print("DONE critical slowing scan")
