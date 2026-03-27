#!/usr/bin/env python3
"""
The Hopf threshold: where does life begin?

1. Confirm the limit cycle at N=200, alpha=8.5 (nonlinear simulation)
2. Find alpha_c(N) for N=200, 500, 1000, 2000, 5000
3. Fit the curve alpha_c(N)
4. Extrapolate: at what N does alpha_c cross 0.3 (biological)?

That N is N_c. The complexity threshold.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
from scipy.optimize import brentq, curve_fit
import sys

def sigmoid(x, a=1.3, theta=4.0):
    arg = -a * (x - theta)
    arg = np.clip(arg, -500, 500)
    return 1.0 / (1.0 + np.exp(arg))

def dsigmoid(x, a=1.3, theta=4.0):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)

def make_balanced_network(N, n_exc, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N - n_exc, replace=False)
    signs[inh_idx] = -1
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs

def find_fixed_point(W, signs, tau_E, tau_I, alpha, P,
                     a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7):
    N = len(signs)
    x = np.ones(N) * 0.3
    for _ in range(3000):
        inputs = alpha * W @ x + P
        x_new = np.zeros(N)
        for i in range(N):
            a_i = a_E if signs[i] > 0 else a_I
            th_i = theta_E if signs[i] > 0 else theta_I
            x_new[i] = sigmoid(inputs[i], a_i, th_i)
        if np.max(np.abs(x_new - x)) < 1e-13:
            break
        x = x_new
    return x

def max_real_eigenvalue(W, signs, tau_E, tau_I, alpha, P):
    """Compute max Re(eigenvalue) of the Jacobian at the fixed point."""
    N = len(signs)
    x_star = find_fixed_point(W, signs, tau_E, tau_I, alpha, P)
    inputs = alpha * W @ x_star + P
    J = np.zeros((N, N))
    for i in range(N):
        a_i = 1.3 if signs[i] > 0 else 2.0
        th_i = 4.0 if signs[i] > 0 else 3.7
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(N):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS / tau_i
    ev = eigvals(J)
    return np.max(ev.real)

def wilson_cowan_rhs(t, x, W, signs, tau_E, tau_I, alpha, P):
    N = len(signs)
    inputs = alpha * W @ x + P
    dxdt = np.zeros(N)
    for i in range(N):
        a_i = 1.3 if signs[i] > 0 else 2.0
        th_i = 4.0 if signs[i] > 0 else 3.7
        tau_i = tau_E if signs[i] > 0 else tau_I
        dxdt[i] = (-x[i] + sigmoid(inputs[i], a_i, th_i)) / tau_i
    return dxdt

tau_E, tau_I = 5.0, 10.0
P = 3.0

# ================================================================
# PART 1: Confirm limit cycle at N=200, alpha=8.5
# ================================================================
print("=" * 70)
print("PART 1: Limit cycle at N=200, alpha=8.5")
print("=" * 70)

N = 200
n_exc = N // 2
W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)
alpha_hopf = 8.5

x_star = find_fixed_point(W, signs, tau_E, tau_I, alpha_hopf, P)
rng = np.random.RandomState(1)
delta = rng.randn(N) * 0.01
x0 = np.clip(x_star + delta, 0, 1)

print(f"  Simulating N={N}, alpha={alpha_hopf}, P={P}...")
print(f"  Fixed point: mean={np.mean(x_star):.4f}, range=[{np.min(x_star):.4f}, {np.max(x_star):.4f}]")
sys.stdout.flush()

sol = solve_ivp(wilson_cowan_rhs, [0, 300], x0,
                args=(W, signs, tau_E, tau_I, alpha_hopf, P),
                method='RK45', max_step=0.5,
                t_eval=np.linspace(0, 300, 3000))

if sol.success:
    # Measure oscillation in three time windows
    for label, t_lo, t_hi in [("early (0-50)", 0, 50),
                               ("mid (100-200)", 100, 200),
                               ("late (200-300)", 200, 300)]:
        mask = (sol.t >= t_lo) & (sol.t <= t_hi)
        if np.sum(mask) > 10:
            segment = sol.y[:, mask]
            var = np.mean(np.var(segment, axis=1))
            amp = np.max(segment) - np.min(segment)
            print(f"  {label}: variance={var:.6e}, amplitude={amp:.4f}")

    # Check if late oscillation is sustained
    late_mask = sol.t > 200
    late_var = np.mean(np.var(sol.y[:, late_mask], axis=1))
    early_var = np.mean(np.var(sol.y[:, sol.t < 50], axis=1))
    if late_var > 0.1 * early_var:
        print(f"\n  >>> LIMIT CYCLE CONFIRMED: late variance {late_var:.2e} "
              f"is {late_var/early_var:.0%} of early")
    else:
        print(f"\n  >>> Oscillation decaying: late/early = {late_var/early_var:.2e}")
else:
    print(f"  Integration failed: {sol.message}")

# ================================================================
# PART 2: Find alpha_c(N) by bisection
# ================================================================
print("\n" + "=" * 70)
print("PART 2: Hopf threshold alpha_c(N)")
print("=" * 70)

alpha_c_data = []  # list of (N, alpha_c) tuples
N_values = [100, 200, 500, 1000, 2000, 5000]

for N in N_values:
    n_exc = N // 2
    print(f"\n  N={N}...", end=" ")
    sys.stdout.flush()

    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

    # First: bracket the Hopf bifurcation
    alpha_lo, alpha_hi = None, None

    for alpha_test in [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]:
        try:
            mr = max_real_eigenvalue(W, signs, tau_E, tau_I, alpha_test, P)
        except Exception:
            break
        if mr < 0:
            alpha_lo = alpha_test
        else:
            alpha_hi = alpha_test
            break

    if alpha_lo is None or alpha_hi is None:
        print(f"could not bracket (last tested: alpha={alpha_test}, max_Re={mr:.4f})")
        continue

    # Bisection to find alpha_c
    for _ in range(30):
        alpha_mid = (alpha_lo + alpha_hi) / 2
        mr = max_real_eigenvalue(W, signs, tau_E, tau_I, alpha_mid, P)
        if mr < 0:
            alpha_lo = alpha_mid
        else:
            alpha_hi = alpha_mid
        if alpha_hi - alpha_lo < 0.001:
            break

    alpha_c = (alpha_lo + alpha_hi) / 2
    alpha_c_data.append((N, alpha_c))
    print(f"alpha_c = {alpha_c:.4f}")

# ================================================================
# PART 3: Fit and extrapolate
# ================================================================
print("\n" + "=" * 70)
print("PART 3: Fit alpha_c(N) and extrapolate to alpha=0.3")
print("=" * 70)

if len(alpha_c_data) >= 3:
    N_arr = np.array([d[0] for d in alpha_c_data], dtype=float)
    ac_arr = np.array([d[1] for d in alpha_c_data], dtype=float)

    print(f"\n  Data points:")
    print(f"  {'N':>8s}  {'alpha_c':>10s}")
    for n, ac in zip(N_arr, ac_arr):
        print(f"  {n:8.0f}  {ac:10.4f}")

    # Try power law fit: alpha_c = A * N^(-beta)
    try:
        def power_law(N, A, beta):
            return A * N**(-beta)

        popt, pcov = curve_fit(power_law, N_arr, ac_arr, p0=[100, 0.5],
                               maxfev=10000)
        A, beta = popt

        print(f"\n  Power law fit: alpha_c(N) = {A:.2f} * N^(-{beta:.4f})")

        # Residuals
        ac_pred = power_law(N_arr, A, beta)
        for n, ac, acp in zip(N_arr, ac_arr, ac_pred):
            print(f"    N={n:6.0f}: measured={ac:.4f}, predicted={acp:.4f}, "
                  f"error={abs(ac-acp)/ac*100:.1f}%")

        # Extrapolate to alpha_c = 0.3
        # A * N_c^(-beta) = 0.3  =>  N_c = (A/0.3)^(1/beta)
        alpha_bio = 0.3
        N_c = (A / alpha_bio) ** (1.0 / beta)

        print(f"\n  >>> Extrapolation: alpha_c = {alpha_bio}")
        print(f"  >>> N_c = {N_c:.0f}")
        print(f"  >>> At N = {N_c:.0f} neurons with balanced E/I,")
        print(f"      biological coupling alpha = {alpha_bio} is sufficient")
        print(f"      for self-sustaining oscillation (Hopf bifurcation).")

        # Sanity checks
        print(f"\n  Sanity checks:")
        for N_check in [100, 500, 1000, 5000, 10000, 50000, 100000]:
            ac_check = power_law(N_check, A, beta)
            print(f"    N={N_check:>7d}: predicted alpha_c = {ac_check:.4f}")

    except Exception as e:
        print(f"\n  Fit failed: {e}")
        # Fallback: log-log linear fit
        log_N = np.log(N_arr)
        log_ac = np.log(ac_arr)
        coeffs = np.polyfit(log_N, log_ac, 1)
        beta = -coeffs[0]
        A = np.exp(coeffs[1])
        print(f"\n  Log-log linear fit: alpha_c = {A:.2f} * N^(-{beta:.4f})")
        N_c = (A / 0.3) ** (1.0 / beta)
        print(f"  N_c = {N_c:.0f}")

else:
    print(f"  Not enough data points for fit ({len(alpha_c_data)})")

# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
