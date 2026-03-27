#!/usr/bin/env python3
"""
The Complexity Threshold: does a critical N_c exist?

At what N does a Wilson-Cowan network sustain oscillation without
external driving? This requires FULL NONLINEAR dynamics (sigmoid),
not the linearized Jacobian.

Test:
  1. Build balanced E/I networks at increasing N
  2. Find the steady state (fixed point of nonlinear dynamics)
  3. Perturb it
  4. Simulate FULL nonlinear Wilson-Cowan (with sigmoid)
  5. Measure: does oscillation amplitude decay or sustain?
  6. Track the Jacobian's max Re(eigenvalue) - does it approach zero?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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
    """Find steady state by iteration."""
    N = len(signs)
    x = np.ones(N) * 0.3
    for _ in range(2000):
        inputs = alpha * W @ x + P
        x_new = np.zeros(N)
        for i in range(N):
            a_i = a_E if signs[i] > 0 else a_I
            th_i = theta_E if signs[i] > 0 else theta_I
            x_new[i] = sigmoid(inputs[i], a_i, th_i)
        if np.max(np.abs(x_new - x)) < 1e-12:
            break
        x = x_new
    return x


def build_jacobian_at(x_star, W, signs, tau_E, tau_I, alpha,
                      a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7):
    """Jacobian of full nonlinear Wilson-Cowan at fixed point."""
    N = len(signs)
    inputs = alpha * W @ x_star + 1.5  # P included
    J = np.zeros((N, N))
    for i in range(N):
        a_i = a_E if signs[i] > 0 else a_I
        th_i = theta_E if signs[i] > 0 else theta_I
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(N):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS / tau_i
    return J


def wilson_cowan_rhs(t, x, W, signs, tau_E, tau_I, alpha, P,
                     a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7):
    """Full nonlinear Wilson-Cowan right-hand side."""
    N = len(signs)
    inputs = alpha * W @ x + P
    dxdt = np.zeros(N)
    for i in range(N):
        a_i = a_E if signs[i] > 0 else a_I
        th_i = theta_E if signs[i] > 0 else theta_I
        tau_i = tau_E if signs[i] > 0 else tau_I
        dxdt[i] = (-x[i] + sigmoid(inputs[i], a_i, th_i)) / tau_i
    return dxdt


def measure_oscillation(sol_y, t_start_frac=0.5):
    """Measure oscillation amplitude in the second half of the trajectory."""
    N, T = sol_y.shape
    t_start = int(T * t_start_frac)
    if t_start >= T - 10:
        return 0

    # Use variance of activity as oscillation measure
    late = sol_y[:, t_start:]
    # Mean variance across neurons
    var_per_neuron = np.var(late, axis=1)
    return np.mean(var_per_neuron)


# ================================================================
# Test 1: Jacobian stability vs N
# ================================================================
print("=" * 70)
print("TEST 1: Does the largest eigenvalue approach zero with N?")
print("(If it crosses zero: Hopf bifurcation, self-sustaining oscillation)")
print("=" * 70)

tau_E, tau_I = 5.0, 10.0

print(f"\n{'N':>5s}  {'alpha':>6s}  {'P':>4s}  {'max_Re':>10s}  {'stable?':>8s}  "
      f"{'n_osc':>6s}  {'max_Im':>10s}")
print("-" * 60)

for N in [10, 20, 50, 100, 200, 500]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

    for alpha, P in [(0.3, 1.5), (1.0, 1.5), (2.0, 1.5), (3.0, 1.5),
                     (1.0, 3.0), (2.0, 3.0), (3.0, 3.0)]:
        x_star = find_fixed_point(W, signs, tau_E, tau_I, alpha, P)
        J = build_jacobian_at(x_star, W, signs, tau_E, tau_I, alpha)
        ev = eigvals(J)

        max_re = np.max(ev.real)
        stable = "YES" if max_re < 0 else "NO (HOPF!)"
        n_osc = np.sum(np.abs(ev.imag) > 1e-6)
        max_im = np.max(np.abs(ev.imag))

        if max_re > -0.01 or N in [10, 100, 500]:
            print(f"  {N:3d}  {alpha:6.1f}  {P:4.1f}  {max_re:10.6f}  "
                  f"{stable:>8s}  {n_osc:6d}  {max_im:10.4f}")


# ================================================================
# Test 2: Full nonlinear simulation at different N
# ================================================================
print("\n" + "=" * 70)
print("TEST 2: Full nonlinear Wilson-Cowan dynamics")
print("Does oscillation sustain or decay?")
print("=" * 70)

for N in [10, 20, 50, 100, 200]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

    for alpha, P in [(1.0, 1.5), (2.0, 3.0), (3.0, 3.0)]:
        x_star = find_fixed_point(W, signs, tau_E, tau_I, alpha, P)

        # Perturbation
        rng = np.random.RandomState(1)
        delta = rng.randn(N) * 0.05
        x0 = np.clip(x_star + delta, 0, 1)

        # Simulate
        T_max = 500.0
        sol = solve_ivp(wilson_cowan_rhs, [0, T_max], x0,
                        args=(W, signs, tau_E, tau_I, alpha, P),
                        method='RK45', max_step=1.0,
                        t_eval=np.linspace(0, T_max, 2000))

        if sol.success:
            # Measure oscillation in early and late parts
            osc_early = measure_oscillation(sol.y, t_start_frac=0.0)
            osc_late = measure_oscillation(sol.y, t_start_frac=0.7)
            ratio = osc_late / osc_early if osc_early > 1e-15 else 0

            status = "SUSTAINED" if ratio > 0.5 else \
                     "DECAYING" if ratio > 0.01 else "DEAD"

            print(f"  N={N:3d}  alpha={alpha:.1f}  P={P:.1f}  "
                  f"early={osc_early:.2e}  late={osc_late:.2e}  "
                  f"ratio={ratio:.4f}  {status}")
        else:
            print(f"  N={N:3d}  alpha={alpha:.1f}  P={P:.1f}  FAILED")


# ================================================================
# Test 3: Sweep alpha at fixed N to find Hopf bifurcation
# ================================================================
print("\n" + "=" * 70)
print("TEST 3: Find the Hopf bifurcation (alpha sweep)")
print("=" * 70)

for N in [50, 100, 200]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)
    P = 3.0

    print(f"\n  N={N}, P={P}:")
    print(f"  {'alpha':>6s}  {'max_Re':>10s}  {'n_osc':>6s}  {'status':>10s}")
    print(f"  {'-'*40}")

    prev_re = None
    for alpha in np.arange(0.5, 10.1, 0.5):
        x_star = find_fixed_point(W, signs, tau_E, tau_I, alpha, P)
        J = build_jacobian_at(x_star, W, signs, tau_E, tau_I, alpha)
        ev = eigvals(J)
        max_re = np.max(ev.real)
        n_osc = np.sum(np.abs(ev.imag) > 1e-6)

        status = "stable" if max_re < 0 else "HOPF!"
        crossed = ""
        if prev_re is not None and prev_re < 0 and max_re >= 0:
            crossed = " <-- CROSSING"
        prev_re = max_re

        print(f"  {alpha:6.1f}  {max_re:10.6f}  {n_osc:6d}  {status:>10s}{crossed}")
