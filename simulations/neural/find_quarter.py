#!/usr/bin/env python3
"""
Find 1/4 in Wilson-Cowan. Six ideas, all tested.

The V-Effect was "not found" until the exact condition was tested.
The same pattern may apply to 1/4.
"""
import numpy as np
from scipy.linalg import eigvals
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


# ================================================================
# IDEA 1: sigma(theta) * (1 - sigma(theta)) = 1/4
# ================================================================
print("=" * 65)
print("IDEA 1: sigma(x) * (1 - sigma(x)) at the inflection point")
print("=" * 65)

print(f"""
For ANY sigmoid sigma(x) = 1/(1+exp(-a(x-theta))):
  sigma(theta) = 1/2  (at inflection point, always)
  sigma(theta) * (1 - sigma(theta)) = 1/2 * 1/2 = 1/4  (always)

This is the neural variance: how uncertain the neuron is.
Maximum uncertainty = 1/4. At the inflection point. Always.

The slope is a * sigma * (1-sigma) = a * 1/4 = a/4.
The VARIANCE (sigma*(1-sigma)) is 1/4 regardless of a.

Connection to quantum: CΨ = purity * coherence.
Neural analog: sigma * (1-sigma) = decided * undecided.
At the fold (maximum sensitivity): sigma*(1-sigma) = 1/4.
""")

# Verify for different sigmoid parameters
print("  Verification:")
for a, theta in [(1.0, 0.0), (1.3, 4.0), (2.0, 3.7), (5.0, 2.0), (0.1, 10.0)]:
    s = sigmoid(theta, a, theta)
    product = s * (1 - s)
    slope = a * product
    print(f"    a={a:4.1f}, theta={theta:4.1f}: sigma(theta)={s:.6f}, "
          f"sigma*(1-sigma)={product:.6f}, slope={slope:.6f}")


# ================================================================
# IDEA 2: V-Effect threshold
# ================================================================
print("\n" + "=" * 65)
print("IDEA 2: What happens at the V-Effect onset?")
print("=" * 65)

from veffect_exact import (build_exact_palindromic_network, build_linear_jacobian,
                            count_freqs, count_corr_freqs, palindrome_residual)

tau_E, tau_I = 5.0, 10.0
N = 20
n_exc = N // 2
W_A, signs_A, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=42)
W_B, signs_B, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=99)

N_c = 2 * N + 1
W_c = np.zeros((N_c, N_c))
signs_c = np.zeros(N_c)
W_c[:N, :N] = W_A; signs_c[:N] = signs_A
W_c[N:2*N, N:2*N] = W_B; signs_c[N:2*N] = signs_B
signs_c[2*N] = 1.0

print(f"\n  Coupling sweep (N={N} per network, exact palindromic):")
print(f"  {'coupling':>10s}  {'K_corr':>7s}  {'residual':>10s}  {'res/0.25':>8s}")
print(f"  {'-'*40}")

for coupling in [0.001, 0.003, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.1]:
    W_test = W_c.copy()
    for offset in [0, N]:
        W_test[2*N, offset] = coupling
        W_test[offset, 2*N] = coupling
        W_test[2*N, offset + N - 1] = coupling
        W_test[offset + N - 1, 2*N] = coupling

    J = build_linear_jacobian(W_test, signs_c, tau_E, tau_I, 0.5)
    ev = np.linalg.eigvals(J)
    K = count_corr_freqs(ev)
    res = palindrome_residual(J, signs_c)
    print(f"  {coupling:10.4f}  {K:7d}  {res:10.6f}  {res/0.25:8.4f}")


# ================================================================
# IDEA 3: CΨ_E * CΨ_I away from balance
# ================================================================
print("\n" + "=" * 65)
print("IDEA 3: Does CΨ_E * CΨ_I stay at 1/4 away from balance?")
print("=" * 65)

from scipy.linalg import expm

N = 20
tau_E, tau_I = 5.0, 10.0

print(f"\n  {'E:I':>6s}  {'CΨ_E(t=5)':>10s}  {'CΨ_I(t=5)':>10s}  {'product':>8s}  {'=1/4?':>6s}")
print(f"  {'-'*50}")

for n_e in [18, 15, 12, 10, 8, 5, 2]:
    n_i = N - n_e
    W, signs = make_balanced_network(N, n_e, density=0.3, seed=42)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, 0.3)

    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]

    rng = np.random.RandomState(1)
    delta = rng.randn(N)
    delta /= np.linalg.norm(delta)

    x_t = expm(J * 5.0) @ delta
    norm_E = np.sum(x_t[e_idx]**2)
    norm_I = np.sum(x_t[i_idx]**2)
    total = norm_E + norm_I
    if total > 1e-30:
        cE = norm_E / total
        cI = norm_I / total
    else:
        cE, cI = 0.5, 0.5

    product = cE * cI
    is_quarter = "YES" if abs(product - 0.25) < 0.01 else "no"
    print(f"  {n_e:2d}:{n_i:<2d}  {cE:10.4f}  {cI:10.4f}  {product:8.4f}  {is_quarter:>6s}")


# ================================================================
# IDEA 4: Normalized midpoint squared
# ================================================================
print("\n" + "=" * 65)
print("IDEA 4: Normalized midpoint m² = 1/4")
print("=" * 65)

print(f"""
  Eigenvalue pairing: mu + mu' = -(1/tau_E + 1/tau_I)
  Normalized: m = mu / (1/tau_E + 1/tau_I), so m + m' = -1
  Midpoint: m = -1/2
  Midpoint squared: m² = 1/4

  This is (0.5)² = 1/4 in disguise. Same as the axiom.
  Not new physics, just the axiom expressed in eigenvalue space.
""")


# ================================================================
# IDEA 5: Exact networks + V-Effect bifurcation
# ================================================================
print("=" * 65)
print("IDEA 5: V-Effect bifurcation in exact networks")
print("=" * 65)

# At what coupling does the FIRST imaginary eigenvalue appear?
N = 20
W_A, signs_A, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=42)
W_B, signs_B, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=99)

N_c = 2 * N + 1
W_c = np.zeros((N_c, N_c))
signs_c = np.zeros(N_c)
W_c[:N, :N] = W_A; signs_c[:N] = signs_A
W_c[N:2*N, N:2*N] = W_B; signs_c[N:2*N] = signs_B
signs_c[2*N] = 1.0

print(f"\n  Fine coupling sweep near V-Effect onset:")
print(f"  {'coupling':>10s}  {'n_complex':>9s}  {'max_Im':>8s}  {'onset?':>6s}")
print(f"  {'-'*40}")

prev_n = 0
onset_coupling = None
for coupling in np.arange(0.0001, 0.02, 0.0005):
    W_test = W_c.copy()
    for offset in [0, N]:
        W_test[2*N, offset] = coupling
        W_test[offset, 2*N] = coupling
        W_test[2*N, offset + N - 1] = coupling
        W_test[offset + N - 1, 2*N] = coupling

    J = build_linear_jacobian(W_test, signs_c, tau_E, tau_I, 0.5)
    ev = np.linalg.eigvals(J)
    n_complex = np.sum(np.abs(ev.imag) > 1e-10)

    if n_complex > 0 and prev_n == 0:
        onset_coupling = coupling
        print(f"  {coupling:10.5f}  {n_complex:9d}  {np.max(np.abs(ev.imag)):8.5f}  ONSET")
    elif coupling < 0.003 or (onset_coupling and coupling < onset_coupling + 0.005):
        print(f"  {coupling:10.5f}  {n_complex:9d}  {np.max(np.abs(ev.imag)):8.5f}")
    prev_n = n_complex

if onset_coupling:
    # Check if onset_coupling is 1/4 of something
    print(f"\n  Onset coupling: {onset_coupling:.5f}")
    print(f"  onset * 4 = {onset_coupling * 4:.5f}")
    print(f"  (Not obviously 1/4 of any known quantity)")


# ================================================================
# IDEA 6: Sigmoid gain at Hopf
# ================================================================
print("\n" + "=" * 65)
print("IDEA 6: What is the sigmoid gain at the Hopf bifurcation?")
print("=" * 65)

print(f"""
  The Hopf occurs when max Re(eigenvalue) crosses zero.
  The Jacobian depends on sigmoid gain g = f'(input) at the
  operating point. The gain is g = a * sigma * (1-sigma).

  At the inflection: g_max = a * 1/4 = a/4.

  Question: at the Hopf bifurcation, is the MEAN gain across
  all neurons equal to 1/4 of g_max? Or some other fraction
  involving 1/4?
""")

# Use approximate networks where Hopf exists
for N in [100, 200]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

    # Find alpha near Hopf (from previous computation)
    # For N=200: Hopf near alpha=5.87
    for alpha in [3.0, 5.0, 5.5, 5.87, 6.0, 7.0, 10.0]:
        x = np.ones(N) * 0.3
        for _ in range(3000):
            inputs = alpha * W @ x + 3.0
            x_new = np.zeros(N)
            for i in range(N):
                a_i = 1.3 if signs[i] > 0 else 2.0
                th_i = 4.0 if signs[i] > 0 else 3.7
                x_new[i] = sigmoid(inputs[i], a_i, th_i)
            if np.max(np.abs(x_new - x)) < 1e-13:
                break
            x = x_new

        # Compute gain at each neuron
        inputs = alpha * W @ x + 3.0
        gains = np.zeros(N)
        variances = np.zeros(N)
        for i in range(N):
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            s = sigmoid(inputs[i], a_i, th_i)
            gains[i] = a_i * s * (1 - s)
            variances[i] = s * (1 - s)  # sigma*(1-sigma), max = 1/4

        J = np.zeros((N, N))
        for i in range(N):
            tau_i = tau_E if signs[i] > 0 else tau_I
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            dS = dsigmoid(inputs[i], a_i, th_i)
            J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
            for j in range(N):
                if i != j:
                    J[i, j] = alpha * W[i, j] * dS / tau_i

        ev = eigvals(J)
        max_re = np.max(ev.real)
        mean_var = np.mean(variances)

        if N == 200:
            status = "HOPF" if max_re > 0 else "stable"
            print(f"  N={N} alpha={alpha:5.2f}: max_Re={max_re:+.4f} "
                  f"mean sigma*(1-sigma)={mean_var:.4f} "
                  f"(1/4={0.25:.4f}) {status}")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "=" * 65)
print("SUMMARY: Where is 1/4?")
print("=" * 65)
