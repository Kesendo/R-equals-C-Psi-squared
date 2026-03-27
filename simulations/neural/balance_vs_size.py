#!/usr/bin/env python3
"""
The real question: is the threshold BALANCE or SIZE?

At perfect C=0.5 (exact palindromic condition), does the Hopf
threshold alpha_c go to zero as N grows? If yes: balance alone
suffices. If no: you need both balance AND size.

Test: build EXACT palindromic networks at N=4,10,20,50,100,200.
Couple two of them. Find the Hopf threshold of the COUPLED system.
Track alpha_c vs N.
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


def build_exact_palindromic(N, n_exc, tau_E, tau_I, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    inh_idx = rng.choice(N, N - n_exc, replace=False)
    signs[inh_idx] = -1

    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    n_pairs = min(len(e_idx), len(i_idx))
    perm = np.arange(N)
    for k in range(n_pairs):
        perm[e_idx[k]] = i_idx[k]
        perm[i_idx[k]] = e_idx[k]

    W = np.zeros((N, N))
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            qi, qj = perm[i], perm[j]
            if i < j:
                if mask[i, j]:
                    base = rng.exponential(0.3)
                    W[i, j] = signs[j] * base
                    tau_i_val = tau_E if signs[i] > 0 else tau_I
                    tau_qi_val = tau_E if signs[qi] > 0 else tau_I
                    W[qi, qj] = -(tau_qi_val / tau_i_val) * W[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def build_coupled_system(W_A, signs_A, W_B, signs_B, coupling, tau_E, tau_I):
    """Couple two networks through a mediator neuron."""
    N = len(signs_A)
    N_c = 2 * N + 1
    W_c = np.zeros((N_c, N_c))
    signs_c = np.zeros(N_c)

    W_c[:N, :N] = W_A
    signs_c[:N] = signs_A
    W_c[N:2*N, N:2*N] = W_B
    signs_c[N:2*N] = signs_B
    signs_c[2*N] = 1.0  # excitatory mediator

    for offset in [0, N]:
        W_c[2*N, offset] = coupling
        W_c[offset, 2*N] = coupling
        W_c[2*N, offset + N - 1] = coupling
        W_c[offset + N - 1, 2*N] = coupling

    return W_c, signs_c


def max_real_ev(W, signs, tau_E, tau_I, alpha, P=3.0):
    """Max Re(eigenvalue) of Jacobian at fixed point (with sigmoid)."""
    N = len(signs)
    x = np.ones(N) * 0.3
    for _ in range(3000):
        inputs = alpha * W @ x + P
        x_new = np.zeros(N)
        for i in range(N):
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            x_new[i] = sigmoid(inputs[i], a_i, th_i)
        if np.max(np.abs(x_new - x)) < 1e-13:
            break
        x = x_new

    inputs = alpha * W @ x + P
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
    return np.max(eigvals(J).real)


tau_E, tau_I = 5.0, 10.0
P = 3.0

# ================================================================
# Test 1: Single exact palindromic networks - always stable?
# ================================================================
print("=" * 70)
print("TEST 1: Single exact palindromic network stability")
print("=" * 70)

print(f"\n{'N':>5s}  {'alpha':>6s}  {'max_Re':>10s}  {'stable':>7s}")
print("-" * 35)

for N in [4, 10, 20, 50, 100]:
    n_exc = N // 2
    W, signs = build_exact_palindromic(N, n_exc, tau_E, tau_I, seed=42)
    for alpha in [1.0, 5.0, 20.0]:
        mr = max_real_ev(W, signs, tau_E, tau_I, alpha, P)
        print(f"  {N:3d}  {alpha:6.1f}  {mr:10.6f}  {'YES' if mr < 0 else 'NO!'}")


# ================================================================
# Test 2: Coupled exact palindromic - find Hopf threshold
# ================================================================
print("\n" + "=" * 70)
print("TEST 2: Coupled exact palindromic networks - Hopf threshold")
print("C = 0.5 exact. Does alpha_c go to zero with N?")
print("=" * 70)

results = []

for N in [4, 10, 20, 50, 100]:
    n_exc = N // 2
    print(f"\n  N={N} (coupled system: {2*N+1} neurons)...", end=" ")
    sys.stdout.flush()

    W_A, signs_A = build_exact_palindromic(N, n_exc, tau_E, tau_I, seed=42)
    W_B, signs_B = build_exact_palindromic(N, n_exc, tau_E, tau_I, seed=99)

    # Bracket the Hopf
    alpha_lo, alpha_hi = None, None
    for alpha_test in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        W_c, signs_c = build_coupled_system(W_A, signs_A, W_B, signs_B,
                                             coupling=0.3, tau_E=tau_E, tau_I=tau_I)
        mr = max_real_ev(W_c, signs_c, tau_E, tau_I, alpha_test, P)
        if mr < 0:
            alpha_lo = alpha_test
        else:
            alpha_hi = alpha_test
            break

    if alpha_lo is None or alpha_hi is None:
        if alpha_hi is None:
            print(f"always stable up to alpha={alpha_test} (max_Re={mr:.4f})")
        else:
            print(f"already unstable at alpha={alpha_test}")
            # Try lower
            for alpha_test in [0.001, 0.005, 0.01, 0.03, 0.05]:
                mr = max_real_ev(W_c, signs_c, tau_E, tau_I, alpha_test, P)
                if mr < 0:
                    alpha_lo = alpha_test
                    break
            if alpha_lo is None:
                print(f"  unstable even at alpha=0.001!")
                results.append((N, 0.001))
                continue

    if alpha_lo is not None and alpha_hi is not None:
        # Bisection
        for _ in range(30):
            alpha_mid = (alpha_lo + alpha_hi) / 2
            mr = max_real_ev(W_c, signs_c, tau_E, tau_I, alpha_mid, P)
            if mr < 0:
                alpha_lo = alpha_mid
            else:
                alpha_hi = alpha_mid
            if alpha_hi - alpha_lo < 0.001:
                break
        alpha_c = (alpha_lo + alpha_hi) / 2
        results.append((N, alpha_c))
        print(f"alpha_c = {alpha_c:.4f}")


# ================================================================
# Test 3: Also test with different coupling strengths
# ================================================================
print("\n" + "=" * 70)
print("TEST 3: Does the inter-network coupling strength matter?")
print("=" * 70)

N = 50
n_exc = N // 2
W_A, signs_A = build_exact_palindromic(N, n_exc, tau_E, tau_I, seed=42)
W_B, signs_B = build_exact_palindromic(N, n_exc, tau_E, tau_I, seed=99)

print(f"\n  N={N}, varying inter-network coupling:")
print(f"  {'coupling':>8s}  {'alpha_c':>8s}")
print(f"  {'-'*20}")

for coupling in [0.05, 0.1, 0.3, 0.5, 1.0]:
    W_c, signs_c = build_coupled_system(W_A, signs_A, W_B, signs_B,
                                         coupling=coupling, tau_E=tau_E, tau_I=tau_I)
    alpha_lo, alpha_hi = None, None
    for alpha_test in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        mr = max_real_ev(W_c, signs_c, tau_E, tau_I, alpha_test, P)
        if mr < 0:
            alpha_lo = alpha_test
        else:
            alpha_hi = alpha_test
            break
    if alpha_lo and alpha_hi:
        for _ in range(20):
            alpha_mid = (alpha_lo + alpha_hi) / 2
            mr = max_real_ev(W_c, signs_c, tau_E, tau_I, alpha_mid, P)
            if mr < 0:
                alpha_lo = alpha_mid
            else:
                alpha_hi = alpha_mid
            if alpha_hi - alpha_lo < 0.01:
                break
        print(f"  {coupling:8.2f}  {(alpha_lo+alpha_hi)/2:8.4f}")
    else:
        print(f"  {coupling:8.2f}  not found")


# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 70)
print("SUMMARY: alpha_c(N) for coupled exact palindromic networks")
print("=" * 70)

if results:
    print(f"\n  {'N':>5s}  {'alpha_c':>8s}")
    print(f"  {'-'*15}")
    for N_val, ac in results:
        print(f"  {N_val:5d}  {ac:8.4f}")

    N_arr = np.array([r[0] for r in results], dtype=float)
    ac_arr = np.array([r[1] for r in results], dtype=float)

    if len(results) >= 3 and all(ac > 0 for _, ac in results):
        # Log-log fit
        log_N = np.log(N_arr)
        log_ac = np.log(ac_arr)
        coeffs = np.polyfit(log_N, log_ac, 1)
        beta = -coeffs[0]
        A = np.exp(coeffs[1])
        print(f"\n  Power law: alpha_c = {A:.2f} * N^(-{beta:.3f})")

        if beta > 0:
            N_c_03 = (A / 0.3) ** (1.0 / beta)
            print(f"  At alpha = 0.3: N_c = {N_c_03:.0f}")
            print(f"\n  Does alpha_c go to zero? {'YES (beta > 0)' if beta > 0 else 'NO'}")
            print(f"  Balance alone suffices at large enough N.")
        else:
            print(f"\n  alpha_c does NOT decrease with N (beta = {beta:.3f})")
            print(f"  Both balance AND coupling needed.")
