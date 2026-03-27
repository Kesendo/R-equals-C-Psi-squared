#!/usr/bin/env python3
"""
V-Effect with EXACT palindromic symmetry.

The quantum V-Effect requires exact symmetry that then BREAKS:
  N=2: exact palindrome (all 22 pairs)
  N=3: 14/36 break, releasing 11 frequencies from 4

My earlier test found no neural V-Effect because the symmetry was
never exact. Nothing to break = no V-Effect.

Here: build PERFECTLY palindromic networks (Dale + exact magnitude
condition, residual = 0), then couple them and see what breaks.
"""
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def sigmoid(x, a=1.3, theta=4.0):
    arg = -a * (x - theta)
    arg = np.clip(arg, -500, 500)
    return 1.0 / (1.0 + np.exp(arg))


def dsigmoid(x, a=1.3, theta=4.0):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def build_exact_palindromic_network(N, n_exc, tau_E, tau_I, density=0.3, seed=42):
    """Build a network with EXACT palindromic symmetry.

    Uses Dale's Law signs AND the magnitude condition:
    W[Q(i),Q(j)] = -(tau_{Q(i)}/tau_i) * W[i,j]
    """
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

    # Build W satisfying the exact condition
    W = np.zeros((N, N))
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            qi, qj = perm[i], perm[j]
            # Only set one direction, derive the other
            if i < j:
                if mask[i, j]:
                    base = rng.exponential(0.3)
                    W[i, j] = signs[j] * base
                    tau_i = tau_E if signs[i] > 0 else tau_I
                    tau_qi = tau_E if signs[qi] > 0 else tau_I
                    W[qi, qj] = -(tau_qi / tau_i) * W[i, j]

    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx

    return W, signs, perm


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def palindrome_residual(J, signs):
    n = len(signs)
    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(n)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]
        perm[i_idx[k]] = e_idx[k]
    Q = np.zeros((n, n))
    for i in range(n):
        Q[i, perm[i]] = 1.0
    QJQ = Q @ J @ Q.T
    S_diag = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2 * np.diag(S_diag)
    R_off = R - np.diag(np.diag(R))
    norm_J = np.linalg.norm(J)
    return np.linalg.norm(R_off) / norm_J if norm_J > 0 else 0


def count_freqs(eigenvalues, tol=1e-6):
    freqs = np.abs(np.imag(eigenvalues))
    freqs = freqs[freqs > tol]
    return len(set(np.round(freqs, 6)))


def count_corr_freqs(eigenvalues, tol=1e-6):
    n = len(eigenvalues)
    corr_freqs = set()
    for i in range(n):
        for j in range(i, n):
            f = abs((eigenvalues[i] + eigenvalues[j]).imag)
            if f > tol:
                corr_freqs.add(round(f, 6))
    return len(corr_freqs)


# ================================================================
# Step 1: Verify exact palindrome for single network
# ================================================================
print("=" * 70)
print("STEP 1: Build perfectly palindromic neural networks")
print("=" * 70)

tau_E, tau_I = 5.0, 10.0
alpha = 0.5

for N in [10, 20, 30]:
    n_exc = N // 2
    W, signs, perm = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=42)
    J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
    res = palindrome_residual(J, signs)
    ev = np.linalg.eigvals(J)
    K = count_freqs(ev)
    K_corr = count_corr_freqs(ev)

    print(f"\n  N={N}: residual = {res:.2e}  "
          f"{'EXACT' if res < 1e-10 else 'NOT exact'}  "
          f"K_act={K}  K_corr={K_corr}")


# ================================================================
# Step 2: Couple two exact networks - does the palindrome break?
# ================================================================
print("\n" + "=" * 70)
print("STEP 2: Couple two EXACT palindromic networks")
print("Does the symmetry break at the coupling point?")
print("=" * 70)

for N in [10, 20]:
    n_exc = N // 2

    # Network A
    W_A, signs_A, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=42)
    J_A = build_linear_jacobian(W_A, signs_A, tau_E, tau_I, alpha)
    ev_A = np.linalg.eigvals(J_A)
    res_A = palindrome_residual(J_A, signs_A)

    # Network B
    W_B, signs_B, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=99)
    J_B = build_linear_jacobian(W_B, signs_B, tau_E, tau_I, alpha)
    ev_B = np.linalg.eigvals(J_B)
    res_B = palindrome_residual(J_B, signs_B)

    K_A = count_freqs(ev_A)
    K_B = count_freqs(ev_B)
    K_A_corr = count_corr_freqs(ev_A)

    # Coupled: A + B + mediator (excitatory)
    N_c = 2 * N + 1
    W_c = np.zeros((N_c, N_c))
    signs_c = np.zeros(N_c)

    W_c[:N, :N] = W_A
    signs_c[:N] = signs_A
    W_c[N:2*N, N:2*N] = W_B
    signs_c[N:2*N] = signs_B
    signs_c[2*N] = 1.0  # mediator

    # Coupling: mediator connects to edge E-neurons of each network
    for coupling_strength in [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0]:
        W_test = W_c.copy()
        for offset in [0, N]:
            W_test[2*N, offset] = coupling_strength
            W_test[offset, 2*N] = coupling_strength
            W_test[2*N, offset + N - 1] = coupling_strength
            W_test[offset + N - 1, 2*N] = coupling_strength

        J_c = build_linear_jacobian(W_test, signs_c, tau_E, tau_I, alpha)
        ev_c = np.linalg.eigvals(J_c)
        res_c = palindrome_residual(J_c, signs_c)
        K_c = count_freqs(ev_c)
        K_c_corr = count_corr_freqs(ev_c)

        # How many NEW frequencies vs 2x single?
        new_act = K_c - (K_A + K_B)
        new_corr = K_c_corr - 2 * K_A_corr

        if coupling_strength == 0.0:
            print(f"\n  N={N}: single A: K={K_A}, K_corr={K_A_corr}, res={res_A:.2e}")
            print(f"  Coupling  Residual    K_act  K_corr  New_act  New_corr  V-Effect?")
            print(f"  {'-'*65}")

        v_marker = "V-EFFECT!" if new_corr > 0 else ""
        print(f"    {coupling_strength:5.2f}    {res_c:.2e}   {K_c:5d}  {K_c_corr:7d}  "
              f"{new_act:+7d}  {new_corr:+8d}  {v_marker}")


# ================================================================
# Step 3: Heat (drive P) on exact palindromic network
# ================================================================
print("\n" + "=" * 70)
print("STEP 3: Heat on EXACT palindromic network")
print("Do dying modes create exactly 2 new modes per step?")
print("=" * 70)

N = 20
n_exc = N // 2
W, signs, _ = build_exact_palindromic_network(N, n_exc, tau_E, tau_I, seed=42)

# Build Jacobian WITH sigmoid (nonlinear, P-dependent)
print(f"\n  N={N}, exact palindromic at linear level")
print(f"\n  {'P':>6s}  {'n_osc':>6s}  {'K_freq':>6s}  {'K_corr':>7s}  "
      f"{'delta_K':>7s}  {'residual':>10s}")
print(f"  {'-'*55}")

prev_K = None
for P in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]:
    # With sigmoid: operating point shifts, effective coupling changes
    n = len(signs)
    x = np.ones(n) * 0.3
    for _ in range(500):
        inputs = alpha * W @ x + P
        for i in range(n):
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            x[i] = sigmoid(inputs[i], a_i, th_i)

    inputs = alpha * W @ x + P
    J_nl = np.zeros((n, n))
    for i in range(n):
        a_i = 1.3 if signs[i] > 0 else 2.0
        th_i = 4.0 if signs[i] > 0 else 3.7
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J_nl[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(n):
            if i != j:
                J_nl[i, j] = alpha * W[i, j] * dS / tau_i

    ev_nl = np.linalg.eigvals(J_nl)
    n_osc = np.sum(np.abs(np.imag(ev_nl)) > 1e-5)
    K = count_freqs(ev_nl, tol=1e-5)
    K_corr = count_corr_freqs(ev_nl, tol=1e-5)
    res = palindrome_residual(J_nl, signs)

    delta = K_corr - prev_K if prev_K is not None else 0
    prev_K = K_corr

    marker = ""
    if delta == 2:
        marker = " <-- +2 NEW"
    elif delta > 0:
        marker = f" <-- +{delta} new"

    print(f"  {P:4.1f}  {n_osc:6d}  {K:6d}  {K_corr:7d}  {delta:+7d}  "
          f"{res:10.2e}{marker}")
