#!/usr/bin/env python3
"""
V-Effect and Thermal Window in Neural Networks

V-Effect: coupling creates new frequencies. In activity space (N dim),
the effect is modest. In correlation space (N^2 dim), it's quadratic:
K frequencies in J become up to K(K+1)/2 in the correlation Liouvillian
L_C = J x I + I x J^T, because eigenvalues of L_C are all pairwise
sums lambda_i + lambda_j.

Thermal window: the external drive P shifts the sigmoid operating point.
Low P: sigmoid is flat, no oscillation. Moderate P: maximum slope,
oscillation. High P: saturation, oscillation dies. The palindromic
quality should peak in the oscillatory window.

No large matrices needed: correlation eigenvalues = pairwise sums of
J-eigenvalues. Works at N=1000 in seconds.
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


def build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, alpha, P,
                                 a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7):
    """Build Jacobian at steady state with sigmoid nonlinearity.

    The operating point depends on P (external drive = temperature).
    """
    n = len(signs)
    # Find steady state iteratively
    x = np.ones(n) * 0.3
    for _ in range(500):
        inputs = alpha * W @ x + P
        for i in range(n):
            a_i = a_E if signs[i] > 0 else a_I
            th_i = theta_E if signs[i] > 0 else theta_I
            x[i] = sigmoid(inputs[i], a_i, th_i)

    # Build Jacobian at steady state
    inputs = alpha * W @ x + P
    J = np.zeros((n, n))
    for i in range(n):
        a_i = a_E if signs[i] > 0 else a_I
        th_i = theta_E if signs[i] > 0 else theta_I
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS / tau_i
    return J, x


def count_distinct_freqs(eigenvalues, tol=1e-4):
    """Count distinct oscillation frequencies."""
    freqs = np.abs(np.imag(eigenvalues))
    freqs = freqs[freqs > tol]
    if len(freqs) == 0:
        return 0
    freqs_rounded = np.round(freqs / tol) * tol
    return len(set(freqs_rounded))


def correlation_freqs(eigenvalues, tol=1e-4):
    """Compute distinct frequencies in correlation space.

    L_C = J x I + I x J^T has eigenvalues lambda_i + lambda_j.
    Frequencies = |Im(lambda_i + lambda_j)| for all i, j.
    """
    n = len(eigenvalues)
    corr_freqs = set()
    for i in range(n):
        for j in range(i, n):
            f = abs((eigenvalues[i] + eigenvalues[j]).imag)
            if f > tol:
                corr_freqs.add(round(f / tol) * tol)
    return len(corr_freqs)


def palindrome_quality(J, signs):
    """Quick palindrome residual."""
    n = len(signs)
    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    n_pairs = min(len(e_idx), len(i_idx))
    perm = np.arange(n)
    for k in range(n_pairs):
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
    return 1.0 - np.linalg.norm(R_off) / norm_J if norm_J > 0 else 0


# ================================================================
# V-EFFECT: Frequency multiplication through coupling
# ================================================================
print("=" * 70)
print("V-EFFECT: Frequency Multiplication")
print("Coupling two networks in correlation space")
print("=" * 70)

tau_E, tau_I = 5.0, 10.0

print(f"\n{'N':>5s}  {'K_single':>8s}  {'K_corr':>8s}  {'K_coupled':>9s}  "
      f"{'K_c_corr':>9s}  {'V-act':>6s}  {'V-corr':>7s}")
print("-" * 65)

for N in [10, 20, 50, 100, 200]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)
    J, _ = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P=1.5)
    ev = np.linalg.eigvals(J)

    K_single = count_distinct_freqs(ev)
    K_corr = correlation_freqs(ev)

    # Coupled: two networks + mediator
    N_c = 2 * N + 1
    W2, signs2 = make_balanced_network(N, n_exc, density=0.3, seed=99)
    W_coupled = np.zeros((N_c, N_c))
    signs_coupled = np.zeros(N_c)
    W_coupled[:N, :N] = W
    signs_coupled[:N] = signs
    W_coupled[N:2*N, N:2*N] = W2
    signs_coupled[N:2*N] = signs2
    signs_coupled[2*N] = 1.0  # mediator is excitatory

    # Couple mediator to edges of both networks
    for offset in [0, N]:
        W_coupled[2*N, offset] = 0.3
        W_coupled[offset, 2*N] = 0.3
        W_coupled[2*N, offset + N - 1] = 0.3
        W_coupled[offset + N - 1, 2*N] = 0.3

    J_c, _ = build_jacobian_with_sigmoid(W_coupled, signs_coupled,
                                          tau_E, tau_I, 0.3, P=1.5)
    ev_c = np.linalg.eigvals(J_c)

    K_coupled = count_distinct_freqs(ev_c)
    K_c_corr = correlation_freqs(ev_c)

    v_act = K_coupled / (2 * K_single) if K_single > 0 else 0
    v_corr = K_c_corr / (2 * K_corr) if K_corr > 0 else 0

    print(f"  {N:3d}  {K_single:8d}  {K_corr:8d}  {K_coupled:9d}  "
          f"{K_c_corr:9d}  {v_act:6.2f}  {v_corr:7.2f}")

print("""
V-act  = K_coupled / (2 * K_single) in activity space
V-corr = K_coupled_corr / (2 * K_single_corr) in correlation space
V > 1 means coupling creates MORE than 2x the frequencies (V-Effect)
""")

# ================================================================
# THERMAL WINDOW: External drive P as temperature
# ================================================================
print("=" * 70)
print("THERMAL WINDOW: External Drive as Temperature")
print("=" * 70)

N = 50
n_exc = 25
W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)

print(f"\nN={N}, sweeping P (external drive)")
print(f"\n{'P':>6s}  {'n_osc':>6s}  {'E_freq':>8s}  {'E_decay':>8s}  "
      f"{'ratio':>6s}  {'K_freq':>6s}  {'K_corr':>7s}  {'pal_q':>6s}")
print("-" * 65)

for P in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0]:
    J, x_star = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P)
    ev = np.linalg.eigvals(J)

    E_freq = np.sum(np.abs(np.imag(ev)))
    E_decay = np.sum(np.abs(np.real(ev)))
    ratio = E_freq / E_decay if E_decay > 0 else 0
    n_osc = np.sum(np.abs(np.imag(ev)) > 1e-4)

    K_freq = count_distinct_freqs(ev)
    K_corr = correlation_freqs(ev)
    pq = palindrome_quality(J, signs)

    marker = ""
    if abs(ratio - 1.0) < 0.1:
        marker = " <-- crossover"

    print(f"  {P:4.1f}  {n_osc:6d}  {E_freq:8.3f}  {E_decay:8.3f}  "
          f"{ratio:6.3f}  {K_freq:6d}  {K_corr:7d}  {pq:6.3f}{marker}")


# ================================================================
# 2x DECAY LAW at large N
# ================================================================
print("\n" + "=" * 70)
print("2x DECAY LAW: Paired vs Unpaired Decay Rates")
print("=" * 70)

print(f"\n{'N':>5s}  {'mean_paired':>12s}  {'mean_unpaired':>14s}  {'ratio':>6s}")
print("-" * 45)

for N in [10, 20, 50, 100]:
    n_exc = N // 2
    W, signs = make_balanced_network(N, n_exc, density=0.3, seed=42)
    J, _ = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P=1.5)
    ev = np.linalg.eigvals(J)

    # Classify: oscillatory (paired) vs purely real (unpaired)
    paired_rates = []
    unpaired_rates = []
    for e in ev:
        rate = -e.real
        if abs(e.imag) > 1e-6:
            paired_rates.append(rate)
        else:
            unpaired_rates.append(rate)

    if paired_rates and unpaired_rates:
        mean_p = np.mean(paired_rates)
        mean_u = np.mean(unpaired_rates)
        ratio = mean_u / mean_p if mean_p > 0 else 0
        print(f"  {N:3d}  {mean_p:12.6f}  {mean_u:14.6f}  {ratio:6.3f}")
    else:
        print(f"  {N:3d}  insufficient data (paired={len(paired_rates)}, "
              f"unpaired={len(unpaired_rates)})")


# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
