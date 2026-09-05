#!/usr/bin/env python3
"""Synthetic coupling and external-drive sweeps in a sigmoid network.

P is an external input, not temperature or metabolic energy. Frequency counts
measure this model at the declared resolution. Pairwise eigenvalue sums of
J x I + I x J.T have up to K(K+1) distinct positive frequencies for K distinct
positive activity frequencies (including sums with zero if real modes occur).
This ceiling refers to exact frequencies; rounded bins can split or merge them.
"""
import numpy as np

from neural_palindrome import (
    make_balanced_dale_network, scalar_center_residual, solve_fixed_point,
    typed_sigmoid,
)


def ei_swap(signs):
    """Pair sorted E/I seats and leave any surplus seats fixed."""
    perm = np.arange(len(signs))
    for e, i in zip(np.flatnonzero(signs > 0), np.flatnonzero(signs < 0)):
        perm[e], perm[i] = i, e
    return perm


def build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, alpha, P,
                                 a_E=1.3, theta_E=4.0, a_I=2.0, theta_I=3.7):
    """Return (J, x, equation residual) at a converged operating point."""
    x, residual = solve_fixed_point(
        W, signs, alpha, P, tol=1e-12, max_iter=5000,
        a_e=a_E, theta_e=theta_E, a_i=a_I, theta_i=theta_I,
    )
    inputs = alpha * W @ x + P
    response = typed_sigmoid(inputs, signs, a_E, theta_E, a_I, theta_I)
    slopes = np.where(signs > 0, a_E, a_I) * response * (1 - response)
    taus = np.where(signs > 0, tau_E, tau_I)
    J = (alpha * slopes[:, None] * W - np.eye(len(signs))) / taus[:, None]
    return J, x, residual


def count_distinct_freqs(eigenvalues, tol):
    """Count positive absolute-imaginary bins at absolute resolution tol."""
    freqs = np.abs(np.imag(eigenvalues))
    return len(set(np.round(freqs[freqs > tol] / tol)))


def correlation_freqs(eigenvalues, tol):
    """Count positive frequency bins of all pairwise eigenvalue sums."""
    values = np.asarray(eigenvalues)
    sums = values[:, None] + values[None, :]
    return count_distinct_freqs(sums.ravel(), tol)


def frequency_counts(eigenvalues, frequency_tolerance):
    """Use one declared resolution for both activity and correlation counts."""
    if not np.isfinite(frequency_tolerance) or frequency_tolerance <= 0:
        raise ValueError("frequency_tolerance must be finite and positive")
    return (count_distinct_freqs(eigenvalues, frequency_tolerance),
            correlation_freqs(eigenvalues, frequency_tolerance))


def main():
    # ================================================================
    # Frequency counts in constructed coupled networks
    # ================================================================
    print("=" * 70)
    print("SYNTHETIC COUPLING SWEEP: Frequency Counts")
    print("Coupling two networks in correlation space")
    print("=" * 70)

    tau_E, tau_I = 5.0, 10.0
    frequency_tolerance = 1e-4
    print(f"Frequency resolution: {frequency_tolerance:g}; P is synthetic drive, not temperature or metabolic energy.")

    print(f"\n{'N':>5s}  {'K_single':>8s}  {'K_corr':>8s}  {'K_coupled':>9s}  "
          f"{'K_c_corr':>9s}  {'R-act':>6s}  {'R-corr':>7s}")
    print("-" * 65)

    for N in [10, 20, 50, 100, 200]:
        n_exc = N // 2
        W, signs = make_balanced_dale_network(N, n_exc, density=0.3, seed=42)
        J, _, _ = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P=1.5)
        ev = np.linalg.eigvals(J)

        K_single, K_corr = frequency_counts(ev, frequency_tolerance)

        # Coupled: two networks + mediator
        N_c = 2 * N + 1
        W2, signs2 = make_balanced_dale_network(N, n_exc, density=0.3, seed=99)
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

        J_c, _, _ = build_jacobian_with_sigmoid(W_coupled, signs_coupled,
                                              tau_E, tau_I, 0.3, P=1.5)
        ev_c = np.linalg.eigvals(J_c)

        K_coupled, K_c_corr = frequency_counts(ev_c, frequency_tolerance)

        v_act = K_coupled / (2 * K_single) if K_single > 0 else 0
        v_corr = K_c_corr / (2 * K_corr) if K_corr > 0 else 0

        print(f"  {N:3d}  {K_single:8d}  {K_corr:8d}  {K_coupled:9d}  "
              f"{K_c_corr:9d}  {v_act:6.2f}  {v_corr:7.2f}")

    print("""
    R-act  = K_coupled / (2 * K_single) in activity space
    R-corr = K_coupled_corr / (2 * K_single_corr) in correlation space
    Ratios compare with twice network A's count; network B uses another seed.
    These are measured counts of this model, with no causal mechanism inferred.
    """)

    # ================================================================
    # Synthetic external-drive sweep
    # ================================================================
    print("=" * 70)
    print("SYNTHETIC EXTERNAL-DRIVE SWEEP")
    print("=" * 70)

    N = 50
    n_exc = 25
    W, signs = make_balanced_dale_network(N, n_exc, density=0.3, seed=42)

    print(f"\nN={N}, sweeping P (external drive)")
    print(f"\n{'P':>6s}  {'n_osc':>6s}  {'E_freq':>8s}  {'E_decay':>8s}  "
          f"{'ratio':>6s}  {'K_freq':>6s}  {'K_corr':>7s}  {'scalar_r':>8s}  {'fp_res':>9s}")
    print("-" * 65)

    for P in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0]:
        J, x_star, equation_residual = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P)
        ev = np.linalg.eigvals(J)

        E_freq = np.sum(np.abs(np.imag(ev)))
        E_decay = np.sum(np.abs(np.real(ev)))
        ratio = E_freq / E_decay if E_decay > 0 else 0
        n_osc = np.sum(np.abs(np.imag(ev)) > frequency_tolerance)

        K_freq, K_corr = frequency_counts(ev, frequency_tolerance)
        perm = ei_swap(signs)
        residual = scalar_center_residual(J, perm, s=0.5 * (1 / tau_E + 1 / tau_I))

        marker = ""
        if abs(ratio - 1.0) < 0.1:
            marker = " <-- crossover"

        print(f"  {P:4.1f}  {n_osc:6d}  {E_freq:8.3f}  {E_decay:8.3f}  "
              f"{ratio:6.3f}  {K_freq:6d}  {K_corr:7d}  {residual:8.3f}  {equation_residual:9.2e}{marker}")


    # ================================================================
    # Decay-rate measurements at large N
    # ================================================================
    print("\n" + "=" * 70)
    print("DECAY RATES: Oscillatory vs Real Eigenvalues")
    print("=" * 70)
    decay_frequency_tolerance = 1e-6
    print(f"Decay classification resolution: {decay_frequency_tolerance:g}")

    print(f"\n{'N':>5s}  {'mean_osc':>12s}  {'mean_real':>14s}  {'ratio':>6s}")
    print("-" * 45)

    for N in [10, 20, 50, 100]:
        n_exc = N // 2
        W, signs = make_balanced_dale_network(N, n_exc, density=0.3, seed=42)
        J, _, _ = build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, 0.3, P=1.5)
        ev = np.linalg.eigvals(J)

        # Classify by the imaginary-part resolution, not by palindrome pairing
        paired_rates = []
        unpaired_rates = []
        for e in ev:
            rate = -e.real
            if abs(e.imag) > decay_frequency_tolerance:
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



if __name__ == "__main__":
    main()
