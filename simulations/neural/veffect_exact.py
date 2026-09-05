#!/usr/bin/env python3
"""Coupling and drive sweeps on constructed networks.

The balanced linear construction satisfies the scalar-center identity. Adding
a fixed excitatory seat does not, including at zero coupling. Frequency counts
are measurements of these specific matrices; no V-effect mechanism is claimed.
"""
import numpy as np

from neural_palindrome import build_exact_weights, build_linear_jacobian, scalar_center_residual
from veffect_and_heat import build_jacobian_with_sigmoid, ei_swap, frequency_counts


def build_exact_palindromic_network(N, n_exc, tau_E, tau_I, density=0.3, seed=42):
    """Compatibility builder; only balanced E/I inputs have the scalar identity."""
    return build_exact_weights(N, n_exc, tau_E, tau_I, density, seed)


def count_freqs(eigenvalues, tol=1e-6):
    """Compatibility activity count at the declared absolute resolution."""
    return frequency_counts(eigenvalues, tol)[0]


def count_corr_freqs(eigenvalues, tol=1e-6):
    """Compatibility correlation count at the declared absolute resolution."""
    return frequency_counts(eigenvalues, tol)[1]


def palindrome_residual(J, signs):
    """Compatibility residual for balanced linear matrices with two leak values.

    This is not a general gate. Infer s only when each sign class has exactly
    one finite real diagonal value and the E/I swap has no fixed seats.
    Odd populations or nonuniform type diagonals require an explicit center
    through scalar_center_residual(J, perm, s).
    """
    J, signs = np.asarray(J), np.asarray(signs)
    message = "use scalar_center_residual(J, perm, s) with an explicit center"
    if (signs.ndim != 1 or signs.size == 0 or J.shape != (signs.size, signs.size)
            or not np.all(np.isin(signs, [-1, 1]))
            or np.count_nonzero(signs > 0) != np.count_nonzero(signs < 0)):
        raise ValueError(message)
    perm = ei_swap(signs)
    diagonal = np.diag(J)
    e_leaks, i_leaks = diagonal[signs > 0], diagonal[signs < 0]
    if (np.any(perm == np.arange(signs.size)) or not np.isrealobj(diagonal)
            or not np.all(np.isfinite(diagonal))
            or not np.all(e_leaks == e_leaks[0]) or not np.all(i_leaks == i_leaks[0])):
        raise ValueError(message)
    s = -0.5 * (e_leaks[0] + i_leaks[0])
    return scalar_center_residual(J, perm, s)


def main():
    # ================================================================
    # Step 1: Verify exact palindrome for single network
    # ================================================================
    print("=" * 70)
    print("STEP 1: Constructed balanced networks: scalar-center residual")
    print("=" * 70)

    tau_E, tau_I = 5.0, 10.0
    alpha = 0.5
    s = 0.5 * (1 / tau_E + 1 / tau_I)
    frequency_tolerance = 1e-6
    print(f"Frequency resolution: {frequency_tolerance:g}; s={s:g}")

    for N in [10, 20, 30]:
        n_exc = N // 2
        W, signs, perm = build_exact_weights(N, n_exc, tau_E, tau_I, seed=42)
        J = build_linear_jacobian(W, signs, tau_E, tau_I, alpha)
        res = scalar_center_residual(J, perm, s)
        ev = np.linalg.eigvals(J)
        K, K_corr = frequency_counts(ev, frequency_tolerance)

        print(f"\n  N={N}: scalar-center residual = {res:.2e}  "
              f"K_act={K}  K_corr={K_corr}")


    # ================================================================
    # Step 2: Couple two exact networks - does the palindrome break?
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Coupling sweep with an extra excitatory seat")
    print("The odd-seat system is not an exact scalar-center palindrome, even at zero coupling.")
    print("=" * 70)

    for N in [10, 20]:
        n_exc = N // 2

        # Network A
        W_A, signs_A, _ = build_exact_weights(N, n_exc, tau_E, tau_I, seed=42)
        J_A = build_linear_jacobian(W_A, signs_A, tau_E, tau_I, alpha)
        ev_A = np.linalg.eigvals(J_A)
        res_A = scalar_center_residual(J_A, ei_swap(signs_A), s)

        # Network B
        W_B, signs_B, _ = build_exact_weights(N, n_exc, tau_E, tau_I, seed=99)
        J_B = build_linear_jacobian(W_B, signs_B, tau_E, tau_I, alpha)
        ev_B = np.linalg.eigvals(J_B)
        res_B = scalar_center_residual(J_B, ei_swap(signs_B), s)

        K_A, K_A_corr = frequency_counts(ev_A, frequency_tolerance)
        K_B, K_B_corr = frequency_counts(ev_B, frequency_tolerance)

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
            res_c = scalar_center_residual(J_c, ei_swap(signs_c), s)
            K_c, K_c_corr = frequency_counts(ev_c, frequency_tolerance)

            # Count differences relative to the sum of the two single-network counts
            new_act = K_c - (K_A + K_B)
            new_corr = K_c_corr - (K_A_corr + K_B_corr)

            if coupling_strength == 0.0:
                print(f"\n  N={N}: single A: K={K_A}, K_corr={K_A_corr}, res={res_A:.2e}")
                print(f"  Coupling  Residual    K_act  K_corr  New_act  New_corr  ")
                print(f"  {'-'*65}")

            print(f"    {coupling_strength:5.2f}    {res_c:.2e}   {K_c:5d}  {K_c_corr:7d}  "
                  f"{new_act:+7d}  {new_corr:+8d}")


    # ================================================================
    # Step 3: Synthetic drive on the constructed network
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Synthetic external-drive sweep")
    print("P is an external input; no thermal or V-effect mechanism is inferred.")
    print("=" * 70)

    N = 20
    n_exc = N // 2
    W, signs, _ = build_exact_weights(N, n_exc, tau_E, tau_I, seed=42)

    # Build Jacobian WITH sigmoid (nonlinear, P-dependent)
    print(f"\n  N={N}, exact palindromic at linear level")
    print(f"\n  {'P':>6s}  {'n_osc':>6s}  {'K_freq':>6s}  {'K_corr':>7s}  "
          f"{'delta_K':>7s}  {'scalar_r':>10s}  {'fp_res':>9s}")
    print(f"  {'-'*55}")

    prev_K = None
    for P in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]:
        J_nl, x, equation_residual = build_jacobian_with_sigmoid(
            W, signs, tau_E, tau_I, alpha, P)
        ev_nl = np.linalg.eigvals(J_nl)
        n_osc = np.sum(np.abs(np.imag(ev_nl)) > frequency_tolerance)
        K, K_corr = frequency_counts(ev_nl, frequency_tolerance)
        res = scalar_center_residual(J_nl, ei_swap(signs), s)

        delta = K_corr - prev_K if prev_K is not None else 0
        prev_K = K_corr

        marker = ""
        if delta == 2:
            marker = " <-- +2 NEW"
        elif delta > 0:
            marker = f" <-- +{delta} new"

        print(f"  {P:4.1f}  {n_osc:6d}  {K:6d}  {K_corr:7d}  {delta:+7d}  "
              f"{res:10.2e}  {equation_residual:9.2e}{marker}")



if __name__ == "__main__":
    main()
