"""Polarity probe #7: stress-test the StandardLindbladPiBalance conjecture.

Conjecture (post-probe-6, 2026-05-26):
    Every L of the form
        L = -i (I ⊗ H − H^T ⊗ I) + Σ_k np.kron(c_k, c_k.conj())
    for ANY (Hermitian or non-Hermitian) H and ANY operators c_k satisfies
        ||M_plus_half||² = ||M_minus_half||²
    in the polarity_coordinates_from_L decomposition.

Probes 1-5 confirmed across 5 candidate-breakers (all factored through this form).
Probe 6 broke the balance with a hand-engineered L OUTSIDE this form.

This probe stress-tests the conjecture across 100+ random configurations:
- N in {2, 3, 4}
- Random complex H (Hermitian or not)
- Random number of jump operators (1 to 4)
- Random complex c_k operators
- Various dephase-style and amplitude-damping-style structures

If asymmetry stays bit-exact 0 across all random configs: conjecture is
strongly empirically supported, proof attempt is justified.
If asymmetry ≠ 0 for any config: edge case found, conjecture refines.
"""

import sys
sys.path.insert(0, 'simulations')

import math
import numpy as np

from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L


def random_operator(d, rng, hermitian=False):
    """Random d x d complex operator, optionally Hermitian."""
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    if hermitian:
        A = (A + A.conj().T) / 2
    return A


def build_L_standard_lindblad(H, jump_ops, gammas):
    """Build L = -i (H ⊗ I − I ⊗ H^T) + Σ_k γ_k · np.kron(c_k, c_k.conj()).

    Standard Lindblad form (CP-trace-preserving when gammas > 0 and c_l, c_l^†
    form a balanced pair; here we allow ANY c_k including non-Hermitian, with
    no CP-preservation guarantee, but the structural form is preserved).
    """
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    # H ⊗ I − I ⊗ H^T in numpy kron convention
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(jump_ops, gammas):
        L = L + g * np.kron(c, c.conj())
    return L


def pauli_basis_transform_inv(N):
    """Inverse of VecToPauliBasisTransform: maps Pauli basis -> vec basis.

    Returns T such that L_pauli = T† · L_vec · T / 2^N
    so L_vec = (2^N) · T · L_pauli · T†.

    We avoid explicit construction; instead, since polarity_coordinates_from_L
    accepts L in Pauli basis directly, we'll convert via the framework helper.
    """
    from framework.pauli import _vec_to_pauli_basis_transform
    return _vec_to_pauli_basis_transform(N)


def L_vec_to_pauli(L_vec, N):
    """Transform L from standard vec basis to Pauli basis."""
    T = pauli_basis_transform_inv(N)  # vec -> pauli
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def run_random_config(N, n_jumps, H_hermitian, jump_hermitian, sigma, rng):
    """Build a random Lindblad-form L, decompose, return asymmetry."""
    d = 2 ** N
    H = random_operator(d, rng, hermitian=H_hermitian)
    jump_ops = [random_operator(d, rng, hermitian=jump_hermitian)
                for _ in range(n_jumps)]
    gammas = rng.uniform(0.01, 0.5, size=n_jumps)

    L_vec = build_L_standard_lindblad(H, jump_ops, gammas)
    L_pauli = L_vec_to_pauli(L_vec, N)

    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    return result


def main():
    rng = np.random.default_rng(seed=2026)

    configs = []
    # Sweep N x H-Hermitian x jump-Hermitian x n_jumps
    for N in [2, 3, 4]:
        for H_herm in [True, False]:
            for jump_herm in [True, False]:
                for n_jumps in [1, 2, 3, 4]:
                    configs.append((N, n_jumps, H_herm, jump_herm))

    print(f"Random Lindblad-form sweep: {len(configs)} configs x 5 random seeds each = {len(configs)*5} trials")
    print()
    print(f"{'N':>2}  {'jumps':>5}  {'H-herm':>7}  {'c-herm':>7}  {'seed':>5}  {'||M||^2':>14}  {'asym':>14}  {'asym/||M||^2':>14}")
    print("-" * 110)

    max_relative_asym = 0.0
    broken_configs = []

    for N, n_jumps, H_herm, jump_herm in configs:
        for trial in range(5):
            sub_rng = np.random.default_rng(seed=2026 + trial * 1000 + N * 10000)
            sigma = sub_rng.uniform(0.0, 0.2)
            try:
                result = run_random_config(
                    N=N, n_jumps=n_jumps,
                    H_hermitian=H_herm, jump_hermitian=jump_herm,
                    sigma=sigma, rng=sub_rng
                )
                ns_M = result['norm_sq']['M']
                asym = result['asymmetry']
                rel_asym = abs(asym) / max(ns_M, 1e-15)
                max_relative_asym = max(max_relative_asym, rel_asym)
                marker = " " if rel_asym < 1e-10 else " *** BROKEN ***"
                print(f"{N:>2}  {n_jumps:>5}  {str(H_herm):>7}  {str(jump_herm):>7}  {trial:>5}  "
                      f"{ns_M:>14.4f}  {asym:>+14.4e}  {rel_asym:>14.4e}{marker}")
                if rel_asym > 1e-10:
                    broken_configs.append((N, n_jumps, H_herm, jump_herm, trial, rel_asym))
            except Exception as e:
                msg = repr(str(e)).encode('ascii', 'backslashreplace').decode('ascii')
                print(f"{N:>2}  {n_jumps:>5}  {str(H_herm):>7}  {str(jump_herm):>7}  {trial:>5}  ERROR: {type(e).__name__}: {msg}")

    print()
    print("=" * 110)
    print("Summary:")
    print("=" * 110)
    print(f"  Total trials:                 {len(configs) * 5}")
    print(f"  Max relative asymmetry:       {max_relative_asym:.4e}")
    print(f"  Broken configs (rel > 1e-10): {len(broken_configs)}")
    print()
    if not broken_configs:
        print("  CONJECTURE STRONGLY SUPPORTED: every random Lindblad-form L gives asymmetry ~ 0.")
        print("  The +/-1/2 balance is propagated by the np.kron(c, c.conj()) structure.")
        print("  Next step: attempt the analytic proof.")
    else:
        print("  CONJECTURE BROKEN in some configs:")
        for N, n_jumps, H_herm, jump_herm, trial, rel_asym in broken_configs[:10]:
            print(f"    N={N}, jumps={n_jumps}, H-herm={H_herm}, c-herm={jump_herm}, "
                  f"trial={trial}, rel_asym={rel_asym:.4e}")
        print()
        print("  Conjecture needs refinement.")


if __name__ == '__main__':
    main()
