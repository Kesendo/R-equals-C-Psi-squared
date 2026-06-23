"""Stress-test Step 5 of F112 proof: ||L_H,+i||^2 = ||L_H,-i||^2 for any H?

Step 5 of the F112 proof (if true) says: for any H (Hermitian or not),
the Hamiltonian-only Liouvillian L_H = -i[H, *] has equal Pi-conjugation
+i and -i Frobenius content. Combined with Steps 1-4 (which removed the
bit_b-homogeneous dissipator from contributing to M_+i / M_-i), this
closes the F112 theorem.

Empirical evidence so far: a few non-Hermitian H tests gave bit-exact 0
asymmetry. But the analytic argument (transpose-symmetry of real Pi)
only directly works for real-entry L_H, which Hermitian H gives but
non-Hermitian H does not.

This probe directly tests Step 5: build L_H for many random H (Hermitian
and non-Hermitian), compute ||L_H,+i||^2 and ||L_H,-i||^2 directly via
Pi eigenspace decomposition. Look for cases where they differ.

If they ALWAYS agree: Step 5 is universal (some structural reason beyond
transpose symmetry).
If they SOMETIMES disagree: Step 5 needs an extra condition on H (not
just any H works); the F112 proof needs to identify it.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np

from framework.symmetry import build_pi_full
from framework.pauli import _vec_to_pauli_basis_transform

PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def ps(letters):
    op = PAULI[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI[l])
    return op


def L2P(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def build_L_H(H, N):
    """L_H = -i[H, *] in Pauli basis (no dissipator)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    return L2P(L_vec, N)


def project_pi_eigenspace(M, Pi, target_eigenvalue):
    """Project M onto Pi-conjugation eigenvalue target_eigenvalue eigenspace.

    For Pi of order 4 with eigenvalues {+1, -1, +i, -i}:
    P_lambda(M) = (1/4) sum_{k=0..3} lambda^{-k} Pi^k M Pi^{-k}
    """
    Pi_inv = Pi.conj().T  # Pi unitary
    result = np.zeros_like(M)
    current_Pi = np.eye(Pi.shape[0], dtype=complex)
    current_Pi_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = 1.0 / (target_eigenvalue ** k) / 4.0
        result = result + coef * (current_Pi @ M @ current_Pi_inv)
        current_Pi = current_Pi @ Pi
        current_Pi_inv = current_Pi_inv @ Pi_inv
    return result


def test_H_Step5(H, N, Pi, label):
    """Compute ||L_H,+i||^2 and ||L_H,-i||^2, report asymmetry."""
    L_H = build_L_H(H, N)
    L_H_plus_i = project_pi_eigenspace(L_H, Pi, 1j)
    L_H_minus_i = project_pi_eigenspace(L_H, Pi, -1j)
    norm_plus_i = float(np.sum(np.abs(L_H_plus_i) ** 2))
    norm_minus_i = float(np.sum(np.abs(L_H_minus_i) ** 2))
    asym = norm_plus_i - norm_minus_i
    rel_asym = abs(asym) / max(norm_plus_i + norm_minus_i, 1e-15)
    is_hermitian = np.allclose(H, H.conj().T)
    marker = "BAL" if rel_asym < 1e-10 else "BREAK"
    print(f"  {label:<70}  Herm={is_hermitian!s:<5}  "
          f"||+i||^2={norm_plus_i:>10.4f}  ||-i||^2={norm_minus_i:>10.4f}  "
          f"rel_asym={rel_asym:>10.4e}  [{marker}]")


def main():
    for N in [2, 3]:
        d = 2 ** N
        Pi = build_pi_full(N)
        print("=" * 130)
        print(f"N = {N}")
        print("=" * 130)

        rng = np.random.default_rng(seed=2026 + N)

        # Test 1: Hermitian H from Pauli sums with real coefficients
        print("\n--- Test 1: Hermitian H (real Pauli coefficients) ---\n")
        from itertools import product
        all_paulis = list(product(['I', 'X', 'Y', 'Z'], repeat=N))
        for trial in range(5):
            coeffs = rng.normal(size=len(all_paulis))
            H = sum(coef * ps(p) for coef, p in zip(coeffs, all_paulis))
            test_H_Step5(H, N, Pi, f"random Hermitian H trial {trial}")

        # Test 2: Non-Hermitian H from Pauli sums with complex coefficients
        print("\n--- Test 2: Non-Hermitian H (complex Pauli coefficients) ---\n")
        for trial in range(5):
            coeffs = rng.normal(size=len(all_paulis)) + 1j * rng.normal(size=len(all_paulis))
            H = sum(coef * ps(p) for coef, p in zip(coeffs, all_paulis))
            test_H_Step5(H, N, Pi, f"random non-Hermitian H (complex Pauli coefs) trial {trial}")

        # Test 3: Fully random complex H (NOT a Pauli sum, breaks Pauli structure)
        print("\n--- Test 3: Random complex matrix H (NOT Pauli-structured) ---\n")
        for trial in range(5):
            H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
            test_H_Step5(H, N, Pi, f"random complex matrix H trial {trial}")


if __name__ == '__main__':
    main()
