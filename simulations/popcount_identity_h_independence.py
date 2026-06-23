"""Falsification test for the claim:

    |Re(lambda_k)| = 2*gamma * <popcount(i XOR j)>_{|rho_k[i,j]|^2}

is an H-INDEPENDENT identity, not a model-specific result. It holds because the
Hermitian part of the Liouvillian is the pure Z-dephasing dissipator (diagonal,
-2*gamma*popcount in the coherence basis) and the Hamiltonian only enters the
anti-Hermitian part, which drops out of Re(lambda) via the Rayleigh quotient.

If the claim is right:
  (1) the identity is bit-exact for Hueckel H,
  (2) the identity is bit-exact for a RANDOM Hermitian H (same dissipator),
  (3) Herm(L) = (L + L^dag)/2 is the SAME matrix for both H's.

If any of these fails, my prose was wrong and I owe a correction.
"""
from __future__ import annotations

import numpy as np

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_op(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def site_op(N, site, letter):
    letters = ["I"] * N
    letters[site] = letter
    return pauli_op(letters)


def two_site_op(N, a, b, la, lb):
    letters = ["I"] * N
    letters[a] = la
    letters[b] = lb
    return pauli_op(letters)


def hueckel_ring_H(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "X", "X") + two_site_op(N, a, b, "Y", "Y")
    for l in range(N):
        H = H + 0.5 * site_op(N, l, "Y")
    return H


def random_hermitian(d, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    return A + A.conj().T


def commutator_superop_vec(H):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T))


def lindbladian_vec(H, c_list, gammas):
    L = commutator_superop_vec(H)
    for c, g in zip(c_list, gammas):
        L = L + g * dissipator_superop_vec(c)
    return L


def unvec_column(v, d):
    return v.reshape(d, d, order="F")


def mean_popcount(rho, N):
    d = 2**N
    num = 0.0
    den = 0.0
    for i in range(d):
        for j in range(d):
            w = abs(rho[i, j]) ** 2
            den += w
            num += bin(i ^ j).count("1") * w
    return num / den if den > 1e-15 else 0.0


def check_identity(label, H, N, gamma):
    d = 2**N
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    L = lindbladian_vec(H, c_holstein, [gamma] * N)
    eigvals, eigvecs = np.linalg.eig(L)
    order = np.argsort(eigvals.real)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    print(f"\n[{label}]  N={N}, gamma={gamma}")
    print(f"  {'k':>3}  {'Re(lambda)':>12}  {'2g*<popcount>':>14}  {'abs diff':>12}")
    max_err = 0.0
    for k in range(min(8, d * d)):
        v = eigvecs[:, k]
        v = v / np.linalg.norm(v)
        rho = unvec_column(v, d)
        pred = -2.0 * gamma * mean_popcount(rho, N)
        re = eigvals[k].real
        err = abs(re - pred)
        max_err = max(max_err, err)
        print(f"  {k:>3}  {re:>+12.8f}  {pred:>+14.8f}  {err:>12.2e}")
    print(f"  --> max abs error over 8 slowest modes: {max_err:.3e}")
    return L


def main():
    N = 4
    gamma = 1.0
    d = 2**N

    H_hueckel = hueckel_ring_H(N)
    H_random = random_hermitian(d, seed=12345)

    L1 = check_identity("Hueckel ring + 0.5*sum Y", H_hueckel, N, gamma)
    L2 = check_identity("RANDOM Hermitian H (nonsense model)", H_random, N, gamma)

    herm1 = (L1 + L1.conj().T) / 2
    herm2 = (L2 + L2.conj().T) / 2
    print(f"\nHerm(L) identical for both H's?  max |Herm(L_hueckel) - Herm(L_random)| = "
          f"{np.max(np.abs(herm1 - herm2)):.3e}")

    diss_only = sum(gamma * dissipator_superop_vec(site_op(N, l, "Z")) for l in range(N))
    print(f"Herm(L_hueckel) == pure dissipator?  max |Herm(L) - D| = "
          f"{np.max(np.abs(herm1 - diss_only)):.3e}")

    diag = np.real(np.diag(herm1))
    popcounts = np.array([bin(i ^ j).count("1") for i in range(d) for j in range(d)], dtype=float)
    print(f"Herm(L) diagonal == -2*gamma*popcount(i XOR j)?  max abs = "
          f"{np.max(np.abs(diag - (-2.0 * gamma * popcounts))):.3e}")
    offdiag = herm1 - np.diag(np.diag(herm1))
    print(f"Herm(L) off-diagonal max abs (should be ~0):  {np.max(np.abs(offdiag)):.3e}")


if __name__ == "__main__":
    main()
