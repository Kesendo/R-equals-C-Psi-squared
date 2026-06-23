"""Algebra test: does F112 polarity balance hold under Peierls-style bond
dephasing on the Hückel benzene ring, even though F1 spectrum palindrome
breaks (the May 2026 result)?

Prediction from docs/carbon/BENZENE_THREE_DEPHASE_LETTERS.md (2026-05-27):
F112 polarity asymmetry should be 0 bit-exact under Peierls bath because
B_b = X_aX_b + Y_aY_b is bit_b-homogeneous (XX bit_b sum = 0, YY bit_b sum
= 2 mod 2 = 0; both terms in the bit_b = 0 class). F112's hypothesis on
bath operators is "each c_k bit_b-homogeneous", which B_b satisfies even
though it is a composite two-site operator.

This script measures F112 asymmetry on cyclobutadiene C₄ and benzene C₆
rings under both Holstein (Z-dephase per site) and Peierls (D[B_b] per
bond) baths. It also reports the F1 spectrum palindrome residual for
comparison with the May 2026 result.

Honest output: whatever the algebra says.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


LETTERS = ["I", "X", "Z", "Y"]
PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}

# Π_Z per-letter action: letter -> (new_letter, phase)
PI_Z_PER_LETTER = {
    "I": ("X", 1.0 + 0j),
    "X": ("I", 1.0 + 0j),
    "Z": ("Y", 1.0j),
    "Y": ("Z", 1.0j),
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


def bond_op_B(N, a, b):
    """B_b = X_a X_b + Y_a Y_b (Hückel hopping for π-electrons in JW basis)."""
    return site_op(N, a, "X") @ site_op(N, b, "X") + site_op(N, a, "Y") @ site_op(N, b, "Y")


def hueckel_ring_H(N):
    """H = Σ_b B_{a, a+1} for cyclic ring (Hermitian, real)."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + bond_op_B(N, a, b)
    return H


def letters_from_flat(k, N):
    """Decompose flat index k into N letters via base-4 little-endian (a + 2·b packing)."""
    letters = []
    ki = k
    for _ in range(N):
        letters.append(LETTERS[ki & 3])
        ki >>= 2
    return letters


def flat_from_letters(letters):
    """Compose letters back to flat index."""
    k = 0
    for i, L in enumerate(letters):
        k += LETTERS.index(L) << (2 * i)
    return k


def vec_to_pauli_transform(N):
    """T (d² × 4^N) with columns vec_F(σ_k) (column-major flatten)."""
    d = 2**N
    n_basis = 4**N
    T = np.zeros((d * d, n_basis), dtype=complex)
    for k in range(n_basis):
        letters = letters_from_flat(k, N)
        sigma = pauli_op(letters)
        T[:, k] = sigma.flatten("F")
    return T


def build_pi_z(N):
    """Π_Z on the 4^N Pauli basis as a signed permutation matrix.

    Π[newK, k] = phase_product, where the new letters and phases come from
    per-letter PI_Z_PER_LETTER applied to the k-th Pauli string's letters.
    """
    n_basis = 4**N
    pi = np.zeros((n_basis, n_basis), dtype=complex)
    for k in range(n_basis):
        letters_in = letters_from_flat(k, N)
        letters_out = []
        phase = 1.0 + 0j
        for L in letters_in:
            new_L, p = PI_Z_PER_LETTER[L]
            letters_out.append(new_L)
            phase = phase * p
        new_k = flat_from_letters(letters_out)
        pi[new_k, k] = phase
    return pi


def commutator_superop_vec(H):
    """L_vec for -i[H, ·]: (-i)(H ⊗ I − I ⊗ H^T)."""
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    """D_vec[c] in vec basis: c ⊗ c* − ½ (c†c ⊗ I + I ⊗ (c†c)^T)."""
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T))


def lindbladian_in_pauli_basis(H, c_list, gammas, N):
    """L = -i[H, ·] + Σ γ_k D[c_k], expressed in the 4^N Pauli basis."""
    d = 2**N
    L_vec = commutator_superop_vec(H)
    for c, gamma in zip(c_list, gammas):
        L_vec = L_vec + gamma * dissipator_superop_vec(c)
    T = vec_to_pauli_transform(N)
    return T.conj().T @ L_vec @ T / d


def palindrome_residual(L_pauli, pi, sigma_total):
    """M = Π L Π⁻¹ + L + 2σ I (the F1 residual matrix)."""
    n = L_pauli.shape[0]
    pi_inv = pi.conj().T
    return pi @ L_pauli @ pi_inv + L_pauli + 2 * sigma_total * np.eye(n, dtype=complex)


def polarity_decompose(M, pi):
    """Return (M_sym, M_+i, M_-i) Π-eigenspace decomposition of M.

    M_sym = (M + Π M Π⁻¹) / 2   (Π-conj +1, -1 even sector)
    M_anti = (M − Π M Π⁻¹) / 2   (Π-conj +i, -i odd sector)
    M_+i  = (M_anti − i Π M_anti Π⁻¹) / 2   (Π-conj +i eigenspace)
    M_-i  = (M_anti + i Π M_anti Π⁻¹) / 2   (Π-conj -i eigenspace)
    """
    pi_inv = pi.conj().T
    pi_M_pi = pi @ M @ pi_inv
    M_sym = (M + pi_M_pi) / 2
    M_anti = (M - pi_M_pi) / 2
    pi_Manti_pi = pi @ M_anti @ pi_inv
    M_plus_i = (M_anti - 1j * pi_Manti_pi) / 2
    M_minus_i = (M_anti + 1j * pi_Manti_pi) / 2
    return M_sym, M_plus_i, M_minus_i


def f1_spectrum_residual(L_pauli, sigma_total):
    """Max |λ_target − nearest_actual| over all eigenvalues, where
    λ_target = -λ - 2σ for each actual λ. Quantifies F1 palindrome break."""
    eigs = np.linalg.eigvals(L_pauli)
    max_res = 0.0
    for lam in eigs:
        target = -lam - 2 * sigma_total
        nearest_dist = np.min(np.abs(eigs - target))
        max_res = max(max_res, nearest_dist)
    return max_res


def run_test(N, gamma=1.0):
    print()
    print(f"N = {N} ring (cyclic, Hückel H = Σ_b (X_a X_b + Y_a Y_b))")
    print("-" * 72)

    H = hueckel_ring_H(N)
    H_hermitian = np.allclose(H, H.conj().T)
    print(f"  H Hermitian: {H_hermitian}")
    print(f"  H shape: {H.shape} (operator space 4^N = {4**N})")

    pi_z = build_pi_z(N)
    print(f"  Π_Z built: shape {pi_z.shape}, Π² order-2 (Π⁴ = I): {np.allclose(np.linalg.matrix_power(pi_z, 4), np.eye(4**N))}")

    # === HOLSTEIN bath: c_l = Z_l per site (single Pauli, bit_b = 1) ===
    print()
    print(f"  Holstein bath: c_l = Z_l, γ = {gamma} per site (N = {N} terms)")
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    gammas_holstein = [gamma] * N
    sigma_holstein = N * gamma

    L_holstein = lindbladian_in_pauli_basis(H, c_holstein, gammas_holstein, N)
    M_holstein = palindrome_residual(L_holstein, pi_z, sigma_holstein)
    _, M_plus_h, M_minus_h = polarity_decompose(M_holstein, pi_z)
    norm_plus_h = np.linalg.norm(M_plus_h) ** 2
    norm_minus_h = np.linalg.norm(M_minus_h) ** 2
    asym_holstein = norm_plus_h - norm_minus_h
    spec_res_holstein = f1_spectrum_residual(L_holstein, sigma_holstein)

    print(f"    F112 asymmetry  = ‖M_+i‖² − ‖M_−i‖² = {asym_holstein:+.3e}")
    print(f"      ‖M_+i‖²  = {norm_plus_h:.6e}")
    print(f"      ‖M_−i‖²  = {norm_minus_h:.6e}")
    print(f"    F1 spec residual = {spec_res_holstein:.3e}")

    # === PEIERLS bath: c_b = B_b per bond (two-site composite, bit_b = 0) ===
    print()
    print(f"  Peierls bath: c_b = B_b = X_a X_b + Y_a Y_b, γ = {gamma} per bond (N = {N} terms)")
    c_peierls = [bond_op_B(N, a, (a + 1) % N) for a in range(N)]
    gammas_peierls = [gamma] * N
    sigma_peierls = N * gamma  # same convention; F1 spectrum residual likely large regardless

    L_peierls = lindbladian_in_pauli_basis(H, c_peierls, gammas_peierls, N)
    M_peierls = palindrome_residual(L_peierls, pi_z, sigma_peierls)
    _, M_plus_p, M_minus_p = polarity_decompose(M_peierls, pi_z)
    norm_plus_p = np.linalg.norm(M_plus_p) ** 2
    norm_minus_p = np.linalg.norm(M_minus_p) ** 2
    asym_peierls = norm_plus_p - norm_minus_p
    spec_res_peierls = f1_spectrum_residual(L_peierls, sigma_peierls)

    print(f"    F112 asymmetry  = ‖M_+i‖² − ‖M_−i‖² = {asym_peierls:+.3e}")
    print(f"      ‖M_+i‖²  = {norm_plus_p:.6e}")
    print(f"      ‖M_−i‖²  = {norm_minus_p:.6e}")
    print(f"    F1 spec residual = {spec_res_peierls:.3e}")

    return {
        "N": N,
        "asym_holstein": asym_holstein,
        "asym_peierls": asym_peierls,
        "spec_res_holstein": spec_res_holstein,
        "spec_res_peierls": spec_res_peierls,
    }


def main():
    print("=" * 72)
    print("F112 polarity balance on Hückel benzene under Peierls bath")
    print("=" * 72)
    print()
    print("Setup: Hermitian H = Σ_b (X_a X_b + Y_a Y_b) on cyclic N-ring.")
    print("Two baths tested:")
    print("  Holstein (c=Z_l, on-site, bit_b=1): F1+F112 both predicted to hold")
    print("  Peierls  (c=B_b, bond,   bit_b=0): F1 broken (May 2026); F112?")
    print()
    print("F112 predicts asymmetry = 0 bit-exact when each c_k is bit_b-homogeneous.")
    print("B_b = XX + YY: both terms bit_b sum = 0, so B IS bit_b-homogeneous.")
    print("Prediction: F112 asymmetry ≈ 0 for both baths.")
    print()

    results = []
    for N in [4, 6]:
        results.append(run_test(N, gamma=1.0))

    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print()
    print(f"{'N':<4} {'bath':<12} {'F112 asymmetry':<22} {'F1 spec residual':<22}")
    print("-" * 72)
    for r in results:
        N = r["N"]
        print(f"{N:<4} {'Holstein':<12} {r['asym_holstein']:+.3e}            {r['spec_res_holstein']:.3e}")
        print(f"{N:<4} {'Peierls':<12} {r['asym_peierls']:+.3e}            {r['spec_res_peierls']:.3e}")

    print()
    print("Read the asymmetry column. If Peierls F112 asymmetry ≈ 1e-10 or smaller,")
    print("the F112 prediction holds. If it's >> 1e-10, the prediction is wrong and")
    print("the doc needs revision (B as composite bath operator does NOT inherit F112).")


if __name__ == "__main__":
    main()
