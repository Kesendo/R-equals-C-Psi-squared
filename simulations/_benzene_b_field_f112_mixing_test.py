"""B-field mixing test on the carbon ring: does F112 polarity balance survive
when H acquires Π²-odd content via a y-direction Zeeman term?

Setup: Hückel ring H_0 = Σ_b (X_a X_b + Y_a Y_b) on N-site cyclic ring,
plus a tunable y-Zeeman perturbation h · Σ_l Y_l. Holstein dephase bath
c_l = Z_l (bit_b-homogeneous, satisfies F112 hypothesis).

H_0 is Π²-truly (all bilinears bit_b sum = 0 even): L_H_0 lives entirely
in the Π²-even (M_sym) sector, so M_anti = 0 trivially. F112 is vacuous
at h = 0.

The y-Zeeman term σ_y has bit_b sum = 1 (odd), so it has Π²-odd content
and L_(h·Σ Y_l) lives in M_anti. As h grows from 0:
  - ‖M_anti‖² should grow → F112 becomes substantive (non-trivial test)
  - F114 ε(H_0 + h·Σ Y_l) becomes Mixed (Hückel: n_Y even / Zeeman_y:
    n_Y odd → coexisting parity classes per term)
  - F112 itself predicts ‖M_+1/2‖² = ‖M_−1/2‖² bit-exact independent of
    F114 ε classification (F112's hypothesis is H Hermitian + c bit_b-
    homogeneous, both still satisfied)

Sweep h ∈ {0, 0.01, 0.1, 0.5, 1.0, 5.0}. Honest output: whatever the
algebra says.
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

PI_Z_PER_LETTER = {
    "I": ("X", 1.0 + 0j),
    "X": ("I", 1.0 + 0j),
    "Z": ("Y", 1.0j),
    "Y": ("Z", 1.0j),
}


def n_y(letters):
    return sum(1 for L in letters if L == "Y")


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
    return site_op(N, a, "X") @ site_op(N, b, "X") + site_op(N, a, "Y") @ site_op(N, b, "Y")


def hueckel_ring_H(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + bond_op_B(N, a, b)
    return H


def zeeman_y_total(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        H = H + site_op(N, l, "Y")
    return H


def letters_from_flat(k, N):
    out = []
    ki = k
    for _ in range(N):
        out.append(LETTERS[ki & 3])
        ki >>= 2
    return out


def flat_from_letters(letters):
    k = 0
    for i, L in enumerate(letters):
        k += LETTERS.index(L) << (2 * i)
    return k


def vec_to_pauli_transform(N):
    d = 2**N
    n_basis = 4**N
    T = np.zeros((d * d, n_basis), dtype=complex)
    for k in range(n_basis):
        sigma = pauli_op(letters_from_flat(k, N))
        T[:, k] = sigma.flatten("F")
    return T


def build_pi_z(N):
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
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T))


def lindbladian_in_pauli_basis(H, c_list, gammas, N):
    d = 2**N
    L_vec = commutator_superop_vec(H)
    for c, gamma in zip(c_list, gammas):
        L_vec = L_vec + gamma * dissipator_superop_vec(c)
    T = vec_to_pauli_transform(N)
    return T.conj().T @ L_vec @ T / d


def polarity_decompose(M, pi):
    pi_inv = pi.conj().T
    pi_M_pi = pi @ M @ pi_inv
    M_sym = (M + pi_M_pi) / 2
    M_anti = (M - pi_M_pi) / 2
    pi_Manti_pi = pi @ M_anti @ pi_inv
    M_plus = (M_anti - 1j * pi_Manti_pi) / 2
    M_minus = (M_anti + 1j * pi_Manti_pi) / 2
    return M_sym, M_anti, M_plus, M_minus


def palindrome_residual(L_pauli, pi, sigma_total):
    n = L_pauli.shape[0]
    pi_inv = pi.conj().T
    return pi @ L_pauli @ pi_inv + L_pauli + 2 * sigma_total * np.eye(n, dtype=complex)


def run_sweep(N, h_values, gamma=1.0):
    print(f"\nN = {N} ring (cyclic Hückel + h · Σ Y_l Zeeman perturbation)")
    print(f"Bath: Holstein c_l = Z_l, γ = {gamma} per site (σ_total = N·γ = {N*gamma})")
    print("-" * 90)

    pi_z = build_pi_z(N)
    sigma_total = N * gamma
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    gammas = [gamma] * N

    print(f"{'h':<8} {'‖M‖²':<14} {'‖M_sym‖²':<14} {'‖M_anti‖²':<14} {'‖M_+1/2‖²':<14} {'‖M_−1/2‖²':<14} {'asymmetry':<14}")
    print("-" * 90)

    H_0 = hueckel_ring_H(N)
    Y_total = zeeman_y_total(N)

    for h in h_values:
        H_total = H_0 + h * Y_total
        L_pauli = lindbladian_in_pauli_basis(H_total, c_holstein, gammas, N)
        M = palindrome_residual(L_pauli, pi_z, sigma_total)
        _, M_anti, M_plus, M_minus = polarity_decompose(M, pi_z)

        norm_M = np.linalg.norm(M) ** 2
        M_sym_part = (M + pi_z @ M @ pi_z.conj().T) / 2
        norm_M_sym = np.linalg.norm(M_sym_part) ** 2
        norm_M_anti = np.linalg.norm(M_anti) ** 2
        norm_plus = np.linalg.norm(M_plus) ** 2
        norm_minus = np.linalg.norm(M_minus) ** 2
        asym = norm_plus - norm_minus

        print(f"{h:<8.3f} {norm_M:<14.6e} {norm_M_sym:<14.6e} {norm_M_anti:<14.6e} {norm_plus:<14.6e} {norm_minus:<14.6e} {asym:+.3e}")


def f114_classification(N, h):
    """Report F114 ε(H_0 + h·Σ Y_l) per the closed form."""
    # Hückel H_0 = Σ_b (X⊗X + Y⊗Y) at sites a, a+1: per-term n_Y is 0 (XX) or 2 (YY) → all even
    # Zeeman_y term Y_l: n_Y = 1 → odd
    # If h = 0: all terms even → ε = -1
    # If h ≠ 0: terms split (Hückel even, Zeeman odd) → Mixed
    if abs(h) < 1e-12:
        return "ε = -1 (all terms n_Y even)"
    return "ε = Mixed (Hückel n_Y even, Zeeman_y n_Y odd → split parity)"


def main():
    print("=" * 90)
    print("B-field mixing test: F112 polarity balance under H_Hückel + h · Σ Y_l Zeeman")
    print("=" * 90)
    print()
    print("Setup:")
    print("  H_0     = Σ_b (X_a X_b + Y_a Y_b)  Hückel ring (Hermitian, all n_Y even)")
    print("  H_pert  = h · Σ_l Y_l               y-direction Zeeman (Hermitian, n_Y odd per term)")
    print("  H_total = H_0 + H_pert              (Hermitian; F114 ε becomes Mixed for h ≠ 0)")
    print("  c_l     = Z_l on each site          Holstein bath, bit_b-homogeneous")
    print()
    print("F112's hypothesis (Hermitian H + each c_k bit_b-homogeneous) is satisfied for")
    print("all h. F112's prediction: ‖M_+1/2‖² = ‖M_−1/2‖² bit-exact independent of h.")
    print()
    print("F114's prediction: ε(H_0 + h·Σ Y_l) is wohl-definiert nur für h = 0 (clean -1);")
    print("für h ≠ 0 wird ε Mixed (Hückel even ↔ Zeeman odd term-parities split).")
    print()
    print("This sweeps h ∈ {0, 0.01, 0.1, 0.5, 1.0, 5.0} on N = 4 and N = 6 rings.")
    print()

    h_values = [0.0, 0.01, 0.1, 0.5, 1.0, 5.0]

    for N in [4, 6]:
        run_sweep(N, h_values, gamma=1.0)
        print()
        print("F114 ε classification per h:")
        for h in h_values:
            print(f"  h = {h:.3f}:  {f114_classification(N, h)}")
        print()

    print("=" * 90)
    print("READING")
    print("=" * 90)
    print()
    print("At h = 0: H_0 is Π²-truly, L_H lives in M_sym only, M_anti = 0 trivially.")
    print("F112 holds vacuously. ‖M_anti‖² should equal 0 bit-exact.")
    print()
    print("At h > 0: the Zeeman_y term carries Π²-odd content. L_H acquires M_anti")
    print("component, so M_anti grows with h. F112 then becomes a substantive test:")
    print("  if asymmetry ≈ 0 (bit-exact) → F112 confirmed under mixed-ε(F114) Hamiltonian")
    print("  if asymmetry > 0 (significant) → F112 hypothesis is more restrictive than")
    print("    the existing typed claim docstring states (would need follow-up).")
    print()
    print("The algebra answers: check the 'asymmetry' column above against the ‖M_anti‖²")
    print("scale. F112 says they should be in different orders of magnitude.")


if __name__ == "__main__":
    main()
