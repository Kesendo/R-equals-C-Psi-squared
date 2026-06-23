"""Algebra test: does F112 polarity balance break when the bath operator
becomes bit_b-mixed?

F112's Tier1Derived hypothesis requires each c_k bit_b-homogeneous (every
Pauli string in c_k shares the same bit_b parity). The previous B-field
mixing test (benzene_b_field_f112_mixing_test.py) confirmed F112 is robust
to F114-Mixed Hamiltonians as long as the bath stays bit_b-homogeneous.

This test inverts the experiment: keep H Hermitian (so F112's H-side
hypothesis holds), but make the bath bit_b-mixed by using amplitude damping
c_l = σ⁻_l = (X_l − i Y_l) / 2. This has X (bit_b=0) and Y (bit_b=1) terms
in the same operator, so c is bit_b-MIXED. F112's hypothesis is now
violated. The algebra answers: does asymmetry become non-zero?

Setup variations:
  (A) H_0 only (Hückel ring) + σ⁻ bath: H is Π²-truly, M_anti=0 trivially
      (vacuous test like the previous Peierls-on-Hückel attempt)
  (B) H_0 + h · Σ Y_l + σ⁻ bath: H has Π²-odd content from Zeeman_y,
      so M_anti ≠ 0, AND bath is bit_b-mixed → genuine F112-violation test

We also include reference cases for comparison:
  (C) H_0 + h · Σ Y_l + Holstein Z bath: F112 hypothesis fully satisfied
      (this is the previous B-field test that confirmed F112 holds bit-exact)
  (D) H_0 + h · Σ Y_l + γ_T1·σ⁻ + γ_Z·Z baths combined: realistic
      "F113-style" setup mixing both bath types

Predicted by F112 docstring (Tier1Derived Hermitian H + each c_k bit_b-
homogeneous): (C) asym = 0 bit-exact; (A), (B), (D) NOT covered by F112's
hypothesis, so asym could be non-zero. Algebra speaks.
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


def zeeman_z_total(N):
    """Z-drive Hamiltonian Σ Z_l: F113 counterexample driver."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        H = H + site_op(N, l, "Z")
    return H


def sigma_minus_site(N, site):
    """σ⁻_l = (X_l − i Y_l) / 2 (amplitude-damping lowering operator).
    bit_b(X) = 0, bit_b(Y) = 1 → bit_b-MIXED."""
    return (site_op(N, site, "X") - 1j * site_op(N, site, "Y")) / 2.0


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
    """D_vec[c] = c ⊗ c* − ½ (c†c ⊗ I + I ⊗ (c†c)^T)."""
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
    M_anti = (M - pi_M_pi) / 2
    pi_Manti_pi = pi @ M_anti @ pi_inv
    M_plus = (M_anti - 1j * pi_Manti_pi) / 2
    M_minus = (M_anti + 1j * pi_Manti_pi) / 2
    return M_anti, M_plus, M_minus


def palindrome_residual(L_pauli, pi, sigma_total):
    n = L_pauli.shape[0]
    pi_inv = pi.conj().T
    return pi @ L_pauli @ pi_inv + L_pauli + 2 * sigma_total * np.eye(n, dtype=complex)


def measure_f112(N, H, c_list, gammas, sigma_total, pi_z):
    L = lindbladian_in_pauli_basis(H, c_list, gammas, N)
    M = palindrome_residual(L, pi_z, sigma_total)
    M_anti, M_plus, M_minus = polarity_decompose(M, pi_z)
    norm_M = np.linalg.norm(M) ** 2
    norm_M_anti = np.linalg.norm(M_anti) ** 2
    norm_plus = np.linalg.norm(M_plus) ** 2
    norm_minus = np.linalg.norm(M_minus) ** 2
    asym = norm_plus - norm_minus
    rel_asym = abs(asym) / max(norm_M, 1e-15)
    return norm_M, norm_M_anti, norm_plus, norm_minus, asym, rel_asym


def run_test(N, h_values, gamma_z=1.0, gamma_t1=1.0):
    print()
    print(f"N = {N} ring (cyclic Hückel + h · Σ Y_l Zeeman)")
    print(f"  Bath options: Holstein γ_Z = {gamma_z}; amplitude damping γ_T1 = {gamma_t1}")
    print("-" * 100)

    pi_z = build_pi_z(N)
    H_0 = hueckel_ring_H(N)
    Y_total = zeeman_y_total(N)

    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    c_amp_damp = [sigma_minus_site(N, l) for l in range(N)]
    sigma_holstein = N * gamma_z
    sigma_amp_damp = N * gamma_t1  # rough convention; F1 centring is offset, only F112 matters here

    print(f"{'h':<6} {'bath':<24} {'‖M_anti‖²':<14} {'‖M_+1/2‖²':<14} {'‖M_−1/2‖²':<14} {'asymmetry':<14} {'rel asym':<12}")
    print("-" * 100)

    for h in h_values:
        H_total = H_0 + h * Y_total

        # (C) Holstein only (F112 hypothesis fully satisfied — reference)
        results_c = measure_f112(N, H_total, c_holstein, [gamma_z] * N, sigma_holstein, pi_z)
        print(f"{h:<6.3f} {'Holstein Z':<24} {results_c[1]:<14.6e} {results_c[2]:<14.6e} {results_c[3]:<14.6e} {results_c[4]:+.3e}     {results_c[5]:<12.3e}")

        # (B) Amplitude damping only (bit_b-mixed bath — F112 hypothesis violated)
        results_b = measure_f112(N, H_total, c_amp_damp, [gamma_t1] * N, sigma_amp_damp, pi_z)
        print(f"{h:<6.3f} {'amplitude damp σ⁻':<24} {results_b[1]:<14.6e} {results_b[2]:<14.6e} {results_b[3]:<14.6e} {results_b[4]:+.3e}     {results_b[5]:<12.3e}")

        # (D) Both baths combined (realistic noisy carbon system)
        c_combined = c_holstein + c_amp_damp
        gammas_combined = [gamma_z] * N + [gamma_t1] * N
        sigma_combined = sigma_holstein + sigma_amp_damp
        results_d = measure_f112(N, H_total, c_combined, gammas_combined, sigma_combined, pi_z)
        print(f"{h:<6.3f} {'Holstein + σ⁻':<24} {results_d[1]:<14.6e} {results_d[2]:<14.6e} {results_d[3]:<14.6e} {results_d[4]:+.3e}     {results_d[5]:<12.3e}")
        print()


def main():
    print("=" * 100)
    print("bit_b-mixed bath test on Hückel carbon ring")
    print("=" * 100)
    print()
    print("F112 hypothesis: H Hermitian AND each c_k bit_b-homogeneous (all Pauli strings in c_k share parity)")
    print("  c_l = Z_l per site         : bit_b = 1 single Pauli, homogeneous ✓")
    print("  c_l = σ⁻_l = (X − iY) / 2  : bit_b(X) = 0, bit_b(Y) = 1, MIXED ✗")
    print()
    print("F112 predicts asym = 0 bit-exact ONLY for the homogeneous case.")
    print("For mixed bath the algebra answers: does asym become non-zero?")
    print()
    print("Three configurations per h:")
    print("  Holstein Z (reference)         : F112 hypothesis satisfied → asym should be 0")
    print("  amplitude damping σ⁻           : F112 hypothesis violated → asym free to grow")
    print("  Holstein + σ⁻ combined         : F112 partially violated → asym from σ⁻ only?")
    print()

    h_values = [0.0, 0.1, 1.0]

    for N in [4, 6]:
        run_test(N, h_values, gamma_z=1.0, gamma_t1=1.0)

    print()
    print("=" * 100)
    print("F113 control test: Z-drive Hamiltonian + σ⁻ bath (the known F112-break setup)")
    print("=" * 100)
    print()
    print("F113 (Welle 4, 2026-05-26) predicts closed-form asymmetry for this exact setup:")
    print("    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)")
    print("With ω = 1.0 on each site, γ_pump = 0, γ_T1 = 1.0 → asymmetry = (4^N / 2) · N · (-1)")
    print("                                                              = -N · 4^N / 2")
    print()

    for N in [2, 3]:  # small N for clear F113 magnitudes
        print(f"N = {N}, H = Σ Z_l (ω = 1), c = σ⁻ at γ_T1 = 1, no Z-dephase:")
        H_zdrive = zeeman_z_total(N)
        c_amp_damp = [sigma_minus_site(N, l) for l in range(N)]
        gammas_t1 = [1.0] * N
        # Heuristic sigma for the F1 residual centring (best-fit; F112 asym is centre-independent)
        sigma_eff = N * 0.5  # σ⁻'s effective dephase contribution
        pi_z = build_pi_z(N)
        results = measure_f112(N, H_zdrive, c_amp_damp, gammas_t1, sigma_eff, pi_z)
        predicted_asym_f113 = -(4**N / 2) * N * 1.0  # ω = 1, (γ_pump − γ_T1) = -1, sum over N sites
        print(f"  ‖M_anti‖²    = {results[1]:.6e}")
        print(f"  ‖M_+1/2‖²    = {results[2]:.6e}")
        print(f"  ‖M_−1/2‖²    = {results[3]:.6e}")
        print(f"  asymmetry    = {results[4]:+.6e}")
        print(f"  F113 predict = {predicted_asym_f113:+.6e}  [from (4^N/2)·Σ ω·(γ_pump - γ_T1)]")
        print(f"  ratio        = {results[4] / predicted_asym_f113:.4f}  (1.0 = exact match)")
        print()

    print()
    print("=" * 100)
    print("READING")
    print("=" * 100)
    print()
    print("rel asym = |‖M_+‖² − ‖M_−‖²| / ‖M‖² is the dimensionless F112-break magnitude.")
    print("  < 1e-12 : F112 holds bit-exact (homogeneous bath confirmation)")
    print("  > 1e-6  : F112 broken substantively (mixed bath / out-of-scope)")
    print()
    print("If σ⁻ rows show rel asym >> Holstein rows, F112 hypothesis is the right gate:")
    print("crossing into bit_b-mixed c regime substantively breaks the polarity balance.")
    print("If σ⁻ rows also show rel asym ≈ 0, F112 is broader than its docstring states.")


if __name__ == "__main__":
    main()
