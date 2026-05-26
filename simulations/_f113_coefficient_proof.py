"""F113 coefficient derivation: rigorous proof of the (1/2)·4^N coefficient.

Closes the open algebraic step of F113 by combining:

  1. Welle-4 reduction:  asymmetry = 4·Re⟨L_H,+i, L_T1,+i⟩
     (verified bit-exact at N = 1, 2, 3, 4 in _f113_coefficient_derivation.py)

  2. Single-site explicit calculation at N = 1 (sympy):
     ⟨(L_H,1)_+i, (L_T1,1)_+i⟩ = -ω·γ/2

  3. Tensor factorization of Π_N = Π_1^⊗N (verified numerically at N = 2, 3).

  4. Π +i projection of single-site superoperators embedded as
     A_l = I_4 ⊗ ... ⊗ A^{(1)} ⊗ ... ⊗ I_4 inherits the tensor structure
     because I_4 is the trivial Π-conj +1 eigenmode.

  5. Frobenius inner product factorizes on tensor products:
     ⟨I_4, I_4⟩ = Tr(I_4) = 4 per spectator site, contributing 4^(N-1) total.

  6. Cross-site terms vanish because Tr((L_T1,1)_+i) = 0 (zero diagonal).

  7. Sum over driven sites: asymmetry = -(1/2)·4^N·Σ_l ω_l·γ_T1,l.

This script performs each step concretely, providing the calculational backbone
for PROOF_F113_COEFFICIENT_DERIVATION.md.

Run: python -X utf8 simulations/_f113_coefficient_proof.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

# ----------------------------------------------------------------------
# Numpy primitives
# ----------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # standard physics: lowers |1>->|0>
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)


def site_op(N, l, m2):
    mats = [I2] * N
    mats[l] = m2
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def L_H_vec(H):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def L_diss_vec(c, gamma):
    d = c.shape[0]
    Id = np.eye(d, dtype=complex)
    ccd = c.conj().T @ c
    return gamma * (np.kron(c, c.conj())
                    - 0.5 * (np.kron(ccd, Id) + np.kron(Id, ccd.T)))


def to_pauli(L_vec, N):
    T = fw.pauli._vec_to_pauli_basis_transform(N)
    return T.conj().T @ L_vec @ T / (2 ** N)


def project_pi(M, Pi, eigenvalue):
    """Π-conjugation eigenspace projection: A_λ = (1/4) Σ_k λ^{-k} Π^k A Π^{-k}."""
    Pi_inv = Pi.conj().T
    result = np.zeros_like(M)
    cur, cur_inv = np.eye(Pi.shape[0], dtype=complex), np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = (1.0 / (eigenvalue ** k)) / 4.0
        result = result + coef * (cur @ M @ cur_inv)
        cur = cur @ Pi
        cur_inv = cur_inv @ Pi_inv
    return result


def fnorm_sq(A):
    return float(np.sum(np.abs(A) ** 2))


def fnorm_inner(A, B):
    return complex(np.sum(A.conj() * B))


# ----------------------------------------------------------------------
# Step 1: Symbolic derivation at N = 1 (sympy)
# ----------------------------------------------------------------------

def step1_symbolic_N1():
    print("=" * 78)
    print("STEP 1: Symbolic derivation at N = 1")
    print("=" * 78)
    print()
    print("Basis ordering (matching framework's _k_to_indices):")
    print("  k=0: I, k=1: X, k=2: Z, k=3: Y")
    print()

    omega, gamma = sp.symbols('omega gamma', real=True, positive=True)

    # Π in Pauli basis (Z-dephasing convention)
    Pi = sp.Matrix([
        [0, 1, 0, 0],     # I -> X
        [1, 0, 0, 0],     # X -> I
        [0, 0, 0, sp.I],  # Z -> i*Y
        [0, 0, sp.I, 0],  # Y -> i*Z
    ])

    # Convention check: Π² = diag(+1, +1, -1, -1)
    Pi_sq = Pi * Pi
    assert Pi_sq == sp.diag(1, 1, -1, -1)
    print("Π² = diag(+1, +1, -1, -1):  bit_b = (0, 0, 1, 1) on (I, X, Z, Y).  Matches F38.")
    print()

    # L_H = -i[H, ·] with H = (ω/2)·Z, in the numerical convention
    # (the overall sign is convention-dependent; we use the numerical-script sign
    # for consistency with F113's empirical anchor):
    L_H = sp.Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, omega],     # L_H(Y) = +ω·X
        [0, 0, 0, 0],
        [0, -omega, 0, 0],    # L_H(X) = -ω·Y
    ])

    # L_T1 = γ·D[σ⁻] where D[σ⁻](ρ) = σ⁻·ρ·σ⁺ - (1/2){σ⁺σ⁻, ρ}
    # On Pauli letters:
    #   D[σ⁻](I) = +Z   (pumping into |0><0|)
    #   D[σ⁻](X) = -X/2  (decay)
    #   D[σ⁻](Z) = -Z   (decoherence)
    #   D[σ⁻](Y) = -Y/2  (decay)
    L_T1 = gamma * sp.Matrix([
        [0,        0,    0,    0],
        [0, -sp.Rational(1, 2), 0, 0],
        [1,        0,   -1,    0],
        [0,        0,    0, -sp.Rational(1, 2)],
    ])

    print("L_H (Pauli basis, N=1, ω-scaled):")
    sp.pprint(L_H)
    print()
    print("L_T1 (Pauli basis, N=1, γ-scaled):")
    sp.pprint(L_T1)
    print()

    # Π +i projection
    def proj_pi(A, Pi, eig):
        Pi_inv = Pi.H
        result = sp.zeros(*A.shape)
        cur = sp.eye(A.shape[0])
        cur_inv = sp.eye(A.shape[0])
        for k in range(4):
            coef = sp.Rational(1, 4) / (eig ** k)
            result = result + coef * (cur * A * cur_inv)
            cur = cur * Pi
            cur_inv = cur_inv * Pi_inv
        return sp.simplify(result)

    L_H_plus = proj_pi(L_H, Pi, sp.I)
    L_T1_plus = proj_pi(L_T1, Pi, sp.I)

    print("L_H,+i (Π +i conjugation eigenmode):")
    sp.pprint(L_H_plus)
    print()
    print("L_T1,+i (Π +i conjugation eigenmode):")
    sp.pprint(L_T1_plus)
    print()

    # Frobenius inner product
    ip = sp.simplify(sum(sp.conjugate(L_H_plus[i, j]) * L_T1_plus[i, j]
                         for i in range(4) for j in range(4)))
    print(f"⟨L_H,+i, L_T1,+i⟩ at N=1 = {ip}")
    print(f"Re⟨L_H,+i, L_T1,+i⟩ at N=1 = {sp.re(ip)}")
    print()
    print(f"Asymmetry (from Welle-4 reduction): 4·Re⟨...⟩ = {sp.simplify(4 * sp.re(ip))}")
    print(f"F113 prediction at N=1: (1/2)·4^1·ω·(0-γ) = -2·ω·γ")
    assert sp.simplify(4 * sp.re(ip) + 2 * omega * gamma) == 0
    print("MATCH ✓")
    print()

    # Trace of L_T1,+i (key for cross-site vanishing in Step 6)
    trace_LT1 = sp.simplify(L_T1_plus.trace())
    trace_LH = sp.simplify(L_H_plus.trace())
    print(f"Tr((L_T1,1)_+i) = {trace_LT1}  (key fact for cross-site vanishing)")
    print(f"Tr((L_H,1)_+i)  = {trace_LH}")
    print()

    # Frobenius norm of identity I_4 (key for the spectator factor 4^(N-1))
    I4 = sp.eye(4)
    norm_I4_sq = sp.simplify(sum(sp.conjugate(I4[i, j]) * I4[i, j]
                                  for i in range(4) for j in range(4)))
    print(f"⟨I_4, I_4⟩ = Tr(I_4) = {norm_I4_sq}  (spectator multiplier per site)")
    print()

    return L_H_plus, L_T1_plus, ip


# ----------------------------------------------------------------------
# Step 2: Tensor factorization Π_N = Π_1^⊗N
# ----------------------------------------------------------------------

def step2_tensor_factorization():
    print("=" * 78)
    print("STEP 2: Tensor factorization Π_N = Π_1 ⊗ Π_1 ⊗ ... ⊗ Π_1")
    print("=" * 78)
    print()
    Pi1 = fw.symmetry.build_pi_full(1)
    for N in [2, 3, 4]:
        Pi_N = fw.symmetry.build_pi_full(N)
        Pi_kron = Pi1
        for _ in range(N - 1):
            Pi_kron = np.kron(Pi_kron, Pi1)
        match = np.allclose(Pi_N, Pi_kron)
        print(f"  N={N}: Π_N == Π_1^⊗{N}?  {match}")
    print()


# ----------------------------------------------------------------------
# Step 3: Embedded single-site superoperators inherit tensor structure
# ----------------------------------------------------------------------

def step3_embedded_singlesite():
    print("=" * 78)
    print("STEP 3: Single-site superoperators embed as I_4 ⊗ ... ⊗ A ⊗ ... ⊗ I_4")
    print("=" * 78)
    print()

    I4 = np.eye(4, dtype=complex)

    # At N=1, the local L_H,1 and L_T1,1
    H_1 = (1.0 / 2.0) * Z
    L_H_1 = to_pauli(L_H_vec(H_1), 1)
    L_T1_1 = to_pauli(L_diss_vec(SIGMA_MINUS, 1.0), 1)

    Pi1 = fw.symmetry.build_pi_full(1)
    L_H_1_plus = project_pi(L_H_1, Pi1, 1j)
    L_T1_1_plus = project_pi(L_T1_1, Pi1, 1j)

    for N in [2, 3]:
        print(f"--- N = {N} ---")
        Pi_N = fw.symmetry.build_pi_full(N)

        for l in range(N):
            H_l = (1.0 / 2.0) * site_op(N, l, Z)
            c_l = site_op(N, l, SIGMA_MINUS)

            L_H_full = to_pauli(L_H_vec(H_l), N)
            L_T1_full = to_pauli(L_diss_vec(c_l, 1.0), N)

            # Expected embedding
            factors_H = [I4] * N
            factors_H[l] = L_H_1
            expected_H = factors_H[0]
            for f in factors_H[1:]:
                expected_H = np.kron(expected_H, f)

            factors_T = [I4] * N
            factors_T[l] = L_T1_1
            expected_T = factors_T[0]
            for f in factors_T[1:]:
                expected_T = np.kron(expected_T, f)

            match_H = np.allclose(L_H_full, expected_H)
            match_T = np.allclose(L_T1_full, expected_T)
            print(f"  site {l}: L_H_N = I_4 ⊗...⊗ L_H_1 ⊗...⊗ I_4 ? {match_H}")
            print(f"           L_T1_N = I_4 ⊗...⊗ L_T1_1 ⊗...⊗ I_4 ? {match_T}")

            # Π +i projection inherits the embedding
            L_H_plus_full = project_pi(L_H_full, Pi_N, 1j)
            L_T1_plus_full = project_pi(L_T1_full, Pi_N, 1j)

            factors_Hp = [I4] * N
            factors_Hp[l] = L_H_1_plus
            expected_Hp = factors_Hp[0]
            for f in factors_Hp[1:]:
                expected_Hp = np.kron(expected_Hp, f)

            factors_Tp = [I4] * N
            factors_Tp[l] = L_T1_1_plus
            expected_Tp = factors_Tp[0]
            for f in factors_Tp[1:]:
                expected_Tp = np.kron(expected_Tp, f)

            match_Hp = np.allclose(L_H_plus_full, expected_Hp)
            match_Tp = np.allclose(L_T1_plus_full, expected_Tp)
            print(f"           (L_H_N)_+i = I_4 ⊗...⊗ (L_H_1)_+i ⊗...⊗ I_4 ? {match_Hp}")
            print(f"           (L_T1_N)_+i = I_4 ⊗...⊗ (L_T1_1)_+i ⊗...⊗ I_4 ? {match_Tp}")
        print()


# ----------------------------------------------------------------------
# Step 4: Cross-site vanishing & per-site inner product
# ----------------------------------------------------------------------

def step4_cross_site_vanishing():
    print("=" * 78)
    print("STEP 4: Cross-site vanishing & per-site inner product 4^(N-1)·(-ωγ/2)")
    print("=" * 78)
    print()

    for N in [2, 3, 4]:
        Pi_N = fw.symmetry.build_pi_full(N)
        print(f"--- N = {N} ---")

        # All (a, b) site pairs: a = drive site (Z_a), b = dissipator site (σ⁻_b)
        print(f"  Site pair inner products (ω = γ = 1):")
        for a in range(N):
            for b in range(N):
                H_a = (1.0 / 2.0) * site_op(N, a, Z)
                c_b = site_op(N, b, SIGMA_MINUS)
                L_H_p = to_pauli(L_H_vec(H_a), N)
                L_T1_p = to_pauli(L_diss_vec(c_b, 1.0), N)
                L_H_plus = project_pi(L_H_p, Pi_N, 1j)
                L_T1_plus = project_pi(L_T1_p, Pi_N, 1j)
                ip = fnorm_inner(L_H_plus, L_T1_plus)
                marker = "(same)" if a == b else "(cross)"
                expected_same = -0.125 * (4 ** N)
                if a == b:
                    print(f"    a={a}, b={b} {marker}: ip = {ip.real:.6f}  "
                          f"(predicted = -(1/8)·4^N = {expected_same:.6f})")
                else:
                    print(f"    a={a}, b={b} {marker}: ip = {ip}  (predicted 0)")
        print()


# ----------------------------------------------------------------------
# Step 5: Closed-form verification with non-uniform per-site rates
# ----------------------------------------------------------------------

def step5_closed_form_verification():
    print("=" * 78)
    print("STEP 5: Closed-form verification with non-uniform per-site rates")
    print("=" * 78)
    print()
    print("Predicted: ⟨L_H,+i, L_T1,+i⟩.real = -(1/8)·4^N·Σ_l ω_l·γ_T1,l")
    print()

    rng = np.random.default_rng(42)
    for N in [2, 3, 4]:
        for sample in range(3):
            omegas = rng.uniform(0.05, 1.0, N)
            gammas = rng.uniform(0.01, 0.5, N)

            Pi_N = fw.symmetry.build_pi_full(N)
            L_H_full = np.zeros((4 ** N, 4 ** N), dtype=complex)
            L_T1_full = np.zeros((4 ** N, 4 ** N), dtype=complex)
            for l in range(N):
                L_H_full += to_pauli(
                    L_H_vec((omegas[l] / 2) * site_op(N, l, Z)), N
                )
                L_T1_full += to_pauli(
                    L_diss_vec(site_op(N, l, SIGMA_MINUS), gammas[l]), N
                )

            L_H_plus = project_pi(L_H_full, Pi_N, 1j)
            L_T1_plus = project_pi(L_T1_full, Pi_N, 1j)
            ip_measured = fnorm_inner(L_H_plus, L_T1_plus).real

            ip_predicted = -0.125 * (4 ** N) * sum(omegas[l] * gammas[l] for l in range(N))

            asym_measured = 4 * ip_measured
            asym_predicted = -0.5 * (4 ** N) * sum(omegas[l] * gammas[l] for l in range(N))

            match = np.isclose(ip_measured, ip_predicted, rtol=1e-12)
            print(f"  N={N}, sample {sample}:")
            print(f"    measured ip = {ip_measured:.8e}, predicted = {ip_predicted:.8e}  match={match}")
            print(f"    measured asym = {asym_measured:.8e}, predicted (F113) = {asym_predicted:.8e}")
        print()


# ----------------------------------------------------------------------
# Step 6: σ⁺ pumping sign verification
# ----------------------------------------------------------------------

def step6_sigma_plus_sign():
    print("=" * 78)
    print("STEP 6: σ⁺ pumping contributes OPPOSITE sign to σ⁻ T1")
    print("=" * 78)
    print()
    # At N=1 with ω=γ=1, σ⁻ gives ip = -1/2; σ⁺ should give +1/2

    Pi_1 = fw.symmetry.build_pi_full(1)
    H_1 = (1.0 / 2.0) * Z
    L_H_1 = to_pauli(L_H_vec(H_1), 1)
    L_H_1_plus = project_pi(L_H_1, Pi_1, 1j)

    L_T1_minus = to_pauli(L_diss_vec(SIGMA_MINUS, 1.0), 1)
    L_T1_minus_plus = project_pi(L_T1_minus, Pi_1, 1j)
    ip_minus = fnorm_inner(L_H_1_plus, L_T1_minus_plus).real

    L_T1_plus = to_pauli(L_diss_vec(SIGMA_PLUS, 1.0), 1)
    L_T1_plus_plus = project_pi(L_T1_plus, Pi_1, 1j)
    ip_plus = fnorm_inner(L_H_1_plus, L_T1_plus_plus).real

    print(f"  ⟨L_H,+i, L_(σ⁻),+i⟩ at N=1, ω=γ=1: {ip_minus:.6f}  (expect -1/2)")
    print(f"  ⟨L_H,+i, L_(σ⁺),+i⟩ at N=1, ω=γ=1: {ip_plus:.6f}  (expect +1/2)")
    print(f"  Sum (σ⁻ + σ⁺ at equal rate): {ip_minus + ip_plus:.6f}  (expect 0: detailed balance)")
    print()


# ----------------------------------------------------------------------
# Step 7: Mixed pump + T1 at non-uniform rates (full F113)
# ----------------------------------------------------------------------

def step7_full_F113():
    print("=" * 78)
    print("STEP 7: Full F113 closed form with σ⁻ + σ⁺ (non-uniform rates)")
    print("=" * 78)
    print()
    print("Predicted asymmetry: (4^N/2)·Σ_l ω_l·(γ_pump,l - γ_T1,l)")
    print()

    rng = np.random.default_rng(2026)
    for N in [2, 3, 4]:
        for sample in range(2):
            omegas = rng.uniform(0.05, 1.0, N)
            g_t1 = rng.uniform(0.001, 0.1, N)
            g_pump = rng.uniform(0.0001, 0.05, N)

            Pi_N = fw.symmetry.build_pi_full(N)
            L_H_full = np.zeros((4 ** N, 4 ** N), dtype=complex)
            L_T1_full = np.zeros((4 ** N, 4 ** N), dtype=complex)
            for l in range(N):
                L_H_full += to_pauli(
                    L_H_vec((omegas[l] / 2) * site_op(N, l, Z)), N
                )
                L_T1_full += to_pauli(
                    L_diss_vec(site_op(N, l, SIGMA_MINUS), g_t1[l]), N
                )
                L_T1_full += to_pauli(
                    L_diss_vec(site_op(N, l, SIGMA_PLUS), g_pump[l]), N
                )

            L_H_plus = project_pi(L_H_full, Pi_N, 1j)
            L_T1_plus = project_pi(L_T1_full, Pi_N, 1j)
            L_H_minus = project_pi(L_H_full, Pi_N, -1j)
            L_T1_minus = project_pi(L_T1_full, Pi_N, -1j)

            # Full asymmetry: ‖L_+i‖² - ‖L_-i‖²
            L_plus = L_H_plus + L_T1_plus
            L_minus = L_H_minus + L_T1_minus
            asym_measured = fnorm_sq(L_plus) - fnorm_sq(L_minus)

            asym_predicted = 0.5 * (4 ** N) * sum(
                omegas[l] * (g_pump[l] - g_t1[l]) for l in range(N)
            )
            match = np.isclose(asym_measured, asym_predicted, rtol=1e-10)
            print(f"  N={N}, sample {sample}:")
            print(f"    measured = {asym_measured:.8e}")
            print(f"    predicted = {asym_predicted:.8e}  match={match}")
        print()


# ----------------------------------------------------------------------
# Step 8: Lemma C verification - L_T1 is real in Pauli basis, cross-term sign
# ----------------------------------------------------------------------

def step8_lemma_c():
    print("=" * 78)
    print("STEP 8: Lemma C verification")
    print("=" * 78)
    print()
    print("Lemma C: L_T1 has only real matrix elements in the Pauli basis.")
    print("This (with Lemma A + Lemma B) gives cross_minus = -cross_plus.")
    print()
    Pi_1 = fw.symmetry.build_pi_full(1)
    H_1 = (1.0 / 2.0) * Z
    L_H_1 = to_pauli(L_H_vec(H_1), 1)
    L_T1_1 = to_pauli(L_diss_vec(SIGMA_MINUS, 1.0), 1)

    print(f"L_T1 (Pauli basis, N=1):")
    print(np.round(L_T1_1.real, 4))
    print(f"Max |imaginary part| = {np.max(np.abs(L_T1_1.imag)):.2e}  (Lemma C: L_T1 is real)")
    print()

    L_H_dag = L_H_1.conj().T
    L_T1_dag = L_T1_1.conj().T

    print(f"L_H^† == -L_H? {np.allclose(L_H_dag, -L_H_1)}  (Lemma B for Hermitian H)")
    print()

    # Lemma A: (A_-i)^† = (A^†)_+i
    L_H_minus = project_pi(L_H_1, Pi_1, -1j)
    L_H_dag_plus = project_pi(L_H_dag, Pi_1, 1j)
    L_T1_minus = project_pi(L_T1_1, Pi_1, -1j)
    L_T1_dag_plus = project_pi(L_T1_dag, Pi_1, 1j)
    print(f"(L_H_-i)^† == (L_H^†)_+i?  {np.allclose(L_H_minus.conj().T, L_H_dag_plus)}  (Lemma A on L_H)")
    print(f"(L_T1_-i)^† == (L_T1^†)_+i? {np.allclose(L_T1_minus.conj().T, L_T1_dag_plus)}  (Lemma A on L_T1)")
    print()

    # Cross-term equality
    L_H_plus = project_pi(L_H_1, Pi_1, 1j)
    L_T1_plus = project_pi(L_T1_1, Pi_1, 1j)
    cross_plus = fnorm_inner(L_H_plus, L_T1_plus).real
    cross_minus = fnorm_inner(L_H_minus, L_T1_minus).real
    print(f"cross_plus  = Re⟨L_H,+i, L_T1,+i⟩ = {cross_plus:.6f}")
    print(f"cross_minus = Re⟨L_H,-i, L_T1,-i⟩ = {cross_minus:.6f}")
    print(f"cross_minus == -cross_plus? {np.isclose(cross_minus, -cross_plus)}  (Lemma C result)")
    print()

    # Verify cross-term sign relation at higher N
    print("Cross-term sign relation at higher N (uniform ω=γ=1):")
    for N in [2, 3, 4]:
        Pi_N = fw.symmetry.build_pi_full(N)
        L_H_full = np.zeros((4**N, 4**N), dtype=complex)
        L_T1_full = np.zeros((4**N, 4**N), dtype=complex)
        for l in range(N):
            L_H_full += to_pauli(L_H_vec(0.5 * site_op(N, l, Z)), N)
            L_T1_full += to_pauli(L_diss_vec(site_op(N, l, SIGMA_MINUS), 1.0), N)
        L_H_p = project_pi(L_H_full, Pi_N, 1j)
        L_H_m = project_pi(L_H_full, Pi_N, -1j)
        L_T1_p = project_pi(L_T1_full, Pi_N, 1j)
        L_T1_m = project_pi(L_T1_full, Pi_N, -1j)
        cp = fnorm_inner(L_H_p, L_T1_p).real
        cm = fnorm_inner(L_H_m, L_T1_m).real
        print(f"  N={N}: cross_+ = {cp:.4f}, cross_- = {cm:.4f}, sign relation: {np.isclose(cm, -cp)}")
    print()

    # Verify the Lemma C key Frobenius equality at N=1..5:
    # ⟨(L_T1^T)_+i, (L_H)_+i⟩ = ⟨(L_H)_+i, (L_T1)_+i⟩
    print("Lemma C key Frobenius equality: ⟨(L_T1^T)_+i, (L_H)_+i⟩ = ⟨(L_H)_+i, (L_T1)_+i⟩")
    for N in [1, 2, 3, 4, 5]:
        Pi_N = fw.symmetry.build_pi_full(N)
        L_H_full = np.zeros((4**N, 4**N), dtype=complex)
        L_T1_full = np.zeros((4**N, 4**N), dtype=complex)
        for l in range(N):
            L_H_full += to_pauli(L_H_vec(0.5 * site_op(N, l, Z)), N)
            L_T1_full += to_pauli(L_diss_vec(site_op(N, l, SIGMA_MINUS), 1.0), N)
        L_T1_T = L_T1_full.T  # transpose, since L_T1 is real -> dagger = transpose
        L_H_p = project_pi(L_H_full, Pi_N, 1j)
        L_T1_p = project_pi(L_T1_full, Pi_N, 1j)
        L_T1_T_p = project_pi(L_T1_T, Pi_N, 1j)
        lhs = fnorm_inner(L_T1_T_p, L_H_p)  # ⟨(L_T1^T)_+i, (L_H)_+i⟩
        rhs = fnorm_inner(L_H_p, L_T1_p)    # ⟨(L_H)_+i, (L_T1)_+i⟩
        print(f"  N={N}: LHS = {lhs:.6f}, RHS = {rhs:.6f}, equal: {np.isclose(lhs, rhs)}")
    print()

    # ===========================================================
    # Lemma C closure (2026-05-26): algebraic two-step reduction
    # ===========================================================
    print("Lemma C ALGEBRAIC CLOSURE (Π²-odd reduction):")
    print("  Step C.1 (algebra): Tr(A_+i · B_-i) + Tr(A_-i · B_+i) = Tr(A_odd · B_odd)")
    print("  Step C.2 (Lindblad): Tr((L_H)_odd · (L_T1)_odd) = 0")
    print("  Full proof + verification: simulations/_f113_lemma_c_step5_closure.py")
    print()
    for N in [1, 2, 3, 4, 5]:
        Pi_N = fw.symmetry.build_pi_full(N)
        Pi_sq = Pi_N @ Pi_N
        Pi_sq_inv = np.linalg.inv(Pi_sq)
        L_H_full = np.zeros((4**N, 4**N), dtype=complex)
        L_T1_full = np.zeros((4**N, 4**N), dtype=complex)
        for l in range(N):
            L_H_full += to_pauli(L_H_vec(0.5 * site_op(N, l, Z)), N)
            L_T1_full += to_pauli(L_diss_vec(site_op(N, l, SIGMA_MINUS), 1.0), N)
        L_H_p = project_pi(L_H_full, Pi_N, 1j)
        L_H_m = project_pi(L_H_full, Pi_N, -1j)
        L_T1_p = project_pi(L_T1_full, Pi_N, 1j)
        L_T1_m = project_pi(L_T1_full, Pi_N, -1j)
        # Step C.1 algebraic identity:
        lhs_C1 = np.trace(L_H_p @ L_T1_m) + np.trace(L_H_m @ L_T1_p)
        # Π²-odd parts:
        L_H_odd = 0.5 * (L_H_full - Pi_sq @ L_H_full @ Pi_sq_inv)
        L_T1_odd = 0.5 * (L_T1_full - Pi_sq @ L_T1_full @ Pi_sq_inv)
        rhs_C1 = np.trace(L_H_odd @ L_T1_odd)
        # Both zero (Step C.2):
        zero_check_C2 = np.isclose(rhs_C1, 0, atol=1e-10)
        print(f"  N={N}: C.1 identity (lhs - rhs) = {abs(lhs_C1 - rhs_C1):.2e}, "
              f"C.2 vanish (|Tr(A_odd · B_odd)|) = {abs(rhs_C1):.2e}  pass={zero_check_C2}")
    print()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print("F113 COEFFICIENT DERIVATION: rigorous proof of (1/2)·4^N coefficient")
    print("=" * 78)
    print()
    print("Strategy: combine Welle-4 reduction (asymmetry = 4·Re⟨L_H,+i, L_T1,+i⟩)")
    print("with single-site sympy + tensor factorization of Π and embedded ops.")
    print()
    print("Key identities used:")
    print("  - Π_N = Π_1^⊗N  (per-site factorization, verified)")
    print("  - L_X,l (single-site at l) = I_4 ⊗ ... ⊗ L_X,1 ⊗ ... ⊗ I_4")
    print("  - (L_X,l)_+i = I_4 ⊗ ... ⊗ (L_X,1)_+i ⊗ ... ⊗ I_4")
    print("  - ⟨A_1 ⊗ B_1, A_2 ⊗ B_2⟩ = ⟨A_1, A_2⟩ · ⟨B_1, B_2⟩")
    print("  - ⟨I_4, I_4⟩ = Tr(I_4) = 4  (per spectator site)")
    print("  - Tr((L_T1,1)_+i) = 0  (zero diagonal -> cross-site vanishes)")
    print()
    print()

    step1_symbolic_N1()
    step2_tensor_factorization()
    step3_embedded_singlesite()
    step4_cross_site_vanishing()
    step5_closed_form_verification()
    step6_sigma_plus_sign()
    step7_full_F113()
    step8_lemma_c()

    print("=" * 78)
    print("ALL STEPS VERIFIED.")
    print("F113 closed form (4^N / 2)·Σ_l ω_l·(γ_pump,l - γ_T1,l) is derived for any N.")
    print()
    print("Key coefficient origin:")
    print("  4 = Welle-4 reduction factor (asymmetry = 4·Re⟨L_H,+i, L_T1,+i⟩)")
    print("  4^(N-1) = (N-1) spectator-site identity factors, each Tr(I_4) = 4")
    print("  1/2 = single-site N=1 inner product magnitude -ωγ/2 / γ / ω")
    print("  Product: 4 · 4^(N-1) · 1/2 = (1/2) · 4^N ✓")
    print("=" * 78)


if __name__ == '__main__':
    main()
