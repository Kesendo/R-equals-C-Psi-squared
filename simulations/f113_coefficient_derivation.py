"""F113 algebraic derivation: where does the (1/2)·4^N coefficient come from?

The Welle 3 derivation showed bit-exact at N=2, 3, 4:
    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

This Welle-4 script attempts to derive the (4^N / 2) coefficient from
Π-eigenspace algebra. Strategy:

  Step 1: Decompose asymmetry algebraically.
    asymmetry = ‖L_total,+i‖² − ‖L_total,−i‖²
             = (‖L_H,+i‖² − ‖L_H,−i‖²)         ← 0 by F112 typed (Hermitian H)
             + (‖L_T1,+i‖² − ‖L_T1,−i‖²)       ← 0 by F112 non-Herm ext (proven N≤4)
             + 2·Re⟨L_H,+i, L_T1,+i⟩
             − 2·Re⟨L_H,−i, L_T1,−i⟩

  So the entire asymmetry sits in the cross-term:
    asymmetry = 2·(Re⟨L_H,+i, L_T1,+i⟩ − Re⟨L_H,−i, L_T1,−i⟩)

  Step 2: Apply Lemma A (dagger maps Π +i ↔ Π −i isometrically) and
  L_H^† = −L_H (anti-Hermitian for Hermitian H) to relate the +i and −i
  cross-products.

  Step 3: Compute the per-site contribution at N=1 explicitly. The N-scaling
  4^N arises from extending to N qubits by tensoring with I^{⊗(N-1)}, which
  multiplies each Frobenius norm² by 4^(N-1) = d²(N-1).

This script numerically verifies each step at small N to confirm the
decomposition is correct, then identifies the algebraic origin of the
coefficient.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)  # lowering (standard)
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)


def site_op(N, l, m2):
    mats = [I2] * N
    mats[l] = m2
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def build_L_H_vec(H):
    """L_H = -i[H, ·] in vec(ρ) basis."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_L_diss_vec(c, gamma):
    """L_diss = γ·(kron(c, c*) - (1/2)(kron(c†c, I) + kron(I, (c†c)^T))) in vec(ρ) basis."""
    d = c.shape[0]
    Id = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    anti = 0.5 * (np.kron(c_dag_c, Id) + np.kron(Id, c_dag_c.T))
    return gamma * (np.kron(c, c.conj()) - anti)


def to_pauli_basis(L_vec, N):
    T = fw.pauli._vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def project_pi_eigenspace(M, Pi, target_eigenvalue):
    """P_λ(M) = (1/4) Σ_k λ^{-k} Π^k M Π^{-k}."""
    Pi_inv = Pi.conj().T
    result = np.zeros_like(M)
    cur_Pi = np.eye(Pi.shape[0], dtype=complex)
    cur_Pi_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = (1.0 / (target_eigenvalue ** k)) / 4.0
        result = result + coef * (cur_Pi @ M @ cur_Pi_inv)
        cur_Pi = cur_Pi @ Pi
        cur_Pi_inv = cur_Pi_inv @ Pi_inv
    return result


def frobenius_inner(A, B):
    """⟨A, B⟩ = Tr(A† B) = sum(A.conj() * B)."""
    return complex(np.sum(A.conj() * B))


def fnorm_sq(A):
    return float(np.sum(np.abs(A) ** 2))


def asymmetry_decomposition(N, omega, gamma_T1, gamma_pump=0.0):
    """Decompose asymmetry into the four terms predicted algebraically.

    Returns dict with:
      total: ‖L,+i‖² − ‖L,−i‖² (full asymmetry, what F113 predicts)
      LH_alone: ‖L_H,+i‖² − ‖L_H,−i‖² (F112 typed, expect 0)
      LT1_alone: ‖L_T1,+i‖² − ‖L_T1,−i‖² (F112 non-Hermitian ext, expect 0)
      cross_plus: 2 Re⟨L_H,+i, L_T1,+i⟩
      cross_minus: 2 Re⟨L_H,−i, L_T1,−i⟩
      cross_contribution: cross_plus − cross_minus (this should equal total)
      predicted_F113: (4^N / 2) · N · ω · (γ_pump − γ_T1)
    """
    # Build full multi-site H + dissipator, but isolate per-site contribution
    # by using H = (ω/2)·Z_l and c = σ⁻_l on site l only (no other sites driven).
    # The per-site additivity from F113 means total = Σ_l contribution_l.
    Pi = fw.symmetry.build_pi_full(N)

    L_H_vec = np.zeros((4**N, 4**N), dtype=complex)
    L_T1_vec = np.zeros((4**N, 4**N), dtype=complex)
    for l in range(N):
        H_l = (omega / 2.0) * site_op(N, l, Z)
        L_H_vec += build_L_H_vec(H_l)
        c_l_minus = site_op(N, l, SIGMA_MINUS)
        L_T1_vec += build_L_diss_vec(c_l_minus, gamma_T1)
        if gamma_pump != 0:
            c_l_plus = site_op(N, l, SIGMA_PLUS)
            L_T1_vec += build_L_diss_vec(c_l_plus, gamma_pump)

    L_H_pauli = to_pauli_basis(L_H_vec, N)
    L_T1_pauli = to_pauli_basis(L_T1_vec, N)
    L_total_pauli = L_H_pauli + L_T1_pauli

    # M = Π·L·Π⁻¹ + L + 2σ·I  --  the F1 palindrome residual; we want its
    # ±i Π-eigenvalue projections. But for the asymmetry computation we just
    # need ‖L,+i‖² etc.; the constants and shift contribute only to M_zero.
    LH_plus = project_pi_eigenspace(L_H_pauli, Pi, 1j)
    LH_minus = project_pi_eigenspace(L_H_pauli, Pi, -1j)
    LT1_plus = project_pi_eigenspace(L_T1_pauli, Pi, 1j)
    LT1_minus = project_pi_eigenspace(L_T1_pauli, Pi, -1j)
    L_plus = LH_plus + LT1_plus
    L_minus = LH_minus + LT1_minus

    LH_alone = fnorm_sq(LH_plus) - fnorm_sq(LH_minus)
    LT1_alone = fnorm_sq(LT1_plus) - fnorm_sq(LT1_minus)
    cross_plus = 2.0 * float(frobenius_inner(LH_plus, LT1_plus).real)
    cross_minus = 2.0 * float(frobenius_inner(LH_minus, LT1_minus).real)
    total = fnorm_sq(L_plus) - fnorm_sq(L_minus)
    cross_contribution = cross_plus - cross_minus

    predicted = 0.5 * (4**N) * N * omega * (gamma_pump - gamma_T1)

    return {
        'N': N, 'omega': omega, 'gamma_T1': gamma_T1, 'gamma_pump': gamma_pump,
        'total': total,
        'LH_alone': LH_alone,
        'LT1_alone': LT1_alone,
        'cross_plus_2Re': cross_plus,
        'cross_minus_2Re': cross_minus,
        'cross_contribution': cross_contribution,
        'predicted_F113': predicted,
    }


def main():
    print("F113 algebraic decomposition: where does (1/2)·4^N come from?")
    print("=" * 78)
    print()
    print("Algebraic claim: asymmetry has three pieces, and the first two vanish")
    print("by existing F112 results. Numerically verify the decomposition:")
    print()
    for N in [1, 2, 3, 4]:
        print(f"=== N = {N} ===")
        for omega, gt1, gp in [(0.13, 0.001, 0.0), (0.5, 0.01, 0.005)]:
            r = asymmetry_decomposition(N, omega, gt1, gp)
            print(f"  ω={omega:.3f}, γ_T1={gt1:.4f}, γ_pump={gp:.4f}")
            print(f"    total asym               = {r['total']:+.6e}")
            print(f"    L_H-alone asym (F112)    = {r['LH_alone']:+.6e}  (expect 0)")
            print(f"    L_T1-alone asym (F112ne) = {r['LT1_alone']:+.6e}  (expect 0)")
            print(f"    cross_+ − cross_−        = {r['cross_contribution']:+.6e}  (should equal total)")
            print(f"    F113 prediction          = {r['predicted_F113']:+.6e}")
            print()

    # Concentrate on the cross-term structure: for Hermitian H,
    # ⟨L_H,+i, L_T1,+i⟩ vs ⟨L_H,-i, L_T1,-i⟩
    print("=" * 78)
    print("Cross-term structure at N=1 (single-qubit minimal case):")
    print("Using H = (ω/2)·Z, c = σ⁻ (lowering), γ_pump=0.")
    print()
    omega = 1.0
    gamma_T1 = 1.0
    r = asymmetry_decomposition(1, omega, gamma_T1, 0.0)
    print(f"  At ω=γ_T1=1, N=1:")
    print(f"    total asym = {r['total']:+.6f}")
    print(f"    cross_plus = {r['cross_plus_2Re']:+.6f}")
    print(f"    cross_minus = {r['cross_minus_2Re']:+.6f}")
    print(f"    cross diff = {r['cross_contribution']:+.6f}")
    print(f"    F113 predicted = {r['predicted_F113']:+.6f}  (= (1/2)·4·1·(0−1) = -2)")
    print()
    print(f"  Cross-term sign relation: cross_minus = {-r['cross_plus_2Re']:+.6f} ?")
    print(f"  cross_plus_actual = {r['cross_plus_2Re']:+.6f}")
    print(f"  cross_minus_actual = {r['cross_minus_2Re']:+.6f}")
    print(f"  Are they equal-magnitude opposite-sign? "
          f"{abs(r['cross_plus_2Re'] + r['cross_minus_2Re']) < 1e-12}")
    print()
    print("If yes (cross_+ = -cross_-): then cross_+ − cross_− = 2·cross_+ = total asymmetry.")
    print("⇒ asymmetry = 2·cross_+ = 4·Re⟨L_H,+i, L_T1,+i⟩.")


if __name__ == '__main__':
    main()
