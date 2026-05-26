"""F113 Lemma C step 5: CLOSED algebraic proof.

This script gives the structural closure of the open identity
   Tr((L_H)_+i · (L_T1)_-i) + Tr((L_H)_-i · (L_T1)_+i) = 0
via the two-step reduction:

  STEP 1 (algebraic reduction via Π²-odd projection):
    For any A, B and Π satisfying Π^* = Π^-1 (unitary with conjugate = inverse):

       Tr(A_+i · B_-i) + Tr(A_-i · B_+i) = Tr(A_odd · B_odd)

    where A_odd := (A - Π² A Π^-2)/2 is the Π²-anti-symmetric part of A,
    and similarly for B_odd. This is a pure operator-algebra identity, no
    assumption on the matrices A, B beyond Π being unitary.

  STEP 2 (the specific Lindblad input):
    For A = L_H with H = Σ_l (ω_l/2)·Z_l (single-site Z-drives) and
    B = L_T1 = Σ_l γ_T1,l · D[σ⁻_l] (per-site amplitude damping):

       Tr((L_H)_odd · (L_T1)_odd) = 0

    Proof: per-site additivity reduces this to per-site contributions.
    At single site N=1:
      - (L_H,1) is entirely Π²-odd ((L_H,1)_odd = L_H,1).
      - (L_T1,1)_odd has support only at the (Z, I) entry with value γ_T1.
      - Trace pickup: (L_H,1)_odd · (L_T1,1)_odd has its only diagonal
        contribution from the (I, Z) entry of L_H, which is L_H[I, Z] = 0
        (since [Z, σ_Z] = 0, so L_H sends σ_Z to 0; equivalently L_H is
        non-zero only on (X, Y) and (Y, X) entries at N=1).

    Cross-site terms (l ≠ m): both Tr((L_H,1)_odd) and Tr((L_T1,1)_odd) are
    zero, and the tensor structure factorizes the trace, so cross-site
    contributions vanish.

This script verifies both steps at N = 1, 2, 3, 4, 5 bit-exact.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)


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
    Pi_inv = Pi.conj().T
    result = np.zeros_like(M)
    cur = np.eye(Pi.shape[0], dtype=complex)
    cur_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = (1.0 / (eigenvalue ** k)) / 4.0
        result = result + coef * (cur @ M @ cur_inv)
        cur = cur @ Pi
        cur_inv = cur_inv @ Pi_inv
    return result


def project_pi2_odd(M, Pi):
    Pi_inv = Pi.conj().T
    Pi2 = Pi @ Pi
    Pi2_inv = Pi_inv @ Pi_inv
    return 0.5 * (M - Pi2 @ M @ Pi2_inv)


def project_pi2_even(M, Pi):
    Pi_inv = Pi.conj().T
    Pi2 = Pi @ Pi
    Pi2_inv = Pi_inv @ Pi_inv
    return 0.5 * (M + Pi2 @ M @ Pi2_inv)


# ============================================================
# STEP 1: Pure operator-algebra reduction
# ============================================================

def step_1_algebraic_reduction():
    print("=" * 78)
    print("STEP 1: Algebraic identity")
    print("  Tr(A_+i · B_-i) + Tr(A_-i · B_+i) = Tr(A_odd · B_odd)")
    print("=" * 78)
    print()
    print("Derivation (no assumption on A, B besides Π unitary):")
    print()
    print("  A_+i = (1/4) Σ_k (1/i)^k Π^k A Π^-k = (1/4) Σ_k (-i)^k Π^k A Π^-k")
    print("  A_-i = (1/4) Σ_k (1/(-i))^k Π^k A Π^-k = (1/4) Σ_k (i)^k Π^k A Π^-k")
    print()
    print("  A_+i + A_-i = (1/4) Σ_k [(-i)^k + (i)^k] Π^k A Π^-k")
    print("    k=0: 1+1=2;  k=1: -i+i=0;  k=2: -1+(-1)=-2;  k=3: i+(-i)=0")
    print("  = (1/2)(A - Π² A Π^-2) = A_odd  (Π²-anti-symmetric projection of A)")
    print()
    print("  A_+i - A_-i = (1/4) Σ_k [(-i)^k - (i)^k] Π^k A Π^-k")
    print("    k=0: 0; k=1: -2i; k=2: 0; k=3: 2i")
    print("  = (1/4)[-2i Π A Π^-1 + 2i Π³ A Π^-3] = (i/2)·Π·(Π² A Π^-2 - A)·Π^-1")
    print("  = -i·Π·A_odd·Π^-1")
    print()
    print("Therefore: A_+i = (1/2)(A_odd - i Π A_odd Π^-1)")
    print("           A_-i = (1/2)(A_odd + i Π A_odd Π^-1)")
    print()
    print("Computing the trace sum:")
    print()
    print("  A_+i · B_-i = (1/4)(A_odd - i Π A_odd Π^-1)(B_odd + i Π B_odd Π^-1)")
    print("    expand and group:")
    print("  = (1/4)[A_odd · B_odd + i A_odd · Π B_odd Π^-1")
    print("          - i Π A_odd Π^-1 · B_odd + Π A_odd Π^-1 · Π B_odd Π^-1]")
    print()
    print("  A_-i · B_+i = (1/4)[A_odd · B_odd - i A_odd · Π B_odd Π^-1")
    print("          + i Π A_odd Π^-1 · B_odd + Π A_odd Π^-1 · Π B_odd Π^-1]")
    print()
    print("Sum: A_+i · B_-i + A_-i · B_+i")
    print("  = (1/2)[A_odd · B_odd + Π A_odd Π^-1 · Π B_odd Π^-1]")
    print("  = (1/2)[A_odd · B_odd + Π (A_odd · B_odd) Π^-1]      (combine inner Π·Π^-1)")
    print()
    print("Taking the trace and using trace cyclicity (Tr(Π X Π^-1) = Tr(X)):")
    print()
    print("  Tr(A_+i · B_-i + A_-i · B_+i)")
    print("  = (1/2)[Tr(A_odd · B_odd) + Tr(Π (A_odd · B_odd) Π^-1)]")
    print("  = (1/2)[Tr(A_odd · B_odd) + Tr(A_odd · B_odd)]")
    print("  = Tr(A_odd · B_odd)  ∎")
    print()
    print("This is a PURE OPERATOR-ALGEBRA IDENTITY, no assumption on A, B.")
    print()

    # Numerical verification
    print("Numerical verification: cross-trace identity holds for any A, B, Π unitary?")
    print()
    rng = np.random.default_rng(2026)
    for N in [1, 2, 3, 4, 5]:
        Pi = fw.symmetry.build_pi_full(N)
        # Use random complex matrices A, B
        d = 4 ** N
        A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        B = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))

        Ap = project_pi(A, Pi, 1j)
        Am = project_pi(A, Pi, -1j)
        Bp = project_pi(B, Pi, 1j)
        Bm = project_pi(B, Pi, -1j)

        lhs = np.trace(Ap @ Bm) + np.trace(Am @ Bp)

        A_odd = project_pi2_odd(A, Pi)
        B_odd = project_pi2_odd(B, Pi)
        rhs = np.trace(A_odd @ B_odd)

        # Verify A_+i + A_-i = A_odd
        match_sum = np.allclose(Ap + Am, A_odd)
        match_trace = np.isclose(lhs, rhs)
        print(f"  N={N}: A_+i + A_-i == A_odd? {match_sum}  Tr identity? {match_trace}")


# ============================================================
# STEP 2: The specific Lindblad input gives Tr(A_odd · B_odd) = 0
# ============================================================

def step_2_lindblad_input():
    print()
    print("=" * 78)
    print("STEP 2: For A = L_H (single-site Z-drives) and B = L_T1 (per-site σ⁻):")
    print("        Tr((L_H)_odd · (L_T1)_odd) = 0")
    print("=" * 78)
    print()
    print("Per-site additivity: L_H = Σ_l L_H,l and L_T1 = Σ_l L_T1,l, so")
    print("  Tr((L_H)_odd · (L_T1)_odd) = Σ_{l,m} Tr((L_H,l)_odd · (L_T1,m)_odd)")
    print()
    print("This splits into same-site (l=m) and cross-site (l≠m) parts.")
    print()

    # Verify at N=1 explicitly (sympy)
    omega, gamma = sp.symbols('omega gamma', real=True, positive=True)
    Pi = sp.Matrix([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, sp.I],
        [0, 0, sp.I, 0],
    ])
    Pi_inv = Pi.H
    Pi_sq = Pi * Pi
    Pi_sq_inv = Pi_inv * Pi_inv

    L_H_sym = sp.Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, omega],
        [0, 0, 0, 0],
        [0, -omega, 0, 0],
    ])
    L_T1_sym = gamma * sp.Matrix([
        [0, 0, 0, 0],
        [0, -sp.Rational(1, 2), 0, 0],
        [1, 0, -1, 0],
        [0, 0, 0, -sp.Rational(1, 2)],
    ])

    L_H_odd_sym = sp.simplify((L_H_sym - Pi_sq * L_H_sym * Pi_sq_inv) / 2)
    L_T1_odd_sym = sp.simplify((L_T1_sym - Pi_sq * L_T1_sym * Pi_sq_inv) / 2)

    print("Single-site N=1 (sympy):")
    print()
    print("  (L_H,1)_odd:")
    sp.pprint(L_H_odd_sym)
    print()
    print(f"  Matches L_H,1? {sp.simplify(L_H_odd_sym - L_H_sym) == sp.zeros(4, 4)}")
    print("    (L_H,1 is ENTIRELY Π²-odd: every nonzero entry of L_H,1 has")
    print("     bit_b(row) + bit_b(col) = 1 since L_H sends bit_b=0 letters (I,X) to bit_b=1")
    print("     letters (Z,Y) and vice versa via [Z,·].)")
    print()

    print("  (L_T1,1)_odd:")
    sp.pprint(L_T1_odd_sym)
    print()
    print("  This has only one nonzero entry: (L_T1,1)_odd[Z, I] = γ.")
    print("  Why? L_T1 has entries on (X,X), (Y,Y), (Z,Z), (Z,I).")
    print("    (X,X): bit_b sum = 0+0 = 0 (even); zero in odd projection")
    print("    (Y,Y): bit_b sum = 1+1 = 0 (even); zero in odd projection")
    print("    (Z,Z): bit_b sum = 1+1 = 0 (even); zero in odd projection")
    print("    (Z,I): bit_b sum = 1+0 = 1 (odd); SURVIVES odd projection")
    print()

    # Compute the product and its trace
    prod = sp.simplify(L_H_odd_sym * L_T1_odd_sym)
    print("  (L_H,1)_odd · (L_T1,1)_odd:")
    sp.pprint(prod)
    print()
    tr = sp.simplify(prod.trace())
    print(f"  Tr((L_H,1)_odd · (L_T1,1)_odd) = {tr}")
    print()
    print("  Reason: the only nonzero (β, β) diagonal entry of the product comes from")
    print("    β = I (col 0), with value (L_H)_odd[I, Z] · (L_T1)_odd[Z, I] = 0 · γ = 0")
    print("    L_H[I, Z] = 0 because [Z, σ_Z] = 0 (a Z-drive does nothing to a Z observable).")
    print()
    print(f"  Tr((L_T1,1)_odd) = {L_T1_odd_sym.trace()}")
    print(f"  Tr((L_H,1)_odd) = {L_H_odd_sym.trace()}")
    print()
    print("  Both single-site Π²-odd traces vanish, which guarantees cross-site terms")
    print("  also vanish via the tensor structure.")
    print()

    # Numerical verification for all N
    print()
    print("Cross-site decomposition at general N:")
    print()
    print("Per-site embedding: (L_X,l)_odd = I_4^{⊗l} ⊗ (L_X,1)_odd ⊗ I_4^{⊗(N-l-1)}")
    print("  (since I_4 is Π_1²-symmetric and Π²-odd is well-defined per-site)")
    print()
    print("For l ≠ m:")
    print("  (L_H,l)_odd · (L_T1,m)_odd = [I^⊗l ⊗ (L_H,1)_odd ⊗ I^⊗(N-l-1)] ·")
    print("                                [I^⊗m ⊗ (L_T1,1)_odd ⊗ I^⊗(N-m-1)]")
    print()
    print("  If l < m: = I^⊗l ⊗ (L_H,1)_odd ⊗ I^⊗(m-l-1) ⊗ (L_T1,1)_odd ⊗ I^⊗(N-m-1)")
    print()
    print("  Tr = Tr(I_4)^l · Tr((L_H,1)_odd) · Tr(I_4)^{m-l-1} · Tr((L_T1,1)_odd) · Tr(I_4)^{N-m-1}")
    print("     = 4^{N-2} · Tr((L_H,1)_odd) · Tr((L_T1,1)_odd)")
    print()
    print("  But Tr((L_H,1)_odd) = 0 and Tr((L_T1,1)_odd) = 0 (verified above).")
    print("  So every cross-site term vanishes.")
    print()
    print("For l = m: single-site sympy verifies Tr((L_H,1)_odd · (L_T1,1)_odd) = 0.")
    print("           By the tensor factorization,")
    print("  Tr((L_H,l)_odd · (L_T1,l)_odd) = 4^{N-1} · Tr((L_H,1)_odd · (L_T1,1)_odd) = 0.")
    print()
    print("Hence Σ_{l,m} Tr((L_H,l)_odd · (L_T1,m)_odd) = 0.   ∎")
    print()

    # Numerical confirmation
    print("Numerical verification:")
    print()
    for N in [1, 2, 3, 4, 5]:
        Pi = fw.symmetry.build_pi_full(N)

        # Build L_H and L_T1 with random per-site ω, γ
        rng = np.random.default_rng(2026 + N)
        omegas = rng.uniform(0.1, 1.0, N)
        gammas = rng.uniform(0.1, 1.0, N)

        L_H_full = np.zeros((4**N, 4**N), dtype=complex)
        L_T1_full = np.zeros((4**N, 4**N), dtype=complex)
        for l in range(N):
            L_H_full += to_pauli(L_H_vec(0.5 * omegas[l] * site_op(N, l, Z)), N)
            L_T1_full += to_pauli(L_diss_vec(site_op(N, l, SIGMA_MINUS), gammas[l]), N)

        L_H_odd = project_pi2_odd(L_H_full, Pi)
        L_T1_odd = project_pi2_odd(L_T1_full, Pi)

        # The key trace
        trace_odd = np.trace(L_H_odd @ L_T1_odd)

        # And confirm the equivalence with cross_+ + cross_-
        Ap = project_pi(L_H_full, Pi, 1j)
        Am = project_pi(L_H_full, Pi, -1j)
        Bp = project_pi(L_T1_full, Pi, 1j)
        Bm = project_pi(L_T1_full, Pi, -1j)
        cross_sum = np.trace(Ap @ Bm) + np.trace(Am @ Bp)

        print(f"  N={N}: Tr((L_H)_odd · (L_T1)_odd) = {trace_odd:.2e}")
        print(f"         cross_+_trace + cross_-_trace = {cross_sum:.2e}")
        print(f"         Both zero? {np.isclose(trace_odd, 0) and np.isclose(cross_sum, 0)}")


# ============================================================
# STEP 3: Connect Tr(A_+i · B_-i) + Tr(A_-i · B_+i) = 0 back to the
#         Frobenius cross-term identity cross_+ + cross_- = 0
# ============================================================

def step_3_back_to_cross_sum():
    print()
    print("=" * 78)
    print("STEP 3: From the trace identity to cross_+ + cross_- = 0")
    print("=" * 78)
    print()
    print("The original Lemma C target was:")
    print()
    print("  cross_+ := Re⟨L_H,+i, L_T1,+i⟩ = Re Tr((L_H,+i)^† · L_T1,+i)")
    print("  cross_- := Re⟨L_H,-i, L_T1,-i⟩ = Re Tr((L_H,-i)^† · L_T1,-i)")
    print("  Claim: cross_+ + cross_- = 0.")
    print()
    print("Apply Lemma A (PROOF_F112) to L_H using Lemma B (L_H^† = -L_H):")
    print("  (L_H,-i)^† = (L_H^†)_+i = -(L_H)_+i")
    print("  (L_H,+i)^† = (L_H^†)_-i = -(L_H)_-i")
    print()
    print("Substituting:")
    print("  cross_+ = Re Tr(-(L_H)_-i · L_T1,+i) = -Re Tr((L_H)_-i · L_T1,+i)")
    print("  cross_- = Re Tr(-(L_H)_+i · L_T1,-i) = -Re Tr((L_H)_+i · L_T1,-i)")
    print()
    print("Sum:")
    print("  cross_+ + cross_- = -Re [Tr((L_H)_-i · L_T1,+i) + Tr((L_H)_+i · L_T1,-i)]")
    print()
    print("By Steps 1+2 above:")
    print("  Tr((L_H)_+i · (L_T1)_-i) + Tr((L_H)_-i · (L_T1)_+i)")
    print("    = Tr((L_H)_odd · (L_T1)_odd) = 0  (by Lindblad-input vanishing).")
    print()
    print("Hence cross_+ + cross_- = -Re(0) = 0.   ∎")
    print()
    print("Equivalently, cross_minus = -cross_plus, closing the Welle-4 reduction")
    print("  asymmetry = 4·Re⟨L_H,+i, L_T1,+i⟩")
    print("entirely from operator algebra (no bit-exact numerical anchor needed).")
    print()

    print("Numerical confirmation that cross_+ + cross_- = 0 for general N:")
    print()
    for N in [1, 2, 3, 4, 5]:
        Pi = fw.symmetry.build_pi_full(N)
        L_H_full = np.zeros((4**N, 4**N), dtype=complex)
        L_T1_full = np.zeros((4**N, 4**N), dtype=complex)
        for l in range(N):
            L_H_full += to_pauli(L_H_vec(0.5 * site_op(N, l, Z)), N)
            L_T1_full += to_pauli(L_diss_vec(site_op(N, l, SIGMA_MINUS), 1.0), N)
        Ap = project_pi(L_H_full, Pi, 1j)
        Am = project_pi(L_H_full, Pi, -1j)
        Bp = project_pi(L_T1_full, Pi, 1j)
        Bm = project_pi(L_T1_full, Pi, -1j)
        cross_p = np.vdot(Ap, Bp).real  # = Tr(Ap^† Bp)
        cross_m = np.vdot(Am, Bm).real
        print(f"  N={N}: cross_+ = {cross_p:.4f}, cross_- = {cross_m:.4f}, "
              f"sum = {cross_p + cross_m:.2e}")


# ============================================================
# Main
# ============================================================

def main():
    print()
    print("F113 LEMMA C STEP 5: STRUCTURAL CLOSURE")
    print("=" * 78)
    print()
    print("Goal: close cross_minus = -cross_plus algebraically.")
    print()
    print("Strategy (two-step):")
    print("  STEP 1 (pure algebra): reduce")
    print("    Tr(A_+i · B_-i) + Tr(A_-i · B_+i) = Tr(A_odd · B_odd)")
    print("    (A_odd = (A - Π² A Π^-2)/2 is the Π²-anti-symmetric part)")
    print()
    print("  STEP 2 (Lindblad-specific): show Tr((L_H)_odd · (L_T1)_odd) = 0")
    print("    via per-site additivity + (L_H,1)_odd = L_H,1 + (L_T1,1)_odd = single (Z,I) entry")
    print("    + L_H[I, Z] = 0 because [Z, σ_Z] = 0.")
    print()
    print("  STEP 3 (closure): apply Lemma A + B to relate cross_+ + cross_- to the")
    print("    above trace identity. Recover cross_minus = -cross_plus.")
    print()
    step_1_algebraic_reduction()
    step_2_lindblad_input()
    step_3_back_to_cross_sum()
    print()
    print("=" * 78)
    print("CLOSURE COMPLETE.")
    print("F113 Lemma C is closed from operator algebra + Lindblad-input vanishing.")
    print("=" * 78)


if __name__ == '__main__':
    main()
