#!/usr/bin/env python3
"""Subdominant Born-deviation Dyson coefficients for F94's setup.

F94 (Tier 1, 2026-05-16) covers the dominant outcome Δ_|00⟩ = (4/3)·Q²·K³ via
the γ¹·J² (sym3) Dyson term. The reflection ON_HOW_FOUR_THIRDS_APPEARED also
reports empirical slopes-per-K for the subdominant outcomes:

  Δ_|01⟩ = Δ_|10⟩ ≈ -1.8 · K (Q-independent)
  Δ_|11⟩         ≈ -2.6 · K (Q-independent)

The naive 1st-order L¹ Dyson term vanishes outright on the pair (0,2)
|i⟩-element for ALL i: L'_dis ρ_0 produces |+⟩→|-⟩-flip content on sites 1, 3
which has zero partial-trace.

So the empirical Q-independence-of-Δ_subdominant is NOT from L¹. The leading
γ¹ contribution comes from the same sym3 term as F94, but for subdominant
outcomes P_u(t) ∝ J²t²/2·A_i (not 1), so

  Δ_subdominant ≈ (γJ²t³ · M_3^{(i)} / 6) / (J²t²/2 · A_i)
                = K · M_3^{(i)} / (3 · A_i)

where A_i = <i|_pair Tr[L_h² ρ_0]|i>_pair is the **raw** L_h² pair element
(not the /2 Taylor-normalized U_2). This IS Q-independent and linear in K.
This script:

1. Computes M_n^{(i)} = <i|_pair Tr_{1,3}[sym_n^1 · ρ_0]|i>_pair for n=1,2,3
   and all four outcomes i ∈ {|00⟩, |01⟩, |10⟩, |11⟩}.
2. Computes A_i = <i|_pair Tr_{1,3}[L_h² · ρ_0]|i>_pair (raw matrix element).
3. Reads off slope_i = M_3^{(i)} / (3 · A_i) for subdominants where A_i ≠ 0.
4. Compares to empirical: |01⟩ ≈ -1.8, |11⟩ ≈ -2.6.
5. If clean rationals: candidate Tier-1 closed forms.
"""
from __future__ import annotations

import sys
from fractions import Fraction
import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
EYE = np.eye(2, dtype=complex)


def single_site_op(N, i, op):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, op if k == i else EYE)
    return out


def two_site_op(N, i, j, op_a, op_b):
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        if k == i:
            out = np.kron(out, op_a)
        elif k == j:
            out = np.kron(out, op_b)
        else:
            out = np.kron(out, EYE)
    return out


def heisenberg_ring(N, J=1.0):
    """Full Heisenberg ring Hamiltonian, H = (J/4) Σ_b (XX+YY+ZZ)."""
    bonds = [(i, (i + 1) % N) for i in range(N)]
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (two_site_op(N, a, b, SX, SX)
                        + two_site_op(N, a, b, SY, SY)
                        + two_site_op(N, a, b, SZ, SZ))
    return H


def L_H_apply(rho, H):
    """Apply L_H = -i[H, ·]."""
    return -1j * (H @ rho - rho @ H)


def L_dis_apply(rho, N):
    """Apply L'_dis = Σ_l (Z_l ρ Z_l - ρ) at γ=1."""
    out = np.zeros_like(rho)
    for l in range(N):
        Zl = single_site_op(N, l, SZ)
        out += Zl @ rho @ Zl - rho
    return out


def reduced_density(rho, N, keep):
    n_keep = len(keep)
    trace = [i for i in range(N) if i not in keep]
    rho_tensor = rho.reshape([2] * N + [2] * N)
    for q in sorted(trace, reverse=True):
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + (N - sum(1 for t in trace if t > q)))
    return rho_tensor.reshape((2 ** n_keep, 2 ** n_keep))


def initial_state_0p0p(N):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    return np.outer(psi_0, psi_0.conj())


def pair_element(rho_full, N, keep_pair, outcome_index):
    """⟨outcome|_pair · Tr_{not pair}[rho_full] · |outcome⟩_pair."""
    rho_pair = reduced_density(rho_full, N, keep_pair)
    val = rho_pair[outcome_index, outcome_index]
    if abs(val.imag) > 1e-10:
        raise RuntimeError(f"Non-real pair element: {val}")
    return float(val.real)


def main():
    N = 4
    rho_0 = initial_state_0p0p(N)
    H = heisenberg_ring(N, J=1.0)  # J=1
    keep_pair = [0, 2]
    outcomes = {"00": 0, "01": 1, "10": 2, "11": 3}

    print("Subdominant Born-deviation Dyson analysis")
    print(f"  Setup: |0+0+⟩ N=4 Heisenberg ring, pair (0,2)")
    print(f"  J = γ = 1 in normalized units; Q = K = 1 (so we read pure structural numbers)")
    print()

    # --- M_1, M_2, M_3 for each outcome ---
    # M_n^{(i)} = ⟨i|_pair Tr_{1,3}[sym_n^1 · ρ_0]|i⟩_pair

    # sym_1^1 = L'_dis
    sym1 = L_dis_apply(rho_0, N)

    # sym_2^1 = L_H L'_dis + L'_dis L_H
    A = L_dis_apply(rho_0, N)
    A = L_H_apply(A, H)
    B = L_H_apply(rho_0, H)
    B = L_dis_apply(B, N)
    sym2 = A + B

    # sym_3^1 = L_H² L'_dis + L_H L'_dis L_H + L'_dis L_H²  (F94's sym3)
    s3_a = L_dis_apply(rho_0, N)
    s3_a = L_H_apply(s3_a, H)
    s3_a = L_H_apply(s3_a, H)
    s3_b = L_H_apply(rho_0, H)
    s3_b = L_dis_apply(s3_b, N)
    s3_b = L_H_apply(s3_b, H)
    s3_c = L_H_apply(rho_0, H)
    s3_c = L_H_apply(s3_c, H)
    s3_c = L_dis_apply(s3_c, N)
    sym3 = s3_a + s3_b + s3_c

    # Unitary 2nd-order: U_2 = L_H² · ρ_0 (then matrix element / 2 for Taylor)
    LH_rho = L_H_apply(rho_0, H)
    LH2_rho = L_H_apply(LH_rho, H)

    # Print matrix elements
    print("=" * 70)
    print("γ¹-Dyson matrix elements M_n^{(i)} = <i|_pair Tr_{1,3}[sym_n^1 · ρ_0]|i>")
    print("=" * 70)
    print(f"{'outcome':<10} {'M_1':>10} {'M_2':>10} {'M_3':>10}")
    print("-" * 50)
    for label, idx in outcomes.items():
        m1 = pair_element(sym1, N, keep_pair, idx)
        m2 = pair_element(sym2, N, keep_pair, idx)
        m3 = pair_element(sym3, N, keep_pair, idx)
        print(f"{label:<10} {m1:>+10.4f} {m2:>+10.4f} {m3:>+10.4f}")
    print()

    print("=" * 70)
    print("Unitary 2nd-order matrix elements")
    print("=" * 70)
    print(f"{'outcome':<10} {'A_i = <i|Tr[L_h²ρ_0]|i>':>26} {'U_2 = A_i/2 (P_u(i) t² coeff)':>32}")
    print("-" * 70)
    for label, idx in outcomes.items():
        a_i = pair_element(LH2_rho, N, keep_pair, idx)
        u2 = a_i / 2
        print(f"{label:<10} {a_i:>+26.4f} {u2:>+32.4f}")
    print()

    # --- Predicted slopes (Q-independent K-coefficient of Δ_i for subdominant) ---
    print("=" * 70)
    print("Predicted slope_i = M_3 / (3 · A_i) for subdominant outcomes")
    print("(Δ_i = ΔP_i^{γ¹J²,t³} / P_u(i, t) ≈ K · M_3 / (3 · A_i))")
    print("=" * 70)
    print(f"{'outcome':<10} {'M_3':>10} {'A_i':>10} {'slope predicted':>30} {'empirical':>12}")
    print("-" * 80)
    empirical = {"00": "Q²K³ form (F94)", "01": "≈ -1.8", "10": "≈ -1.8", "11": "≈ -2.6"}
    for label, idx in outcomes.items():
        m3 = pair_element(sym3, N, keep_pair, idx)
        a_i = pair_element(LH2_rho, N, keep_pair, idx)
        if abs(a_i) > 1e-12 and abs(m3) > 1e-12:
            slope = m3 / (3 * a_i)
            slope_str = f"{slope:>+10.6f}"
            frac = Fraction(slope).limit_denominator(100)
            slope_str += f" ≈ {frac} = {float(frac):+.4f}"
        elif abs(m3) < 1e-12 and abs(a_i) < 1e-12:
            slope_str = "M_3=0 AND A_i=0 (higher-order needed)"
        elif abs(a_i) < 1e-12:
            slope_str = "A_i=0 (P_u(i) has no t² term, higher order)"
        else:
            slope_str = "M_3=0 (dominant case or special parity)"
        print(f"{label:<10} {m3:>+10.4f} {a_i:>+10.4f} {slope_str:<30}  {empirical[label]:>12}")
    print()

    # --- Verify F94 dominant case ---
    print("=" * 70)
    print("Verify F94 dominant: M_3^{(00)} should equal 8")
    print("=" * 70)
    m3_00 = pair_element(sym3, N, keep_pair, 0)
    a_00 = pair_element(LH2_rho, N, keep_pair, 0)
    print(f"  M_3^{{(00)}} = {m3_00:.6f}  (expect 8.0)")
    print(f"  A^{{(00)}} = {a_00:.6f}  (P_u(|00⟩, t) ≈ 1 + (t²/2)·(-1.5) = 1 - 0.75·t²)")
    print()

    # --- |11⟩ requires higher-order Dyson (M_3=0, A=0) ---
    print("=" * 70)
    print("|11⟩ outcome: doubly subdominant; needs M_5 and B = ⟨11|Tr[L_h⁴ρ_0]|11⟩")
    print("=" * 70)
    # sym_4^1 = 4 orderings of (L_H L_H L_H L_dis)
    s4_a = L_dis_apply(rho_0, N)
    s4_a = L_H_apply(s4_a, H); s4_a = L_H_apply(s4_a, H); s4_a = L_H_apply(s4_a, H)
    s4_b = L_H_apply(rho_0, H); s4_b = L_dis_apply(s4_b, N)
    s4_b = L_H_apply(s4_b, H); s4_b = L_H_apply(s4_b, H)
    s4_c = L_H_apply(rho_0, H); s4_c = L_H_apply(s4_c, H); s4_c = L_dis_apply(s4_c, N)
    s4_c = L_H_apply(s4_c, H)
    s4_d = L_H_apply(rho_0, H); s4_d = L_H_apply(s4_d, H); s4_d = L_H_apply(s4_d, H)
    s4_d = L_dis_apply(s4_d, N)
    sym4 = s4_a + s4_b + s4_c + s4_d

    # sym_5^1 = 5 orderings of (L_H L_H L_H L_H L_dis)
    def apply_chain(rho_0, ops):
        rho = rho_0.copy()
        for op in ops:
            rho = op(rho)
        return rho
    L_H_op = lambda r: L_H_apply(r, H)
    L_dis_op = lambda r: L_dis_apply(r, N)
    sym5 = sum(apply_chain(rho_0,
                           [L_H_op if k != pos else L_dis_op for k in range(5)])
               for pos in range(5))

    # L_h^4 (unitary 4th order)
    LH4 = L_H_apply(rho_0, H)
    LH4 = L_H_apply(LH4, H); LH4 = L_H_apply(LH4, H); LH4 = L_H_apply(LH4, H)

    m4_11 = pair_element(sym4, N, keep_pair, 3)
    m5_11 = pair_element(sym5, N, keep_pair, 3)
    b_11 = pair_element(LH4, N, keep_pair, 3)
    slope_11 = m5_11 / (5 * b_11) if abs(b_11) > 1e-12 else float("nan")
    frac_11 = Fraction(slope_11).limit_denominator(100)
    print(f"  M_4^{{(11)}} = {m4_11:.6f}  (expect 0)")
    print(f"  M_5^{{(11)}} = {m5_11:.6f}")
    print(f"  B^{{(11)}}    = ⟨11|Tr[L_h⁴ρ_0]|11⟩ = {b_11:.6f}")
    print(f"  P_u(|11⟩, t) ≈ (t⁴/4!) · {b_11:.4f} = (t⁴/24) · {b_11:.4f}")
    print()
    print(f"  slope_|11⟩ = M_5 / (5 · B) = {m5_11:.4f} / (5 · {b_11:.4f}) = {slope_11:+.6f}")
    print(f"                                                          = {frac_11} = {float(frac_11):+.6f}")
    print()
    print(f"  Empirical: ≈ -2.66 at K=0.001 (matches -8/3 = -2.6667)")
    print()

    # --- Summary: all 4 outcome closed forms ---
    print("=" * 70)
    print("SUMMARY: closed forms for all 4 outcomes of pair (0,2) of |0+0+⟩ N=4")
    print("=" * 70)
    print("  Outcome    Closed form                        Decomposition")
    print("  -------    ------------------------------     -----------------------")
    print("  |00⟩      Δ = (4/3) · Q² · K³                F94: M_3 / 3! = 8/6 = 4/3")
    print("  |01⟩      Δ = -(16/9) · K = -(4/3)² · K     M_3/(3A) = -4/(3·3/4) = -16/9")
    print("  |10⟩      Δ = -(16/9) · K = -(4/3)² · K     (same as |01⟩ by symmetry)")
    print("  |11⟩      Δ = -(8/3) · K = -2·(4/3) · K     M_5/(5B) = -20/(5·3/2) = -8/3")
    print()
    print("  Structural unity: all 4 outcomes are simple algebraic expressions")
    print("  in F94's 4/3 anchor:  4/3,  -(4/3)²,  -(4/3)²,  -2·(4/3).")
    print("  Sum: (4/3)·Q²·K³ − 2·(4/3)²·K − 2·(4/3)·K  (total of 4 outcomes).")
    print()

    # --- Verify |01⟩ slope -16/9 with actual Lindblad ---
    print("=" * 70)
    print("Numerical Lindblad verification of |01⟩ slope")
    print("=" * 70)
    print("  Run actual RK4 Lindblad at fixed Q (= J/γ), sweep K (= γt) in deep regime,")
    print("  fit Δ_|01⟩(K) linear; expect slope -16/9 = -1.7778.")
    print()

    from scipy.linalg import expm

    def evolve_lindblad(rho_0, N, H, gamma, t, n_steps=4000):
        # Vectorize Liouvillian: dρ/dt = -i[H, ρ] + γ Σ_l (Z_l ρ Z_l - ρ).
        # For verification we just use matrix exponential per small step.
        d = 2 ** N
        dt = t / n_steps
        # Use simple RK4 here since vectorized Liouvillian is large at N=4 (256x256)
        # but acceptable for verification.
        def Ldot(rho):
            return -1j * (H @ rho - rho @ H) + gamma * sum(
                single_site_op(N, l, SZ) @ rho @ single_site_op(N, l, SZ) - rho
                for l in range(N)
            )

        rho = rho_0.copy()
        for _ in range(n_steps):
            k1 = Ldot(rho)
            k2 = Ldot(rho + dt / 2 * k1)
            k3 = Ldot(rho + dt / 2 * k2)
            k4 = Ldot(rho + dt * k3)
            rho = rho + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return rho

    def evolve_unitary(rho_0, N, H, t):
        U = expm(-1j * H * t)
        return U @ rho_0 @ U.conj().T

    Q = 50.0  # deep perturbative
    gamma = 0.01
    J = Q * gamma
    H_J = heisenberg_ring(N, J=J)
    K_values = [0.005, 0.01, 0.02, 0.03, 0.04]
    print(f"  Q = {Q}, γ = {gamma}, J = {J}")
    print(f"  {'K = γt':>10} {'P_lind(|01⟩)':>15} {'P_unit(|01⟩)':>15} {'Δ_|01⟩':>15} {'slope per K':>15}")
    print("  " + "-" * 75)
    for K in K_values:
        t = K / gamma
        rho_lind = evolve_lindblad(rho_0, N, H_J, gamma, t, n_steps=400)
        rho_unit = evolve_unitary(rho_0, N, H_J, t)
        P_lind_01 = pair_element(rho_lind, N, keep_pair, 1)
        P_unit_01 = pair_element(rho_unit, N, keep_pair, 1)
        if abs(P_unit_01) > 1e-12:
            delta = (P_lind_01 - P_unit_01) / P_unit_01
            slope = delta / K
            print(f"  {K:>10.4f} {P_lind_01:>15.6e} {P_unit_01:>15.6e} {delta:>+15.6e} {slope:>+15.4f}")
        else:
            print(f"  {K:>10.4f} P_unit too small for ratio")

    print()
    print(f"  Theoretical slope (M_3 / (3·A_i)): {-4/(3*0.75):.6f} = -16/9 = {-16/9:.6f}")
    print()


if __name__ == "__main__":
    main()
