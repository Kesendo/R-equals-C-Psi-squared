#!/usr/bin/env python3
"""Tier-1 derivation attempt: evaluate the Dyson sym3 = L_H²L'_dis + L_HL'_disL_H + L'_disL_H²
acting on ρ_0 directly, then extract c via partial trace.

If c_derived = 4/3 bit-exact, the Tier-1 closed form is:
  Δ_|00⟩(Q, K) = (4/3) · Q² · K³ + O(Q³K⁴)
for the dominant outcome of |0+0+⟩ pair (0,2) under Heisenberg ring + Z-dephasing.
"""
from __future__ import annotations

import sys
from pathlib import Path
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
    bonds = [(i, (i + 1) % N) for i in range(N)]
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (two_site_op(N, a, b, SX, SX)
                       + two_site_op(N, a, b, SY, SY)
                       + two_site_op(N, a, b, SZ, SZ))
    return H


def L_H(rho, H):
    return -1j * (H @ rho - rho @ H)


def L_dis(rho, N):
    """Z-dephasing dissipator at γ=1: Σ_l (Z_l ρ Z_l - ρ)."""
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


def main():
    N = 4
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    rho_0 = np.outer(psi_0, psi_0.conj())

    print("Tier-1 derivation attempt for Δ_|00⟩")
    print(f"  Setup: |0+0+⟩ N=4 Heisenberg ring, pair (0,2), |00⟩ outcome")
    print()

    # P_u(0) = ⟨00|_pair Tr_{1,3}[ρ_0] |00⟩_pair
    rho_pair_0 = reduced_density(rho_0, N, keep=[0, 2])
    P_u0 = float(np.real(rho_pair_0[0, 0]))
    print(f"  P_u(t=0) at |00⟩-pair-outcome = {P_u0}")
    print()

    # Build H at J=1 (so L_H carries one factor J=1)
    H = heisenberg_ring(N, J=1.0)

    # Apply each ordering of {L_H, L_H, L'_dis} to ρ_0
    # Ordering 1: L_H² · L'_dis
    A1 = L_dis(rho_0, N)
    A1 = L_H(A1, H)
    A1 = L_H(A1, H)

    # Ordering 2: L_H · L'_dis · L_H
    A2 = L_H(rho_0, H)
    A2 = L_dis(A2, N)
    A2 = L_H(A2, H)

    # Ordering 3: L'_dis · L_H²
    A3 = L_H(rho_0, H)
    A3 = L_H(A3, H)
    A3 = L_dis(A3, N)

    sym3 = A1 + A2 + A3

    # Extract ⟨00|_pair Tr_{1,3}[sym3] |00⟩_pair
    sym3_pair = reduced_density(sym3, N, keep=[0, 2])
    sym3_element = float(np.real(sym3_pair[0, 0]))

    print(f"  sym3 = L_H²L'_dis + L_HL'_disL_H + L'_disL_H²")
    print(f"  ⟨00|_pair Tr_{{1,3}}[sym3 ρ_0] |00⟩_pair = {sym3_element}")
    print()

    # Δ contribution at γ¹J² order:
    #   ΔP_|00⟩ ≈ (γt³/6) · sym3_element  with J=γ=1 baked in
    #          = (J² γ t³ / 6) · sym3_element / (J²γ) at general (J, γ)
    #          = (1/6) · sym3_element · J²γt³
    # And Δ = ΔP / P_u0:
    #   Δ = (sym3_element / (6 · P_u0)) · J²γt³ = (sym3_element / (6·P_u0)) · Q²K³

    c_derived = sym3_element / (6.0 * P_u0)
    print(f"  Predicted coefficient c = sym3_element / (6 · P_u0)")
    print(f"                          = {sym3_element} / ({6 * P_u0})")
    print(f"                          = {c_derived}")
    print()

    # Compare with 4/3
    four_thirds = 4.0 / 3.0
    print(f"  4/3 = {four_thirds}")
    print(f"  c_derived - 4/3 = {c_derived - four_thirds:+.10e}")
    print()

    # Try to identify a clean rational
    print(f"  As fraction (approximation):")
    frac = Fraction(c_derived).limit_denominator(100)
    print(f"    c ≈ {frac} = {float(frac)}")
    print()

    # Cross-check numerically
    print(f"  Numerical extraction (from previous run, mean over 16 samples): c ≈ 1.32992")
    print(f"  Symbolic prediction:                                            c = {c_derived:.5f}")
    print(f"  Match within numerical precision: {abs(c_derived - 1.32992) < 0.01}")


if __name__ == "__main__":
    main()
