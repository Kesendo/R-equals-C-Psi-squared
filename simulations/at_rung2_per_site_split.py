#!/usr/bin/env python3
"""The rung-2 four splits per-site: |00⟩⟨11| is an exact N=2 eigenmode at −2(γ₁+γ₂).

The Absorption Theorem's rung 2 (rate 4γ, two absorption quanta of 2γ each) is the
single object behind every dynamical 4γ₀ in the repository: t_peak = 1/(4γ₀)
(TPeakLaw), the clock's Takt pin 4γ₀, the L_eff mirror axis −4γ₀, the
approach-family carrier e^{−4γt}, and the F25 Bell+ rate. This script certifies the
per-site split that makes "two quanta" literal: for the N=2 XY chain
H = J(XX + YY) under site-dependent Z-dephasing (γ₁, γ₂), the coherence |00⟩⟨11|
is an EXACT Liouvillian eigenmode,

    L(|00⟩⟨11|) = −2(γ₁ + γ₂)·|00⟩⟨11|,    residual 0.0,

because (XX+YY)|00⟩ = (XX+YY)|11⟩ = 0 (the Hamiltonian part vanishes on this mode,
for any J) and each site contributes its own absorption quantum 2γ_l (Theorem 2,
the vector form: the bra and ket labels differ at both sites, Δ₁ = Δ₂ = 1). The
uniform case γ₁ = γ₂ = γ recovers the rung-2 rate 4γ.

Counterfactual arithmetic (the "two different fours" disambiguation, see
PROOF_ABSORPTION_THEOREM.md, Remark (two different fours)): the rung-2 four is
CENTRE-based, the midpoint of the HD=(1,3) channel pair, (2γ+6γ)/2 = 4γ. The
discriminant four (a₋₁ = d² on the Pi2 dyadic ladder, the EP discriminant 4γ₀²)
is HALF-GAP-SQUARED-based, ((6γ−2γ)/2)² = (2γ)² = 4γ², a quantum squared. They
coincide only at HD=(1,3), glued by the dyadic root d² − 2d = 0. At the
counterfactual HD=(1,5) pair the two genealogies separate: centre-based
(2γ+10γ)/2 = 6γ vs discriminant-based ((10γ−2γ)/2)² = 16γ². 6 ≠ 16: same number
"4" at home, two different objects abroad.

Self-validating: asserts residual < 1e-12 for two non-uniform γ profiles (each at
two J values) plus the uniform 4γ case, and the counterfactual 6-vs-16 arithmetic.
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def liouvillian_n2_xy(J: float, gamma1: float, gamma2: float) -> np.ndarray:
    """N=2 XY-chain Liouvillian, column-stacking convention vec(AρB) = (Bᵀ⊗A)vec(ρ)."""
    H = J * (np.kron(X, X) + np.kron(Y, Y))
    Id = np.eye(4, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for gamma, Zk in ((gamma1, np.kron(Z, I2)), (gamma2, np.kron(I2, Z))):
        # D(ρ) = ZρZ − ρ; Z real symmetric, so vec(ZρZ) = (Z⊗Z)vec(ρ).
        L += gamma * (np.kron(Zk, Zk) - np.eye(16, dtype=complex))
    return L


def mode_00_11() -> np.ndarray:
    """vec(|00⟩⟨11|), column-stacking: ρ[0, 3] = 1 → index 0 + 4·3 = 12."""
    v = np.zeros(16, dtype=complex)
    v[12] = 1.0
    return v


def check(J: float, gamma1: float, gamma2: float) -> float:
    L = liouvillian_n2_xy(J, gamma1, gamma2)
    v = mode_00_11()
    lam = -2.0 * (gamma1 + gamma2)
    residual = float(np.max(np.abs(L @ v - lam * v)))
    print(f"  J={J:<4}  γ=({gamma1}, {gamma2})   λ_predicted = {lam:+.6f}"
          f"   residual = {residual:.3e}")
    assert residual < 1e-12, f"|00⟩⟨11| not exact at γ=({gamma1},{gamma2}): {residual}"
    return residual


def main() -> None:
    print("=" * 78)
    print("  rung-2 per-site split: |00⟩⟨11| exact at −2(γ₁+γ₂) on the N=2 XY chain")
    print("=" * 78)

    # Non-uniform profiles: the 4 splits into two per-site absorption quanta.
    for J in (1.0, 0.7):
        for g1, g2 in ((0.05, 0.13), (0.02, 0.31)):
            check(J, g1, g2)

    # Uniform case: rate = 4γ, the rung-2 four itself.
    gamma = 0.05
    L = liouvillian_n2_xy(1.0, gamma, gamma)
    v = mode_00_11()
    residual = float(np.max(np.abs(L @ v - (-4.0 * gamma) * v)))
    print(f"  J=1.0   γ=({gamma}, {gamma})   λ_predicted = {-4.0 * gamma:+.6f}"
          f"   residual = {residual:.3e}   (uniform: rate = 4γ)")
    assert residual < 1e-12, f"uniform 4γ case not exact: {residual}"

    # Counterfactual arithmetic: the two fours coincide at HD=(1,3), separate at (1,5).
    centre_13, halfgap_sq_13 = (2 + 6) / 2, ((6 - 2) / 2) ** 2
    centre_15, halfgap_sq_15 = (2 + 10) / 2, ((10 - 2) / 2) ** 2
    print(f"  HD=(1,3): centre = {centre_13:g}γ, half-gap² = {halfgap_sq_13:g}γ²  (both '4')")
    print(f"  HD=(1,5): centre = {centre_15:g}γ, half-gap² = {halfgap_sq_15:g}γ²  (6 ≠ 16)")
    assert centre_13 == halfgap_sq_13 == 4.0
    assert (centre_15, halfgap_sq_15) == (6.0, 16.0)

    print("-" * 78)
    print("  PASS: the dynamical 4γ₀ is rung 2 (two per-site quanta 2γ_l), exact;")
    print("        the discriminant 4 (a₋₁ = d², half-gap squared) is a different four,")
    print("        glued to it only at HD=(1,3) by d² − 2d = 0.")


if __name__ == "__main__":
    main()
