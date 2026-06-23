#!/usr/bin/env python3
"""Attempt Tier-1 coefficient for Δ_|00⟩ ≈ c · Q² · K³ at dominant outcome.

Strategy:
  1. Sample Δ_|00⟩ on a fine grid of very small (Q, K) — deep in perturbative regime.
  2. Check Δ_|00⟩ / (Q² · K³) converges to a constant c.
  3. Cross-check by computing the Dyson-series leading term symbolically:
       ΔP_|00⟩ ≈ γ · t³/6 · ⟨00_pair| Tr_{1,3}[ sym{L_H, L_H, L'_dis} ρ_0 ] |00⟩_pair
     where the symmetric ordering sym{L_H, L_H, L'_dis} = L_H²L'_dis + L_HL'_disL_H +
     L'_disL_H². This is the leading non-vanishing γ-contribution in the time-Taylor
     expansion of e^{Lt} ρ_0.

  4. Compare numerical c with the symbolic derivation. If they match, that IS the
     Tier-1 coefficient (provided the Dyson term is the leading one).

This is the Februar-style 'try to derive, see if it lands' move.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
EYE = np.eye(2, dtype=complex)


def heisenberg_ring(N, J=1.0):
    bonds = [(i, (i + 1) % N) for i in range(N)]

    def op_at(op_a, op_b, a, b):
        out = np.array([[1]], dtype=complex)
        for k in range(N):
            if k == a:
                out = np.kron(out, op_a)
            elif k == b:
                out = np.kron(out, op_b)
            else:
                out = np.kron(out, EYE)
        return out

    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (op_at(SX, SX, a, b) + op_at(SY, SY, a, b) + op_at(SZ, SZ, a, b))
    return H


def reduced_density(rho, N, keep):
    n_keep = len(keep)
    trace = [i for i in range(N) if i not in keep]
    rho_tensor = rho.reshape([2] * N + [2] * N)
    for q in sorted(trace, reverse=True):
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + (N - sum(1 for t in trace if t > q)))
    return rho_tensor.reshape((2 ** n_keep, 2 ** n_keep))


def P_00_pair(J, gamma, t, N=4, use_lindblad=True):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    rho_0 = np.outer(psi_0, psi_0.conj())

    H = heisenberg_ring(N, J)
    if use_lindblad and gamma > 0:
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        rho_vec_t = expm(L * t) @ rho_0.flatten()
        rho_t = rho_vec_t.reshape(2 ** N, 2 ** N)
    else:
        U = expm(-1j * H * t)
        rho_t = U @ rho_0 @ U.conj().T

    rho_pair = reduced_density(rho_t, N, keep=[0, 2])
    return float(np.real(rho_pair[0, 0]))


def delta_00(J, gamma, t, N=4):
    P_u = P_00_pair(J, gamma, t, N=N, use_lindblad=False)
    P_l = P_00_pair(J, gamma, t, N=N, use_lindblad=True)
    if P_u < 1e-14:
        return np.nan
    return P_l / P_u - 1


def main():
    print("Δ_|00⟩ coefficient extraction — small-(Q, K) limit")
    print("  Hypothesis: Δ_|00⟩ ≈ c · Q² · K³ at small (Q, K)")
    print()

    # Deep perturbative regime: small J · t and small γ · t
    # Sample multiple (J, γ, t) with various Q²·K³ values
    samples = []
    print(f"{'γ':>6s}  {'J':>6s}  {'t':>6s}    {'Q':>6s}  {'K':>8s}  {'Q²K³':>10s}    "
          f"{'Δ_|00⟩':>12s}    {'c=Δ/(Q²K³)':>13s}")
    print('-' * 100)
    test_configs = [
        # γ,    J,     t
        (0.05, 1.0, 0.05),
        (0.05, 1.0, 0.10),
        (0.05, 1.0, 0.20),
        (0.05, 2.0, 0.05),
        (0.05, 2.0, 0.10),
        (0.05, 0.5, 0.05),
        (0.05, 0.5, 0.10),
        (0.05, 0.5, 0.20),
        (0.05, 4.0, 0.05),
        (0.05, 4.0, 0.025),
        (0.025, 1.0, 0.05),
        (0.025, 1.0, 0.10),
        (0.025, 2.0, 0.05),
        (0.1, 1.0, 0.05),
        (0.1, 1.0, 0.10),
        (0.1, 0.5, 0.05),
    ]
    for gamma, J, t in test_configs:
        Q = J / gamma
        K = gamma * t
        Q2K3 = Q ** 2 * K ** 3
        d = delta_00(J, gamma, t)
        if d is not None and not np.isnan(d) and Q2K3 > 1e-12:
            c = d / Q2K3
            samples.append((Q, K, Q2K3, d, c))
            print(f"{gamma:>6.3f}  {J:>6.3f}  {t:>6.3f}    "
                  f"{Q:>6.1f}  {K:>8.5f}  {Q2K3:>10.5g}    "
                  f"{d:>12.5g}    {c:>13.5g}")

    # Filter for very small Q²K³ (deepest perturbative regime)
    small = [s for s in samples if s[2] < 0.01]
    if small:
        c_values = [s[4] for s in small]
        c_mean = float(np.mean(c_values))
        c_std = float(np.std(c_values))
        print()
        print(f"Coefficient c at small Q²K³ (< 0.01):")
        print(f"  mean  = {c_mean:.5f}")
        print(f"  std   = {c_std:.5f}")
        print(f"  range = [{min(c_values):.5f}, {max(c_values):.5f}]")
        print()
        print(f"Empirical (Tier-2):  Δ_|00⟩ ≈ {c_mean:.4f} · Q² · K³  in the deep")
        print(f"                                                       perturbative regime.")
    print()

    # Sketch of the derivation pathway:
    print("Derivation pathway (Tier-1 attempt):")
    print()
    print("  The Liouvillian is L = L_H + γ · L'_dis where")
    print("    L_H[ρ]    = -i[H, ρ]                    (∝ J)")
    print("    L'_dis[ρ] = Σ_l (Z_l ρ Z_l - ρ)         (γ-free operator)")
    print()
    print("  ρ(t) = e^{Lt} ρ_0 = ρ_0 + Lt·ρ_0 + L²t²/2·ρ_0 + L³t³/6·ρ_0 + ...")
    print()
    print("  Expanding L³ at γ-coefficient (1 γ-vertex, 2 H-vertices):")
    print()
    print("    L³|_{γ-coef} = L_H² L'_dis  +  L_H L'_dis L_H  +  L'_dis L_H²")
    print()
    print("  Each L_H carries one factor J; L'_dis carries no J. So the contribution")
    print("  to ρ(t) at γ¹ J² is t³/6 · γ · (sum of 3 orderings) acting on ρ_0.")
    print()
    print("  The dimensionful prefactor: γ · J² · t³ = (J/γ)² · (γt)³ = Q² · K³.")
    print("  ✓ matches the empirical scaling.")
    print()
    print("  The coefficient c is the structural constant:")
    print()
    print("    c = (1/6) · ⟨00|_pair Tr_{1,3}[ sym{L_H, L_H, L'_dis}/J²/γ · ρ_0 ] |00⟩_pair / P_u(t→0)")
    print()
    print("  Need to evaluate this matrix element with H = Heisenberg ring (N=4),")
    print("  Z-dephasing on each site, |0+0+⟩ initial state. Tractable in principle —")
    print("  this script verifies c numerically; the explicit symbolic computation")
    print("  (yielding e.g. c = 3/2 or similar rational/integer-combinatoric form) is")
    print("  the next step toward Tier-1.")


if __name__ == "__main__":
    main()
