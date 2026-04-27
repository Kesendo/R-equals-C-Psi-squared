#!/usr/bin/env python3
"""EQ-020 hardware prediction from the chiral mirror law.

Claim (computationally verifiable, hardware-testable):

For an OBC uniform XY chain at N qubits with Z-dephasing, prepared in
the bonding state |φ_k⟩ = (|vac⟩ + |ψ_k⟩)/√2 where |ψ_k⟩ is the k-th
single-excitation sine mode, the per-site Bloch components evolve such
that the K_1-paired state |φ_{N+1-k}⟩ gives:

    ⟨X_i⟩(t)|_{ψ_{N+1-k}} = +(−1)^i · ⟨X_i⟩(t)|_{ψ_k}
    ⟨Y_i⟩(t)|_{ψ_{N+1-k}} = −(−1)^i · ⟨Y_i⟩(t)|_{ψ_k}    (extra sign!)
    ⟨Z_i⟩(t)|_{ψ_{N+1-k}} = +⟨Z_i⟩(t)|_{ψ_k}     (identical)

and consequently:

    P_i(t)|_{ψ_{N+1-k}} = P_i(t)|_{ψ_k}    (per-site purity identical)

Derivation. ψ_k is real for OBC sine modes, so off-diagonal ρ_i[0,1](t) =
e^{iE_k t} ψ_k(i) / 2 (in single-excitation bonding-mode dynamics). For
K_1-paired states, ψ_{N+1-k}(i) = (−1)^i ψ_k(i) AND E_{N+1-k} = −E_k.
The time-evolved ρ_i[0,1] gets:
  - sign (-1)^i from the wave-function mirror
  - complex conjugation from the energy reversal
⟨X_i⟩ = 2 Re(ρ_i[0,1]) is invariant under conjugation → only (−1)^i.
⟨Y_i⟩ = -2 Im(ρ_i[0,1]) flips under conjugation → −(−1)^i.

The Y EXTRA-SIGN is a SHARP test, distinguishing the chiral mirror from
naive sublattice-Z₂ symmetry alone.

Hardware test: prepare both states, do single-qubit tomography at several
times, verify the predicted sign relations. Backend-independent prediction:
follows from K_1-symmetry of H + Z-dephasing. Any deviation indicates
state-prep errors, K_1-symmetry-breaking noise (coherent, non-sublattice-
symmetric), or measurement-axis miscalibration.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sps

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from _eq014_chiral_mirror_multi_N import (
    build_xy_chain_H, build_L, sine_mode_state, rk4_step,
    site_op_sparse, X, Y,
)


def bonding_density(N, k, d):
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = sine_mode_state(N, k, d)
    phi = (vac + psi) / np.sqrt(2)
    rho = np.outer(phi, phi.conj())
    return rho.flatten()


def site_z_op(N, i):
    Z2 = sps.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
    full = sps.eye(1, dtype=complex, format='csr')
    I2 = sps.eye(2, dtype=complex, format='csr')
    for j in range(N):
        full = sps.kron(full, Z2 if j == i else I2, format='csr')
    return full


def expectation_per_site(rho_flat, N):
    """Returns (⟨X_i⟩, ⟨Y_i⟩, ⟨Z_i⟩, P_i) per site, len-N each."""
    d = 2 ** N
    rho = rho_flat.reshape(d, d)
    Xexp = np.zeros(N)
    Yexp = np.zeros(N)
    Zexp = np.zeros(N)
    Pe = np.zeros(N)
    for i in range(N):
        Xi = site_op_sparse(X, i, N).toarray()
        Yi = site_op_sparse(Y, i, N).toarray()
        Zi = site_z_op(N, i).toarray()
        Xexp[i] = float(np.real(np.trace(rho @ Xi)))
        Yexp[i] = float(np.real(np.trace(rho @ Yi)))
        Zexp[i] = float(np.real(np.trace(rho @ Zi)))
        # P_i = Tr(ρ_i²) = (1 + ⟨X⟩² + ⟨Y⟩² + ⟨Z⟩²) / 2
        Pe[i] = 0.5 * (1.0 + Xexp[i] ** 2 + Yexp[i] ** 2 + Zexp[i] ** 2)
    return Xexp, Yexp, Zexp, Pe


def evolve_to_t(L, rho_flat, t, dt_small=0.02):
    n_steps = max(1, int(np.ceil(t / dt_small)))
    dt = t / n_steps
    current = rho_flat.astype(complex).copy()
    for _ in range(n_steps):
        current = rk4_step(L, current, dt)
    return current


def main():
    N = 5
    k_pairs = [(2, 4)]   # k=2 ↔ k=4 mirror partners (k=1↔5 also valid; pick interior pair)
    J = 1.0
    gamma = 0.05
    d = 2 ** N

    print("=" * 80)
    print("EQ-020 hardware-testable chiral mirror prediction")
    print("=" * 80)
    print(f"  N = {N}, J = {J}, γ = {gamma}, OBC XY chain + Z-dephasing")
    print()

    H = build_xy_chain_H(N, J)
    L = build_L(H, N, gamma)

    times = [0.0, 1.0, 2.0, 5.0, 10.0]

    for k, k_m in k_pairs:
        print(f"--- K_1-paired states: k = {k} ↔ k = {k_m} (N + 1 − k) ---")
        rho0_k = bonding_density(N, k, d)
        rho0_km = bonding_density(N, k_m, d)

        for t in times:
            rho_k_t = evolve_to_t(L, rho0_k, t)
            rho_km_t = evolve_to_t(L, rho0_km, t)
            X_k, Y_k, Z_k, P_k = expectation_per_site(rho_k_t, N)
            X_m, Y_m, Z_m, P_m = expectation_per_site(rho_km_t, N)

            # Predicted:
            #   X_m = +(-1)^i · X_k    (sublattice sign)
            #   Y_m = -(-1)^i · Y_k    (extra conjugation sign from E flip)
            #   Z_m = Z_k              (identical)
            #   P_m = P_k              (identical)
            sign = np.array([(-1) ** i for i in range(N)])
            X_diff = float(np.max(np.abs(X_m - (+sign) * X_k)))
            Y_diff = float(np.max(np.abs(Y_m - (-sign) * Y_k)))
            Z_diff = float(np.max(np.abs(Z_m - Z_k)))
            P_diff = float(np.max(np.abs(P_m - P_k)))

            print(f"  t = {t:5.1f}:")
            print(f"    max|⟨X_i⟩_m − (+(−1)^i)·⟨X_i⟩_k| = {X_diff:.2e}")
            print(f"    max|⟨Y_i⟩_m − (−(−1)^i)·⟨Y_i⟩_k| = {Y_diff:.2e}")
            print(f"    max|⟨Z_i⟩_m − ⟨Z_i⟩_k|           = {Z_diff:.2e}")
            print(f"    max|P_i_m − P_i_k|              = {P_diff:.2e}")


if __name__ == "__main__":
    main()
