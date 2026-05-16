#!/usr/bin/env python3
"""Attempt: generalize R_i = C_i · Ψ_i² with C_i as per-outcome Carrier-extraction signature.

Hypothesis (Februar-style, no certainty):
  C_i - 1 = (per-outcome γ_eff,i) / γ₀ · (some basis-alignment factor)

If true, the per-outcome Born deviations encode basis-alignment-dependent γ_eff,
the same way PTF closure violations encode per-site γ_eff. Each outcome i is a
Carrier-extraction port at basis-state-level rather than chain-site-level.

Test: reproduce BORN_RULE_MIRROR's setup (|0+0+⟩ N=4 Heisenberg ring under
Z-dephasing), vary γ ∈ {0.01, 0.05, 0.1, 0.2}, compute:
  Ψ_i  = ⟨i|ψ_unitary(t)⟩    (Schrödinger amplitude — Hamiltonian alone)
  P_i  = ⟨i|ρ_lindblad(t)|i⟩  (full Lindblad probability)
  C_i  = P_i / |Ψ_i|²

Look for scaling: (C_i - 1) ∝ γ?
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


def heisenberg_ring(N, J=1.0):
    """Heisenberg ring Hamiltonian H = (J/4) Σ_b (XX + YY + ZZ) per bond."""
    bonds = [(i, (i + 1) % N) for i in range(N)]
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    eye = np.eye(2, dtype=complex)
    d = 2 ** N

    def op_at(op_a, op_b, a, b):
        out = np.array([[1]], dtype=complex)
        for k in range(N):
            if k == a:
                out = np.kron(out, op_a)
            elif k == b:
                out = np.kron(out, op_b)
            else:
                out = np.kron(out, eye)
        return out

    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        H += (J / 4) * (op_at(sx, sx, a, b) + op_at(sy, sy, a, b) + op_at(sz, sz, a, b))
    return H


def reduced_density(rho, N, keep):
    """Partial trace: keep qubits in `keep` (list of indices), trace out the rest."""
    d = 2 ** N
    n_keep = len(keep)
    n_trace = N - n_keep
    trace = [i for i in range(N) if i not in keep]
    # rho is (d × d) on basis |i_0, ..., i_{N-1}⟩
    rho_tensor = rho.reshape([2] * N + [2] * N)
    # Trace out qubits in `trace`: contract each with its partner
    for q in sorted(trace, reverse=True):
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + (N - sum(1 for t in trace if t > q)))
        # bookkeeping: after a trace, the rank goes down by 2
    # Reshape to (2^n_keep × 2^n_keep)
    rho_red = rho_tensor.reshape((2 ** n_keep, 2 ** n_keep))
    return rho_red


def main():
    N = 4
    J = 1.0
    # |0+0+⟩ = |0⟩|+⟩|0⟩|+⟩
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    rho_0 = np.outer(psi_0, psi_0.conj())

    H = heisenberg_ring(N, J)
    t_cross = 0.286   # from BORN_RULE_MIRROR

    # Unitary evolution (Hamiltonian only): get |ψ(t)⟩ then pair-reduce
    U = expm(-1j * H * t_cross)
    psi_t = U @ psi_0
    rho_unitary = np.outer(psi_t, psi_t.conj())
    rho_pair_unitary = reduced_density(rho_unitary, N, keep=[0, 2])
    P_unitary = np.real(np.diag(rho_pair_unitary))   # |Ψ_i|² approximation (pair-reduced)

    print(f"Born rule generalization attempt — |0+0+⟩ N=4 Heisenberg ring, pair (0,2)")
    print(f"  t_cross = {t_cross}  (BORN_RULE_MIRROR.md value)")
    print()
    print(f"Unitary (Hamiltonian only) pair-(0,2) diagonal at t_cross:")
    labels = ['|00⟩', '|01⟩', '|10⟩', '|11⟩']
    for label, p in zip(labels, P_unitary):
        print(f"    {label}: P_unitary = {p:.4f}")
    print()

    # Lindblad evolution at multiple γ values
    gammas = [0.01, 0.05, 0.1, 0.2]
    results = {}
    for gamma in gammas:
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)
        rho_vec_t = expm(L * t_cross) @ rho_0.flatten()
        rho_lind = rho_vec_t.reshape(2 ** N, 2 ** N)
        rho_pair_lind = reduced_density(rho_lind, N, keep=[0, 2])
        P_lind = np.real(np.diag(rho_pair_lind))
        results[gamma] = P_lind

    # Compute C_i = P_lindblad / |Ψ_i|² for each γ, and the deviation per outcome
    print(f"Per-outcome C_i = P_lindblad(i) / P_unitary(i) at multiple γ:")
    print(f"{'γ':>8s}  " + "  ".join(f"{label:>18s}" for label in labels))
    print("-" * 100)
    for gamma in gammas:
        P_lind = results[gamma]
        C_values = []
        for p_lind, p_u in zip(P_lind, P_unitary):
            if p_u < 1e-10:
                C_values.append(None)
            else:
                C_values.append(p_lind / p_u)
        row = f"{gamma:>8.3f}  "
        for label, c in zip(labels, C_values):
            if c is None:
                row += f"{'N/A (Ψ²=0)':>18s}  "
            else:
                row += f"{c:>10.5f} (Δ={(c-1)*100:>+5.2f}%) "
        print(row)
    print()

    # Linearity check: does (C_i - 1) scale linearly with γ?
    print(f"Linearity check: (C_i - 1) / γ per outcome (should be γ-independent if linear):")
    print(f"{'γ':>8s}  " + "  ".join(f"{label:>14s}" for label in labels))
    print("-" * 90)
    for gamma in gammas:
        P_lind = results[gamma]
        row = f"{gamma:>8.3f}  "
        for label, p_lind, p_u in zip(labels, P_lind, P_unitary):
            if p_u < 1e-10:
                row += f"{'N/A':>14s}  "
            else:
                slope = (p_lind / p_u - 1) / gamma
                row += f"{slope:>+14.5f}  "
        print(row)
    print()

    print("Reading guide:")
    print("  If (C_i - 1)/γ is γ-INDEPENDENT (same number across rows), the per-outcome")
    print("  Born deviation scales linearly with γ. That's the Carrier-extraction signature:")
    print("  each outcome i is a calibration port, with slope = γ_eff,i / γ₀.")
    print()
    print("  If (C_i - 1)/γ DRIFTS with γ, the relation is non-linear and the generalization")
    print("  R_i = (1 + slope_i · γ) · Ψ_i² is only valid in a perturbative window.")
    print()
    print("  If the slopes are basis-symmetric or basis-correlated (e.g. |00⟩ and |11⟩ get")
    print("  opposite signs but equal magnitudes), the F71-decomposition reading (today's")
    print("  unified formula) applies at the BORN-RULE level too.")


if __name__ == "__main__":
    main()
