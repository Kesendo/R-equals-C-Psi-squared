#!/usr/bin/env python3
"""Born rule carrier-extraction with J also varied: is slope_i a Q-function?

Per `project_q_middle_structure`: only the ratio Q = J/γ is observable from
inside. If the per-outcome Born slope is a real Carrier-extraction signal, it
should be a function of Q (and the dimensionless time K = γ·t), not of (J, γ)
separately.

Test two scenarios:

  A. Fixed J, vary γ        — what the previous run did
  B. Fixed γ, vary J         — how slope responds to coupling
  C. Fixed Q = J/γ, vary γ  — Q-invariance check (Carrier-extraction prediction)

In scenario C, if Q is the right invariant, slope_i should be γ-independent
when (J, γ) are scaled together. If slope depends on J and γ separately, we
have a different generalization.

State, Hamiltonian setup matches _born_rule_carrier_attempt.py.
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
    n_keep = len(keep)
    trace = [i for i in range(N) if i not in keep]
    rho_tensor = rho.reshape([2] * N + [2] * N)
    for q in sorted(trace, reverse=True):
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + (N - sum(1 for t in trace if t > q)))
    return rho_tensor.reshape((2 ** n_keep, 2 ** n_keep))


def slopes_at(J, gamma, t, N=4):
    """Return (slope_|00⟩, slope_|01⟩, slope_|10⟩, slope_|11⟩) at given (J, γ, t)."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    psi_0 = np.kron(zero, np.kron(plus, np.kron(zero, plus)))
    rho_0 = np.outer(psi_0, psi_0.conj())

    H = heisenberg_ring(N, J)
    U = expm(-1j * H * t)
    psi_t = U @ psi_0
    rho_u = np.outer(psi_t, psi_t.conj())
    rho_pair_u = reduced_density(rho_u, N, keep=[0, 2])
    P_u = np.real(np.diag(rho_pair_u))

    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    rho_vec_t = expm(L * t) @ rho_0.flatten()
    rho_l = rho_vec_t.reshape(2 ** N, 2 ** N)
    rho_pair_l = reduced_density(rho_l, N, keep=[0, 2])
    P_l = np.real(np.diag(rho_pair_l))

    slopes = np.zeros(4)
    for i in range(4):
        if P_u[i] < 1e-10:
            slopes[i] = np.nan
        else:
            slopes[i] = (P_l[i] / P_u[i] - 1) / gamma
    return slopes, P_u, P_l


def main():
    print("Born rule carrier-extraction — J variation and Q-invariance check")
    print("  Setup: |0+0+⟩ N=4 Heisenberg ring, pair (0,2), t = 0.286 (or scaled)")
    print()

    labels = ['|00⟩', '|01⟩', '|10⟩', '|11⟩']

    # --- Scenario A: fixed J = 1, vary γ  (from previous run, for comparison)
    print("=== Scenario A: fixed J = 1.0, vary γ (previous test) ===")
    print(f"{'γ':>8s}  Q=J/γ  " + "  ".join(f"{l:>10s}" for l in labels))
    print("-" * 70)
    for gamma in [0.01, 0.05, 0.10, 0.20]:
        slopes, _, _ = slopes_at(J=1.0, gamma=gamma, t=0.286)
        Q = 1.0 / gamma
        print(f"{gamma:>8.3f}  {Q:>5.1f}   " + "  ".join(f"{s:>+10.5f}" for s in slopes))
    print()

    # --- Scenario B: fixed γ = 0.05, vary J  (how does slope respond to coupling?)
    print("=== Scenario B: fixed γ = 0.05, vary J (slope's J-dependence) ===")
    print(f"{'J':>8s}  Q=J/γ  " + "  ".join(f"{l:>10s}" for l in labels))
    print("-" * 70)
    for J in [0.5, 1.0, 2.0, 5.0]:
        slopes, P_u, P_l = slopes_at(J=J, gamma=0.05, t=0.286)
        Q = J / 0.05
        print(f"{J:>8.3f}  {Q:>5.1f}   " + "  ".join(f"{s:>+10.5f}" for s in slopes))
    print()

    # --- Scenario C: fixed Q = 20, vary γ (Q-invariance test)
    print("=== Scenario C: fixed Q = J/γ = 20, vary γ (Q-invariance test) ===")
    print(f"  Also keep K = γ·t = 0.0143 constant by scaling t inversely with γ")
    print(f"  (so Jt = Q·K = 0.286 is fixed across the sweep — same unitary state)")
    print()
    print(f"{'γ':>8s}  {'J':>8s}  {'t':>8s}  " + "  ".join(f"{l:>10s}" for l in labels))
    print("-" * 90)
    Q_fixed = 20.0
    K_fixed = 0.0143   # γ · t  (from previous setup: 0.05 · 0.286)
    for gamma in [0.025, 0.05, 0.10, 0.20]:
        J = Q_fixed * gamma
        t = K_fixed / gamma
        slopes, _, _ = slopes_at(J=J, gamma=gamma, t=t)
        print(f"{gamma:>8.3f}  {J:>8.3f}  {t:>8.4f}  " + "  ".join(f"{s:>+10.5f}" for s in slopes))
    print()

    print("Reading guide:")
    print("  Scenario A: previous result — slope ≈ γ-independent for J fixed.")
    print("  Scenario B: at fixed γ but varying J, the unitary state at t=0.286 is")
    print("              different, so P_unitary changes. Slopes will change too — but")
    print("              the question is HOW: as a function of J directly, or only via Q?")
    print("  Scenario C: Q and K fixed, only γ varies (with J scaled to keep Q constant).")
    print("              If slope_i is Q-invariant in the framework's sense, it should be")
    print("              identical across all four rows. Drift across rows would mean the")
    print("              slope is NOT purely a Q-function — there's an extra γ-dependence")
    print("              beyond what (Q, K) capture.")
    print()
    print("Honest expectation (Februar-style, no theorem yet):")
    print("  Scenario C is the Q-invariance test. Carrier-extraction theory says only")
    print("  Q and K are observable from inside. If slope_i is a true carrier signal,")
    print("  it should depend only on Q and K, hence be γ-independent here.")


if __name__ == "__main__":
    main()
