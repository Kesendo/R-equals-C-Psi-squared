#!/usr/bin/env python3
"""Quick (Q, K) map of the per-outcome Born deviation Δ_i at γ = 0.05 fixed.

γ fixed at 0.05. Then Q = J/γ = 20·J and K = γt = 0.05·t. Sweep (J, t),
tabulate Δ_i for each basis state of pair (0,2) of |0+0+⟩ N=4 Heisenberg.

For fun: see the shape of Δ_i(Q, K).
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


def delta_at(J, gamma, t, N=4):
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

    Delta = np.array([np.nan if P_u[i] < 1e-10 else (P_l[i] / P_u[i] - 1)
                       for i in range(4)])
    return Delta, P_u, P_l


def main():
    GAMMA = 0.05
    print(f"Δ_i(Q, K) map at γ = {GAMMA} fixed, |0+0+⟩ N=4 Heisenberg ring, pair (0,2)")
    print(f"  Q = J/γ = {1/GAMMA:.0f}·J ;  K = γt = {GAMMA}·t")
    print()

    # Reasonably-spaced grid
    J_values = [0.5, 1.0, 2.0, 5.0, 10.0]
    t_values = [0.1, 0.286, 0.5, 1.0, 2.0]

    for outcome_idx, label in enumerate(['|00⟩', '|01⟩', '|10⟩', '|11⟩']):
        print(f"━━━ Δ_{label} (in %) ━━━")
        print(f"{'':>8s} | " + "  ".join(f"K={GAMMA*t:>6.4f}" for t in t_values))
        print(f"{'':>8s}   " + "  ".join(f"  (t={t:>4.2f})" for t in t_values))
        print("-" * 90)
        for J in J_values:
            Q = J / GAMMA
            row = f"Q={Q:>5.0f} | "
            for t in t_values:
                Delta, _, _ = delta_at(J=J, gamma=GAMMA, t=t)
                d = Delta[outcome_idx]
                if np.isnan(d):
                    row += f"{'  Ψ²≈0':>8s}  "
                else:
                    row += f"{d*100:>+7.3f}%  "
            print(row)
        print()

    print("Reading:")
    print("  Each cell = Δ_i = (P_lindblad / P_unitary − 1) at the given (Q, K).")
    print("  Positive = decoherence FAVORS this outcome at (Q, K).")
    print("  Negative = decoherence SUPPRESSES this outcome at (Q, K).")
    print()
    print("  Look for: is Δ_i monotonic in K?  In Q?  Are there sign-flips?")
    print("  Does Δ_|00⟩ = −Δ_|11⟩ (F71-anti-symmetric)?")


if __name__ == "__main__":
    main()
