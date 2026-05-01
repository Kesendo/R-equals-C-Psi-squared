"""Klein-class trajectories of (C, Ψ) for Bell+ under three pure Hamiltonians.

The seed-to-fruit question: how does the fold-catastrophe parameter space
(C, Ψ) project onto the Π² Klein-class trichotomy?

Hypothesis: each pure Klein class of H drives a qualitatively distinct
trajectory shape in (C, Ψ)-space.

  truly  XX+YY (Π²-even matched, M=0): F1-algebra closes, CΨ follows
                                       the idealized fold path.
  odd    XY+YX (pure Π²-odd):          M_anti = L_{H_odd}, dynamics IS
                                       recirculation, CΨ re-pumped.
  even   YZ+ZY (Π²-even non-truly):    M_sym ≠ 0, standing-wave hold,
                                       no drive but no recovery.

Setup: 2-qubit Bell+ initial state, J=1, γ_Z=0.1, t in [0, 5].
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing


def concurrence(rho):
    """Wootters formula: C = max(0, λ₁ − λ₂ − λ₃ − λ₄) with λᵢ ordered."""
    sy = np.array([[0, -1j], [1j, 0]])
    SySy = np.kron(sy, sy)
    rho_tilde = SySy @ rho.conj() @ SySy
    R = rho @ rho_tilde
    eigs = np.linalg.eigvals(R)
    eigs = np.real(eigs)
    eigs[eigs < 0] = 0
    sqrt_eigs = np.sort(np.sqrt(eigs))[::-1]
    return float(max(0.0, sqrt_eigs[0] - sqrt_eigs[1] - sqrt_eigs[2] - sqrt_eigs[3]))


def psi_l1(rho):
    """Ψ = ℓ₁(ρ) / (d² − 1) where ℓ₁ = sum of |off-diagonal entries|."""
    d = rho.shape[0]
    off = rho - np.diag(np.diag(rho))
    return float(np.sum(np.abs(off))) / (d * d - 1)


def trajectory(rho_init, terms, t_grid, J=1.0, gamma_Z=0.1, N=2):
    bonds = [(0, 1)]
    H = _build_bilinear(N, bonds, terms) * J
    L = lindbladian_z_dephasing(H, [gamma_Z] * N)
    rho_vec_init = rho_init.flatten('F')
    Cs, Psis, CPsis = [], [], []
    for t in t_grid:
        rho_t_vec = expm(L * t) @ rho_vec_init
        rho_t = rho_t_vec.reshape(2 ** N, 2 ** N, order='F')
        rho_t = 0.5 * (rho_t + rho_t.conj().T)
        C = concurrence(rho_t)
        Psi = psi_l1(rho_t)
        Cs.append(C)
        Psis.append(Psi)
        CPsis.append(C * Psi)
    return np.array(Cs), np.array(Psis), np.array(CPsis)


def topology_signature(C_traj, Psi_traj, t_grid):
    """Classify the trajectory topology (monotone / oscillatory / pump-collapse)."""
    C = np.array(C_traj)
    Psi = np.array(Psi_traj)
    n_revivals = int(np.sum((np.diff(np.sign(np.diff(C))) > 0)))
    Psi_pump = float(Psi.max() / max(Psi[0], 1e-9))
    C_zero_hits = int(np.sum(C < 1e-6))
    return {
        'n_C_revivals': n_revivals,
        'Psi_max_over_init': Psi_pump,
        'n_C_zero_hits': C_zero_hits,
        'C_max_after_init': float(np.max(C[1:])),
        'Psi_max_overall': float(Psi.max()),
        'CPsi_max': float(np.max(C * Psi)),
    }


def main():
    cases = [
        ('truly  (XX+YY)', [('X', 'X', 1.0), ('Y', 'Y', 1.0)]),
        ('odd    (XY+YX)', [('X', 'Y', 1.0), ('Y', 'X', 1.0)]),
        ('even   (YZ+ZY)', [('Y', 'Z', 1.0), ('Z', 'Y', 1.0)]),
    ]

    plus = np.array([1, 1]) / np.sqrt(2)
    minus = np.array([1, -1]) / np.sqrt(2)
    zero = np.array([1, 0])
    one = np.array([0, 1])

    initial_states = {
        'Bell+ (|00⟩+|11⟩)/√2':   np.array([1, 0, 0, 1]) / np.sqrt(2),
        'Bell- (|00⟩-|11⟩)/√2':   np.array([1, 0, 0, -1]) / np.sqrt(2),
        'Ψ+    (|01⟩+|10⟩)/√2':   np.array([0, 1, 1, 0]) / np.sqrt(2),
        'Ψ-    (|01⟩-|10⟩)/√2':   np.array([0, 1, -1, 0]) / np.sqrt(2),
        '|+,+⟩ (product, C=0)':   np.kron(plus, plus),
        '|+,-⟩ (product, C=0)':   np.kron(plus, minus),
        '|0,+⟩ (asym, C=0)':      np.kron(zero, plus),
    }

    t_grid = np.linspace(0, 5.0, 101)

    print(f'Klein-class trajectory topology across initial states (J=1, γ_Z=0.1, N=2)')
    print('=' * 100)
    print()
    print(f'{"initial":<25s} | {"H class":<14s} | {"C revivals":>10s} | '
          f'{"Ψ_max/Ψ_0":>10s} | {"C zero hits":>11s} | {"CΨ_max":>7s}')
    print('-' * 100)

    for state_label, psi in initial_states.items():
        psi = psi.astype(complex)
        if abs(np.linalg.norm(psi)) < 1e-9:
            continue
        psi /= np.linalg.norm(psi)
        rho_init = np.outer(psi, psi.conj())
        for case_label, terms in cases:
            C_traj, Psi_traj, _ = trajectory(rho_init, terms, t_grid)
            sig = topology_signature(C_traj, Psi_traj, t_grid)
            psi_pump = sig['Psi_max_over_init'] if sig['Psi_max_over_init'] < 1e3 else float('inf')
            print(f'{state_label:<25s} | {case_label:<14s} | '
                  f'{sig["n_C_revivals"]:>10d} | '
                  f'{psi_pump:>10.2f} | '
                  f'{sig["n_C_zero_hits"]:>11d} | '
                  f'{sig["CPsi_max"]:>7.4f}')
        print()


if __name__ == '__main__':
    main()
