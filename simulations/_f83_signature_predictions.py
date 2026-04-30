"""F83 anti-fraction signature predictions for hardware testing.

Designs a 4-Hamiltonian discriminator test that operationally distinguishes
F83's anti-fraction predictions across the trichotomy + Π²-class:

    Hamiltonian       F77 class            F83 anti-fraction r=hen/hod
    ----------------  -------------------  ----------------- ---------
    XX+YY              truly                n/a (M=0)           —
    XY+YX              soft (pure Π²-odd)    1/2               0.0
    YZ+ZY              soft (pure Π²-even)   0                 ∞
    XY+YZ              mixed                  1/6               1.0

For each Hamiltonian, this script computes the framework's predicted ⟨P⟩(t)
signature for all 16 (q0, q2) Pauli expectations on chain N=3, |+,−,+⟩
initial state, t=0.8, n_trotter=3 (matching the soft_break Marrakesh
2026-04-26 setup for direct comparison).

Two prediction models per Hamiltonian:
  - Continuous Lindblad (γ_Z=0.1, no T1)
  - Trotter n=3 (matching actual circuit, γ_Z=0.1, no T1)

The Trotter model fit Marrakesh's ⟨X₀Z₂⟩ to within 0.0014 for the existing
3 Hamiltonians (XX+YY, XY+YX, XX+XY); we expect the same precision for
the new 4th Hamiltonian (XY+YZ) plus the alternative case (YZ+ZY).

Hardware test design summary:
  - 4 Hamiltonians × 9 Pauli bases × 4096 shots = 36 circuits
  - ~3-5 minutes QPU on ibm_marrakesh path [48, 49, 50]
  - QPU budget: well within 180 min/year (about 2-3 % of annual budget)
  - Expected discrimination: F83 anti-fraction maps to specific
    Pauli-expectation patterns; differences between r=0, 1/6, 1/2 should
    be observable at >>30σ given the soft_break run's ±0.015 statistical
    error per ⟨P⟩ at 4096 shots.

Pipeline integration (when Tom decides to submit):
  Edit `D:\...\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\
       run_soft_break.py` CATEGORIES list to the 4 Hamiltonians below;
  run `python run_soft_break.py --hardware --backend ibm_marrakesh
       --path 48,49,50 --shots 4096`. Wait ~3 min, get JSON output. Compare
  with predictions in this script.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))

import framework as fw
from framework.lindblad import lindbladian_z_dephasing
from framework.pauli import _build_bilinear, ur_pauli, site_op


N = 3
J = 1.0
T_EVAL = 0.8
N_TROTTER = 3
GAMMA_Z = 0.1

# 4 Hamiltonians for F83 signature discrimination
CATEGORIES = [
    ('truly_unbroken',       [('X', 'X'), ('Y', 'Y')], 'truly',                 None),
    ('pi2_odd_pure',         [('X', 'Y'), ('Y', 'X')], 'soft (pure Π²-odd)',    0.5),
    ('pi2_even_nontruly',    [('Y', 'Z'), ('Z', 'Y')], 'soft (Π²-even non-truly)', 0.0),
    ('mixed_anti_one_sixth', [('X', 'Y'), ('Y', 'Z')], 'mixed (Π²-odd + Π²-even non-truly)', 1.0/6),
]

PAULIS = ['I', 'X', 'Y', 'Z']


def vec_F(M):
    return M.flatten('F')


def unvec_F(v, d):
    return v.reshape((d, d), order='F')


def trotter_lindblad_propagate(N, bonds, terms, n_trot, t, gamma_z, rho_0):
    """First-order Trotter circuit with per-step Z-dephasing."""
    delta_t = t / n_trot
    Id = np.eye(2 ** N, dtype=complex)

    # Per-step Trotter unitary
    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q, c) in terms:
        for (l, m) in bonds:
            ops = [ur_pauli('I')] * N
            ops[l] = ur_pauli(P); ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * c * delta_t * op_full) @ U_step

    # Per-step Z-dephasing channel (vec basis)
    L_deph = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        L_deph += gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    M_deph = expm(L_deph * delta_t)

    rho = rho_0
    for _ in range(n_trot):
        rho = U_step @ rho @ U_step.conj().T
        rho = unvec_F(M_deph @ vec_F(rho), 2 ** N)
    return rho


def two_qubit_pauli_obs_q0_q2(N, p0, p2):
    """Operator P_q0 ⊗ I_q1 ⊗ P_q2 on N=3 chain."""
    return np.kron(np.kron(ur_pauli(p0), np.eye(2)), ur_pauli(p2))


def predict_observables(category, terms, mode='trotter'):
    """Return dict of {(p0, p2): predicted ⟨P_q0 ⊗ I ⊗ P_q2⟩ at t=T_EVAL}."""
    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    bilinear = [(a, b, J) for (a, b) in terms]

    # Initial state |+,-,+>
    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi0 = np.kron(np.kron(ket_p, ket_m), ket_p)
    rho_0 = np.outer(psi0, psi0.conj())

    if mode == 'trotter':
        rho_t = trotter_lindblad_propagate(N, bonds, bilinear, N_TROTTER, T_EVAL, GAMMA_Z, rho_0)
    elif mode == 'continuous':
        H = _build_bilinear(N, bonds, bilinear)
        L = lindbladian_z_dephasing(H, [GAMMA_Z] * N)
        rho_t = unvec_F(expm(L * T_EVAL) @ vec_F(rho_0), 2 ** N)
    else:
        raise ValueError(f"unknown mode: {mode}")

    expectations = {}
    for p0 in PAULIS:
        for p2 in PAULIS:
            obs = two_qubit_pauli_obs_q0_q2(N, p0, p2)
            expectations[(p0, p2)] = float(np.real(np.trace(rho_t @ obs)))
    return expectations


def main():
    print("=" * 78)
    print("F83 SIGNATURE TEST — Hardware Prediction Table")
    print("=" * 78)
    print(f"  Setup: N={N} chain, J={J}, t={T_EVAL}, n_trotter={N_TROTTER}, γ_Z={GAMMA_Z}")
    print(f"  Initial state: |+,−,+⟩")
    print(f"  Target backend: ibm_marrakesh path [4, 5, 6] (top-ranked by 2026-04-30 calibration; rank 1/223)")
    print(f"  Continuity alternative: path [48, 49, 50] (April 26 soft_break, now rank 171/223)")
    print()

    # Compute predictions (Trotter, our anchor model from Marrakesh fit) for all 4
    print("Predicted ⟨P_q0 ⊗ I_q1 ⊗ P_q2⟩ (Trotter n=3, γ_Z=0.1):")
    print()
    header = f'  {"Pauli":<6}'
    for cat, _, _, _ in CATEGORIES:
        header += f' | {cat[:14]:>14}'
    print(header)
    print('  ' + '-' * (8 + 17 * len(CATEGORIES)))

    predictions = {cat: predict_observables(cat, terms, 'trotter')
                   for cat, terms, _, _ in CATEGORIES}

    for p0 in PAULIS:
        for p2 in PAULIS:
            if p0 == 'I' and p2 == 'I':
                continue
            row = f'  {p0},{p2:<5}'
            for cat, _, _, _ in CATEGORIES:
                v = predictions[cat][(p0, p2)]
                row += f' | {v:>+14.4f}'
            print(row)

    # F83 structural summary
    print()
    print("=" * 78)
    print("F83 anti-fraction predictions (closed form, no propagation needed)")
    print("=" * 78)
    print()
    chain = fw.ChainSystem(N=N)
    print(f'  {"Hamiltonian":<24} {"F77 class":<32} {"F83 anti-frac":>15} {"r":>8}')
    print('  ' + '-' * 80)
    for cat, terms, descr, _ in CATEGORIES:
        d = chain.predict_pi_decomposition(terms)
        anti = d['anti_fraction']
        r = d['r']
        r_str = f'{r:.3f}' if r != float('inf') else '∞'
        anti_str = 'n/a (M=0)' if d['M_sq'] < 1e-12 else f'{anti:.4f}'
        print(f'  {cat:<24} {descr:<32} {anti_str:>15} {r_str:>8}')

    # Discrimination signal
    print()
    print("=" * 78)
    print("Discriminating Pauli expectations (where the 4 Hamiltonians differ most)")
    print("=" * 78)
    print()
    # For each (p0, p2), compute spread = max(|pred|) − min(|pred|) across 4 H
    spreads = []
    for p0 in PAULIS:
        for p2 in PAULIS:
            if p0 == 'I' and p2 == 'I':
                continue
            vals = [predictions[cat][(p0, p2)] for cat, _, _, _ in CATEGORIES]
            spread = max(vals) - min(vals)
            spreads.append((p0, p2, spread, vals))

    spreads.sort(key=lambda x: -x[2])
    print(f'  Top 6 discriminating Paulis (largest spread across 4 H):')
    print(f'  {"Pauli":<6} {"truly":>10} {"odd":>10} {"even-nt":>10} {"mixed":>10} {"spread":>10}')
    print('  ' + '-' * 70)
    for p0, p2, spread, vals in spreads[:6]:
        labels = ['truly', 'odd', 'even-nt', 'mixed']
        line = f'  {p0},{p2:<5}'
        for v in vals:
            line += f' {v:>10.4f}'
        line += f' {spread:>10.4f}'
        print(line)

    print()
    print("=" * 78)
    print("Hardware test recipe")
    print("=" * 78)
    print(f"""
  Edit:
    D:\\Entwicklung\\Projekte\\.NET Projekte\\AIEvolution\\AIEvolution.UI\\
        experiments\\ibm_quantum_tomography\\run_soft_break.py

  Replace CATEGORIES (line ~107) with:

      CATEGORIES = [
          ('truly_unbroken',       [('X', 'X'), ('Y', 'Y')]),  # M = 0
          ('pi2_odd_pure',         [('X', 'Y'), ('Y', 'X')]),  # F83 anti = 1/2
          ('pi2_even_nontruly',    [('Y', 'Z'), ('Z', 'Y')]),  # F83 anti = 0
          ('mixed_anti_one_sixth', [('X', 'Y'), ('Y', 'Z')]),  # F83 anti = 1/6
      ]

  Then:
    cd "D:\\Entwicklung\\Projekte\\.NET Projekte\\AIEvolution\\AIEvolution.UI\\
        experiments\\ibm_quantum_tomography"
    python run_soft_break.py --hardware --backend ibm_marrakesh \\
                             --path 4,5,6 --shots 4096
        (or --path 48,49,50 for continuity with April 26 soft_break)

  Total: 36 circuits, ~3-5 min QPU. Output JSON in results/.

  Then compare measured ⟨P_q0 P_q2⟩ values against the prediction table
  above. Categories should be discriminable at the discrimination spread
  shown (typically 0.5+ for top observables, statistical error 0.015 at
  4096 shots → >>30σ separation).

  Validates: F83 anti-fraction structurally observable on hardware via
  characteristic Pauli-expectation patterns per Π²-class.
""")


if __name__ == "__main__":
    main()
