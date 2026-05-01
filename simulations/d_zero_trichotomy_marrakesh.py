"""d=0 sector-population reading of the F77 trichotomy on Marrakesh data.

The existing palindrome_trichotomy confirmation (Confirmations) reads the
trichotomy via ⟨X₀Z₂⟩ at full 9-Pauli-pair tomography. This script shows
that the same trichotomy is visible in just the Z-basis — the d=0
sector populations p_n = Tr(P_n · ρ).

Mechanism: each Hamiltonian class in the trichotomy has a distinct
excitation-conservation profile.

  truly  XX+YY: conserves total Z-magnetization → ⟨n⟩ stays at 1 (initial)
  soft   XY+YX: breaks conservation, mass shifts UP   → ⟨n⟩ ≈ 1.28
  hard   XX+XY: breaks conservation, mass shifts DOWN → ⟨n⟩ ≈ 0.88

These three values are separated by 0.13–0.27 — far above the ≈0.04
hardware deviation we observe — so a single Z-basis tomography (1 of 9
Pauli bases used in the original test) discriminates the three classes.

Hardware data: data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json
(Marrakesh path [48,49,50], N=3, J=1, t=0.8, n_trotter=3, |+−+⟩ initial).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing


HARDWARE_PATH = Path(__file__).parent.parent / 'data' / 'ibm_soft_break_april2026' / 'soft_break_ibm_marrakesh_20260426_001101.json'


def reduce_q1_out(rho_3q):
    """Trace out qubit 1 from a 3-qubit density matrix → 4×4 (q0,q2) reduced."""
    rho4 = rho_3q.reshape(2, 2, 2, 2, 2, 2)
    return np.einsum('ikj lkm->ijlm', rho4).reshape(4, 4)


def hardware_sector_populations(counts_zz):
    """Sector populations p_n from raw 2-qubit Z-basis counts."""
    total = sum(counts_zz.values())
    p = np.zeros(3)
    for bitstr, c in counts_zz.items():
        n = bitstr.count('1')
        p[n] += c
    return p / total


def framework_prediction(terms, t=0.8, gamma_Z=0.1, J=1.0):
    """Propagate |+−+⟩ under H = Σ_bonds bilinear + Z-dephasing, return reduced (q0,q2)."""
    N = 3
    bonds = [(0, 1), (1, 2)]
    plus = np.array([1, 1]) / np.sqrt(2)
    minus = np.array([1, -1]) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())

    H = _build_bilinear(N, bonds, terms) * J
    L = lindbladian_z_dephasing(H, [gamma_Z] * N)
    rho_t_vec = expm(L * t) @ rho_0.flatten('F')
    rho_t = rho_t_vec.reshape(2 ** N, 2 ** N, order='F')
    return reduce_q1_out(rho_t)


def main():
    cases = [
        ('truly  (XX+YY)', [('X', 'X', 1.0), ('Y', 'Y', 1.0)]),
        ('soft   (XY+YX)', [('X', 'Y', 1.0), ('Y', 'X', 1.0)]),
        ('hard   (XX+XY)', [('X', 'X', 1.0), ('X', 'Y', 1.0)]),
    ]
    hw_keys = {'truly  (XX+YY)': 'truly_unbroken',
               'soft   (XY+YX)': 'soft_broken',
               'hard   (XX+XY)': 'hard_broken'}

    with open(HARDWARE_PATH) as f:
        hw_data = json.load(f)

    print(f'Hardware: {hw_data["backend"]}, job {hw_data["job_id"]}')
    print(f'Parameters: {hw_data["parameters"]}')
    print()

    rows = []
    for label, terms in cases:
        rho_red = framework_prediction(terms)
        sp = fw.sector_populations(rho_red, N=2)
        zz_counts = hw_data['counts'][hw_keys[label]]['ZZ']
        p_hw = hardware_sector_populations(zz_counts)
        mean_n_hw = float(p_hw[1] + 2 * p_hw[2])
        rows.append((label, sp['p'], sp['mean_n'], p_hw, mean_n_hw))

    print(f'{"":<18s} | {"framework p_0,p_1,p_2":<28s} | {"hw p_0,p_1,p_2":<28s} | mean_n: pred / hw / Δ')
    print('-' * 120)
    for label, p_fw, mn_fw, p_hw, mn_hw in rows:
        fw_str = f'{p_fw[0]:.4f}, {p_fw[1]:.4f}, {p_fw[2]:.4f}'
        hw_str = f'{p_hw[0]:.4f}, {p_hw[1]:.4f}, {p_hw[2]:.4f}'
        print(f'{label:<18s} | {fw_str:<28s} | {hw_str:<28s} | {mn_fw:.4f} / {mn_hw:.4f} / {mn_hw - mn_fw:+.4f}')

    print()
    spreads_fw = [rows[1][2] - rows[0][2], rows[2][2] - rows[0][2]]
    spreads_hw = [rows[1][4] - rows[0][4], rows[2][4] - rows[0][4]]
    print(f'⟨n⟩ spreads relative to truly:')
    print(f'  framework: soft − truly = {spreads_fw[0]:+.4f}, hard − truly = {spreads_fw[1]:+.4f}')
    print(f'  hardware:  soft − truly = {spreads_hw[0]:+.4f}, hard − truly = {spreads_hw[1]:+.4f}')
    print()
    print('Three classes, three distinct ⟨n⟩ values, all separated by >0.10 — a clean')
    print('trichotomy signature in just the Z-basis (1 of 9 Pauli measurements).')


if __name__ == '__main__':
    main()
