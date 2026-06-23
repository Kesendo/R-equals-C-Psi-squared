#!/usr/bin/env python3
"""Hardware test design: distinguishing 3 truly-unbroken from 19 soft-broken
Hamiltonians on IBM hardware via tomography.

This is the FIRST framework-based experimental design. All Hamiltonians,
operators, and predicted observables are constructed via framework.py
primitives (ur_pauli, ur_xyz, lindbladian_z_dephasing). The pre-framework
scripts in `D:\\...\\ibm_quantum_tomography\\` use ad-hoc Pauli definitions
and don't call framework. This script does.

Goal
----
On IBM hardware, distinguish three Hamiltonian categories:
  - 1 truly-unbroken: XX + YY  (both-parity-even, palindrome holds exactly)
  - 1 soft-broken:    XY + YX  (eigenvalues pair, eigenvectors don't)
  - 1 hard-broken:    XX + XY  (eigenvalues don't pair)

Setup
-----
  N = 3 chain (qubits 0, 1, 2)
  Two bonds: (0,1) and (1,2), same Hamiltonian per bond
  γ_l = γ per site (Z-dephasing from hardware native sources)
  Initial state: prepared via gates (TBD: |+−+⟩ probes the X-basis Néel sector)
  Trotter evolution for time t, then 9-Pauli tomography on (q0, q2)

Discriminating observable
-------------------------
The framework predicts that the time-resolved 2-qubit reduced density matrix
ρ_{0,2}(t) carries different signatures across the three categories:

  - Truly-unbroken: Π-symmetric structure preserved at operator level.
    Specific matrix elements of ρ_{0,2} (those connecting Π-paired
    eigenvectors of L) follow the palindrome relation predictably.

  - Soft-broken: same EIGENVALUES but scrambled EIGENVECTORS. The
    matrix elements that "should" follow Π-symmetry don't. Specifically,
    the off-diagonal coupling between weight sectors of ρ_{0,2} differs
    from the truly-unbroken case.

  - Hard-broken: spectrum itself non-palindromic. Different decay
    rates for different modes; the entire decoherence curve has a
    different shape.

This script computes ρ_{0,2}(t) predictions for each category via
framework.lindbladian_z_dephasing + matrix exponential, and identifies
which 2-qubit Pauli expectation values most strongly differ between
the three categories.

Output: a measurement protocol with predicted hardware observables and
their discriminating-power ranking.
"""
import math
import sys
from pathlib import Path

import numpy as np

import framework as fw

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ----------------------------------------------------------------------
# Test Hamiltonians (one from each category)
# ----------------------------------------------------------------------

CATEGORIES = [
    ('truly_unbroken', [('X', 'X'), ('Y', 'Y')]),  # XX + YY
    ('soft_broken',    [('X', 'Y'), ('Y', 'X')]),  # XY + YX
    ('hard_broken',    [('X', 'X'), ('X', 'Y')]),  # XX + XY
]


def build_test_hamiltonian(N, bonds, terms, J=1.0):
    """Build H = J Σ_bond Σ_term σ_term[0]_i σ_term[1]_j on each bond.

    Uses framework._build_bilinear under the hood.
    """
    bilinear_terms = [(t[0], t[1], J) for t in terms]
    return fw._build_bilinear(N, bonds, bilinear_terms)


def evolve_density_matrix(rho_0, L_vec, t, d):
    """Evolve ρ_0 under L for time t. Returns ρ(t) as 2^N × 2^N matrix.

    L_vec is the column-stack vec Liouvillian. We exponentiate it.
    """
    rho_0_vec = rho_0.flatten('F')  # column-stack
    U = np.linalg.matrix_power
    # exp(L_vec * t) via eigendecomposition (sufficient for small d)
    expL = np_expm(L_vec * t)
    rho_t_vec = expL @ rho_0_vec
    rho_t = rho_t_vec.reshape((d, d), order='F')
    return rho_t


def np_expm(M):
    """Matrix exponential via scipy."""
    from scipy.linalg import expm
    return expm(M)


def reduced_density_matrix(rho, sites_to_keep, N):
    """Trace out all sites except sites_to_keep."""
    sites_to_keep = sorted(sites_to_keep)
    out = rho.reshape([2] * (2 * N))
    ket_axes = list(range(N))
    bra_axes = list(range(N, 2 * N))
    for j in range(N - 1, -1, -1):
        if j in sites_to_keep:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        lo, hi = sorted((a_k, a_b))
        for k in range(N):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    k = len(sites_to_keep)
    d_sub = 2 ** k
    remaining_axes_ket = [ket_axes[i] for i in sites_to_keep]
    remaining_axes_bra = [bra_axes[i] for i in sites_to_keep]
    current_order = remaining_axes_ket + remaining_axes_bra
    perm = np.argsort(current_order)
    out = np.transpose(out, perm)
    return out.reshape(d_sub, d_sub)


def two_qubit_pauli_expectations(rho_2):
    """Compute all 16 Pauli expectations ⟨P_a P_b⟩ for 2-qubit ρ.

    Returns dict {('I','I'): val, ('I','X'): val, ..., ('Z','Z'): val}.
    """
    out = {}
    for label_a in ['I', 'X', 'Y', 'Z']:
        for label_b in ['I', 'X', 'Y', 'Z']:
            P_a = fw.ur_pauli(label_a)
            P_b = fw.ur_pauli(label_b)
            P_ab = np.kron(P_a, P_b)
            val = float(np.trace(P_ab @ rho_2).real)
            out[(label_a, label_b)] = val
    return out


def main():
    N = 3
    gamma = 0.1
    gamma_l = [gamma] * N
    Sigma_gamma = sum(gamma_l)
    bonds = [(0, 1), (1, 2)]
    J = 1.0
    t_eval = 0.8  # similar to receiver-engineering experiment
    d = 2 ** N

    print("=" * 90)
    print("FIRST framework-based hardware test design")
    print(f"N={N}, J={J}, γ={gamma} per site, t_eval={t_eval}")
    print("Three Hamiltonians, three categories. Initial state: |+−+⟩ X-Néel.")
    print("=" * 90)

    # Initial state: X-Néel |+−+⟩ (probes the multi-exc sector, similar to Z⊗N partnership)
    psi_init = fw.ur_xneel(N, sign_pattern=[+1, -1, +1])
    rho_init = np.outer(psi_init, psi_init.conj())
    print(f"\nInitial state ‖ψ‖² = {np.linalg.norm(psi_init)**2:.6f}, "
          f"Tr(ρ) = {np.trace(rho_init).real:.6f}")

    # Compute the three Hamiltonians and evolutions
    results = {}
    for category, terms in CATEGORIES:
        H = build_test_hamiltonian(N, bonds, terms, J=J)
        L = fw.lindbladian_z_dephasing(H, gamma_l)
        rho_t = evolve_density_matrix(rho_init, L, t_eval, d)
        rho_02 = reduced_density_matrix(rho_t, [0, 2], N)
        expectations = two_qubit_pauli_expectations(rho_02)
        # Also compute palindrome residual norm
        M_residual = fw.palindrome_residual(L, Sigma_gamma, N)
        residual_norm = float(np.linalg.norm(M_residual))
        results[category] = {
            'H_terms': terms,
            'rho_02': rho_02,
            'expectations': expectations,
            'palindrome_residual': residual_norm,
        }

    print(f"\n{'Category':>16s}  {'H terms':>14s}  {'palindrome residual':>22s}")
    for category, data in results.items():
        terms_str = '+'.join(t[0] + t[1] for t in data['H_terms'])
        print(f"{category:>16s}  {terms_str:>14s}  {data['palindrome_residual']:>22.4e}")

    # Find discriminating observables: 2-qubit Pauli expectations that differ
    # most across the three categories.
    print("\n" + "=" * 90)
    print(f"2-qubit Pauli expectations ⟨P_a^0 P_b^2⟩ at t={t_eval}, on sites (0, 2):")
    print("=" * 90)
    print(f"\n{'Pauli (P0, P2)':>16s}  {'truly_unbroken':>16s}  {'soft_broken':>14s}  "
          f"{'hard_broken':>14s}  {'max-min spread':>16s}")

    pauli_keys = [(a, b) for a in ['I', 'X', 'Y', 'Z'] for b in ['I', 'X', 'Y', 'Z']]
    discriminators = []
    for key in pauli_keys:
        vals = [results[cat]['expectations'][key] for cat, _ in CATEGORIES]
        spread = max(vals) - min(vals)
        if abs(spread) > 1e-6:
            discriminators.append((key, vals, spread))

    discriminators.sort(key=lambda x: -abs(x[2]))
    for (key, vals, spread) in discriminators[:10]:
        key_str = f"({key[0]}, {key[1]})"
        print(f"{key_str:>16s}  {vals[0]:>16.4f}  {vals[1]:>14.4f}  "
              f"{vals[2]:>14.4f}  {spread:>16.4f}")

    # Truly vs soft discrimination
    print("\n" + "=" * 90)
    print("Specifically: truly_unbroken vs soft_broken (the framework's NEW prediction)")
    print("=" * 90)
    print(f"\n{'Pauli (P0, P2)':>16s}  {'truly_unbroken':>16s}  {'soft_broken':>14s}  "
          f"{'difference':>14s}")
    truly_vs_soft = []
    for key in pauli_keys:
        v_truly = results['truly_unbroken']['expectations'][key]
        v_soft = results['soft_broken']['expectations'][key]
        diff = v_truly - v_soft
        if abs(diff) > 1e-6:
            truly_vs_soft.append((key, v_truly, v_soft, diff))
    truly_vs_soft.sort(key=lambda x: -abs(x[3]))
    for (key, v_t, v_s, d_) in truly_vs_soft[:8]:
        key_str = f"({key[0]}, {key[1]})"
        print(f"{key_str:>16s}  {v_t:>16.4f}  {v_s:>14.4f}  {d_:>14.4f}")

    print()
    print("=" * 90)
    print("Hardware protocol (proposed):")
    print("=" * 90)
    print("""
    1. Backend: ibm_marrakesh (Heron r2) or ibm_kingston (Heron r2).
    2. Path: linear 3-qubit chain via rank_paths (top-rated).
    3. Initial state: prepare |+−+⟩ via Hadamard on qubit 0,
       Hadamard + Z on qubit 1, Hadamard on qubit 2 (3 SQ gates).
    4. Trotter steps: 3 (first-order Trotter, dt = 0.267).
    5. For each of 3 Hamiltonians:
       - Apply Trotter circuit
       - Run 9-Pauli tomography on (q0, q2)
       - 8192 shots × 9 bases = 73728 shots per Hamiltonian
    6. Total: 3 × 9 = 27 circuits, ~3-5 QPU minutes.

    Discriminating observable (predicted from this simulation):
    The top discriminators between truly_unbroken and soft_broken should
    show difference ≥ 0.05 in absolute Pauli expectation value (above
    typical IBM tomography noise of ~0.02 at 8192 shots).
    """)


if __name__ == "__main__":
    main()
