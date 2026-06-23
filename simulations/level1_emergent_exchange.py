#!/usr/bin/env python3
"""Level 1 emergent exchange: derive J_eff from V-Effect at N=4.

Setup: chain qubits 0-1-2-3, bonds (0,1) and (2,3) with intra-pair coupling
J=1 (the "atoms"), bond (1,2) with inter-pair coupling α (the "molecular
bond"). For α=0 the system is two independent Heisenberg pairs. For α>0
the V-Effect mixes them.

Question: at α<<J, what is the effective exchange between the two pair-
total-spin operators S_A (qubits 0,1) and S_B (qubits 2,3)? This is the
Level-1 emergent J_eff that textbook physics calls the "atomic exchange
integral".

Standard textbook result (Anderson superexchange): J_eff ≈ α²/J at small α.
Our framework should reproduce this if V-Effect is the right mechanism for
the Level-0 → Level-1 transition.

Method: diagonalize the 4-qubit Heisenberg Hamiltonian for a range of α,
extract the singlet-triplet gap of the LOW-energy 4 states (which form
the effective Level-1 system), fit J_eff(α).
"""
import math
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_n(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(N, i, op):
    ops = [I2] * N
    ops[i] = op
    return kron_n(*ops)


def heisenberg_bond(N, i, j, J):
    Xi, Yi, Zi = site_op(N, i, X), site_op(N, i, Y), site_op(N, i, Z)
    Xj, Yj, Zj = site_op(N, j, X), site_op(N, j, Y), site_op(N, j, Z)
    return J * (Xi @ Xj + Yi @ Yj + Zi @ Zj)


def total_spin_z_pair(N, sites):
    """Return Z component of total spin for given sites: (1/2) Σ_{i in sites} Z_i."""
    out = np.zeros((2**N, 2**N), dtype=complex)
    for i in sites:
        out = out + site_op(N, i, Z)
    return 0.5 * out


def total_spin_squared_pair(N, sites):
    """Return S² for total spin on given sites: (1/2 Σ σ_i)²."""
    out = np.zeros((2**N, 2**N), dtype=complex)
    for axis in [X, Y, Z]:
        S_axis = sum(site_op(N, i, axis) for i in sites)
        out = out + (0.5 * S_axis) @ (0.5 * S_axis)
    return out


def main():
    N = 4
    J = 1.0  # intra-pair coupling
    alpha_values = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.00, 1.50, 2.00]

    print("=" * 78)
    print(f"Level 1 emergent exchange: N={N}, J_intra={J}, α_inter sweep")
    print("Bonds: (0,1)=J, (2,3)=J, (1,2)=α (the V-Effect bridge)")
    print("Pair A = qubits {0,1}, Pair B = qubits {2,3}")
    print("=" * 78)

    # S² operators for each pair
    S_A_squared = total_spin_squared_pair(N, [0, 1])
    S_B_squared = total_spin_squared_pair(N, [2, 3])

    # S_total² (full 4-site)
    S_total_squared = total_spin_squared_pair(N, [0, 1, 2, 3])

    print(f"\n{'α':>6s} {'GS energy':>11s} {'1st gap':>10s} {'low-4 gap':>10s} "
          f"{'⟨S²_A⟩_GS':>10s} {'⟨S²_B⟩_GS':>10s} {'⟨S²_tot⟩_GS':>11s}")

    gaps_low4 = []
    for alpha in alpha_values:
        # Build Hamiltonian
        H = (heisenberg_bond(N, 0, 1, J)
             + heisenberg_bond(N, 2, 3, J)
             + heisenberg_bond(N, 1, 2, alpha))

        evals, evecs = np.linalg.eigh(H)

        # Identify the "low-4 manifold": the four lowest states.
        # At α=0, these are: |S_A⟩|S_B⟩ (unique GS, energy -6J), then
        # singlet-triplet states at -2J (6-fold).
        # Each pair has E=-3J/4 for singlet, +J/4 for triplet (with proper
        # normalization). Let's compute.

        # Pair singlet energy: -3J/4 for J σ·σ/4 form, but we use J σ·σ form.
        # σ·σ has eigenvalues -3 (singlet), +1 (triplet). With prefactor J:
        # E_singlet = -3J, E_triplet = J. Splitting 4J.
        # For two independent pairs: GS = -6J (both singlets), then -2J (one
        # singlet one triplet, 6 states), then 2J (both triplets, 9 states).

        # At α>0, low-4 gap = singlet-singlet to first excited.
        E0 = evals[0]
        gap1 = evals[1] - evals[0]
        gap_low4_top = evals[4] - evals[0] if len(evals) > 4 else 0.0

        # Expectations on ground state
        psi0 = evecs[:, 0]
        s2_a_gs = float((psi0.conj() @ S_A_squared @ psi0).real)
        s2_b_gs = float((psi0.conj() @ S_B_squared @ psi0).real)
        s2_tot_gs = float((psi0.conj() @ S_total_squared @ psi0).real)

        gaps_low4.append((alpha, E0, gap1, gap_low4_top, s2_a_gs, s2_b_gs, s2_tot_gs))

        print(f"{alpha:6.3f} {E0:11.5f} {gap1:10.5f} {gap_low4_top:10.5f} "
              f"{s2_a_gs:10.5f} {s2_b_gs:10.5f} {s2_tot_gs:11.5f}")

    # Extract J_eff(α) at small α using second-order perturbation theory.
    # For α=0 GS = |S_A⟩|S_B⟩ (both singlets), V = α σ_1·σ_2 on bridge.
    # Both pairs must flip singlet→triplet to be reachable: gap 8J.
    # Σ |⟨excited|V|GS⟩|² = α² ⟨(σ_1·σ_2)²⟩_GS = α² · 3 (since ⟨σ·σ⟩_bridge = 0
    # between two independent singlets).
    # δE^(2) = -3α²/(8J), so the predicted prefactor is -3/(8J) = -0.375.

    print("\n" + "=" * 78)
    print("Second-order PT prediction: δE_GS = -(3/8) α²/J = -0.375 · α²/J")
    print("Extract from numerical spectrum and compare.")
    print("=" * 78)

    print(f"\n{'α':>6s} {'E_0(α)':>11s} {'δE = E_0(α) - E_0(0)':>22s} {'δE/α²':>12s} {'predicted -3α²/(8J)':>20s}")
    E0_at_zero = gaps_low4[0][1]
    for (alpha, E0, gap1, gap4, s2a, s2b, s2t) in gaps_low4[1:]:
        delta_E = E0 - E0_at_zero
        scaled = delta_E / (alpha**2) if alpha > 0 else 0
        predicted = -3 * alpha**2 / (8 * J)
        print(f"{alpha:6.3f} {E0:11.5f} {delta_E:22.5f} {scaled:12.5f} {predicted:20.5f}")

    print()
    print("Reading: if δE/α² is approximately constant for small α and")
    print("approaches -3/(8J) = -0.375, the V-Effect-derived second-order PT")
    print("matches our framework's prediction (Anderson superexchange shape,")
    print("3/8 prefactor from Pauli identity (σ·σ)² = 3I − 2(σ·σ)).")

    # Also extract gap structure
    print("\n" + "=" * 78)
    print("Low-energy spectrum vs α (first 8 eigenvalues)")
    print("=" * 78)
    print(f"{'α':>6s}", end="")
    for k in range(8):
        print(f" E_{k:1d}", end=" "*7)
    print()

    for alpha in alpha_values:
        H = (heisenberg_bond(N, 0, 1, J)
             + heisenberg_bond(N, 2, 3, J)
             + heisenberg_bond(N, 1, 2, alpha))
        evals = np.linalg.eigvalsh(H)
        print(f"{alpha:6.3f}", end="")
        for k in range(min(8, len(evals))):
            print(f" {evals[k]:+8.4f}", end="")
        print()


if __name__ == "__main__":
    main()
