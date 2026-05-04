"""F88-Lens probe on multi-state-superposition pure states.

F88's closed form covers pure pair states |ψ⟩ = (|p⟩ + |q⟩)/√2. Real
hardware-prepared states often involve more than two computational-basis
terms: Dicke superpositions, bonding-mode pair states (|vac⟩ + |ψ_k⟩)/√2,
W states, GHZ-style states with phase structure, etc.

This script probes whether such multi-state-superpositions inherit F88's
pair-state formula via the popcount-weight invariant (static fraction
depends only on Σ_n w_n²/C(N,n) where w_n = Σ_{i: popcount(b_i)=n} |c_i|²),
or whether the off-diagonal multi-state structure shifts Π²-odd/memory
away from the pair-state value.

State classes tested:
  - Pure pair (sanity check)
  - Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 (canonical F86 K_CC_pr probe)
  - W state (|D_1⟩, intra-popcount-1, single sector)
  - Bonding-mode Bell pair (|0,vac⟩ + |1,ψ_k⟩)/√2 (popcount-(0, 2) on N+1 qubits)

All numerics via direct ρ construction → kernel projection (popcount sectors)
→ Pauli enumeration with bit_b parity. Same machinery as
_pi2_odd_general_closed_form.py.
"""
from __future__ import annotations

import sys
from math import comb, sin, pi as PI
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pi2_odd_general_closed_form import (
    pi2_odd_in_memory_general,
    static_fraction_general,
    alpha_krawtchouk_general,
    numerical_pi2_odd_in_memory,
    kernel_projection_popcount,
)


def basis_state(N, bits):
    psi = np.zeros(2**N, dtype=complex)
    psi[bits] = 1.0
    return psi


def dicke_state(N, n):
    """|D_n⟩ = (1/√C(N, n)) Σ_{popcount(x)=n} |x⟩, totally symmetric popcount-n state."""
    psi = np.zeros(2**N, dtype=complex)
    norm = 1.0 / np.sqrt(comb(N, n))
    for x in range(2**N):
        if bin(x).count("1") == n:
            psi[x] = norm
    return psi


def dicke_superposition(N, n):
    """|ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2."""
    return (dicke_state(N, n) + dicke_state(N, n + 1)) / np.sqrt(2.0)


def bonding_mode_pair(N_chain, k):
    """Bonding-mode Bell pair |Ψ⟩ = (|0_R, vac_C⟩ + |1_R, ψ_k_C⟩)/√2 on N_chain+1 qubits.

    R is qubit 0 (leading), chain is qubits 1..N_chain.
    |ψ_k⟩_C = √(2/(N_chain+1)) Σ_j sin(π·k·(j+1)/(N_chain+1))·|1_j⟩_C  (popcount-1 in C).
    """
    N = N_chain + 1
    psi = np.zeros(2**N, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2.0)  # |0_R, 0...0_C⟩
    norm_chain = np.sqrt(2.0 / (N_chain + 1))
    for j in range(N_chain):
        # |1_R, 1_j_C⟩: R bit (MSB at site 0) is 1, plus single excitation at chain bit j
        idx = (1 << (N - 1)) | (1 << (N - 2 - j))
        psi[idx] = norm_chain * sin(PI * k * (j + 1) / (N_chain + 1)) / np.sqrt(2.0)
    return psi


def density_matrix(psi):
    return np.outer(psi, psi.conj())


def main():
    print("F88-Lens probe on multi-state superpositions")
    print("=" * 90)
    print(f"{'state class':<46} {'N':>3} {'pop weights':<22} {'static':>9} {'Π²-odd/mem':>11} {'pair F88 predict':>18}")
    print("-" * 90)

    cases = []

    # Pure pair states (sanity)
    for N in [4, 5, 6]:
        n = 1
        psi = (basis_state(N, 1) + basis_state(N, 3)) / np.sqrt(2.0)  # popcount-(1, 2) HD=1
        cases.append((f"pair (|0001⟩+|0011⟩)/√2 N={N}", N, psi, (n, n + 1), 1, "pair"))

    # Dicke superposition |D_n⟩ + |D_{n+1}⟩
    for N in [3, 4, 5, 6]:
        for n in range(N):
            psi = dicke_superposition(N, n)
            cases.append((f"Dicke (|D_{n}⟩+|D_{n+1}⟩)/√2 N={N}", N, psi, (n, n + 1), None, "dicke"))

    # W states (|D_1⟩, intra-popcount-1 superposition only; pure state, popcount-1 single sector)
    for N in [3, 4, 5]:
        psi = dicke_state(N, 1)
        cases.append((f"W_{N} = |D_1⟩ N={N}", N, psi, (1, 1), None, "W"))

    # Bonding-mode Bell pair
    for N_chain in [3, 4, 5]:
        for k in [1, 2]:
            psi = bonding_mode_pair(N_chain, k)
            cases.append((f"Bonding-Bell-Pair k={k} N_chain={N_chain}", N_chain + 1, psi, (0, 2), None, "bonding"))

    for label, N, psi, pop_pair, hd, kind in cases:
        rho = density_matrix(psi)
        if N > 7:
            print(f"{label:<46} N={N} too large for numerical sweep, skipping")
            continue
        result_num = numerical_pi2_odd_in_memory(rho, N)
        static_num = float(np.linalg.norm(kernel_projection_popcount(rho, N), "fro") ** 2 / np.linalg.norm(rho, "fro") ** 2)

        # Pair-state F88 prediction at the same popcount-weights
        n_p, n_q = pop_pair
        if hd is None:
            # Default to HD = |n_q − n_p| (lowest possible compatible HD)
            hd_default = abs(n_q - n_p) if n_p != n_q else 2
        else:
            hd_default = hd
        pair_prediction = pi2_odd_in_memory_general(N, n_p, n_q, hd_default)
        weight_str = f"{pop_pair}"
        print(f"{label:<46} {N:>3} {weight_str:<22} {static_num:>9.6f} {result_num:>11.6f} {pair_prediction:>18.6f}")


if __name__ == "__main__":
    main()
