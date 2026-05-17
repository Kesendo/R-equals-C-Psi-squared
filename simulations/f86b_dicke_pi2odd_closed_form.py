#!/usr/bin/env python3
"""F86b Dicke superposition Π²-odd closed-form (Tier 1 derived 2026-05-17).

Closes the open analytical proof of PROOF_F86B_UNIVERSAL_SHAPE.md:129
("The K-intermediate 3/8 anchor ... analytical proof remains open").

THEOREM. For |ψ⟩ = (|D_n⟩+|D_{n+1}⟩)/√2 on N qubits, the total Π²-odd
Frobenius² of ρ = |ψ⟩⟨ψ| equals:

    α_total = (1 − ⟨ψ|X⊗N|ψ⟩²) / 2
            = 2·|c_+|²·|c_-|²
            = (1+γ)(1-γ)/2     where γ = ⟨ψ|X⊗N|ψ⟩

Three structural cases:
  ⟨ψ|X⊗N|ψ⟩ = 1   →  α = 0    (Dicke-mirror, N odd, 2n+1=N)
  ⟨ψ|X⊗N|ψ⟩ = 1/2 →  α = 3/8  (Dicke-K-intermediate, N even, n ∈ {N/2-1, N/2})
  ⟨ψ|X⊗N|ψ⟩ = 0   →  α = 1/2  (generic)

PROOF (4 lines):
  1. Decompose ψ = c_+ ψ_+ + c_- ψ_- where ψ_± are X⊗N-eigenstates (±1 eigenvalue).
  2. ρ = |ψ⟩⟨ψ| has Π²-EVEN diagonal blocks |c_±|²|ψ_±⟩⟨ψ_±| + Π²-ODD off-diagonal
     c_+c_-*|ψ_+⟩⟨ψ_-| + h.c. (X⊗N-conjugation eigenvalue tracks bit_b parity).
  3. Π²-ODD Frobenius² = 2|c_+c_-|² (the two off-diagonal terms are HS-orthogonal).
  4. |c_±|² = (1 ± γ)/2 with γ = ⟨ψ|X⊗N|ψ⟩, so α_total = 2·(1+γ)(1-γ)/4 = (1-γ²)/2.

The three-anchor structure follows from X⊗N action on Dicke pairs: X⊗N maps
|D_k⟩ → |D_{N-k}⟩, so X⊗N|ψ⟩ = (|D_{N-n}⟩+|D_{N-n-1}⟩)/√2. The overlap
⟨ψ|X⊗N|ψ⟩ = (1/2)(δ_{n,N-n} + δ_{n,N-n-1} + δ_{n+1,N-n} + δ_{n+1,N-n-1})
hits {0, 1/2, 1} depending on which Kronecker deltas align with the (n, n+1) pair.

This is THE F50-style mechanism — orthogonal-symmetry decomposition replacing
opaque Krawtchouk computation — applied to F86b. Same approach as F50 max-spin
Dicke endpoint ladder rungs (commit 5523171).
"""
from __future__ import annotations

import sys
import numpy as np
from itertools import combinations, product

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

PAULI = {
    "I": np.eye(2, dtype=np.complex128),
    "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
}


def pauli_string(N, letters):
    op = np.array([[1]], dtype=np.complex128)
    for L in letters:
        op = np.kron(op, PAULI[L])
    return op


def dicke_state(N, k):
    d = 2 ** N
    v = np.zeros(d, dtype=np.complex128)
    for positions in combinations(range(N), k):
        idx = 0
        for p in positions:
            idx |= 1 << (N - 1 - p)
        v[idx] = 1.0
    return v / np.linalg.norm(v)


def x_total(N):
    op = np.array([[1]], dtype=np.complex128)
    for _ in range(N):
        op = np.kron(op, PAULI["X"])
    return op


def alpha_total_closed_form(N, n):
    """Closed-form: α_total = (1 − ⟨ψ|X⊗N|ψ⟩²) / 2."""
    D_n = dicke_state(N, n)
    D_n1 = dicke_state(N, n + 1)
    psi = (D_n + D_n1) / np.sqrt(2)
    X_N = x_total(N)
    gamma = np.real(np.conj(psi) @ (X_N @ psi))
    return (1 - gamma ** 2) / 2


def alpha_total_brute(N, n):
    """Brute-force: Π²-odd Frobenius² via Pauli enumeration."""
    D_n = dicke_state(N, n)
    D_n1 = dicke_state(N, n + 1)
    psi = (D_n + D_n1) / np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    d = 2 ** N
    alpha = 0.0
    for L in product("IXYZ", repeat=N):
        op = pauli_string(N, list(L))
        c = np.trace(op.conj().T @ rho)
        bit_b = sum(1 for cc in L if cc in "YZ")
        if all(cc == "I" for cc in L):
            continue
        if bit_b % 2 == 1:
            alpha += abs(c) ** 2 / d
    return alpha


def overlap_gamma(N, n):
    """⟨ψ|X⊗N|ψ⟩ for Dicke superposition."""
    D_n = dicke_state(N, n)
    D_n1 = dicke_state(N, n + 1)
    psi = (D_n + D_n1) / np.sqrt(2)
    X_N = x_total(N)
    return np.real(np.conj(psi) @ (X_N @ psi))


def main():
    print("=" * 78)
    print("F86b Dicke superposition Π²-odd closed-form verification")
    print("Theorem: α_total = (1 − ⟨ψ|X⊗N|ψ⟩²) / 2  for ψ = (|D_n⟩+|D_{n+1}⟩)/√2")
    print("=" * 78)
    print()
    print(f"  {'N':>3} {'n':>3} {'category':>15} {'⟨ψ|X⊗N|ψ⟩':>13} {'α_closed':>10} {'α_brute':>10} {'ok':>3}")
    print("  " + "-" * 70)

    all_pass = True
    for N in [3, 4, 5, 6, 7]:
        for n in range(N):
            gamma = overlap_gamma(N, n)
            alpha_cf = alpha_total_closed_form(N, n)
            alpha_br = alpha_total_brute(N, n)
            if 2 * n + 1 == N:
                cat = "Dicke-mirror"
            elif (N % 2 == 0) and n in [N // 2 - 1, N // 2]:
                cat = "Dicke-K-int"
            else:
                cat = "generic"
            ok = abs(alpha_cf - alpha_br) < 1e-10
            if not ok:
                all_pass = False
            marker = "✓" if ok else "✗"
            print(f"  {N:>3} {n:>3} {cat:>15} {gamma:>13.6f} {alpha_cf:>10.6f} {alpha_br:>10.6f} {marker:>3}")
    print()
    print(f"  ALL CHECKS PASS: {all_pass}")
    print()
    print("Closed-form values at the three anchors:")
    print("  γ = 0   → α = (1−0)/2 = 1/2  (generic)")
    print("  γ = 1/2 → α = (1−1/4)/2 = 3/8  (K-intermediate)")
    print("  γ = 1   → α = (1−1)/2 = 0     (mirror, ψ is X⊗N-eigenstate)")


if __name__ == "__main__":
    main()
