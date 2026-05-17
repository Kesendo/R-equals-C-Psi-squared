#!/usr/bin/env python3
"""F86b Dicke superposition Π²-odd closed-form: independent brute-force verifier.

Canonical statement and proof: docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md §Statement 2.

For ψ = (|D_n⟩+|D_{n+1}⟩)/√2: α_total = (1 − γ²)/2 with γ = ⟨ψ|X⊗N|ψ⟩,
yielding the three anchors 0 (mirror) / 3/8 (K-intermediate) / 1/2 (generic).

This script computes α_total two ways and checks bit-exact agreement:
(1) the closed form via X⊗N overlap; (2) a deliberate independent Pauli-enumeration
brute force, kept so the verification doesn't share code paths with the proof.
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
    """Independent witness: Π²-odd Frobenius² via direct 4^N Pauli enumeration.

    Deliberately shares no code path with alpha_total_closed_form so the
    verification is honest (an error in the closed form cannot mask itself).
    """
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
    for N in [3, 4, 5, 6]:
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
