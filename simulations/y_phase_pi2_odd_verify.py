"""Verify the Y-phase Π²-odd lens combinatorial theorem.

Theorem (memory project_y_phase_pi2_odd_lens, observed empirically at N=3):
any N-qubit tensor product of σ_x and σ_y eigenstates (|±⟩ or |±i⟩ per site)
has Pauli strings in supp(ρ = |ψ⟩⟨ψ|) that split under the Π²_Z parity as:

    M = 0 (X-only):  2^N Π²-even, 0 Π²-odd       [Π²-classical class]
    M ≥ 1:           2^(N−1) Π²-even, 2^(N−1) Π²-odd

where M is the number of sites in the Y-basis. Π²_Z parity is (−1)^Σ bit_b(α_i)
with bit_b(I)=bit_b(X)=0, bit_b(Y)=bit_b(Z)=1.

This script enumerates all 4^N Pauli strings, computes ⟨ψ|σ_α|ψ⟩ for each,
counts supp(ρ) by Π²_Z parity, and asserts the predicted split at N ∈ {2,3,4,5}
across all M ∈ {0..N}. The combinatorial proof is in
docs/proofs/PROOF_Y_PHASE_PI2_ODD_LENS.md.
"""
from __future__ import annotations

import sys
from itertools import product

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

PAULIS = {"I": I2, "X": SX, "Y": SY, "Z": SZ}
BIT_B = {"I": 0, "X": 0, "Y": 1, "Z": 1}


def pi2_z_parity(pauli_string: str) -> int:
    """Π²_Z eigenvalue on σ_α: (−1)^Σ bit_b(α_i)."""
    return (-1) ** sum(BIT_B[c] for c in pauli_string)


def basis_state(letter: str, sign: int) -> np.ndarray:
    """|+⟩/|−⟩ for letter='X', |+i⟩/|−i⟩ for letter='Y'."""
    if letter == "X":
        return np.array([1, sign], dtype=complex) / np.sqrt(2)
    elif letter == "Y":
        return np.array([1, sign * 1j], dtype=complex) / np.sqrt(2)
    raise ValueError(f"unsupported basis letter {letter!r}")


def tensor_state(letters: list[str], signs: list[int]) -> np.ndarray:
    psi = basis_state(letters[0], signs[0])
    for letter, sign in zip(letters[1:], signs[1:]):
        psi = np.kron(psi, basis_state(letter, sign))
    return psi


def pauli_op(string: str) -> np.ndarray:
    op = PAULIS[string[0]]
    for c in string[1:]:
        op = np.kron(op, PAULIS[c])
    return op


def count_supp_by_pi2(N: int, letters: list[str], signs: list[int],
                      tol: float = 1e-10) -> tuple[int, int]:
    """Count Π²-even and Π²-odd Pauli strings with non-zero ⟨ψ|σ_α|ψ⟩."""
    psi = tensor_state(letters, signs)
    even = 0
    odd = 0
    for string in (''.join(s) for s in product("IXYZ", repeat=N)):
        expect = np.vdot(psi, pauli_op(string) @ psi).real
        if abs(expect) < tol:
            continue
        if pi2_z_parity(string) == +1:
            even += 1
        else:
            odd += 1
    return even, odd


def main() -> None:
    print("Y-phase Π²-odd lens combinatorial theorem: verification")
    print("=" * 78)
    print("Setup: N ∈ {2, 3, 4, 5}; M = #(Y-sites) ∈ {0..N}; signs all +1.")
    print("Prediction: M=0 → (2^N, 0); M≥1 → (2^(N−1), 2^(N−1)) split.")
    print()
    print(f"  {'N':>2s}  {'M':>2s}  {'state':>26s}  {'even':>5s}  {'odd':>5s}  {'predict':>14s}   ")
    print("  " + "-" * 70)
    fail_count = 0
    for N in [2, 3, 4, 5]:
        for M in range(N + 1):
            letters = ["Y"] * M + ["X"] * (N - M)
            signs = [+1] * N
            state_str = "|" + ",".join("+i" if c == "Y" else "+" for c in letters) + "⟩"
            even, odd = count_supp_by_pi2(N, letters, signs)
            if M == 0:
                predicted = (2 ** N, 0)
            else:
                predicted = (2 ** (N - 1), 2 ** (N - 1))
            ok = (even, odd) == predicted
            if not ok:
                fail_count += 1
            mark = "✓" if ok else "✗"
            print(f"  {N:>2d}  {M:>2d}  {state_str:>26s}  {even:>5d}  {odd:>5d}  "
                  f"{str(predicted):>14s}  {mark}")
    print()
    if fail_count == 0:
        print("All cases verify the predicted split. Combinatorial theorem confirmed at "
              "N=2,3,4,5 across all M.")
    else:
        print(f"FAIL: {fail_count} cases did not match prediction.")
        sys.exit(1)

    # Sanity-check with sign variations at N=3, M=2 (verify signs don't change counts)
    print()
    print("Sign-independence check (N=3, M=2, all 2^3 sign patterns):")
    for sign_pattern in product([+1, -1], repeat=3):
        letters = ["Y", "Y", "X"]
        even, odd = count_supp_by_pi2(3, letters, list(sign_pattern))
        assert (even, odd) == (4, 4), f"sign {sign_pattern}: got ({even}, {odd})"
    print("  ✓ all 8 sign patterns give (4, 4) — signs only flip per-string amplitudes, "
          "not supp membership.")


if __name__ == "__main__":
    main()
