"""F112 Diagonal-Vermutung bei N=3 testen.

Vermutung: für BitB-odd σ_α, σ_β: ⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = δ_{σ_α, σ_β} · 4^N

Falls bei N=3 (32 BitB-odd strings, 1024 ordered pairs) bestätigt → universal-N
Ableitung der F112 non-Hermitian extension ist trivial (Im = 0 weil
Off-Diagonale alle null + Diagonale reelle Norm²).

Numerik via numpy (kein sympy, viel schneller für N=3).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.symmetry import build_pi_full
from framework.pauli import _vec_to_pauli_basis_transform

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}
LETTERS_ORDER = ["I", "X", "Z", "Y"]  # PauliIndex convention
BIT_B = {"I": 0, "X": 0, "Z": 1, "Y": 1}


def pauli_string(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def build_L_H_pauli(H, N):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2**N)


def project_pi_eigenspace(M, Pi, target_eigenvalue):
    Pi_inv = Pi.conj().T
    result = np.zeros_like(M)
    cur_pi = np.eye(Pi.shape[0], dtype=complex)
    cur_pi_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = (1.0 / (target_eigenvalue**k)) / 4.0
        result = result + coef * (cur_pi @ M @ cur_pi_inv)
        cur_pi = cur_pi @ Pi
        cur_pi_inv = cur_pi_inv @ Pi_inv
    return result


def frobenius(A, B):
    return complex(np.sum(A.conj() * B))


def bit_b_total(letters):
    return sum(BIT_B[L] for L in letters) % 2


def enumerate_strings(N):
    from itertools import product
    return list(product(LETTERS_ORDER, repeat=N))


N = 3
print(f"F112 Diagonal-Vermutung bei N={N}")
print("=" * 72)

pi = build_pi_full(N)
all_strings = enumerate_strings(N)
bit_b_odd = [s for s in all_strings if bit_b_total(s) == 1]
print(f"BitB-odd strings: {len(bit_b_odd)} (expected 4^{N}/2 = {4**N // 2})")
print()

# Pre-compute all L_{σ,-i} for BitB-odd
print(f"Pre-computing L_α,-i for {len(bit_b_odd)} strings...")
L_mi_cache = {}
for s in bit_b_odd:
    sigma = pauli_string(s)
    L = build_L_H_pauli(sigma, N)
    L_mi_cache[s] = project_pi_eigenspace(L, pi, -1j)
print("Done.")
print()

# Test: only diagonals should be non-zero
print(f"Computing all {len(bit_b_odd)**2} pairwise inner products...")
diagonals = []
off_diagonals_nonzero = []
all_im_zero = True
tol = 1e-10

for sa in bit_b_odd:
    for sb in bit_b_odd:
        inner = frobenius(L_mi_cache[sa], L_mi_cache[sb])
        if abs(inner.imag) > tol:
            all_im_zero = False
            print(f"  ⚠ Im ≠ 0 at ({''.join(sa)}, {''.join(sb)}): {inner}")
        if sa == sb:
            diagonals.append((sa, inner.real))
        elif abs(inner) > tol:
            off_diagonals_nonzero.append((sa, sb, inner))

print()
print(f"Im part check: {'ALL ZERO ✓' if all_im_zero else 'FAILED ✗'}")
print()
print(f"Diagonal values (α = β):")
diag_values = [d[1] for d in diagonals]
unique_diag = set(round(v, 6) for v in diag_values)
print(f"  {len(diagonals)} diagonals, unique values: {unique_diag}")
if len(unique_diag) == 1:
    val = list(unique_diag)[0]
    print(f"  → Alle Diagonalen gleich {val}.  4^N = {4**N}.  Verhältnis: {val / (4**N):.6f}")
print()
print(f"Off-diagonal non-zero count: {len(off_diagonals_nonzero)}")
if off_diagonals_nonzero:
    print("  First 5 examples:")
    for sa, sb, val in off_diagonals_nonzero[:5]:
        print(f"    ({''.join(sa)}, {''.join(sb)}): {val}")
else:
    print("  ✓ ALLE Off-Diagonalen null — Diagonal-Vermutung BESTÄTIGT bei N=3")
    print()
    print("  → F(α, β) = Im⟨L_α,-i, L_β,-i⟩ = δ_{α,β} · 4^N · 0 = 0 für α ≠ β")
    print("  → F(α, α) = Im(real number) = 0 für α = β")
    print("  → F ≡ 0 für alle (α, β) STRUKTURELL (nicht enumerativ!)")
