"""F112 non-Hermitian extension: basis-enumeration proof attempt for the open identity.

Background:
  F112 typed Tier1Derived covers Hermitian H + bit_b-homogeneous c.
  F112 non-Hermitian Tier1Candidate extends to non-Hermitian H = H_re + i·H_im
  (with H_re, H_im both Hermitian). The proof reduces to the open identity:

      F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0  for any Hermitian H_re, H_im.

  My algebraic reduction (2026-05-26): F is real-bilinear in (H_re, H_im) and
  antisymmetric under H_re ↔ H_im exchange (since ⟨X, Y⟩^* = ⟨Y, X⟩ flips
  Im sign). So F is determined by its values on a basis of pairs of Hermitian
  operators. The Hermitian operator space is spanned by Pauli strings with
  real coefficients (4^N of them).

  Constructive proof strategy: if F(σ_α, σ_β) = 0 bit-exact for every pair of
  Pauli strings (σ_α, σ_β) at some N, then F ≡ 0 on the entire Hermitian
  operator space at that N (by basis spanning + bilinearity). This script
  enumerates the 4^N × 4^N pair grid at N=2 (16 × 16 = 256 pairs, of which
  136 are upper-triangular distinct) and reports.

  If F bit-exact 0 across all pairs → F112 non-Hermitian extension is proven
  at N=2 (or whatever N this runs at). The pattern would lift to general N
  by an inductive argument or by recognizing the structural symmetry that
  makes F vanish.

  If F has non-zero entries → counterexample found; the conjecture is wrong
  (or our derivation has a bug).

Method:
  For each ordered pair (α, β) of Pauli string indices:
    - Build σ_α, σ_β as 2^N × 2^N Hermitian matrices
    - Compute L_α := -i[σ_α, ·] and L_β := -i[σ_β, ·] in Pauli basis
    - Project each onto Π-conjugation -i eigenspace via the standard projector
    - Compute Frobenius inner product ⟨L_α,-i, L_β,-i⟩
    - Report Im part; flag if |Im| > 1e-10

Reports:
  - Max |Im|, mean |Im|, count of non-zero entries
  - Structural breakdown: by (bit_b parity of α, bit_b parity of β)
  - Conclusion: bit-exact 0 → identity holds at this N
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.symmetry import build_pi_full
from framework.pauli import _vec_to_pauli_basis_transform

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
LETTERS = ['I', 'X', 'Y', 'Z']


def pauli_string(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def build_L_H_pauli(H, N):
    """L = -i[H, ·] in Pauli basis."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def project_pi_eigenspace(M, Pi, target_eigenvalue):
    """Pi-conjugation eigenvalue projection: P_λ(M) = (1/4) Σ_k λ^{-k} Π^k M Π^{-k}."""
    Pi_inv = Pi.conj().T
    result = np.zeros_like(M)
    cur_Pi = np.eye(Pi.shape[0], dtype=complex)
    cur_Pi_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = (1.0 / (target_eigenvalue ** k)) / 4.0
        result = result + coef * (cur_Pi @ M @ cur_Pi_inv)
        cur_Pi = cur_Pi @ Pi
        cur_Pi_inv = cur_Pi_inv @ Pi_inv
    return result


def frobenius_inner(A, B):
    """⟨A, B⟩ = Tr(A† B)."""
    return complex(np.sum(A.conj() * B))


def bit_b_of_string(letters):
    return (sum(1 for L in letters if L == 'Y')
            + sum(1 for L in letters if L == 'Z')) % 2


def enumerate_F_on_basis(N, tol=1e-10):
    """Compute F(σ_α, σ_β) for every ordered pair of Pauli strings at chain length N.

    Returns (max_im, mean_im, count_nonzero, total, by_bitb_class_dict)."""
    Pi = build_pi_full(N)
    strings = list(product(LETTERS, repeat=N))
    print(f"N = {N}: enumerating {len(strings)}^2 = {len(strings)**2} ordered pairs ({len(strings)*(len(strings)+1)//2} distinct unordered)")

    # Pre-compute L_α,-i for every α
    print(f"  Pre-computing L_α,-i for {len(strings)} Pauli strings...")
    L_minus_i_per_alpha = {}
    for alpha in strings:
        sigma = pauli_string(alpha)
        L = build_L_H_pauli(sigma, N)
        L_minus_i_per_alpha[alpha] = project_pi_eigenspace(L, Pi, -1j)
    print(f"  Done. Computing F(α, β) for all pairs...")

    max_im = 0.0
    sum_im_abs = 0.0
    count_nonzero = 0
    total = 0
    nonzero_examples = []
    by_bitb_class = {(0, 0): [], (0, 1): [], (1, 0): [], (1, 1): []}

    # Enumerate UNORDERED pairs (alpha ≤ beta indexwise) to avoid F antisymmetry double count
    alpha_idx_map = {a: i for i, a in enumerate(strings)}
    for alpha in strings:
        for beta in strings:
            if alpha_idx_map[alpha] > alpha_idx_map[beta]:
                continue  # upper triangle only
            f_val = frobenius_inner(L_minus_i_per_alpha[alpha], L_minus_i_per_alpha[beta])
            im_val = float(f_val.imag)
            abs_im = abs(im_val)
            sum_im_abs += abs_im
            if abs_im > max_im:
                max_im = abs_im
            if abs_im > tol:
                count_nonzero += 1
                if len(nonzero_examples) < 10:
                    nonzero_examples.append((alpha, beta, im_val))
            total += 1
            bbk = (bit_b_of_string(alpha), bit_b_of_string(beta))
            by_bitb_class[bbk].append(abs_im)

    mean_im = sum_im_abs / max(total, 1)

    print(f"\n  --- N = {N} basis enumeration result ---")
    print(f"  Total distinct ordered pairs (upper triangle): {total}")
    print(f"  Max |Im⟨L_α,-i, L_β,-i⟩|: {max_im:.4e}")
    print(f"  Mean |Im|: {mean_im:.4e}")
    print(f"  Pairs with |Im| > {tol:.0e}: {count_nonzero}")
    if nonzero_examples:
        print(f"  First non-zero examples (BAD: F112 conjecture violated):")
        for a, b, im in nonzero_examples:
            print(f"    F({''.join(a)}, {''.join(b)}) = {im:+.4e}")
    else:
        print(f"  *** All {total} pair F-values are bit-exact 0 (< {tol:.0e}) ***")
        print(f"  → F ≡ 0 on Pauli-basis at N={N} (basis-spanning proof of the open identity at this N)")

    print(f"\n  Breakdown by bit_b parity class (alpha_bit_b, beta_bit_b):")
    for bbk in sorted(by_bitb_class.keys()):
        vals = by_bitb_class[bbk]
        if not vals:
            continue
        max_v = max(vals)
        print(f"    {bbk}: {len(vals)} pairs, max |Im| = {max_v:.4e}")

    return max_im, mean_im, count_nonzero, total


def main():
    print("F112 NON-HERMITIAN EXTENSION: basis-enumeration proof attempt")
    print("=" * 80)
    print(f"""
The open identity to prove or refute:

    F(H_re, H_im) := Im⟨L_{{H_re,-i}}, L_{{H_im,-i}}⟩ = 0   for any Hermitian H_re, H_im

F is real-bilinear in (H_re, H_im), antisymmetric under H_re ↔ H_im exchange
(since ⟨X, Y⟩^* = ⟨Y, X⟩ flips Im sign).

By bilinearity + basis spanning: if F(σ_α, σ_β) = 0 bit-exact for every pair of
Pauli strings (σ_α, σ_β) at some N, then F ≡ 0 on the entire Hermitian
operator space at that N.
""")
    for N in [2, 3]:
        enumerate_F_on_basis(N, tol=1e-10)
        print()

    print("=" * 80)
    print("If both N=2 and N=3 enumerations give max |Im| bit-exact 0:")
    print("  → F112 non-Hermitian extension proven at N=2 AND N=3 via basis-spanning argument")
    print("  → For general N, the pattern needs inductive lift OR a structural reason")
    print()
    print("If any non-zero entries surface:")
    print("  → The conjecture is wrong; F112 non-Hermitian extension fails on a constructive case")


if __name__ == '__main__':
    main()
