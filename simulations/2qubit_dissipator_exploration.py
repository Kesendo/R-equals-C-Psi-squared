"""Exploration of 2-qubit dissipators and their F81 violation structure.

This is a Tier 1 + Tier 2 research script (not a closed framework theorem yet).
Counterpart to F82/F84 (single-site amplitude damping).

Findings (commit f9820d4 or later):

  1. **Pauli-Channel Cancellation Lemma extends to 2-qubit Pauli channels.**
     For any Pauli pair c, c' ∈ {I, X, Y, Z}, the Hermitian dissipator
     D[c ⊗ c'] is Π²-symmetric and contributes zero to f81_violation.
     Verified all 9 = 3×3 pure-Pauli combinations (X⊗X, X⊗Y, ..., Z⊗Z) at N=3.

  2. **σ-channels (σ⁻⊗σ⁻, σ⁺⊗σ⁺, σ⁻⊗σ⁺, σ⁺⊗σ⁻) break F81 with closed form.**
     For a single bond (i, j) at chain N ≥ 3:

         ‖D_2qubit_σ_odd‖²_F  =  γ²  ·  6  ·  4^(N-3)

     Verified bit-exact at N = 3 (v² = 0.06 for γ = 0.1), N = 4 (v² = 0.24),
     N = 5 (v² = 0.96). Linear in γ; v/γ = √(6 · 4^(N-3)) = √6 · 2^(N-3).

  3. **Detailed-balance pairs (per bond):**
     (σ⁻⊗σ⁻, σ⁺⊗σ⁺)   correlated cooling/heating  → cancel exactly
     (σ⁻⊗σ⁺, σ⁺⊗σ⁻)   exchange directionality      → cancel exactly
     The two pairs are orthogonal (Pythagoras): cooling + exchange gives
     v² = 2 · single-channel v², i.e. independent contributions.

  4. **Non-overlapping bonds: additive.**
     For two bonds with no shared qubit (e.g., (0,1) and (2,3) at N=4),
     v²_total  =  v²_bond1  +  v²_bond2  exactly.

  5. **Overlapping bonds (shared qubit): non-additive, constructive.**
     For two bonds sharing one qubit (e.g., (0,1) and (1,2) at N=4),
     v²_total exceeds the independent-additive prediction. At γ=0.1, N=4:
     non-overlapping (0,1)+(2,3): v² = 0.48 = 2 · 0.24
     overlapping (0,1)+(1,2):     v² = 0.64 > 0.48
     Excess 0.16 reflects constructive cross-term from shared qubit's
     Π²-anti-symmetric matrix elements.

This is the partial F86 picture. The single-bond closed form and
non-overlapping bond additivity are clean. The overlapping-bond
cross-term structure is open (would need to derive how shared-qubit
σ⁻ contributions overlap in Pauli basis).

Promote to F86 framework theorem when overlap structure is closed.
For now, this script + ANALYTICAL_FORMULAS draft are the documentation.
"""
from __future__ import annotations

import sys
from pathlib import Path
from itertools import product

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from framework.lindblad import lindbladian_general, palindrome_residual
from framework.pauli import _vec_to_pauli_basis_transform, _build_bilinear, ur_pauli
from framework.symmetry import build_pi_full


SIGMA_MINUS = np.array([[0, 1], [0, 0]], dtype=complex)
SIGMA_PLUS = np.array([[0, 0], [1, 0]], dtype=complex)
ID2 = np.eye(2, dtype=complex)


def two_qubit_op(N, i, j, op_a, op_b):
    """Build a 2-qubit operator op_a on site i, op_b on site j, identity elsewhere."""
    ops = [ID2] * N
    ops[i] = op_a
    ops[j] = op_b
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def f81_violation_2qubit(N, c_op_list_with_gamma):
    """Compute f81_violation for soft H = J(XY+YX) plus given c_ops with rates."""
    bonds = [(i, i + 1) for i in range(N - 1)]
    H = _build_bilinear(N, bonds, [('X', 'Y', 1.0), ('Y', 'X', 1.0)])
    c_ops = [np.sqrt(g) * c for (c, g) in c_op_list_with_gamma]
    L = lindbladian_general(H, c_ops)
    M = palindrome_residual(L, 0.0, N)
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)
    PiMPi = Pi @ M @ Pi_inv
    M_anti = (M - PiMPi) / 2
    T = _vec_to_pauli_basis_transform(N)
    d = 2 ** N
    L_vec = -1j * (np.kron(H, np.eye(d, dtype=complex)) -
                   np.kron(np.eye(d, dtype=complex), H.T))
    L_H = (T.conj().T @ L_vec @ T) / d
    return float(np.linalg.norm(M_anti - L_H))


def main():
    print('=' * 72)
    print('2-qubit Dissipator Exploration (J / partial F86)')
    print('=' * 72)
    print()

    # 1. Pauli-Channel Cancellation Lemma extension
    print('1. Pauli-Channel Cancellation extends to 2-qubit Pauli channels')
    print()
    N = 3
    print(f'   N = {N}, bond (0,1):')
    for la, lb in product('XYZ', repeat=2):
        op = two_qubit_op(N, 0, 1, ur_pauli(la), ur_pauli(lb))
        v = f81_violation_2qubit(N, [(op, 0.1)])
        mark = '✓ Π²-sym' if v < 1e-9 else f'✗ violation = {v:.4e}'
        print(f'     D[{la}⊗{lb}]: {mark}')
    print()

    # 2. Single-bond σ-channel closed form
    print('2. Single-bond σ-channel f81_violation closed form')
    print(f'   Predicted: v² = γ² · 6 · 4^(N-3)')
    print()
    print(f'   {"N":>3} {"channel":<24} {"v":>10} {"v²":>10} {"predicted v²":>14}')
    print('   ' + '-' * 70)
    for N in [3, 4, 5]:
        for label, op_a, op_b in [
            ('σ⁻⊗σ⁻ (cooling)', SIGMA_MINUS, SIGMA_MINUS),
            ('σ⁺⊗σ⁺ (heating)', SIGMA_PLUS, SIGMA_PLUS),
            ('σ⁻⊗σ⁺ (exchange)', SIGMA_MINUS, SIGMA_PLUS),
        ]:
            op = two_qubit_op(N, 0, 1, op_a, op_b)
            v = f81_violation_2qubit(N, [(op, 0.1)])
            pred = 0.01 * 6 * 4 ** (N - 3)
            print(f'   {N:>3} {label:<24} {v:>10.4f} {v*v:>10.6f} {pred:>14.6f}')
    print()

    # 3. Detailed-balance pairs
    print('3. Detailed balance: cooling + heating cancels per bond')
    print()
    N = 3
    op_dd = two_qubit_op(N, 0, 1, SIGMA_MINUS, SIGMA_MINUS)
    op_uu = two_qubit_op(N, 0, 1, SIGMA_PLUS, SIGMA_PLUS)
    op_du = two_qubit_op(N, 0, 1, SIGMA_MINUS, SIGMA_PLUS)
    op_ud = two_qubit_op(N, 0, 1, SIGMA_PLUS, SIGMA_MINUS)
    cases = [
        ('σ⁻⊗σ⁻ + σ⁺⊗σ⁺ (γ=0.1 each)', [(op_dd, 0.1), (op_uu, 0.1)]),
        ('σ⁻⊗σ⁺ + σ⁺⊗σ⁻ (γ=0.1 each)', [(op_du, 0.1), (op_ud, 0.1)]),
        ('all 4 σ-channels (γ=0.1 each)', [(op_dd, 0.1), (op_uu, 0.1), (op_du, 0.1), (op_ud, 0.1)]),
    ]
    for label, c_ops in cases:
        v = f81_violation_2qubit(N, c_ops)
        mark = '✓ balanced' if abs(v) < 1e-10 else f'violation = {v:.4f}'
        print(f'   {label:<40} {mark}')
    print()

    # 4. Non-overlapping vs overlapping bonds
    print('4. Multi-bond: non-overlapping additive, overlapping non-additive')
    print()
    N = 4
    op01 = two_qubit_op(N, 0, 1, SIGMA_MINUS, SIGMA_MINUS)
    op12 = two_qubit_op(N, 1, 2, SIGMA_MINUS, SIGMA_MINUS)
    op23 = two_qubit_op(N, 2, 3, SIGMA_MINUS, SIGMA_MINUS)
    v_01_alone = f81_violation_2qubit(N, [(op01, 0.1)])
    v_disjoint = f81_violation_2qubit(N, [(op01, 0.1), (op23, 0.1)])
    v_overlap = f81_violation_2qubit(N, [(op01, 0.1), (op12, 0.1)])
    print(f'   N=4, γ=0.1, σ⁻⊗σ⁻ on each bond:')
    print(f'     (0,1) only:                v={v_01_alone:.4f}, v²={v_01_alone**2:.4f}')
    print(f'     (0,1) + (2,3) non-overlap: v={v_disjoint:.4f}, v²={v_disjoint**2:.4f}')
    print(f'     additivity check: 2·v_single² = {2*v_01_alone**2:.4f}  (matches?)')
    print(f'     (0,1) + (1,2) overlap:     v={v_overlap:.4f}, v²={v_overlap**2:.4f}')
    print(f'     excess over additive:       {v_overlap**2 - 2*v_01_alone**2:.4f}')
    print()

    print('=' * 72)
    print('Summary')
    print('=' * 72)
    print("""
Pauli-Channel Cancellation generalizes verbatim to 2-qubit Pauli
channels D[c ⊗ c'] for c, c' ∈ {I, X, Y, Z}: all give zero F81 violation.

The 2-qubit σ-channels (σ⁻⊗σ⁻, σ⁺⊗σ⁺, σ⁻⊗σ⁺, σ⁺⊗σ⁻) break F81 with
single-bond closed form ‖D_odd‖² = 6·γ²·4^(N-3). They form two
detailed-balance pairs that cancel at thermal equilibrium per bond.

Non-overlapping bonds are additive; overlapping bonds (shared qubit)
have constructive cross-terms not captured by the per-bond formula.
The overlap structure is open and would need separate analysis
(e.g., via Pauli-decomposition of σ⁻⊗σ⁻ on adjacent bonds).

For hardware diagnostics: this would extend F84's f81_violation reading
to also detect *correlated* amplitude-damping channels (vacuum-fluctuation
contribution from pair emission/absorption, distinct from single-site σ±).
Without hardware data showing correlated decay, the F86 promotion is
deferred until the overlap structure is closed and an operational use
case appears.
""")


if __name__ == "__main__":
    main()
