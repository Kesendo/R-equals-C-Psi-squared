#!/usr/bin/env python3
"""F50 max-spin closed-form verification: Dicke endpoint ladder rungs.

Canonical statement and proof: docs/proofs/PROOF_WEIGHT1_DEGENERACY.md
§Spin-isotypic decomposition / max-spin closed-form.

Checks that the conjectured Dicke endpoint basis spans the empirically-extracted
pure-weight-w sym-supported subspace for N = 2..5 across all w. Two routes:
(1) brute SVD nullspace of (I - P_sym⊗P_sym) on weight-w Pauli ops; (2) the
proposed closed-form basis {|D_0⟩⟨D_w| ± h.c., |D_{N-w}⟩⟨D_N| ± h.c.}.
"""
from __future__ import annotations

import sys
import numpy as np
from itertools import product, combinations

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string_op(N, letters):
    op = np.array([[1]], complex)
    for L in letters:
        op = np.kron(op, PAULI[L])
    return op


def weight(letters):
    return sum(1 for c in letters if c in "XY")


def dicke_state(N, k):
    d = 2 ** N
    vec = np.zeros(d, complex)
    for positions in combinations(range(N), k):
        idx = 0
        for p in positions:
            idx |= 1 << (N - 1 - p)
        vec[idx] += 1.0
    return vec / np.linalg.norm(vec)


def sym_projector(N):
    d = 2 ** N
    P = np.zeros((d, d), complex)
    for k in range(N + 1):
        D_k = dicke_state(N, k)
        P += np.outer(D_k, D_k.conj())
    return P


def empirical_pure_weight_basis(N, target_w, P_sym, tol=1e-9):
    """Compute empirical pure-weight-w ops supported in sym (max-spin block)."""
    d = 2 ** N
    all_letters = list(product("IXYZ", repeat=N))
    weight_ops = [pauli_string_op(N, L) for L in all_letters if weight(L) == target_w]
    M_w = np.column_stack([op.flatten("F") for op in weight_ops])
    n_w = M_w.shape[1]
    T = np.zeros((d * d, n_w), complex)
    for k in range(n_w):
        A = M_w[:, k].reshape(d, d, order="F")
        T[:, k] = (A - P_sym @ A @ P_sym).flatten("F")
    U, S, Vh = np.linalg.svd(T)
    null_dim = int(np.sum(S < tol))
    null_basis = Vh[-null_dim:, :].conj().T if null_dim > 0 else np.zeros((n_w, 0), complex)
    ops = []
    for j in range(null_dim):
        c = null_basis[:, j]
        A = (M_w @ c).reshape(d, d, order="F")
        ops.append(A)
    return ops


def closed_form_basis(N, w):
    """Closed-form: Dicke endpoint ladder rungs."""
    d = 2 ** N
    D0 = dicke_state(N, 0)
    DN = dicke_state(N, N)
    ops = []
    if w == 0:
        ops.append(np.outer(D0, D0.conj()))
        ops.append(np.outer(DN, DN.conj()))
    elif w == N:
        ops.append(np.outer(D0, DN.conj()) + np.outer(DN, D0.conj()))
        ops.append(1j * (np.outer(D0, DN.conj()) - np.outer(DN, D0.conj())))
    else:
        Dw = dicke_state(N, w)
        DNmw = dicke_state(N, N - w)
        ops.append(np.outer(D0, Dw.conj()) + np.outer(Dw, D0.conj()))
        ops.append(1j * (np.outer(D0, Dw.conj()) - np.outer(Dw, D0.conj())))
        ops.append(np.outer(DNmw, DN.conj()) + np.outer(DN, DNmw.conj()))
        ops.append(1j * (np.outer(DNmw, DN.conj()) - np.outer(DN, DNmw.conj())))
    return ops


def verify_endpoint_projector_formulas(N):
    """Verify |D_0⟩⟨D_0| = (1/2^N) Σ_k e_k(Z), |D_N⟩⟨D_N| = (1/2^N) Σ_k (-1)^k e_k(Z)."""
    d = 2 ** N
    weight0_letters = list(product("IZ", repeat=N))
    sum_ek = np.zeros((d, d), complex)
    sum_alt_ek = np.zeros((d, d), complex)
    for L in weight0_letters:
        z_count = sum(1 for c in L if c == "Z")
        op = pauli_string_op(N, list(L))
        sum_ek += op
        sum_alt_ek += (-1) ** z_count * op
    D0 = dicke_state(N, 0)
    DN = dicke_state(N, N)
    diff_1 = np.linalg.norm(sum_ek / (2 ** N) - np.outer(D0, D0.conj()))
    diff_2 = np.linalg.norm(sum_alt_ek / (2 ** N) - np.outer(DN, DN.conj()))
    return diff_1, diff_2


def main():
    print("=" * 78)
    print("F50 max-spin closed-form verification (Dicke endpoint ladder rungs)")
    print("=" * 78)
    print()

    print("(1) Endpoint projector formulas |D_0⟩⟨D_0|, |D_N⟩⟨D_N|:")
    for N in [2, 3, 4, 5]:
        d1, d2 = verify_endpoint_projector_formulas(N)
        marker1 = "✓" if d1 < 1e-10 else "✗"
        marker2 = "✓" if d2 < 1e-10 else "✗"
        print(f"  N={N}: ||Π(I+Z)/2^N − |D_0⟩⟨D_0|||  = {d1:.2e} {marker1}")
        print(f"        ||Π(I-Z)/2^N − |D_N⟩⟨D_N|||  = {d2:.2e} {marker2}")
    print()

    print("(2) Full max-spin pure-weight basis = endpoint ladder rungs:")
    print(f"  {'N':>3} {'w':>3} {'empirical':>9} {'closed-form':>11} {'spans match?':>13}")
    print("  " + "-" * 50)
    all_pass = True
    for N in [2, 3, 4, 5]:
        P_sym = sym_projector(N)
        for w in range(N + 1):
            emp_ops = empirical_pure_weight_basis(N, w, P_sym)
            cf_ops = closed_form_basis(N, w)
            M_emp = np.column_stack([op.flatten("F") for op in emp_ops]) if emp_ops else np.zeros((4 ** N, 0), complex)
            M_cf = np.column_stack([op.flatten("F") for op in cf_ops]) if cf_ops else np.zeros((4 ** N, 0), complex)
            r_emp = np.linalg.matrix_rank(M_emp, tol=1e-8) if M_emp.size else 0
            r_cf = np.linalg.matrix_rank(M_cf, tol=1e-8) if M_cf.size else 0
            combined = np.column_stack([M_emp, M_cf]) if M_emp.size and M_cf.size else M_emp
            r_comb = np.linalg.matrix_rank(combined, tol=1e-8) if combined.size else 0
            span_match = (r_comb == r_emp == r_cf)
            marker = "✓" if span_match else "✗"
            if not span_match:
                all_pass = False
            print(f"  {N:>3} {w:>3} {r_emp:>9} {r_cf:>11} {marker:>13}")
    print()
    print(f"  ALL CHECKS PASS: {all_pass}")
    print()

    print("(3) Summary: max-spin block of M(N+1) decomposes as")
    print("  - pure-weight ops (4N total): 2 at edges (w=0, N), 4 at each interior weight")
    print("  - multi-weight ops ((N-1)² total): middle Dicke transitions |D_k⟩⟨D_l|, k,l ∈ {1,...,N-1}")
    print("  - Total: (N+1)² operators in M(N+1)")


if __name__ == "__main__":
    main()
