"""Symmetry operators and panels: Π conjugation, K-chiral, F71 chain-mirror, Y-parity.

Public API:
  Π conjugation:  pi_action, pi_squared_eigenvalue, build_pi_full
  Parity selection: respects_bit_a_parity, respects_bit_b_parity, is_both_parity_even
  Chain-mirror (F71): chain_mirror_state, f71_symmetric_projector,
                       f71_antisymmetric_projector, f71_eigenstate_class,
                       receiver_engineering_signature, bond_mirror_basis
  Chiral panel: chiral_K_full, k_classify_pauli, k_classify_hamiltonian, chiral_panel
  Y-parity panel: y_parity_panel
"""
from __future__ import annotations

import numpy as np

from .pauli import (
    PAULI_LABELS, _PAULI_MATRICES,
    _resolve, bit_a, bit_b,
    _k_to_indices, _indices_to_k,
    _vec_to_pauli_basis_transform, pauli_basis_vector,
)


# ----------------------------------------------------------------------
# Π conjugation (per site and full Pauli-string basis)
# ----------------------------------------------------------------------
# Π per site: I↔X (sign 1), Y↔Z (sign i). bit_a flips, bit_b preserved,
# phase = i if bit_b=1.

def pi_action(idx):
    """Return ((new_a, new_b), phase) for Π acting on Pauli at idx."""
    a, b = _resolve(idx)
    return (1 - a, b), (1j if b == 1 else 1)


def pi_squared_eigenvalue(indices):
    """Π² eigenvalue on a Pauli string = (-1)^{total_bit_b}."""
    return (-1) ** (sum(bit_b(idx) for idx in indices) % 2)


def build_pi_full(N):
    """Π in the 4^N Pauli-string basis: 4^N × 4^N matrix."""
    d2 = 4 ** N
    Pi = np.zeros((d2, d2), dtype=complex)
    for k in range(d2):
        indices = _k_to_indices(k, N)
        new_indices = []
        sign = 1
        for idx in indices:
            (na, nb), phase = pi_action(idx)
            new_indices.append((na, nb))
            sign *= phase
        new_k = _indices_to_k(new_indices)
        Pi[new_k, k] = sign
    return Pi


# ----------------------------------------------------------------------
# Parity selection on bilinears (bit_a / bit_b parities)
# ----------------------------------------------------------------------

def respects_bit_a_parity(idx_pair):
    return (bit_a(idx_pair[0]) + bit_a(idx_pair[1])) % 2 == 0


def respects_bit_b_parity(idx_pair):
    return (bit_b(idx_pair[0]) + bit_b(idx_pair[1])) % 2 == 0


def is_both_parity_even(idx_pair):
    return respects_bit_a_parity(idx_pair) and respects_bit_b_parity(idx_pair)


# ----------------------------------------------------------------------
# F71: chain-mirror spatial symmetry (Section 17)
# ----------------------------------------------------------------------
# F71 is the chain-reflection symmetry of uniform-J Liouvillian:
# R: |b₀b₁…b_{N-1}⟩ → |b_{N-1}…b₁b₀⟩ commutes with H, D, hence L. Proved
# in docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md.
#
# Bond-input parity: R induces a permutation J_b ↔ J_{N-2-b} on the (N-1)
# bonds. Symmetric/antisymmetric subspaces have sizes:
#   odd N (even bond count):  balanced k+k split.
#   even N (odd bond count):  unbalanced (k+1)+k split (one self-mirror bond).
#
# Empirical (EQ-024 2026-04-28): F71-eigenstate receivers are
# capacity-optimal at odd N, capacity-suboptimal at even N. See
# experiments/J_BLIND_RECEIVER_CLASSES.md.

def chain_mirror_state(N):
    """Chain-mirror operator R on 2^N Hilbert space.

    R|b₀b₁…b_{N-1}⟩ = |b_{N-1}…b₁b₀⟩ (bit-reversal permutation).
    Real, symmetric, involutory (R² = I).
    """
    dim = 2 ** N
    R = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        rev = 0
        x = i
        for _ in range(N):
            rev = (rev << 1) | (x & 1)
            x >>= 1
        R[rev, i] = 1.0
    return R


def f71_symmetric_projector(N):
    """P_sym = (I + R)/2 onto the F71-symmetric (+1) eigenspace."""
    return (np.eye(2 ** N, dtype=complex) + chain_mirror_state(N)) / 2.0


def f71_antisymmetric_projector(N):
    """P_anti = (I − R)/2 onto the F71-antisymmetric (−1) eigenspace."""
    return (np.eye(2 ** N, dtype=complex) - chain_mirror_state(N)) / 2.0


def f71_eigenstate_class(psi, tol=1e-6):
    """Classify state psi as F71-eigenstate.

    Returns +1, -1, or None (F71-mixed/breaking).
    F71-eigenstates of either sign trigger the bond-block decomposition
    of the J-Jacobian.
    """
    N = int(round(np.log2(len(psi))))
    if 2 ** N != len(psi):
        raise ValueError(f"psi length {len(psi)} is not a power of 2")
    R = chain_mirror_state(N)
    Rpsi = R @ psi
    norm_psi = np.linalg.norm(psi)
    if norm_psi < tol:
        raise ValueError("psi has zero norm")
    overlap = np.vdot(psi, Rpsi) / (norm_psi ** 2)
    if abs(overlap.imag) > tol:
        return None
    o_real = float(overlap.real)
    if abs(o_real - 1.0) < tol:
        return +1
    if abs(o_real + 1.0) < tol:
        return -1
    return None


def bond_mirror_basis(N):
    """Symmetric/antisymmetric basis vectors for the bond-input space.

    R̄: J_b → J_{N-2-b} acts on R^{N-1}. Returns (sym_basis, asym_basis)
    as orthonormal arrays.

    Dimensional structure (verified at N=3..8):
      odd N: balanced k+k split, n_sym = n_asym = (N-1)/2
      even N: unbalanced (k+1)+k split, self-mirror bond at (N-2)/2
    """
    n_bonds = N - 1
    sym, asym = [], []
    used = [False] * n_bonds
    for b in range(n_bonds):
        if used[b]:
            continue
        b_mirror = n_bonds - 1 - b
        if b == b_mirror:
            v = np.zeros(n_bonds)
            v[b] = 1.0
            sym.append(v)
            used[b] = True
        else:
            v_s = np.zeros(n_bonds); v_s[b] = 1 / np.sqrt(2); v_s[b_mirror] = 1 / np.sqrt(2)
            v_a = np.zeros(n_bonds); v_a[b] = 1 / np.sqrt(2); v_a[b_mirror] = -1 / np.sqrt(2)
            sym.append(v_s)
            asym.append(v_a)
            used[b] = True
            used[b_mirror] = True
    return np.array(sym), np.array(asym)


def receiver_engineering_signature(psi):
    """F71-based receiver-engineering favorability forecast.

    Returns dict with N, f71_eigenvalue (+1/-1/None), bond_block_dims,
    bond_block_balanced, prediction (str). Combines f71_eigenstate_class
    with bond_mirror_basis.

    Tier 2: verified at N=5 (balanced, F71-eigenstate optimal) and N=6
    (unbalanced, F71-breaking wins). See J_BLIND_RECEIVER_CLASSES.md.
    """
    N = int(round(np.log2(len(psi))))
    eig = f71_eigenstate_class(psi)
    sym, asym = bond_mirror_basis(N)
    n_sym, n_asym = len(sym), len(asym)
    balanced = (n_sym == n_asym)

    if eig is None:
        prediction = 'no-prediction (F71-mixed state — no block structure)'
    elif balanced:
        prediction = 'capacity-optimal (balanced split, F71-eigenstate)'
    else:
        prediction = 'capacity-suboptimal (unbalanced split, F71-eigenstate)'

    return {
        'N': N,
        'f71_eigenvalue': eig,
        'bond_block_dims': (n_sym, n_asym),
        'bond_block_balanced': balanced,
        'prediction': prediction,
    }


# ----------------------------------------------------------------------
# Chiral (K-sublattice) panel
# ----------------------------------------------------------------------
# K_full = ⊗_{odd i} Z_i. K anti-commutes with H_xy = (X_a X_{a+1} +
# Y_a Y_{a+1})/2: K H K = -H. Z-dephasing commutes with K trivially.
# When K anti-commutes with H, K is a super-operator symmetry of L.

def chiral_K_full(N):
    """K_full = ⊗_{odd i} Z_i — chain chiral / sublattice operator."""
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    out = np.array([[1.0]], dtype=complex)
    for i in range(N):
        out = np.kron(out, Z if (i % 2 == 1) else I2)
    return out


def k_classify_pauli(N, K_full=None):
    """For each Pauli string α, return its K-eigenvalue (+1, -1, or 0=mixed)."""
    if K_full is None:
        K_full = chiral_K_full(N)
    out = {}
    for alpha in range(4 ** N):
        idx_tuple = _k_to_indices(alpha, N)
        P = np.array([[1.0]], dtype=complex)
        for idx in idx_tuple:
            P = np.kron(P, _PAULI_MATRICES[idx])
        KPK = K_full @ P @ K_full
        if np.allclose(KPK, P, atol=1e-10):
            out[alpha] = +1
        elif np.allclose(KPK, -P, atol=1e-10):
            out[alpha] = -1
        else:
            out[alpha] = 0
    return out


def k_classify_hamiltonian(H, N, K_full=None):
    """Classify H as 'K-even', 'K-odd', or 'K-mixed'."""
    if K_full is None:
        K_full = chiral_K_full(N)
    KHK = K_full @ H @ K_full
    if np.allclose(KHK, H, atol=1e-10):
        return 'K-even'
    if np.allclose(KHK, -H, atol=1e-10):
        return 'K-odd'
    return 'K-mixed'


def chiral_panel(H, rho_0, N, K_full=None):
    """Chiral structure panel: K-symmetry status of L + ρ_0's K-decomposition.

    Returns dict with K_status, K_symmetric_L, rho_decomposition (w_plus,
    w_minus, rho_plus, rho_minus), rho_is_K_eigenstate, pauli_classification
    (n_K_even, n_K_odd), static_zero_at_t0, chiral_protected_dynamic.

    Dynamical K-protection requires: K is a symmetry of L AND ρ_0 is a
    K-eigenstate. Both are reported.
    """
    if K_full is None:
        K_full = chiral_K_full(N)

    K_status = k_classify_hamiltonian(H, N, K_full=K_full)
    K_symmetric_L = K_status in ('K-even', 'K-odd')

    rho_KKt = K_full @ rho_0 @ K_full
    rho_plus = (rho_0 + rho_KKt) / 2
    rho_minus = (rho_0 - rho_KKt) / 2
    w_plus = float(np.real(np.trace(rho_plus @ rho_plus)))
    w_minus = float(np.real(np.trace(rho_minus @ rho_minus)))
    rho_is_K_eigenstate = (w_plus < 1e-12) or (w_minus < 1e-12)

    k_class = k_classify_pauli(N, K_full=K_full)
    n_even = sum(1 for v in k_class.values() if v == +1)
    n_odd = sum(1 for v in k_class.values() if v == -1)

    static_zero = []
    for alpha in range(1, 4 ** N):
        P_alpha_idx = _k_to_indices(alpha, N)
        P_alpha = np.array([[1.0]], dtype=complex)
        for idx in P_alpha_idx:
            P_alpha = np.kron(P_alpha, _PAULI_MATRICES[idx])
        if k_class[alpha] == +1:
            exp_val = float(np.real(np.trace(P_alpha @ rho_plus)))
        elif k_class[alpha] == -1:
            exp_val = float(np.real(np.trace(P_alpha @ rho_minus)))
        else:
            exp_val = 1.0
        if abs(exp_val) < 1e-10:
            static_zero.append(alpha)

    chiral_protected_dynamic = []
    if K_symmetric_L and rho_is_K_eigenstate:
        for alpha in range(1, 4 ** N):
            if (k_class[alpha] == +1 and w_minus < 1e-12) or \
               (k_class[alpha] == -1 and w_plus < 1e-12):
                chiral_protected_dynamic.append(alpha)

    static_zero_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(a, N))
        for a in static_zero
    ]
    chiral_dynamic_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(a, N))
        for a in chiral_protected_dynamic
    ]

    return {
        'K_status': K_status,
        'K_symmetric_L': K_symmetric_L,
        'rho_decomposition': {
            'w_plus': w_plus, 'w_minus': w_minus,
            'rho_plus': rho_plus, 'rho_minus': rho_minus,
        },
        'rho_is_K_eigenstate': rho_is_K_eigenstate,
        'pauli_classification': {'n_K_even': n_even, 'n_K_odd': n_odd},
        'static_zero_at_t0': static_zero_labels,
        'chiral_protected_dynamic': chiral_dynamic_labels,
    }


# ----------------------------------------------------------------------
# Y-parity panel: Z₂ grading on Pauli-Y count
# ----------------------------------------------------------------------

def _y_parity_classify_paulis(N):
    """Return Y-parity (0 or 1) for each Pauli index α = 0..4^N-1."""
    out = np.zeros(4 ** N, dtype=int)
    for alpha in range(4 ** N):
        idx_tuple = _k_to_indices(alpha, N)
        n_y = sum(1 for idx in idx_tuple if idx == (1, 1))
        out[alpha] = n_y % 2
    return out


def y_parity_panel(H, gamma_l, rho_0, N, gamma_t1_l=None):
    """Y-parity panel: tests whether L preserves Y-parity, reports protected Y-odd Paulis.

    Returns dict with L_preserves_Y_parity, L_offdiag_weight,
    rho_0_Y_decomposition, pauli_classification, Y_odd_observables,
    Y_parity_protected.

    Y-parity-protected observables stay 0 forever iff L preserves Y-parity
    AND ρ_0 has no Y-odd content (or vice versa).
    """
    from .lindblad import lindbladian_z_dephasing, lindbladian_z_plus_t1
    if gamma_t1_l is None:
        gamma_t1_l = [0.0] * N

    if any(g != 0 for g in gamma_t1_l):
        L = lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = lindbladian_z_dephasing(H, gamma_l)

    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

    y_parity = _y_parity_classify_paulis(N)
    even_idx = np.where(y_parity == 0)[0]
    odd_idx = np.where(y_parity == 1)[0]
    n_even = len(even_idx)
    n_odd = len(odd_idx)

    L_oe = L_pauli[np.ix_(odd_idx, even_idx)]
    L_eo = L_pauli[np.ix_(even_idx, odd_idx)]
    offdiag_norm = float(np.linalg.norm(L_oe)) + float(np.linalg.norm(L_eo))
    diag_norm = float(np.linalg.norm(L_pauli[np.ix_(even_idx, even_idx)])) \
                 + float(np.linalg.norm(L_pauli[np.ix_(odd_idx, odd_idx)]))
    relative_offdiag = offdiag_norm / (diag_norm + 1e-15)
    L_preserves_Y_parity = (relative_offdiag < 1e-10)

    rho_pauli_vec = pauli_basis_vector(rho_0, N)
    rho_even_pauli = rho_pauli_vec.copy()
    rho_even_pauli[odd_idx] = 0
    rho_odd_pauli = rho_pauli_vec.copy()
    rho_odd_pauli[even_idx] = 0
    w_even = float(np.linalg.norm(rho_even_pauli))
    w_odd = float(np.linalg.norm(rho_odd_pauli))

    y_odd_labels = [
        ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
        for a in odd_idx
    ]

    Y_parity_protected = []
    if L_preserves_Y_parity and w_odd < 1e-10:
        Y_parity_protected = list(y_odd_labels)
    elif L_preserves_Y_parity and w_even < 1e-10:
        even_labels = [
            ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(int(a), N))
            for a in even_idx if a != 0
        ]
        Y_parity_protected = even_labels

    return {
        'L_preserves_Y_parity': L_preserves_Y_parity,
        'L_offdiag_weight': offdiag_norm,
        'L_relative_offdiag': relative_offdiag,
        'rho_0_Y_decomposition': {'w_even': w_even, 'w_odd': w_odd},
        'pauli_classification': {'n_Y_even': int(n_even), 'n_Y_odd': int(n_odd)},
        'Y_odd_observables': y_odd_labels,
        'Y_parity_protected': Y_parity_protected,
    }
