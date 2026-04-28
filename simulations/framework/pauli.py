"""Pauli matrices, multi-qubit strings, basis transforms, bilinear builder.

Contains the bit_a / bit_b indexed Pauli primitives plus the basis-transform
infrastructure used by every other framework module.

Public API:
  ur_pauli, pauli_string, site_op, _build_bilinear
  pauli_matrix(letter)         (simple letter-indexed matrix)
  PAULI_LABELS, LABEL_TO_INDEX
  bit_a, bit_b, total_bit_a, total_bit_b_parity
  _k_to_indices, _indices_to_k
  _vec_to_pauli_basis_transform, pauli_basis_vector, _pauli_label
"""
from __future__ import annotations

from itertools import product as iproduct

import numpy as np

# ----------------------------------------------------------------------
# (bit_a, bit_b) indexed Pauli matrices
# ----------------------------------------------------------------------
#  (0,0) = I    (1,0) = X    (0,1) = Z    (1,1) = Y

_PAULI_MATRICES = {
    (0, 0): np.eye(2, dtype=complex),
    (1, 0): np.array([[0, 1], [1, 0]], dtype=complex),
    (0, 1): np.array([[1, 0], [0, -1]], dtype=complex),
    (1, 1): np.array([[0, -1j], [1j, 0]], dtype=complex),
}

PAULI_LABELS = {(0, 0): 'I', (1, 0): 'X', (0, 1): 'Z', (1, 1): 'Y'}
LABEL_TO_INDEX = {v: k for k, v in PAULI_LABELS.items()}


def ur_pauli(a, b=None):
    """Return the 2×2 Pauli at index (a, b) or by label string ('I','X','Y','Z')."""
    if b is None and isinstance(a, str):
        return _PAULI_MATRICES[LABEL_TO_INDEX[a]]
    return _PAULI_MATRICES[(a, b)]


def pauli_matrix(letter):
    """Letter-indexed Pauli matrix. Equivalent to ur_pauli(letter)."""
    return ur_pauli(letter)


def bit_a(idx):
    """bit_a (n_XY) of a Pauli (idx is (a,b) tuple or label)."""
    if isinstance(idx, str):
        idx = LABEL_TO_INDEX[idx]
    return idx[0]


def bit_b(idx):
    """bit_b (n_YZ) of a Pauli (idx is (a,b) tuple or label)."""
    if isinstance(idx, str):
        idx = LABEL_TO_INDEX[idx]
    return idx[1]


def total_bit_a(indices):
    """XY-weight (sum of bit_a) of a Pauli string."""
    return sum(bit_a(idx) for idx in indices)


def total_bit_b_parity(indices):
    """Π²-parity (sum of bit_b mod 2) of a Pauli string."""
    return sum(bit_b(idx) for idx in indices) % 2


# ----------------------------------------------------------------------
# Multi-qubit Pauli strings
# ----------------------------------------------------------------------

def _resolve(idx):
    """Convert label string or (a,b) tuple into (a,b) tuple."""
    if isinstance(idx, str):
        return LABEL_TO_INDEX[idx]
    return idx


def pauli_string(indices_or_labels):
    """Build a 2^N × 2^N Pauli string from a list of N (a,b) tuples or labels."""
    out = ur_pauli(*_resolve(indices_or_labels[0]))
    for idx in indices_or_labels[1:]:
        out = np.kron(out, ur_pauli(*_resolve(idx)))
    return out


def site_op(N, site, idx_or_label):
    """N-qubit operator with single-Pauli idx on the given site, I elsewhere."""
    ops = [ur_pauli('I')] * N
    ops[site] = ur_pauli(*_resolve(idx_or_label))
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


# ----------------------------------------------------------------------
# Index encoding for the 4^N Pauli-string basis
# ----------------------------------------------------------------------

def _k_to_indices(k, N):
    """Convert flat index k ∈ [0, 4^N) to N-tuple of (a,b) tuples."""
    out = []
    kk = k
    for _ in range(N):
        i = kk % 4
        kk //= 4
        out.append((i & 1, (i >> 1) & 1))
    return tuple(reversed(out))


def _indices_to_k(indices):
    """Inverse of _k_to_indices."""
    k = 0
    for (a, b) in indices:
        i = a + 2 * b
        k = k * 4 + i
    return k


def _pauli_label(k, N):
    """Convert flat Pauli index k to label string like 'XIZ' (left-most-first)."""
    return ''.join(PAULI_LABELS[idx] for idx in _k_to_indices(k, N))


# ----------------------------------------------------------------------
# Basis transforms (Pauli ↔ vec)
# ----------------------------------------------------------------------

def _vec_to_pauli_basis_transform(N):
    """Transformation matrix from Pauli-string basis to column-stack vec.

    M (d² × 4^N) such that vec(σ_α) = M[:, α]. M†M = 2^N · I (orthogonality).
    """
    d2 = 4 ** N
    d = 2 ** N
    M = np.zeros((d * d, d2), dtype=complex)
    for k in range(d2):
        indices = _k_to_indices(k, N)
        sigma = pauli_string(list(indices))
        M[:, k] = sigma.flatten('F')
    return M


def pauli_basis_vector(rho, N):
    """Express ρ as a 4^N Pauli-basis coefficient vector.

    vec[k] = (1/2^N) · Tr(σ_k · ρ), so that ρ = (1/2^N) · Σ_k vec[k] · σ_k.
    """
    d2 = 4 ** N
    vec = np.zeros(d2, dtype=complex)
    for k in range(d2):
        sigma_k = pauli_string(list(_k_to_indices(k, N)))
        vec[k] = np.trace(sigma_k @ rho) / (2 ** N)
    return vec


# ----------------------------------------------------------------------
# Hamiltonian builder from Pauli-pair bilinears
# ----------------------------------------------------------------------

def _build_bilinear(N, bonds, terms):
    """Build H = Σ_bond Σ_term coeff · σ_a^l σ_b^{l+1}.

    Args:
        N: number of qubits
        bonds: list of (i, j) bond tuples
        terms: list of (label1, label2, coeff) triples

    Returns:
        2^N × 2^N Hermitian matrix.
    """
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for (la, lb, coeff) in terms:
            H = H + coeff * site_op(N, i, la) @ site_op(N, j, lb)
    return H


def _site_op_kron(op, site, N):
    """Helper: 2×2 op on `site`, identity elsewhere (different signature from site_op)."""
    I2 = np.eye(2, dtype=complex)
    factors = [I2] * N
    factors[site] = op
    out = factors[0]
    for f in factors[1:]:
        out = np.kron(out, f)
    return out
