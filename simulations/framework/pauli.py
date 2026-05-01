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
  bonding_mode_state(N, k), bonding_mode_pair_state(N, k)
  polarity_state(N, signs)
"""
from __future__ import annotations

from functools import lru_cache
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

@lru_cache(maxsize=None)
def _vec_to_pauli_basis_transform(N):
    """Transformation matrix from Pauli-string basis to column-stack vec.

    M (d² × 4^N) such that vec(σ_α) = M[:, α]. M†M = 2^N · I (orthogonality).
    Cached by N; callers should treat the returned array as read-only.
    """
    d2 = 4 ** N
    d = 2 ** N
    M = np.zeros((d * d, d2), dtype=complex)
    for k in range(d2):
        indices = _k_to_indices(k, N)
        sigma = pauli_string(list(indices))
        M[:, k] = sigma.flatten('F')
    M.flags.writeable = False
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


def _build_kbody_chain(N, terms):
    """Build H = Σ_l Σ_term coeff · σ_{a_1}^l σ_{a_2}^{l+1} ... σ_{a_k}^{l+k-1}
    over a chain (sliding window of consecutive sites).

    F85 generalization of `_build_bilinear` to k-body Pauli terms on a chain.
    For each term (a_1, a_2, ..., a_k, coeff), sums the k-body Pauli string
    over all chain positions l = 0, 1, ..., N - k.

    Args:
        N: number of qubits
        terms: list of (a_1, a_2, ..., a_k, coeff) tuples where the k Pauli
            letters are followed by the coefficient. Mixed body counts are
            allowed (different terms can have different k).

    Returns:
        2^N × 2^N Hermitian matrix.
    """
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for term in terms:
        letters = term[:-1]
        coeff = term[-1]
        k = len(letters)
        if k > N:
            raise ValueError(
                f"k-body term of length {k} cannot fit on N={N} chain"
            )
        for l in range(N - k + 1):
            op = site_op(N, l, letters[0])
            for i in range(1, k):
                op = op @ site_op(N, l + i, letters[i])
            H = H + coeff * op
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


# ----------------------------------------------------------------------
# Single-excitation bonding-mode states (F65)
# ----------------------------------------------------------------------
# ψ_k(j) = √(2/(N+1)) · sin(πk(j+1)/(N+1)) is the j-th site amplitude of
# the k-th single-excitation eigenmode of the open uniform XX chain. The
# state vector embeds these amplitudes at single-excitation basis kets
# |1_j⟩ — site j → bit position N−1−j in the 2^N index (MSB convention,
# consistent with site_op's Kronecker order: site 0 = leftmost factor).

def bonding_mode_state(N, k):
    """F65 single-excitation bonding mode |ψ_k⟩ = Σ_j ψ_k(j) |1_j⟩.

    The handshake-algebra k coordinate (HANDSHAKE_ALGEBRA.md): selecting
    a value of k IS the per-side handshake. Two callers who independently
    construct `bonding_mode_state(N, k)` with the same k have agreed on
    the receiver — no exchange step needed.

    K-partner: bonding_mode_state(N, k_partner(N, k)) gives identical
    mirror-pair |·|²-observables under bipartite NN-hopping with real H.

    Args:
        N: chain length.
        k: bonding-mode index, 1 ≤ k ≤ N.

    Returns:
        Length-2^N normalized complex state vector.
    """
    if not 1 <= k <= N:
        raise ValueError(f"k={k} outside [1, N={N}]")
    psi = np.zeros(2 ** N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for j in range(N):
        psi[2 ** (N - 1 - j)] = norm * np.sin(np.pi * k * (j + 1) / (N + 1))
    return psi


def bonding_mode_pair_state(N, k):
    """Bonding-mode pair state (|vac⟩ + |ψ_k⟩) / √2.

    The canonical PTF / handshake initial state: a coherent superposition
    of the vacuum and the k-th bonding mode. Its single-qubit purity
    trajectories are the painter's canvas in PERSPECTIVAL_TIME_FIELD.md.
    """
    if not 1 <= k <= N:
        raise ValueError(f"k={k} outside [1, N={N}]")
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0
    norm = np.sqrt(2.0 / (N + 1))
    for j in range(N):
        psi[2 ** (N - 1 - j)] += norm * np.sin(np.pi * k * (j + 1) / (N + 1))
    psi /= np.linalg.norm(psi)
    return psi


# ----------------------------------------------------------------------
# Polarity-axis states (the +0 / 0 / −0 layer)
# ----------------------------------------------------------------------
# X-basis eigenstates carry the polarity-layer projection of the qubit:
# |+⟩ = (|0⟩+|1⟩)/√2 is the +0 polarity; |−⟩ = (|0⟩−|1⟩)/√2 is the −0
# polarity. Tensor products over N sites — |s_0, s_1, ..., s_{N-1}⟩ with
# s_j ∈ {+, −} — span the full polarity sublattice. THE_POLARITY_LAYER.md
# names this the bit_a-axis projection of the qubit; the d=0 axis lives
# at the center, between +0 and −0.

_X_PLUS = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)
_X_MINUS = np.array([1.0, -1.0], dtype=complex) / np.sqrt(2.0)


def polarity_state(N, signs):
    """N-site X-basis polarity state on the +0/−0 axis.

    Each site is in |+⟩ (sign +1, +0 polarity) or |−⟩ (sign −1, −0 polarity).
    The state is a tensor product following site_op's MSB convention:
    the signs argument is left-to-right = site 0 to site N-1 = leftmost
    Kronecker factor first.

    Args:
        N: chain length.
        signs: scalar ±1 (uniform polarity, |+⟩^N or |−⟩^N) or sequence
            of N values in {+1, −1}, one per site.

    Returns:
        Length-2^N normalized complex state vector.

    Examples:
        polarity_state(3, +1)            → |+,+,+⟩  (uniform +0)
        polarity_state(3, [+1,-1,+1])    → |+,−,+⟩  (X-basis Néel)
        polarity_state(3, [-1,-1,-1])    → |−,−,−⟩  (uniform −0)
    """
    if np.isscalar(signs):
        signs_list = [int(signs)] * N
    else:
        signs_list = list(signs)
    if len(signs_list) != N:
        raise ValueError(
            f"signs has length {len(signs_list)}, expected N={N}"
        )
    if not all(s in (1, -1) for s in signs_list):
        raise ValueError(f"signs must each be +1 or -1; got {signs_list}")
    factors = [_X_PLUS if s == 1 else _X_MINUS for s in signs_list]
    state = factors[0]
    for f in factors[1:]:
        state = np.kron(state, f)
    return state
