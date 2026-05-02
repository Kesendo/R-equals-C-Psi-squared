"""(n, n+1) U(1)-coherence block primitives for chromaticity / Q-scale analysis.

The Liouvillian L of an open XY chain under Z-dephasing preserves popcount,
so for any (n, n+1) coherence subspace L acts as a finite block of size
C(N, n) · C(N, n+1) — far smaller than the full d² = 4^N. These primitives
construct the block-restricted L, expose the J-coefficient per bond (for
J-perturbation observables like Q_SCALE / F86's K_CC_pr), and provide the
HD-channel-uniform projectors used in the chromaticity decomposition.

References:
- F74: chromaticity c(n, N) = min(n, N−1−n) + 1
- F73: c=1 spatial-sum coherence closure (the (0, 1)-block special case)
- F86: Q_peak chromaticity-specific N-invariant constants
- Q_SCALE_THREE_BANDS: per-block band structure

Convention: site 0 = MSB (big-endian); H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1})
matches `ChainSystem(H_type='xy')`. For Heisenberg, the diagonal Z·Z part
adds to D within each popcount sector but the U(1)-block structure is
identical.

Public API:
  popcount_states(N, n)
  block_basis(N, n)                       returns (P_n, P_np1)
  chromaticity(N, n)                      returns c = min(n, N−1−n) + 1
  block_L_split_xy(N, n, gamma_0)         returns (D, [M_H_b], P_n, P_np1)
  hd_channel_basis(N, n)                  returns (P, [HD_values])
  dicke_block_probe(N, n)                 returns flat (Mn·Mnp1)-vector
  spatial_sum_coherence_kernel(N, n)      returns Mtot×Mtot quadratic form

F86 EP-derived constants:
  t_peak(gamma_0)                         = 1/(4γ₀)        universal EP time

(Earlier closed-form conjectures `q_peak_endpoint(N) = csc(π/(N+1))` and
`Q_PEAK_INTERIOR_C3_ANCHOR = csc(π/5)` were retracted after extended-N
data showed both were N=7-specific coincidence-matches; see F86 in
ANALYTICAL_FORMULAS.md and PROOF_F86_QPEAK.md for the rollback details.)
"""
from __future__ import annotations

from math import comb, pi, sin

import numpy as np


def t_peak(gamma_0):
    """E-folding time of the EP-degenerate Liouvillian eigenvalue at Q_peak.
    Universal across c, N, n, bond position: t_peak = 1/(4γ₀).

    The exceptional point between adjacent pure-rate channels (gap
    Δ = 4γ₀) has degenerate eigenvalue Re(λ) = −4γ₀ at J·g_eff = 2γ₀.
    The probe's amplitude on that mode decays as exp(−4γ₀·t), so the
    natural timescale of the J-derivative peak is 1/(4γ₀). See F86.
    """
    return 1.0 / (4.0 * gamma_0)


def popcount_states(N, n):
    """List of computational basis states (as ints) with popcount n.

    Big-endian: site 0 = MSB. Sorted by integer value.
    """
    return [x for x in range(2 ** N) if bin(x).count("1") == n]


def block_basis(N, n):
    """Basis state lists (P_n, P_np1) for the (n, n+1) coherence block."""
    return popcount_states(N, n), popcount_states(N, n + 1)


def _block_index_maps(N, n):
    """(P_n, P_np1, p_to_idx, q_to_idx, Mnp1, Mtot) — the boilerplate every
    block-level constructor needs. Big-endian flat index is
    `p_to_idx[p] * Mnp1 + q_to_idx[q]`.
    """
    P_n = popcount_states(N, n)
    P_np1 = popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    return (P_n, P_np1,
            {p: i for i, p in enumerate(P_n)},
            {q: i for i, q in enumerate(P_np1)},
            Mnp1, len(P_n) * Mnp1)


def _bond_flip_targets(state, N):
    """Yield (bond, flipped_state) for every bond where the state's two
    adjacent bits are opposite (i.e. an XX+YY swap is non-trivial).
    Big-endian: site 0 = MSB.
    """
    for bond in range(N - 1):
        bit_hi = (state >> (N - 1 - bond)) & 1
        bit_lo = (state >> (N - 2 - bond)) & 1
        if bit_hi != bit_lo:
            mask = (1 << (N - 1 - bond)) | (1 << (N - 2 - bond))
            yield bond, state ^ mask


def chromaticity(N, n):
    """F74 chromaticity: number of distinct pure dephasing rates in the
    (n, n+1) coherence block at J = 0.

        c(n, N) = min(n, N − 1 − n) + 1

    Hamming-distance values HD ∈ {1, 3, ..., 2c−1} are realised; pure rates
    are 2γ₀·HD.
    """
    return min(n, N - 1 - n) + 1


def block_L_split_xy(N, n, gamma_0):
    """Construct the (n, n+1) coherence block of L for the uniform XY chain
    under Z-dephasing, split into a J-independent dephasing part and per-bond
    H coefficients.

    For uniform J across all bonds:  L = D + J · sum(M_H_per_bond).

    Args:
        N: chain length.
        n: lower popcount of the coherence block. Must satisfy 0 ≤ n ≤ N−1.
        gamma_0: uniform Z-dephasing rate per site.

    Returns:
        D: (Mtot, Mtot) diagonal complex matrix, D_{i,i} = −2γ₀ · HD(p, q),
           Mtot = C(N, n) · C(N, n+1).
        M_H_per_bond: list of N−1 (Mtot, Mtot) complex matrices. M_H_per_bond[b]
                      is the J-coefficient contribution to L from bond b.
                      Entries are ±1j on the (p ↔ p XOR mask_b, q) and
                      (p, q ↔ q XOR mask_b) couplings where the bond's bits
                      are opposite.
        P_n: list of popcount-n basis states.
        P_np1: list of popcount-(n+1) basis states.

    Convention: H_b = (J/2)·(X_b X_{b+1} + Y_b Y_{b+1}) on adjacent sites,
    big-endian. H_b|...01...⟩ = J|...10...⟩ when the bond's bits are opposite.
    """
    if not 0 <= n <= N - 1:
        raise ValueError(f"n must be in [0, N-1]; got n={n}, N={N}")

    P_n, P_np1, p_to_idx, q_to_idx, Mnp1, Mtot = _block_index_maps(N, n)

    D = np.zeros((Mtot, Mtot), dtype=complex)
    M_H_per_bond = [np.zeros((Mtot, Mtot), dtype=complex) for _ in range(N - 1)]

    p_flips = {p: list(_bond_flip_targets(p, N)) for p in P_n}
    q_flips = {q: list(_bond_flip_targets(q, N)) for q in P_np1}

    for p in P_n:
        for q in P_np1:
            i = p_to_idx[p] * Mnp1 + q_to_idx[q]
            D[i, i] = -2.0 * gamma_0 * bin(p ^ q).count("1")

            # H acts -i[H, |p⟩⟨q|]: -i·H|p⟩⟨q| (left, sign -1j) and
            # +i·|p⟩⟨q|·H (right, sign +1j). Each bond flips a single
            # adjacent-bit pair when the two bits are opposite.
            for bond, p_flipped in p_flips[p]:
                j = p_to_idx[p_flipped] * Mnp1 + q_to_idx[q]
                M_H_per_bond[bond][j, i] += -1j
            for bond, q_flipped in q_flips[q]:
                j = p_to_idx[p] * Mnp1 + q_to_idx[q_flipped]
                M_H_per_bond[bond][j, i] += 1j

    return D, M_H_per_bond, P_n, P_np1


def hd_channel_basis(N, n):
    """Orthonormalised HD-channel-uniform projectors for the (n, n+1) block.

    Each channel k ∈ 0..c−1 collects all basis pairs (p, q) with
    HD(p, q) = 2k+1 and assigns them equal weight. The result |c_k⟩ is unit-
    normalised. These projectors diagonalise M_H_eff in the c-dim subspace:
    ⟨c_j | L_H | c_k⟩ = 0 for j ≠ k, a structural extension of F73
    (c = 1 case at (0, 1) block) to all c.

    Returns:
        P: (Mtot, c) complex matrix, columns are the orthonormal channel
           vectors.
        HDs: list of c HD values [1, 3, 5, ..., 2c−1].
    """
    P_n, P_np1, p_to_idx, q_to_idx, Mnp1, Mtot = _block_index_maps(N, n)
    c = chromaticity(N, n)
    HDs = [2 * k + 1 for k in range(c)]

    P = np.zeros((Mtot, c), dtype=complex)
    for k, hd in enumerate(HDs):
        count = 0
        for p in P_n:
            for q in P_np1:
                if bin(p ^ q).count("1") == hd:
                    P[p_to_idx[p] * Mnp1 + q_to_idx[q], k] = 1.0
                    count += 1
        if count > 0:
            P[:, k] /= np.sqrt(count)
    return P, HDs


def dicke_block_probe(N, n):
    """Dicke coherence probe |S_n⟩⟨S_{n+1}| / 2 expressed as a flat coefficient
    vector in the (n, n+1)-block basis.

    |S_n⟩ = (1/√C(N, n))·Σ_{popcount(p)=n} |p⟩, hence the Dicke probe has
    every basis-element coefficient equal to 1/(2·√(C(N,n)·C(N,n+1))).
    Used as initial state ρ₀ for Q_SCALE / F86 K_CC_pr observables.
    """
    Mn = comb(N, n)
    Mnp1 = comb(N, n + 1)
    Mtot = Mn * Mnp1
    coeff = 1.0 / (2.0 * np.sqrt(Mn * Mnp1))
    return np.full(Mtot, coeff, dtype=complex)


def spatial_sum_coherence_kernel(N, n):
    """F73-style spatial-sum coherence kernel as a quadratic form on the
    flat (n, n+1)-block coefficient vector.

    S(t) = Σᵢ 2 · |(ρᵢ(t))_{0,1}|²

    can be written as ρ(t)† · S_kernel · ρ(t) with S_kernel constructed here.
    For each site i, the (0, 1) element of the i-th reduced density matrix
    is a linear functional of the block coefficients (the (p, q) pairs that
    differ ONLY at site i and agree elsewhere). |·|² becomes A†A in matrix
    form, summed over sites.

    Returns:
        S_kernel: (Mtot, Mtot) Hermitian positive semi-definite complex matrix.
    """
    P_n, _, p_to_idx, q_to_idx, Mnp1, Mtot = _block_index_maps(N, n)

    S = np.zeros((Mtot, Mtot), dtype=complex)
    for i_site in range(N):
        A = np.zeros(Mtot, dtype=complex)
        mask_i = 1 << (N - 1 - i_site)
        for p in P_n:
            if (p >> (N - 1 - i_site)) & 1:
                continue
            q = p | mask_i
            if bin(q).count("1") != n + 1:
                continue
            A[p_to_idx[p] * Mnp1 + q_to_idx[q]] = 1.0
        S += 2.0 * np.outer(A.conj(), A)
    return S
