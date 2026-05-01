"""Slow-mode lens analysis: extract optimal single-excitation receiver amplitudes.

The Lindbladian L's slowest non-stationary modes carry the longest-lived
coherences. For chains, the dominant slow modes typically project onto
the single-excitation (SE) sector. The "lens state" is the single-qubit
amplitude pattern that maximally aligns with the slowest SE-accessible
mode — the optimal receiver/sensor for that decoherence regime.

Python parallel of the C# `LensAnalysis.cs` pipeline (compute/RCPsiSquared.Compute)
for N ≤ 6 where dense numpy eigendecomposition is feasible. The C# version
handles N = 7-8 with native memory and ILP64 LAPACK.

Used by experiments/CUSP_LENS_CONNECTION.md and several Python lens scripts
(slow_mode_lens_analysis.py, slow_mode_shapes_vs_J.py, n7_central_defect_check.py,
observer_time_slow_mode_analysis.py, kingston_gamma0_test.py, etc.).

Public API:
  slow_modes(chain, n_slow=3, ...)
  lens_pipeline(chain, n_slow=3, se_threshold=0.01)
"""
from __future__ import annotations

import numpy as np


def _se_indices(N):
    """Single-excitation basis indices in the d=2^N computational basis (MSB convention).

    Qubit k (k=0..N-1) flipped: index 2^(N-1-k). All other qubits = 0.
    """
    return [1 << (N - 1 - k) for k in range(N)]


def _extract_se_block(eigvec_vec, N):
    """Extract the N×N single-excitation block from a vec'd d² operator.

    Args:
        eigvec_vec: vec-form eigenvector of L (length d² = 4^N), row-major
            convention (vec[i*d + j] = M[i, j]).
        N: chain length.

    Returns:
        (block, frob_ratio): N×N complex matrix and ratio
        ‖block‖_F / ‖eigvec‖_F (how much of the eigenvec lives in SE).
    """
    d = 2 ** N
    se_idx = _se_indices(N)
    block = np.zeros((N, N), dtype=complex)
    for a in range(N):
        for b in range(N):
            block[a, b] = eigvec_vec[se_idx[a] * d + se_idx[b]]
    full_norm = float(np.linalg.norm(eigvec_vec))
    block_norm = float(np.linalg.norm(block))
    frob_ratio = block_norm / full_norm if full_norm > 1e-30 else 0.0
    return block, frob_ratio


def _lens_amplitudes(se_block):
    """Optimal single-qubit amplitudes that diagonalize the Hermitized SE block.

    Construct M = (Q + Q†) / 2 with Q = se_block.T, find the eigenvector
    with largest |eigenvalue|. The resulting amplitudes are the receiver
    pattern that maximally projects onto the slow mode.

    Returns:
        (amplitudes, eigenvalue, projection): real positive amplitudes
        (length N, normalized), the dominant eigenvalue (real), and
        ⟨a|Q|a⟩ projection magnitude.
    """
    Q = se_block.T  # match C# LensAnalysis convention: Q[i,j] = block[j,i]
    M = (Q + Q.conj().T) / 2.0
    evals, evecs = np.linalg.eigh(M)
    # eigh returns ascending real evals; we want largest |eigenvalue|
    abs_evals = np.abs(evals)
    best = int(np.argmax(abs_evals))
    amps = np.abs(evecs[:, best])
    norm = float(np.linalg.norm(amps))
    if norm > 1e-30:
        amps = amps / norm
    a_vec = amps.astype(complex)
    proj = a_vec.conj() @ Q @ a_vec
    return amps, float(evals[best].real), float(np.abs(proj))


def _default_n_slow(N):
    """Default number of slow modes to scan: min(d²/4, 50), matching the C# lens pipeline."""
    return min(4 ** N // 4, 50)


def slow_modes(chain, n_slow=None, exclude_stationary=True, stationary_tol=1e-10,
               L=None):
    """Find n_slow slowest non-stationary Liouvillian modes.

    Eigendecomposes the Liouvillian (L override, else chain.L) in vec form.
    Modes are sorted by |Re(eigenvalue)| ascending; the stationary mode
    (Re ≈ 0) is excluded by default.

    Args:
        chain: ChainSystem (provides N for shape consistency).
        n_slow: number of slow modes to return. Default min(d²/4, 50)
            matches the C# `lens` pipeline scan width.
        exclude_stationary: drop modes with |Re(λ)| < stationary_tol.
        stationary_tol: threshold for "stationary".
        L: optional 4^N × 4^N Liouvillian override. If None, uses chain.L
            (the cached uniform-γ Heisenberg-on-chain Liouvillian).
            Pass a custom L to study non-uniform γ profiles, transverse
            fields, T1 noise, etc.

    Returns:
        dict with keys:
          'eigenvalues': complex array of length n_slow (slowest first)
          'right_eigvecs': d² × n_slow array, columns are right eigenvectors of L
          'left_covecs':  n_slow × d² array, rows are left covectors (rows of R⁻¹)
          'indices': original indices into eig(L)'s output
          'rates': real array, -Re(eigvals[i]) (the decay rates, all > 0)
    """
    if n_slow is None:
        n_slow = _default_n_slow(chain.N)
    if L is None:
        L = chain.L
    evals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    re_parts = -np.real(evals)  # decay rates: positive for stable modes
    candidates = np.arange(len(evals))
    if exclude_stationary:
        candidates = candidates[np.abs(re_parts[candidates]) > stationary_tol]
    order = candidates[np.argsort(re_parts[candidates])]
    sel = order[:n_slow]
    return {
        'eigenvalues': evals[sel],
        'right_eigvecs': R[:, sel],
        'left_covecs': R_inv[sel, :],
        'indices': sel,
        'rates': re_parts[sel],
    }


def lens_pipeline(chain, n_slow=None, se_threshold=0.01, L=None):
    """Full lens analysis: find slowest SE-accessible mode and its amplitudes.

    Walks the slow modes from slowest upward; for each, extracts the SE-block
    of its left covector (the "drainage direction" of the mode in operator
    space) and computes the SE Frobenius ratio. The first mode with
    ratio > se_threshold is taken as the lens-accessible slow mode. The
    next non-SE mode is reported as inaccessible.

    Args:
        chain: ChainSystem.
        n_slow: how deep to scan.
        se_threshold: minimum SE-block Frobenius ratio to qualify as
            "lens-accessible".
        L: optional Liouvillian override (see `slow_modes`).

    Returns:
        dict with:
          'lens_amplitudes': real array length N, normalized (None if no
              SE-accessible mode found within n_slow).
          'lens_eigenvalue': dominant SE-block eigenvalue (real).
          'lens_projection': ⟨a|Q|a⟩ magnitude.
          'lens_rate': -Re(eigval) of the lens-accessible mode.
          'lens_se_ratio': SE Frobenius ratio of the chosen mode.
          'second_mode_accessible': True if a second slow mode is also
              SE-accessible.
          'second_mode_rate': decay rate of the second slow mode (NaN if
              none scanned).
          'slow_modes': raw output of `slow_modes(chain, n_slow)`.
    """
    sm = slow_modes(chain, n_slow=n_slow, L=L)
    N = chain.N
    lens_amps = None
    lens_evalue = float('nan')
    lens_proj = float('nan')
    lens_rate = float('nan')
    lens_ratio = 0.0
    second_accessible = False
    second_rate = float('nan')
    found_lens = False
    for k in range(len(sm['eigenvalues'])):
        left_covec = sm['left_covecs'][k]
        block, ratio = _extract_se_block(left_covec, N)
        if ratio > se_threshold:
            if not found_lens:
                lens_amps, lens_evalue, lens_proj = _lens_amplitudes(block)
                lens_rate = float(sm['rates'][k])
                lens_ratio = ratio
                found_lens = True
            else:
                second_accessible = True
                second_rate = float(sm['rates'][k])
                break
        elif found_lens:
            second_rate = float(sm['rates'][k])
            break
    return {
        'lens_amplitudes': lens_amps,
        'lens_eigenvalue': lens_evalue,
        'lens_projection': lens_proj,
        'lens_rate': lens_rate,
        'lens_se_ratio': lens_ratio,
        'second_mode_accessible': second_accessible,
        'second_mode_rate': second_rate,
        'slow_modes': sm,
    }
