"""F82 closed-form prediction and inversion of F81 violation under T1 amplitude damping."""
from __future__ import annotations

import numpy as np


def predict_T1_dissipator_violation(chain, gamma_t1_l):
    """F82 closed form: predict the F81 violation from per-site T1 rates.

    Theorem F82 (proved in PROOF_F82_T1_DISSIPATOR_CORRECTION):

        ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1)

    where D_{T1, odd} is the Π²-anti-symmetric part of the T1 amplitude-
    damping dissipator, and ‖·‖_F is the Frobenius norm in the framework's
    Pauli basis. This equals the f81_violation reported by
    `pi_decompose_M(gamma_t1=...)`. The result is Hamiltonian-independent
    and γ_z-independent (Master Lemma propagates through F82).

    Args:
        gamma_t1_l: per-site T1 amplitude-damping rates (list of length N
            or scalar for uniform).

    Returns:
        Predicted ‖D_{T1, odd}‖_F as float. For uniform γ_T1 across all N
        sites this simplifies to γ_T1 · √N · 2^(N−1).

    Use case: predict the F81 violation expected for a given hardware T1
    profile, then compare with `pi_decompose_M` (or with measured-then-
    fitted Lindblad data) to validate the F82 model.
    """
    if np.isscalar(gamma_t1_l):
        gt1 = [float(gamma_t1_l)] * chain.N
    else:
        gt1 = list(gamma_t1_l)
        if len(gt1) != chain.N:
            raise ValueError(
                f"gamma_t1_l must have length {chain.N}, got {len(gt1)}"
            )
    sum_sq = sum(g * g for g in gt1)
    return float(np.sqrt(sum_sq) * (2 ** (chain.N - 1)))


def estimate_T1_from_violation(chain, f81_violation):
    """F82 inverse closed form: extract RMS γ_T1 from a measured F81 violation.

    Given a measured (or fitted-Lindblad-implied) f81_violation, F82
    inverts the closed-form scaling to recover the RMS per-site T1 rate:

        γ_T1, RMS = √((Σ_l γ²_T1_l) / N) = f81_violation / (√N · 2^(N−1)).

    For uniform per-site T1, this returns the actual γ_T1. For non-uniform
    rates, the result is the root-mean-square across sites; per-site rates
    cannot be disentangled from the F81 violation alone (additional probes
    are needed to break the degeneracy).

    Args:
        f81_violation: ‖M_anti − L_{H_odd}‖_F as returned by `pi_decompose_M`
            (or computed from process-tomography data).

    Returns:
        RMS γ_T1 across the N sites, as float.

    Use case: hardware T1-rate readout from a noise-channel-blind L fit.
    For the Marrakesh dataset (joint fit gives γ_T1 ≈ 0): F82 inverse
    recovers γ_T1, RMS ≈ 0, consistent with the empirical refutation of
    the T1 amplification hypothesis.
    """
    if f81_violation < 0:
        raise ValueError(f"f81_violation must be non-negative, got {f81_violation}")
    return float(f81_violation / (np.sqrt(chain.N) * (2 ** (chain.N - 1))))
