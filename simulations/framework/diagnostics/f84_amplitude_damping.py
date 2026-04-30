"""F84 closed-form prediction and inversion of F81 violation under thermal amplitude damping."""
from __future__ import annotations

import numpy as np


def predict_amplitude_damping_violation(chain, gamma_t1_l, gamma_pump_l=None):
    """F84 closed form: predict the F81 violation from cooling + heating rates.

    Theorem F84 (proved in PROOF_F84_AMPLITUDE_DAMPING) generalizes F82
    to thermal amplitude damping with both cooling (σ⁻) and heating (σ⁺):

        ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)

    The F81 violation depends only on the *net* cooling rate Δγ_l =
    γ_↓_l − γ_↑_l. At thermal equilibrium γ_↓ = γ_↑ (detailed balance),
    the violation is zero even though both channels are active. The
    violation isolates the vacuum-fluctuation contribution.

    Pauli-channel dissipators D[Z], D[X], D[Y] do not contribute to
    f81_violation (Pauli-Channel Cancellation Lemma in PROOF_F84): they
    are Π²-symmetric. F84 violation is exclusive to σ⁻/σ⁺ channels.

    Args:
        gamma_t1_l: per-site cooling rates γ_↓_l (σ⁻ amplitude damping).
            Scalar (uniform) or list of length N.
        gamma_pump_l: per-site heating rates γ_↑_l (σ⁺ amplitude damping).
            Scalar (uniform), list of length N, or None (= no heating, pure F82).

    Returns:
        Predicted ‖D_{AmplDamp, odd}‖_F as float. Generalizes F82's
        `predict_T1_dissipator_violation` (which is the special case
        gamma_pump_l = None / 0).

    Use case: predict the F81 violation expected for a given thermal
    amplitude-damping profile (cooling + heating rates per site), then
    compare with `pi_decompose_M(gamma_t1, gamma_pump)` or with
    process-tomography-derived L on hardware. At any temperature, the
    violation reads out the temperature-independent vacuum component.
    """
    if np.isscalar(gamma_t1_l):
        gt1 = [float(gamma_t1_l)] * chain.N
    else:
        gt1 = list(gamma_t1_l)
        if len(gt1) != chain.N:
            raise ValueError(
                f"gamma_t1_l must have length {chain.N}, got {len(gt1)}"
            )
    if gamma_pump_l is None:
        gp = [0.0] * chain.N
    elif np.isscalar(gamma_pump_l):
        gp = [float(gamma_pump_l)] * chain.N
    else:
        gp = list(gamma_pump_l)
        if len(gp) != chain.N:
            raise ValueError(
                f"gamma_pump_l must have length {chain.N}, got {len(gp)}"
            )
    sum_sq = sum((d - u) ** 2 for d, u in zip(gt1, gp))
    return float(np.sqrt(sum_sq) * (2 ** (chain.N - 1)))


def estimate_net_cooling_from_violation(chain, f81_violation):
    """F84 inverse closed form: extract RMS net cooling rate from f81_violation.

    Inverts F84:

        |Δγ|_RMS = √(Σ_l (γ_↓_l − γ_↑_l)² / N) = f81_violation / (√N · 2^(N−1))

    For uniform γ_↓, γ_↑: returns the actual Δγ = γ_↓ − γ_↑. For
    non-uniform rates: returns the RMS net rate; per-site rates and
    the cooling/heating split require additional probes.

    At any temperature, this recovers the vacuum-fluctuation amplitude:
    thermal photon balance (γ_↓ ↔ γ_↑) cancels under the F81 anti-symmetric
    projection, leaving only the spontaneous-emission rate.

    Equivalent to `estimate_T1_from_violation` when the bath is at T=0
    (γ_↑ = 0, hence Δγ = γ_↓).

    Args:
        f81_violation: ‖M_anti − L_{H_odd}‖_F as returned by `pi_decompose_M`.

    Returns:
        RMS |γ_↓ − γ_↑| across the N sites, as float ≥ 0.
    """
    if f81_violation < 0:
        raise ValueError(f"f81_violation must be non-negative, got {f81_violation}")
    return float(f81_violation / (np.sqrt(chain.N) * (2 ** (chain.N - 1))))
