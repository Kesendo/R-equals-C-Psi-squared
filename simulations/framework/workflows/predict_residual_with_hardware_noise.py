"""Predict ‖M‖² for an IBM-like hardware noise model in closed form."""
from __future__ import annotations

from ..lindblad import (
    HARDWARE_DISSIPATORS,
    HARDWARE_DISSIPATOR_D1,
)


def predict_residual_with_hardware_noise(chain, terms=None,
                                          T1_l=None, Tphi_l=None,
                                          T1pump_l=None,
                                          Xnoise_l=None, Ynoise_l=None):
    """Predict ‖M‖² for an IBM-like hardware noise model in closed form.

    Aggregates the per-class Frobenius forms and the cross-class d1, d2
    terms from HARDWARE_DISSIPATORS / HARDWARE_DISSIPATOR_D1. Returns a
    dict with the per-class contributions plus the total.

    All rate arguments are per-site lists (length N) of floats. Pass None
    to omit a class. Pass terms=None for purely-dissipative ‖M‖².

    σ_offset = 0 (raw Lindbladian L; do NOT subtract Π-trivial Z-dephasing
    offset, since that's already absorbed into c1=0 for Tphi/Ynoise).

    Args:
        terms: optional Pauli-pair Hamiltonian terms; if None uses no H.
        T1_l, Tphi_l, T1pump_l, Xnoise_l, Ynoise_l: per-site rate lists or
            None (= zero). Lists must be length N.

    Returns:
        dict with keys 'hamiltonian', 'per_class' (sub-dict), 'cross'
        (sub-dict), and 'total' (the sum).
    """
    from ..diagnostics.f49_frobenius_scaling import (
        predict_residual_norm_squared_from_terms,
    )
    rate_args = {
        'T1':     T1_l,
        'T1pump': T1pump_l,
        'Tphi':   Tphi_l,
        'Xnoise': Xnoise_l,
        'Ynoise': Ynoise_l,
    }
    rates = {}
    for name, gl in rate_args.items():
        if gl is None or all(g == 0 for g in gl):
            continue
        if len(gl) != chain.N:
            raise ValueError(
                f"{name}_l length {len(gl)} != N {chain.N}"
            )
        rates[name] = list(gl)

    factor = 4 ** (chain.N - 1)

    # Hamiltonian palindrome part — reuse predict_residual_norm_squared_from_terms
    if terms:
        h_part = predict_residual_norm_squared_from_terms(chain, terms)
    else:
        h_part = 0.0

    # Per-class contributions
    per_class = {}
    for name, gl in rates.items():
        spec = HARDWARE_DISSIPATORS[name]
        sum_g_sq = sum(g * g for g in gl)
        sum_g = sum(gl)
        contrib = factor * (spec['c1'] * sum_g_sq + spec['c2'] * sum_g * sum_g)
        per_class[name] = contrib

    # Cross-class contributions
    cross = {}
    active = list(rates.keys())
    for i, k1 in enumerate(active):
        for j, k2 in enumerate(active):
            if i >= j:
                continue
            key = (k1, k2) if (k1, k2) in HARDWARE_DISSIPATOR_D1 else (k2, k1)
            d1 = HARDWARE_DISSIPATOR_D1[key]
            p1 = HARDWARE_DISSIPATORS[k1]['pauli']
            p2 = HARDWARE_DISSIPATORS[k2]['pauli']
            d2 = 32.0 * (sum(abs(x)**2 for x in p1)) * (sum(abs(x)**2 for x in p2))
            gl1, gl2 = rates[k1], rates[k2]
            sum_g1_g2 = sum(g1 * g2 for g1, g2 in zip(gl1, gl2))
            prod_sums = sum(gl1) * sum(gl2)
            contrib = factor * (d1 * sum_g1_g2 + d2 * prod_sums)
            cross[(k1, k2)] = contrib

    total = h_part + sum(per_class.values()) + sum(cross.values())
    return {
        'hamiltonian': h_part,
        'per_class':   per_class,
        'cross':       cross,
        'total':       total,
    }
