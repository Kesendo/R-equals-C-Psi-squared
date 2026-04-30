"""Optimal γ-Sensing setup + γ-extraction from CΨ measurement (Bell+, multi-axis channels)."""
from __future__ import annotations

import numpy as np


def gamma_probe_setup(chain, gamma_assumed=None, target_precision=0.01,
                       channel='Z'):
    """Optimal γ-Sensing-Parameter via Cusp-nahe CΨ-Probe (Bell+).

    Generalizes F25 → F26 (multi-axis Pauli channels). The cusp at CΨ=1/4
    is a good probe region but not the Fisher-Information optimum; the
    optimum lies post-cusp.

    Channels:
        'Z'             — pure Z-dephasing only (γ_z = γ; F25 special case)
        'X'             — pure X-noise (γ_x = γ)
        'Y'             — pure Y-noise (γ_y = γ; symmetric with X for Bell+)
        'depolarizing'  — γ/3 on each axis

    K = γ·t is γ-invariant; channel determines K_optimal numerically.
    Reference K_cusp values per channel (F27): Z=0.0374, X=0.0867,
    Y=0.0867, depolarizing=0.0440.

    Args:
        gamma_assumed: prior γ-estimate. If None, uses self.gamma_0.
        target_precision: relative precision Δγ/γ desired.
        channel: noise channel — one of 'Z', 'X', 'Y', 'depolarizing'.

    Returns:
        dict with t_optimal, cpsi_target, K_optimal, fisher_per_shot,
        shots_needed, plus cusp parameters and the channel itself.
    """
    from scipy.optimize import minimize_scalar, brentq
    from ..lindblad import cpsi_bell_plus

    if gamma_assumed is None:
        gamma_assumed = chain.gamma_0
    gamma = float(gamma_assumed)

    # Map channel → axis assignment of γ
    def gammas_for_channel(g):
        if channel == 'Z':
            return (0.0, 0.0, g)
        if channel == 'X':
            return (g, 0.0, 0.0)
        if channel == 'Y':
            return (0.0, g, 0.0)
        if channel == 'depolarizing':
            return (g/3, g/3, g/3)
        raise ValueError(f"channel must be 'Z'/'X'/'Y'/'depolarizing'; got {channel!r}")

    gx, gy, gz = gammas_for_channel(gamma)

    def cpsi(t):
        return cpsi_bell_plus(gx, gy, gz, t)

    def dcpsi_dgamma(t):
        # Numerical derivative — analytic form is messy for general channel.
        eps = max(1e-7, 1e-6 * gamma)
        gx_p, gy_p, gz_p = gammas_for_channel(gamma + eps)
        gx_m, gy_m, gz_m = gammas_for_channel(gamma - eps)
        return (cpsi_bell_plus(gx_p, gy_p, gz_p, t)
                - cpsi_bell_plus(gx_m, gy_m, gz_m, t)) / (2 * eps)

    def neg_fisher(t):
        c = cpsi(t)
        dc = dcpsi_dgamma(t)
        return -(dc * dc / max(1 - c * c, 1e-12))

    res = minimize_scalar(neg_fisher,
                          bounds=(0.001 / gamma, 5.0 / gamma),
                          method='bounded')
    t_opt = float(res.x)
    cpsi_opt = float(cpsi(t_opt))
    fisher_opt = -float(res.fun)

    # Cusp t (where CΨ = 1/4)
    try:
        t_cusp = brentq(lambda t: cpsi(t) - 0.25, 1e-6, 1000.0 / gamma)
        cpsi_cusp = 0.25
    except ValueError:
        t_cusp = float('nan')
        cpsi_cusp = float('nan')

    delta_gamma = float(target_precision) * gamma
    shots = int(np.ceil(1.0 / (fisher_opt * delta_gamma * delta_gamma)))

    return {
        'channel':          channel,
        't_optimal':        t_opt,
        'cpsi_target':      cpsi_opt,
        'K_optimal':        gamma * t_opt,
        'fisher_per_shot':  fisher_opt,
        'shots_needed':     shots,
        'target_precision': float(target_precision),
        't_cusp':           float(t_cusp),
        'K_cusp':           gamma * float(t_cusp) if not np.isnan(t_cusp) else float('nan'),
        'cpsi_cusp':        cpsi_cusp,
        'gamma_assumed':    gamma,
    }


def estimate_gamma_from_cpsi(chain, cpsi_measured, t, channel='Z'):
    """Invert F26 to extract γ from CΨ at time t (Bell+, multi-axis channels).

    For 'Z' channel uses F25 cubic root (closed form). For 'X', 'Y',
    'depolarizing' uses numerical inversion of F26.

    Bell+ initial CΨ(t=0) = 1/3, so cpsi_measured must satisfy
    0 < cpsi_measured < 1/3.

    Args:
        cpsi_measured: measured CΨ value (e.g. from hardware tomography).
        t: probe time at which CΨ was measured.
        channel: 'Z' (default), 'X', 'Y', or 'depolarizing'.

    Returns:
        Estimated γ as float.
    """
    from scipy.optimize import brentq
    from ..lindblad import cpsi_bell_plus
    cpsi = float(cpsi_measured)
    if cpsi >= 1.0 / 3.0:
        raise ValueError(
            f"cpsi_measured = {cpsi} ≥ 1/3 = CΨ(t=0); cannot invert "
            f"(Bell+ never starts above this)."
        )
    if cpsi <= 0:
        raise ValueError(
            f"cpsi_measured = {cpsi} ≤ 0; physically invalid for Bell+."
        )
    if t <= 0:
        raise ValueError(f"t must be > 0; got {t}")

    if channel == 'Z':
        f = brentq(lambda fv: fv * (1 + fv * fv) / 6.0 - cpsi, 1e-12, 1.0 - 1e-12)
        return -float(np.log(f)) / (4.0 * float(t))

    # General channel: numerical inversion of F26
    def gammas_for_channel(g):
        if channel == 'X':            return (g, 0.0, 0.0)
        if channel == 'Y':            return (0.0, g, 0.0)
        if channel == 'depolarizing': return (g/3, g/3, g/3)
        raise ValueError(f"channel must be 'Z'/'X'/'Y'/'depolarizing'; got {channel!r}")

    def f(gamma):
        gx, gy, gz = gammas_for_channel(gamma)
        return cpsi_bell_plus(gx, gy, gz, t) - cpsi

    # CΨ monotonically decreases in γ ⇒ use bisection
    gamma_est = brentq(f, 1e-12, 1000.0 / float(t))
    return float(gamma_est)
