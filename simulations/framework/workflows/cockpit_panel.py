"""Full Lebensader cockpit panel orchestrator (skeleton + trace + cusp + chiral + Y-parity)."""
from __future__ import annotations

import numpy as np

from ..pauli import _build_bilinear
from ..lebensader import cockpit_panel as _cockpit_panel_primitive


def cockpit_panel(chain, receiver, terms=None, gamma_t1=None,
                  t_max=10.0, dt=0.005,
                  threshold=1e-9, cluster_tol=1e-8):
    """Full Lebensader cockpit panel: skeleton + trace + cusp + chiral + Y-parity.

    Args:
        receiver: Receiver instance providing ρ_0.
        terms: optional Pauli-pair Hamiltonian terms; if None uses chain's default H.
        gamma_t1: scalar or list of length N for T1 amplitude-damping rates.
        t_max, dt: time grid.
    """
    if receiver.N != chain.N:
        raise ValueError(
            f"receiver.N ({receiver.N}) does not match chain.N ({chain.N})"
        )

    if terms is not None:
        bilinear = [(t[0], t[1], chain.J) for t in terms]
        H = _build_bilinear(chain.N, chain.bonds, bilinear)
    else:
        H = chain.H

    gamma_l = [chain.gamma_0] * chain.N
    if gamma_t1 is None:
        gamma_t1_l = [0.0] * chain.N
    elif np.isscalar(gamma_t1):
        gamma_t1_l = [float(gamma_t1)] * chain.N
    else:
        gamma_t1_l = list(gamma_t1)
        if len(gamma_t1_l) != chain.N:
            raise ValueError(
                f"gamma_t1 list length {len(gamma_t1_l)} != N {chain.N}"
            )

    return _cockpit_panel_primitive(
        H, gamma_l, receiver.rho, chain.N,
        gamma_t1_l=gamma_t1_l,
        t_max=t_max, dt=dt,
        threshold=threshold, cluster_tol=cluster_tol,
    )
