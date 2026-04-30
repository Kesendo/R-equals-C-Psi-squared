"""Tests for cpsi_bell_plus formula and the F27 cusp K-values per channel."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_cpsi_bell_plus_recovers_f25_for_pure_z():
    """F26 with γ_x=γ_y=0 must reduce to F25: CΨ = f·(1+f²)/6 with f=exp(-4γz·t)."""
    gz = 0.05
    t = 2.0
    f = np.exp(-4 * gz * t)
    cpsi_f25 = f * (1 + f**2) / 6.0
    cpsi_f26 = fw.cpsi_bell_plus(0.0, 0.0, gz, t)
    assert abs(cpsi_f26 - cpsi_f25) < 1e-12


def test_cpsi_bell_plus_at_t0_gives_one_third():
    """Bell+ at t=0: CΨ = 1·(1+1+1+1)/12 = 4/12 = 1/3."""
    cpsi = fw.cpsi_bell_plus(0.05, 0.07, 0.03, 0.0)
    assert abs(cpsi - 1.0/3.0) < 1e-12


def test_cpsi_bell_plus_monotonic_decay():
    """CΨ monotonically decreases with t for any nonzero noise (F26 corollary)."""
    cpsi_values = [fw.cpsi_bell_plus(0.05, 0.0, 0.05, t) for t in [0.0, 1.0, 2.0, 5.0, 10.0]]
    for i in range(len(cpsi_values) - 1):
        assert cpsi_values[i+1] < cpsi_values[i]


def test_cpsi_cusp_K_per_channel_matches_F27():
    """F27 K-values: Z=0.0374, X=Y=0.0867, depol=0.0440. Cross-check via cusp finder."""
    from scipy.optimize import brentq
    for channel, K_expected in fw.CPSI_CUSP_K_PER_CHANNEL.items():
        if channel == 'Z':
            gx, gy, gz = 0.0, 0.0, 1.0
        elif channel == 'X':
            gx, gy, gz = 1.0, 0.0, 0.0
        elif channel == 'Y':
            gx, gy, gz = 0.0, 1.0, 0.0
        elif channel == 'depolarizing':
            gx, gy, gz = 1/3, 1/3, 1/3
        # Solve CΨ(t) = 1/4 with γ = 1
        t_cusp = brentq(lambda t: fw.cpsi_bell_plus(gx, gy, gz, t) - 0.25, 1e-6, 100)
        K_computed = 1.0 * t_cusp
        assert abs(K_computed - K_expected) < 0.001, \
            f"channel {channel}: K computed {K_computed} vs F27 {K_expected}"
