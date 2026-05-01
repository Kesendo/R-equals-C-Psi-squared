"""Tests for bridge_dynamics: per-site Bloch trajectories + polarity crossings.

The bridge is geometrically the closed parametric curve on each site's
Bloch ball indexed by t. These primitives expose that curve and locate
the moments where it crosses the polarity boundary.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_bloch_trajectory_shape_and_initial_state():
    """Output shape (N, len(t_grid), 3); at t=0 matches polarity_diagnostic."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.01)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 1.0, 50)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    assert traj.shape == (N, 50, 3)
    # At t=0: |+⟩^N has X=+1, Y=Z=0 per site
    assert np.allclose(traj[:, 0, 0], 1.0, atol=1e-10)
    assert np.allclose(traj[:, 0, 1], 0.0, atol=1e-10)
    assert np.allclose(traj[:, 0, 2], 0.0, atol=1e-10)


def test_bloch_trajectory_y_zero_for_real_amplitude_state():
    """Real H + real initial state ⟹ ρ(t) real ⟹ ⟨Y_i⟩ = 0 for all t.

    Tom 2026-05-01: "Wenn kann es nur Z sein" — Y requires complex
    amplitudes, which real-amplitude polarity/bonding modes don't have.
    """
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.01)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 5.0, 100)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    assert np.allclose(traj[:, :, 1], 0.0, atol=1e-10), \
        f"Y component non-zero (max {np.max(np.abs(traj[:, :, 1])):.2e})"


def test_bloch_trajectory_polarity_state_oscillates_through_boundary():
    """|+⟩^N under XY+Z-deph oscillates: ⟨X_i⟩ does cross 0 for at least one site."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.005)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 5.0, 200)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    # Some site must have changed sign at some point
    x_min_per_site = traj[:, :, 0].min(axis=1)
    assert (x_min_per_site < 0).any(), \
        f"no site crossed polarity boundary; min ⟨X⟩ per site: {x_min_per_site}"


def test_bloch_trajectory_static_initial_state_no_crossings():
    """If chain.H = 0 (J=0), |+⟩^N is stationary up to dephasing → no X crossings.

    Z-dephasing damps |+⟩ toward I/2 (X→0 monotonically) but doesn't flip sign.
    """
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05, J=0.0)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 10.0, 100)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    events = fw.polarity_crossings(traj, t_grid)
    assert len(events) == 0, \
        f"unexpected crossings under J=0: {events}"


def test_polarity_crossings_returns_events():
    """polarity_crossings finds the boundary events for an oscillating state."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.005)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 5.0, 200)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    events = fw.polarity_crossings(traj, t_grid)
    assert len(events) > 0
    for ev in events:
        assert 0 <= ev['site'] < N
        assert ev['t_cross'] >= 0
        assert ev['direction'] in ('+→−', '−→+')
        assert abs(ev['bloch_at_cross'][0]) < 1e-3, \
            "bloch_at_cross[0] should be ~0 at the polarity boundary"


def test_polarity_crossings_y_axis():
    """axis_index=1 reports Y-axis crossings (different boundary on the Bloch ball)."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.005)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 2.0, 100)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    # For real H + real ρ, Y is identically 0; no Y crossings should fire.
    events_y = fw.polarity_crossings(traj, t_grid, axis_index=1)
    assert len(events_y) == 0


def test_polarity_crossings_invalid_axis_raises():
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    psi = fw.polarity_state(N, +1)
    t_grid = np.linspace(0, 1.0, 10)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    with pytest.raises(ValueError):
        fw.polarity_crossings(traj, t_grid, axis_index=5)


def test_bridge_reflection_signature_synthetic_rotation():
    """Synthetic +→−, −→+ pair with explicit YZ vectors gives the
    expected rotation angle.
    """
    events = [
        {'site': 0, 't_cross': 1.0, 'bloch_at_cross': (0.0, 0.6, 0.8),
         'direction': '+→−'},
        {'site': 0, 't_cross': 2.0, 'bloch_at_cross': (0.0, 0.8, 0.6),
         'direction': '−→+'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert len(sig['cycles']) == 1
    cyc = sig['cycles'][0]
    assert cyc['site'] == 0
    assert cyc['cycle_duration'] == 1.0
    # angle = arccos((0.6·0.8 + 0.8·0.6) / (1 · 1)) = arccos(0.96)
    expected = float(np.degrees(np.arccos(0.96)))
    assert abs(cyc['rotation_angle_deg'] - expected) < 1e-9
    assert sig['incomplete_sites'] == []


def test_bridge_reflection_signature_90_degree_rotation():
    """YZ_in = (1, 0), YZ_out = (0, 1) gives the canonical 90° signature."""
    events = [
        {'site': 0, 't_cross': 0.5, 'bloch_at_cross': (0.0, 1.0, 0.0),
         'direction': '+→−'},
        {'site': 0, 't_cross': 1.5, 'bloch_at_cross': (0.0, 0.0, 1.0),
         'direction': '−→+'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert abs(sig['cycles'][0]['rotation_angle_deg'] - 90.0) < 1e-9


def test_bridge_reflection_signature_180_degree_rotation():
    """YZ_in = (1, 0), YZ_out = (-1, 0) gives 180° (direct return)."""
    events = [
        {'site': 0, 't_cross': 0.5, 'bloch_at_cross': (0.0, 1.0, 0.0),
         'direction': '+→−'},
        {'site': 0, 't_cross': 1.5, 'bloch_at_cross': (0.0, -1.0, 0.0),
         'direction': '−→+'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert abs(sig['cycles'][0]['rotation_angle_deg'] - 180.0) < 1e-9


def test_bridge_reflection_signature_degenerate_yz_origin():
    """When either YZ vector is at the origin, rotation is undefined (None)."""
    events = [
        {'site': 0, 't_cross': 1.0, 'bloch_at_cross': (0.0, 0.0, 0.0),
         'direction': '+→−'},
        {'site': 0, 't_cross': 2.0, 'bloch_at_cross': (0.0, 0.5, 0.5),
         'direction': '−→+'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert sig['cycles'][0]['rotation_angle_deg'] is None
    assert sig['cycles'][0]['rotation_angle_rad'] is None


def test_bridge_reflection_signature_incomplete_cycle_flagged():
    """A +→− entry with no following −→+ return is flagged as incomplete."""
    events = [
        {'site': 0, 't_cross': 1.0, 'bloch_at_cross': (0.0, 0.5, 0.5),
         'direction': '+→−'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert len(sig['cycles']) == 0
    assert (0, '+→−') in sig['incomplete_sites']


def test_bridge_reflection_signature_per_site_independence():
    """Events from different sites are paired independently."""
    events = [
        {'site': 0, 't_cross': 1.0, 'bloch_at_cross': (0.0, 1.0, 0.0),
         'direction': '+→−'},
        {'site': 1, 't_cross': 1.5, 'bloch_at_cross': (0.0, 0.0, 1.0),
         'direction': '+→−'},
        {'site': 0, 't_cross': 2.0, 'bloch_at_cross': (0.0, 0.0, 1.0),
         'direction': '−→+'},
        {'site': 1, 't_cross': 2.5, 'bloch_at_cross': (0.0, 1.0, 0.0),
         'direction': '−→+'},
    ]
    sig = fw.bridge_reflection_signature(events)
    assert len(sig['cycles']) == 2
    sites = {c['site'] for c in sig['cycles']}
    assert sites == {0, 1}
    # Both sites have a 90° rotation
    for cyc in sig['cycles']:
        assert abs(cyc['rotation_angle_deg'] - 90.0) < 1e-9


def test_bridge_reflection_signature_real_trajectory_runs():
    """Integration: real chain trajectory → events → signature without error.

    For a permutation-symmetric initial state the YZ at crossings is at
    the origin (Tom-degenerate); the function returns None angle, which
    is the correct flagging.
    """
    chain = fw.ChainSystem(N=3, H_type='xy', gamma_0=0.005)
    psi = fw.polarity_state(3, +1)
    t_grid = np.linspace(0, 5.0, 200)
    traj = fw.bloch_trajectory(chain, psi, t_grid)
    events = fw.polarity_crossings(traj, t_grid)
    sig = fw.bridge_reflection_signature(events)
    # function runs, returns the right shape
    assert isinstance(sig, dict)
    assert 'cycles' in sig
    assert 'incomplete_sites' in sig
    # |+⟩^N is permutation-symmetric → YZ-at-crossing is degenerate
    for cyc in sig['cycles']:
        assert cyc['rotation_angle_deg'] is None


def test_bloch_trajectory_accepts_density_matrix_input():
    """Both |ψ⟩ and ρ inputs supported."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.01)
    psi = fw.polarity_state(N, +1)
    rho = np.outer(psi, psi.conj())
    t_grid = np.linspace(0, 1.0, 20)
    traj_psi = fw.bloch_trajectory(chain, psi, t_grid)
    traj_rho = fw.bloch_trajectory(chain, rho, t_grid)
    assert np.allclose(traj_psi, traj_rho, atol=1e-12)
