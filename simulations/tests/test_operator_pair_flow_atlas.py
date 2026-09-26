"""Independent checks for the finite operator-pair/readout atlas pilot."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from framework.chain_system import ChainSystem
from framework.lindblad import lindbladian_z_dephasing
from framework.pauli import site_op
from operator_pair_flow_atlas import (
    exact_heisenberg_se_h,
    exact_readout_ode,
    exact_readout_rank,
    one_excitation_basis,
    single_excitation_flow,
)


def _complex_se_density(n: int) -> np.ndarray:
    states = one_excitation_basis(n)
    amplitudes = np.array([1.0, 1.0j, 0.5 + 0.25j] + [0.0] * (n - 3))
    amplitudes /= np.linalg.norm(amplitudes)
    psi = np.zeros(1 << n, dtype=complex)
    psi[list(states)] = amplitudes
    rho = np.outer(psi, psi.conj())
    assert not np.allclose(rho, rho.T)
    return rho


def _direct_rhs(h: np.ndarray, gamma: list[float], rho: np.ndarray) -> np.ndarray:
    """State-level equation, without any Liouvillian or pair-flow builder."""
    n = len(gamma)
    d = 1 << n
    rhs = -1j * (h @ rho - rho @ h)
    for site, rate in enumerate(gamma):
        bit = 1 << (n - 1 - site)
        signs = np.array([1 if (state & bit) == 0 else -1 for state in range(d)])
        rhs += rate * (signs[:, None] * rho * signs[None, :] - rho)
    return rhs


def test_one_excitation_basis_is_big_endian():
    assert one_excitation_basis(3) == (4, 2, 1)
    assert one_excitation_basis(5) == (16, 8, 4, 2, 1)


@pytest.mark.parametrize("n", [3, 4, 5])
def test_pair_flow_matches_full_generator_and_direct_rhs(n: int):
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    gamma = [0.13 + 0.07 * site for site in range(n)]
    flow = single_excitation_flow(h, gamma)
    full = lindbladian_z_dephasing(h, gamma)
    flat = list(flow.flat_indices)
    np.testing.assert_allclose(flow.generator, full[np.ix_(flat, flat)], atol=1e-12)
    assert np.max(np.abs(full[np.ix_([i for i in range(full.shape[0]) if i not in flat], flat)])) < 1e-12

    rho = _complex_se_density(n)
    direct = _direct_rhs(h, gamma, rho)
    sector_rho = rho[np.ix_(flow.states, flow.states)]
    sector_rhs = direct[np.ix_(flow.states, flow.states)]
    np.testing.assert_allclose(
        (flow.generator @ sector_rho.flatten()).reshape(n, n),
        sector_rhs,
        atol=1e-12,
    )


def test_complex_observable_signal_matches_full_state_at_finite_time():
    n = 3
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    gamma = [0.0, 0.7, 0.2]
    flow = single_excitation_flow(h, gamma)
    rho = _complex_se_density(n)
    observable = np.zeros_like(rho)
    a, b, c = flow.states
    observable[a, b] = 1j
    observable[b, a] = -1j
    observable[c, c] = 0.5
    t = 0.37

    full = lindbladian_z_dephasing(h, gamma)
    full_rho = (expm(t * full) @ rho.flatten()).reshape(1 << n, 1 << n)
    sector_rho = rho[np.ix_(flow.states, flow.states)]
    sector_obs = observable[np.ix_(flow.states, flow.states)]
    pair_rho = (expm(t * flow.generator) @ sector_rho.flatten()).reshape(n, n)
    assert np.trace(observable @ full_rho) == pytest.approx(
        np.trace(sector_obs @ pair_rho), abs=1e-11
    )


def test_local_x_field_rejects_false_single_excitation_closure():
    h = ChainSystem(3, J=1.0, H_type="heisenberg").H
    h = h + 0.3 * site_op(3, 1, "X")
    with pytest.raises(ValueError, match="not invariant"):
        single_excitation_flow(h, [0.1, 0.2, 0.3])


def test_exact_readout_hankel_rank_depends_on_preparation_and_observable():
    h = exact_heisenberg_se_h(3, 1)
    assert h == sp.Matrix([[0, 2, 0], [2, -2, 2], [0, 2, 0]])
    end = exact_readout_rank(h, [0, 1, 0], preparation=0, readout=0)
    centre = exact_readout_rank(h, [0, 1, 0], preparation=1, readout=1)
    trace = exact_readout_rank(h, [0, 1, 0], preparation=0, readout="trace")
    assert end.dimension == centre.dimension == trace.dimension == 9
    assert (end.rank, centre.rank, trace.rank) == (8, 4, 1)
    assert end.moments[:3] == (1, 0, -8)
    assert centre.moments[:3] == (1, 0, -16)
    assert trace.moments[:5] == (1, 0, 0, 0, 0)


def test_exact_moments_agree_with_independent_complex_pair_generator():
    h = ChainSystem(3, J=1.0, H_type="heisenberg").H
    flow = single_excitation_flow(h, [0, 1, 0])
    exact_h = exact_heisenberg_se_h(3, 1)
    for site in (0, 1):
        expected = exact_readout_rank(
            exact_h, [0, 1, 0], preparation=site, readout=site
        ).moments[:7]
        initial = np.zeros(9, dtype=complex)
        initial[site * 3 + site] = 1
        measurement = np.zeros(9, dtype=complex)
        measurement[site * 3 + site] = 1
        state = initial
        observed = []
        for _ in expected:
            observed.append(measurement @ state)
            state = flow.generator @ state
        np.testing.assert_allclose(observed, expected, atol=1e-9)


def test_centre_signal_has_exact_fourth_order_ode_and_same_time_curve():
    h = exact_heisenberg_se_h(3, 1)
    centre = exact_readout_rank(h, [0, 1, 0], preparation=1, readout=1)
    coefficients = exact_readout_ode(centre)
    assert coefficients == (0, 64, 40, 4)

    r = centre.rank
    companion = np.zeros((r, r))
    companion[np.arange(r - 1), np.arange(1, r)] = 1
    companion[-1, :] = -np.array(coefficients, dtype=float)
    derivatives_at_zero = np.array(centre.moments[:r], dtype=float)

    chain_h = ChainSystem(3, J=1.0, H_type="heisenberg").H
    flow = single_excitation_flow(chain_h, [0, 1, 0])
    initial = np.zeros(9, dtype=complex)
    initial[4] = 1  # centre population, row-major (a,b)=(1,1)
    for t in (0.1, 0.4, 1.2):
        reduced = (expm(t * companion) @ derivatives_at_zero)[0]
        full_pair = (expm(t * flow.generator) @ initial)[4]
        assert reduced == pytest.approx(full_pair.real, abs=1e-10)
        assert abs(full_pair.imag) < 1e-12


def test_exact_readout_rejects_float_input():
    h = exact_heisenberg_se_h(3, 1)
    with pytest.raises(ValueError, match="integer"):
        exact_readout_rank(h, [0, 0.3, 0], preparation=1, readout=1)
