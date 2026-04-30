"""Tests for ChainSystem.propagate_with_hardware_noise (Marrakesh-confirmed values)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_propagate_with_hardware_noise_no_evolution():
    """t=0 returns ρ_0 (with symmetrization)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(rho_0, t=0.0)
    np.testing.assert_allclose(rho_t, rho_0, atol=1e-12)


def test_propagate_with_hardware_noise_no_dissipation_unitary():
    """No T1/Tphi/etc. → unitary evolution, trace preserved exactly."""
    chain = fw.ChainSystem(N=3, gamma_0=0.0)  # no Z-dephasing
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    assert abs(np.trace(rho_t).real - 1.0) < 1e-10
    # Purity preserved (unitary on pure state)
    assert abs(np.trace(rho_t @ rho_t).real - 1.0) < 1e-9


def test_propagate_with_hardware_noise_matches_marrakesh_xz_for_truly():
    """Marrakesh April 26: truly-unbroken (XX+YY) Hamiltonian on |+−+⟩ at t=0.8.

    With γ_Z=0.1 (idealized), <X_0 Z_2> should be near 0 (matches HW within readout noise)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')])
    # Trace out middle qubit, compute <X_0 Z_2>
    rho_3q = rho_t.reshape(2,2,2,2,2,2)
    rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    X0Z2 = np.kron(X, Z)
    exp = float(np.real(np.trace(rho_02 @ X0Z2)))
    # Idealized framework prediction: 0.000 (per Marrakesh README); HW measured +0.011.
    # Our propagation should match the framework idealized number ~0.
    assert abs(exp) < 0.01


def test_propagate_with_hardware_noise_matches_marrakesh_xz_for_soft():
    """Marrakesh: soft (XY+YX) under γ_Z=0.1 should give <X_0 Z_2> ≈ -0.62 (matches README)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')])
    rho_3q = rho_t.reshape(2,2,2,2,2,2)
    rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    exp = float(np.real(np.trace(rho_02 @ np.kron(X, Z))))
    # README idealized: -0.623; HW: -0.711. Our propagation should match idealized.
    assert abs(exp - (-0.623)) < 0.02


def test_propagate_with_hardware_noise_zz_crosstalk_changes_dynamics():
    """J_zz adds Z-Z Hamiltonian correction; non-trivial effect on Pauli expectations."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_no_zz = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')])
    rho_with_zz = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')], J_zz=0.3)

    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    Y0Z2 = np.kron(Y, Z)
    def y0z2(rho):
        rho_3q = rho.reshape(2,2,2,2,2,2)
        rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
        return float(np.real(np.trace(rho_02 @ Y0Z2)))
    e_no = y0z2(rho_no_zz)
    e_with = y0z2(rho_with_zz)
    # ZZ-crosstalk should shift <Y_0 Z_2> by > 0.1
    assert abs(e_with - e_no) > 0.1, \
        f"J_zz=0.3 should shift <Y_0 Z_2> by >0.1; got Δ={abs(e_with - e_no):.4f}"


def test_propagate_with_hardware_noise_t1_has_nontrivial_effect():
    """T1 changes <X_0 Z_2> measurably (sanity that the T1 channel actually fires)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_id = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')])
    rho_t1 = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')],
        T1_l=[0.05]*3)  # 10× larger T1 to ensure visible effect

    def x0z2(rho):
        rho_3q = rho.reshape(2,2,2,2,2,2)
        rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
        X = np.array([[0,1],[1,0]], dtype=complex)
        Z = np.array([[1,0],[0,-1]], dtype=complex)
        return float(np.real(np.trace(rho_02 @ np.kron(X, Z))))

    e_id = x0z2(rho_id)
    e_t1 = x0z2(rho_t1)
    assert abs(e_t1 - e_id) > 0.01, \
        f"T1=0.05 must shift <X_0 Z_2> by >0.01; got Δ={abs(e_t1 - e_id):.4f}"
