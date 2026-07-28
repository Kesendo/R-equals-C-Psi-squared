"""The vec convention of the Liouvillian builders, pinned on a non-symmetric state.

Every other propagation test in this suite starts from a real-symmetric rho
(Bell+, computational basis, real bonding modes), and for those flatten() and
flatten('F') return the same vector. That is why a row-stack Liouvillian paired
with column-stack callers went unnoticed: the helpers propagated the CONJUGATE
state and no test could see it. These tests use a state with complex
off-diagonals, where the two conventions part company.
"""
import numpy as np
import pytest
from scipy.linalg import expm

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from framework.chain_system import ChainSystem
from framework.lindblad import lindbladian_z_dephasing
from framework.pauli import site_paulis
from framework.workflows._propagation import (
    propagation_setup,
    per_site_bloch_trajectory,
)


def _non_symmetric_state(N):
    """A pure state whose density matrix has complex off-diagonals, rho != rho.T."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[1] = 0.5j
    psi[d - 1] = 0.5
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    assert not np.allclose(rho, rho.T), "the separating state must not be symmetric"
    return rho


def test_liouvillian_acts_on_row_stacked_vec():
    """L @ rho.flatten() is the RHS; L @ rho.flatten('F') is not."""
    N = 2
    H = ChainSystem(N=N).H
    gammas = [0.2, 0.3]
    L = lindbladian_z_dephasing(H, gammas)
    rho = _non_symmetric_state(N)

    Z = np.diag([1.0, -1.0]).astype(complex)

    def site_z(l):
        op = np.array([[1.0]], dtype=complex)
        for k in range(N):
            op = np.kron(op, Z if k == l else np.eye(2))
        return op

    rhs = -1j * (H @ rho - rho @ H)
    for l in range(N):
        Zl = site_z(l)
        rhs += gammas[l] * (Zl @ rho @ Zl - rho)

    assert np.allclose(L @ rho.flatten(), rhs.flatten())
    assert not np.allclose(L @ rho.flatten('F'), rhs.flatten('F'))


@pytest.mark.parametrize("t", [0.5, 1.5, 3.0])
def test_bloch_trajectory_matches_expm_on_a_non_symmetric_state(t):
    """The spectral helper reproduces expm(L t) to machine precision."""
    N = 3
    H = ChainSystem(N=N).H
    L = lindbladian_z_dephasing(H, [0.2, 0.3, 0.4])
    rho = _non_symmetric_state(N)
    d = 2 ** N
    paulis = site_paulis(N)

    evals, R, _, c0 = propagation_setup(L, rho)
    got = np.asarray(per_site_bloch_trajectory(evals, R, c0, np.array([t]), paulis))

    rho_t = (expm(L * t) @ rho.flatten()).reshape(d, d)
    rho_t = (rho_t + rho_t.conj().T) / 2

    for site in range(N):
        for comp in range(3):
            expected = float(np.real(np.trace(paulis[site][comp] @ rho_t)))
            assert got[site, 0, comp] == pytest.approx(expected, abs=1e-10)
