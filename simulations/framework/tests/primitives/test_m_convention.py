"""The Pauli-basis transform's stacking convention, pinned.

Sibling of tests/workflows/test_vec_convention.py, which pins the Liouvillian's
row-stack convention. That one deliberately left the Pauli basis transform alone;
this is the missing half.

Everything here fails on a state with COMPLEX off-diagonals and passes on a real
symmetric one, because flatten() and flatten('F') agree exactly when ρ = ρᵀ.
Every ρ the framework builds by default (bonding modes, polarity states, Bell,
computational basis) is real, which is why a green suite proved nothing about
this for so long.
"""
import numpy as np
import pytest
from scipy.linalg import expm

import framework as fw
from framework.pauli import (
    _vec_to_pauli_basis_transform, pauli_basis_vector, _pauli_label,
)


def _complex_rho(N, seed=11):
    """Hermitian, unit-trace, with genuinely complex off-diagonals (ρ != ρᵀ)."""
    rng = np.random.default_rng(seed)
    d = 2 ** N
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = A @ A.conj().T
    return rho / np.trace(rho).real


def _y_parity_signs(N):
    return np.array([(-1) ** _pauli_label(k, N).count('Y') for k in range(4 ** N)],
                    dtype=complex)


@pytest.mark.parametrize("N", [1, 2, 3])
def test_the_two_stackings_differ_by_the_y_parity_sign(N):
    """M_F = M_C · C with C = diag((−1)^n_Y). This is the whole discrepancy."""
    M_F = np.asarray(_vec_to_pauli_basis_transform(N, order='F'))
    M_C = np.asarray(_vec_to_pauli_basis_transform(N, order='C'))
    c = _y_parity_signs(N)
    assert np.allclose(M_F, M_C * c)
    for M in (M_F, M_C):
        assert np.allclose(M.conj().T @ M, (2 ** N) * np.eye(4 ** N))


@pytest.mark.parametrize("N", [1, 2, 3])
def test_row_stack_transform_inverts_pauli_basis_vector(N):
    """pauli_basis_vector returns Tr(σρ)/2^N, which is row-stack-consistent.

    The column-stack transform reconstructs ρᵀ from those same coefficients: that
    is the trap, and on a Hermitian ρ it hands back the conjugate state.
    """
    rho = _complex_rho(N)
    a = pauli_basis_vector(rho, N)
    M_C = np.asarray(_vec_to_pauli_basis_transform(N, order='C'))
    M_F = np.asarray(_vec_to_pauli_basis_transform(N, order='F'))
    assert np.allclose(M_C @ a, rho.flatten())
    assert np.allclose(M_F @ a, rho.T.flatten())
    assert not np.allclose(M_F @ a, rho.flatten())


@pytest.mark.parametrize("N", [2, 3])
def test_pauli_basis_generator_propagates_rho_not_its_transpose(N):
    """A generator built for coefficient propagation must use the row-stack transform.

    Build it column-stack and exp(L_pauli t) transports the coefficients of ρᵀ,
    so the trajectory is that of the conjugate state. Checked against expm on the
    row-stack Liouvillian, which needs no Pauli basis at all.
    """
    rho = _complex_rho(N)
    H = fw._build_bilinear(N, [(i, i + 1) for i in range(N - 1)],
                           [('X', 'X', 1.0), ('Y', 'Y', 1.0), ('Z', 'Z', 1.0)])
    L = fw.lindbladian_z_dephasing(H, [0.05] * N)
    d = 2 ** N
    a = pauli_basis_vector(rho, N)

    for order, expect_match in (('C', True), ('F', False)):
        M = np.asarray(_vec_to_pauli_basis_transform(N, order=order))
        G = M.conj().T @ L @ M / (2 ** N)
        for t in (0.3, 1.7):
            coeffs = expm(G * t) @ a
            # Reconstruction from coefficients is convention-free: ρ = Σ a_k σ_k.
            got = sum(coeffs[k] * fw.pauli_string(list(fw._k_to_indices(k, N)))
                      for k in range(4 ** N))
            want = (expm(L * t) @ rho.flatten()).reshape(d, d)
            assert bool(np.allclose(got, want)) == expect_match, (
                f"order={order!r} at t={t}: expected match={expect_match}, "
                f"max|Δ| = {np.max(np.abs(got - want)):.3e}")


@pytest.mark.parametrize("N", [2, 3])
def test_real_symmetric_rho_cannot_see_any_of_this(N):
    """The blindness itself, pinned: on ρ = ρᵀ both stackings agree exactly.

    Guards against 'fixing' the convention by switching test inputs to real states
    and declaring the discrepancy gone.
    """
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[-1] = 1 / np.sqrt(2)
    rho = np.outer(psi, psi.conj())          # Bell-like, real symmetric
    assert np.allclose(rho, rho.T)
    a = pauli_basis_vector(rho, N)
    M_C = np.asarray(_vec_to_pauli_basis_transform(N, order='C'))
    M_F = np.asarray(_vec_to_pauli_basis_transform(N, order='F'))
    assert np.allclose(M_C @ a, M_F @ a)


# ----------------------------------------------------------------------
# The call sites. The pins above hold the primitive; these hold the four
# places that were actually wrong, so that reverting any of them goes red.
# ----------------------------------------------------------------------

def test_cockpit_trajectory_is_rho_not_its_conjugate():
    """lebensader: the CΨ trace must follow ρ, not ρᵀ.

    Oracle is expm on the row-stack Liouvillian, which needs no Pauli basis, so
    this cannot agree with the generator by construction.
    """
    from framework.lebensader import cockpit_panel, _purity_psi_norm_cpsi
    N = 2
    rho = _complex_rho(N, seed=5)
    H = fw._build_bilinear(N, [(0, 1)], [('X', 'X', 1.0), ('Y', 'Y', 1.0), ('Z', 'Z', 1.0)])
    gamma_l = [0.05] * N
    out = cockpit_panel(H, gamma_l, rho, N, t_max=2.0, dt=0.02)
    times = out['_trajectory_for_inspection']['times']
    cpsi = out['_trajectory_for_inspection']['cpsi']

    L = fw.lindbladian_z_dephasing(H, gamma_l)
    d = 2 ** N
    for t_want in (0.5, 1.0, 2.0):
        i = int(np.argmin(np.abs(times - t_want)))
        rho_t = (expm(L * times[i]) @ rho.flatten()).reshape(d, d)
        assert np.isclose(cpsi[i], _purity_psi_norm_cpsi(rho_t)[2], atol=1e-9), (
            f"CΨ at t={times[i]:.3f}: cockpit {cpsi[i]:.9f} vs expm "
            f"{_purity_psi_norm_cpsi(rho_t)[2]:.9f}")


def test_protected_observables_read_rho_not_its_conjugate():
    """observables: the cluster contributions must be those of ρ.

    Under the column-stack generator they are ρᵀ's, so the two states' answers
    are swapped; a complex ρ separates them.
    """
    from framework.observables import pi_protected_observables
    N = 2
    rho = _complex_rho(N, seed=5)
    H = fw._build_bilinear(N, [(0, 1)], [('X', 'Z', 1.0), ('Y', 'Z', 1.0)])
    gamma_l = [0.05] * N
    got = {e['k']: e['max_cluster_contribution']
           for e in pi_protected_observables(H, gamma_l, rho, N)['active']}
    swapped = {e['k']: e['max_cluster_contribution']
               for e in pi_protected_observables(H, gamma_l, rho.T, N)['active']}
    shared = set(got) & set(swapped)
    assert shared, "no active observables in common; pick a different H"
    assert any(not np.isclose(got[k], swapped[k]) for k in shared), (
        "ρ and ρᵀ gave identical contributions, so this test cannot see the "
        "convention; choose a state further from symmetric")

    # The reference: same pipeline, generator built explicitly row-stack.
    M_C = np.asarray(_vec_to_pauli_basis_transform(N, order='C'))
    L = fw.lindbladian_z_dephasing(H, gamma_l)
    G = M_C.conj().T @ L @ M_C / (2 ** N)
    evals, V = np.linalg.eig(G)
    c = np.linalg.inv(V) @ pauli_basis_vector(rho, N)
    for k, want in got.items():
        ref = max(abs(V[k, j] * c[j]) for j in range(4 ** N))
        assert np.isclose(want, ref, rtol=1e-7), f"observable {k}: {want} vs {ref}"


def test_d_zero_pauli_decomposition_matches_its_own_docstring():
    """d_zero: the returned coefficients must be Tr(σ_k · mode)/2^N, as documented.

    Uses a Y-dephasing kernel, which lives in {I, Y}^N and so carries the odd-n_Y
    strings the column-stack transform sign-flips. Basis-independent: each mode is
    checked against its own reconstruction.
    """
    from framework.diagnostics.d_zero import stationary_modes
    N = 2
    d = 2 ** N
    chain = fw.ChainSystem(N=N)
    L = fw.lindbladian_general(np.zeros((d, d), dtype=complex),
                               [np.sqrt(0.05) * fw.site_op(N, l, 'Y') for l in range(N)])
    out = stationary_modes(chain, L=L)
    assert out['kernel_dimension'] >= 2
    carries_y = False
    for ii in range(out['kernel_dimension']):
        mode = out['right_eigvecs'][:, ii].reshape(d, d)
        want = pauli_basis_vector(mode, N)
        if any(abs(want[k]) > 1e-9 and _pauli_label(k, N).count('Y') % 2 for k in range(4 ** N)):
            carries_y = True
        assert np.allclose(out['pauli_decomposition'][ii], want, atol=1e-10)
    assert carries_y, "kernel carried no odd-Y string; the test would be blind"
