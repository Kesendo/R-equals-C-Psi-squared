"""The two continuous-crossover pairs (XZ+YZ, ZX+ZY) are LOCAL, not non-local.

These tests pin the correction of the old "genuinely non-local / entangled Π" verdict:
a single uniform per-site unitary M (M² = −I), the closed-form crossover map, mirrors the
full Liouvillian to machine precision at N = 2, 3, 4, while the discrete signed-permutation
crossover P1 fails. The mirror is a product, hence local.
"""
import numpy as np
import pytest

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import framework as fw  # noqa: E402

XZ_YZ = (('X', 'Z'), ('Y', 'Z'))
ZX_ZY = (('Z', 'X'), ('Z', 'Y'))


def test_crossover_map_is_unitary_involution():
    """M is unitary with M² = −I (an order-4 element, eigenvalues ±i)."""
    M = fw.crossover_map()
    assert np.allclose(M @ M.conj().T, np.eye(4), atol=1e-12)
    assert np.allclose(M @ M, -np.eye(4), atol=1e-12)


def test_crossover_map_is_continuous_not_permutation():
    """M mixes X and Y: 8 nonzeros (a signed permutation would have 4)."""
    M = fw.crossover_map()
    assert int(np.sum(np.abs(M) > 1e-9)) == 8


@pytest.mark.parametrize("combo", [XZ_YZ, ZX_ZY])
@pytest.mark.parametrize("N", [2, 3, 4])
def test_crossover_product_pi_is_local(combo, N):
    """Π = M^⊗N mirrors the full Liouvillian to machine precision: the case is LOCAL."""
    resid = fw.verify_crossover_local(combo, N, gamma=0.5)
    assert resid < 1e-7, f"{combo} N={N}: residual {resid:.2e} should be ~0"


@pytest.mark.parametrize("combo", [XZ_YZ, ZX_ZY])
def test_discrete_crossover_fails(combo):
    """The discrete signed-permutation crossover P1^⊗N does NOT close these cases.

    This is the contrast that the old discrete-only search saw and read as 'non-local'.
    """
    # P1 in framework order [I,X,Z,Y]: I->X, X->I, Z->iY, Y->iZ.
    P1 = np.zeros((4, 4), dtype=complex)
    P1[1, 0] = 1; P1[0, 1] = 1; P1[3, 2] = 1j; P1[2, 3] = 1j
    N = 3
    resid = fw.product_pi_residual(combo, [P1] * N, N, gamma=0.5)
    assert resid > 1.0, f"discrete P1^⊗N should fail for {combo}; got {resid:.2e}"


def test_is_crossover_pair():
    assert fw.is_crossover_pair('XZ', 'YZ')
    assert fw.is_crossover_pair('ZY', 'ZX')        # order-independent
    assert fw.is_crossover_pair(('X', 'Z'), ('Y', 'Z'))
    assert not fw.is_crossover_pair('XZ', 'ZY')    # different sites -> discrete-local
    assert not fw.is_crossover_pair('XX', 'YY')    # truly


def test_verify_rejects_non_crossover():
    with pytest.raises(ValueError):
        fw.verify_crossover_local((('X', 'Z'), ('Z', 'Y')), 2)
