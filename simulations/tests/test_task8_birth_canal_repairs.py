"""Small controls for the Task-8 birth-canal producer repairs."""

import numpy as np

from simulations import gamma_profile_shortcut as shortcut
from simulations import light_content


def test_tolerance_cluster_separates_edge_from_projector_mean():
    # H=0 and nearly equal site rates put two distinct real rates in the default tolerance cluster.
    profile = [0.5, 0.50000025]
    h_zero = np.zeros((4, 4), dtype=complex)

    edge, mean, basis, cluster_size = light_content.slow_subspace(2, 1.0, profile, h_zero)
    per_site, _ = light_content.projector_light(basis, 2)
    absorption = 2.0 * np.dot(profile, per_site)

    assert edge == 1.0
    assert np.isclose(mean, 1.00000025, atol=1e-12)
    assert mean != edge
    assert cluster_size == 8
    assert np.isclose(absorption, mean, atol=1e-12)


def test_shortcut_diagnostics_are_invariant_under_subspace_basis_rotation():
    rng = np.random.default_rng(404)
    raw = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
    basis, _ = np.linalg.qr(raw)
    rotation, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    rotated = basis @ rotation
    h = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)

    assert np.isclose(shortcut.subspace_overlap(basis, rotated), 1.0, atol=1e-12)
    assert np.isclose(
        shortcut.commutator_residual(basis, h, 1),
        shortcut.commutator_residual(rotated, h, 1),
        atol=1e-12,
    )
