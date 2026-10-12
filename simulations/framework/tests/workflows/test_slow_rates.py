"""slow_rates: instruments reading one DephasedSystem, checked against each other."""
import numpy as np

import framework as fw
from framework.workflows.slow_rates import slow_rates, ExactBlocks, SecularRoute, SparseBlocks


def _quiet_spokes_star(scale=1.0):
    star = fw.ChainSystem(N=5, J=0.25, topology='star', H_type='heisenberg')
    cov = np.zeros((5, 5))
    cov[1, 2] = cov[2, 1] = 0.015
    g = np.array([0.08, 0.02, 0.03, 0.04, 0.05])
    return fw.DephasedSystem.from_chain(star, gamma=scale * g, covariance=scale * cov)


def test_exact_blocks_meet_the_full_liouvillian():
    # Two routes to the slowest rate: the blocks of the workflow and one 4^N eigensolve.
    s = _quiet_spokes_star()
    full = -np.linalg.eigvals(s.lindbladian()).real
    full_slowest = full[full > 1e-9].min()
    got = slow_rates(s, [ExactBlocks()]).slowest()
    assert abs(got.rate - full_slowest) < 1e-10
    assert got.block in [(1, 1), (4, 4)]


def test_translation_invariant_ring_has_uniform_site_weights():
    ring = fw.ChainSystem(N=5, J=0.25, topology='ring', H_type='heisenberg', gamma_0=0.05)
    res = slow_rates(fw.DephasedSystem.from_chain(ring))
    for block, rows, by in res.rows:
        for r in by.values():
            if hasattr(r, 'weights'):
                assert np.ptp(r.weights) < 1e-9, (block, r.instrument, r.weights)


def test_route_error_is_a_cubic_law_not_a_threshold():
    # The secular route is first order in the price; its second order is purely imaginary,
    # so the rate error is third order. Gate the LAW: (route - exact) / s^3 constant when
    # every rate and covariance is scaled by s, over three halvings.
    devs = []
    for s in (1.0, 0.5, 0.25):
        sys_ = _quiet_spokes_star(s)
        ex = ExactBlocks().slowest(sys_, (1, 1)).rate
        ro = SecularRoute().slowest(sys_, (1, 1)).rate
        devs.append((ro - ex) / s ** 3)
    assert np.ptp(devs) / abs(np.mean(devs)) < 0.05, devs


def _profiled_chain(N=7):
    g = np.round(np.random.default_rng(2026).uniform(0.02, 0.08, N), 3)
    return fw.DephasedSystem.from_bonds(N, fw.DephasedSystem.topology_bonds(N, 'chain'), J=0.25, gamma=g)


def test_sparse_propagator_arnoldi_meets_exact_where_lr_missed():
    # ARPACK which='LR' on L returned 0.35478 on this block, exact 0.30933: Krylov favours the
    # hull of the spectrum and the slowest mode sits inside it. Arnoldi on e^{Lt} targets the
    # modulus e^{-rate t} itself; the Ritz residual is the error certificate it reports.
    s = _profiled_chain()
    ex = ExactBlocks().slowest(s, (2, 5))
    sx = SparseBlocks().slowest(s, (2, 5))
    res = float(sx.note.split()[2].rstrip(';'))
    assert res < 1e-12, sx.note
    assert abs(sx.rate - ex.rate) < 1e-9 * ex.rate, (sx.rate, ex.rate)
    assert abs(sx.omega - ex.omega) < 1e-9 * max(1.0, ex.omega)
    assert np.abs(sx.weights - ex.weights).max() < 1e-9


def test_route_refuses_a_cluster_it_cannot_hold():
    # The N = 30 star puts 143 699 of the 189 225 (2, 2) coherences in one Bohr cluster: the
    # route must refuse with a reason instead of asking for hundreds of GB.
    s = fw.DephasedSystem.from_bonds(30, fw.DephasedSystem.topology_bonds(30, 'star'), J=0.25)
    ok, why = SecularRoute().covers(s, (2, 2), 435 * 435)
    assert not ok and why.startswith("cluster 143699"), why
