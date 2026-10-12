"""DephasedSystem: the system every instrument reads (N, H, Kossakowski a)."""
import numpy as np
import pytest

import framework as fw


def _star_system(gamma, cov=None):
    chain = fw.ChainSystem(N=4, J=1.0, topology='star', H_type='heisenberg')
    return fw.DephasedSystem.from_chain(chain, gamma=gamma, covariance=cov)


def test_from_chain_defaults_to_the_chains_uniform_rate():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    s = fw.DephasedSystem.from_chain(chain)
    assert np.array_equal(s.a, np.diag([0.05, 0.05, 0.05]))
    assert s.sigma == pytest.approx(0.15)


def test_covariance_goes_off_diagonal_only():
    cov = np.full((4, 4), 0.0625)
    s = _star_system([0.5, 0.25, 0.125, 0.375], cov)
    assert np.array_equal(np.diag(s.a), [0.5, 0.25, 0.125, 0.375])
    assert s.a[0, 1] == 0.0625 and s.a[2, 3] == 0.0625


def test_refuses_an_invalid_watching():
    chain = fw.ChainSystem(N=2)
    with pytest.raises(ValueError):
        fw.DephasedSystem(chain.H, [[0.1, 0.2], [0.0, 0.1]])       # not symmetric
    with pytest.raises(ValueError):
        fw.DephasedSystem(chain.H, [[0.1, 0.5], [0.5, 0.1]])       # not PSD


def test_block_liouvillian_is_the_full_liouvillian_restricted_entry_for_entry():
    # Dyadic rates and covariances (multiples of 1/16) keep every product and sum exact, so
    # the two assemblies (Kronecker sums in lindblad.py, the price law here) must agree to
    # the bit; a nonzero residual would be a construction error, not rounding.
    cov = np.array([[0, 1, 0, 2], [1, 0, -1, 0], [0, -1, 0, 1], [2, 0, 1, 0]]) / 16
    s = _star_system(np.array([8, 4, 2, 6]) / 16, cov)
    L = s.lindbladian()
    assert s.is_number_conserving()
    pc = np.array([bin(i).count('1') for i in range(2 ** s.N)])
    for (p, q), rows in s.blocks():
        I, J = np.flatnonzero(pc == p), np.flatnonzero(pc == q)
        flat = (I[:, None] * 2 ** s.N + J[None, :]).ravel()
        assert np.array_equal(s.block_liouvillian(p, q), L[np.ix_(flat, flat)]), (p, q)


def test_blocks_carry_the_whole_liouvillian():
    s = _star_system(np.array([8, 4, 2, 6]) / 16)
    assert sum(rows for _, rows in s.blocks()) == 4 ** s.N


def test_uniform_independent_price_is_two_gamma_hamming():
    s = fw.DephasedSystem.from_chain(fw.ChainSystem(N=3, gamma_0=0.25))
    I = np.arange(8)
    z = 1 - 2 * ((I[:, None] >> (2 - np.arange(3))[None, :]) & 1)
    R = s.coherence_price(z, z)
    hamming = np.array([[bin(i ^ j).count('1') for j in I] for i in I])
    assert np.array_equal(R, 2 * 0.25 * hamming)


def test_bonds_build_the_same_hamiltonian_as_chain_system():
    # Two builders of one H: ChainSystem's dense Pauli sums and the sector hop from a bond
    # list. Every entry is 0, +-J or 2J sums of exact terms, so they must agree to the bit.
    for topo in ('chain', 'ring', 'star', 'complete'):
        chain = fw.ChainSystem(N=5, J=0.25, topology=topo, H_type='heisenberg')
        b = fw.DephasedSystem.from_bonds(5, fw.DephasedSystem.topology_bonds(5, topo), J=0.25, delta=1.0)
        assert np.array_equal(b.dense_hamiltonian(), np.real(chain.H)), topo
        # ChainSystem's 'xy' carries J/2 in front of XX + YY (its docstring), 'heisenberg' J
        chain = fw.ChainSystem(N=5, J=0.25, topology=topo, H_type='xy')
        b = fw.DephasedSystem.from_bonds(5, fw.DephasedSystem.topology_bonds(5, topo), J=0.125, delta=0.0)
        assert np.array_equal(b.dense_hamiltonian(), np.real(chain.H)), topo


def test_sectors_agree_between_the_dense_and_the_bond_form():
    chain = fw.ChainSystem(N=6, J=0.25, topology='star', H_type='heisenberg')
    dense = fw.DephasedSystem.from_chain(chain)
    bonds = fw.DephasedSystem.from_bonds(6, fw.DephasedSystem.topology_bonds(6, 'star'), J=0.25)
    for p in range(7):
        zd, Hd = dense.sector(p)
        zb, Hb = bonds.sector(p)
        assert np.array_equal(zd, zb) and np.array_equal(Hd, Hb), p


def test_a_sector_past_the_wall_is_built_without_the_dense_matrix():
    s = fw.DephasedSystem.from_bonds(100, fw.DephasedSystem.topology_bonds(100, 'chain'), J=0.25, delta=0.0)
    z, H1 = s.sector(1)
    assert H1.shape == (100, 100)
    # the open-chain single-excitation band: 2 * (2J) cos(k pi / (N + 1))
    k = np.arange(1, 101)
    assert np.allclose(np.sort(np.linalg.eigvalsh(H1)), np.sort(4 * 0.25 * np.cos(k * np.pi / 101)))


def test_block_meets_the_independent_weight_coherence_builder_entry_for_entry():
    # A second builder, written for the F89 pencil (little-endian, ascending mask, gamma
    # profile indexed by bit): map our big-endian sector states to its masks and compare
    # every entry. Dyadic rates and couplings keep both sides exact, so equality is to the bit.
    from framework.weight_coherence_block import weight_block_build, weight_block_configs
    N, J, delta = 7, 0.25, 0.5
    gamma = np.array([3, 1, 4, 1, 5, 2, 6]) / 64
    s = fw.DephasedSystem.from_bonds(N, fw.DephasedSystem.topology_bonds(N, 'chain'),
                                     J=J, delta=delta, gamma=gamma)
    for (p, q) in [(0, 1), (1, 1), (1, 2), (2, 3), (3, 3)]:
        theirs = weight_block_build(N, p, q, J, gamma=gamma, delta=delta)
        # their bit s is site s, so the profile is passed as it stands
        masks = {m: k for k, m in enumerate(weight_block_configs(N, p))}
        maskq = {m: k for k, m in enumerate(weight_block_configs(N, q))}
        to_mask = lambda z: sum(1 << l for l in range(N) if z[l] == -1)
        zp, zq = s.sector(p)[0], s.sector(q)[0]
        perm = [masks[to_mask(a)] * len(maskq) + maskq[to_mask(b)] for a in zp for b in zq]
        assert np.array_equal(s.block_liouvillian(p, q), theirs[np.ix_(perm, perm)]), (p, q)


def test_complement_carries_block_over_entry_for_entry_without_field():
    # X^N maps the coherence |i><j| to |~i><~j|: popcounts (p, q) -> (N-p, N-q), H unchanged
    # when [H, X^N] = 0, the price even in Delta. Dyadic inputs: exact, not approximate.
    N = 6
    g = np.array([0.0625, 0.03125, 0.125, 0.09375, 0.046875, 0.078125])
    cov = np.zeros((N, N))
    cov[1, 4] = cov[4, 1] = 0.015625
    cov[0, 2] = cov[2, 0] = -0.0234375
    for topo in ('chain', 'star'):
        for h, symmetric in ((None, True), (0.125, False)):
            s = fw.DephasedSystem.from_bonds(N, fw.DephasedSystem.topology_bonds(N, topo), J=0.25,
                                             delta=0.5, h=h, gamma=g, covariance=cov)
            assert s.has_xn_symmetry() == symmetric
            assert fw.DephasedSystem(s.dense_hamiltonian(), s.a).has_xn_symmetry() == symmetric
            worst = 0.0
            for p in range(N + 1):
                for q in range(N + 1):
                    zp, zq = s.sector(p)[0], s.sector(q)[0]
                    ip = {tuple(r): k for k, r in enumerate(s.sector(N - p)[0])}
                    iq = {tuple(r): k for k, r in enumerate(s.sector(N - q)[0])}
                    perm = np.array([ip[tuple(-a)] * len(zq) + iq[tuple(-b)] for a in zp for b in zq])
                    B = s.block_liouvillian(N - p, N - q)
                    worst = max(worst, np.abs(s.block_liouvillian(p, q) - B[np.ix_(perm, perm)]).max())
            assert (worst == 0.0) == symmetric, (topo, h, worst)
