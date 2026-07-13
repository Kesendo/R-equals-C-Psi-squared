"""Parity pin: the four frozen (1,2)-pencil transcriptions == the framework canonical.

The F89 de-monolith found the joint-popcount pencil L(q) = A + qC independently transcribed
four times in committed engines of record. Those engines stay byte-untouched (their value is
stability; committed proofs and Reproduce lines import them). This gate is the consolidation:
it pins every transcription ENTRY-EXACT against the framework mirror
`framework/weight_coherence_block.py` of the C# canonical
`compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs`, through each script's documented
adapter. The adapter table below IS the documentation of the four conventions; the gate stays
discriminating forever because the legacy builders remain independent transcriptions.

Adapter table (legs = (w_ket, w_bra); order = basis enumeration; sign/object = what the script
returns in terms of the module's A (AT diagonal, gamma=1) and K (real unit-hop difference)):

| script.builder                            | legs  | order                       | object returned              | q-book  |
|-------------------------------------------|-------|-----------------------------|------------------------------|---------|
| seed_existence_nullity_check.build(N)     | (2,1) | combinations x combinations | (A, C) with C = -1j*K        | carrier |
| resonant_n_twinning.blocks(N)             | (2,1) | combinations x combinations | K66/K26 sub-blocks of MINUS K, t6 | (none; K only) |
| cross_triple_orthogonality.block_space(N) | (2,1) | combinations x plain-int bra| H = hop SUM, K = MINUS K, J/nu/six | (none; K, H only) |
| f89_why_diabolic_probe.build_L(J, gamma)  | (1,2) | ascending SE x DE_PAIRS     | full L (module sign), octic q = J | octic |

Float caveat (the one non-exact row): f89's diagonal is the two-term program -6*gamma + 4*gamma*pov,
the module's the one-term -2*gamma*n_diff; equal in value, but different floating-point programs,
so at generic gamma (e.g. 0.3) the overlap entries differ in the last ULP. The gate is therefore
entry-exact for all hop entries at every gamma, entry-exact for the diagonal at gamma in {1, 0},
and <= 1e-15 for the diagonal at generic gamma. All other rows are entry-exact (integer-valued
constructions, exact IEEE scalings).

Import note: f89_why_diabolic_probe mutates process-global state at import
(sys.stdout.reconfigure, np.set_printoptions); the gate snapshots and restores np printoptions.

Run: python simulations/weight_coherence_parity.py   (~seconds; all gates must print PASS)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import framework as fw
import seed_existence_nullity_check as seed
import resonant_n_twinning as twin
import cross_triple_orthogonality as ct

_printoptions = np.get_printoptions()
import f89_why_diabolic_probe as f89  # mutates np printoptions + stdout encoding at import
np.set_printoptions(**_printoptions)


def _index_map_21(n):
    """Legacy (2,1) combinations-order index -> module ascending-mask index."""
    perm2 = fw.combinations_order_permutation(n, 2)
    d1 = n  # C(n,1); w=1 orders coincide (identity permutation)
    return np.array([p2 * d1 + i1 for p2 in perm2 for i1 in range(d1)], dtype=int)


def gate_seed(n_range=range(3, 10)):
    for n in n_range:
        a_s, c_s = seed.build(n)
        a_m, k_m = fw.weight_block_pencil(n, 2, 1, gamma=1.0)
        m = _index_map_21(n)
        assert np.array_equal(a_s, a_m[m]), f"seed A mismatch at N={n}"
        assert np.array_equal(c_s, (-1j * k_m)[np.ix_(m, m)]), f"seed C mismatch at N={n}"
    print(f"PASS seed_existence_nullity_check.build == pencil(2,1), C = -1j*K, "
          f"combinations order, N={n_range.start}..{n_range.stop - 1} (entry-exact)")


def gate_twinning(n_range=range(3, 10)):
    for n in n_range:
        k66, k26, t6 = twin.blocks(n)
        a_m, k_m = fw.weight_block_pencil(n, 2, 1, gamma=1.0)
        m = _index_map_21(n)
        k_legacy_order = (-k_m)[np.ix_(m, m)]          # the script's sign: K = MINUS(module K)
        a_legacy_order = a_m[m]
        m2 = a_legacy_order == -2.0
        m6 = a_legacy_order == -6.0
        assert np.array_equal(k66, k_legacy_order[np.ix_(m6, m6)]), f"twinning K66 mismatch at N={n}"
        assert np.array_equal(k26, k_legacy_order[np.ix_(m2, m6)]), f"twinning K26 mismatch at N={n}"
        # the bipartite sign, re-derived from module masks: t = (-1)^(sum ket sites + bra site)
        kets = fw.weight_block_configs(n, 2)
        bras = fw.weight_block_configs(n, 1)
        t_m = np.array([(-1.0) ** (sum(fw.sites_of_mask(k, n)) + fw.sites_of_mask(b, n)[0])
                        for k in kets for b in bras])
        assert np.array_equal(t6, t_m[m][m6]), f"twinning t6 mismatch at N={n}"
    print(f"PASS resonant_n_twinning.blocks == MINUS(K) sub-blocks + bipartite sign, "
          f"N={n_range.start}..{n_range.stop - 1} (entry-exact)")


def gate_cross_triple(n_range=range(3, 8)):
    for n in n_range:
        basis, h_ct, k_ct, j_ct, nu_ct, six_ct = ct.block_space(n)
        a_m, k_m = fw.weight_block_pencil(n, 2, 1, gamma=1.0)
        h_m = (np.kron(fw.weight_block_hop(n, 2), np.eye(n))
               + np.kron(np.eye(len(fw.weight_block_configs(n, 2))), fw.weight_block_hop(n, 1)))
        m = _index_map_21(n)
        assert np.array_equal(k_ct, (-k_m)[np.ix_(m, m)]), f"cross_triple K mismatch at N={n}"
        assert np.array_equal(h_ct, h_m[np.ix_(m, m)]), f"cross_triple H mismatch at N={n}"
        # basis-derived vectors, re-derived from module masks (bra = plain site int in the script)
        kets = fw.weight_block_configs(n, 2)
        bras = fw.weight_block_configs(n, 1)
        j_m = np.diag([(-1.0) ** fw.sites_of_mask(b, n)[0] for k in kets for b in bras])
        nu_m = np.array([sum(1 for x in fw.sites_of_mask(k, n) if x < fw.sites_of_mask(b, n)[0])
                         for k in kets for b in bras], float)
        six_m = np.array([not ((k >> fw.sites_of_mask(b, n)[0]) & 1) for k in kets for b in bras])
        assert np.array_equal(j_ct, j_m[np.ix_(m, m)]), f"cross_triple J mismatch at N={n}"
        assert np.array_equal(nu_ct, nu_m[m]), f"cross_triple nu mismatch at N={n}"
        assert np.array_equal(six_ct, six_m[m]), f"cross_triple six mismatch at N={n}"
    print(f"PASS cross_triple_orthogonality.block_space == (H sum, MINUS K, J, nu, six), "
          f"plain-int bra, N={n_range.start}..{n_range.stop - 1} (entry-exact)")


def gate_f89():
    n = 4
    perm_bra = fw.combinations_order_permutation(n, 2)   # DE_PAIRS -> ascending: [3,5,9,6,10,12]
    m = np.array([i * 6 + perm_bra[j] for i in range(4) for j in range(6)], dtype=int)
    off = ~np.eye(24, dtype=bool)
    for q in (1.5, 0.658983):
        for gamma in (1.0, 0.3, 0.0):
            l_f, *_ = f89.build_L(q, gamma)
            l_m = fw.weight_block_build(n, 1, 2, q, gamma=gamma, book="octic")[np.ix_(m, m)]
            assert np.array_equal(l_f[off], l_m[off]), \
                f"f89 hop mismatch at q={q}, gamma={gamma}"
            d_f, d_m = np.diag(l_f), np.diag(l_m)
            if gamma in (1.0, 0.0):
                assert np.array_equal(d_f, d_m), f"f89 diagonal mismatch at q={q}, gamma={gamma}"
            else:
                # two-term (-6g + 4g*pov) vs one-term (-2g*n_diff) float programs: last-ULP only
                assert np.abs(d_f - d_m).max() <= 1e-15, \
                    f"f89 diagonal beyond ULP at q={q}, gamma={gamma}"
    print("PASS f89_why_diabolic_probe.build_L == build(4,1,2,q,'octic'), DE_PAIRS permutation, "
          "q in {1.5, q_EP}, gamma in {1, 0.3, 0} (hops entry-exact; diagonal exact at gamma in "
          "{1,0}, <= 1e-15 at gamma=0.3)")


if __name__ == "__main__":
    gate_seed()
    gate_twinning()
    gate_cross_triple()
    gate_f89()
    print("ALL PARITY GATES PASS: the four frozen transcriptions match the one pencil "
          "(three rows entry-exact; the f89 diagonal at generic gamma to <= 1e-15).")
