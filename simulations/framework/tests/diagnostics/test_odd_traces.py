"""The odd centred power sums modulo primes (framework.diagnostics.odd_traces), the Python twin of the sextant's
OddTraceObstruction: the same verdicts as compute/MirrorWorld.Tests/SextantTests.cs on the same settings, and the
seventh against the closed form of experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import odd_traces as ot  # noqa: E402


def chain(n, bond, jumps):
    """bond: {letters: coefficient} on every edge of the open chain; jumps: [(site, matrix, rate)]."""
    terms = [(lab, (b, b + 1), c) for b in range(n - 1) for lab, c in bond.items()]
    return ct.Setting(n, terms, [((l,), A, r) for l, A, r in jumps])


def test_the_warm_dna_model_is_obstructed_at_the_third_power():
    # SextantTests: X fields, ZZ, loss 3, gain 1 and Z at 1 on every site, N = 3: the third trace obstructs
    n = 3
    s = chain(n, {'ZZ': 2}, [j for l in range(n) for j in ((l, ct.LOWER, 3), (l, ct.RAISE, 1), (l, ct.PAULI['Z'], 1))])
    s.terms += [('X', (l,), -1) for l in range(n)]
    assert ot.odd_trace_obstruction(s, 9)[0] == 3


def test_dm_bonds_under_loss_have_every_odd_trace_zero():
    # F137 pairs them; every odd trace up to the dimension vanishes modulo both primes
    s = chain(3, {'XY': 1, 'YX': -1}, [(l, ct.LOWER, 1) for l in range(3)])
    assert ot.odd_trace_obstruction(s) is None


def test_the_cooled_chain_breaks_at_the_seventh_and_reads_the_closed_form():
    # XX + 2XY + ZZ under loss: Δ = 1, S = 1, A = 2, W = 5; Tr((2L + c)^7) = -7·2^17·Δ·S·A·W·(γ1γ2 - γ2γ3) at N = 3,
    # and M = L + c/2, so Tr(M^7) is that over 2^7
    bond = {'XX': 1, 'XY': 2, 'ZZ': 1}
    broken = chain(3, bond, [(l, ct.LOWER, g) for l, g in enumerate((1, 2, 3))])
    assert ot.odd_trace_obstruction(broken)[0] == 7
    law = sp.Rational(-7 * 2 ** 17 * 1 * 1 * 2 * 5 * (1 * 2 - 2 * 3), 2 ** 7)
    for p in ot.PRIMES:
        assert ot.odd_traces_mod_p(broken, [7], p)[7] == int(law.p) * pow(int(law.q), -1, p) % p
    kept = chain(3, bond, [(l, ct.LOWER, g) for l, g in enumerate((1, 2, 1))])
    assert ot.odd_trace_obstruction(kept) is None


def test_the_low_odd_traces_vanish_and_a_bond_with_a_mixed_letter_breaks_at_the_fifth():
    # the class: p1 = p3 = p5 = 0 at every profile; a mixed letter (XZ) breaks earlier
    s = chain(3, {'XX': 2, 'XY': -1, 'YX': 3, 'YY': 1, 'ZZ': -2}, [(l, ct.LOWER, g) for l, g in enumerate((2, 1, 3))])
    t = ot.odd_traces_mod_p(s, [1, 3, 5, 7], ot.PRIMES[0])
    assert t[1] == t[3] == t[5] == 0 and t[7] != 0
    mixed = chain(2, {'XZ': 1}, [(l, ct.LOWER, g) for l, g in enumerate((1, 2))])
    assert ot.odd_trace_obstruction(mixed, 9)[0] == 5
