"""Tests for framework.diagnostics.cube_reflector (the Theodolite). Found answers are checked by an independent route:
the full generator, built row-major from the Lindblad formula and only then permuted to per-site order, must
satisfy E·L·E⁻¹ = 𝒯(L) exactly. NONE rows are spectra that do not pair (an exact odd trace, or F137's −656
triangle); slow rows are marked and deselected by default with -m 'not slow'."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import cube_reflector as cr  # noqa: E402

R = sp.Rational
DOWN = [R(1, 2), R(1), R(3, 2)]
UP = [R(1, 8), R(1, 3), R(1, 5)]


def thermal_jumps(down, up):
    return [((l,), ct.LOWER, d) for l, d in enumerate(down)] + [((l,), ct.RAISE, u) for l, u in enumerate(up)]


def bond(a, b, J=1, D=0, Dz=0):
    terms = []
    if J:
        terms += [('XX', (a, b), J), ('YY', (a, b), J)]
    if D:
        terms += [('XY', (a, b), D), ('YX', (a, b), -D)]
    if Dz:
        terms += [('ZZ', (a, b), Dz)]
    return terms


def full_generator(setting):
    """Row-major Lindblad generator on N sites, then permuted to per-site order (an independent route)."""
    n = setting.n
    d = 2 ** n

    def on(op_by_site):
        return sp.kronecker_product(*[op_by_site.get(l, sp.eye(2)) for l in range(n)])

    H = sp.zeros(d, d)
    for letters, sites, coef in setting.terms:
        H += coef * on({s: ct.PAULI[c] for s, c in zip(sites, letters)})
    Id = sp.eye(d)
    L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for sites, A, rate in setting.jumps:
        J = on({sites[0]: sp.Matrix(A)})
        JdJ = J.H * J
        L += rate * (sp.kronecker_product(J, J.conjugate())
                     - R(1, 2) * (sp.kronecker_product(JdJ, Id) + sp.kronecker_product(Id, JdJ.T)))
    P = ct._per_site_perm(n)
    return (P * L * P.T).applyfunc(sp.expand)


def handshake(setting, relation, point):
    E = sp.kronecker_product(*[point[l] for l in range(setting.n)])
    L = full_generator(setting)
    c = cr.constant_c(setting)
    rhs = (cr.target(L, relation) - c * sp.eye(L.shape[0])) * E
    return (E * L - rhs).applyfunc(sp.expand).is_zero_matrix


def heisenberg_gradient():
    return ct.Setting(3, bond(0, 1, 1, 0, 1) + bond(1, 2, 1, 0, 1), thermal_jumps(DOWN, UP))


def test_site_family_dimensions():
    thermal = ct.site_table(ct.Setting(1, [], thermal_jumps([R(1, 2)], [R(1, 8)])), 0)
    zdeph = ct.site_table(ct.Setting(1, [], [((0,), ct.PAULI['Z'], 1)]), 0)
    empty = ct.site_table(ct.Setting(1), 0)
    for rel in cr.RELATIONS:
        assert len(cr.site_family(thermal, rel)) == 6
        assert len(cr.site_family(zdeph, rel)) == 8
        assert len(cr.site_family(empty, rel)) == 16


def test_forty_families_and_the_candidates_per_site():
    assert len(cr._families()) == 40                     # 24 weight-after + 24 weight-before, 8 coincide
    thermal = ct.site_table(ct.Setting(1, [], thermal_jumps([R(1, 2)], [R(1, 8)])), 0)
    zdeph = ct.site_table(ct.Setting(1, [], [((0,), ct.PAULI['Z'], 1)]), 0)
    t1 = ct.site_table(ct.Setting(1, [], [((0,), ct.LOWER, 1)]), 0)
    assert [len(cr.cube_candidates(thermal, r)) for r in cr.RELATIONS] == [2, 2]
    assert [len(cr.cube_candidates(zdeph, r)) for r in cr.RELATIONS] == [4, 4]
    assert len(cr.cube_candidates(t1, 'pi')) == 0            # cooling alone: no cube move for the Pi-reflection


def test_precheck_on_the_gradient_chain():
    s = heisenberg_gradient()
    assert cr.precheck(s, 'refl')[0] == 'PASS'          # E^-1(I) ∝ Z^⊗N commutes with every number-conserving H
    assert cr.precheck(s, 'pi')[0] == 'NONE'            # the identity row ⊗(Z + κ) fails across a gradient
    one_ratio = ct.Setting(3, s.terms, thermal_jumps(DOWN, [d / 4 for d in DOWN]))
    assert cr.precheck(one_ratio, 'pi')[0] == 'PASS'


def test_answer_on_the_gradient_chain_finds_the_reflection_and_checks_it():
    s = heisenberg_gradient()
    v = cr.answer(s, 'refl', timeout=None)
    assert v.kind == 'CUBE MOVES' and v.c == sum(DOWN) + sum(UP)
    assert v.point is not None and handshake(s, 'refl', v.point)
    w = cr.answer(s, 'pi', timeout=None)
    assert w.kind == 'NONE' and w.note.startswith('pre-check')


def test_dm_chain_reflection_is_a_cube_move():
    s = ct.Setting(3, bond(0, 1, 4, 3) + bond(1, 2, 4, 3), thermal_jumps(DOWN, UP))
    v = cr.answer(s, 'refl', timeout=None)
    assert v.kind == 'CUBE MOVES'
    if v.point is not None:
        assert handshake(s, 'refl', v.point)


def test_m_beta_on_the_flux_triangle_at_one_ratio():
    down = [R(1), R(2), R(3)]
    terms = bond(0, 1) + bond(1, 2) + bond(2, 0, 0, 1)            # flux pi/2: J + iD = i on bond (2, 0)
    s = ct.Setting(3, terms, thermal_jumps(down, [d / 4 for d in down]))
    v = cr.answer(s, 'pi', timeout=None)
    assert v.kind == 'CUBE MOVES'
    if v.point is not None:
        assert handshake(s, 'pi', v.point)


def test_the_minus_656_triangle_has_none():
    terms = bond(0, 1) + bond(1, 2) + bond(2, 0, 0, 1)
    s = ct.Setting(3, terms, thermal_jumps(DOWN, UP))
    for rel in cr.RELATIONS:
        v = cr.answer(s, rel, timeout=None)
        assert v.kind == 'NONE'


def test_the_continuous_rotation_rows_are_cube_moves():
    for t1, t2 in (('XZ', 'YZ'), ('ZX', 'ZY')):
        terms = [(t, (a, a + 1), 1) for a in range(2) for t in (t1, t2)]
        s = ct.Setting(3, terms, [((l,), ct.PAULI['Z'], 1) for l in range(3)])
        assert cr.cube_search(s, 'pi')[0] == 'FOUND'


def test_golden_router_is_found_as_cube_moves():
    terms = [(lab, (w, w + 1, w + 2), c) for w in range(2) for lab, c in (('XZX', 3), ('XZY', 2), ('YZX', 2))]
    s = ct.Setting(4, terms, [((l,), ct.PAULI['Z'], R(1, 1 + l)) for l in range(4)])
    status, hit = cr.cube_search(s, 'pi')
    assert status == 'FOUND'
    assert all(name == hit['names'][0] for name in hit['names'].values())   # one move on every site, X·rho


def test_every_n_kinds():
    heis = cr.Cell(jumps=[(ct.LOWER, 1), (ct.RAISE, R(1, 4))],
                   windows=[('XX', (0, 1), 1), ('YY', (0, 1), 1), ('ZZ', (0, 1), 1)])
    assert cr.every_n(heis, 'refl', rate_free=True) == 'UNIFORM'
    dm = cr.Cell(jumps=[(ct.LOWER, 1), (ct.RAISE, R(1, 4))],
                 windows=[('XX', (0, 1), 4), ('YY', (0, 1), 4), ('XY', (0, 1), 3), ('YX', (0, 1), -3)])
    assert cr.every_n(dm, 'refl', rate_free=True) == 'TURNED'
    for t2 in ('XY', 'XZ'):
        lead = 'XX' if t2 == 'XY' else 'XY'
        hard = cr.Cell(jumps=[(ct.PAULI['Z'], 1)], windows=[(lead, (0, 1), 1), (t2, (0, 1), 1)])
        assert cr.every_n(hard, 'pi') is None       # hard rows pair at N = 2 but are refused for every N


def test_t1_alone_needs_the_machine_for_its_none():
    # cooling only: the pre-check passes (⊗(Z - I) is an eigenprojector of H), the machine proves NONE
    for down in ([R(1), R(1), R(1)], DOWN):
        s = ct.Setting(3, bond(0, 1, 1, 0, 1) + bond(1, 2, 1, 0, 1), [((l,), ct.LOWER, d) for l, d in enumerate(down)])
        assert cr.precheck(s, 'pi')[0] == 'PASS'
        v = cr.answer(s, 'pi', timeout=None)
        assert v.kind == 'NONE' and v.note.startswith('machine')
        assert cr.answer(s, 'refl', timeout=None).kind == 'CUBE MOVES'


@pytest.mark.slow
def test_every_n_router_is_periodic_four_and_the_small_n_artifact_is_refused():
    router = cr.Cell(jumps=[(ct.PAULI['Z'], 1)],
                     windows=[('XZX', (0, 1, 2), 3), ('XZY', (0, 1, 2), 2), ('YZX', (0, 1, 2), 2)])
    assert cr.every_n(router, 'pi') == ('PERIODIC', 4)
    artifact = cr.Cell(jumps=[(ct.PAULI['Z'], 1)],
                       windows=[('XZX', (0, 1, 2), 3), ('XZY', (0, 1, 2), 2), ('YZX', (0, 1, 2), 1)])
    assert cr.every_n(artifact, 'pi') is None


def test_xz_plus_zxz_has_none():
    # its spectrum does not pair (Tr((L - mu)^5) = 46080 at these rates); the pre-check alone decides it
    terms = [('XZ', (a, a + 1), 1) for a in range(3)] + [('ZXZ', (a, a + 1, a + 2), 1) for a in range(2)]
    s = ct.Setting(4, terms, [((l,), ct.PAULI['Z'], R(1, 1 + l)) for l in range(4)])
    for rel in cr.RELATIONS:
        assert cr.answer(s, rel, timeout=None).kind == 'NONE'


def test_identity_letters_inside_a_term_do_not_change_the_answer():
    # a term is traceless only where its letter is not I; padded with I it must land in its true support
    jumps = [((l,), ct.PAULI['Z'], 1) for l in range(2)]
    as_fields = ct.Setting(2, [('X', (0,), 1), ('Y', (0,), 1)], jumps)
    padded = ct.Setting(2, [('XI', (0, 1), 1), ('Y', (0,), 1)], jumps)
    for rel in cr.RELATIONS:
        assert cr.answer(as_fields, rel, timeout=None).kind == 'CUBE MOVES'
        assert cr.answer(padded, rel, timeout=None).kind == 'CUBE MOVES'


def test_every_n_reads_windows_at_any_offset():
    for off in ((0, 1), (1, 2)):
        hard = cr.Cell(jumps=[(ct.PAULI['Z'], 1)], windows=[('XX', off, 1), ('XY', off, 1)])
        assert cr.every_n(hard, 'pi') is None
    heis = [('XX', (0, 1), 1), ('YY', (0, 1), 1), ('ZZ', (0, 1), 1)]
    by_field = cr.Cell(jumps=[(ct.PAULI['Z'], 1)], fields=[('X', 1)], windows=heis)
    by_window = cr.Cell(jumps=[(ct.PAULI['Z'], 1)], windows=heis + [('X', (1,), 1)])
    assert cr.every_n(by_field, 'pi') == cr.every_n(by_window, 'pi')


def test_outside_is_reported():
    s = ct.Setting(4, [('XXXX', (0, 1, 2, 3), 1)])
    assert cr.answer(s, 'pi').kind == 'OUTSIDE'
