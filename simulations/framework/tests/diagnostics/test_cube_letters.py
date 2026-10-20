"""Tests for framework.diagnostics.cube_letters, the letter layer of the Theodolite. Read elements are checked by an
independent route: the full generator built row-major from the Lindblad formula, only then permuted to per-site order,
must satisfy E·L·E⁻¹ = 𝒯(L) - c exactly. Silent rows have a reason to be silent (no site-local element at all, or a
spectrum that does not pair at every N); slow rows are marked."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics import cube_letters as cl  # noqa: E402
from framework.diagnostics import cube_reflector as cr  # noqa: E402

R = sp.Rational
X, Y, Z = ct.PAULI['X'], ct.PAULI['Y'], ct.PAULI['Z']
HEIS = [('XX', (0, 1), 1), ('YY', (0, 1), 1), ('ZZ', (0, 1), 1)]
DOWN = [R(1, 2), R(1), R(3, 2)]


def bond(a, b, J=1, D=0, Dz=0):
    t = [('XX', (a, b), J), ('YY', (a, b), J)] if J else []
    t += [('XY', (a, b), D), ('YX', (a, b), -D)] if D else []
    return t + ([('ZZ', (a, b), Dz)] if Dz else [])


def heisenberg(n):
    return [t for a in range(n - 1) for t in bond(a, a + 1, 1, 0, 1)]


def full_generator(setting):
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


def handshake(setting, relation, maps):
    E = sp.kronecker_product(*maps)
    L = full_generator(setting)
    rhs = (cr.target(L, relation) - cr.constant_c(setting) * sp.eye(L.shape[0])) * E
    return (E * L - rhs).applyfunc(lambda e: sp.simplify(sp.expand(e))).is_zero_matrix


def read_and_check(setting, relation):
    rd = cl.read(setting, relation)
    assert rd is not None
    assert cl.verify(setting, relation, rd.maps)
    assert handshake(setting, relation, rd.maps)
    return rd


# ---------------------------------------------------------------- the binomial solver

def test_solve_binomials_takes_roots_and_frees_the_gauge():
    x, y, z = sp.symbols('x y z')
    sol = cl.solve_binomials([({x: 1, y: -1}, -3, 'a'), ({x: 1, y: 1}, -1, 'b')], [x, y])
    assert sp.simplify(sol[x] ** 2 - 3) == 0 and sp.simplify(sol[x] * sol[y] + 1) == 0
    sol = cl.solve_binomials([({x: 1, y: -1}, 1, 'a')], [x, y, z])
    assert sol == {x: 1, y: 1, z: 1}


def test_solve_binomials_closes_cycles_exactly():
    x, y, z = sp.symbols('x y z')
    ring = [({x: 1, y: -1}, sp.I, 'a'), ({y: 1, z: -1}, sp.I, 'b')]
    assert cl.solve_binomials(ring + [({z: 1, x: -1}, -1, 'c')], [x, y, z]) is not None
    assert cl.solve_binomials(ring + [({z: 1, x: -1}, 1, 'c')], [x, y, z]) is None
    assert cl.solve_binomials([({x: 1}, 2, 'a'), ({x: 1}, 3, 'b')], [x]) is None


# ---------------------------------------------------------------- reading settings

def test_reads_pi_under_z_dephasing():
    s = ct.Setting(3, heisenberg(3), [((l,), Z, g) for l, g in enumerate([R(1, 5), R(2, 7), R(1, 3)])])
    rd = read_and_check(s, 'pi')
    assert rd.transposed and all(k == ('d', 'a') for k in rd.kinds)


def test_reads_the_thermal_element_weighted_by_the_steady_state():
    s = ct.Setting(3, heisenberg(3), [((l,), ct.LOWER, 1) for l in range(3)] + [((l,), X, R(1, 2)) for l in range(3)])
    rd = read_and_check(s, 'pi')
    for x, y in zip(rd.x, rd.y):            # x/y = −(γ↓ + γ⊥)/(γ↑ + γ⊥) = −3 and x·y = −1: x² = 3
        assert sp.simplify(x / y + 3) == 0 and sp.simplify(x * y + 1) == 0


def test_the_gradient_is_silent_for_pi_and_reads_s_xy_for_the_reflection():
    s = ct.Setting(3, heisenberg(3), [((l,), ct.LOWER, d) for l, d in enumerate(DOWN)]
                   + [((l,), ct.RAISE, u) for l, u in enumerate([R(1, 8), R(1, 3), R(1, 5)])])
    assert cl.read(s, 'pi') is None            # nothing site-local exists (the pre-check proves it)
    assert cr.precheck(s, 'pi')[0] == 'NONE'
    read_and_check(s, 'refl')


def test_reads_the_golden_bond_and_the_c0_window():
    golden = ct.Setting(3, [(p, (a, a + 1), 1) for a in range(2) for p in ('XX', 'XY', 'YX')],
                        [((l,), Z, 1) for l in range(3)])
    read_and_check(golden, 'pi')
    window = ct.Setting(4, [(p, (w, w + 1, w + 2), 1) for w in range(2) for p in ('XZY', 'YZX')],
                        [((l,), Z, 1) for l in range(4)])
    read_and_check(window, 'pi')


def test_reads_f137s_chain_on_the_field_triangle_and_refuses_the_controls():
    def triangle(delta, J2, D2):
        terms = bond(0, 1, 1, 0, delta) + bond(1, 2, 1, 0, delta) + bond(2, 0, J2, D2, delta)
        terms += [('Z', (l,), h) for l, h in enumerate((R(3, 10), R(-1, 5), R(1, 2)))]
        return ct.Setting(3, terms, [((l,), ct.LOWER, d) for l, d in enumerate(DOWN)]
                          + [((l,), ct.RAISE, d / 4) for l, d in enumerate(DOWN)])
    rd = read_and_check(triangle(0, 0, 1), 'pi')
    assert not rd.transposed and all(k == ('a', 'a') for k in rd.kinds)
    for x, y in zip(rd.x, rd.y):               # the spin flip with the ratio of the arrows: x/y = −γ↓/γ↑
        assert sp.simplify(x / y + 4) == 0
    assert cl.read(triangle(R(7, 10), 0, 1), 'pi') is None
    assert cl.read(triangle(0, 1, 0), 'pi') is None


def test_identity_letters_inside_a_term_do_not_change_the_reading():
    jumps = [((l,), Z, 1) for l in range(2)]
    for s in (ct.Setting(2, [('X', (0,), 1), ('Y', (0,), 1)], jumps),
              ct.Setting(2, [('XI', (0, 1), 1), ('Y', (0,), 1)], jumps)):
        for rel in cl.cr.RELATIONS:
            read_and_check(s, rel)


def test_answer_reads_the_letters_first():
    s = ct.Setting(3, heisenberg(3), [((l,), ct.LOWER, 1) for l in range(3)] + [((l,), X, R(1, 2)) for l in range(3)])
    v = cr.answer(s, 'pi', timeout=None)
    assert v.kind == 'CUBE MOVES' and v.note.startswith('read from the letters')
    assert handshake(s, 'pi', v.point)


# ---------------------------------------------------------------- every N

def test_every_n_uniform_rows():
    assert cl.every_n(cr.Cell(jumps=[(Z, 1)], windows=HEIS), 'pi')[0] == 'UNIFORM'
    assert cl.every_n(cr.Cell(jumps=[(ct.LOWER, 1), (X, R(1, 2))], windows=HEIS), 'pi')[0] == 'UNIFORM'
    golden = cr.Cell(jumps=[(Z, 1)], windows=[('XX', (0, 1), 1), ('XY', (0, 1), 1), ('YX', (0, 1), 1)])
    assert cl.every_n(golden, 'pi')[0] == 'UNIFORM'
    c0 = cr.Cell(jumps=[(Z, 1)], windows=[('XZY', (0, 1, 2), 1), ('YZX', (0, 1, 2), 1)])
    assert cl.every_n(c0, 'pi')[0] == 'UNIFORM'


def test_every_n_heisenberg_dm_turns_by_arctan_d_over_j():
    dm = cr.Cell(jumps=[(Z, 1)], windows=HEIS + [('XY', (0, 1), R(1, 2)), ('YX', (0, 1), -R(1, 2))])
    out = cl.every_n(dm, 'pi')
    assert out[0] == 'TURNED'
    z = out[1]['turn']                         # e^{i·arctan(1/2)} = (2 + i)/√5, up to a turn by π
    assert sp.simplify(z ** 2 - ((2 + sp.I) / sp.sqrt(5)) ** 2) == 0
    pure = cr.Cell(jumps=[(Z, 1)], windows=[('XY', (0, 1), 1), ('YX', (0, 1), -1)])
    out = cl.every_n(pure, 'pi')
    assert out[0] == 'TURNED' and sp.simplify(out[1]['turn'] ** 2 + 1) == 0     # a quarter turn


def test_every_n_hard_rows_are_silent():
    for t1, t2 in (('XX', 'XY'), ('XY', 'XZ')):
        assert cl.every_n(cr.Cell(jumps=[(Z, 1)], windows=[(t1, (0, 1), 1), (t2, (0, 1), 1)]), 'pi') is None


def test_every_n_rate_free():
    thermal = [(ct.LOWER, 1), (ct.RAISE, R(1, 4))]
    assert cl.every_n(cr.Cell(jumps=thermal, windows=HEIS), 'refl', rate_free=True)[0] == 'UNIFORM'
    dm = cr.Cell(jumps=thermal, windows=[('XX', (0, 1), 4), ('YY', (0, 1), 4), ('XY', (0, 1), 3), ('YX', (0, 1), -3)])
    assert cl.every_n(dm, 'refl', rate_free=True)[0] == 'TURNED'


@pytest.mark.slow
def test_every_n_golden_window_is_periodic_four():
    golden = cr.Cell(jumps=[(Z, 1)], windows=[('XZX', (0, 1, 2), 1), ('XZY', (0, 1, 2), 1), ('YZX', (0, 1, 2), 1)])
    out = cl.every_n(golden, 'pi')
    assert out[0] == 'PERIODIC' and out[1] == 4
