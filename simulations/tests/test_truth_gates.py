"""Truth gates: executable checks of TRUE physics and arithmetic that the repository states.

Each gate asserts a fact from below: an exact rational or matrix identity (sympy,
fractions), a counterexample that must stay live, or the behaviour of a producer
function copied by source into a scratch module. No producer is imported into the
test process or run as a script; the one import-safety gate below imports a producer
in a separate interpreter and guards its result file. A gate that carries a mutation
shows that it can fail.

Three groups are not physics but protect a true statement where it is written:
  * F47 path-metric labels: a one-dimensional path metric has no intrinsic curvature,
    so the four hosts of the Bures path coefficient must name the coordinate-shape
    second derivative and must not call it a Gaussian or Bures curvature.
  * Calibration-producer labels: the r = T2/(2 T1) < R* band is the free-|+>
    normalized-purity proxy, not a CΨ = 1/4 crossing (the free-|+> CΨ falls from 1
    to 0 and so crosses 1/4 at every r; F24), and the Kingston F97 lens is a radial
    diagnostic, not the F97 parameter.
  * Import safety of simulations/veffect_cavity_modes.py: importing it in a
    subprocess must print nothing and must leave its tracked result file byte-identical.

Run with explicit paths only (python -m pytest simulations/tests/test_truth_gates.py);
collecting simulations/ as a whole imports one-shot producers that write at import.
"""
from __future__ import annotations

import ast
import os
import runpy
import struct
import subprocess
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

ROOT = Path(os.environ.get("RCPSI_REPO_ROOT") or Path(__file__).resolve().parents[2])


def read_host(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8-sig").replace("\r\n", "\n").replace("\r", "\n")


def source_copy_definitions(source, names, tmp_path, *, prelude="import numpy as np\n", assignments=()):
    """Run only the named top-level functions (and named constant assignments) of a
    producer source in a scratch module. Module-level plotting, file opens and
    repository imports of the producer are never evaluated."""
    tree = ast.parse(source)
    selected = [node for node in tree.body
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names]
    assert {node.name for node in selected} == set(names)
    constants = [node for node in tree.body
                 if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id in assignments for target in node.targets)]
    assert {target.id for node in constants for target in node.targets} >= set(assignments)
    scratch = tmp_path / "isolated_source.py"
    scratch.write_text(prelude + "\n" + ast.unparse(ast.Module(body=constants + selected, type_ignores=[])),
                       encoding="utf-8", newline="\n")
    return runpy.run_path(str(scratch))


# ----------------------------------------------------------------------------- F1
def f1_partner(z, sigma):
    """Independent holomorphic F1 oracle, not a production implementation."""
    return -z - 2 * sigma


def maximum_cardinality_f1_match(values, sigma, tolerance):
    values = np.asarray(values, dtype=complex)
    graph = abs(values[None, :] - f1_partner(values[:, None], sigma)) < tolerance
    return int(np.count_nonzero(maximum_bipartite_matching(csr_matrix(graph)) >= 0))


F1_CASES = (
    ((2 + 3j, -2 - 3j), 0, 1e-12, 2),
    ((2 + 3j, -8 - 3j), 3, 1e-12, 2),
    ((1, 1, 1, -1), 0, 1e-12, 2),
    (1.9e-7 * np.array((2 - 7j, -5 + 3j, 1j, 3 + 3j)), 0, 1e-6, 4),
    ((2.0 ** -21,), 0, 2.0 ** -19, 1),
    ((2.0 ** -21,), 0, 2.0 ** -20, 0),
    ((-2 + 2.0 ** -21,), 2, 2.0 ** -20, 0),
    ((1e-6 / 2,), 0, 1e-6, 0),
)


@pytest.mark.parametrize("values,sigma,tolerance,expected", F1_CASES)
def test_f1_independent_multiset_matching(values, sigma, tolerance, expected):
    assert maximum_cardinality_f1_match(values, sigma, tolerance) == expected


def test_f1_partner_and_strict_boundary_are_exact():
    d = 2.0 ** -20
    assert f1_partner(2 + 3j, 0) == -2 - 3j
    assert f1_partner(2 + 3j, 3) == -8 - 3j
    z = -2 + d / 2
    assert f1_partner(z, 2) == -2 - d / 2
    assert abs(z - f1_partner(z, 2)) == d
    assert max(1, abs(f1_partner(z, 2))) == 2 + d / 2
    assert abs(d / 2 - f1_partner(d / 2, 0)) == d
    assert abs(1e-6 / 2 - f1_partner(1e-6 / 2, 0)).hex() == (1e-6).hex()


def test_f1_production_matcher_preserves_multiplicity_without_recursion(tmp_path):
    source = read_host("simulations/pairing_structure.py")
    copied = source_copy_definitions(
        source, {"f1_partner", "maximum_cardinality_f1_match"}, tmp_path,
        prelude="import numpy as np\nfrom scipy.sparse import csr_matrix\n"
                "from scipy.sparse.csgraph import maximum_bipartite_matching\n")
    match = copied["maximum_cardinality_f1_match"]
    for values, sigma, tolerance, expected in F1_CASES:
        assert match(values, sigma, tolerance) == expected
    recursion_limit = sys.getrecursionlimit()
    assert match(np.zeros(1024), 0, 1e-12) == 1024
    assert sys.getrecursionlimit() == recursion_limit


# ------------------------------------------------- CPsi monotonicity and absorbers
def _outer(ket, bra=None):
    bra = ket if bra is None else bra
    return ket * bra.conjugate().T


def _dissipator(jump, rho):
    adjoint = jump.conjugate().T
    return sp.simplify(jump * rho * adjoint - (adjoint * jump * rho + rho * adjoint * jump) / 2)


def _rhs(rho, hamiltonian, jumps):
    result = -sp.I * (hamiltonian * rho - rho * hamiltonian)
    for jump in jumps:
        result += _dissipator(jump, rho)
    return sp.simplify(result)


def _l1(rho):
    return sp.simplify(sum(sp.Abs(rho[i, j]) for i in range(rho.rows) for j in range(rho.rows) if i != j))


def _cpsi(rho):
    assert rho.rows == 4
    return sp.simplify(sp.trace(rho * rho) * _l1(rho) / 3)


def _assert_density(rho):
    assert rho == rho.conjugate().T
    assert sp.simplify(sp.trace(rho)) == 1
    assert all(eigenvalue.is_nonnegative is True for eigenvalue in rho.eigenvals())


def _directional_tuple(rho, rho_dot):
    """Exact (dP/dt, dl1/dt, dCPsi/dt) at t = 0 along rho_dot."""
    purity = sp.simplify(sp.trace(rho * rho))
    l1 = _l1(rho)
    purity_dot = sp.simplify(2 * sp.re(sp.trace(rho * rho_dot)))
    l1_dot = 0
    for i in range(rho.rows):
        for j in range(rho.rows):
            if i == j:
                continue
            if rho[i, j] == 0:
                # l1 stays differentiable: none of these exact germs creates a new coherence at t = 0.
                assert sp.simplify(rho_dot[i, j]) == 0
            else:
                l1_dot += sp.re(sp.conjugate(rho[i, j]) * rho_dot[i, j]) / sp.Abs(rho[i, j])
    l1_dot = sp.simplify(l1_dot)
    return purity_dot, l1_dot, sp.simplify((purity_dot * l1 + purity * l1_dot) / 3)


def _pointwise_fixture():
    identity = sp.eye(2)
    x = sp.Matrix([[0, 1], [1, 0]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    z = sp.diag(1, -1)
    projector0 = _outer(sp.Matrix([1, 0]))
    rho = sp.kronecker_product((identity + x / 2 + z / 2) / 2, projector0)
    hamiltonian = sp.kronecker_product(y, identity)
    z_first = sp.kronecker_product(z, identity)
    dephasing_jump = z_first / 2  # rate 1/4 on the first site
    return identity, x, z, projector0, rho, hamiltonian, sp.zeros(4), z_first, dephasing_jump


def _assert_positive_pointwise_germ(rho, hamiltonian, dephasing_jump):
    assert _directional_tuple(rho, _rhs(rho, hamiltonian, (dephasing_jump,))) == (
        sp.Rational(-1, 8), sp.Rational(3, 4), sp.Rational(1, 6))


def test_monotonicity_exact_pointwise_counterexample_and_h_zero_mutation():
    _, _, _, _, rho, hamiltonian, zero, _, jump = _pointwise_fixture()
    _assert_density(rho)
    assert _cpsi(rho) == sp.Rational(1, 8)
    _assert_positive_pointwise_germ(rho, hamiltonian, jump)
    assert _directional_tuple(rho, _rhs(rho, zero, (jump,))) == (
        sp.Rational(-1, 8), sp.Rational(-1, 4), sp.Rational(-1, 12))
    with pytest.raises(AssertionError):
        _assert_positive_pointwise_germ(rho, zero, jump)


def test_monotonicity_hadamard_change_and_same_lab_z_pulse_are_exact():
    identity, x, z, projector0, rho, hamiltonian, _, z_first, jump = _pointwise_fixture()
    rho00 = sp.kronecker_product(projector0, projector0)
    local_hadamard = sp.kronecker_product((x + z) / sp.sqrt(2), identity)
    after_hadamard = sp.simplify(local_hadamard * rho00 * local_hadamard.conjugate().T)
    _assert_density(after_hadamard)
    assert _cpsi(rho00) == 0
    assert _cpsi(after_hadamard) == sp.Rational(1, 3)

    pulsed = sp.simplify(z_first * rho * z_first)
    assert _cpsi(pulsed) == sp.Rational(1, 8)
    assert _directional_tuple(pulsed, _rhs(pulsed, hamiltonian, (jump,))) == (
        sp.Rational(-1, 8), sp.Rational(-5, 4), sp.Rational(-1, 3))
    passive_hamiltonian = sp.simplify(z_first * hamiltonian * z_first)
    assert _directional_tuple(pulsed, _rhs(pulsed, passive_hamiltonian, (jump,))) == (
        sp.Rational(-1, 8), sp.Rational(3, 4), sp.Rational(1, 6))


def _liouvillian_nullity(dimension, jumps):
    columns = []
    for i in range(dimension):
        for j in range(dimension):
            basis = sp.zeros(dimension)
            basis[i, j] = 1
            columns.append(sp.Matrix(list(_rhs(basis, sp.zeros(dimension), jumps))))
    return len(sp.Matrix.hstack(*columns).nullspace())


def _assert_unique_stationary_target(first_jump, second_jump, target, expected_cpsi):
    _assert_density(target)
    assert _dissipator(first_jump, target) == sp.zeros(4)
    assert _dissipator(second_jump, target) == sp.zeros(4)
    assert _liouvillian_nullity(4, (first_jump, second_jump)) == 1
    assert _cpsi(target) == expected_cpsi


def test_absorber_local_markov_stationary_state_and_axis_mutation_are_exact():
    identity = sp.eye(2)
    ket0, ket1 = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    plus, minus = (ket0 + ket1) / sp.sqrt(2), (ket0 - ket1) / sp.sqrt(2)
    projector0 = _outer(ket0)
    coherent_axis_jump = sp.kronecker_product(_outer(plus, minus), identity)
    second_jump = sp.kronecker_product(identity, _outer(ket0, ket1))
    rho_star = sp.kronecker_product(_outer(plus), projector0)
    _assert_unique_stationary_target(coherent_axis_jump, second_jump, rho_star, sp.Rational(1, 3))

    computational_axis_jump = sp.kronecker_product(_outer(ket0, ket1), identity)
    rho00 = sp.kronecker_product(projector0, projector0)
    assert _rhs(rho00, sp.zeros(4), (computational_axis_jump, second_jump)) == sp.zeros(4)
    assert _cpsi(rho00) == 0
    with pytest.raises(AssertionError):
        _assert_unique_stationary_target(computational_axis_jump, second_jump, rho00, sp.Rational(1, 3))


def _rank_one_damping_state(target, source, initial, q):
    assert sp.simplify((target.conjugate().T * target)[0]) == 1
    assert sp.simplify((source.conjugate().T * source)[0]) == 1
    assert sp.simplify((target.conjugate().T * source)[0]) == 0
    aa = sp.simplify((target.conjugate().T * initial * target)[0])
    ab = sp.simplify((target.conjugate().T * initial * source)[0])
    ba = sp.simplify((source.conjugate().T * initial * target)[0])
    bb = sp.simplify((source.conjugate().T * initial * source)[0])
    return sp.simplify((aa + (1 - q) * bb) * _outer(target) + q * bb * _outer(source)
                       + sp.sqrt(q) * ab * _outer(target, source) + sp.sqrt(q) * ba * _outer(source, target))


def _cpsi_on_open_unit_interval(rho, q):
    l1 = sp.simplify(sum(sp.refine(sp.Abs(rho[i, j]), sp.Q.lt(q, 1))
                         for i in range(rho.rows) for j in range(rho.cols) if i != j))
    return sp.factor(sp.trace(rho * rho) * l1 / 3)


def _assert_exact_upward_crossing_for_axis(target, source):
    # q = e^(-t), so d/dt = -q d/dq for the unit-rate local semigroup.
    t = sp.symbols("t", nonnegative=True)
    assert sp.diff(sp.exp(-t), t) == -sp.exp(-t)
    q = sp.symbols("q", positive=True)
    identity = sp.eye(2)
    ket0, ket1 = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    projector0 = _outer(ket0)
    rho00 = sp.kronecker_product(projector0, projector0)
    first_jump = sp.kronecker_product(_outer(target, source), identity)
    second_jump = sp.kronecker_product(identity, _outer(ket0, ket1))
    rho_q = sp.kronecker_product(_rank_one_damping_state(target, source, projector0, q), projector0)
    assert sp.simplify(rho_q.subs(q, 1) - rho00) == sp.zeros(4)
    _assert_density(rho_q.subs(q, sp.Rational(1, 8)))
    assert sp.simplify(-q * rho_q.diff(q) - _rhs(rho_q, sp.zeros(4), (first_jump, second_jump))) == sp.zeros(4)

    cpsi_q = _cpsi_on_open_unit_interval(rho_q, q)
    assert sp.simplify(cpsi_q - (1 - q) * (q ** 2 - q + 2) / 6) == 0
    assert sp.factor(-q * sp.diff(cpsi_q, q)) == q * (3 * q ** 2 - 4 * q + 3) / 6
    assert sp.discriminant(3 * q ** 2 - 4 * q + 3, q) == -20
    assert cpsi_q.subs(q, 1) == 0
    assert cpsi_q.subs(q, sp.Rational(1, 8)) == sp.Rational(847, 3072)
    assert sp.Rational(847, 3072) - sp.Rational(1, 4) == sp.Rational(79, 3072)
    assert cpsi_q.subs(q, 0) == sp.Rational(1, 3)
    return cpsi_q


def test_absorber_exact_upward_crossing_from_below_under_fixed_local_semigroup():
    q = sp.symbols("q", positive=True)
    ket0, ket1 = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    plus, minus = (ket0 + ket1) / sp.sqrt(2), (ket0 - ket1) / sp.sqrt(2)
    projector0 = _outer(ket0)
    assert _assert_exact_upward_crossing_for_axis(plus, minus).subs(q, sp.Rational(1, 8)) == sp.Rational(847, 3072)

    # The computational axis generates the stationary |00><00| curve itself; it never crosses upward.
    mutated_rho = sp.kronecker_product(_rank_one_damping_state(ket0, ket1, projector0, q), projector0)
    mutated_curve = _cpsi_on_open_unit_interval(mutated_rho, q)
    assert mutated_rho == sp.kronecker_product(projector0, projector0)
    assert mutated_curve == 0
    assert mutated_curve.subs(q, sp.Rational(1, 8)) < sp.Rational(1, 4)
    with pytest.raises(AssertionError):
        _assert_exact_upward_crossing_for_axis(ket0, ket1)


# ----------------------------------------------------------- Bell+ and the F86 toy
def test_bellplus_spin_flip_eigenvalues_and_concurrence_are_analytic():
    f = sp.symbols("f", real=True, nonnegative=True)
    rho = sp.Matrix([[sp.Rational(1, 2), 0, 0, f / 2], [0, 0, 0, 0],
                     [0, 0, 0, 0], [f / 2, 0, 0, sp.Rational(1, 2)]])
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    flip = sp.kronecker_product(y, y) * sp.conjugate(rho) * sp.kronecker_product(y, y)
    assert flip == rho
    lam = sp.Symbol("lambda")
    assert (rho * flip).charpoly().as_expr().factor() == (
        lam ** 2 * (lam - (1 + f) ** 2 / 4) * (lam - (1 - f) ** 2 / 4)).expand().factor()
    # 0 <= f <= 1 orders the nonnegative square roots without a float eigensolver.
    roots = ((1 + f) / 2, (1 - f) / 2, 0, 0)
    assert sp.simplify(roots[0] - sum(roots[1:])) == f
    gamma, t = sp.symbols("gamma t", positive=True, finite=True)
    assert sp.exp(-4 * gamma * t).is_positive
    radius = f * (1 + f * f) / 6
    assert radius.subs(f, sp.Rational(1, 2)) == sp.Rational(5, 48) < sp.Rational(1, 4)
    assert (roots[0] - sum(roots[1:])).subs(f, sp.Rational(1, 2)) == sp.Rational(1, 2)
    assert sp.diff(radius, f) == f * f / 3 + (f * f + 1) / 6
    crossing_f = brentq(lambda x: x ** 3 + x - 1.5, 0, 1, xtol=1e-14)
    assert abs(crossing_f - 0.8612240997395736) < 1e-14
    assert 0 < crossing_f < 1
    assert abs(crossing_f * (1 + crossing_f ** 2) / 6 - .25) < 1e-14
    assert -np.log(crossing_f) / 4 > 0


def test_f86_toy_pair_lifetime_halves():
    gamma = sp.symbols("gamma", positive=True)
    j = sp.symbols("J", nonnegative=True)
    # The k=1 toy pair has rates -4 gamma +/- sqrt(4 gamma^2 - J^2).
    slow = -4 * gamma + sp.sqrt(4 * gamma ** 2 - j ** 2)
    assert sp.simplify(-1 / slow.subs(j, 0)) == 1 / (2 * gamma)
    assert sp.simplify(-1 / slow.subs(j, 2 * gamma)) == 1 / (4 * gamma)
    assert slow.subs(j, 2 * gamma) != 0


# --------------------------------------------------- dephasing channels and F14
def test_collective_and_local_z_jump_dissipators_are_distinct_exactly():
    z = sp.diag(1, -1)
    a, b = sp.kronecker_product(z, sp.eye(2)), sp.kronecker_product(sp.eye(2), z)
    rho = sp.Matrix([[sp.Rational(1, 2), 0, 0, sp.Rational(1, 2)], [0, 0, 0, 0],
                     [0, 0, 0, 0], [sp.Rational(1, 2), 0, 0, sp.Rational(1, 2)]])

    def dissipator(jump):
        square = jump.H * jump
        return jump * rho * jump.H - (square * rho + rho * square) / 2

    product, collective, local = dissipator(a * b), dissipator(a + b), dissipator(a) + dissipator(b)
    assert product == sp.zeros(4)
    assert collective != local
    assert [sum(abs(x) ** 2 for x in matrix) for matrix in (product, collective, local)] == [0, 32, 8]


def test_markovian_log_coherence_can_curve_exactly():
    x, y, z = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)
    rho = (sp.eye(2) + x / 2 + z / 2) / 2

    def derivatives(h):
        def generator(state):
            return -sp.I * (h * state - state * h) + (z * state * z - state) / 4

        first, second = generator(rho), generator(generator(rho))
        l, dl, ddl = 2 * rho[0, 1], 2 * first[0, 1], 2 * second[0, 1]
        return dl / l, sp.simplify(ddl / l - (dl / l) ** 2)

    assert derivatives(y) == (sp.Rational(3, 2), -7)
    assert derivatives(sp.zeros(2)) == (-sp.Rational(1, 2), 0)


def _nonuniform_eigenstate_control():
    """Exact 32x32 matrices, no trajectory rerun and no eigensolver."""
    letters = (sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1))

    def site(letter, index):
        return sp.kronecker_product(*(letter if j == index else sp.eye(2) for j in range(5)))

    h = sp.zeros(32)
    for j in range(4):
        for letter in letters:
            h += site(letter, j) * site(letter, j + 1)
    rho = sp.ones(32) / 32
    d = sp.zeros(32)
    for j, rate in enumerate((5, 1, 1, 1, 1)):
        z = site(letters[2], j)
        d += rate * (z * rho * z - rho)
    return h * rho - rho * h, h * d - d * h


def test_nonuniform_dephasing_leaves_the_initial_h_eigenspace():
    initial, later = _nonuniform_eigenstate_control()
    assert initial == sp.zeros(32)
    assert max(abs(x) for x in later) == sp.Rational(1, 2)
    assert sum(abs(x) ** 2 for x in later) == 128
    assert max(abs(x) / 20 for x in later) == sp.Rational(1, 40)
    assert sum(abs(x / 20) ** 2 for x in later) == sp.Rational(8, 25)
    with pytest.raises(AssertionError):
        assert later == sp.zeros(32)  # "purely dephasing-driven" would need this commutator to vanish
    u = brentq(lambda u: u ** 3 + u - 0.5, 0., 1.)
    independent = [-np.log(u) / rate for rate in (.5, .1, .1, .1, .1)]
    assert np.allclose(independent, [1.716733, 8.583667, 8.583667, 8.583667, 8.583667], rtol=0, atol=1e-6)
    assert not np.allclose(independent, [4.71, 5.39, 5.27, 5.00, 4.44], rtol=0, atol=.005)


def test_markovian_dynamics_can_transiently_create_coherence():
    x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    hamiltonian = x / 2.0

    def generator(rho):
        return -1j * (hamiltonian @ rho - rho @ hamiltonian) + 0.1 * (z @ rho @ z - rho)

    basis = []
    for row in range(2):
        for column in range(2):
            matrix = np.zeros((2, 2), dtype=complex)
            matrix[row, column] = 1.0
            basis.append(matrix)
    liouvillian = np.column_stack([generator(matrix).reshape(-1) for matrix in basis])
    rho_initial = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
    rho_later = (expm(0.1 * liouvillian) @ rho_initial.reshape(-1)).reshape(2, 2)
    assert abs(rho_initial[0, 1]) == 0.0
    assert abs(rho_later[0, 1]) > 0.049


# ------------------------------------------------------------ Born shadow, channels
def test_born_shadow_identity_is_linearity_not_a_born_derivation():
    past = sp.Matrix([[sp.Rational(2, 5), sp.Rational(1, 7)], [sp.Rational(1, 7), sp.Rational(1, 10)]])
    future = sp.Matrix([[sp.Rational(1, 10), -sp.Rational(1, 7)], [-sp.Rational(1, 7), sp.Rational(2, 5)]])
    rho = past + future
    assert rho.diagonal() == past.diagonal() + future.diagonal()
    assert sp.trace(rho * rho) == sp.trace(past * past) + sp.trace(future * future) + 2 * sp.trace(past * future)


def test_primitive_replacement_channel_counterexample_is_executable(tmp_path):
    source = read_host("simulations/subsystem_crossing.py")
    names = {"purity", "psi_norm", "cpsi", "primitive_replacement_fixed_point"}
    copied = source_copy_definitions(source, names, tmp_path)
    sigma = copied["primitive_replacement_fixed_point"]()
    np.testing.assert_allclose(np.linalg.eigvalsh(sigma), [0.0125, 0.0125, 0.0125, 0.9625], rtol=0, atol=4e-16)
    assert copied["cpsi"](sigma) == pytest.approx(28177 / 96000, rel=0, abs=3e-16)
    assert copied["cpsi"](sigma) > 0.25

    mutant = source.replace("0.95 * bell_projector", "0.05 * bell_projector", 1)
    assert mutant != source
    mutant_copy = source_copy_definitions(mutant, names, tmp_path)
    assert mutant_copy["cpsi"](mutant_copy["primitive_replacement_fixed_point"]()) != pytest.approx(
        28177 / 96000, rel=0, abs=3e-16)


def test_no_signalling_verdict_uses_reduced_state_not_scalar(tmp_path):
    source = read_host("simulations/test2_no_signalling.py")
    namespace = source_copy_definitions(source, {"quarter_band", "no_signalling_reading"}, tmp_path)
    assert namespace["quarter_band"](0.30) == "above 1/4"
    assert namespace["quarter_band"](0.25) == "equal to 1/4"
    assert namespace["quarter_band"](0.20) == "below 1/4"
    before, after = np.diag([0.7, 0.3]), np.diag([0.3, 0.7])
    assert namespace["no_signalling_reading"](before, after)[1] is False

    eigenvalue_only = source.replace(
        "rho_A_before, rho_A_after, atol=tolerance",
        "np.linalg.eigvalsh(rho_A_before), np.linalg.eigvalsh(rho_A_after), atol=tolerance", 1)
    assert eigenvalue_only != source
    mutated = source_copy_definitions(eigenvalue_only, {"no_signalling_reading"}, tmp_path)["no_signalling_reading"]
    assert mutated(before, after)[1] is True


# ------------------------------------------------------------- F94-F97 arithmetic
def test_f94_f96_exact_rational_controls_are_distinct():
    ring_sym3 = tuple(map(Fraction, (8, -4, -4, 0)))
    chain_sym3 = tuple(map(Fraction, (5, -4, -1, 0)))
    assert sum(ring_sym3) == sum(chain_sym3) == 0
    assert ring_sym3[0] / 6 == Fraction(4, 3)
    assert ring_sym3[1] / (3 * Fraction(3, 4)) == Fraction(-16, 9)
    assert ring_sym3[2] / (3 * Fraction(3, 4)) == Fraction(-16, 9)
    assert Fraction(-20) / (5 * Fraction(3, 2)) == Fraction(-8, 3)
    assert chain_sym3[2] / (3 * Fraction(1, 4)) == Fraction(-4, 3)
    assert Fraction(40, 27) != Fraction(4, 3)


def test_qutrit_closed_form_is_40_over_27(tmp_path):
    source = read_host("simulations/f94_qutrit_born_mirror.py")
    prelude = "from fractions import Fraction\n"
    copied = source_copy_definitions(source, {"c_closed"}, tmp_path, prelude=prelude)
    assert copied["c_closed"](2) == Fraction(4, 3)
    assert copied["c_closed"](3) == Fraction(40, 27)
    mutant = source.replace("3 * d * d", "4 * d * d", 1)
    assert mutant != source
    assert source_copy_definitions(mutant, {"c_closed"}, tmp_path, prelude=prelude)["c_closed"](3) != Fraction(40, 27)


def test_formula_four_qutrit_refutes_ancestry(tmp_path):
    source = read_host("simulations/formula_four_classification.py")
    prelude = "from fractions import Fraction\n"
    function = source_copy_definitions(source, {"f94_qudit_coefficient"}, tmp_path, prelude=prelude)[
        "f94_qudit_coefficient"]
    assert function(2) == Fraction(4, 3)
    assert function(3) == Fraction(40, 27)
    assert function(3) != Fraction(3, 1)
    mutant = source.replace("3 * d * d", "4 * d * d", 1)
    assert mutant != source
    changed = source_copy_definitions(mutant, {"f94_qudit_coefficient"}, tmp_path, prelude=prelude)[
        "f94_qudit_coefficient"]
    assert changed(3) != Fraction(40, 27)


@pytest.mark.parametrize("b,c,expected", ((2.0, 5.0, np.arctan(0.5)), (1.0, 4 / 3, np.pi / 6), (7.0, 49.0, 0.0)))
def test_f95_positive_b_quadratic_angle(b, c, expected):
    actual = np.arctan(np.sqrt(c / b ** 2 - 1)) if c >= b ** 2 else np.nan
    assert actual == pytest.approx(expected, abs=1e-14)
    assert (np.arctan(np.sqrt((9 * c) / (3 * b) ** 2 - 1)) if c >= b ** 2 else np.nan) == pytest.approx(actual, abs=1e-14)
    with np.errstate(invalid="ignore"):
        assert np.isnan(np.sqrt(-1.0))


def test_f96_ring_slopes_and_topology_counterexample_are_distinct():
    assert (Fraction(-16, 9), Fraction(-16, 9), Fraction(-8, 3)) == (Fraction(-16, 9), Fraction(-16, 9), Fraction(-8, 3))
    assert Fraction(-4, 3) not in {Fraction(-16, 9), Fraction(-8, 3)}
    assert sum((8, -4, -4, 0)) == 0
    assert sum((5, -4, -1, 0)) == 0


def test_f97_period_one_cardioid_is_not_the_quarter_circle():
    for phi in (0.0, np.pi / 2, np.pi, 2 * np.pi - 1e-12):
        z = np.exp(1j * phi) / 2
        c = z - z * z
        roots = np.roots((1, -1, c))
        assert min(abs(root - z) for root in roots) < 2e-8
        assert abs(2 * z) == pytest.approx(1.0, abs=1e-14)
    assert abs((0.5j) - (0.5j) ** 2) ** 2 == pytest.approx(5 / 16)
    assert abs((0.5 + 0j) - (0.5 + 0j) ** 2) ** 2 == pytest.approx(1 / 16)
    assert tuple(np.roots((1, -1, 0))) == pytest.approx((1, 0))


def test_f97_selected_and_other_roots_are_distinct(tmp_path):
    source = read_host("simulations/cardioid_parametrization_tier1.py")
    names = {"period_one_cardioid_point", "fixed_points_at", "fixed_point_at", "wrapped_angle"}
    functions = source_copy_definitions(source, names, tmp_path)
    radial_magnitudes = []
    for phi in (0.0, np.pi / 2, np.pi, 2 * np.pi - 1e-12):
        z_star, c = functions["period_one_cardioid_point"](phi)
        radial_magnitudes.append(abs(c))
        assert z_star * z_star - z_star + c == pytest.approx(0j, abs=1e-14)
        assert abs(2 * z_star) == pytest.approx(1.0, rel=0, abs=1e-14)
        assert min(abs(root - z_star) for root in functions["fixed_points_at"](c)) < 2e-8
        assert functions["wrapped_angle"](z_star) == pytest.approx(phi % (2 * np.pi), rel=0, abs=2e-12)
    z_star, c = functions["period_one_cardioid_point"](np.pi / 2)
    assert abs(2 * (1 - z_star)) == pytest.approx(np.sqrt(5), rel=0, abs=1e-14)
    assert functions["fixed_points_at"](0j) == pytest.approx((0j, 1 + 0j))
    assert functions["fixed_point_at"](0j, 0.1 + 0j) == 0j
    assert functions["fixed_point_at"](0j, 0.9 + 0j) == 1 + 0j
    assert max(radial_magnitudes) - min(radial_magnitudes) == pytest.approx(0.5, rel=0, abs=1e-14)

    for old, new in (("c = z_star - z_star * z_star", "c = z_star + z_star * z_star"),
                     ("return (angle + 2 * np.pi) % (2 * np.pi)", "return angle")):
        mutant = source.replace(old, new, 1)
        assert mutant != source, old
        changed = source_copy_definitions(mutant, names, tmp_path)
        mz, mc = changed["period_one_cardioid_point"](np.pi / 2)
        checks = (abs(mz * mz - mz + mc) < 1e-14,
                  abs(changed["wrapped_angle"](changed["period_one_cardioid_point"](2 * np.pi - 1e-12)[0])
                      - (2 * np.pi - 1e-12)) < 2e-12)
        assert checks != (True, True), old

    selector_mutant = source.replace(
        "return min(fixed_points_at(c), key=lambda root: abs(root - target))", "return target", 1)
    assert selector_mutant != source
    changed = source_copy_definitions(selector_mutant, {"fixed_points_at", "fixed_point_at"}, tmp_path)
    assert changed["fixed_point_at"](0j, 0.1 + 0j) != 0j


# ------------------------------------------------------------------ R* and F47
def test_rstar_exact_double_bits_and_independent_stationarity_root():
    expected = 0.21275477982200533
    assert struct.unpack(">Q", struct.pack(">d", expected))[0] == 0x3FCB3B8C72B39BD9

    def y_stationary(r):
        return brentq(lambda y: -2 + 2 * y + (1 / r) * y ** (1 / r - 1), 0.5, 0.8)

    root = brentq(lambda r: (1 - 2 * y_stationary(r) + y_stationary(r) ** 2 + y_stationary(r) ** (1 / r)) - 0.25,
                  0.2127547798220052, 0.2127547798220055, xtol=5e-16)
    assert root == pytest.approx(expected, abs=5e-16)
    y = y_stationary(0.4)
    assert 1 - 2 * y + y * y + y ** (1 / 0.4) == pytest.approx(0.425334, abs=1e-6)


def test_f47_coordinate_shape_is_not_intrinsic_curvature():
    coordinate = sp.symbols("coordinate", positive=True)

    def coordinate_shape(metric_coefficient):
        return sp.simplify(-sp.diff(sp.log(metric_coefficient), coordinate, 2) / (2 * metric_coefficient))

    assert coordinate_shape(sp.Integer(1)) == 0
    reparameterized_flat_line = coordinate_shape(1 / (4 * coordinate))
    assert sp.simplify(reparameterized_flat_line + 2 / coordinate) == 0
    assert reparameterized_flat_line.subs(coordinate, 1) == -2
    # An intrinsic curvature would keep the flat line's zero under reparameterization;
    # the exact counterexample must stay live.
    assert reparameterized_flat_line.subs(coordinate, 1) != 0


# The one-dimensional Bures path coefficient g_path(CPsi) and its coordinate-shape
# second derivative, per host: (phrases the host states, curvature labels it must not use).
F47_HOSTS = {
    "experiments/COCKPIT_UNIVERSALITY.md": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative", "not Gaussian or intrinsic curvature"),
        ("Bures curvature", "curvature grows with system size"),
    ),
    "simulations/cockpit_n5.py": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative", "not intrinsic curvature",
         "g_path(CPsi)", "S_CPsi"),
        ("BURES CURVATURE",),
    ),
    "simulations/results/cockpit_n5.txt": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative", "not intrinsic curvature",
         "g_path(CPsi)", "S_CPsi"),
        ("K_Gauss",),
    ),
    "review/OPEN_QUESTIONS_INDEX.md": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative"),
        ("Gaussian curvature",),
    ),
    "experiments/INFORMATION_GEOMETRY.md": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative", "not a curvature"),
        ("Gaussian curvature", "hyperbolic"),
    ),
    "simulations/information_geometry.py": (
        ("Bures path-metric coefficient", "coordinate-shape second derivative", "not intrinsic curvature"),
        ("GAUSSIAN CURVATURE", "K (Gauss)"),
    ),
    "simulations/results/information_geometry.txt": (
        ("BURES PATH-METRIC COEFFICIENT", "coordinate-shape second derivative", "not intrinsic curvature"),
        ("GAUSSIAN CURVATURE", "K (Gauss)"),
    ),
}


def f47_findings(path, source):
    required, forbidden = F47_HOSTS[path]
    normalized = " ".join(source.split())
    found = [f"{path}: missing {phrase}" for phrase in required if phrase not in normalized]
    found += [f"{path}: curvature label {phrase}" for phrase in forbidden if phrase.lower() in normalized.lower()]
    return found


@pytest.mark.parametrize("path", sorted(F47_HOSTS))
def test_f47_host_names_the_one_dimensional_path_object(path):
    assert f47_findings(path, read_host(path)) == []


@pytest.mark.parametrize("path", sorted(F47_HOSTS))
def test_f47_curvature_label_restoration_fails_on_each_host(path):
    source = read_host(path)
    for stale in F47_HOSTS[path][1]:
        mutant = source.replace("coordinate-shape second derivative", stale, 1)
        assert mutant != source
        assert any(stale in item for item in f47_findings(path, mutant)), (path, stale)


def test_bures_geodesic_holds_for_the_commuting_bell_family_and_fails_where_h_moves_the_state():
    # Bell+ under equal local Z-dephasing stays p|Phi+><Phi+| + (1-p)|Phi-><Phi-|, a commuting
    # family. There the Bures angle is the Hellinger angle: with sqrt(p) = cos(a), the angle
    # between p and q is exactly |a(p) - a(q)|, so angles add along the monotone path and the
    # path is the geodesic (INFORMATION_GEOMETRY Phase 3, F46).
    a1, a2 = sp.symbols("a1 a2", real=True)
    fidelity_root = sp.cos(a1) * sp.cos(a2) + sp.sin(a1) * sp.sin(a2)
    assert sp.simplify(fidelity_root - sp.cos(a1 - a2)) == 0
    # ... and Bell+ really stays in that family: exactly, over the integers. The Heisenberg
    # bond acts as +1 on both Bell projectors, and local Z-dephasing maps each projector into
    # the span of the two (Z Phi+ Z = Phi-).
    sx = sp.Matrix([[0, 1], [1, 0]]); sy = sp.Matrix([[0, -sp.I], [sp.I, 0]]); sz = sp.diag(1, -1)
    i2 = sp.eye(2)
    kron = sp.kronecker_product
    heis = kron(sx, sx) + kron(sy, sy) + kron(sz, sz)
    phi_p = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2)
    phi_m = sp.Matrix([1, 0, 0, -1]) / sp.sqrt(2)
    proj_p, proj_m = phi_p * phi_p.T, phi_m * phi_m.T
    for proj in (proj_p, proj_m):
        assert heis * proj - proj * heis == sp.zeros(4, 4)
        for z in (kron(sz, i2), kron(i2, sz)):
            assert sp.simplify(z * proj * z - (proj_m if proj == proj_p else proj_p)) == sp.zeros(4, 4)
    # The counterexample that must stay live: |+0>, which the Heisenberg bond moves, runs a
    # path several times longer than the endpoint angle (7.68 in the producer).
    paulis = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0]).astype(complex)]
    ham = sum(np.kron(p, p) for p in paulis)
    ident = np.eye(4)
    gen = -1j * (np.kron(ham, ident) - np.kron(ident, ham.T))
    for z in (np.kron(paulis[2], np.eye(2)), np.kron(np.eye(2), paulis[2])):
        gen += 0.05 * (np.kron(z, z.conj()) - np.eye(16))

    def angle(r, q):
        w, v = np.linalg.eigh(r)
        root = v @ np.diag(np.sqrt(np.maximum(w, 0))) @ v.conj().T
        ev = np.linalg.eigvalsh(root @ q @ root)
        return float(np.arccos(min(1.0, np.sum(np.sqrt(np.maximum(ev, 0))))))

    psi = np.kron([1, 1], [1, 0]).astype(complex) / np.sqrt(2)
    states = [np.outer(psi, psi.conj())]
    step = expm(gen * (1.5 / 599))
    for _ in range(599):
        states.append((step @ states[-1].reshape(-1)).reshape(4, 4))
    length = sum(angle(states[k], states[k + 1]) for k in range(599))
    assert length / angle(states[0], states[-1]) > 5


# The geodesic and susceptibility labels INFORMATION_GEOMETRY once carried: a one-dimensional
# geodesic equation that holds for any monotone curve, and d2CPsi/dgamma2 called a Fisher
# susceptibility. Each host must not say either again.
GEODESIC_HOSTS = ("experiments/INFORMATION_GEOMETRY.md", "simulations/information_geometry.py",
                  "simulations/results/information_geometry.txt")
GEODESIC_STALE = ("approximately geodesic", "geodesic deviation", "Fisher susceptibility")


@pytest.mark.parametrize("path", GEODESIC_HOSTS)
def test_information_geometry_hosts_do_not_restore_the_one_dimensional_geodesic(path):
    text = " ".join(read_host(path).split()).lower()
    assert [w for w in GEODESIC_STALE if w.lower() in text] == []
    mutant = text + " the lindblad trajectory is approximately geodesic"
    assert [w for w in GEODESIC_STALE if w.lower() in mutant] == ["approximately geodesic"]


# ---------------------------------------------------------- V-Effect N=3 census
def test_v_effect_n3_distinct_pair_census_is_exact():
    simulations = str(ROOT / "simulations")
    added = simulations not in sys.path
    if added:
        sys.path.insert(0, simulations)
    try:
        import framework as fw
    finally:
        if added:
            sys.path.remove(simulations)
    labels = tuple(a + b for a in "XYZ" for b in "XYZ")
    routed = {f"{left}+{right}": fw.classify_two_term_palindrome(left, right)["fate"]
              for left, right in combinations(labels, 2)}
    assert len(routed) == 36
    assert {key for key, fate in routed.items() if fate == "truly"} == {"XX+YY", "XX+ZZ", "YY+ZZ"}
    assert [*routed.values()].count("truly") == 3
    assert [*routed.values()].count("soft") == 19
    assert [*routed.values()].count("hard") == 14


# ------------------------------------------ calibration producers: the R* proxy band
CALIBRATION_CONTRACTS = {
    "data/ibm_history/ibm_history_analysis.py": {
        "required": ("R_STAR = 0.21275477982200533", "normalized-purity proxy", "proxy_band",
                     "normalized_purity_proxy_analysis.png", "normalized_purity_proxy_at_time",
                     "find_proxy_quarter_crossing", "proxy = 1/4"),
        "forbidden": ("Purity crosses C*Psi = 1/4", "NEVER reaches 1/4", "r* (kritisch)",
                      "quarter_boundary_analysis.png", "C·Ψ = ¼", "below_rstar = r_vals < R_STAR",
                      "pct = 100*(rv<R_STAR).mean()"),
    },
    "simulations/ptf_clock_field.py": {
        "required": ("normalized-purity proxy", "both proxy bands"),
        "forbidden": ("Q3 (the quarter-boundary field)", "Q3 quarter boundary"),
    },
    "simulations/marrakesh_quarter_boundary_review.py": {
        "required": ("R_STAR = 0.21275477982200533", "is_below_rstar", "below-R*", "at-or-above-R*"),
        "forbidden": ("R_THRESHOLD = 0.213", "quantum-side", "classical-side", "pure-quantum",
                      "specific physical threshold", "flip_to_quantum", "stable_quantum"),
    },
    "simulations/marrakesh_uniform_quantum_chain.py": {
        "required": ("below_rstar_stable", "below_fraction", "below %", "at-or-above-R*"),
        "forbidden": ("quantum_stable", "uniform-quantum", "uniform-classical", "quantum-side", "\"crossing\":",
                      "cross %"),
    },
    "simulations/marrakesh_kingston_fez_compare.py": {
        "required": ("below_rstar_rows", "below-R*", "at-or-above-R*"),
        "forbidden": ("stable_below_rstar", "r_threshold", "stable_quantum", "uniform-quantum", "uniform-classical"),
    },
    "simulations/qubit_biography.py": {
        "required": ("R_STAR = 0.21275477982200533", "below_flags = rs < R_STAR", "stable-below",
                     "stable-at-or-above", "r = t2 / (2.0 * t1)"),
        "forbidden": ("np.sign(rs - R_STAR)", "pulse-stable", "silent-stable", "classic-stable", "quantum-side",
                      "classical-side"),
    },
    "simulations/marrakesh_may05_preflight.py": {
        "required": ("empirical composite heuristic", "volatility penalty", "variable-near penalty"),
        "forbidden": ("volatility veto", "twitch penalty", "4. RECOMMENDATION", "uniform-quantum",
                      "uniform-classical"),
    },
    "simulations/marrakesh_path_biography.py": {
        "required": ("normalized-purity proxy", "R* bands", "below_fraction", "below %"),
        "forbidden": ("CΨ regime identity", "crossing =", "cross %"),
    },
    "simulations/f88b_lens_ibm_kingston_uniform_quantum.py": {
        "required": ("predominantly below-R*", "at-or-above-R*", "confounded association"),
        "forbidden": ("all below-R*", "all below R*", "all-below-R*", "uniform-quantum", "uniform-classical",
                      "quantum-side dephasing dominates"),
    },
    "simulations/kingston_f97_lens.py": {
        "required": ("radial diagnostic", "not the F97 parameter", "c=+1/4"),
        "forbidden": ("cross the cardioid boundary", "cardioid attractor", "F97-anchor crossings", "Quarter cusp",
                      "Half cardioid"),
    },
}


def calibration_findings(path, text):
    contract = CALIBRATION_CONTRACTS[path]
    lowered = text.lower()
    found = [f"{path}: missing positive control: {phrase}"
             for phrase in contract["required"] if phrase.lower() not in lowered]
    found += [f"{path}: restored unsupported label: {phrase}"
              for phrase in contract["forbidden"] if phrase.lower() in lowered]
    return tuple(found)


@pytest.mark.parametrize("path", sorted(CALIBRATION_CONTRACTS))
def test_calibration_producer_labels_are_current(path):
    source = read_host(path)
    ast.parse(source, filename=path)
    assert calibration_findings(path, source) == ()


@pytest.mark.parametrize("path", sorted(CALIBRATION_CONTRACTS))
def test_calibration_producer_label_restoration_fails(path):
    original = read_host(path)
    for phrase in CALIBRATION_CONTRACTS[path]["forbidden"]:
        mutant = original + "\n# " + phrase + "\n"
        assert any(phrase in item for item in calibration_findings(path, mutant)), (path, phrase)


HISTORY_PATH = "data/ibm_history/ibm_history_analysis.py"
HISTORY_LEGACY_KEYS = {"r_param", "cpsi_min", "distance_from_quarter", "crosses_quarter", "t_star_us",
                       "t_star_over_T2", "t_coherence_crossing_us", "t_coh_over_T2"}
HISTORY_NEUTRAL_KEYS = {"r", "proxy_band", "normalized_purity_proxy_min", "proxy_offset_from_quarter",
                        "is_below_rstar", "proxy_quarter_time_us", "proxy_quarter_time_over_T2"}


def _function(tree, name):
    return next((node for node in ast.walk(tree)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def history_schema_findings(text):
    tree = ast.parse(text)
    found = []
    compute = _function(tree, "compute_qubit_record")
    if compute is None:
        found.append("missing compute_qubit_record")
    else:
        keys = {key.value for candidate in ast.walk(compute) if isinstance(candidate, ast.Dict)
                for key in candidate.keys if isinstance(key, ast.Constant) and isinstance(key.value, str)}
        if HISTORY_NEUTRAL_KEYS - keys:
            found.append(f"neutral output keys missing: {sorted(HISTORY_NEUTRAL_KEYS - keys)}")
        if HISTORY_LEGACY_KEYS & keys:
            found.append(f"legacy output keys remain active: {sorted(HISTORY_LEGACY_KEYS & keys)}")
    names = {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    retired = {"cpsi_from_purity", "cpsi", "find_quarter_crossing", "normalized_purity_proxy",
               "find_proxy_quarter_time"} & names
    if retired:
        found.append(f"legacy active functions remain: {sorted(retired)}")
    adapter = _function(tree, "_adapt_legacy_calibration_schema")
    if adapter is None:
        found.append("missing legacy calibration schema adapter")
    elif "LEGACY-CALIBRATION-SCHEMA" not in (ast.get_docstring(adapter) or ""):
        found.append("legacy adapter is not visibly marked")
    for node in tree.body:
        if (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == "mkdir"):
            found.append("import-time directory creation remains")
    return tuple(found)


def test_history_writer_uses_neutral_output_schema_and_marked_legacy_adapter():
    assert history_schema_findings(read_host(HISTORY_PATH)) == ()


def test_history_schema_gate_kills_active_legacy_output_and_unmarked_adapter():
    text = read_host(HISTORY_PATH)
    active_legacy = text.replace("'proxy_band': band,", "'crosses_quarter': band,", 1)
    assert active_legacy != text
    assert any("legacy output keys" in item for item in history_schema_findings(active_legacy))
    unmarked = text.replace("LEGACY-CALIBRATION-SCHEMA", "legacy calibration schema", 1)
    assert "legacy adapter is not visibly marked" in history_schema_findings(unmarked)
    mkdir = text + "\nOUTPUT_DIR.mkdir(exist_ok=True)\n"
    assert "import-time directory creation remains" in history_schema_findings(mkdir)


EXACT_RSTAR_HOSTS = ("data/ibm_history/ibm_history_analysis.py", "simulations/ptf_clock_field.py",
                     "simulations/marrakesh_quarter_boundary_review.py", "simulations/qubit_biography.py")

ACTIVE_RSTAR_BOUNDARIES = {
    "data/ibm_history/ibm_history_analysis.py": ("proxy_band", 'return "below-R*" if r < R_STAR else "at-or-above-R*"'),
    "simulations/marrakesh_quarter_boundary_review.py": ("is_below_rstar", "return r_param(q) < R_STAR"),
    "simulations/marrakesh_uniform_quantum_chain.py": ("main", "below_flags = rs < R_STAR"),
    "simulations/marrakesh_kingston_fez_compare.py": ("below_rstar_rows", "(q.t2_us / (2 * q.t1_us)) < R_STAR"),
    "simulations/qubit_biography.py": ("archetype_from_series", "below_flags = rs < R_STAR"),
    "simulations/marrakesh_may05_preflight.py": ("band_count", "by_id[qid].t2_us / (2.0 * by_id[qid].t1_us) < R_STAR"),
    "simulations/marrakesh_path_biography.py": ("main", "below_flags = rs < R_STAR"),
}


def assigned_numeric_constant(text, name):
    for node in ast.parse(text).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, (int, float)):
                return float(node.value.value)
    return None


def active_rstar_comparisons(text, function_name):
    owner = _function(ast.parse(text), function_name)
    if owner is None:
        return ()
    comparisons = []
    for node in ast.walk(owner):
        if isinstance(node, ast.Compare) and len(node.ops) == 1 and any(
                isinstance(item, ast.Name) and item.id == "R_STAR" for item in (node.left, *node.comparators)):
            comparisons.extend(type(op) for op in node.ops)
    return tuple(comparisons)


@pytest.mark.parametrize("path", EXACT_RSTAR_HOSTS)
def test_rstar_constant_is_the_exact_root_and_mutation_sensitive(path):
    text = read_host(path)
    assert assigned_numeric_constant(text, "R_STAR") == 0.21275477982200533
    mutant = text.replace("0.21275477982200533", "0.212754779822", 1)
    assert mutant != text
    assert assigned_numeric_constant(mutant, "R_STAR") != 0.21275477982200533


@pytest.mark.parametrize("path", sorted(ACTIVE_RSTAR_BOUNDARIES))
def test_rstar_band_is_strictly_below_and_rejects_equality_mutation(path):
    text = read_host(path)
    function_name, fragment = ACTIVE_RSTAR_BOUNDARIES[path]
    assert fragment in text, path
    ops = active_rstar_comparisons(text, function_name)
    assert ast.Lt in ops, path
    assert ast.LtE not in ops, path
    mutant = text.replace(fragment, fragment.replace("< R_STAR", "<= R_STAR", 1), 1)
    assert mutant != text
    assert ast.LtE in active_rstar_comparisons(mutant, function_name), path


def test_history_record_at_exactly_rstar_is_at_or_above_and_legacy_flag_is_carried(tmp_path):
    source = read_host(HISTORY_PATH)
    names = {"purity_single_qubit", "normalized_purity_proxy_from_purity", "normalized_purity_proxy_at_time",
             "find_proxy_minimum", "find_proxy_quarter_crossing", "coherence_one_quarter_time", "r_parameter",
             "proxy_band", "compute_qubit_record", "_adapt_legacy_calibration_schema"}
    module = source_copy_definitions(
        source, names, tmp_path, assignments=("R_STAR", "LN4"),
        prelude="import numpy as np\nfrom scipy.optimize import brentq, minimize_scalar\n")
    r_star = module["R_STAR"]
    assert r_star == 0.21275477982200533
    record = module["compute_qubit_record"]("fixture", 1, 0.5, r_star)
    assert record["r"] == r_star  # T2/(2 T1) with T1 = 1/2 is T2 itself: exact
    assert record["proxy_band"] == "at-or-above-R*"
    assert record["is_below_rstar"] is False
    assert HISTORY_LEGACY_KEYS.isdisjoint(record)

    legacy = {"date": "fixture", "qubit": "1", "T1_us": "0.5", "T2_us": repr(r_star), "frequency_GHz": "5.0",
              "r_param": "0.212755", "cpsi_min": "0.25", "distance_from_quarter": "0.0", "crosses_quarter": "True",
              "t_star_us": "0.1", "t_star_over_T2": "0.2", "t_coherence_crossing_us": "0.3",
              "t_coh_over_T2": "1.386294"}
    adapted = module["_adapt_legacy_calibration_schema"](legacy)
    # The archive rounded T1/T2 after computing its flag, so the adapter carries the
    # archived raw-derived flag instead of reclassifying the rounded display columns.
    assert adapted["proxy_band"] == "below-R*"
    assert adapted["is_below_rstar"] is True
    assert adapted["band_provenance"] == "legacy-raw-derived-flag"
    assert HISTORY_LEGACY_KEYS.isdisjoint(adapted)


# ------------------------------------------------ import safety of a producer
VEFFECT_PRODUCER = "simulations/veffect_cavity_modes.py"
VEFFECT_RESULT = "simulations/results/veffect_cavity_modes.txt"
IMPORT_PROBE = (
    "import importlib.util, sys\n"
    "spec = importlib.util.spec_from_file_location('import_probe', sys.argv[1])\n"
    "module = importlib.util.module_from_spec(spec)\n"
    "spec.loader.exec_module(module)\n"
    "sys.stdout.write('MAIN ' + str(callable(getattr(module, 'main', None))))\n"
)


def import_in_subprocess(script: Path, watched: Path, cwd: Path):
    """Import `script` (not as __main__) in a fresh interpreter and report whether the
    watched file kept its bytes. A changed file is restored before returning, so a
    firing gate leaves the tracked result as it found it."""
    before = watched.read_bytes()
    try:
        completed = subprocess.run(
            [sys.executable, "-c", IMPORT_PROBE, str(script)], cwd=cwd, capture_output=True, text=True,
            encoding="utf-8", timeout=300, check=False,
            env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONDONTWRITEBYTECODE": "1"})
    finally:
        unchanged = watched.read_bytes() == before
        if not unchanged:
            watched.write_bytes(before)
    return completed, unchanged


def unguarded_main_calls(source):
    """Module-level calls of main() outside an `if __name__ == "__main__":` block."""
    guards = ("__name__ == '__main__'", "'__main__' == __name__")
    calls = []
    for node in ast.parse(source).body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if isinstance(node, ast.If) and ast.unparse(node.test) in guards:
            continue
        calls += [ast.unparse(sub) for sub in ast.walk(node)
                  if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name) and sub.func.id == "main"]
    return calls


def test_veffect_producer_calls_main_only_under_the_main_guard():
    source = read_host(VEFFECT_PRODUCER)
    assert unguarded_main_calls(source) == []
    assert 'if __name__ == "__main__":\n    main()' in source
    unguarded = source.replace('if __name__ == "__main__":\n    main()', "main()", 1)
    assert unguarded != source
    assert unguarded_main_calls(unguarded) == ["main()"]


def test_veffect_import_is_silent_and_leaves_its_result_bytes(tmp_path):
    completed, unchanged = import_in_subprocess(ROOT / VEFFECT_PRODUCER, ROOT / VEFFECT_RESULT, tmp_path)
    assert completed.returncode == 0, completed.stderr[-3000:]
    assert unchanged, f"importing {VEFFECT_PRODUCER} rewrote {VEFFECT_RESULT} (restored by the gate)"
    assert completed.stdout == "MAIN True", completed.stdout[-500:]
    assert completed.stderr == ""
    assert list(tmp_path.iterdir()) == []


def test_import_probe_catches_a_writing_and_printing_import(tmp_path):
    workdir = tmp_path / "cwd"
    workdir.mkdir()
    result = tmp_path / "fake_result.txt"
    result.write_bytes(b"original\n")
    fake = tmp_path / "fake_producer.py"
    fake.write_text("from pathlib import Path\n"
                    "Path(__file__).with_name('fake_result.txt').write_text('overwritten\\n')\n"
                    "print('loud')\n"
                    "def main():\n    pass\n", encoding="utf-8")
    completed, unchanged = import_in_subprocess(fake, result, workdir)
    assert completed.returncode == 0, completed.stderr
    assert not unchanged
    assert result.read_bytes() == b"original\n"
    assert completed.stdout != "MAIN True"
