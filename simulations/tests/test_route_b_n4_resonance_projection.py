"""What each half of the generator supplies at the N=4 end-weight profile.

The certificate is read, never re-executed: the numbers below are re-derived
from parts(4) and compared against the committed artifact. Each check carries
the mutation or the control that makes its own instrument read the other
answer, because an instrument that cannot report the other answer reports
nothing.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'simulations'))

import pytest
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM

from framework.weight_coherence_block import weight_block_configs
from route_b_artifact_provenance import artifact_provenance, canonical_text_sha256
from route_b_other_n_unfolding import parts

SCRIPT = ROOT/'simulations/route_b_n4_resonance_projection.py'
ARTIFACT = ROOT/'simulations/results/route_b_n4_resonance_projection.json'
DIMENSION = 24
RESONANT = 5


@pytest.fixture(scope='module')
def certificate():
    return json.loads(ARTIFACT.read_text(encoding='utf-8'))


@pytest.fixture(scope='module')
def block():
    """The (1,2) coherence block rebuilt independently of the producer."""
    field = s.QQ.algebraic_field(s.sqrt(3))
    d, c, left, right, _ = parts(4)

    def integer(a):
        # The producer's own guard, kept: a silent truncation here would make
        # the independent rebuild weaker than the thing it rebuilds.
        assert all(v.imag == 0 and int(v.real) == v.real for row in a for v in row)
        return s.Matrix([[int(v.real) for v in row] for row in a])

    matrix = lambda a: DM.from_Matrix(integer(a)).convert_to(field)
    bulk, ends = integer((c-left-right)/1j), integer((left+right)/1j)
    dissipator = matrix(d)
    hop = DM.from_Matrix(bulk+s.sqrt(3)/2*ends).convert_to(field).to_dense()
    identity = DM.eye(DIMENSION, field).to_dense()
    resonant = (hop-identity).nullspace().transpose().to_dense()
    labels = [(a, b) for a in weight_block_configs(4, 1) for b in weight_block_configs(4, 2)]
    return dict(field=field, dissipator=dissipator, bulk=bulk, ends=ends, hop=hop,
                identity=identity, resonant=resonant,
                hamming=[bin(a ^ b).count('1') for a, b in labels])


def compress(space, operator):
    """The compression onto a subspace, in that subspace's own basis."""
    columns = space.to_Matrix() if hasattr(space, 'to_Matrix') else space
    acting = operator.to_Matrix() if hasattr(operator, 'to_Matrix') else operator
    dual = (columns.T*columns).inv()*columns.T
    return s.simplify(dual*acting*columns)


def multiplicities(operator, field, size):
    """Eigenvalue multiplicities as exact ranks: no eigensolver, no tolerance.

    The sum can only UNDERCOUNT, so reaching `size` proves both that the
    scanned window holds the whole spectrum and that the operator is
    semisimple. There is no path to a false positive.
    """
    identity = DM.eye(size, field).to_dense()
    found = {k: size-(operator-identity.scalarmul(field.convert(k))).rank()
             for k in range(-size, size+1)}
    return {k: m for k, m in found.items() if m}


def test_certificate_provenance_is_canonical_and_survives_checkout_line_endings(certificate):
    assert certificate['script_sha256'] == canonical_text_sha256(SCRIPT)
    expected = artifact_provenance(SCRIPT, SCRIPT)
    expected.pop('source_sha256')
    # The dependency closure is pinned too: the labelling this certificate
    # rests on comes from framework.weight_coherence_block, not from the
    # script alone.
    assert 'simulations/framework/weight_coherence_block.py' in expected['dependencies']
    assert 'simulations/route_b_other_n_unfolding.py' in expected['dependencies']
    for field, value in expected.items():
        assert certificate[field] == value, f'provenance mismatch: {field}'
    lf = SCRIPT.read_bytes().replace(b'\r\n', b'\n').replace(b'\r', b'\n')
    assert b'\n' in lf
    raw = set()
    for ending in (b'\n', b'\r\n', b'\r'):
        checkout = lf.replace(b'\n', ending)
        canonical = checkout.decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')
        assert hashlib.sha256(canonical.encode('utf-8')).hexdigest() == certificate['script_sha256']
        raw.add(hashlib.sha256(checkout).hexdigest())
    # The mutation: a raw-bytes book gives three different answers for one
    # file, so a hash kept that way is unreachable on any checkout but the one
    # that wrote it.
    assert len(raw) == 3


def test_dephasing_half_is_the_absorption_diagonal_and_its_two_rungs_are_mirror_partners(block, certificate):
    hamming = block['hamming']
    diagonal = s.Matrix(block['dissipator'].to_Matrix())
    assert diagonal == s.diag(*[-2*k for k in hamming])
    strata = sorted(set(hamming))
    assert strata == [1, 3]
    assert certificate['hamming_strata'] == {str(k): hamming.count(k) for k in strata}
    assert certificate['absorption_rate_per_stratum'] == {str(k): -2*k for k in strata}
    # The two rungs pair about −4γ because their POPULATIONS balance, and the
    # control is the same block one site up, where the values are still 1 and 3
    # but the populations are 20 and 30, so the same reading gives no pairing.
    assert hamming.count(1) == hamming.count(3)
    fifth = [bin(a ^ b).count('1') for a in weight_block_configs(5, 1)
             for b in weight_block_configs(5, 2)]
    assert sorted(set(fifth)) == [1, 3]
    assert [fifth.count(k) for k in (1, 3)] == [20, 30] != [25, 25]


def test_coherent_half_is_symmetric_at_every_end_weight(block, certificate):
    weight = s.symbols('w', real=True)
    symbolic = block['bulk']+weight*block['ends']
    assert (symbolic-symbolic.T).is_zero_matrix
    assert certificate['coherent_half_symmetric_at_symbolic_weight'] is True
    # The control: the same reading on a perturbation that breaks symmetry
    # reports it, so the emptiness above is a measurement.
    broken = symbolic.copy()
    broken[0, 1] += 1
    assert not (broken-broken.T).is_zero_matrix


def test_the_degeneracy_is_free_fermion_and_the_equal_spacing_only_sets_its_size(block, certificate):
    weight = s.symbols('w', real=True)
    factors = s.factor_list((block['bulk']+weight*block['ends']).charpoly().as_expr())[1]
    # Over Q(w), so at EVERY end weight and not at sampled ones: two quadratic
    # factors to the fourth power.
    observed = sorted(((s.Poly(f, s.Symbol('lambda')).degree(), m) for f, m in factors), reverse=True)
    assert [list(pair) for pair in observed] == certificate['free_fermion_factor_multiplicities_over_Qw']
    assert [m for _, m in observed] == [4, 4, 1, 1, 1, 1]
    # The control: a matrix with no repeated factor reads all ones, so a four
    # is a degeneracy and not the instrument's habit.
    simple = s.diag(1, 2, 3)
    assert [m for _, m in s.factor_list(simple.charpoly().as_expr())[1]] == [1, 1, 1]
    comb = multiplicities(block['hop'], block['field'], DIMENSION)
    assert {str(k): m for k, m in sorted(comb.items())} == certificate['coherent_half_frequency_comb']
    assert sum(comb.values()) == DIMENSION == certificate['coherent_half_geometric_total']
    assert comb[1] == RESONANT == max(comb.values())
    # Four at every weight, five only here: the value +1 lies on the
    # multiplicity-four factor at this end weight alone, and there a simple
    # factor passes through it as well.
    here = s.sqrt(3)/2
    at_resonance = [(m, s.solve(s.Eq(s.simplify(f.subs(s.Symbol('lambda'), 1)), 0), weight))
                    for f, m in factors]
    assert sum(1 for m, roots in at_resonance if m == 4 and here in roots) == 1
    assert sum(1 for m, roots in at_resonance if m == 1 and here in roots) == 1
    assert sorted({str(root) for m, roots in at_resonance if m == 4
                   for root in roots if root.is_positive})         == certificate['end_weight_lifting_the_multiplet_to_five'] == ['sqrt(3)/2']
    assert certificate['resonance_is_four_plus_one'] is True
    # The control: the same instrument on a Jordan block falls short of the
    # dimension, so reaching it is semisimplicity and not arithmetic.
    defective = DM.from_Matrix(s.Matrix([[1, 1], [0, 1]])).convert_to(block['field']).to_dense()
    assert sum(multiplicities(defective, block['field'], 2).values()) == 1


def test_the_resonant_eigenspace_meets_each_stratum_which_straddling_alone_would_not_give(block, certificate):
    hamming, strata, field = block['hamming'], [1, 3], block['field']
    columns = block['resonant'].to_Matrix()
    ranks = [columns[[i for i in range(DIMENSION) if hamming[i] == k], :].rank() for k in strata]
    assert ranks == [4, 4]
    assert certificate['resonance_rank_per_stratum'] == {str(k): r for k, r in zip(strata, ranks)}
    meets = {str(k): len(columns[[i for i in range(DIMENSION) if hamming[i] != k], :].nullspace())
             for k in strata}
    assert meets == {'1': 1, '3': 1} == certificate['resonance_stratum_intersections']
    # The control on that instrument: a subspace sitting inside stratum 1 meets
    # stratum 3 in dimension zero, so a one is a measurement.
    inside = s.zeros(DIMENSION, 3)
    for column, row in enumerate([i for i in range(DIMENSION) if hamming[i] == 1][:3]):
        inside[row, column] = 1
    assert [len(inside[[i for i in range(DIMENSION) if hamming[i] != k], :].nullspace())
            for k in strata] == [3, 0]
    # The control that matters: straddling both strata is NOT by itself
    # non-scalarity. This plane straddles and compresses to exactly −4·I.
    balanced = s.zeros(DIMENSION, 2)
    for column, (one, three) in enumerate(zip([i for i in range(DIMENSION) if hamming[i] == 1][:2],
                                              [i for i in range(DIMENSION) if hamming[i] == 3][:2])):
        balanced[one, column] = balanced[three, column] = 1
    plane = DM.from_Matrix(balanced).convert_to(field).to_dense()
    assert [balanced[[i for i in range(DIMENSION) if hamming[i] == k], :].rank() for k in strata] == [2, 2]
    assert compress(plane, block['dissipator']) == -4*s.eye(2)
    rates = compress(block['resonant'], block['dissipator']).eigenvals()
    assert len(rates) == sum(rates.values()) == RESONANT
    assert {s.sympify(k) for k in certificate['restricted_dissipator_rates']} == set(rates)
    assert set(rates) >= {s.Integer(-2), s.Integer(-6)}, "the strata's extreme rates themselves"
    assert s.simplify(sum(k*m for k, m in rates.items())/RESONANT) == certificate['restricted_dissipator_mean'] == -4


def test_every_multiplet_splits_evenly_between_the_strata_so_minus_four_is_the_blocks_own(block, certificate):
    field, recorded = block['field'], certificate['restricted_mean_per_multiplet']
    hamming = block['hamming']
    projector = s.diag(*[1 if k == 1 else 0 for k in hamming])
    means, weights = {}, {}
    for k in [int(key) for key in recorded]:
        space = (block['hop']-block['identity'].scalarmul(field.convert(k))).nullspace().transpose().to_dense()
        assert space.shape[1] > 0
        means[str(k)] = s.simplify(compress(space, block['dissipator']).trace()/space.shape[1])
        weights[str(k)] = s.simplify(compress(space, projector).trace()/space.shape[1])
    assert means == {k: s.Integer(v) for k, v in recorded.items()}
    assert set(means.values()) == {s.Integer(-4)}
    # The reason, and it is the block's and not this eigenspace's: every
    # multiplet carries exactly half its weight in each stratum.
    assert weights == {k: s.sympify(v) for k, v in certificate['stratum_weight_per_multiplet'].items()}
    assert set(weights.values()) == {s.Rational(1, 2)}
    # The control on the weight instrument itself: a subspace that does NOT
    # split evenly reads its own fraction, so one half is a measurement.
    lopsided = s.zeros(DIMENSION, 4)
    for column, row in enumerate([i for i in range(DIMENSION) if hamming[i] == 1][:3]
                                 +[i for i in range(DIMENSION) if hamming[i] == 3][:1]):
        lopsided[row, column] = 1
    tilted_space = DM.from_Matrix(lopsided).convert_to(field).to_dense()
    assert s.simplify(compress(tilted_space, projector).trace()/4) == s.Rational(3, 4)
    # The control: a dephasing that is NOT a function of Hamming breaks the
    # constancy, so the same reading does discriminate. A Hamming-indexed one
    # would not, whatever its two values, which is the point.
    tilted = s.Matrix(block['dissipator'].to_Matrix())
    tilted[0, 0] -= 7
    drifting = set()
    for k in [int(key) for key in recorded]:
        space = (block['hop']-block['identity'].scalarmul(field.convert(k))).nullspace().transpose().to_dense()
        drifting.add(s.simplify(compress(space, tilted).trace()/space.shape[1]))
    assert len(drifting) > 1


def test_flattening_the_dephasing_un_resolves_the_coincidence_rather_than_removing_it(block, certificate):
    field = block['field']
    hopping = compress(block['resonant'], block['ends'])
    u, v = s.symbols('u v', real=True)
    flattened = s.factor((s.I*v*s.eye(RESONANT)-s.I*u/(4*s.sqrt(3))*hopping).det(method='domain-ge'))
    assert s.expand(flattened-s.I*(u-4*v)**4*(u+4*v)/1024) == 0
    assert str(flattened) == certificate['scalar_dephasing_determinant']
    shift = (hopping-s.sqrt(3)*s.eye(RESONANT)).applyfunc(s.simplify)
    geometric = RESONANT-DM.from_Matrix(shift).convert_to(field).rank()
    algebraic = RESONANT-DM.from_Matrix((shift**RESONANT).applyfunc(s.simplify)).convert_to(field).rank()
    # Fourfold and SEMISIMPLE: the flattened generator keeps a degenerate line,
    # it does not lose the coincidence.
    assert algebraic == geometric == 4
    assert certificate['scalar_dephasing_line'] == dict(eigenvalue='sqrt(3)', algebraic=algebraic, geometric=geometric)
    # The control: the same pair of ranks on a Jordan block separates, so
    # equality is semisimplicity and not an identity of the instrument.
    jordan = s.Matrix([[1, 1], [0, 1]])
    step = (jordan-s.eye(2))
    assert 2-DM.from_Matrix(step).convert_to(field).rank() == 1
    assert 2-DM.from_Matrix(step**2).convert_to(field).rank() == 2


def test_the_four_points_carry_a_double_root_of_geometric_multiplicity_one(block, certificate):
    dephasing = compress(block['resonant'], block['dissipator'])
    hopping = compress(block['resonant'], block['ends'])
    u, v = s.symbols('u v', real=True)
    effective = dephasing+s.I*u/(4*s.sqrt(3))*hopping
    polynomial = s.factor(((-4+s.I*v)*s.eye(RESONANT)-effective).det(method='domain-ge'))
    residual = (u+4*v)*(u-4*v)**2+40*v-8*u
    # The DOUBLE root lives in this factorization; rank four below is only its
    # geometric multiplicity. Neither alone is the Jordan block.
    assert s.expand(polynomial-s.I*((u-4*v)**2+64)*residual/1024) == 0
    assert str(polynomial) == certificate['effective_determinant']
    for row in certificate['effective_ranks']:
        uu, vv = s.sympify(row['u']), s.sympify(row['v'])
        numberfield = s.QQ.algebraic_field(s.I, uu, s.sqrt(3))
        shifted = ((-4+s.I*vv)*s.eye(RESONANT)-effective.subs(u, uu)).applyfunc(s.simplify)
        assert DM.from_Matrix(shifted).convert_to(numberfield).rank() == row['rank'] == RESONANT-1
        plane = s.Matrix.hstack(*(shifted*shifted).applyfunc(s.simplify).nullspace()).applyfunc(s.simplify)
        assert plane.shape == (RESONANT, 2), 'the eigenvalue is not semisimple'
        # Exact rank, never simplify(det) != 0: a structural inequality reads
        # an unreduced zero as nonzero.
        control = (s.I*vv*s.eye(RESONANT)-s.I*uu/(4*s.sqrt(3))*hopping).applyfunc(s.simplify)
        assert DM.from_Matrix(control).convert_to(numberfield).rank() == row['scalar_dephasing_control_rank'] == RESONANT


def test_the_coalescing_plane_straddles_two_multiplets_so_the_descent_premise_fails(block, certificate):
    dephasing = compress(block['resonant'], block['dissipator'])
    hopping = compress(block['resonant'], block['ends'])
    u = s.symbols('u', real=True)
    effective = dephasing+s.I*u/(4*s.sqrt(3))*hopping
    # The lemma needs both coalescing directions to descend from ONE multiplet.
    # T' sees the four-plus-one split as its own spectrum, so containment in
    # one of its eigenspaces is exactly that premise.
    slopes = s.simplify(hopping).eigenvals()
    assert slopes == {s.sqrt(3): 4, -s.sqrt(3): 1}
    assert {str(k): m for k, m in sorted(slopes.items(), key=str)} == certificate['end_weight_term_slopes']
    eigenspaces = {str(k): s.Matrix.hstack(*(hopping-k*s.eye(RESONANT)).applyfunc(s.simplify).nullspace())
                   for k in slopes}
    recorded = {row['u']: row for row in certificate['effective_term_plane_ranks']}
    for row in certificate['effective_ranks']:
        uu, vv = s.sympify(row['u']), s.sympify(row['v'])
        numberfield = s.QQ.algebraic_field(s.I, uu, s.sqrt(3))
        shifted = ((-4+s.I*vv)*s.eye(RESONANT)-effective.subs(u, uu)).applyfunc(s.simplify)
        plane = s.Matrix.hstack(*(shifted*shifted).applyfunc(s.simplify).nullspace()).applyfunc(s.simplify)

        def joint(operator):
            columns = s.Matrix.hstack(plane, (operator*plane).applyfunc(s.simplify))
            return DM.from_Matrix(columns).convert_to(numberfield).rank()

        half = (s.I*uu/(4*s.sqrt(3))*hopping).applyfunc(s.simplify)
        # Rank above two is the plane leaving itself: neither term of the
        # first-order operator preserves it.
        assert joint(half) == recorded[row['u']]['end_weight_term'] == 3
        assert joint(dephasing) == recorded[row['u']]['dephasing_term'] == 3
        # The premise itself: a plane inside a multiplet would join that
        # multiplet at its own dimension. Both joints exceed it.
        for slope, space in eigenspaces.items():
            columns = s.Matrix.hstack(space, plane)
            rank = DM.from_Matrix(columns).convert_to(numberfield).rank()
            assert [rank, space.cols] == recorded[row['u']]['multiplet_containment'][slope]
            assert rank > space.cols
        # The control: a plane genuinely inside the fourfold multiplet joins it
        # at four, so exceeding is straddling and not the instrument's floor.
        inside = eigenspaces['sqrt(3)'][:, :2]
        columns = s.Matrix.hstack(eigenspaces['sqrt(3)'], inside)
        assert DM.from_Matrix(columns).convert_to(numberfield).rank() == 4
        # The control: a plane the operator DOES preserve reads two, so three
        # is non-invariance and not the instrument's floor.
        vectors = [vector for _, _, basis in dephasing.eigenvects() for vector in basis][:2]
        invariant = s.Matrix.hstack(*vectors).applyfunc(s.simplify)
        columns = s.Matrix.hstack(invariant, (dephasing*invariant).applyfunc(s.simplify))
        assert DM.from_Matrix(columns).convert_to(numberfield).rank() == 2
