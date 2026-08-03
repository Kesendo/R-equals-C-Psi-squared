"""F138 clause 2, decided over several primes: no eigensolver, no tolerance.

F138 states the boundary law of the dephasing palindrome, and grades its two
clauses differently: clause 1 exhaustively swept, clause 2 spot-checked. This
sweeps clause 2.

THE ROUTE. The palindrome is spec(L) = {-lambda - s}, so with p(x) = det(xI-L)
and 4^N even it is exactly the polynomial identity p(x) == p(-x - s). Every
entry of L is a Gaussian rational here (gamma = 1/20, magnitudes 30/100,
22/100, 41/100, 17/100, 35/100, J and delta integers), so L maps into GF(p)
with i a square root of -1 and the identity is tested by comparing
determinants. An earlier version of this sweep instead ran an eigensolver and
built a calibrated tolerance on the premise that no exact route existed. The
premise was false, and CLAUDE.md's rounding rule asks for exactly the question
that would have caught it: why is there anything to tolerate?

WHAT THE VERDICTS ARE WORTH. Asymmetric, and the asymmetry is the point:

  BREAKS ARE PROVED. One witness x with p(x) != p(-x-s) mod any prime rules the
  identity out over Q(i) as well.

  HOLDS ARE CERTIFIED, NOT DECIDED. Two escape channels: Schwartz-Zippel (two
  distinct degree-n polynomials agree at at most n of the p elements) and the
  modulus itself (q = p(x) - p(-x-s) can be nonzero over Q(i) yet vanish
  identically mod one prime, which no number of points can see; the witness is
  guard 3 of f138_exact_palindrome_test.py). Hence THREE independent primes and
  two points each: a false HOLDS must escape in EVERY prime, and to escape
  structurally it would need every coefficient of q divisible by all three.
  Per prime the Schwartz-Zippel bound is (n/p)^2, so with the smallest prime
  used it is 1.9e-14 at n = 64 and 3.0e-13 at n = 256, and a row survives all
  three only below 6e-42 and 3e-38 respectively. Stated per n because the
  dimension differs by stage: Stage B runs at n = 256, not 64.

  The asymmetry runs the safe way for the headline. The load-bearing result is
  'no false positives', i.e. no row where the predicate says HOLD and the
  spectrum breaks. Those are the PROVED verdicts.

SCOPE HELD FIXED, and named so it can be read: J = 1 and delta = 1 on every
bond and gamma equal on every dephasing site, except in Stage E, which repeats
the N=3 grid with non-uniform rational J and gamma because F138 quantifies
over those and a law tested at one point of an axis is not tested on it.
"""

import itertools
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('\\', 1)[0].rsplit('/', 1)[0])

from f138_exact_palindrome_test import (PRIMES, det_mod_np_p,  # noqa: E402
                                     inv_p, sqrt_minus_one)

LETTERS = (0, 1, 2, 3)
LETTER_NAME = {0: '.', 1: 'X', 2: 'Y', 3: 'Z'}
POINTS = 2

GAMMA = (1, 20)
FIELD_NUM, FIELD_DEN = (30, 22, 41, 17, 35), 100

rng = np.random.default_rng(20260803)


def pauli(p):
    i = sqrt_minus_one(p)
    return {1: np.array([[0, 1], [1, 0]], dtype=np.int64),
            2: np.array([[0, (p - i) % p], [i, 0]], dtype=np.int64),
            3: np.array([[1, 0], [0, p - 1]], dtype=np.int64)}, i


def site_op(op, site, n, p):
    return bond_op(op, site, None, n, p)


def bond_op(op, a, c, n, p):
    """P_op at sites a (and c, if given), identity elsewhere, as ONE kron.

    Built as a single tensor product rather than as site_op(a) @ site_op(c).
    That matmul was safe only because a Pauli tensor product is monomial, so
    its unreduced accumulator ever held one product; with 16 terms at N=4 and
    entries near p, a non-monomial operand would overflow int64 silently. A
    correctness that depends on an unstated property of the operands is a trap
    rather than a proof, and the kron form does not need the property.
    """
    out = np.ones((1, 1), dtype=np.int64)
    eye = np.eye(2, dtype=np.int64)
    for k in range(n):
        out = np.kron(out, op if k in (a, c) else eye) % p
    return out


def build_L(n, edges, deph, field_letters, p, signs=None,
            bond_terms=(1, 2, 3), j_num=None, gamma_num=None):
    """(L mod p, shift mod p). j_num / gamma_num give the non-uniform variants."""
    PA, i_img = pauli(p)
    dim = 2 ** n
    H = np.zeros((dim, dim), dtype=np.int64)
    for b, (a, c) in enumerate(edges):
        jj = 1 if j_num is None else (j_num[b] * inv_p(FIELD_DEN, p)) % p
        for t in bond_terms:
            H = (H + jj * bond_op(PA[t], a, c, n, p)) % p
    for s in range(n):
        if field_letters[s] == 0:
            continue
        mag = FIELD_NUM[s] * (1 if signs is None else signs[s])
        h = (mag * inv_p(FIELD_DEN, p)) % p
        H = (H + h * site_op(PA[field_letters[s]], s, n, p)) % p
    ident = np.eye(dim, dtype=np.int64)
    minus_i = (p - i_img) % p
    L = (minus_i * ((np.kron(H, ident) - np.kron(ident, H.T)) % p)) % p
    shift = 0
    for s in range(n):
        if deph[s] == 0:
            continue
        g = ((GAMMA[0] if gamma_num is None else gamma_num[s])
             * inv_p(GAMMA[1] if gamma_num is None else 100, p)) % p
        A = site_op(PA[deph[s]], s, n, p)
        # reduce the kron BEFORE scaling: unreduced entries reach ~1e18 and the
        # scalar is ~1e9, so the product overflows int64 silently. This was
        # wrong until the object-dtype reference disagreed with it.
        L = (L + g * (np.kron(A, A.T) % p)) % p
        L = (L - g * np.eye(dim * dim, dtype=np.int64)) % p
        shift = (shift + g) % p
    return L, (2 * shift) % p


def palindromic(n, edges, deph, fld, signs=None, bond_terms=(1, 2, 3),
                j_num=None, gamma_num=None):
    """False is a PROOF of a break; True is a certificate, see the header."""
    for p in PRIMES:
        L, shift = build_L(n, edges, deph, fld, p, signs=signs,
                           bond_terms=bond_terms, j_num=j_num,
                           gamma_num=gamma_num)
        eye = np.eye(L.shape[0], dtype=np.int64)
        for _ in range(POINTS):
            x = int(rng.integers(0, p))
            y = (-x - shift) % p
            if (det_mod_np_p((x * eye - L) % p, p)
                    != det_mod_np_p((y * eye - L) % p, p)):
                return False
    return True


def components(n, edges):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i, j in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    groups = {}
    for s in range(n):
        groups.setdefault(find(s), set()).add(s)
    return [frozenset(v) for v in groups.values()]


def predicate(n, edges, deph, field, bond_terms=(1, 2, 3)):
    ceiling = 2 if len(bond_terms) >= 2 else 3
    for comp in components(n, edges):
        d = {deph[s] for s in comp if deph[s] != 0}
        if not d:
            continue
        if len(d) > ceiling:
            return False, 'clause1'
        f = {field[s] for s in comp if field[s] != 0}
        if len(f) > 1:
            return False, 'clause2-not-common'
        if f and next(iter(f)) in d:
            return False, 'clause2-not-orthogonal'
    return True, 'holds'


# At N=3 the four named topologies are only TWO graphs: chain [(0,1),(1,2)] and
# star [(0,1),(0,2)] are both the path P3, and ring and complete are the SAME
# edge set, the triangle K3. Running all four would re-verify the same physics
# up to four times and report it as four confirmations.
GRAPHS_N3 = (('P3 (chain=star)', [(0, 1), (1, 2)]),
             ('K3 (ring=complete)', [(0, 1), (1, 2), (2, 0)]),
             ('bond+isolate', [(0, 1)]))

BOND_SETS = (((1, 2, 3), 'XX+YY+ZZ'), ((1, 2), 'XX+YY'), ((1, 3), 'XX+ZZ'),
             ((2, 3), 'YY+ZZ'), ((1,), 'XX'), ((2,), 'YY'), ((3,), 'ZZ'))

N4_DEPHASING = ((3, 3, 3, 3), (1, 1, 1, 1), (2, 2, 2, 2), (3, 1, 3, 1),
                (3, 3, 1, 1), (3, 1, 1, 3), (3, 3, 2, 2), (1, 2, 1, 2),
                (2, 3, 3, 2), (3, 3, 3, 0), (3, 3, 0, 0), (3, 0, 0, 0),
                (3, 1, 0, 0), (0, 3, 1, 0), (3, 1, 2, 0), (3, 1, 2, 3),
                (1, 2, 3, 1), (0, 0, 0, 0))


def sign_patterns(fld):
    """Sign patterns that are DISTINCT configurations.

    A sign on a site with no field term changes nothing, so enumerating all
    2^N would re-run the same operator and report the repeats as coverage.
    Only the sites actually carrying a field vary.
    """
    free = [s for s in range(len(fld)) if fld[s] != 0]
    if not free:
        yield None
        return
    for combo in itertools.product((1, -1), repeat=len(free)):
        signs = [1] * len(fld)
        for s, v in zip(free, combo):
            signs[s] = v
        yield tuple(signs)


def score(n, edges, cases, bond_terms=(1, 2, 3), j_num=None, gamma_num=None):
    fp, fn, seen = [], [], 0
    for deph, fld, signs in cases:
        measured = palindromic(n, edges, deph, fld, signs=signs,
                               bond_terms=bond_terms, j_num=j_num,
                               gamma_num=gamma_num)
        predicted, reason = predicate(n, edges, deph, fld,
                                      bond_terms=bond_terms)
        seen += 1
        if predicted and not measured:
            fp.append((deph, fld, signs, reason))
        elif measured and not predicted:
            fn.append((deph, fld, signs, reason))
    return seen, fp, fn


def show(label, seen, fp, fn):
    print(f'  {label:34s} rows={seen:6d}  false+={len(fp):5d}  '
          f'false-={len(fn):5d}')
    for tag, items in (('FALSE POSITIVE', fp), ('FALSE NEGATIVE', fn)):
        for deph, fld, signs, reason in items[:4]:
            d = ''.join(LETTER_NAME[a] for a in deph)
            f = ''.join(LETTER_NAME[a] for a in fld)
            s = 'all+' if signs is None else ''.join(
                '+' if v > 0 else '-' for v in signs)
            print(f'      {tag}: deph={d} field={f} signs={s} ({reason})')
        if len(items) > 4:
            print(f'      ... and {len(items) - 4} more of the same kind')


def grid3():
    return [(d, f, s)
            for d in itertools.product(LETTERS, repeat=3)
            for f in itertools.product(LETTERS, repeat=3)
            for s in sign_patterns(f)]


def main():
    print('F138 clause 2, decided over several primes')
    print('=' * 74)
    print(__doc__.split('\n\n', 1)[1].strip())
    print()
    print('## Stage A: N=3, exhaustive, full Heisenberg bond')
    print()
    print('  Every 4^3 dephasing assignment x every 4^3 field assignment x')
    print('  every DISTINCT sign pattern, on each distinct graph. Signs are')
    print('  varied only where a field term exists, so these are distinct')
    print('  configurations rather than re-runs.')
    print()
    for name, edges in GRAPHS_N3:
        show(name, *score(3, edges, grid3()))
    print()
    print('## Stage B: N=4, stated cross-section, full Heisenberg bond')
    print()
    print('  18 named dephasing patterns x every 4^4 field assignment. NOT')
    print('  exhaustive and not claimed to be; signs are exhaustive at N=3.')
    print()
    for name, edges in (('P4 chain', [(0, 1), (1, 2), (2, 3)]),
                        ('C4 ring', [(0, 1), (1, 2), (2, 3), (3, 0)]),
                        ('two bonds', [(0, 1), (2, 3)])):
        cases = [(d, f, None) for d in N4_DEPHASING
                 for f in itertools.product(LETTERS, repeat=4)]
        show(name, *score(4, edges, cases))
    print()
    print('## Stage E: the coupling axis F138 quantifies over')
    print()
    print('  The N=3 grid again with non-uniform rational J per bond and')
    print('  non-uniform rational gamma per site. Uniform couplings ADD')
    print('  symmetry, which biases toward the palindrome, i.e. toward hiding')
    print('  exactly the false positives the headline rests on.')
    print()
    jn, gn = [43, 55, 13], [3, 7, 11]
    for name, edges in GRAPHS_N3:
        cases = [(d, f, None)
                 for d in itertools.product(LETTERS, repeat=3)
                 for f in itertools.product(LETTERS, repeat=3)]
        show(f'{name} non-uniform J,gamma',
             *score(3, edges, cases, j_num=jn[:len(edges)], gamma_num=gn))
    print()
    print('## Stage C: tightness against bond content, INSIDE the domain')
    print()
    print('  F138 reads "For H a sum of bond terms ... holds exactly when",')
    print('  and its two-term proviso conditions only the CEILING in clause 1,')
    print('  not the reach of the law. So a two-term bond such as XX+YY is')
    print('  squarely inside the stated domain, and a row there where the')
    print('  predicate forbids and the spectrum pairs is a counterexample to')
    print('  the "only when" direction, not an out-of-scope curiosity.')
    print()
    for bond_terms, bname in BOND_SETS:
        for name, edges in GRAPHS_N3:
            cases = [(d, f, None)
                     for d in itertools.product(LETTERS, repeat=3)
                     for f in itertools.product(LETTERS, repeat=3)]
            show(f'{bname} / {name}',
                 *score(3, edges, cases, bond_terms=bond_terms))
        print()


if __name__ == '__main__':
    main()
