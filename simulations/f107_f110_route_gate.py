"""Gate: the routes of F107, F109 and F110 at every body count, and what F108's mirrors reach.

F107 (truly => y_par = 0), F109 (y_par-homogeneous mother-cell soft => y_par = 1) and F110 Aspect A (hard
Klein-homogeneous pairs only in the diagonal Klein cell) are statements about Pauli strings of any body count. Their routes, checked here:

  F107  The proof derives the Z criterion from Pi_Z = R.D and carries it to X and Y by Pi_Y = Pi_Z^-1 and the
        Hadamard transport of Pi_Z. Per site this is one fact, checked in G1: each canonical palindromizer Pi_D
        (framework.symmetry.pi_action, imported, not copied) carries left multiplication by a letter a to
        +-(right multiplication by a) and back, with the minus signs on the two letters other than the mirror's flip
        letter A(D), one on each side:
            Pi_Z: left Z, right Y      Pi_X: left X, right Y      Pi_Y: left Y, right Z
        For a string both maps are products over sites, so Pi_D L_P Pi_D^-1 = -i(eps_l(P) r_P - eps_r(P) l_P) and
        M_P = 0 iff both products are +1: both letters other than A(D) occur an even number of times. Every
        criterion contains "#Y even" (A(D) is X or Z). G2 compares this with the directly computed M_P.
  F109, F110  Pi_D^2 is the turn by A(D)^N (checked per site in G1), so a string is Pi^2-D-even exactly when it
        commutes with A(D)^N; A(D)^N
        anticommutes with every jump D_l. The one-sided multiplication rho -> rho A(D)^N therefore keeps L_H and
        reflects the dissipator, R L R^-1 = -L^dag - 2 sigma, for EVERY Pi^2-D-even H, and by F158's sufficiency step
        (PROOF_PALINDROME_TWO_END_COUNT section (e)) the spectrum pairs about -sigma. The third letter's cell is
        coloured by its own letter; only the diagonal cell can be hard (G4).
  F108  Its mirrors (Parts 1-3, defined below as their proofs define them) are canonical mirrors composed with
        conjugation by Y^N, applied first: Pi5(Z) = Pi_Z.Ad_Y, Pi5(Y) = Pi_Y.Ad_Y, Pi5(X) = Pi_X^-1.Ad_Y (G1). They
        carry the same per-site form with other signs (G1), and cover exactly the Pi^2-D-even strings of EVEN weight
        (an even number of non-identity letters), on any number of sites. The mother cell's non-truly strings have
        all three counts odd, hence odd weight, so no F108 mirror covers any of them (G3). The operator-space Klein
        V4 (D the transpose, H the X <-> Z letter swap fixing I and Y, Q_zx = H.D) commutes with Ad_Y and permutes
        Pi_Z^{+-1}, Pi_X^{+-1} simply transitively; on F108's three (the X one being Pi_X^-1.Ad_Y) it swaps Z<->Y,
        Z<->X and X<->Y and sends the third variant to Pi_X.Ad_Y (G1).

Every comparison is exact: Pauli label maps carry phases in {+-1, +-i} and superoperators are small-integer complex
matrices, compared with ==. Expected values are literals with their derivation in a comment.
Run: python simulations/f107_f110_route_gate.py   (about ten seconds)"""
import itertools
import os
import random
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from framework.symmetry import pi_action                  # the canonical mirrors, imported, not copied

LET = 'IXYZ'
PM = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], complex),
      'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1]).astype(complex)}
KLEIN = {'I': (0, 0), 'X': (1, 0), 'Z': (0, 1), 'Y': (1, 1)}
NAME = {v: k for k, v in KLEIN.items()}
MUL = {}
for a, b in itertools.product(LET, repeat=2):
    m = PM[a] @ PM[b]
    for c in LET:
        t = np.trace(PM[c].conj().T @ m) / 2
        if t != 0:
            MUL[(a, b)] = (complex(t), c)
A = {'Z': 'X', 'X': 'Z', 'Y': 'X'}                         # the canonical mirror's flip letter (F103 section 8)
F108 = {                                                  # PROOF_F108_PART1/2/3, "The ... operator" sections
    'Z': {'I': ('X', 1), 'X': ('I', -1), 'Y': ('Z', 1j), 'Z': ('Y', -1j)},
    'X': {'I': ('Z', 1), 'Z': ('I', -1), 'X': ('Y', -1j), 'Y': ('X', 1j)},
    'Y': {'I': ('X', 1), 'X': ('I', -1), 'Y': ('Z', -1j), 'Z': ('Y', 1j)},
}
fails = 0


def check(name, ok):
    global fails
    fails += (not ok)
    print('  %-112s %s' % (name, 'PASS' if ok else 'FAIL'))


def canonical(D):
    def f(c):
        (na, nb), ph = pi_action(c, dephase_letter=D)
        return NAME[(na, nb)], ph
    return f


def site_matrix(f):                                       # a per-site letter map as a 4x4 matrix on labels
    M = np.zeros((4, 4), complex)
    for j, c in enumerate(LET):
        c2, ph = f(c)
        M[LET.index(c2), j] = ph
    return M


def lmul(a):
    return site_matrix(lambda s: (MUL[(a, s)][1], MUL[(a, s)][0]))


def rmul(a):
    return site_matrix(lambda s: (MUL[(s, a)][1], MUL[(s, a)][0]))


def inverse(M):                                           # a signed permutation: its inverse is its adjoint
    Mi = M.conj().T
    assert np.array_equal(Mi @ M, np.eye(M.shape[0]))
    return Mi


def sign_table(f):
    """{letter: (eps_l, eps_r)} with M l_a M^-1 = eps_l r_a and M r_a M^-1 = eps_r l_a; None where not of that form."""
    M = site_matrix(f)
    Mi = inverse(M)
    out = {}
    for a in LET:
        cl, cr = M @ lmul(a) @ Mi, M @ rmul(a) @ Mi
        out[a] = (next((e for e in (1, -1) if np.array_equal(cl, e * rmul(a))), None),
                  next((e for e in (1, -1) if np.array_equal(cr, e * lmul(a))), None))
    return out


def minus_letters(T):
    return ''.join(a for a in LET if T[a][0] == -1), ''.join(a for a in LET if T[a][1] == -1)


print('G1  F107: the per-site sign table of each canonical mirror (literal expectations)')
EXPECT = {'Z': ('Z', 'Y'), 'X': ('X', 'Y'), 'Y': ('Y', 'Z')}  # (left-side minus letter, right-side minus letter)
for D in 'ZXY':
    T = sign_table(canonical(D))
    ok_form = all(v[0] is not None and v[1] is not None for v in T.values())
    check('Pi_%s: every letter is carried to +-the other side; minus signs (left, right) = %s, the two letters other '
          'than A = %s' % (D, minus_letters(T), A[D]), ok_form and minus_letters(T) == EXPECT[D])
check('control through the same door: F108 Part 1\'s mirror has a different table (%s, %s), so the table reads the '
      'phases' % minus_letters(sign_table(lambda c: F108['Z'][c])), minus_letters(sign_table(lambda c: F108['Z'][c])) != EXPECT['Z'])
# F108's three mirrors have the same per-site form with other signs (literal expectations); with the product argument
# this fixes their reach at every N. String signs (left, right): Z variant ((-1)^#X, (-1)^weight), X variant
# ((-1)^weight, (-1)^#Z), Y variant ((-1)^weight, (-1)^#X); both +1 is a Pi^2-D-even string of even weight
EXPECT5 = {'Z': ('X', 'XYZ'), 'X': ('XYZ', 'Z'), 'Y': ('XYZ', 'X')}
for D in 'ZXY':
    T5 = sign_table(lambda c, D=D: F108[D][c])
    check('F108 %s variant: every letter carried to +-the other side; minus signs (left, right) = %s'
          % (D, minus_letters(T5)), all(None not in v for v in T5.values()) and minus_letters(T5) == EXPECT5[D])
# the two transports the proof uses, per site and exact: Pi_Y = Pi_Z^-1, and Pi_X = Ad_Had . Pi_Z . Ad_Had
had = site_matrix(lambda c: {'I': ('I', 1), 'X': ('Z', 1), 'Y': ('Y', -1), 'Z': ('X', 1)}[c])
PZ, PX, PY = (site_matrix(canonical(D)) for D in 'ZXY')
check('transports: Pi_Y = Pi_Z^-1 and Pi_X = Ad_Had.Pi_Z.Ad_Had per site (Ad_Had: X <-> Z, Y -> -Y), exact',
      np.array_equal(PY @ PZ, np.eye(4)) and np.array_equal(had @ PZ @ had, PX) and not np.array_equal(PX, PZ))


def ad(letter):                                           # per-site conjugation by a letter: the sign of commuting
    return site_matrix(lambda c: (c, 1 if c in ('I', letter) else -1))


# the square of each canonical mirror is the turn by its flip letter, Pi_D^2 = Ad_A(D) per site
check('Pi_D^2 = Ad_A(D) per site for D = Z, X, Y (the turn by the flip letter), exact',
      all(np.array_equal(P @ P, ad(A[D])) for D, P in (('Z', PZ), ('X', PX), ('Y', PY))))
# F108's mirrors are canonical mirrors composed with conjugation by Y^N (applied first), exact per site
P5 = {D: site_matrix(lambda c, D=D: F108[D][c]) for D in 'ZXY'}
check('F108 mirrors: Pi5(Z) = Pi_Z.Ad_Y, Pi5(Y) = Pi_Y.Ad_Y, Pi5(X) = Pi_X^-1.Ad_Y per site, exact',
      np.array_equal(P5['Z'], PZ @ ad('Y')) and np.array_equal(P5['Y'], PY @ ad('Y'))
      and np.array_equal(P5['X'], inverse(PX) @ ad('Y')) and not np.array_equal(P5['X'], PX @ ad('Y')))
# the operator-space Klein V4 on the mirrors: D the transpose (Y -> -Y), H the X <-> Z letter swap fixing I and Y,
# Q_zx = H.D = Ad_Had. Each commutes with Ad_Y, so it acts on F108's mirrors as on the canonical ones; on the four
# mirrors Pi_Z, Pi_Y = Pi_Z^-1, Pi_X, Pi_X^-1 it swaps two pairs (literal tables), and F108's set holds the X mirror
# of the other orientation, Pi_X^-1.Ad_Y, its fourth being Pi_X.Ad_Y
dsw = site_matrix(lambda c: (c, -1 if c == 'Y' else 1))
hsw = site_matrix(lambda c: {'I': ('I', 1), 'X': ('Z', 1), 'Y': ('Y', 1), 'Z': ('X', 1)}[c])
V4 = {'D': dsw, 'H': hsw, 'Qzx': had}


def acts(G, fam):                                         # where G.P.G^-1 lands in fam, exactly; '' when outside
    return {k: next((n for n, W in fam.items() if np.array_equal(G @ V @ inverse(G), W)), '') for k, V in fam.items()}


CANON4 = {'Z': PZ, 'Y': PY, 'X': PX, 'X-1': inverse(PX)}
F108_4 = {'Z': P5['Z'], 'Y': P5['Y'], 'X': P5['X'], 'X+': PX @ ad('Y')}
EXPECT_C = {'D': {'Z': 'Y', 'Y': 'Z', 'X': 'X-1', 'X-1': 'X'}, 'H': {'Z': 'X-1', 'Y': 'X', 'X': 'Y', 'X-1': 'Z'},
            'Qzx': {'Z': 'X', 'Y': 'X-1', 'X': 'Z', 'X-1': 'Y'}}
EXPECT_F = {'D': {'Z': 'Y', 'Y': 'Z', 'X': 'X+', 'X+': 'X'}, 'H': {'Z': 'X', 'X': 'Z', 'Y': 'X+', 'X+': 'Y'},
            'Qzx': {'Z': 'X+', 'X+': 'Z', 'X': 'Y', 'Y': 'X'}}
check('Klein V4: H.D = Q_zx, each element commutes with Ad_Y; on Pi_Z, Pi_Y, Pi_X, Pi_X^-1: D swaps Z<->Y, Q_zx Z<->X, '
      'H Y<->X, each sending the third to Pi_X^-1, exact',
      np.array_equal(hsw @ dsw, had) and all(np.array_equal(G @ ad('Y'), ad('Y') @ G) for G in V4.values())
      and all(acts(G, CANON4) == EXPECT_C[g] for g, G in V4.items()))
check('Klein V4 on F108\'s three: D swaps Z<->Y, H Z<->X, Q_zx X<->Y, each sending the third to Pi_X.Ad_Y; '
      'Q_zx.Pi5(Z).Q_zx and H.Pi5(Y).H are not +-Pi5(X), exact',
      all(acts(G, F108_4) == EXPECT_F[g] for g, G in V4.items())
      and not any(np.array_equal(had @ P5['Z'] @ had, s * P5['X']) or np.array_equal(hsw @ P5['Y'] @ hsw, s * P5['X'])
                  for s in (1, -1)))

N = 4
labels = [''.join(t) for t in itertools.product(LET, repeat=N)]
idx = {s: i for i, s in enumerate(labels)}
strings = [s for s in labels if set(s) != {'I'}]


def full_map(f):
    M = np.zeros((4 ** N, 4 ** N), complex)
    for s in labels:
        ph, out = 1, ''
        for c in s:
            c2, p = f(c)
            ph *= p
            out += c2
        M[idx[out], idx[s]] = ph
    return M


def LH(p):                                                # -i [p, .] on the label basis, built from the product table
    M = np.zeros((4 ** N, 4 ** N), complex)
    for s in labels:
        ph, out, anti = 1, '', 0
        for x, y in zip(p, s):
            t, c = MUL[(x, y)]
            ph *= t
            out += c
            anti ^= (x != 'I' and y != 'I' and x != y)
        if anti:
            M[idx[out], idx[s]] = -2j * ph
    return M


def LD(D):                                                # sum_l D[D_l] at unit rates on the label basis
    diag = [-2 * sum(c not in ('I', D) for c in s) for s in labels]
    return np.diag(diag).astype(complex)


def weight(s):
    return len(s) - s.count('I')


def klein(s):
    a = b = 0
    for c in s:
        a ^= KLEIN[c][0]
        b ^= KLEIN[c][1]
    return (a, b)


def pi2_even(s, D):                                       # Pi_D^2 parity: bit_b for Z and Y, bit_a for X
    return klein(s)[0 if D == 'X' else 1] == 0


def commutes(s, t):
    return sum(x != 'I' and y != 'I' and x != y for x, y in zip(s, t)) % 2 == 0


LHs = {p: LH(p) for p in strings}


def truly_set(f):
    Q = full_map(f)
    Qi = inverse(Q)
    return {p for p in strings if not np.any(Q @ LHs[p] @ Qi + LHs[p])}


def product_rule(T, p):
    if any(None in T[c] for c in p):                      # a letter not carried to +-the other side: no product form
        return False
    el = er = 1
    for c in p:
        el *= T[c][0]
        er *= T[c][1]
    return el == 1 and er == 1


print('G2  F107 at N = 4, every non-identity string (body counts 1-4): the direct M against the sign product')
for D in 'ZXY':
    Q = full_map(canonical(D))
    check('Pi_%s palindromizes the dissipator exactly: Pi L_D Pi^-1 + L_D + 2 sigma = 0' % D,
          not np.any(Q @ LD(D) @ inverse(Q) + LD(D) + 2 * N * np.eye(4 ** N)))
    tru = truly_set(canonical(D))
    T = sign_table(canonical(D))
    rule = {p for p in strings if all(p.count(c) % 2 == 0 for c in 'XYZ' if c != A[D])}
    # 71 = (4^4 + 2^4 + 2^4 + 0^4)/4 - 1: words with two given letters each even, minus the identity
    check('Pi_%s: M_P = 0 on exactly the strings the sign product admits = both letters other than %s even, %d strings'
          % (D, A[D], len(tru)), tru == {p for p in strings if product_rule(T, p)} == rule and len(tru) == 71)
    check('Pi_%s: every truly string has #Y even (y_par = 0)' % D, bool(tru) and all(p.count('Y') % 2 == 0 for p in tru))

print('G3  F108: its three mirrors cover exactly the Pi^2-D-even strings of even weight')
for D in 'ZXY':
    f = (lambda c, D=D: F108[D][c])
    Q = full_map(f)
    check('F108 %s mirror palindromizes the %s dissipator exactly' % (D, D),
          not np.any(Q @ LD(D) @ inverse(Q) + LD(D) + 2 * N * np.eye(4 ** N)))
    tru = truly_set(f)
    even_w = {p for p in strings if pi2_even(p, D) and weight(p) % 2 == 0}
    missed = [p for p in strings if pi2_even(p, D) and p not in tru]
    by_w = {w: sum(weight(p) == w for p in missed) for w in (1, 2, 3, 4)}
    # missed: the 4 single-site fields A(D)_l, and the three-body strings, (4 choices of three sites) x (the 13
    # Pi^2-even words of three letters from XYZ: A(D)A(D)A(D) and one A(D) with two of the other two letters, 1 + 3 x 4)
    n3 = sum(1 for w in itertools.product('XYZ', repeat=3) if pi2_even(''.join(w), D))
    check('F108 %s: covered = Pi^2-even AND even weight (%d strings, its own bilinears included); missed by weight %s'
          % (D, len(tru), by_w),
          tru == even_w and by_w == {1: 4, 2: 0, 3: 4 * n3, 4: 0} and n3 == 13
          and all(p in tru for p in strings if weight(p) == 2 and pi2_even(p, D)))
    mother_nontruly = [p for p in strings if klein(p) == (0, 0)
                       and not all(p.count(c) % 2 == 0 for c in 'XYZ' if c != A[D])]
    # all three counts odd with at most four letters: the 24 arrangements of X, Y, Z and one I
    check('F108 %s: the mother cell\'s %d non-truly strings all have odd weight and none is covered'
          % (D, len(mother_nontruly)),
          len(mother_nontruly) == 24 and all(weight(p) % 2 == 1 and p not in tru for p in mother_nontruly))

print('G4  F109, F110: Pi^2-D-even = commutes with A(D)^N; the colouring by A(D)^N reflects L for every such H')
check('for every string at N = 4 and every D: Pi^2-D-even <=> commutes with A(D)^N',
      all(pi2_even(p, D) == commutes(p, A[D] * N) for p in strings for D in 'ZXY'))
d = 2 ** N


def op(s):
    m = np.array([[1.0 + 0j]])
    for c in s:
        m = np.kron(m, PM[c])
    return m


def lr(X, Y):                                             # rho -> X rho Y, row-stacking vec
    return np.kron(X, Y.T)


Id = np.eye(d)


def lindbladian(H, D, g):
    out = -1j * (lr(H, Id) - lr(Id, H))
    for l in range(N):
        J = op(''.join(D if m == l else 'I' for m in range(N)))
        out = out + g[l] * (lr(J, J) - lr(Id, Id))
    return out


def colour_residual(Lv, letter, g):
    K = op(letter * N)
    return Lv.conj().T + lr(Id, K) @ Lv @ lr(Id, K) + 2 * sum(g) * np.eye(d * d)


random.seed(20261004)
for D in 'ZXY':
    third = NAME[(KLEIN[D][0] ^ KLEIN[A[D]][0], KLEIN[D][1] ^ KLEIN[A[D]][1])]
    even = [p for p in strings if pi2_even(p, D)]
    ok, kinds = True, set()
    for _ in range(60):
        chosen = random.sample(even, random.randint(2, 6))
        kinds.add((len({klein(p) for p in chosen}) > 1, len({p.count('Y') % 2 for p in chosen}) > 1,
                   any(weight(p) % 2 for p in chosen)))
        H = sum(random.choice([1, 2, -1, 3, -2]) * op(p) for p in chosen)
        g = [random.randint(0, 4) for _ in range(N)]
        ok &= not np.any(colour_residual(lindbladian(H, D, g), A[D], g))
    check('D=%s: 60 random Pi^2-even H, 2-6 strings of any body count, rates 0-4, exact; mixed cells, mixed y_par and '
          'odd weights all occur' % D, ok and (True, True, True) in kinds)
    # what neither mirror reaches: a pair of the A(D) cell across both y_par halves, and a mother string
    cell_a = [p for p in strings if klein(p) == KLEIN[A[D]]]
    y0 = next((p for p in cell_a if p.count('Y') % 2 == 0 and weight(p) == 3), None)
    y1 = next((p for p in cell_a if p.count('Y') % 2 == 1 and weight(p) == 2), None)
    mom = next((p for p in strings if klein(p) == (0, 0) and weight(p) == 3), None)
    if None in (y0, y1, mom):                             # a wrong flip-letter table reports, it does not crash
        check('D=%s: the %s cell holds a weight-3 string with #Y even and a weight-2 one with #Y odd, the mother '
              'cell a weight-3 string' % (D, A[D]), False)
    else:
        Qc, Q5 = full_map(canonical(D)), full_map(lambda c, D=D: F108[D][c])
        for name, plist in (('%s + %s (the %s cell, both y_par)' % (y0, y1, A[D]), (y0, y1)),
                            ('%s (mother)' % mom, (mom,))):
            Hs = sum(LHs[p] for p in plist)
            neither = bool(np.any(Qc @ Hs @ inverse(Qc) + Hs)) and bool(np.any(Q5 @ Hs @ inverse(Q5) + Hs))
            g = [1, 2, 3, 4]
            col = not np.any(colour_residual(lindbladian(sum(op(p) for p in plist), D, g), A[D], g))
            check('D=%s, H = %s: neither Pi_%s nor F108\'s mirror has M = 0, the colouring is exact' % (D, name, D),
                  neither and col)
    ok3 = True
    cell3 = [p for p in strings if klein(p) == KLEIN[third]]
    for _ in range(40):
        H = sum(random.choice([1, 2, -1, 3]) * op(p) for p in random.sample(cell3, random.randint(1, 5)))
        g = [random.randint(0, 4) for _ in range(N)]
        ok3 &= not np.any(colour_residual(lindbladian(H, D, g), third, g))
    check('D=%s: the third letter\'s cell (%s) is coloured by %s^N, 40 random H exact' % (D, third, third), ok3)
    okm = True                                            # the Mother cell takes either lit letter, the third one too
    cellm = [p for p in strings if klein(p) == (0, 0)]
    for _ in range(20):
        H = sum(random.choice([1, 2, -1, 3]) * op(p) for p in random.sample(cellm, random.randint(1, 4)))
        g = [random.randint(0, 4) for _ in range(N)]
        okm &= not np.any(colour_residual(lindbladian(H, D, g), third, g))
    check('D=%s: the Mother cell is coloured by the third letter %s^N as well, 20 random H exact' % (D, third), okm)
    # controls through the same door: a Pi^2-odd string breaks the A(D) colouring; the diagonal cell takes no lit letter
    ctrl = True
    for _ in range(10):
        g = [random.randint(1, 4) for _ in range(N)]
        Hodd = op(random.choice([p for p in strings if not pi2_even(p, D)])) + op(random.choice(even))
        ctrl &= bool(np.any(colour_residual(lindbladian(Hodd, D, g), A[D], g)))
        Hd = sum(op(p) for p in random.sample([p for p in strings if klein(p) == KLEIN[D]], 3))
        ctrl &= all(bool(np.any(colour_residual(lindbladian(Hd, D, g), lit, g))) for lit in (A[D], third))
    check('D=%s controls: one Pi^2-odd string breaks it; no lit letter colours the diagonal cell' % D, ctrl)

print('\nALL PASS' if fails == 0 else '\n%d FAIL' % fails)
sys.exit(0 if fails == 0 else 1)
