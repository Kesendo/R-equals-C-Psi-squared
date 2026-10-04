"""Gate: F103's open question 3 (Pattern B vs Pattern C in the off-diagonal soft cells), read on the letter cube.

Rule. The canonical palindromizer for dephasing letter D (framework.symmetry.pi_action) flips every site by one letter
A(D): Z-dephasing -> X, X-dephasing -> Z, Y-dephasing -> X. Its truly strings (M = Pi L_H Pi^-1 + L_H = 0) are those in
which both letters other than A occur an even number of times (F85's syntactic "#Y and #Z even" under Z, rotated), and
a pair is truly iff both strings are (M sends each string to the single label P xor sigma, so different strings never
cancel). With K the pair's Klein letter (the xor of its letters):
  K = I mother cell, K = D diagonal cell, K = A Pattern C (y_par = 0 pairs truly, y_par = 1 soft),
  K = D xor A, the third letter, Pattern B (nothing truly, the whole cell soft, so its shape is the cell's own
  enumeration: 10/6 strings at k = 3, 32/32 at k = 4).
No off-diagonal cell is ever hard: for a cell letter L outside {I, D}, W(rho) = D^N rho L^N D^N palindromizes every
Hamiltonian of that cell (strings of net letter L commute with L^N and anticommute with D^N; the right factor L^N
reflects the dissipator, the turn by D^N flips H). The assignment follows the choice of mirror: the quarter turn about
D carries A to D xor A and gives another palindromizer of the same dissipator, under which B and C swap; and the
phases count too (F108 Part 1's Pi_5bilinear flips by X like the canonical Pi_Z, and against it the X cell's truly half
is y_par = 1).

Routes, all exact (Gaussian-integer superoperators in the Pauli basis, comparisons with ==):
  G1  N = 4, k = 3, the 294 pairs: M computed directly for the canonical mirror and for the turned mirror, both
      checked to palindromize the dissipator exactly; truly pair by pair = the syntactic rule; the canonical counts
      reproduce F103 section 4's truly table.
  G2  k = 4: the rule's truly / soft counts in every mother and off-diagonal cell, and no truly in the diagonal ones,
      reproduce F106's measured grid (simulations/results/f87_z2cubed_split_n4_k4_counts.json, from the C# classifier).
  G3  the colouring rho -> rho K^N (K a letter outside {I, D}, lit on every site) gives -L^dag - 2 sigma exactly for
      the mother and both off-diagonal cells, and W = D^N rho K^N D^N gives -L - 2 sigma in the off-diagonal cells:
      N = 4, random Hamiltonians of 1 to 5 strings of any body count, random per-site rates, all three letters;
      controls: no lit letter colours the diagonal cell, W fails on the mother cell.
  G4  Pi_5bilinear: palindromizes the Z dissipator exactly, its truly X-cell pairs are the y_par = 1 half, and XZZ
      and the mother cell's XYZ permutations are not truly against it (F108 Part 1's mirror reaches the Pi^2-even
      strings of even weight only; simulations/f107_f110_route_gate.py G3 counts them).
  G5  the 256 phase choices on the X flip all palindromize the Z dissipator; their X-cell truly sets run from
      empty to all sixteen strings, so the phases decide which strings are truly, not merely which half.
Run: python simulations/f103_pattern_bc_gate.py   (a few seconds)"""
import itertools
import json
import os
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

N, d = 4, 16
P = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], complex),
     'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1]).astype(complex)}
MUL = {}
for a, b in itertools.product('IXYZ', repeat=2):
    m = P[a] @ P[b]
    for c in 'IXYZ':
        t = np.trace(P[c].conj().T @ m) / 2
        if t != 0: MUL[(a, b)] = (t, c)
def string(s):
    m = np.array([[1.0 + 0j]])
    for c in s: m = np.kron(m, P[c])
    return m
labels = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
idx = {s: i for i, s in enumerate(labels)}
V = np.array([string(s).reshape(-1) for s in labels]).T
def to_pauli(S): return V.conj().T @ S @ V / d
def lr(A, B): return np.kron(A, B.T)
I = np.eye(d)
KLEIN = {'I': (0, 0), 'X': (1, 0), 'Z': (0, 1), 'Y': (1, 1)}
NAME = {v: k for k, v in KLEIN.items()}
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from framework.symmetry import pi_action                  # the canonical Pi, imported, not copied


def pi_letter(c, D):
    (a, b), phase = pi_action(c, dephase_letter=D)
    return NAME[(a, b)], phase
def letter_map_matrix(f):                  # f(letter) -> (letter, phase), applied sitewise
    M = np.zeros((4 ** N, 4 ** N), complex)
    for s in labels:
        ph, out = 1, ''
        for c in s:
            c2, p = f(c); ph *= p; out += c2
        M[idx[out], idx[s]] = ph
    return M
def quarter_turn(D):                       # Q -> i Q D for Q not in {I, D}
    def f(c):
        if c in ('I', D): return c, 1
        t, c2 = MUL[(c, D)]
        return c2, 1j * t
    return f
fails = 0
def check(name, ok):
    global fails
    fails += (not ok)
    print('  %-100s %s' % (name, 'PASS' if ok else 'FAIL'))
def place(t, l): return ''.join(t[i - l] if l <= i < l + 3 else 'I' for i in range(N))
window = [t for t in itertools.product('IXYZ', repeat=3) if set(t) != {'I'}]
LH = {}
for t in window:
    H = sum(string(place(t, l)) for l in range(N - 2))
    LH[t] = to_pauli(-1j * (lr(H, I) - lr(I, H)))
def klein(t):
    a = b = 0
    for c in t: a ^= KLEIN[c][0]; b ^= KLEIN[c][1]
    return (a, b)
pairs, seen = [], set()
for t1, t2 in itertools.product(window, repeat=2):
    if klein(t1) != klein(t2) or t1.count('Y') % 2 != t2.count('Y') % 2: continue
    key = tuple(sorted([t1, t2]))
    if key in seen: continue
    seen.add(key); pairs.append((t1, t2, klein(t1), t1.count('Y') % 2))
print('G1  N = 4, k = 3, both mirrors, direct M')
print('  pairs: %d' % len(pairs))
F103 = {'Z': {(0,0): (45,0), (0,1): (0,0), (1,0): (55,0), (1,1): (0,0)},
        'X': {(0,0): (45,0), (0,1): (55,0), (1,0): (0,0), (1,1): (0,0)},
        'Y': {(0,0): (45,0), (0,1): (0,0), (1,0): (55,0), (1,1): (0,0)}}
A = {'Z': 'X', 'X': 'Z', 'Y': 'X'}
def syntactic(t, Aletter): return all(t.count(B) % 2 == 0 for B in 'XYZ' if B != Aletter)
for D in 'ZXY':
    Pi = letter_map_matrix(lambda c: pi_letter(c, D))
    U = letter_map_matrix(quarter_turn(D))
    Pi2 = U @ Pi @ U.conj().T
    LD = sum(to_pauli(lr(string(''.join(D if m == l else 'I' for m in range(N))), string(''.join(D if m == l else 'I' for m in range(N))).conj().T) - lr(I, I)) for l in range(N))
    Aprime = [c for c in 'XYZ' if c not in (D, A[D])][0]
    # U carries A(D) to A'(D): the mirror letter of Pi' is the image of A under the quarter turn
    img = quarter_turn(D)(A[D])[0]
    check('D=%s: quarter turn about D carries A=%s to %s = D xor A' % (D, A[D], img), img == Aprime)
    for name, Pm in (('canonical', Pi), ('other', Pi2)):
        res = Pm @ LD @ Pm.conj().T + LD + 2 * N * np.eye(4 ** N)
        check('D=%s %s mirror palindromizes the dissipator exactly (Pi L_D Pi^-1 + L_D + 2 sigma = 0)' % (D, name),
              np.array_equal(res, np.zeros_like(res)))
        M = {t: Pm @ LH[t] @ Pm.conj().T + LH[t] for t in window}
        counts = {K: [0, 0] for K in [(0,0), (0,1), (1,0), (1,1)]}
        agree = True
        for t1, t2, K, y in pairs:
            Mt = M[t1] + M[t2]
            tr = not np.any(Mt)
            counts[K][y] += tr
            pred = syntactic(t1, A[D] if name == 'canonical' else Aprime) and syntactic(t2, A[D] if name == 'canonical' else Aprime)
            agree &= (tr == pred)
        if name == 'canonical':
            check('D=%s canonical: truly counts %s reproduce F103 section 4' % (D, {k: tuple(v) for k, v in counts.items()}),
                  all(tuple(counts[K]) == F103[D][K] for K in counts))
        check('D=%s %s: direct M = 0 agrees pair by pair with the syntactic rule (letters other than %s even); counts %s'
              % (D, name, A[D] if name == 'canonical' else Aprime, {k: tuple(v) for k, v in counts.items()}), agree)

print('G2  k = 4, F106 measured grid')
data = json.load(open(os.path.join(HERE, 'results', 'f87_z2cubed_split_n4_k4_counts.json'), encoding='utf-8'))
meas = {(e['dephase'], (e['klein_a'], e['klein_b']), e['y_par'], e['trichotomy']): e['count'] for e in data['grid']}
strs4 = [s for s in itertools.product('IXYZ', repeat=4) if set(s) != {'I'}]
for D in 'ZXY':
    ok = True
    for K in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        if NAME[K] == D:                                  # diagonal: the rule gives no truly, the grid has none
            ok &= all(meas.get((D, K, y, 'truly'), 0) == 0 for y in (0, 1))
            continue
        for y in (0, 1):
            cell = [s for s in strs4 if klein(s) == K and s.count('Y') % 2 == y]
            nt = sum(syntactic(s, A[D]) for s in cell)
            tp, allp = nt * (nt + 1) // 2, len(cell) * (len(cell) + 1) // 2
            ok &= (tp, allp - tp, 0) == (meas.get((D, K, y, 'truly'), 0), meas.get((D, K, y, 'soft'), 0),
                                          meas.get((D, K, y, 'hard'), 0))
    check('D=%s: every mother and off-diagonal cell of F106 (k=4) = the rule (truly, soft, hard = 0)' % D, ok)
print('G3  every non-diagonal cell is coloured, at any body count; W is the colouring as an operator identity')
import random


def dissipator(J):
    JdJ = J.conj().T @ J
    return lr(J, J.conj().T) - 0.5 * (lr(JdJ, I) + lr(I, JdJ))


random.seed(20261004)
allstr = [''.join(s) for s in itertools.product('IXYZ', repeat=N) if set(s) != {'I'}]


def lindbladian(H, D, g):
    return -1j * (lr(H, I) - lr(I, H)) + sum(
        g[l] * dissipator(string(''.join(D if m == l else 'I' for m in range(N)))) for l in range(N))


def random_cell_H(cell):
    return sum(random.choice([1, 2, -1, 3]) * string(s) for s in random.sample(cell, random.randint(1, 5)))


for D in 'ZXY':
    third = NAME[(KLEIN[D][0] ^ KLEIN[A[D]][0], KLEIN[D][1] ^ KLEIN[A[D]][1])]
    cells = {Kl: [s for s in allstr if klein(s) == KLEIN[Kl]] for Kl in ('I', D, A[D], third)}
    for Kl, lit in (('I', A[D]), (A[D], A[D]), (third, third)):
        ok_col = ok_w = True
        for _ in range(40):
            H = random_cell_H(cells[Kl])
            g = [random.randint(1, 4) for _ in range(N)]
            Lv = lindbladian(H, D, g)
            Kn = string(lit * N)
            Rc, Rcinv = lr(I, Kn), lr(I, Kn.conj().T)                  # the colouring rho -> rho K^N
            ok_col &= not np.any(Rc @ Lv @ Rcinv + Lv.conj().T + 2 * sum(g) * np.eye(d * d))
            if Kl != 'I':
                Dn, Bn = string(D * N), Kn @ string(D * N)
                W, Winv = lr(Dn, Bn), lr(Dn.conj().T, Bn.conj().T)       # D^N rho K^N D^N, unitary
                ok_w &= not np.any(W @ Lv @ Winv + Lv + 2 * sum(g) * np.eye(d * d))
        check('D=%s, %s cell: the colouring rho -> rho %s^N gives -L^dag - 2 sigma exactly (40 random H, any body count, '
              'random per-site rates)' % (D, 'mother' if Kl == 'I' else 'letter-' + Kl, lit), ok_col)
        if Kl != 'I':
            check('D=%s, letter-%s cell: W(rho) = D^N rho %s^N D^N gives W L W^-1 = -L - 2 sigma exactly, same 40 H'
                  % (D, Kl, Kl), ok_w)
    fails_ok = True                                                    # controls: what must NOT hold
    for _ in range(10):
        g = [random.randint(1, 4) for _ in range(N)]
        Hd = random_cell_H(cells[D])
        Lv = lindbladian(Hd, D, g)
        for lit in (A[D], third):
            Kn = string(lit * N)
            fails_ok &= bool(np.any(lr(I, Kn) @ Lv @ lr(I, Kn.conj().T) + Lv.conj().T + 2 * sum(g) * np.eye(d * d)))
        Hm = random_cell_H(cells['I'])
        Lm = lindbladian(Hm, D, g)
        Dn, Bn = string(D * N), string(A[D] * N) @ string(D * N)
        fails_ok &= bool(np.any(lr(Dn, Bn) @ Lm @ lr(Dn.conj().T, Bn.conj().T) + Lm + 2 * sum(g) * np.eye(d * d)))
    check('D=%s controls: no lit letter colours the diagonal cell, and W fails on the mother cell' % D, fails_ok)

print('G4  the phases matter too: F108 Part 1\'s Pi_5bilinear flips by X like the canonical Pi_Z')
pi5 = {'I': ('X', 1), 'X': ('I', -1), 'Y': ('Z', 1j), 'Z': ('Y', -1j)}
Pm = letter_map_matrix(lambda c: pi5[c])
LDz = sum(to_pauli(lr(string(''.join('Z' if m == l else 'I' for m in range(N))),
                      string(''.join('Z' if m == l else 'I' for m in range(N))).conj().T) - lr(I, I)) for l in range(N))
check('Pi_5bilinear palindromizes the Z dissipator exactly', not np.any(Pm @ LDz @ Pm.conj().T + LDz + 2 * N * np.eye(4 ** N)))
M5 = {t: Pm @ LH[t] @ Pm.conj().T + LH[t] for t in window}
c5 = {}
for t1, t2, K, y in pairs:
    c5[(K, y)] = c5.get((K, y), 0) + (not np.any(M5[t1] + M5[t2]))
check('against Pi_5bilinear the X cell\'s truly pairs are its y_par = 1 half (21; XZZ is not truly) and the mother keeps'
      ' (45, 0), its XYZ permutations not truly: %s' % sorted((k, v) for k, v in c5.items() if v),
      c5.get(((1, 0), 1)) == 21 and c5.get(((1, 0), 0)) == 0 and c5.get(((0, 0), 0)) == 45
      and c5.get(((0, 0), 1), 0) == 0 and sum(v for (K, y), v in c5.items() if K not in ((0, 0), (1, 0))) == 0
      and np.any(M5[('X', 'Z', 'Z')]) and all(np.any(M5[p]) for p in itertools.permutations('XYZ')))

print('G5  all 256 phase choices on the X flip palindromize the Z dissipator; their X-cell truly sets run from none to all')
xcell = [t for t in window if klein(t) == KLEIN['X']]
sets, all_pal = set(), True
for ph in itertools.product([1, 1j, -1, -1j], repeat=4):
    mp = {'I': ('X', ph[0]), 'X': ('I', ph[1]), 'Y': ('Z', ph[2]), 'Z': ('Y', ph[3])}
    Pv = letter_map_matrix(lambda c: mp[c])
    all_pal &= not np.any(Pv @ LDz @ Pv.conj().T + LDz + 2 * N * np.eye(4 ** N))
    tru = tuple(sorted(t for t in xcell if not np.any(Pv @ LH[t] @ Pv.conj().T + LH[t])))
    sets.add((sum(t.count('Y') % 2 == 0 for t in tru), sum(t.count('Y') % 2 == 1 for t in tru)))
check('256 maps, all palindromize D[Z]; X-cell truly (y0 of 10, y1 of 6) takes the values %s, empty and all 16 among them'
      % sorted(sets), all_pal and (0, 0) in sets and (10, 6) in sets and len(sets) > 3)

print('\nALL PASS' if fails == 0 else '\n%d FAIL' % fails)
sys.exit(0 if fails == 0 else 1)
