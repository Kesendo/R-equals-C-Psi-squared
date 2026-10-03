"""Gate for F63 turned on the letter cube (docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md).

S_P(rho) = P^N rho P^N for P in {X, Y, Z}; on a Pauli string S_P is the sign (-1)^{k_P}, k_P the string's
letter-cube coordinate (docs/THE_ONE_SQUARE.md section 7; fw.cube_coords). Pi^2 of F63 is S_X.

G1  [L, S_P] = 0 iff every Hamiltonian term has even k_P, whatever the Pauli-string jumps (N = 3: every single
    term and every pair of terms, three letters, four Pauli jump sets). Exact: all entries are small Gaussian
    integers, the commutator is compared with == 0.0.
G2  sigma^- = (X - iY)/2 on one site: breaks S_X and S_Y for every Hamiltonian, keeps S_Z exactly when the
    Hamiltonian is Z-even (its two letters share k_Z = 1 but differ in k_X and k_Y).
G1b targeted controls: a field along X keeps only the turn by X^N (the four sectors need F61 as well), the
    Dzyaloshinskii-Moriya bond YZ - ZY keeps the turn by X^N, D[I + Y] = D[Y].
G2b the jump condition is sufficient, not necessary: sigma^- with sigma^+ at equal rates and the pair {X+Y, X-Y}
    keep S_X; unequal rates break it.
G3  dim ker L per S_X sector (Pauli basis, w_YZ = k_X parity), exact integer rank (python-flint), Z-dephasing on
    one site: on the open chain, N = 2..5, every seat, XXX and XY, the kernel is (floor(N/2)+1, ceil(N/2)), the
    e_d count, exactly at the seats where F157 counts no blind state, and larger at the others (G3a); the cells
    of the document's table for chain, ring and star (G3b).
G4  the Re = 0 class, exactly: the largest L_H-invariant subspace inside ker D (Hermitian jumps), an integer
    rank. It equals the kernel, N + 1, at the sighted seats checked (end seats N = 2..5, interior seats N = 4, 5);
    the pairs (Re = 0 class, kernel) at the blind and graph-symmetric seats are pinned. Requires python-flint.

Run: python simulations/f63_three_turns_gate.py
"""
import itertools
import sys
import numpy as np
import flint

PAULI = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
SIGMA_MINUS = np.array([[0, 0], [1, 0]])
FAILS = []


def check(name, ok, detail=''):
    print(('PASS ' if ok else 'FAIL ') + name + ((': ' + detail) if detail else ''))
    if not ok:
        FAILS.append(name)


def site_op(N, ops):
    m = np.array([[1.0 + 0j]])
    for l in range(N):
        o = ops.get(l, 'I')
        m = np.kron(m, PAULI[o] if isinstance(o, str) else o)
    return m


def string(s):
    return site_op(len(s), dict(enumerate(s)))


def k(s, p):
    return sum(1 for c in s if c not in ('I', p))


def liouvillian(H, jumps):
    d = H.shape[0]
    I = np.eye(d)
    L = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for J in jumps:
        JJ = J.conj().T @ J
        L = L + np.kron(J, J.conj()) - 0.5 * (np.kron(JJ, I) + np.kron(I, JJ.T))
    return L


def conj_super(U):
    return np.kron(U, U.conj())


# G1 and G2 -------------------------------------------------------------------------------------------------
N = 3
strings = [''.join(s) for s in itertools.product('IXYZ', repeat=N) if set(s) != {'I'}]
jumpsets = {'Z0': ['ZII'], 'X0,Y2': ['XII', 'IIY'], 'Z0,X1,Y2': ['ZII', 'IXI', 'IIY'], 'XYZ on 1': ['IXI', 'IYI', 'IZI']}
sm = site_op(N, {0: SIGMA_MINUS})
termsets = [[t] for t in strings] + [list(c) for c in itertools.combinations(strings, 2)]
for p in 'XYZ':
    S = conj_super(string(p * N))
    agree = total = 0
    sm_keep = sm_even = sm_wrong = 0
    for terms in termsets:
        H = sum(string(t) for t in terms)
        even = all(k(t, p) % 2 == 0 for t in terms)
        for js in jumpsets.values():
            L = liouvillian(H, [string(j) for j in js])
            total += 1
            agree += ((np.abs(L @ S - S @ L).max() == 0.0) == even)
        L = liouvillian(H, [sm])
        keeps = np.abs(L @ S - S @ L).max() == 0.0
        sm_keep += keeps
        sm_even += even
        sm_wrong += keeps != (even and p == 'Z')
    check('G1 S_%s: commutes iff every term has even k_%s (%d term sets x %d jump sets)' % (p, p, len(termsets), len(jumpsets)),
          agree == total, '%d/%d agree' % (agree, total))
    check('G2 S_%s under sigma^-: kept %d of %d %s-even Hamiltonians' % (p, sm_keep, sm_even, p),
          sm_wrong == 0 and sm_keep == (sm_even if p == 'Z' else 0))

# G2b: the jump condition is sufficient, not necessary. Conjugation can permute a set of dissipators whose jumps
# each mix parities: sigma^- with sigma^+ at EQUAL rates (= half of D[X] + D[Y]), and the pair {X + Y, X - Y}.
# Unequal rates pick a direction and break it. Heisenberg XXX on the chain, jumps on site 0.
Hh = sum(site_op(N, {a: c, a + 1: c}) for a in range(N - 1) for c in 'XYZ')
SX = conj_super(string('X' * N))
sp_ = site_op(N, {0: SIGMA_MINUS.T})
# Expected (turn by X^N, by Z^N, by Y^N), the rows of the document's table.
for name, js, want in [('sigma^- and sigma^+ at equal rates, one site', [sm, sp_], (True, True, True)),
                       ('sigma^- and sigma^+ at rates 1 and 1/4', [sm, 0.5 * sp_], (False, True, False)),
                       ('{X + Y, X - Y}', [string('XII') + string('YII'), string('XII') - string('YII')], (True, True, True))]:
    L = liouvillian(Hh, js)
    got = []
    for p in 'XZY':
        S = conj_super(string(p * N))
        got.append(np.abs(L @ S - S @ L).max() == 0.0)
    check('G2b %s: turns by X, Z, Y kept = %s' % (name, want), tuple(got) == want, 'got %s' % (tuple(got),))
# G1b, three targeted controls (handed over by an outside review, Codex, and recomputed here): at N = 2 with one Z
# jump, a field along X keeps the turn by X^N and breaks the other two, so the four sectors need F61 as well;
# the Dzyaloshinskii-Moriya bond YZ - ZY (along x) keeps the turn by X^N; and D[I + Y] = D[Y].
Z2 = [string('ZI')]
for name, H2, want in [('field along X', string('XI'), (True, False, False)),
                       ('DM bond YZ - ZY', string('YZ') - string('ZY'), (True, False, False))]:
    L = liouvillian(H2, Z2)
    got = tuple(bool(np.abs(L @ conj_super(string(p * 2)) - conj_super(string(p * 2)) @ L).max() == 0.0) for p in 'XZY')
    check('G1b N=2 %s, Z jump: turns by X, Z, Y kept = %s' % (name, want), got == want, 'got %s' % (got,))
L1, L2 = liouvillian(np.zeros((4, 4)), [string('II') + string('YI')]), liouvillian(np.zeros((4, 4)), [string('YI')])
check('G1b D[I + Y] = D[Y] exactly', np.abs(L1 - L2).max() == 0.0, 'max |L1 - L2| = %g' % np.abs(L1 - L2).max())
# An identity part hides a mixture: D[Z + I] = D[Z], so the single-jump rule needs a traceless jump.
L1, L2 = liouvillian(Hh, [string('ZII') + string('III')]), liouvillian(Hh, [string('ZII')])
check('G2b D[Z + I] = D[Z] exactly', np.abs(L1 - L2).max() == 0.0, 'max |L1 - L2| = %g' % np.abs(L1 - L2).max())


# G3 --------------------------------------------------------------------------------------------------------
def kernel_per_sector(N, bonds, letters, seat):
    """dim ker L in the w_YZ-even and w_YZ-odd sectors (L is block diagonal there by G1)."""
    H = sum(site_op(N, {a: c, b: c}) for a, b in bonds for c in letters)
    L = liouvillian(H, [site_op(N, {seat: 'Z'})])
    labels = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
    V = np.array([string(s).reshape(-1) for s in labels]).T
    M = (V.conj().T @ L @ V) / 2 ** N  # V has entries 0, +-1, +-i; dividing by a power of two stays exact
    R = np.rint(M.real * 2)
    # The Pauli-basis matrix is real with half-integer entries for these models; anything else is a defect.
    if not (np.abs(M.imag).max() == 0.0 and np.abs(M.real * 2 - R).max() == 0.0):
        raise SystemExit('Pauli-basis matrix not exactly real half-integer at N=%d' % N)
    out = []
    for parity in (0, 1):
        idx = [i for i, s in enumerate(labels) if k(s, 'X') % 2 == parity]
        B = flint.fmpz_mat(R[np.ix_(idx, idx)].astype(int).tolist())
        out.append(B.nrows() - B.rank())
    return tuple(out)


def chain(N):
    return [(l, l + 1) for l in range(N - 1)]


def blind(N, seat, letters):
    """F157 on the uniform chain, seats counted from 0 (compute/MirrorWorld/BlindSeat.cs)."""
    from math import gcd
    return (gcd(2 * seat + 1, N) - 1) // 2 if letters == 'XYZ' else gcd(seat + 1, N + 1) - 1


# G3a: on the open chain the kernel is the e_d count, (floor(N/2)+1, ceil(N/2)) per sector, at exactly the seats
# where F157 counts no blind state; at a blind seat it is larger. Every seat, N = 2..5, both books.
for letters, tag in (('XYZ', 'XXX'), ('XY', 'XY')):
    for N in range(2, 6):
        rows = []
        for seat in range(N):
            got = kernel_per_sector(N, chain(N), letters, seat)
            plain = got == (N // 2 + 1, (N + 1) // 2)
            rows.append((seat, got, blind(N, seat, letters), plain))
        ok = all(plain == (b == 0) and (plain or sum(got) > N + 1) for _, got, b, plain in rows)
        check('G3a %s chain N=%d: kernel = e_d count iff F157 blind(seat) = 0' % (tag, N), ok,
              '; '.join('seat %d %s blind %d' % (s, g, b) for s, g, b, _ in rows))

# G3b: the cells of the document's table, exact.
cases = [('chain centre', 3, chain(3), 1, (3, 3), (6, 6)), ('chain centre', 5, chain(5), 2, (6, 6), (12, 12)),
         ('chain seat 1', 5, chain(5), 1, (3, 3), (10, 10)),
         ('ring', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], 0, (7, 6), (12, 10)),
         ('star hub', 4, [(0, 1), (0, 2), (0, 3)], 0, (11, 6), (16, 11)),
         ('star leaf', 4, [(0, 1), (0, 2), (0, 3)], 1, (5, 3), (7, 5))]
for name, N, bonds, seat, want_xxx, want_xy in cases:
    for letters, tag, want in (('XYZ', 'XXX', want_xxx), ('XY', 'XY', want_xy)):
        got = kernel_per_sector(N, bonds, letters, seat)
        check('G3b %s %s N=%d: kernel per sector' % (tag, name, N), got == want, '%s, want %s' % (got, want))


# G4: the Re = 0 class against the kernel, by an exact route. With Hermitian jumps L = L_H + D in the Pauli basis,
# L_H antisymmetric and D diagonal (-2 gamma on strings with X or Y at the dephased seat, 0 otherwise), so
# Re(lambda) = <v, D v> / <v, v> for an eigenvector v: the undamped modes span the largest L_H-invariant subspace
# inside ker D, W = {v : P L_H^j v = 0 for all j}, P the coordinates of the damped strings. dim W is an exact integer
# rank (rows stacked until the rank stops growing). The eigensolver count of |Re| below 1e-8 is printed beside it
# as a second route, not gated.
def undamped_dim(N, bonds, letters, seat):
    H = sum(site_op(N, {a: c, b: c}) for a, b in bonds for c in letters)
    LH = liouvillian(H, [])
    labels = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
    V = np.array([string(s).reshape(-1) for s in labels]).T
    M = (V.conj().T @ LH @ V) / 2 ** N
    R = np.rint(M.real)
    if not (np.abs(M.imag).max() == 0.0 and np.abs(M.real - R).max() == 0.0):
        raise SystemExit('L_H not an exact integer matrix in the Pauli basis at N=%d' % N)
    damped = [i for i, s in enumerate(labels) if s[seat] in 'XY']
    Mi = R.astype(np.int64)
    block = np.eye(len(labels), dtype=np.int64)[damped]
    rows, rank = None, -1
    while True:
        rows = block if rows is None else np.vstack([rows, block])
        F = flint.fmpz_mat(rows.tolist())
        r = F.rank()
        if r == rank:
            break
        rank = r
        E = F.rref()[0]
        nz = [i for i in range(E.nrows()) if any(E[i, j] != 0 for j in range(E.ncols()))]
        rows = np.array([[int(E[i, j]) for j in range(E.ncols())] for i in nz], dtype=object)
        block = (rows @ Mi.astype(object))
    L = liouvillian(H, [site_op(N, {seat: 'Z'})])
    re = np.abs(np.linalg.eigvals(L).real)
    return len(labels) - rank, int((re < 1e-8).sum())


# Sighted seats: Re = 0 class = kernel = N + 1. Elsewhere the exact pairs (Re = 0 class, kernel) are pinned; a
# blind seat can carry undamped oscillating modes beyond the kernel or none.
g4 = [('chain end', n, chain(n), n - 1, l, None) for n in (2, 3, 4, 5) for l in ('XYZ', 'XY')] + \
     [('chain seat 1', 4, chain(4), 1, l, None) for l in ('XYZ', 'XY')] + \
     [('chain seat 2', 4, chain(4), 2, l, None) for l in ('XYZ', 'XY')] + \
     [('chain seat 1', 5, chain(5), 1, 'XYZ', None), ('chain seat 1', 5, chain(5), 1, 'XY', (20, 20)),
      ('chain centre', 3, chain(3), 1, 'XYZ', (10, 6)), ('chain centre', 5, chain(5), 2, 'XY', (64, 24)),
      ('ring', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], 0, 'XYZ', (21, 13))]
for name, N, bonds, seat, letters, want in g4:
    exact, eig = undamped_dim(N, bonds, letters, seat)
    kern = sum(kernel_per_sector(N, bonds, letters, seat))
    tag = 'XXX' if letters == 'XYZ' else 'XY'
    if want is None:
        ok, label = exact == kern == N + 1, 'Re = 0 class = kernel = N + 1 (sighted seat)'
    else:
        ok, label = (exact, kern) == want, '(Re = 0 class, kernel) = %s' % (want,)
    check('G4 %s %s N=%d: %s' % (tag, name, N, label), ok, 'exact %d (eigensolver %d), kernel %d' % (exact, eig, kern))

print('\n%d FAIL' % len(FAILS) if FAILS else '\nALL PASS')
sys.exit(1 if FAILS else 0)
