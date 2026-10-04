"""Gate: which part of a dissipator is odd under the three turns (F82's open generalization, read on the letter cube).

A jump J = sum_i c_i P_i (Pauli strings) gives D[J](rho) = sum_{i,j} c_i conj(c_j) (P_i rho P_j^dag - 1/2 {P_j^dag P_i, rho}).
The (i, j) term carries a Pauli string sigma to a multiple of P_i sigma P_j^dag, so it moves the string's count k_Q
(Q = X, Y or Z: the number of letters that anticommute with Q) by k_Q(P_i) + k_Q(P_j) mod 2. The turn by Q^N acts on a
string as the sign (-1)^{k_Q}; for Q = X this is F81/F82's Pi^2 with bit_b = k_X mod 2. Hence:
  the Q-odd part of D[J] is exactly the sum of its (i, j) terms with k_Q(P_i) + k_Q(P_j) odd.
A single Pauli string jump has only i = j: odd part 0 (X-, Y-, Z-noise, ZZ-dephasing; depolarizing as their sum), and
F81's identity Pi M Pi^-1 = M - 2 L_{H_odd} holds exactly under such noise. T1 (sigma- = (X + iY)/2) pairs X (k_X = 0)
with Y (k_X = 1); its odd part is F82's single (Z, I) entry per site, 4^(N-1) entries of size gamma_l.

Exact throughout: integer and Gaussian-integer data, superoperators on row-major vec(rho), every comparison with ==.
Checks:
  G1  N=3, all 63 single-string jumps: odd part 0 under each of the three turns.
  G2  N=3, 60 random jumps per turn (2-4 strings, Gaussian-integer coefficients): odd part = the odd-pair formula.
  G3  T1 at N=2, 3, 4, three rate profiles: ||D_odd||_F^2 = 4^(N-1) sum_l gamma_l^2 (F82's closed form, squared).
  G4  F81's identity under Z-dephasing plus X-noise, Y-noise, ZZ-dephasing or depolarizing, N=3 chain, H = XY + YX
      (Pi^2-odd) and H = XX + XY (mixed): Pi M Pi^-1 - M = -2 L_{H_odd} exactly; with T1 added it fails by exactly
      -2 D_{T1,odd}.
  G5  T1's odd pairs split into sandwich and anticommutator: each alone moves I -> Z and Z -> I, the sum only I -> Z.
  G6  which single jumps register, both sides of the line: X + Z (Hermitian, unital: odd part = its X/Z cross terms,
      the sqrt2 gamma the F84 proof pins at N = 1, sqrt2 gamma 2^(N-1) on one site of N) and X + Y register; Y + Z (no
      Pauli axis, both letters bit_b 1), I + Z and the dephasing jump (I - Z)/2 (identity components), and decay in the
      x basis (Z - iY)/2 do not; relaxation toward a y eigenstate (Z + iX)/2 does; sigma- and sigma+ at equal rates
      each register and cancel as a pair.
  G7  pair decay sigma- x sigma- and sigma+ x sigma+ on one pair, N = 2, 3, 4: ||D_odd||^2 = (3/2) 4^(N-2), the
      6 4^(N-3) of simulations/2qubit_dissipator_exploration.py, every odd entry a one-site I <-> Z move on that pair.
Run: python simulations/f82_dissipator_odd_part_gate.py   (under a second)"""
import itertools
import random
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

P = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], complex),
     'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1]).astype(complex)}
K_AXIS = {'Z': 0, 'X': 1, 'Y': 2}
CORNER = {'I': (0, 0, 0), 'X': (1, 0, 1), 'Y': (1, 1, 0), 'Z': (0, 1, 1)}   # (k_Z, k_X, k_Y) of one letter


def string(s):
    m = np.array([[1.0 + 0j]])
    for c in s:
        m = np.kron(m, P[c])
    return m


def lr(A, B):
    """rho -> A rho B on row-major vec(rho)."""
    return np.kron(A, B.T)


def dissipator(J):
    d = J.shape[0]
    I = np.eye(d)
    JdJ = J.conj().T @ J
    return lr(J, J.conj().T) - 0.5 * (lr(JdJ, I) + lr(I, JdJ))


def turn(N, Q):
    U = string(Q * N)
    return lr(U, U.conj().T)


def k(s, Q):
    return sum(CORNER[c][K_AXIS[Q]] for c in s)


def site(N, l, c):
    return ''.join(c if m == l else 'I' for m in range(N))


FAILS = 0


def check(name, ok):
    global FAILS
    FAILS += (not ok)
    print('  %-104s %s' % (name, 'PASS' if ok else 'FAIL'))


print('G1  single-string jumps')
N = 3
strings = [''.join(s) for s in itertools.product('IXYZ', repeat=N) if set(s) != {'I'}]
for Q in 'XYZ':
    T = turn(N, Q)
    bad = [s for s in strings if np.any(T @ dissipator(string(s)) @ T - dissipator(string(s)))]
    check('N=3, all %d single-string jumps, turn by %s^N: odd part exactly 0 (exceptions: %d)' % (len(strings), Q, len(bad)),
          not bad)

print('G2  sums of strings: the odd part is the odd-pair terms')
random.seed(20261004)
for Q in 'XYZ':
    T = turn(N, Q)
    I = np.eye(2 ** N)
    ok = True
    for _ in range(60):
        comps = random.sample(strings + ['I' * N], random.randint(2, 4))
        cs = [complex(random.randint(-3, 3), random.randint(-3, 3)) for _ in comps]
        J = sum(c * string(s) for c, s in zip(cs, comps))
        D = dissipator(J)
        odd = (D - T @ D @ T) / 2
        pred = np.zeros_like(D)
        for (ci, si), (cj, sj) in itertools.product(zip(cs, comps), repeat=2):
            if (k(si, Q) + k(sj, Q)) % 2 == 1:
                Pi_, Pj = string(si), string(sj)
                PjdPi = Pj.conj().T @ Pi_
                pred += ci * np.conj(cj) * (lr(Pi_, Pj.conj().T) - 0.5 * (lr(PjdPi, I) + lr(I, PjdPi)))
        ok &= not np.any(odd - pred)
    check('N=3, 60 random jumps, turn by %s^N: odd part = sum of the (i, j) terms with k_%s(P_i) + k_%s(P_j) odd' % (Q, Q, Q),
          ok)

print('G3  T1: one move per site, F82 closed form')
sm = np.array([[0, 1], [0, 0]], complex)              # framework convention sigma- = (X + iY)/2
for N in (2, 3, 4):
    T = turn(N, 'X')
    for rates in ([1] * N, [3, 0, 1, 2][:N], [2, 5, 0, 1][:N]):
        D = sum(r * dissipator(np.kron(np.kron(np.eye(2 ** l), sm), np.eye(2 ** (N - l - 1))))
                for l, r in enumerate(rates))
        odd = (D - T @ D @ T) / 2
        lhs = float((odd.real ** 2 + odd.imag ** 2).sum())
        rhs = 4 ** (N - 1) * sum(r * r for r in rates)
        check('N=%d, T1 rates %s: ||D_odd||^2 = %g = 4^(N-1) sum gamma_l^2 = %d' % (N, rates, lhs, rhs), lhs == rhs)

print('G4  F81 under Pauli-channel noise, F82 under T1')
KLEIN = {'I': (0, 0), 'X': (1, 0), 'Z': (0, 1), 'Y': (1, 1)}
NAME = {v: c for c, v in KLEIN.items()}
N = 3
d = 2 ** N
labels = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
idx = {s: i for i, s in enumerate(labels)}
V = np.array([string(s).reshape(-1) for s in labels]).T


def to_pauli(S):
    return V.conj().T @ S @ V / d


Pi = np.zeros((4 ** N, 4 ** N), complex)               # framework.symmetry.pi_action, dephase letter Z
for s in labels:
    ph, out = 1, ''
    for c in s:
        a, b = KLEIN[c]
        out += NAME[(1 - a, b)]
        ph *= (1j if b == 1 else 1)
    Pi[idx[out], idx[s]] = ph
I = np.eye(d)


def bond_H(letters):
    return sum(string(''.join(letters[0] if m == l else (letters[1] if m == l + 1 else 'I') for m in range(N)))
               for l in range(N - 1))


def LH(H):
    return -1j * (lr(H, I) - lr(I, H))


Dz = sum(dissipator(string(site(N, l, 'Z'))) for l in range(N))
noises = {
    'X-noise': sum(dissipator(string(site(N, l, 'X'))) for l in range(N)),
    'Y-noise': sum(dissipator(string(site(N, l, 'Y'))) for l in range(N)),
    'ZZ-dephasing': sum(dissipator(string(''.join('Z' if m in (l, l + 1) else 'I' for m in range(N)))) for l in range(N - 1)),
    'depolarizing': sum(dissipator(string(site(N, l, c))) for l in range(N) for c in 'XYZ'),
}
for hname, terms in (('XY + YX', [('X', 'Y'), ('Y', 'X')]), ('XX + XY', [('X', 'X'), ('X', 'Y')])):
    H = sum(bond_H(t) for t in terms)
    Hodd = sum(bond_H(t) for t in terms if (KLEIN[t[0]][1] + KLEIN[t[1]][1]) % 2 == 1)
    for nname, Dn in noises.items():
        L = to_pauli(LH(H) + Dz + 2 * Dn)
        M = Pi @ L @ Pi.conj().T + L
        lhs = Pi @ M @ Pi.conj().T - M
        check('H = %s, Z-dephasing + %s: Pi M Pi^-1 - M = -2 L_{H_odd} exactly' % (hname, nname),
              not np.any(lhs + 2 * to_pauli(LH(Hodd))))
    T1 = sum(3 * dissipator(np.kron(np.kron(np.eye(2 ** l), sm), np.eye(2 ** (N - l - 1)))) for l in range(N))
    T = turn(N, 'X')
    T1odd = (T1 - T @ T1 @ T) / 2
    L = to_pauli(LH(H) + Dz + T1)
    M = Pi @ L @ Pi.conj().T + L
    lhs = Pi @ M @ Pi.conj().T - M
    check('H = %s, Z-dephasing + T1: Pi M Pi^-1 - M = -2 L_{H_odd} - 2 D_{T1,odd} exactly (F82)' % hname,
          not np.any(lhs + 2 * to_pauli(LH(Hodd)) + 2 * to_pauli(T1odd)))

print('G5  the two pieces of T1\'s odd pairs: each moves both ways, only their sum is the single (Z, I) entry')
basis1 = ['I', 'X', 'Y', 'Z']
V1 = np.array([P[c].reshape(-1) for c in basis1]).T


def entries_1(S):
    Mp = V1.conj().T @ S @ V1 / 2
    return {(basis1[i], basis1[j]) for i, j in zip(*np.nonzero(Mp))}


cX, cY = 0.5, 0.5j                                     # sigma- = (X + iY)/2
sandwich = np.zeros((4, 4), complex)
anti = np.zeros((4, 4), complex)
for (ci, si), (cj, sj) in (((cX, 'X'), (cY, 'Y')), ((cY, 'Y'), (cX, 'X'))):
    Pi_, Pj = P[si], P[sj]
    PjdPi = Pj.conj().T @ Pi_
    sandwich += ci * np.conj(cj) * lr(Pi_, Pj.conj().T)
    anti += -0.5 * ci * np.conj(cj) * (lr(PjdPi, np.eye(2)) + lr(np.eye(2), PjdPi))
check('sandwich alone: entries %s (output, input)' % sorted(entries_1(sandwich)), entries_1(sandwich) == {('Z', 'I'), ('I', 'Z')})
check('anticommutator alone: entries %s' % sorted(entries_1(anti)), entries_1(anti) == {('Z', 'I'), ('I', 'Z')})
check('their sum: entries %s, the single move I -> Z' % sorted(entries_1(sandwich + anti)),
      entries_1(sandwich + anti) == {('Z', 'I')})

print('G6  which single jumps register: the parity rule draws the line, nothing else does')
T1x = lr(P['X'], P['X'])
D = dissipator(P['X'] + P['Z'])                        # D[(X + Z)/sqrt2] at rate gamma is (gamma/2) D[X + Z]
odd = (D - T1x @ D @ T1x) / 2
pair = sum(lr(P[a], P[b]) - 0.5 * (lr(P[b] @ P[a], np.eye(2)) + lr(np.eye(2), P[b] @ P[a])) for a, b in (('X', 'Z'), ('Z', 'X')))
check('N=1, J = X + Z: odd part = its X/Z cross terms, ||odd||^2 = %g = 8, i.e. sqrt2 gamma for (X+Z)/sqrt2 at rate gamma'
      % float((odd.real ** 2 + odd.imag ** 2).sum()), not np.any(odd - pair) and float((odd.real ** 2 + odd.imag ** 2).sum()) == 8.0)


def odd_norm2(J, N=1):
    T = turn(N, 'X')
    D = dissipator(J)
    o = (D - T @ D @ T) / 2
    return float((o.real ** 2 + o.imag ** 2).sum())


cases = [                                              # (name, jump, ||odd||^2 expected)
    ('X + Y (both an X part and a Y part)', P['X'] + P['Y'], 8.0),
    ('Y + Z (two letters of bit_b 1, no Pauli axis)', P['Y'] + P['Z'], 0.0),
    ('I + Z (an identity component, real)', P['I'] + P['Z'], 0.0),
    ('n = (I - Z)/2, the textbook dephasing jump', (P['I'] - P['Z']) / 2, 0.0),
    ('(Z - iY)/2 = |+><-|, decay in the x basis', (P['Z'] - 1j * P['Y']) / 2, 0.0),
    ('(Z + iX)/2, relaxation toward a y eigenstate', (P['Z'] + 1j * P['X']) / 2, 1.0),
]
for name, J, expect in cases:
    v = odd_norm2(J)
    check('N=1, %s: ||odd||^2 = %g (expected %g)' % (name, v, expect), v == expect)
for N in (2, 3):
    J = np.kron(P['X'] + P['Z'], np.eye(2 ** (N - 1)))
    v = odd_norm2(J, N)
    check('N=%d, X + Z on one site: ||odd||^2 = %g = 8 * 4^(N-1), i.e. sqrt2 gamma 2^(N-1) for the normalized axis' % (N, v),
          v == 8.0 * 4 ** (N - 1))
smp = sm.conj().T
D = dissipator(sm) + dissipator(smp)
o = (D - turn(1, 'X') @ D @ turn(1, 'X')) / 2
check('N=1, sigma- and sigma+ at equal rates: the two odd parts cancel exactly (each jump mixes, the pair does not)',
      not np.any(o) and odd_norm2(sm) == 1.0 and odd_norm2(smp) == 1.0)

print('G7  pair decay sigma- x sigma- (and sigma+ x sigma+)')
for N in (2, 3, 4):
    T = turn(N, 'X')
    labs = [''.join(s) for s in itertools.product('IXYZ', repeat=N)]
    VN = np.array([string(s).reshape(-1) for s in labs]).T
    ok_norm = ok_move = True
    for l, op in itertools.product(range(N - 1), (sm, sm.conj().T)):
        J = np.kron(np.kron(np.eye(2 ** l), np.kron(op, op)), np.eye(2 ** (N - l - 2)))
        D = dissipator(J)
        odd = (D - T @ D @ T) / 2
        ok_norm &= float((odd.real ** 2 + odd.imag ** 2).sum()) == 1.5 * 4 ** (N - 2)
        Mp = VN.conj().T @ odd @ VN / 2 ** N
        for i, j in zip(*np.nonzero(Mp)):
            diff = [m for m in range(N) if labs[i][m] != labs[j][m]]
            ok_move &= len(diff) == 1 and {labs[i][diff[0]], labs[j][diff[0]]} == {'I', 'Z'} and diff[0] in (l, l + 1)
    check('N=%d, every pair of neighbours, both jumps: ||D_odd||^2 = (3/2) 4^(N-2) = 6 4^(N-3) exactly' % N, ok_norm)
    check('N=%d, every odd entry is a one-site I <-> Z move on the jump\'s own pair' % N, ok_move)

print('\nALL PASS' if FAILS == 0 else '\n%d FAIL' % FAILS)
sys.exit(0 if FAILS == 0 else 1)
