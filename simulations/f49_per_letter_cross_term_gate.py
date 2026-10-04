"""Gate: the F49 cross term when the light comes along X, Y and Z, per site (depolarizing: equal letters per site).

Setting: Pauli dephasing with rate gamma_{l,C} >= 0 for letter C in {X, Y, Z} on site l, L_D = sum gamma_{l,C} D[C_l],
centred L_Dc = L_D + Gamma I with Gamma = sum_{l,C} gamma_{l,C}, and any Hamiltonian H = sum_tau c_tau T_tau of Pauli
strings T_tau (fields, two-site terms, k-site terms, any graph). Prediction, read on the letter cube:
  the centred rate of a letter P on site l is the character sum e_l(P) = sum_C gamma_{l,C} chi_C(P)
  (chi_C(P) = +1 if P commutes with C, -1 if not); a term whose letter on site l is a moves P to aP there, and
  e_l(P) + e_l(aP) = 2 gamma_{l,a} chi_a(P): only the light of the moving letter survives. Hence
    ||{L_H, L_Dc}||^2 = sum_tau ||L_H^tau||^2 [ 4 sum_{m not in supp tau} g_m^2 + B_tau ],   g_m^2 = sum_C gamma_{m,C}^2,
    B_tau = 4 gamma_{i,a}^2 (field a_i), 4 (gamma_{i,a} - gamma_{j,b})^2 (two-site a_i b_j), 4 sum_l gamma_{l,a_l}^2 (k >= 3),
    ||L_Dc||^2 = 4^N sum_m g_m^2,    ||L_H^tau||^2 = 2 4^N |c_tau|^2.
Special cases: F49 (uniform Z; every term balanced), F49b, F49c (X_i Z_j under Z: N-2 -> N-1), F49d (per-site Z:
(gamma_i - gamma_j)^2 only from ZZ terms), and uniform depolarizing (one rate for every letter on every site): every
two-site term has B = 0, so F49's constant sqrt((N-2)/(N 4^(N-1))) holds for EVERY two-site coupling, the nine letter
pairs included; per-site depolarizing rates p_l give every two-site term B = 4 (p_i - p_j)^2 instead.

Route: the anticommutator built as a superoperator on row-major vec(rho) from the Lindblad form itself (no Pauli-basis
shortcut), Frobenius norms as sums of re^2 + im^2 (unitary-invariant, so equal to the Pauli-basis norms), integer data,
every comparison with ==.
Checks:
  G1  random term sets (fields, two-site, k >= 3 strings) with random integer rates per site and letter, N = 3, 4.
  G2  the special cases: F49 (XXZ, XY, DM chains, uniform Z), F49c (X_i Z_j chain), F49d (per-site Z, XXZ chain),
      mixed bonds and the crossing term under per-site Z rates, uniform depolarizing (all nine two-site letter pairs,
      chain and ring, N = 3, 4, 5; every two-site term at N = 2), and the controls that do NOT give F49's constant:
      a field under depolarizing (B = 4 gamma_a^2 > 0) and per-site depolarizing (N = 3, p = (1, 2, 3): 1/42).
  G3  the moving-letter identity e(P) + e(aP) = 2 gamma_a chi_a(P) for every letter pair, random rates.
  G4  the commutator [L_H, L_Dc] loses the spectators: for a single string, its norm does not change when the rates
      off the string's support change (e(aP) - e(P) carries only the lights the move flips).
Run: python simulations/f49_per_letter_cross_term_gate.py   (about half a minute)"""
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


def string(s):
    m = np.array([[1.0 + 0j]])
    for c in s:
        m = np.kron(m, P[c])
    return m


def lr(A, B):
    return np.kron(A, B.T)


def dissipator(J):
    d = J.shape[0]
    I = np.eye(d)
    JdJ = J.conj().T @ J
    return lr(J, J.conj().T) - 0.5 * (lr(JdJ, I) + lr(I, JdJ))


def chi(C, a):
    return 1 if a in ('I', C) else -1


def norm2(A):
    return float((A.real ** 2 + A.imag ** 2).sum())


def build(N, terms, g):
    """terms: list of (coef, string); g[l][C]: integer rate of letter C on site l."""
    d = 2 ** N
    I = np.eye(d)
    H = sum(c * string(s) for c, s in terms)
    LH = -1j * (lr(H, I) - lr(I, H))
    LD = sum(g[l][C] * dissipator(string(''.join(C if m == l else 'I' for m in range(N))))
             for l in range(N) for C in 'XYZ' if g[l][C])
    Gamma = sum(g[l][C] for l in range(N) for C in 'XYZ')
    LDc = LD + Gamma * np.eye(d * d)
    return LH, LDc


def predicted(N, terms, g):
    gm2 = [sum(g[m][C] ** 2 for C in 'XYZ') for m in range(N)]
    total = 0
    for c, s in terms:
        supp = [l for l in range(N) if s[l] != 'I']
        spect = 4 * sum(gm2[m] for m in range(N) if m not in supp)
        if len(supp) == 1:
            B = 4 * g[supp[0]][s[supp[0]]] ** 2
        elif len(supp) == 2:
            B = 4 * (g[supp[0]][s[supp[0]]] - g[supp[1]][s[supp[1]]]) ** 2
        else:
            B = 4 * sum(g[l][s[l]] ** 2 for l in supp)
        total += 2 * 4 ** N * abs(c) ** 2 * (spect + B)
    return total, 4 ** N * sum(gm2)


FAILS = 0


def check(name, ok):
    global FAILS
    FAILS += (not ok)
    print('  %-104s %s' % (name, 'PASS' if ok else 'FAIL'))


def terms_on(N, pairs_letters, bonds, coef=1):
    out = []
    for (i, j) in bonds:
        for a, b in pairs_letters:
            out.append((coef, ''.join(a if m == i else (b if m == j else 'I') for m in range(N))))
    return out


def merge(terms):
    acc = {}
    for c, s in terms:
        acc[s] = acc.get(s, 0) + c
    return [(c, s) for s, c in acc.items() if c != 0]


print('G1  random terms, random rates per site and letter')
random.seed(20261004)
for N, trials in ((3, 40), (4, 12)):
    strings = [''.join(s) for s in itertools.product('IXYZ', repeat=N) if set(s) != {'I'}]
    ok_x = ok_d = True
    for t in range(trials):
        g = [{C: random.randint(0, 4) for C in 'XYZ'} for _ in range(N)]
        terms = [(random.choice([1, -1, 2, 3]), s) for s in random.sample(strings, random.randint(2, 7))]
        LH, LDc = build(N, terms, g)
        A = LH @ LDc + LDc @ LH
        rhs, rhsD = predicted(N, terms, g)
        ok_x &= norm2(A) == rhs
        ok_d &= norm2(LDc) == rhsD
    check('N=%d, %d random term sets: ||{L_H, L_Dc}||^2 = the per-letter formula' % (N, trials), ok_x)
    check('N=%d, same sets: ||L_Dc||^2 = 4^N sum_m g_m^2' % N, ok_d)

print('G2  special cases')


def ratio_matches(N, terms, g, num, den):
    """R^2 = ||A||^2 / (||L_H||^2 ||L_Dc||^2) equals num/den, compared cross-multiplied."""
    LH, LDc = build(N, terms, g)
    A = LH @ LDc + LDc @ LH
    return norm2(A) * den == num * norm2(LH) * norm2(LDc)


chain = lambda N: [(l, l + 1) for l in range(N - 1)]
ring = lambda N: [(l, (l + 1) % N) for l in range(N)]
for N in (3, 4, 5):
    gz = [{'X': 0, 'Y': 0, 'Z': 2} for _ in range(N)]
    ok = all(ratio_matches(N, merge(terms_on(N, L, chain(N))), gz, N - 2, N * 4 ** (N - 1))
             for L in ([('X', 'X'), ('Y', 'Y'), ('Z', 'Z'), ('Z', 'Z')], [('X', 'X'), ('Y', 'Y')], [('X', 'Y')]))
    ok &= ratio_matches(N, [(c, s) for c, s in [(1, x) for _, x in terms_on(N, [('X', 'Y')], chain(N))]] +
                        [(-1, x) for _, x in terms_on(N, [('Y', 'X')], chain(N))], gz, N - 2, N * 4 ** (N - 1))
    check('F49, N=%d, uniform Z: XXZ (Delta=2), XY, XY bond alone, DM XY-YX: R^2 = (N-2)/(N 4^(N-1))' % N, ok)
    check('F49c, N=%d, uniform Z, X_i Z_j chain: R^2 = (N-1)/(N 4^(N-1))' % N,
          ratio_matches(N, terms_on(N, [('X', 'Z')], chain(N)), gz, N - 1, N * 4 ** (N - 1)))
for N in (3, 4):
    gs = [{'X': 0, 'Y': 0, 'Z': z} for z in [1, 3, 2, 5][:N]]
    terms = merge(terms_on(N, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z'), ('Z', 'Z')], chain(N)))
    LH, LDc = build(N, terms, gs)
    A = LH @ LDc + LDc @ LH
    gam = [x['Z'] for x in gs]
    spect = sum(2 * 4 ** N * abs(c) ** 2 * 4 * sum(gam[m] ** 2 for m in range(N) if s[m] == 'I') for c, s in terms)
    asym = sum(2 * 4 ** N * abs(c) ** 2 * 4 * (gam[i] - gam[i + 1]) ** 2
               for c, s in terms for i in range(N - 1) if s[i] == 'Z' and s[i + 1] == 'Z')
    check('F49d, N=%d, per-site Z rates %s, XXZ chain: spectator part + (gamma_i - gamma_j)^2 from ZZ only' % (N, gam),
          norm2(A) == spect + asym)
for N in (3, 4):
    gz = [{'X': 0, 'Y': 0, 'Z': 2} for _ in range(N)]
    mixed = merge(terms_on(N, [('X', 'X'), ('Y', 'Y')], chain(N), 2) + terms_on(N, [('X', 'Z')], chain(N), 3))
    LH, LDc = build(N, mixed, gz)
    A = LH @ LDc + LDc @ LH
    rhs, _ = predicted(N, mixed, gz)
    bal = sum(2 * 4 ** N * abs(c) ** 2 for c, s in mixed if 'Z' not in s)
    crs = sum(2 * 4 ** N * abs(c) ** 2 for c, s in mixed if 'Z' in s)
    check('mixed bond, N=%d, uniform Z: 2(XX+YY) + 3 XZ on each bond = (N-2) part + (N-1) part, term by term' % N,
          norm2(A) == rhs == 4 * 4 * ((N - 2) * bal + (N - 1) * crs))
    gs = [{'X': 0, 'Y': 0, 'Z': z} for z in [1, 3, 2, 5][:N]]
    crossing = terms_on(N, [('X', 'Z')], chain(N))
    LH, LDc = build(N, crossing, gs)
    A = LH @ LDc + LDc @ LH
    gam = [x['Z'] for x in gs]
    pred = sum(2 * 4 ** N * 4 * (sum(gam[m] ** 2 for m in range(N) if s[m] == 'I') + gam[s.index('Z')] ** 2)
               for c, s in crossing)
    check('per-site Z rates %s, X_i Z_j chain (F49d\'s out-of-scope case): bond term 4 gamma_j^2' % gam, norm2(A) == pred)
gd2 = [{'X': 1, 'Y': 1, 'Z': 1} for _ in range(2)]
ok = all(norm2((lambda LH, LDc: LH @ LDc + LDc @ LH)(*build(2, [(1, a + b)], gd2))) == 0.0
         for a, b in itertools.product('XYZ', repeat=2))
check('N=2, depolarizing: every one of the nine two-site terms has {L_H, L_Dc} = 0 exactly (F48 for every bond)', ok)
for N in (3, 4, 5):
    gd = [{'X': 1, 'Y': 1, 'Z': 1} for _ in range(N)]
    for topo, bonds in (('chain', chain(N)), ('ring', ring(N))):
        if topo == 'ring' and N < 3:
            continue
        ok = all(ratio_matches(N, terms_on(N, [(a, b)], bonds), gd, N - 2, N * 4 ** (N - 1))
                 for a, b in itertools.product('XYZ', repeat=2))
        mix = merge(terms_on(N, [('X', 'Z')], bonds, 2) + terms_on(N, [('Y', 'Y')], bonds, 1) +
                    terms_on(N, [('Z', 'X')], bonds, -1))
        ok &= ratio_matches(N, mix, gd, N - 2, N * 4 ** (N - 1))
        check('depolarizing (equal letter rates), N=%d %s: all nine two-site letter pairs and a mix give F49\'s constant'
              % (N, topo), ok)
    fields = [(1, ''.join('X' if m == l else 'I' for m in range(N))) for l in range(N)]
    check('control, N=%d, depolarizing: an X field does NOT give F49\'s constant (B = 4 gamma_X^2 > 0)' % N,
          not ratio_matches(N, fields, gd, N - 2, N * 4 ** (N - 1)))
gp = [{C: p for C in 'XYZ'} for p in (1, 2, 3)]
ok = all(ratio_matches(3, terms_on(3, [L], chain(3)), gp, 1, 42) and not ratio_matches(3, terms_on(3, [L], chain(3)), gp, 1, 48)
         for L in (('X', 'X'), ('X', 'Z'), ('Z', 'Z')))
check('control, N=3, per-site depolarizing p = (1, 2, 3): XX, XZ and ZZ chains all give R^2 = 1/42, not F49\'s 1/48', ok)

print('G3  the moving-letter identity')
ok = True
for _ in range(200):
    g = {C: random.randint(0, 9) for C in 'XYZ'}
    e = {Pl: sum(g[C] * chi(C, Pl) for C in 'XYZ') for Pl in 'IXYZ'}
    for Pl, a in itertools.product('IXYZ', 'XYZ'):
        prod = string(a) @ string(Pl)
        aP = [c for c in 'IXYZ' if abs(np.trace(string(c).conj().T @ prod)) == 2][0]
        ok &= e[Pl] + e[aP] == 2 * g[a] * chi(a, Pl)
check('200 random rate triples: e(P) + e(aP) = 2 gamma_a chi_a(P) for every letter P and moving letter a', ok)

print('G4  the commutator is the complement: it loses the spectators')
ok = True
for _ in range(30):
    N = 4
    s = random.choice([x for x in (''.join(t) for t in itertools.product('IXYZ', repeat=N)) if 1 <= sum(c != 'I' for c in x) <= 3])
    g1 = [{C: random.randint(0, 4) for C in 'XYZ'} for _ in range(N)]
    g2 = [dict(g1[m]) if s[m] != 'I' else {C: random.randint(0, 4) for C in 'XYZ'} for m in range(N)]
    vals = []
    for g in (g1, g2):
        LH, LDc = build(N, [(1, s)], g)
        vals.append(norm2(LH @ LDc - LDc @ LH))
    ok &= vals[0] == vals[1]
check('30 random strings at N=4: ||[L_H, L_Dc]||^2 does not change when the rates off the support change', ok)

print('\nALL PASS' if FAILS == 0 else '\n%d FAIL' % FAILS)
sys.exit(0 if FAILS == 0 else 1)
