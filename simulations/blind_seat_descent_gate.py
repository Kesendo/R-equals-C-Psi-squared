"""Gate: the descent from W's top component (PROOF_ONE_DEPHASED_SEAT, section on blind dephased seats, bullet
"At N = 4, a seat that is not blind" under the paragraph "The sheet and blindness, at N = 3 and N = 4").

Setting: Heisenberg bonds sum J_ij sigma_i.sigma_j on a connected graph (integer couplings, signs free), the jump Z on
seat d = 0, planar integer fields on the sites != d, none on d; omega_j = f_x + i f_y, bit 1 = excited,
F+ = sum_j omega_j sigma^-_j raises the number of excitations, S+ = sum_j sigma^-_j. Exact or certified, no threshold.

1. The fence. For each sector n, the vectors with d excited in every term at the vacuum's energy E0 (the sum of the
   couplings), solved exactly (sympy nullspace). Predicted empty at n = 1, 2 and N - 1 where d is not blind, on drawn
   rows at N = 4, 5, 6; the other sectors are read. Beside it, [S+, H_b] = 0 exactly (a check of the bit-action
   bond matrix; on Heisenberg input it cannot fail). Controls: at a constructed blind row at each N, with a vector of
   ker L_J vanishing at d, the fence at N - 1 is not empty (the n = 1 fence cannot fail on a connected graph, and no
   control here reaches a nonempty fence at n = 2); the constructed N = 5 row of the bullet is not blind and its fence
   at n = 3 is one-dimensional, spanned exactly by X^N applied to |0>_d (x) s_12 (x) s_34 (s the singlet).
2. The hand steps at N = 4, symbolically, every operator built from the bit action, none typed: the determinant
   -2 conj(w_u w_v w_w) of F- on two excitations of the three sites other than d and its kernel at w_w = 0; the hop of
   d's bonds carrying |1>_d (x) u' to |0>_d (x) sum 2(J_dk u'_j + J_dj u'_k)|jk>; the |uv> entry
   2(J_vw conj(w_v) - J_uw conj(w_u)) of the bonds on y0 = conj(w_v)|uw> - conj(w_u)|vw>; and step (i)'s two maps,
   S+ on |1>_d (x) u' and S- on |0>_d (x) y, as 3 x 3 matrices with nonzero determinant.
3. Condition (I) at N = 4: the cyclic space alg(H_b, Z_d, F+)|vac> has dimension 16, by a rank modulo two primes
   (a rank modulo p is at most the rational rank, so 16 there certifies it), at drawn rows where d is not blind: the
   fields off one line with all three nonzero, with one of them zero, at rows whose two-excitation sector is reducible
   (commutant dimension above one, exactly; how many the draw finds is read), and class A with all three nonzero. Far elements homogeneous under Z^4
   counted by parity by a rank modulo p (a kernel modulo p is at least the rational one, so a zero is exact):
   predicted none off one line; on one line at most one by the rank and at least one by class A's colouring
   (m.sigma)^4, checked exactly, of even height, and none of odd height. Control: at a blind row with d bonded to two
   sites the vector |1>_0 (x) (J_01 e_1 - J_02 e_2) is exactly an eigenvector of H_b killed by F-, so the cyclic space
   has dimension at most 15 there, and the rank modulo p reads exactly 15.
4. At N = 5: drawn rows where d is not blind, conditions (I) and (II) and the far counts (on one line the 1 bounded
   below by class A's colouring, exactly); where both conditions hold the descent decides the row; the
   constructed row, where (II) fails at n = 3, at four field sets, two of them antiparallel within the pairs, and
   there the vacuum equation alone in the multiple of its fence vector, tau and chi, exactly: one-dimensional with the
   multiple nonzero at the two antiparallel sets, only zero at a third.
5. N = 4 on one line with a field zero: the matrix entries the steps (v) and (vi) use, symbolically from the bit
   action, and the steps' conclusion on an exhaustive grid (couplings in {-1, 0, 1, 2}, connected, d not blind,
   fields on the x line with a field zero on each site in turn): (I) by a full rank modulo two primes off the
   antiparallel locus, a short one on it, and on every locus row of the grid the cyclic dimension 14 and the far
   elements (even, odd) = (2, 0); on the locus, symbolically in the
   couplings and the field, the subspace Q = t (x) span{|01>, |10>} of 0 and 3 kept by the bonds and killed by F-, and
   the element h = (X0 X3 + Y0 Y3) (x) |t><t| commuting with H and anticommuting with Z_0, with a control that breaks
   it; at three locus rows the cyclic dimension 14 and the far elements (even, odd) = (2, 0), bounded below by class
   A's colouring and h (rank 2, so independent); the blind locus row at J03 = 0 read; controls off the locus.

Run: python simulations/blind_seat_descent_gate.py   (about a minute and a quarter; writes results/blind_seat_descent_gate.txt)
"""
import itertools
import random
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import pauli_product

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")   # the console's default code page cannot print every character below

OUT, FAILS = [], []


def say(s=""):
    print(s)
    OUT.append(s)


def check(cond, what):
    say(("  PASS  " if cond else "  FAIL  ") + what)
    if not cond:
        FAILS.append(what)


def _is_prime(q):
    return q > 1 and all(q % r for r in range(2, int(q ** 0.5) + 1))


def _prime_1mod4(start):
    q = start
    while not (_is_prime(q) and q % 4 == 1):
        q += 1
    return q


PRIMES = [_prime_1mod4(10 ** 9), _prime_1mod4(10 ** 9 + 100)]


def sqrt_m1(p):
    g = next(g for g in range(2, 200) if pow(g, (p - 1) // 2, p) == p - 1)
    return pow(g, (p - 1) // 4, p)


def rank_mod(M, p):
    M = np.array(M, dtype=object) % p
    M = np.array(M, dtype=np.int64)
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r


# ---- rows --------------------------------------------------------------------------------------------------------------
def connected(n, bonds):
    adj = {i: set() for i in range(n)}
    for (i, j), J in bonds.items():
        if J:
            adj[i].add(j); adj[j].add(i)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x] - seen:
            seen.add(y); st.append(y)
    return len(seen) == n


def laplacian(n, bonds):
    L = sp.zeros(n, n)
    for (i, j), J in bonds.items():
        L[i, j] -= J; L[j, i] -= J; L[i, i] += J; L[j, j] += J
    return L


def blind_dim(n, bonds, d=0):
    L = laplacian(n, bonds)
    v = sp.zeros(n, 1); v[d] = 1
    K = [v]
    for _ in range(n - 1):
        K.append(L * K[-1])
    return n - sp.Matrix.hstack(*K).rank()


def on_one_line(fields):
    vs = [v for v in fields.values() if v != (0, 0)]
    return all(a[0] * b[1] - a[1] * b[0] == 0 for a in vs for b in vs)


CPL = (-3, -2, -1, 1, 2, 3)
FLD = [(1, 0), (0, 1), (-1, 0), (3, 4), (4, -3), (1, 1), (-1, 2), (2, 0), (0, 2), (5, 1)]


def draw_bonds(n, rng):
    edges = list(itertools.combinations(range(n), 2))
    while True:
        b = {e: rng.choice(CPL) for e in edges if rng.random() < 0.6}
        if connected(n, b):
            return b


def key(n, bonds, fields=None):
    return (n, tuple(sorted(bonds.items())), tuple(sorted((fields or {}).items())))


# ---- the state space ---------------------------------------------------------------------------------------------------
def hb(n, bonds):
    D = 2 ** n
    M = sp.zeros(D, D)
    for x in range(D):
        for (i, j), J in bonds.items():
            bi, bj = (x >> i) & 1, (x >> j) & 1
            M[x, x] += J * (1 if bi == bj else -1)
            if bi != bj:
                M[x ^ (1 << i) ^ (1 << j), x] += 2 * J
    return M


def fplus(n, fields):
    D = 2 ** n
    M = sp.zeros(D, D)
    for j, (fx, fy) in fields.items():
        for x in range(D):
            if not (x >> j) & 1:
                M[x | (1 << j), x] += fx + sp.I * fy
    return M


def splus(n):
    return fplus(n, {j: (1, 0) for j in range(n)})


def zd(n, d=0):
    return sp.diag(*[(-1 if (x >> d) & 1 else 1) for x in range(2 ** n)])


def sector(n, k):
    return [x for x in range(2 ** n) if bin(x).count('1') == k]


def fence(n, bonds, H=None, d=0):
    """Basis of {v in sector k : H_b v = E0 v, bit d = 1 on the support}, for k = 1 .. N - 1, exactly."""
    H = hb(n, bonds) if H is None else H
    E0 = sum(bonds.values())
    out = {}
    for k in range(1, n):
        idx = [x for x in sector(n, k) if (x >> d) & 1]
        A = sp.Matrix([[H[a, b] - (E0 if a == b else 0) for b in idx] for a in range(2 ** n)])
        out[k] = [(idx, v) for v in A.nullspace()]
    return out


def cyclic_dim(n, bonds, fields):
    """dim alg(H_b, Z_d, F+)|vac> by a rank modulo each prime (a lower bound on the rational dimension); the max."""
    H = np.array(hb(n, bonds).tolist(), dtype=object)
    Z = np.array(zd(n).tolist(), dtype=object)
    best = 0
    for p in PRIMES:
        ip = sqrt_m1(p)
        Fp = np.zeros((2 ** n, 2 ** n), dtype=object)
        for j, (fx, fy) in fields.items():
            for x in range(2 ** n):
                if not (x >> j) & 1:
                    Fp[x | (1 << j), x] = (Fp[x | (1 << j), x] + fx + fy * ip) % p
        gens = [H % p, Z % p, Fp]
        basis, piv = [], []
        queue = [np.array([1] + [0] * (2 ** n - 1), dtype=object)]
        while queue:
            v = queue.pop() % p
            r = v.copy()
            for b, c in zip(basis, piv):
                if r[c]:
                    r = (r - r[c] * b) % p
            nz = [i for i in range(len(r)) if r[i]]
            if not nz:
                continue
            c = nz[0]
            r = (r * pow(int(r[c]), p - 2, p)) % p
            for t in range(len(basis)):
                if basis[t][c]:
                    basis[t] = (basis[t] - basis[t][c] * r) % p
            basis.append(r); piv.append(c)
            for G in gens:
                queue.append(G.dot(v) % p)
        best = max(best, len(basis))
    return best


# ---- far elements by parity (Pauli strings, rank modulo p) --------------------------------------------------------------
P = PRIMES[0]
I_P = sqrt_m1(P)
IPOW = [1, I_P, P - 1, P - I_P]


def anticomm(s, t):
    return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1


def h_terms(n, bonds, fields):
    t = {}
    for (i, j), J in bonds.items():
        for L in 'XYZ':
            s = ['I'] * n; s[i] = L; s[j] = L
            t[''.join(s)] = (J, 0)
    for j, (fx, fy) in fields.items():
        for L, c in (('X', fx), ('Y', fy)):
            if c:
                s = ['I'] * n; s[j] = L
                t[''.join(s)] = (c, 0)
    return t


def far_count(n, bonds, fields, parity, d=0):
    terms = h_terms(n, bonds, fields)
    span = [''.join(w) for w in itertools.product('IXYZ', repeat=n)
            if w[d] in 'XY' and sum(1 for a in w if a in 'XY') % 2 == parity]
    rowidx, entries = {}, {}
    for col, s in enumerate(span):
        for Pp, (re, im) in terms.items():
            if not anticomm(Pp, s):
                continue
            e, u = pauli_product(Pp, s)
            r = rowidx.setdefault(u, len(rowidx))
            entries[(r, col)] = (entries.get((r, col), 0) + 2 * (re + im * I_P) * IPOW[e]) % P
    M = np.zeros((len(rowidx), len(span)), dtype=np.int64)
    for (r, c), v in entries.items():
        M[r, c] = v
    return len(span) - rank_mod(M, P)


def colouring_exact(n, bonds, fields):
    """(m.sigma)^N for the common axis m: commutes with H, anticommutes with Z_d, exactly."""
    m = next(v for v in fields.values() if v != (0, 0))
    ms = sp.Matrix([[0, m[0] - sp.I * m[1]], [m[0] + sp.I * m[1], 0]])
    W = ms
    for _ in range(n - 1):
        W = sp.kronecker_product(W, ms)
    H = hb(n, bonds) + fplus(n, fields) + fplus(n, fields).H
    return all(sp.expand(e) == 0 for e in (H * W - W * H)) and all(sp.expand(e) == 0 for e in (zd(n) * W + W * zd(n)))


# ---- part 1: the fence ---------------------------------------------------------------------------------------------------
say("Part 1. The fence: vectors with d excited at the vacuum's energy, sector by sector, exactly")
rng = random.Random(20261118)
for n, count in ((4, 30), (5, 20), (6, 6)):
    seen, nb, readings, comm_ok = set(), 0, {}, True
    S = splus(n)
    while nb < count:
        b = draw_bonds(n, rng)
        if key(n, b) in seen:
            continue
        seen.add(key(n, b))
        H = hb(n, b)
        comm_ok &= (S * H - H * S).is_zero_matrix
        bl = blind_dim(n, b)
        if bl:
            continue
        nb += 1
        f = fence(n, b, H)
        dims = tuple(len(f[k]) for k in range(1, n))
        readings[dims] = readings.get(dims, 0) + 1
    # a constructed control: the triangle (1, 2, 3) on the reciprocal locus e2 = 0 (3, 6, -2), every other site a leaf
    # on site 1; the triangle's second kernel vector, shifted to vanish at 1, extends by zeros and vanishes at d
    b = {(1, 2): 3, (1, 3): 6, (2, 3): -2, (0, 1): 1}
    b.update({(1, j): 1 for j in range(4, n)})
    f = fence(n, b)
    check(len(laplacian(n, b).nullspace()) >= 2 and len(f[n - 1]) >= 1,
          f"N = {n}: control, the constructed blind row {sorted(b.items())}: ker L_J holds a vector vanishing at d, "
          f"and the fence at N - 1 is not empty ({len(f[n - 1])})")
    check(comm_ok, f"N = {n}: [S+, H_b] = 0 exactly on every drawn row")
    zero_at = [1, 2, n - 1]
    check(all(all(d[k - 1] == 0 for k in zero_at) for d in readings),
          f"N = {n}: {nb} distinct rows where d is not blind, fence empty at n = {sorted(set(zero_at))} on every one")
    say(f"  read: fence dimensions by sector n = 1..{n - 1} -> rows: {readings}")

B5 = {(0, 1): 1, (0, 2): 1, (0, 3): 2, (0, 4): 2, (1, 3): 1, (2, 4): 2, (2, 3): 3, (1, 2): -1, (3, 4): -2}
check(blind_dim(5, B5) == 0, "the constructed N = 5 row 0-1:1, 0-2:1, 0-3:2, 0-4:2, 1-3:1, 2-4:2, 2-3:3, 1-2:-1, 3-4:-2: "
                              "d not blind")
f5 = fence(5, B5)
check(tuple(len(f5[k]) for k in range(1, 5)) == (0, 0, 1, 0), "its fence by sector is (0, 0, 1, 0), exactly")
# X^5 (|0>_0 (x) s_12 (x) s_34), s = |01> - |10> on (a, b) meaning a unexcited b excited minus the reverse
y = {}
for (a1, b1, s1), (a2, b2, s2) in itertools.product(((1, 2, 1), (2, 1, -1)), ((3, 4, 1), (4, 3, -1))):
    x = (1 << b1) | (1 << b2)                     # the excited sites of the two singlet terms, d unexcited
    y[x ^ 31] = s1 * s2
idx5, v5 = f5[3][0]
vec = {x: v5[i] for i, x in enumerate(idx5) if v5[i] != 0}
ratio = next(iter(vec.values())) / y[next(iter(vec))]
check(set(vec) == set(y) and all(vec[x] == ratio * y[x] for x in y),
      "its fence vector is X^5 applied to |0>_0 (x) s_12 (x) s_34")
say(f"  read: E0 = {sum(B5.values())}, -3 J_12 - 3 J_34 = {-3 * B5[(1, 2)] - 3 * B5[(3, 4)]}")

# ---- part 2: the hand steps at N = 4, symbolically ----------------------------------------------------------------------
say()
say("Part 2. The hand steps of the case analysis at N = 4, symbolically (d = 0; u, v, w = 1, 2, 3)")
wu, wv, ww = sp.symbols('wu wv ww')
cu, cv, cw = sp.conjugate(wu), sp.conjugate(wv), sp.conjugate(ww)
# F- = sum conj(w_j) sigma+_j (lowering) from the bit action, cut to two excitations of {1, 2, 3} -> one, d unexcited
cw_of = {1: cu, 2: cv, 3: cw}
pair_states = [(1 << 1) | (1 << 2), (1 << 1) | (1 << 3), (1 << 2) | (1 << 3)]       # y_12, y_13, y_23
one_states = [1 << 1, 1 << 2, 1 << 3]                                              # |1>, |2>, |3>
Fm = sp.zeros(3, 3)
for c, x in enumerate(pair_states):
    for j in (1, 2, 3):
        if (x >> j) & 1:
            Fm[one_states.index(x ^ (1 << j)), c] += cw_of[j]
check(sp.simplify(Fm.det() + 2 * cu * cv * cw) == 0, "det of F- on two excitations of u, v, w is -2 conj(w_u w_v w_w)")
Fm0 = Fm.subs(ww, 0)
ns = Fm0.nullspace()
check(len(ns) == 1 and sp.simplify(ns[0][0]) == 0 and sp.simplify(ns[0][1] * (-cu) - ns[0][2] * cv) == 0,
      "with w_w = 0 its kernel is spanned by y0 = conj(w_v)|uw> - conj(w_u)|vw>")
Js = {e: sp.Symbol(f'J{e[0]}{e[1]}') for e in itertools.combinations(range(4), 2)}
H4 = hb(4, Js)
pairs = [(1, 2), (1, 3), (2, 3)]
ok = True
for j in (1, 2, 3):
    col = 1 | (1 << j)                               # |1>_0 (x) e_j
    for (a, b) in pairs:
        row = (1 << a) | (1 << b)                    # |0>_0 (x) |ab>
        pred = 2 * (Js[(0, b)] if a == j else 0) + 2 * (Js[(0, a)] if b == j else 0)
        ok &= sp.simplify(H4[row, col] - pred) == 0
check(ok, "the hop of d's bonds: <0_d, jk| H_b |1_d, e_i> = 2 J_dk [i = j] + 2 J_dj [i = k], from the full bond matrix")
y0 = {(1 << 1) | (1 << 3): cv, (1 << 2) | (1 << 3): -cu}   # d unexcited
comp = sum(H4[(1 << 1) | (1 << 2), x] * c for x, c in y0.items())
check(sp.simplify(comp - 2 * (Js[(2, 3)] * cv - Js[(1, 3)] * cu)) == 0,
      "the |uv> entry of H_b on |0>_d (x) y0 is 2 (J_vw conj(w_v) - J_uw conj(w_u))")
Sp4 = splus(4)
exc1 = [1 | (1 << j) for j in (1, 2, 3)]                                           # |1>_0 (x) e_j
exc2 = [1 | x for x in pair_states]                                                # |1>_0 (x) |jk>
Mp = sp.Matrix([[Sp4[r, c] for c in exc1] for r in exc2])
leak = [x for x in range(16) if not (x & 1)]
check(Mp.det() != 0 and all(Sp4[r, c] == 0 for r in leak for c in exc1),
      f"step (i): S+ keeps |1>_d (x) u' at d excited and acts on u' as a 3 x 3 matrix of determinant {Mp.det()}, "
      "so S+ (|1>_d (x) u') = 0 only at u' = 0")
Sm4 = Sp4.T
Mm = sp.Matrix([[Sm4[r, c] for c in pair_states] for r in one_states])
leak = [x for x in range(16) if x & 1]
check(Mm.det() != 0 and all(Sm4[r, c] == 0 for r in leak for c in pair_states),
      f"step (i): S- keeps |0>_d (x) y at d unexcited and acts on y as a 3 x 3 matrix of determinant {Mm.det()}, "
      "so S- (|0>_d (x) y) = 0 only at y = 0")

# ---- part 3: condition (I) and the far counts at N = 4 --------------------------------------------------------------------
say()
say("Part 3. Condition (I) and the far elements at N = 4, at seats that are not blind")


def draw_fields(n, rng, zero=None, line=False):
    if line:
        m = rng.choice([(1, 0), (3, 4), (1, 1)])
        f = {}
        for j in range(1, n):
            k = rng.choice([-2, -1, 1, 2])
            f[j] = (k * m[0], k * m[1])
        return f
    f = {j: rng.choice(FLD) for j in range(1, n)}
    if zero is not None:
        f[zero] = (0, 0)
    return f


def commutant_dim_modp(n, bonds, k):
    """An upper bound on the rational commutant dimension (a kernel modulo p is at least the rational one)."""
    H = hb(n, bonds)
    idx = sector(n, k)
    A = np.array([[int(H[a, b]) for b in idx] for a in idx], dtype=np.int64)
    Z = np.diag([(-1 if x & 1 else 1) for x in idx]).astype(np.int64)
    I = np.eye(len(idx), dtype=np.int64)
    rows = np.vstack([np.kron(I, A) - np.kron(A.T, I), np.kron(I, Z) - np.kron(Z.T, I)])
    return len(idx) ** 2 - rank_mod(rows, P)


def commutant_dim_exact(n, bonds, k):
    H = hb(n, bonds)
    idx = sector(n, k)
    A = sp.Matrix([[H[a, b] for b in idx] for a in idx])
    Z = sp.diag(*[(-1 if x & 1 else 1) for x in idx])
    m = len(idx)
    I = sp.eye(m)
    rows = sp.Matrix.vstack(sp.kronecker_product(I, A) - sp.kronecker_product(A.T, I),
                            sp.kronecker_product(I, Z) - sp.kronecker_product(Z.T, I))
    from sympy.polys.matrices import DomainMatrix
    return m * m - DomainMatrix.from_Matrix(rows).convert_to(sp.QQ).rank()   # exact rank over Q, not Matrix.rank


rng = random.Random(4118)
cells = {'off one line, all three fields nonzero': dict(), 'off one line, one field zero': dict(zero=True),
         'on one line, all three nonzero (class A)': dict(line=True)}
for name, opt in cells.items():
    seen, rows, okI, okfar = set(), 0, True, True
    target = 60
    while rows < target:
        b = draw_bonds(4, rng)
        f = draw_fields(4, rng, zero=(rng.choice([1, 2, 3]) if opt.get('zero') else None), line=opt.get('line', False))
        if key(4, b, f) in seen or blind_dim(4, b):
            continue
        if opt.get('line'):
            if not on_one_line(f):
                continue
        elif on_one_line(f):
            continue
        seen.add(key(4, b, f))
        rows += 1
        okI &= cyclic_dim(4, b, f) == 16
        fe, fo = far_count(4, b, f, 0), far_count(4, b, f, 1)
        if opt.get('line'):
            okfar &= (fe, fo) == (1, 0) and colouring_exact(4, b, f)
        else:
            okfar &= (fe, fo) == (0, 0)
    check(okI, f"{name}: {rows} distinct rows, the cyclic space has dimension 16 at every one (I)")
    check(okfar, f"{name}: far elements homogeneous under Z^4 (even, odd) = "
                 + ("(1, 0) at every row, the 1 being class A's colouring, exact" if opt.get('line') else "(0, 0) at every row"))

# rows whose two-excitation sector is reducible while d is not blind
rng = random.Random(77)
red, tried = [], 0
while tried < 2500 and len(red) < 6:
    b = draw_bonds(4, rng)
    tried += 1
    if blind_dim(4, b):
        continue
    if commutant_dim_modp(4, b, 2) > 1 and commutant_dim_exact(4, b, 2) > 1:
        red.append(b)
okI, okfar, nsets = True, True, 0
for b in red:
    for z in (None, 1, 2, 3):
        for _ in range(3):
            f = draw_fields(4, rng, zero=z)
            if on_one_line(f):
                continue
            nsets += 1
            okI &= cyclic_dim(4, b, f) == 16
            okfar &= (far_count(4, b, f, 0), far_count(4, b, f, 1)) == (0, 0)
check(len(red) >= 1, f"rows with d not blind and the two-excitation sector reducible (commutant > 1, exact) exist "
                     f"in the draw")
say(f"  read: {len(red)} such rows in {tried} draws")
check(okI and okfar, f"at them, {nsets} field sets off one line (some with a field zero): the cyclic space has "
                     f"dimension 16 and no far element homogeneous under Z^4")

# control: a blind row where (I) fails
Bc = {(0, 1): 1, (0, 2): 2, (1, 2): -2, (1, 3): 4, (2, 3): 2}
Fc = {1: (2, 0), 2: (1, 0), 3: (0, 1)}
v = sp.zeros(16, 1)
v[1 | 2] = Bc[(0, 1)]; v[1 | 4] = -Bc[(0, 2)]
Hc = hb(4, Bc)
lam = (Hc * v)[1 | 2] / v[1 | 2]
Fm_c = fplus(4, Fc).H
check(blind_dim(4, Bc) == 1 and (Hc * v - lam * v).is_zero_matrix and (Fm_c * v).is_zero_matrix and
      (zd(4) * v + v).is_zero_matrix,
      "control, the blind row 0-1:1, 0-2:2, 1-2:-2, 1-3:4, 2-3:2 with fields (2,0), (1,0), (0,1): "
      "|1>_0 (x) (J_01 e_1 - J_02 e_2) is an eigenvector of H_b with d excited killed by F-, exactly")
cd = cyclic_dim(4, Bc, Fc)
check(cd == 15, f"so the cyclic space there has dimension at most 15, and modulo p it reads exactly 15 ({cd})")

# ---- part 4: measured at N = 5 --------------------------------------------------------------------------------------------
say()
say("Part 4. At N = 5: the conditions at drawn rows, and the constructed row where (II) fails")
rng = random.Random(5118)
for name, line in (('off one line', False), ('on one line, all four nonzero (class A)', True)):
    seen, rows, okI, okII, okfar = set(), 0, True, True, True
    while rows < 25:
        b = draw_bonds(5, rng)
        f = draw_fields(5, rng, line=line)
        if key(5, b, f) in seen or blind_dim(5, b) or on_one_line(f) != line:
            continue
        seen.add(key(5, b, f))
        rows += 1
        okI &= cyclic_dim(5, b, f) == 32
        fz = fence(5, b)
        okII &= all(len(fz[k]) == 0 for k in fz)
        fo, fe = far_count(5, b, f, 1), far_count(5, b, f, 0)
        okfar &= (fo, fe) == ((1, 0) if line else (0, 0)) and (not line or colouring_exact(5, b, f))
    check(okI and okII, f"{name}: {rows} distinct rows, (I) and (II) hold at every one")
    check(okfar, f"{name}: far elements homogeneous under Z^5 (odd, even) = " +
          ("(1, 0), the 1 being class A's colouring, exact" if line else "(0, 0)"))
for f in ({1: (1, 0), 2: (-1, 0), 3: (0, 1), 4: (0, -1)}, {1: (3, 4), 2: (-3, -4), 3: (1, 0), 4: (-1, 0)},
          {1: (5, 1), 2: (4, -3), 3: (1, 1), 4: (0, 2)}, {1: (1, 0), 2: (2, 0), 3: (3, 4), 4: (1, 0)}):
    cdim = cyclic_dim(5, B5, f)
    fo, fe = far_count(5, B5, f, 1), far_count(5, B5, f, 0)
    check(cdim == 32 and (fo, fe) == (0, 0),
          f"the constructed row, fields {f}: (I) holds, (II) fails at n = 3, and no far element homogeneous under Z^5")

# the vacuum equation alone at the antiparallel fields: alpha F+ v - tau Xi F+|vac> + chi F-|top> = 0
Hd5 = hb(5, B5)
vfen = sp.zeros(32, 1)
for i, x in enumerate(idx5):
    vfen[x] = v5[i]
for f in ({1: (1, 0), 2: (-1, 0), 3: (0, 1), 4: (0, -1)}, {1: (3, 4), 2: (-3, -4), 3: (1, 0), 4: (-1, 0)},
          {1: (5, 1), 2: (4, -3), 3: (1, 1), 4: (0, 2)}):
    Fp = fplus(5, f)
    vac = sp.zeros(32, 1); vac[0] = 1
    top = sp.zeros(32, 1); top[31] = 1
    Xi = sp.zeros(32, 32)
    for x in range(32):
        Xi[x ^ 31, x] = 1
    M = sp.Matrix.hstack(Fp * vfen, -Xi * Fp * vac, Fp.H * top)
    ker = M.nullspace()
    anti = f[2] == (-f[1][0], -f[1][1])
    good = len(ker) == 1 and ker[0][0] != 0 if anti else len(ker) == 0
    check(good, f"the constructed row, fields {f}: the vacuum equation alone in (alpha, tau, chi) has "
                + ("a one-dimensional solution with alpha != 0 (fields antiparallel within the pairs)" if anti
                   else "only the zero solution"))

# ---- part 5: N = 4 on one line with a field zero -----------------------------------------------------------------------
say()
say("Part 5. N = 4 on one line with a field zero: steps (v) to (vii) and the antiparallel locus")
st = lambda *s: sum(1 << i for i in s)
J01, J02, J03, J12, J13, J23, r1, r2 = sp.symbols('J01 J02 J03 J12 J13 J23 r1 r2')
u1, u2, u3 = sp.symbols('u1 u2 u3')
Hs = hb(4, {(0, 1): J01, (0, 2): J02, (0, 3): J03, (1, 2): J12, (1, 3): J13, (2, 3): J23})
e23 = sp.zeros(16, 1); e23[st(2, 3)] = 1
check(sp.expand((Hs * e23)[st(1, 3)] - 2 * J12) == 0 and sp.expand((Hs * e23)[st(1, 2)] - 2 * J13) == 0,
      "step (v): the bonds carry |0>_0 (x) |23> to |13> and |12> with 2 J12 and 2 J13")
uu = sp.zeros(16, 1); uu[st(0, 1)] = u1; uu[st(0, 2)] = u2; uu[st(0, 3)] = u3
Hu = Hs * uu
check(all(sp.expand(Hu[st(*pq)] - val) == 0 for pq, val in (((1, 2), 2 * (J02 * u1 + J01 * u2)),
                                                            ((1, 3), 2 * (J03 * u1 + J01 * u3)),
                                                            ((2, 3), 2 * (J03 * u2 + J02 * u3)))),
      "the hop of 0's bonds on |1>_0 (x) u: entries 2(J02 u1 + J01 u2), 2(J03 u1 + J01 u3), 2(J03 u2 + J02 u3)")
alpha = sp.expand((r1 * Hu[st(1, 3)] + r2 * Hu[st(2, 3)]).subs(u1, -r2 * u2 / r1))
check(sp.expand(alpha - 2 * u3 * (r1 * J01 + r2 * J02)) == 0,
      "step (alpha): with r1 u1 + r2 u2 = 0 the ratio r2 : -r1 of the |13> and |23> entries asks (r1 J01 + r2 J02) u3 = 0")
r0 = sp.zeros(16, 1); r0[st(1, 3)] = r2; r0[st(2, 3)] = -r1
Hr = Hs * r0
check([sp.expand(Hr[st(0, k)]) for k in (1, 2, 3)] == [sp.expand(2 * r2 * J03), sp.expand(-2 * r1 * J03),
                                                       sp.expand(2 * (r2 * J01 - r1 * J02))],
      "step (beta): the hop carries |0>_0 (x) (r2|13> - r1|23>) to |1>_0 (x) 2(r2 J03, -r1 J03, r2 J01 - r1 J02)")
gam2 = sp.expand(Hr[st(1, 3)] * (-r1) - Hr[st(2, 3)] * r2)
check(sp.expand(Hr[st(1, 2)] - 2 * (r2 * J23 - r1 * J13)) == 0 and
      sp.expand(gam2 + 2 * (r1 * r2 * (J13 - J23 - J01 + J02) - J12 * (r1 ** 2 - r2 ** 2))) == 0,
      "step (gamma): the |12> entry 2(r2 J23 - r1 J13), and the ratio of the |13> and |23> entries "
      "J13 - J23 - J01 + J02 = J12 (r1^2 - r2^2)/(r1 r2)")
e1v = sp.zeros(16, 1); e1v[st(0, 1)] = 1
e2v = sp.zeros(16, 1); e2v[st(0, 2)] = 1
check(sp.expand((Hs * e1v)[st(0, 3)] - 2 * J13) == 0 and sp.expand((Hs * e2v)[st(0, 3)] - 2 * J23) == 0,
      "step (delta): within the excited 0 the bonds carry e1 and e2 to e3 with 2 J13 and 2 J23")

# the exhaustive grid: (I) fails exactly on the antiparallel locus (real fields, so primes need no square root of -1)
GP = (1000003, 1000033)


def cyclic_dim_real4(H, F):
    best = 0
    Zm = np.diag([(-1 if x & 1 else 1) for x in range(16)]).astype(np.int64)
    for p in GP:
        gens = [H % p, Zm % p, F % p]
        basis, piv, queue = [], [], [np.eye(16, dtype=np.int64)[0]]
        while queue:
            v = queue.pop() % p
            r = v.copy()
            for bv, c in zip(basis, piv):
                if r[c]:
                    r = (r - r[c] * bv) % p
            nz = np.nonzero(r)[0]
            if len(nz) == 0:
                continue
            c = nz[0]
            r = (r * pow(int(r[c]), p - 2, p)) % p
            for t in range(len(basis)):
                if basis[t][c]:
                    basis[t] = (basis[t] - basis[t][c] * r) % p
            basis.append(r); piv.append(c)
            for G in gens:
                queue.append(G.dot(v) % p)
        best = max(best, len(basis))
    return best


def on_locus(b, fl, z):
    if len(fl) != 2:
        return False
    (i, wi), (j, wj) = sorted(fl.items())
    g = lambda x, y: b.get(tuple(sorted((x, y))), 0)
    return wi == (-wj[0], -wj[1]) and g(0, i) == g(z, i) == -g(0, j) == -g(z, j)


edges4 = list(itertools.combinations(range(4), 2))
FSETS = []
for z in (1, 2, 3):
    i, j = [s for s in (1, 2, 3) if s != z]
    for a, bb in ((1, -1), (1, 2), (2, -1), (1, 1), (1, 0)):
        fl = {s: (v, 0) for s, v in ((i, a), (j, bb)) if v}
        FSETS.append((z, fl, np.array(fplus(4, fl).tolist(), dtype=np.int64)))
counts, miss, locus_rows, locus_ok = {}, 0, 0, True
for Js in itertools.product((-1, 0, 1, 2), repeat=6):
    b = {e: J for e, J in zip(edges4, Js) if J}
    if not connected(4, b) or blind_dim(4, b):
        continue
    H = np.array(hb(4, b).tolist(), dtype=np.int64)
    for z, fl, F in FSETS:
        full = cyclic_dim_real4(H, F) == 16
        pred = on_locus(b, fl, z)
        counts[(len(fl), full, pred)] = counts.get((len(fl), full, pred), 0) + 1
        miss += full == pred
        if pred:
            locus_rows += 1
            locus_ok &= (far_count(4, b, fl, 0), far_count(4, b, fl, 1)) == (2, 0)
check(miss == 0 and counts.get((2, False, True), 0) > 0 and counts.get((1, True, False), 0) > 0,
      "grid exhaustive in the couplings, {-1, 0, 1, 2} on the six bonds, connected, 0 not blind, five pairs of fields "
      "on the x line with a field zero on each site in turn: the rank modulo p reaches 16 at every row off the antiparallel locus, which "
      "certifies (I) there, and falls short at every row on it, where Q (below, with the labels relabelled) bounds it")
check(locus_rows > 0 and locus_ok,
      f"at the grid's {locus_rows} locus rows the far elements (even, odd) = (2, 0) modulo p, at most 2 and 0, "
      "and at least class A's colouring and h (below)")
say(f"  read: (nonzero fields, cyclic = 16, on the locus) -> cases: {dict(sorted(counts.items()))}")

# the locus: Q, the element h and the two dimensions, exactly
c_, e_, f_, fx, fy = sp.symbols('c e f fx fy')
BL = {(0, 1): -c_, (0, 2): c_, (1, 3): -c_, (2, 3): c_, (0, 3): e_, (1, 2): f_}
HL = hb(4, BL)
FL = fplus(4, {1: (fx, fy), 2: (-fx, -fy)})
HLt = HL + FL + FL.H
q0 = sp.zeros(16, 1); q0[st(1, 3)] = 1; q0[st(2, 3)] = 1        # |0>_0 |1>_3 (x) (|10> + |01>) on 1, 2
q1 = sp.zeros(16, 1); q1[st(0, 1)] = 1; q1[st(0, 2)] = 1        # |1>_0 |0>_3 (x) (|10> + |01>) on 1, 2
check(all(sp.expand(x) == 0 for x in (HL * q0 - (f_ - e_) * q0 - 2 * e_ * q1)) and
      all(sp.expand(x) == 0 for x in (HL * q1 - (f_ - e_) * q1 - 2 * e_ * q0)) and
      all(sp.expand(x) == 0 for x in (FL.H * q0)) and all(sp.expand(x) == 0 for x in (FL.H * q1)),
      "on the locus, symbolically in c, J03 = e, J12 = f and the field (fx, fy): Q = t (x) span{|01>, |10>} on 0, 3 is "
      "kept by the bonds, which act as J12 - J03 + 2 J03 SWAP_03, and killed by F-")
hh = 2 * (q0 * q1.T + q1 * q0.T) / 2                            # (X0X3 + Y0Y3) (x) |t><t|, |t> = (q-part)/sqrt 2
check(all(sp.expand(x) == 0 for x in (HLt * hh - hh * HLt)) and all(x == 0 for x in (zd(4) * hh + hh * zd(4))),
      "h = (X0 X3 + Y0 Y3) (x) |t><t| commutes with the whole H and anticommutes with Z_0, symbolically")
BLd = dict(BL); BLd[(0, 1)] = -c_ + 1
HLd = hb(4, BLd) + FL + FL.H
check(any(sp.expand(x) != 0 for x in (HLd * hh - hh * HLd)),
      "control: with J01 moved off -c the same h no longer commutes with H")
for (cc, ee, ff), m in (((1, 2, 3), (1, 0)), ((2, -1, 3), (3, 4)), ((-1, 3, 0), (1, 1))):
    bl = {k: v.subs({c_: cc, e_: ee, f_: ff}) for k, v in BL.items()}
    bl = {k: int(v) for k, v in bl.items() if v != 0}
    fl = {1: m, 2: (-m[0], -m[1])}
    cd = cyclic_dim(4, bl, fl)
    fe, fo = far_count(4, bl, fl, 0), far_count(4, bl, fl, 1)
    hnum = hh.subs({c_: cc, e_: ee, f_: ff})
    ms = sp.Matrix([[0, m[0] - sp.I * m[1]], [m[0] + sp.I * m[1], 0]])
    A4 = sp.kronecker_product(ms, ms, ms, ms)
    indep = hnum.rank() == 2 and all(sp.expand(x) == 0 for x in (A4 * A4 - (m[0] ** 2 + m[1] ** 2) ** 4 * sp.eye(16)))   # h of rank 2, A invertible
    check(blind_dim(4, bl) == 0 and cd == 14 and (fe, fo) == (2, 0) and colouring_exact(4, bl, fl) and indep,
          f"locus row c, J03, J12 = {cc}, {ee}, {ff}, field {m} on 1 and its negative on 2: 0 not blind, cyclic "
          f"dimension 14 (at least 14 by the rank, at most 14 with Q orthogonal to it), far elements (even, odd) = "
          f"(2, 0), the 2 reached by class A's colouring and h, independent, exactly")
bl0 = {(0, 1): -1, (0, 2): 1, (1, 3): -1, (2, 3): 1, (1, 2): 3}
say(f"  read, the locus at J03 = 0: blind dimension {blind_dim(4, bl0)}, far elements (even, odd) = "
    f"({far_count(4, bl0, {1: (1, 0), 2: (-1, 0)}, 0)}, {far_count(4, bl0, {1: (1, 0), 2: (-1, 0)}, 1)}) modulo p")
for fl, what in (({1: (1, 0), 2: (1, 0)}, "parallel"), ({1: (1, 0), 2: (-2, 0)}, "opposite, unequal")):
    bl = {(0, 1): -1, (0, 2): 1, (1, 3): -1, (2, 3): 1, (0, 3): 2, (1, 2): 3}
    check(cyclic_dim(4, bl, fl) == 16 and (far_count(4, bl, fl, 0), far_count(4, bl, fl, 1)) == (1, 0),
          f"control, the locus's bonds with the fields {what}: cyclic dimension 16 and far elements (even, odd) = (1, 0)")

say()
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} FAILED: " + "; ".join(FAILS))
(Path(__file__).parent / "results" / "blind_seat_descent_gate.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
sys.exit(1 if FAILS else 0)
