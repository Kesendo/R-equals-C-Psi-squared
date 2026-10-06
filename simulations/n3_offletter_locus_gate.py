"""Gate for the N = 3 triangle with the jump Z on one site d alone: fields off the letters, and XX + YY bonds.

PROOF_PALINDROME_COMPLEMENT_CONNECTION decides the triangle's palindromic locus for Heisenberg bonds and letter fields
(classes A, B, C). This gate reads its neighbours. A row pairs exactly when near = far >= 1 (F158: near and far are the
kernels of [H, .] on the strings whose d-letter is in {I, Z} and in {X, Y}); the sweeps rank them modulo two primes that
must agree, the families are checked as exact elements. Classes for fields as vectors: A every nonzero field along one
common axis perpendicular to z ((n.sigma)^3); B and C as the proposition states them; D below.

G1  (exact, symbolic) Family D, d at an end of the path d - v - u (J_du = 0), f_d = p X_d, f_u = q Y_u: under
    Heisenberg bonds with f_v = -(f_d + f_u), and under XX + YY bonds with f_v = 0, the element
    W = J_dv q YYY + p q XYI + J_uv p XXX (sites d, u, v) commutes with H, anticommutes with Z_d and squares to
    (J_dv^2 q^2 + J_uv^2 p^2 + p^2 q^2) I. So the proposition read for vector fields is false: a Heisenberg row of D
    with p, q != 0 lies outside A, B and C. Controls (this W fails): f_v changed, a Z part on f_d, an X part on f_u, J_du != 0.
G2  (exact, symbolic) XX + YY bonds on the path u - d - v (J_uv = 0), f_d = 0, f_u = a X_u, f_v = b Y_v:
    W = J_du X_d X_u + J_dv Y_d Y_v commutes with H, anticommutes with Z_d and squares to (J_du^2 + J_dv^2) I. At equal
    weights this is the census row of THE_PALINDROME_AS_A_COLOURING (element XXI + IYY); the weights are free here.
    Controls: a ZZ part on a bond, J_uv != 0, a field on d.
G3  (measured) On a seeded sweep of distinct Heisenberg rows with integer vector fields, every palindromic row lies in
    A, B, C or D, and every row of A, B, C pairs; counted separately, the rows no turn about z makes letter rows. With C
    left out of the rule palindromic rows are reported outside.
G4  (measured) XX + YY bonds: three rows of C outside A break, and on the path u - d - v with f_d = 0, fields of equal
    magnitude perpendicular to z pair at six angles and three couplings.
G5  (exact, one row) Locus E, rows J = (J_du, J_dv, J_uv) with all couplings nonzero and f_d = 0: the row J = (1, 2, 6),
    f_u = 8X, f_v = -5X + 3Y pairs through an eighteen-string W with W^2 = 4500 I; two perturbations break it; the row
    satisfies the three E relations, and its dephased seat is not blind.
G6  (measured) Points on the three E relations, r and f_v solved from J_du, J_dv and f_u (f_v's Y part a square
    root, taken modulo primes where it exists): each pairs with near = far = 1 at two primes, and each breaks when f_v
    moves by X/10; points that keep the second and third relations at J_uv + 1, off the first, break.
Output: simulations/results/n3_offletter_locus_gate.txt. Run: python simulations/n3_offletter_locus_gate.py"""
import itertools, os, random, sys
from fractions import Fraction as F
import numpy as np
import sympy as sp

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass
OUT = []
_print = print
def print(*a, **k):
    _print(*a, **k); OUT.append(" ".join(str(x) for x in a))
FAILS = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

L = "IXYZ"
MUL = {}
for a in L:
    for b in L:
        if a == "I": MUL[(a, b)] = (1, b)
        elif b == "I": MUL[(a, b)] = (1, a)
        elif a == b: MUL[(a, b)] = (1, "I")
for a, b, c in (("X", "Y", "Z"), ("Y", "Z", "X"), ("Z", "X", "Y")):
    MUL[(a, b)] = (sp.I, c); MUL[(b, a)] = (-sp.I, c)
def smul(s, t):
    ph = sp.Integer(1); out = []
    for x, y in zip(s, t):
        p, c = MUL[(x, y)]; ph *= p; out.append(c)
    return ph, "".join(out)
def pmul(A, B):
    out = {}
    for s, a in A.items():
        for t, b in B.items():
            ph, w = smul(s, t); out[w] = sp.expand(out.get(w, 0) + ph * a * b)
    return {w: c for w, c in out.items() if c != 0}
def padd(A, B, sign=1):
    out = dict(A)
    for w, c in B.items(): out[w] = sp.expand(out.get(w, 0) + sign * c)
    return {w: c for w, c in out.items() if c != 0}
def comm(A, B): return padd(pmul(A, B), pmul(B, A), -1)
def acomm(A, B): return padd(pmul(A, B), pmul(B, A), 1)
def term(letters_by_site, c):
    s = ["I"] * 3
    for k, v in letters_by_site.items(): s[k] = v
    return {"".join(s): c}
WORDS = ["".join(w) for w in itertools.product(L, repeat=3)]
IDX = {w: k for k, w in enumerate(WORDS)}
D = 0  # the dephased site; sites 1 = u, 2 = v

def hamiltonian(J, f, bond_letters):
    H = {}
    for (a, b), j in {(0, 1): J[0], (0, 2): J[1], (1, 2): J[2]}.items():
        if j:
            for P in bond_letters: H = padd(H, term({a: P, b: P}, sp.Integer(j)))
    for l in range(3):
        for k, P in enumerate("XYZ"):
            if f[l][k]: H = padd(H, term({l: P}, sp.Integer(f[l][k])))
    return H

# ---------------------------------------------------------------- ends modulo primes
def sqrtm1(p):
    for g in range(2, 400):
        x = pow(g, (p - 1) // 4, p)
        if x * x % p == p - 1: return x
PRIMES = [998244353, 469762049]
NEAR = [w for w in WORDS if w[D] in "IZ"]
FAR = [w for w in WORDS if w[D] in "XY"]
def nullity_mod(H, cols, p):
    ip = sqrtm1(p); M = np.zeros((64, len(cols)), dtype=np.int64)
    for j, P in enumerate(cols):
        for h, c in H.items():
            ph, w = smul(h, P); ph2, _ = smul(P, h)
            if ph == ph2: continue
            z = sp.expand(2 * ph * c); re, im = int(sp.re(z)), int(sp.im(z))
            M[IDX[w], j] = (M[IDX[w], j] + re + im * ip) % p
    r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, 64) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]; M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in range(64):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
    return len(cols) - r
def verdict(H):
    reads = [(nullity_mod(H, NEAR, p), nullity_mod(H, FAR, p)) for p in PRIMES]
    n, fa = reads[0]
    return reads[0] == reads[1], n == fa and fa >= 1

# ---------------------------------------------------------------- classes, fields as vectors
def dot(a, b): return sum(F(x) * F(y) for x, y in zip(a, b))
def zero(v): return all(x == 0 for x in v)
def class_A(f):
    nz = [v for v in f if not zero(v)]
    if not nz: return True
    if any(v[2] for v in nz): return False
    return all(nz[0][0] * v[1] - nz[0][1] * v[0] == 0 for v in nz)
def rot_pi(m, v):
    k = 2 * dot(m, v) / dot(m, m); return tuple(k * F(mi) - F(vi) for mi, vi in zip(m, v))
def class_B(J, f):
    if J[0] != J[1]: return False
    f0, f1, f2 = f
    cands = [(f0[0], f0[1], 0)] if (f0[0] or f0[1]) else []
    if (f1[0] + f2[0], f1[1] + f2[1]) != (0, 0): cands.append((f1[0] + f2[0], f1[1] + f2[1], 0))
    cands += [(-v[1], v[0], 0) for v in (f1, f2) if v[0] or v[1]] + [(1, 0, 0), (0, 1, 0)]
    return any(rot_pi(m, f0) == tuple(map(F, f0)) and rot_pi(m, f1) == tuple(map(F, f2)) for m in cands)
def class_C(J, f):
    e2 = J[0] * J[1] + J[0] * J[2] + J[1] * J[2]
    return e2 == 0 and zero(f[0]) and f[1][2] == 0 and f[2][2] == 0 and dot(f[1], f[1]) == dot(f[2], f[2])
def class_D(J, f, heis):
    """d at an end: one of J_du, J_dv zero; f_d and the far field nonzero, perpendicular, perpendicular to z;
    the middle field -(f_d + f_far) under Heisenberg bonds, 0 under XX + YY."""
    if (J[0] == 0) == (J[1] == 0) or J[2] == 0: return False
    far, mid = (1, 2) if J[0] == 0 else (2, 1)
    fd, ff, fm = f[0], f[far], f[mid]
    if zero(fd) or zero(ff) or fd[2] or ff[2] or dot(fd, ff) != 0: return False
    want = tuple(-(x + y) for x, y in zip(fd, ff)) if heis else (0, 0, 0)
    return tuple(fm) == want
def letter_turnable(f):
    """some turn about z makes every field a letter field"""
    xy = []
    for v in f:
        if zero(v): continue
        if v[2] and (v[0] or v[1]): return False
        if not v[2]: xy.append((v[0], v[1]))
    return all(a[0] * b[1] - a[1] * b[0] == 0 or a[0] * b[0] + a[1] * b[1] == 0 for a, b in itertools.combinations(xy, 2))

# ---------------------------------------------------------------- G1
j1, j2, p_, q_, a, b, Jdu, Jdv, Juv, g, e = sp.symbols("J_dv J_uv p q a b J_du J_dv J_uv g e", real=True)  # G1: j1, j2; G2: Jdu, Jdv, Juv
def sym_bond(i, j, c, letters):
    H = {}
    for P in letters: H = padd(H, term({i: P, j: P}, c))
    return H
Zd = term({0: "Z"}, sp.Integer(1))
WD = padd(padd(term({0: "Y", 1: "Y", 2: "Y"}, j1 * q_), term({0: "X", 1: "Y"}, p_ * q_)), term({0: "X", 1: "X", 2: "X"}, j2 * p_))
def H_D(letters, mid):
    H = padd(sym_bond(0, 2, j1, letters), sym_bond(1, 2, j2, letters))
    H = padd(H, padd(term({0: "X"}, p_), term({1: "Y"}, q_)))
    return padd(H, mid)
mid_heis = padd(term({2: "X"}, -p_), term({2: "Y"}, -q_))
sqD = {"III": sp.expand(j1 ** 2 * q_ ** 2 + j2 ** 2 * p_ ** 2 + p_ ** 2 * q_ ** 2)}
for letters, mid, label in (("XYZ", mid_heis, "Heisenberg, f_v = -(f_d + f_u)"), ("XY", {}, "XX + YY, f_v = 0")):
    H = H_D(letters, mid)
    check(f"G1 family D, {label}: [H, W] = 0, {{Z_d, W}} = 0, W^2 = (J_dv^2 q^2 + J_uv^2 p^2 + p^2 q^2) I",
          comm(H, WD) == {} and acomm(Zd, WD) == {} and pmul(WD, WD) == sqD)
Hh = H_D("XYZ", mid_heis)
check("G1 control: f_v changed by e X_v breaks [H, W] = 0", comm(padd(Hh, term({2: "X"}, e)), WD) != {})
check("G1 control: a Z part on f_d breaks [H, W] = 0", comm(padd(Hh, term({0: "Z"}, e)), WD) != {})
check("G1 control: an X part on f_u (no longer perpendicular to f_d) breaks [H, W] = 0", comm(padd(Hh, term({1: "X"}, e)), WD) != {})
check("G1 control: J_du != 0 (a Heisenberg bond d - u) breaks [H, W] = 0", comm(padd(Hh, sym_bond(0, 1, e, "XYZ")), WD) != {})
rowD = ((0, 1, 3), ((2, 0, 0), (0, 5, 0), (-2, -5, 0)))
agree, pal = verdict(hamiltonian(*rowD, "XYZ"))
check("G1 the Heisenberg row J = (0, 1, 3), f_d = 2X, f_u = 5Y, f_v = -2X - 5Y pairs by the ends, lies outside A, B, C, "
      "in D, and no turn about z makes it a letter row",
      agree and pal and not (class_A(rowD[1]) or class_B(*rowD) or class_C(*rowD)) and class_D(*rowD, True)
      and not letter_turnable(rowD[1]))

# ---------------------------------------------------------------- G2
def xy_bond(i, j, c): return sym_bond(i, j, c, "XY")
H2 = padd(padd(xy_bond(0, 1, Jdu), xy_bond(0, 2, Jdv)), padd(term({1: "X"}, a), term({2: "Y"}, b)))
W2 = padd(term({0: "X", 1: "X"}, Jdu), term({0: "Y", 2: "Y"}, Jdv))
check("G2 [H, W] = 0, {Z_d, W} = 0, W^2 = (J_du^2 + J_dv^2) I for every J_du, J_dv, a, b",
      comm(H2, W2) == {} and acomm(Zd, W2) == {} and pmul(W2, W2) == {"III": sp.expand(Jdu ** 2 + Jdv ** 2)})
check("G2 control: a ZZ part on the bond d - u breaks [H, W] = 0", comm(padd(H2, term({0: "Z", 1: "Z"}, g)), W2) != {})
check("G2 control: J_uv (XX + YY) != 0 breaks [H, W] = 0", comm(padd(H2, xy_bond(1, 2, Juv)), W2) != {})
check("G2 control: a field on d breaks [H, W] = 0 (X, Y and Z each)", all(comm(padd(H2, term({0: P}, e)), W2) != {} for P in "XYZ"))
row2 = ((2, 3, 0), ((0, 0, 0), (1, 0, 0), (0, 5, 0)))
agree, pal = verdict(hamiltonian(*row2, "XY"))
check("G2 the row J = (2, 3, 0), f_u = X, f_v = 5Y pairs by the ends and lies outside A, B and C",
      agree and pal and not (class_A(row2[1]) or class_B(*row2) or class_C(*row2)))

# ---------------------------------------------------------------- G3
COUPLINGS = [(1, 1, 1), (1, 2, 3), (2, 2, 1), (2, 2, -3), (3, 6, -2), (1, -2, -2), (-1, 2, 2), (2, -1, 2), (1, 1, 0), (1, 0, 2),
             (2, 3, 0), (0, 1, 3), (0, 2, -1)]
VECS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 2, 0), (3, 4, 0), (4, 3, 0), (5, 0, 0), (0, 5, 0),
        (-3, 4, 0), (1, 0, 1), (1, 1, 1), (0, 1, -1), (3, 0, 4), (2, 1, 0), (-1, 0, 0), (-2, -5, 0), (-1, -1, 0)]
def g3_rows(seed, n):
    rnd = random.Random(seed); rows = set()
    while len(rows) < n:
        J = rnd.choice(COUPLINGS); mode = rnd.random()
        if mode < 0.2:
            m = rnd.choice([(1, 1, 0), (1, 2, 0), (3, 4, 0), (2, -1, 0)])
            f = tuple(tuple(k * x for x in m) for k in (rnd.choice([0, 1, -2, 3]) for _ in range(3)))
        elif mode < 0.4:
            J = (J[0], J[0], J[2]) if J[0] else (1, 1, J[2])
            m = rnd.choice([(1, 1, 0), (1, 2, 0), (2, 1, 0), (3, 4, 0)])
            f1 = rnd.choice(VECS); f2 = rot_pi(m, f1)
            if any(x.denominator != 1 for x in f2): continue
            f = (tuple(rnd.choice([0, 1, -1]) * x for x in m), f1, tuple(int(x) for x in f2))
        elif mode < 0.6:
            J = rnd.choice([(3, 6, -2), (1, -2, -2), (-1, 2, 2), (2, -1, 2)])
            fa, fb = rnd.choice([((3, 4, 0), (5, 0, 0)), ((3, 4, 0), (0, 5, 0)), ((1, 1, 0), (1, -1, 0)), ((4, 3, 0), (-3, 4, 0)),
                                 ((1, 2, 0), (2, -1, 0))])
            f = ((0, 0, 0), fa, fb)
        elif mode < 0.7:
            J = rnd.choice([(0, 1, 3), (0, 2, -1), (0, 1, 1)])
            pd, qu = rnd.choice([((2, 0, 0), (0, 5, 0)), ((1, 0, 0), (0, 1, 0)), ((3, 4, 0), (-4, 3, 0))])
            f = (pd, qu, tuple(-(x + y) for x, y in zip(pd, qu))) if rnd.random() < 0.7 else (pd, qu, rnd.choice(VECS))
        else:
            f = tuple(rnd.choice(VECS) for _ in range(3))
            if rnd.random() < 0.5: f = (f[0], f[1], tuple(-x for x in f[1]) if rnd.random() < 0.5 else f[1])
        rows.add((J, f))
    return sorted(rows)
st = dict(rows=0, new=0, pal=0, pal_new=0, outside=0, class_breaks=0, prime_split=0, D_rows=0, outside_without_C=0)
split = dict(B=0, C=0, D=0, D_planted=0)
for J, f in g3_rows(20261006, 1200):
    agree, pal = verdict(hamiltonian(J, f, "XYZ"))
    inA, inB, inC, inD = class_A(f), class_B(J, f), class_C(J, f), class_D(J, f, True)
    new = not letter_turnable(f)
    st["rows"] += 1; st["new"] += new; st["pal"] += pal; st["pal_new"] += pal and new; st["prime_split"] += not agree
    st["D_rows"] += inD
    if pal and new:
        for k, v in (("B", inB), ("C", inC), ("D", inD)): split[k] += v
        split["D_planted"] += inD and J in [(0, 1, 3), (0, 2, -1), (0, 1, 1)]
    if pal and not (inA or inB or inC or inD): st["outside"] += 1
    if (inA or inB or inC) and not pal: st["class_breaks"] += 1
    if pal and not (inA or inB or inD): st["outside_without_C"] += 1
print("G3 sweep (distinct rows):", st)
print("G3 palindromic rows no turn about z makes letter rows, by class:", split)
check("G3 the two primes agree on every row", st["prime_split"] == 0)
check("G3 every palindromic row lies in A, B, C or D", st["outside"] == 0, f"{st['outside']} outside")
check("G3 every row of A, B or C pairs", st["class_breaks"] == 0, f"{st['class_breaks']} break")
check("G3 the sweep holds rows of D and palindromic rows no turn about z makes letter rows", st["D_rows"] > 0 and st["pal_new"] > 0)
check("G3 those rows split 31 B, 8 C, 9 D with no overlap, all D rows on the planted couplings",
      (split["B"], split["C"], split["D"], split["D_planted"]) == (31, 8, 9, 9) and sum(split[k] for k in "BCD") == st["pal_new"])
check("G3 control: with C left out, palindromic rows are reported outside", st["outside_without_C"] > 0, f"{st['outside_without_C']}")

# ---------------------------------------------------------------- G5
JE = (1, 2, 6); fE = ((0, 0, 0), (8, 0, 0), (-5, 3, 0))
WE = {k: sp.Integer(v) for k, v in {"XII": 12, "XIX": 30, "XIY": -10, "XXI": 15, "XXX": -30, "XXY": -2, "XYI": 5, "XYX": -8,
      "XZZ": -12, "YII": -4, "YIX": -10, "YIY": 30, "YXI": 5, "YXY": 24, "YYI": 15, "YYX": 6, "YYY": 10, "YZZ": 4}.items()}
HE = hamiltonian(JE, fE, "XYZ")
check("G5 locus E, row J = (1, 2, 6), f_u = 8X, f_v = -5X + 3Y: the eighteen-string W commutes with H, anticommutes with "
      "Z_d and squares to 4500 I (exact), and the row lies outside A, B, C and D",
      comm(HE, WE) == {} and acomm(Zd, WE) == {} and pmul(WE, WE) == {"III": sp.Integer(4500)}
      and not (class_A(fE) or class_B(JE, fE) or class_C(JE, fE) or class_D(JE, fE, True)))
cE = [verdict(hamiltonian(JE, ((0, 0, 0), (8, 0, 0), (-5, 4, 0)), "XYZ")), verdict(hamiltonian((1, 2, 7), fE, "XYZ"))]
check("G5 control: f_v = -5X + 4Y and J_uv = 7 each break the row (two primes agreeing)", all(ag and not pl for ag, pl in cE))
s_, t_, r_ = JE; e2E = s_ * t_ + s_ * r_ + t_ * r_
check("G5 the row satisfies the three E relations exactly (6 = 6, f_u.f_v = -40 = -2 e2, 448 + 272 - 720 = 0)",
      (r_ * (s_ - t_) ** 2, s_ * t_ * (s_ + t_), dot(fE[1], fE[2]), e2E,
       (s_ + r_) * dot(fE[1], fE[1]) + (t_ + r_) * dot(fE[2], fE[2]) + 2 * (s_ + t_ + r_) * dot(fE[1], fE[2])) == (6, 6, -40, 20, 0))
check("G5 the dephased seat's Krylov determinant -(J_du - J_dv) e2 is 20 at the row (not blind)", -(s_ - t_) * e2E == 20)

# ---------------------------------------------------------------- G6
def ham_mod(J, f, p):
    """Heisenberg H with coefficients already reduced mod p (fields may hold square roots mod p)"""
    H = {}
    for (a_, b_), j in {(0, 1): J[0], (0, 2): J[1], (1, 2): J[2]}.items():
        if j % p:
            for P in "XYZ":
                s = ["I"] * 3; s[a_] = P; s[b_] = P; H["".join(s)] = j % p
    for l in range(3):
        for k, P in enumerate("XYZ"):
            if f[l][k] % p:
                s = ["I"] * 3; s[l] = P; w = "".join(s); H[w] = (H.get(w, 0) + f[l][k]) % p
    return H
def nullity_modp(Hm, cols, p):
    ip = sqrtm1(p); M = np.zeros((64, len(cols)), dtype=object)
    for j, P in enumerate(cols):
        for h, c in Hm.items():
            ph, w = smul(h, P); ph2, _ = smul(P, h)
            if ph == ph2: continue
            z = sp.expand(2 * ph); re, im = int(sp.re(z)), int(sp.im(z))
            M[IDX[w], j] = (M[IDX[w], j] + c * (re + im * ip)) % p
    r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, 64) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]; inv = pow(int(M[r, c]), p - 2, p); M[r] = [(x * inv) % p for x in M[r]]
        for i in range(64):
            if i != r and M[i, c]: M[i] = [(x - M[i, c] * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return len(cols) - r
def sqrt_modp(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    return sp.ntheory.residue_ntheory.sqrt_mod(a, p)
def e_point(s, t, a_, r=None):
    """the three E relations solved for r and f_v = (x, y, 0) given J_du = s, J_dv = t, f_u = (a, 0, 0); y^2 returned.
    With r given, only the second and third relations are solved (the control off the first)."""
    if r is None: r = s * t * (s + t) / (s - t) ** 2
    e2 = s * t + s * r + t * r; p1 = s + t + r
    x = -2 * e2 / a_
    return r, x, (4 * p1 * e2 - (s + r) * a_ * a_) / (t + r) - x * x
def e_verdict(J, fu, x, y2, p, dx=F(0)):
    y = sqrt_modp(y2.numerator * pow(y2.denominator, p - 2, p), p)
    if y is None: return None
    red = lambda q: q.numerator * pow(q.denominator, p - 2, p) % p
    Hm = ham_mod(tuple(red(F(j)) for j in J), ((0, 0, 0), (red(fu), 0, 0), (red(x + dx), y, 0)), p)
    return nullity_modp(Hm, NEAR, p), nullity_modp(Hm, FAR, p)
E_PRIMES = [998244353, 469762049, 167772161, 754974721, 1004535809, 2013265921]
e_ok = e_tot = e_ctrl = e_off_tot = e_off_break = 0
for s, t, a_ in [(1, 2, 8), (1, 2, 7), (1, 3, 5), (2, 1, 9), (1, -2, 3), (3, 1, 4), (1, 4, 6), (-1, 2, 5), (2, 5, 11), (3, -1, 7), (1, 5, 13), (2, 3, 9)]:
    s, t, a_ = F(s), F(t), F(a_)
    r, x, y2 = e_point(s, t, a_)
    if y2 <= 0: continue
    reads = [v for v in (e_verdict((s, t, r), a_, x, y2, p) for p in E_PRIMES) if v is not None][:2]
    if len(reads) < 2: continue
    ctrl = [v for v in (e_verdict((s, t, r), a_, x, y2, p, F(1, 10)) for p in E_PRIMES) if v is not None][:2]
    e_tot += 1; e_ok += reads[0] == reads[1] == (1, 1); e_ctrl += len(ctrl) == 2 and all(c[0] != c[1] or c[1] == 0 for c in ctrl)
    r1, x1, y21 = e_point(s, t, a_, r + 1)
    if y21 > 0:
        off = [v for v in (e_verdict((s, t, r1), a_, x1, y21, p) for p in E_PRIMES) if v is not None][:2]
        if len(off) == 2: e_off_tot += 1; e_off_break += all(c[0] != c[1] or c[1] == 0 for c in off)
print(f"G6 points on the three E relations read modulo two primes where y^2 has a root: {e_ok}/{e_tot} pair with near = far = 1; "
      f"{e_ctrl}/{e_tot} break when f_v moves by X/10")
check("G6 every point read on the three E relations pairs (near = far = 1 at two primes), and each breaks when f_v moves",
      e_tot >= 6 and e_ok == e_tot and e_ctrl == e_tot, f"{e_ok}/{e_tot}, controls {e_ctrl}/{e_tot}")
check("G6 control: points that keep the second and third relations at J_uv + 1 (off the first) break",
      e_off_tot >= 4 and e_off_break == e_off_tot, f"{e_off_break}/{e_off_tot}")

# ---------------------------------------------------------------- G4
c_rows = [((3, 6, -2), ((0, 0, 0), (1, 0, 0), (0, 1, 0))), ((1, -2, -2), ((0, 0, 0), (1, 0, 0), (0, 1, 0))),
          ((3, 6, -2), ((0, 0, 0), (3, 4, 0), (5, 0, 0)))]
def c_row_reads(J, f):
    h, x = verdict(hamiltonian(J, f, "XYZ")), verdict(hamiltonian(J, f, "XY"))
    return class_C(J, f) and not class_A(f) and h[0] and h[1] and x[0] and not x[1]
check("G4 XX + YY bonds: three rows of C outside A pair under Heisenberg bonds and break under XX + YY (two primes agreeing)",
      all(c_row_reads(J, f) for J, f in c_rows))
pairs = [((3, 4, 0), (-3, 4, 0)), ((3, 4, 0), (5, 0, 0)), ((1, 0, 0), (0, 1, 0)), ((4, 3, 0), (0, 5, 0)), ((1, 1, 0), (1, -1, 0)),
         ((5, 0, 0), (-4, 3, 0))]
res = [verdict(hamiltonian(J, ((0, 0, 0), fu, fv), "XY")) for J in ((2, 3, 0), (1, -4, 0), (5, 7, 0)) for fu, fv in pairs]
check("G4 XX + YY on the path u - d - v, f_d = 0, equal magnitudes perpendicular to z at six angles: every row of the grid pairs",
      all(ag and pl for ag, pl in res), f"{sum(pl for _, pl in res)}/{len(res)}")
unequal = verdict(hamiltonian((2, 3, 0), ((0, 0, 0), (3, 4, 0), (1, 0, 0)), "XY"))
check("G4 control: unequal magnitudes at a non-right angle break", unequal[0] and not unequal[1])

print("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
os.makedirs("simulations/results", exist_ok=True)
open("simulations/results/n3_offletter_locus_gate.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
