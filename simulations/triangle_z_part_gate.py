"""Gate: the Heisenberg triangle with fields that have a z part (PROOF_ONE_DEPHASED_SEAT, section
"The triangle with fields off the plane").

Setting: sites (d, u, v) = (0, 1, 2), Heisenberg bonds s·σ_d·σ_u + t·σ_d·σ_v + r·σ_u·σ_v with real couplings, at most
one zero, the jump Z on d alone, fields f_d, f_u, f_v vectors of R^3. q_u, q_v here are the z parts f_u^z, f_v^z of the
proof, and (i), (ii), (iii) its three relations. Exact throughout (sympy over Q(i), integer ranks), no threshold. Five
parts.

1. A z part on d. Tr(H·Z_d) = 8·f_d^z identically, the shortest odd word of F158 §(f5): no row with a z part on d
   pairs. Read exactly at rows (far end zero there).
2. The twin's spectrum. H_00 and H_11 (d's bit 0 and 1) are the dimer r·σ_u·σ_v + a·σ_u + b·σ_v with a = f_u + s·ẑ,
   b = f_v + t·ẑ and a' = f_u - s·ẑ, b' = f_v - t·ẑ (when f_d has no z part). The dimer's characteristic polynomial is
      x^4 - 2(3r^2 + |a|^2 + |b|^2) x^2 - 8r(a·b - r^2) x + (|a|^2 - |b|^2)^2 + 2r^2(|a|^2 + |b|^2) + 8r^2·a·b - 3r^4,
   so H_00 and H_11 are isospectral exactly when
      (i): s·q_u + t·q_v = 0,   (ii) (r != 0): t·q_u + s·q_v = 0,   (iii): (s·q_u - t·q_v)(|f_u|^2 + s^2 - |f_v|^2 - t^2) = 0,
   all checked symbolically; on the stratum s = t the hop space's dimension from the discriminant of H_01^dag·H_01.
3. The words. Tr(Z_d·H^2) = 16·(s·q_u + t·q_v) and Tr(Z_d·H^3) = 16·(3r(t·q_u + s·q_v) - s^2·q_u - t^2·q_v) for every
   field on d without a z part, symbolic, and the three-jump identity 3·Tr(Z_d·H^3) + Tr(Z_d·H·Z_d·H·Z_d·H) =
   192·r·(t·q_u + s·q_v), so r·(ii) is a combination of odd words. The determinants of the proof's case logic, taken from these
   polynomials: r != 0 with t = -s and q_u = q_v gives the word -32·s^2·q; r = 0 gives s·t·(s - t). Planted rows on
   the branch s·t + 3r(s + t) = 0, where Tr(Z_d·H^2) and Tr(Z_d·H^3) vanish, one of them with (i) and (iii) holding as well, so
   that (ii) alone decides, read exactly.
4. J_du = J_dv with a field on d. J_du = J_dv = t, f_d = (p, 0, 0). The hop space (I_d and Z_d parts of [H, W]) does not
   depend on f_u, f_v or r; the previous section's basis W0, W1', W_H, W2' at s = t is checked to solve the hop
   equations with z parts present, each element with a string the others lack. The X_d and Y_d parts on it are that
   section's twelve rows, unchanged, and eight rows that q_u·Z_u + q_v·Z_v adds, all pinned as the proof displays
   them; with a z part the kernel is exactly f_v = R_x·f_u (B off the plane with its axis along f_d, 2·X_d (x) SWAP·X_uX_v) and f_u = f_v (the
   singlet row S), each one-dimensional. Step 5's element SWAP_uv·(n·σ)^(x)3 commutes with H, anticommutes with Z_d
   and squares to 1 symbolically on the whole class (n = (c, sn, 0) rational on the circle, f_d = λ·n, f_v = R_n·f_u).
5. An exact census by a second route, the full commutator on the 64 Pauli strings: a grid of rows with z parts on u
   or v, rows planted on the three isospectral strata and on the word branch, and a seeded census over R^3 fields on
   all three sites against the combined rule: F158 §(f5) for a z part on d, this section for a z part on u or v, and
   otherwise the two sections with planar fields, which hold the letter section's rows without a z part. A far end read as zero at one prime is zero over Q (the rank modulo p is at most the rank over Q), so a
   row read with far = 0 is broken exactly, its near count an upper bound; every row with a nonzero far end at both
   primes is recounted with integer ranks (Bareiss, every division asserted exact). Four mutated rules (B without its
   axis, B without the reversed z part, the rule without B, a z part on d read as planar) are run against the same
   reads and must disagree. The far
   elements without pairing are sorted into their families and printed (measured), one row of them by height. The
   witness's pinned rows and its grid off the plane are read here too, the second route for its tests.

Run: python simulations/triangle_z_part_gate.py   (about three and a half minutes; writes simulations/results/triangle_z_part_gate.txt)
"""
import itertools
import random
import sys
from collections import Counter

import numpy as np
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = []
FAILS = []


def say(msg=""):
    print(msg)
    OUT.append(msg)


def check(cond, what):
    say(("   PASS  " if cond else "   FAIL  ") + what)
    if not cond:
        FAILS.append(what)


# ---------------------------------------------------------------------------------------------------------------
# Pauli algebra on (d, u, v)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
I2 = sp.eye(2)
PAULI = [X, Y, Z]
PM = [I2, X, Y, Z]


def op(k, A, n=3):
    m = sp.Matrix([[1]])
    for j in range(n):
        m = sp.kronecker_product(m, A if j == k else I2)
    return m


def smat(st):
    m = sp.Matrix([[1]])
    for a in st:
        m = sp.kronecker_product(m, PM[a])
    return m


def triangle(J, fd, fu, fv):
    H = sp.zeros(8)
    for c, (a, b) in zip(J, [(0, 1), (0, 2), (1, 2)]):
        H += c * sum((op(a, Q) * op(b, Q) for Q in PAULI), sp.zeros(8))
    for k, f in enumerate((fd, fu, fv)):
        H += sum((f[i] * op(k, PAULI[i]) for i in range(3)), sp.zeros(8))
    return H


s, t, r, p = sp.symbols('s t r p', real=True)
p1, p2, hz = sp.symbols('p1 p2 hz', real=True)
u1, u2, qu, v1, v2, qv = sp.symbols('u1 u2 q_u v1 v2 q_v', real=True)
Zd = op(0, Z)

# ---------------------------------------------------------------------------------------------------------------
say("1. A z part on d: the shortest odd word of F158 §(f5)")
H_full = triangle((s, t, r), (p1, p2, hz), (u1, u2, qu), (v1, v2, qv))
w1 = sp.expand((Zd * H_full).trace())
check(w1 == 8 * hz, f"Tr(H·Z_d) = 8·f_d^z identically in every coupling and field (read {w1})")

# ---------------------------------------------------------------------------------------------------------------
say()
say("2. The twin's spectrum: H_00 and H_11 are the dimer in the conditioned fields")
a1, a2, a3, b1, b2, b3, x = sp.symbols('a1 a2 a3 b1 b2 b3 x')
Hdim = r * sum((sp.kronecker_product(Q, Q) for Q in PAULI), sp.zeros(4)) \
    + sum((av * sp.kronecker_product(PAULI[i], I2) for i, av in enumerate((a1, a2, a3))), sp.zeros(4)) \
    + sum((bv * sp.kronecker_product(I2, PAULI[i]) for i, bv in enumerate((b1, b2, b3))), sp.zeros(4))
A2, B2, AB = a1**2 + a2**2 + a3**2, b1**2 + b2**2 + b3**2, a1 * b1 + a2 * b2 + a3 * b3
cp = sp.expand((x * sp.eye(4) - Hdim).det())
cp_closed = sp.expand(x**4 - 2 * (3 * r**2 + A2 + B2) * x**2 - 8 * r * (AB - r**2) * x
                      + (A2 - B2)**2 + 2 * r**2 * (A2 + B2) + 8 * r**2 * AB - 3 * r**4)
check(sp.expand(cp - cp_closed) == 0, "the dimer's characteristic polynomial is the closed form (three invariants)")

H_planar_d = triangle((s, t, r), (p1, p2, 0), (u1, u2, qu), (v1, v2, qv))
blk = lambda M, i, j: M.extract([4 * i + k for k in range(4)], [4 * j + k for k in range(4)])
cond = {}
for bit, sg in ((0, 1), (1, -1)):
    fa = (u1, u2, qu + sg * s)
    fb = (v1, v2, qv + sg * t)
    D = Hdim.subs({a1: fa[0], a2: fa[1], a3: fa[2], b1: fb[0], b2: fb[1], b3: fb[2]})
    cond[bit] = (fa, fb)
    check(sp.expand(blk(H_planar_d, bit, bit) - D) == sp.zeros(4),
          f"H_{bit}{bit} is the dimer with a = f_u {'+' if sg > 0 else '-'} s·ẑ, b = f_v {'+' if sg > 0 else '-'} t·ẑ")
inv = lambda fa, fb: (sum(c**2 for c in fa), sum(c**2 for c in fb), sum(fa[i] * fb[i] for i in range(3)))
(sa, ta, ua), (sb, tb, ub) = inv(*cond[0]), inv(*cond[1])
FU2, FV2 = u1**2 + u2**2 + qu**2, v1**2 + v2**2 + qv**2
check(sp.expand((sa + ta) - (sb + tb) - 4 * (s * qu + t * qv)) == 0, "|a|^2 + |b|^2 - |a'|^2 - |b'|^2 = 4·(s·q_u + t·q_v)  (i)")
check(sp.expand(ua - ub - 2 * (t * qu + s * qv)) == 0, "a·b - a'·b' = 2·(t·q_u + s·q_v)  (ii)")
check(sp.expand((sa - ta)**2 - (sb - tb)**2 - 8 * (s * qu - t * qv) * (FU2 + s**2 - FV2 - t**2)) == 0,
      "(|a|^2 - |b|^2)^2 - (|a'|^2 - |b'|^2)^2 = 8·(s·q_u - t·q_v)(|f_u|^2 + s^2 - |f_v|^2 - t^2)  (iii)")
say("   so with r != 0 isospectral <=> (i), (ii), (iii); with r = 0 (the x^1 coefficient and the r^2 terms drop) <=> (i), (iii)")
# the hop space on the whole stratum s = t: H_01^dag H_01 has a simple spectrum
pp_, tt_ = sp.symbols('pp tt', positive=True)
sgm = sp.Matrix([[0, 0], [1, 0]])                     # sigma^- = |1><0|
H01 = pp_ * sp.eye(4) + 2 * tt_ * (sp.kronecker_product(sgm, I2) + sp.kronecker_product(I2, sgm))
lam = sp.symbols('lam')
disc = sp.factor(sp.discriminant(sp.expand((lam * sp.eye(4) - H01.H * H01).det()), lam))
check(sp.expand(disc - 2**24 * pp_**6 * tt_**14 * (4 * pp_**4 + 13 * pp_**2 * tt_**2 + 32 * tt_**4)) == 0,
      f"s = t: disc of H_01^dag·H_01 = 2^24·p^6·t^14·(4p^4 + 13p^2·t^2 + 32t^4), each factor positive at p, t != 0 (read {disc})")

# ---------------------------------------------------------------------------------------------------------------
say()
say("3. The words Tr(Z_d·H^2) and Tr(Z_d·H^3), any field on d without a z part")
H2 = sp.expand(H_planar_d * H_planar_d)
w2 = sp.expand((Zd * H2).trace())
w3 = sp.expand((Zd * H2 * H_planar_d).trace())
check(sp.expand(w2 - 16 * (s * qu + t * qv)) == 0, "Tr(Z_d·H^2) = 16·(s·q_u + t·q_v), 16 times (i)")
W3 = 16 * (3 * r * (t * qu + s * qv) - s**2 * qu - t**2 * qv)
check(sp.expand(w3 - W3) == 0, "Tr(Z_d·H^3) = 16·(3r(t·q_u + s·q_v) - s^2·q_u - t^2·q_v)")
w111 = sp.expand((Zd * H_planar_d * Zd * H_planar_d * Zd * H_planar_d).trace())
check(sp.expand(3 * w3 + w111 - 192 * r * (t * qu + s * qv)) == 0,
      "3·Tr(Z_d·H^3) + Tr(Z_d·H·Z_d·H·Z_d·H) = 192·r·(t·q_u + s·q_v): r·(ii) is a combination of odd words, one with three jumps")
q = sp.symbols('q')
check(sp.expand(w3.subs({t: -s, qu: q, qv: q})) == -32 * s**2 * q,
      "r != 0, t = -s, (i) => q_u = q_v = q: the word reads -32·s^2·q (nonzero)")
coef = lambda P: [sp.expand(P).coeff(qu), sp.expand(P).coeff(qv)]
E1p, E2p = s * qu + t * qv, t * qu + s * qv
check(sp.factor(sp.Matrix([coef(w2 / 16), coef(w3.subs(r, 0) / 16)]).det()) == sp.factor(-s * t * (t - s)),
      "r = 0: Tr(Z_d·H^2) and Tr(Z_d·H^3) as forms in (q_u, q_v), determinant -s·t·(t - s)")
check(sp.factor(sp.Matrix([coef(E1p), coef(sp.expand(ua - ub) / 2)]).det()) == sp.factor((s - t) * (s + t)),
      "r != 0: (i) and (ii) as forms in (q_u, q_v), determinant s^2 - t^2")
check(sp.factor(sp.Matrix([coef(w2 / 16), coef(w3 / 16)]).det()) == sp.factor((s - t) * (s * t + 3 * r * (s + t))),
      "the two one-jump words alone: determinant (s - t)(s·t + 3r(s + t)), so they leave a branch")
say("   case logic (argued in the proof): r != 0 => s = +-t; s = -t dies by the word; s = t => q_v = -q_u and by (iii)")
say("   |f_u| = |f_v|; r = 0 => s = t (both nonzero there) => the same.")
Jb = (6, -3, 2)
check(Jb[0] * Jb[1] + 3 * Jb[2] * (Jb[0] + Jb[1]) == 0, "J = (6, -3, 2) lies on s·t + 3r(s + t) = 0")
for (fu_, fv_) in (((0, 0, 1), (0, 0, 2)), ((1, 0, 2), (4, 0, 4))):
    sub = {s: 6, t: -3, r: 2, u1: fu_[0], u2: fu_[1], qu: fu_[2], v1: fv_[0], v2: fv_[1], qv: fv_[2]}
    iii = ((s * qu - t * qv) * (FU2 + s**2 - FV2 - t**2)).subs(sub)
    check(w2.subs(sub) == 0 and w3.subs(sub) == 0 and E2p.subs(sub) != 0,
          f"f_u = {fu_}, f_v = {fv_}: both one-jump words vanish, (ii) = {E2p.subs(sub)}, (iii) = {iii}, "
          f"Tr(Z_d·H·Z_d·H·Z_d·H) = {w111.subs(sub)}")
check(((s * qu - t * qv) * (FU2 + s**2 - FV2 - t**2)).subs({s: 6, t: -3, r: 2, u1: 1, u2: 0, qu: 2, v1: 4, v2: 0, qv: 4}) == 0,
      "at f_u = (1, 0, 2), f_v = (4, 0, 4) (i) and (iii) hold too: (ii) alone decides that row")

# ---------------------------------------------------------------------------------------------------------------
say()
say("4. J_du = J_dv with a field on d: the previous section's hop space and its table, with z parts")
pn = sp.symbols('p', nonzero=True)
tn = sp.symbols('t', nonzero=True)
S3 = list(itertools.product(range(4), repeat=3))
AM = {st: smat(st) for st in S3}
FAR = [st for st in S3 if st[0] in (1, 2)]
FM = [AM[st] for st in FAR]
lbl = 'IXYZ'
NAME = {''.join(lbl[a] for a in st): st for st in S3}


def comps(M):
    return {st: sp.expand((AM[st] * M).trace() / 8) for st in S3}


def S_(name):
    return AM[NAME[name]]


# the committed basis at s = t (PROOF_ONE_DEPHASED_SEAT, the triangle with a planar field on d, step 4)
W0 = S_('XXX')
W1p = S_('XYY') + S_('XZZ') - S_('XII')
WH = pn * S_('XII') + tn * (S_('XXI') + S_('YYI') + S_('XIX') + S_('YIY'))
W2p = pn * (S_('XXI') + S_('XIX')) + tn * (S_('XYY') + S_('YXY') + S_('YYX'))
basis = [W0, W1p, WH, W2p]
Hd = tn * sum((op(0, Q) * op(1, Q) + op(0, Q) * op(2, Q) for Q in PAULI), sp.zeros(8)) + pn * op(0, X)
Huv = r * sum((op(1, Q) * op(2, Q) for Q in PAULI), sp.zeros(8)) \
    + sum((f * op(1, PAULI[i]) for i, f in enumerate((u1, u2, qu))), sp.zeros(8)) \
    + sum((f * op(2, PAULI[i]) for i, f in enumerate((v1, v2, qv))), sp.zeros(8))
H4 = Hd + Huv
hop_rows = [st for st in S3 if st[0] in (0, 3)]
check(all(comps(H4 * B - B * H4)[st] == 0 for B in basis for st in hop_rows),
      "W0, W1', W_H, W2' solve the hop equations (the parts along I_d and Z_d vanish) with z parts present")
private = ['XXX', 'XZZ', 'YIY', 'YXY']
coef = lambda B, name: sp.expand((S_(name) * B).trace() / 8)
check(all((coef(basis[k], private[k]) != 0) and all(coef(basis[j], private[k]) == 0 for j in range(4) if j != k)
          for k in range(4)), "each has a string the other three lack (XXX, XZZ, YIY, YXY): independent at t, p != 0")
ns_dim = len(sp.Matrix([[comps(Hd * B - B * Hd)[st] for B in FM] for st in hop_rows]).nullspace())
check(ns_dim == 4, f"the hop space at s = t is four-dimensional at generic t, p ({ns_dim}; on the whole stratum by the "
      "previous section's discriminant)")
let_rows = [st for st in S3 if st[0] in (1, 2)]
cc = sp.symbols('c0 c1 cH c2')
tabl = {}
for st in let_rows:
    row = [sp.factor(comps(H4 * B - B * H4)[st] / (2 * sp.I)) for B in basis]
    if any(e != 0 for e in row):
        tabl[''.join(lbl[a] for a in st)] = row
check(all(not e.has(r) for row in tabl.values() for e in row), "J_uv enters none of the rows")
lin = lambda row: sp.expand(sum(e * c for e, c in zip(row, cc)))
c0, c1, cH, c2 = cc
ux, uy, vx, vy = u1, u2, v1, v2
planar = {  # the twelve rows of the previous section's step 4, unchanged by z parts
    'YIZ': tn * (pn + vx) * cH, 'YZI': tn * (pn + ux) * cH,
    'YXZ': tn**2 * cH + tn * (pn + vx) * c2, 'YZX': tn**2 * cH + tn * (pn + ux) * c2,
    'YYZ': -tn * vy * c2, 'YZY': -tn * uy * c2,
    'XIZ': -vy * (tn * cH + pn * c2), 'XZI': -uy * (tn * cH + pn * c2),
    'XXZ': -vy * c0 + uy * c1, 'XZX': -uy * c0 + vy * c1,
    'XYZ': (vx - ux) * c1 - tn**2 * cH + tn * vx * c2, 'XZY': (ux - vx) * c1 - tn**2 * cH + tn * ux * c2}
zrows = {  # the eight rows q_u·Z_u + q_v·Z_v adds
    'YXI': -tn * qu * cH, 'YIX': -tn * qv * cH,
    'YXX': -tn * (qu + qv) * c2, 'YYY': tn * (qu + qv) * c2,
    'XYI': qu * (tn * cH + pn * c2), 'XIY': qv * (tn * cH + pn * c2),
    'XXY': qv * c0 - qu * c1 - tn * qu * c2, 'XYX': qu * c0 - qv * c1 - tn * qv * c2}
check(set(tabl) == set(planar) | set(zrows), f"the table has exactly the twelve planar rows and eight z rows ({len(tabl)} rows)")
check(all(sp.expand(lin(tabl[k]) - v) == 0 for k, v in planar.items()), "the twelve planar rows read as displayed there")
check(all(sp.expand(lin(tabl[k]) - v) == 0 for k, v in zrows.items()), "the eight z rows read as displayed")
for k, v in zrows.items():
    say(f"     {k[0]}_d·{k[1:]}: {v}")
say("   elimination (argued in the proof): Y_d·XI, Y_d·IX give c_H = 0; X_d·YI, X_d·IY then c2 = 0 (p != 0); X_d·XY,")
say("   X_d·YX on (c0, c1) have determinant q_u^2 - q_v^2; q_u = -q_v gives c1 = -c0 and with X_d·XZ, X_d·YZ the axis,")
say("   q_u = q_v gives c1 = c0 and f_u = f_v (S).")
blk2 = sp.Matrix([[sp.diff(zrows['XXY'], c0), sp.diff(zrows['XXY'], c1)], [sp.diff(zrows['XYX'], c0), sp.diff(zrows['XYX'], c1)]])
check(sp.factor(blk2.det()) == sp.factor((qu - qv) * (qu + qv)), f"X_d·XY, X_d·YX on (c0, c1): determinant {sp.factor(blk2.det())}")
TAB = sp.Matrix([[sp.diff(lin(row), c) for c in cc] for row in tabl.values()])


def kernel_at(sub):
    return TAB.subs(sub).applyfunc(sp.simplify).nullspace()


check(len(TAB.nullspace()) == 0, "at generic fields the kernel is zero")
for name, sub, want in (
        ("q_v = -q_u, f_v = R_x f_u (B off the plane, axis along f_d)", {qv: -qu, v1: u1, v2: -u2}, 1),
        ("q_v = q_u, f_v = f_u (the singlet row S)", {qv: qu, v1: u1, v2: u2}, 1),
        ("q_v = -q_u, planar parts equal, off the axis", {qv: -qu, v1: u1, v2: u2}, 0),
        ("q_v = -q_u, the mirror image turned by a quarter", {qv: -qu, v1: u2, v2: u1}, 0),
        ("q_v = q_u, planar parts mirrored", {qv: qu, v1: u1, v2: -u2}, 0)):
    kn = kernel_at(sub)
    check(len(kn) == want, f"{name}: kernel dimension {len(kn)} at generic remaining parameters")
kB = kernel_at({qv: -qu, v1: u1, v2: -u2})[0]
check(kB[2] == 0 and kB[3] == 0 and sp.simplify(kB[0] + kB[1]) == 0, "B off the plane: c_H = c2 = 0, c1 = -c0")
WB = W0 - W1p
SWAP = sp.Matrix(8, 8, lambda i, j: 1 if (i >> 2) == (j >> 2) and (((i >> 1) & 1) == (j & 1)) and ((i & 1) == ((j >> 1) & 1)) else 0)
check(WB == 2 * SWAP * op(0, X) * op(1, X) * op(2, X), "W0 - W1' = 2·X_d (x) SWAP_uv·X_uX_v, the element of B off the plane")
kS = kernel_at({qv: qu, v1: u1, v2: u2})[0]
check(kS[2] == 0 and kS[3] == 0 and sp.simplify(kS[0] - kS[1]) == 0, "S: c_H = c2 = 0, c1 = c0")
check(W0 + W1p == -4 * op(0, X) * sp.kronecker_product(sp.eye(2), (sp.eye(4) - sum((sp.kronecker_product(Q, Q) for Q in PAULI), sp.zeros(4))) / 4),
      "W0 + W1' = -4·X_d (x) P_s, P_s = (1 - σ_u·σ_v)/4 the singlet projector of u and v")
say("   step 5 on the whole class: n = (c, sn, 0) with c = (1 - m^2)/(1 + m^2), sn = 2m/(1 + m^2), f_d = λ·n, f_v = R_n·f_u")
mm, lm = sp.symbols('m lambda', real=True)
cn_, sn_ = (1 - mm**2) / (1 + mm**2), 2 * mm / (1 + mm**2)
nvec = (cn_, sn_, 0)
fuv = (u1, u2, qu)
dn = cn_ * u1 + sn_ * u2
fvR = (2 * dn * cn_ - u1, 2 * dn * sn_ - u2, -qu)
nsig = lambda k: cn_ * op(k, X) + sn_ * op(k, Y)
W5 = SWAP * nsig(0) * nsig(1) * nsig(2)
H5 = triangle((s, s, r), (lm * cn_, lm * sn_, 0), fuv, fvR)
check((H5 * W5 - W5 * H5).applyfunc(sp.simplify) == sp.zeros(8), "[H, SWAP_uv·(n·σ)^(x)3] = 0 identically")
check((Zd * W5 + W5 * Zd).applyfunc(sp.simplify) == sp.zeros(8), "{Z_d, SWAP_uv·(n·σ)^(x)3} = 0")
check((W5 * W5).applyfunc(sp.simplify) == sp.eye(8), "(SWAP_uv·(n·σ)^(x)3)^2 = 1, so it is invertible")

# ---------------------------------------------------------------------------------------------------------------
say()
say("5. The census: the full commutator on Pauli strings, exact ranks, against the rule")


def _mul1(a, b):
    if a == 0:
        return 0, b
    if b == 0:
        return 0, a
    if a == b:
        return 0, 0
    return (1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - a - b


IDX = {st: k for k, st in enumerate(S3)}
NEAR = [st for st in S3 if st[0] in (0, 3)]


def comm_cols(terms, basis_):
    cols = []
    for b in basis_:
        col = [0] * 64
        for (a, h) in terms:
            ph, out = 0, []
            for x_, y_ in zip(a, b):
                qq, cc = _mul1(x_, y_)
                ph += qq
                out.append(cc)
            if ph % 2 == 1:
                col[IDX[tuple(out)]] += (1 if (ph - 1) % 4 == 0 else -1) * h
        cols.append(col)
    return cols


def rank_int(cols):
    M = [list(c) for c in cols]
    n, m = len(M), len(M[0])
    rk, prev = 0, 1
    for cc in range(m):
        piv = next((i for i in range(rk, n) if M[i][cc] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(rk + 1, n):
            row = []
            for j in range(m):
                num = M[rk][cc] * M[i][j] - M[i][cc] * M[rk][j]
                if num % prev:
                    raise ArithmeticError("Bareiss division not exact")
                row.append(num // prev)
            M[i] = row
        prev = M[rk][cc]
        rk += 1
        if rk == n:
            break
    return rk


PRIMES = (1000003, 998244353)


def rank_p(cols, pr):
    A = np.array(cols, dtype=np.int64) % pr
    n, m = A.shape
    rk = 0
    for c in range(m):
        piv = np.nonzero(A[rk:, c])[0]
        if piv.size == 0:
            continue
        i = rk + piv[0]
        A[[rk, i]] = A[[i, rk]]
        A[rk] = (A[rk] * pow(int(A[rk, c]), pr - 2, pr)) % pr
        f = A[:, c].copy()
        f[rk] = 0
        A = (A - np.outer(f, A[rk])) % pr
        rk += 1
        if rk == n:
            break
    return rk


def terms_of(J, fd, fu, fv):
    T = []
    for (i, j, c) in ((0, 1, J[0]), (0, 2, J[1]), (1, 2, J[2])):
        if c:
            for l in (1, 2, 3):
                a = [0, 0, 0]
                a[i] = a[j] = l
                T.append((tuple(a), c))
    for site, f in ((0, fd), (1, fu), (2, fv)):
        for l in (1, 2, 3):
            if f[l - 1]:
                a = [0, 0, 0]
                a[site] = l
                T.append((tuple(a), f[l - 1]))
    return T


def ends(J, fd, fu, fv):
    """(near, far, exact): far = 0 at one prime is exact (rank mod p <= rank over Q); far > 0 at both primes is
    recounted with integer ranks, and then both counts are exact."""
    T = terms_of(J, fd, fu, fv)
    cn, cf = comm_cols(T, NEAR), comm_cols(T, FAR)
    farp = [32 - rank_p(cf, pr) for pr in PRIMES]
    if min(farp) == 0:
        return min(32 - rank_p(cn, pr) for pr in PRIMES), 0, False
    return 32 - rank_int(cn), 32 - rank_int(cf), True


def n2(f):
    return sum(c * c for c in f)


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot2(a, b):
    return a[0] * b[0] + a[1] * b[1]


def mirror(f, axis):
    """the pi turn about the planar axis (axis[0], axis[1], 0), scaled by |axis|^2 to stay integer"""
    nn = dot2(axis, axis)
    d = dot2(f, axis)
    return (2 * d * axis[0] - nn * f[0], 2 * d * axis[1] - nn * f[1], -nn * f[2]), nn


def Bz(J, fd, fu, fv):
    """J_du = J_dv and f_v = R_n f_u for a pi turn about an axis n ⊥ z, with f_d on n (any n if f_d = 0)"""
    if J[0] != J[1] or fu[2] != -fv[2]:
        return False
    if fd[0] or fd[1]:
        g, nn = mirror(fu, fd)
        return g == (nn * fv[0], nn * fv[1], nn * fv[2])
    return fu[0]**2 + fu[1]**2 == fv[0]**2 + fv[1]**2


def rule(J, fd, fu, fv):
    """the combined rule (F158 §(f5), this section, the two planar sections): the class a row pairs in, or None"""
    sJ, tJ, rJ = J
    if fd[2]:
        return None                                        # a z part on d: F158 §(f5)
    if fu[2] or fv[2]:
        return 'Boff' if Bz(J, fd, fu, fv) else None       # this section: B off the plane
    if fd[0] or fd[1]:                                     # a planar field on d
        if cross2(fd, fu) == 0 and cross2(fd, fv) == 0:
            return 'A'
        if sJ == tJ and Bz(J, fd, fu, fv):
            return 'B'
        if sJ == 0 and tJ and rJ and dot2(fd, fu) == 0 and fv[:2] == (-(fd[0] + fu[0]), -(fd[1] + fu[1])):
            return 'D'
        if tJ == 0 and sJ and rJ and dot2(fd, fv) == 0 and fu[:2] == (-(fd[0] + fv[0]), -(fd[1] + fv[1])):
            return "D'"
        return None
    if cross2(fu, fv) == 0:                                # planar fields, none on d
        return 'A'
    if sJ == tJ and n2(fu) == n2(fv):
        return 'B'
    e2 = sJ * tJ + sJ * rJ + tJ * rJ
    if e2 == 0 and n2(fu) == n2(fv):
        return 'C'
    R1 = rJ * (sJ - tJ)**2 - sJ * tJ * (sJ + tJ)
    R2 = dot2(fu, fv) + 2 * e2
    R3 = (sJ + rJ) * n2(fu) + (tJ + rJ) * n2(fv) + 2 * (sJ + tJ + rJ) * dot2(fu, fv)
    if R1 == 0 and R2 == 0 and R3 == 0:
        return 'E'
    return None


tab = Counter()
bad = []
seen = set()
nexact = [0]
ROWS = []          # (label, key, near, far, exact)


def read(label, J, fd, fu, fv):
    if list(J).count(0) > 1:
        return
    key = (tuple(J), tuple(fd), tuple(fu), tuple(fv))
    near, far, exact = ends(J, fd, fu, fv)
    nexact[0] += exact
    seen.add(key)
    k = rule(J, fd, fu, fv)
    pairs = far > 0 and near == far
    ROWS.append((label, key, near, far, exact))
    tab[(label, k or '-', (near, far) if exact else ('<=' + str(near), 0))] += 1
    if pairs != (k is not None):
        bad.append((label, key, near, far, k))


# (a) the grid: z parts on u or v, entries in {-1, 0, 1}, f_d in {0, x}
Jset = [-1, 1, 2, 0]
F3 = list(itertools.product((-1, 0, 1), repeat=3))
for J in itertools.product(Jset, repeat=3):
    for fd in ((0, 0, 0), (1, 0, 0)):
        for fu in F3:
            for fv in F3:
                if fu[2] or fv[2]:
                    read('grid', J, fd, fu, fv)
say(f"   grid done ({len(seen)} distinct rows so far)")

# (b) planted rows on the strata and the branch
random.seed(11)
V = list(itertools.product(range(-3, 4), repeat=2))
for (sJ, tJ) in [(1, 2), (2, 1), (1, -2), (-1, 3), (3, 1), (2, -3)]:            # r = 0 alone
    for k in (1, -1, 2):
        need = (k * k - 1) * (sJ * sJ - tJ * tJ)
        sols = [(a, b) for a in V for b in V if n2(a) - n2(b) == need]
        random.shuffle(sols)
        for a, b in sols[:6]:
            for fd in ((0, 0, 0), (1, 0, 0), (2, 1, 0)):
                read('r = 0 alone', (sJ, tJ, 0), fd, (a[0], a[1], k * tJ), (b[0], b[1], -k * sJ))
for sJ in (1, 2, -1):                                                             # s = -t, equal z parts
    for rJ in (1, -1, 2, 0):
        for qq in (1, 2, -1):
            sols = [(a, b) for a in V for b in V if n2(a) == n2(b) and a != (0, 0)]
            random.shuffle(sols)
            for a, b in sols[:4]:
                for fd in ((0, 0, 0), (1, 0, 0), (2, 1, 0)):
                    read('s = -t, equal z parts', (sJ, -sJ, rJ), fd, (a[0], a[1], qq), (b[0], b[1], qq))
for sJ in (1, 2, -1):                                                             # s = t, on and off the axis
    for rJ in (1, -1, 2, 3, 0):
        for qq in (1, 2, -1):
            for _ in range(3):
                pu = random.choice([v for v in V if v != (0, 0)])
                for fd in ((1, 0, 0), (2, 1, 0), (-1, 2, 0)):
                    nn = dot2(fd, fd)
                    fu = (nn * pu[0], nn * pu[1], nn * qq)
                    fvB, _ = mirror(fu, fd)
                    fvB = tuple(c // nn for c in fvB)
                    read('s = t, on the axis', (sJ, sJ, rJ), fd, fu, fvB)
                    read('s = t, planar parts equal', (sJ, sJ, rJ), fd, fu, (fu[0], fu[1], -fu[2]))
                    read('s = t, turned off the axis', (sJ, sJ, rJ), fd, fu, (-fvB[1], fvB[0], fvB[2]))
                read('s = t, no field on d', (sJ, sJ, rJ), (0, 0, 0), (pu[0], pu[1], qq), (pu[1], pu[0], -qq))
for qq in (1, 2, -1):                                                             # the word branch
    for fd in ((0, 0, 0), (1, 0, 0), (2, 1, 0)):
        for a, b in [((0, 0), (0, 0)), ((1, 0), (2, 1)), ((1, 2), (0, 3)), ((3, 0), (1, -2))]:
            read('word branch (6,-3,2)', (6, -3, 2), fd, (a[0], a[1], qq), (b[0], b[1], 2 * qq))
            read('word branch (6,-3,2)', (-6, 3, -2), fd, (a[0], a[1], qq), (b[0], b[1], 2 * qq))

for fd in ((0, 0, 0), (1, 0, 0), (2, 1, 0)):                                       # the branch, (ii) alone
    for Jb_ in ((6, -3, 2), (-6, 3, -2)):
        read('word branch, (ii) alone', Jb_, fd, (1, 0, 2), (4, 0, 4))

# (c) the combined rule over R^3: seeded random rows and planted rows of every class, z parts toggled
JJ = [-3, -2, -1, 1, 2, 3]
for _ in range(1500):
    J = (random.choice(JJ + [0]), random.choice(JJ), random.choice(JJ))
    f = [tuple(random.randint(-2, 2) for _ in range(3)) for _ in range(3)]
    for k in range(3):
        if random.random() < 0.5:
            f[k] = (f[k][0], f[k][1], 0)
    read('random R^3', J, *f)
for _ in range(60):
    sJ, rJ = random.choice(JJ), random.choice(JJ)
    m = random.choice([(1, 0), (2, 1), (1, -1), (0, 1)])
    al, be, ga = random.randint(-2, 2), random.randint(-2, 2), random.randint(-2, 2)
    read('A planted', (sJ, random.choice(JJ), rJ), (ga * m[0], ga * m[1], 0), (al * m[0], al * m[1], 0), (be * m[0], be * m[1], 0))
    read('A, z on d', (sJ, random.choice(JJ), rJ), (m[0], m[1], 1), (al * m[0], al * m[1], 0), (be * m[0], be * m[1], 0))
    pu = random.choice([v for v in V if v != (0, 0)])
    read('B planted', (sJ, sJ, rJ), (0, 0, 0), (pu[0], pu[1], 0), (pu[1], pu[0], 0))
    read('B, z on d', (sJ, sJ, rJ), (0, 0, 2), (pu[0], pu[1], 0), (pu[1], pu[0], 0))
    qq = random.choice([1, 2, -1, 3])
    read('B off the plane planted', (sJ, sJ, rJ), (0, 0, 0), (pu[0], pu[1], qq), (pu[1], pu[0], -qq))
    read('B off the plane, z on d', (sJ, sJ, rJ), (0, 0, 1), (pu[0], pu[1], qq), (pu[1], pu[0], -qq))
    read('B off the plane, q_v = q_u', (sJ, sJ, rJ), (0, 0, 0), (pu[0], pu[1], qq), (pu[1], pu[0], qq))
    tc = random.choice([2, 3, -2, -3, 4])                       # C: the chart (1 + t, t(1 + t), -t)
    Jc = (1 + tc, tc * (1 + tc), -tc)
    read('C planted', Jc, (0, 0, 0), (pu[0], pu[1], 0), (pu[0], -pu[1], 0))
    read('C, z on u', Jc, (0, 0, 0), (pu[0], pu[1], 1), (pu[0], -pu[1], 0))
    pd = random.choice([(1, 0), (2, 1), (-1, 3)])                 # D: J_du = 0, f_u ⊥ f_d, sum zero
    kk_ = random.choice([1, -1, 2])
    fuD = (-kk_ * pd[1], kk_ * pd[0], 0)
    fvD = (-(pd[0] + fuD[0]), -(pd[1] + fuD[1]), 0)
    tD, rD = random.choice(JJ), random.choice(JJ)
    read('D planted', (0, tD, rD), (pd[0], pd[1], 0), fuD, fvD)
    read("D' planted", (tD, 0, rD), (pd[0], pd[1], 0), fvD, fuD)
    read("D', z on u", (tD, 0, rD), (pd[0], pd[1], 0), (fvD[0], fvD[1], 1), fuD)
    read('D, z on v', (0, tD, rD), (pd[0], pd[1], 0), fuD, (fvD[0], fvD[1], 1))
for Jr, fu, fv in [((1, 2, 6), (8, 0, 0), (-5, 3, 0))]:           # E, the proposition's row
    read('E planted', Jr, (0, 0, 0), fu, fv)
    for qq in (1, -1, 2):
        read('E, z on u', Jr, (0, 0, 0), (fu[0], fu[1], qq), fv)
        read('E, z on both', Jr, (0, 0, 0), (fu[0], fu[1], qq), (fv[0], fv[1], -qq))

n_rows = sum(tab.values())
say(f"   {n_rows} reads of {len(seen)} distinct rows, {nexact[0]} recounted with integer ranks (a nonzero far end at both primes);")
say("   (near, far) by planted label and class of the combined rule ('-' for none; '<=n' a near count from the primes,")
say("   an upper bound, on rows whose far end is zero exactly):")
labels = sorted({(lab, k) for (lab, k, _) in tab})
for (lab, k) in labels:
    cells = {nf: c for (l2, k2, nf), c in tab.items() if (l2, k2) == (lab, k)}
    say(f"     {lab:30s} {k:3s} {dict(sorted(cells.items(), key=str))}")
check(not bad, f"every row: pairs exactly when the combined rule names a class ({len(bad)} disagree)")
for b in bad[:10]:
    say(f"     disagrees: {b}")
classes_met = {k for (_, k, _) in tab}
check({'A', 'B', 'C', 'D', "D'", 'E', 'Boff', '-'} <= classes_met, "each class of the combined rule is met")
zd_pairs = sum(1 for (lab, key, near, far, ex) in ROWS if key[1][2] and far > 0 and near == far)
zd_rows = sum(1 for (lab, key, near, far, ex) in ROWS if key[1][2])
check(zd_rows > 0 and zd_pairs == 0, f"the ends: none of the {zd_rows} reads with a z part on d pairs")


def disagree(rl):
    return sum(1 for (lab, key, near, far, ex) in ROWS if (far > 0 and near == far) != (rl(*key) is not None))


def Bz_any_axis(J, fd, fu, fv):
    return J[0] == J[1] and fu[2] == -fv[2] and fu[0]**2 + fu[1]**2 == fv[0]**2 + fv[1]**2


def Bz_any_sign(J, fd, fu, fv):
    return Bz(J, fd, fu, fv) or Bz(J, fd, fu, (fv[0], fv[1], -fv[2]))


for name, mutant in (
        ("B off the plane without its axis (any n when f_d != 0)",
         lambda J, fd, fu, fv: ('x' if Bz_any_axis(J, fd, fu, fv) else None) if (not fd[2] and (fu[2] or fv[2])) else rule(J, fd, fu, fv)),
        ("B off the plane without the reversed z part",
         lambda J, fd, fu, fv: ('x' if Bz_any_sign(J, fd, fu, fv) else None) if (not fd[2] and (fu[2] or fv[2])) else rule(J, fd, fu, fv)),
        ("the rule without B off the plane",
         lambda J, fd, fu, fv: None if (not fd[2] and (fu[2] or fv[2])) else rule(J, fd, fu, fv)),
        ("a z part on d read as if it were planar",
         lambda J, fd, fu, fv: rule(J, (fd[0], fd[1], 0), fu, fv))):
    nd = disagree(mutant)
    check(nd > 0, f"mutated rule, {name}: {nd} reads disagree")

say()
say("   far elements without pairing, sorted into families (measured on these reads):")
fam = Counter()
for (lab, key, near, far, ex) in ROWS:
    J, fd, fu, fv = key
    if far == 0 or near == far or fd[2] or not (fu[2] or fv[2]):
        continue
    if J[0] == J[1] and fu == fv:
        k = 'S: J_du = J_dv, f_u = f_v'
    elif J[0] == -J[1] and fu == tuple(-c for c in fv):
        k = 'J_dv = -J_du, f_v = -f_u'
    elif fd[:2] == (0, 0) and fu[:2] == (0, 0) and fv[:2] == (0, 0):
        k = 'every field along z'
    else:
        k = 'other'
    fam[('no field on d' if fd[:2] == (0, 0) else 'planar field on d', k)] += 1
for k in sorted(fam):
    say(f"     {k[0]:18s} {k[1]:28s} {fam[k]}")
check(fam and all(k[1] != 'other' for k in fam) and all(k[1].startswith('S') for k in fam if k[0] == 'planar field on d'),
      "with a z part on u or v, every far element without pairing read here lies in S (with a planar field on d: only S), "
      "J_dv = -J_du with f_v = -f_u, or a row with every field along z")
FARODD = [x for x in FAR if sum(1 for a in x if a in (1, 2)) % 2 == 1]
FAREVEN = [x for x in FAR if sum(1 for a in x if a in (1, 2)) % 2 == 0]
Teh = terms_of((-1, 1, -1), (0, 0, 0), (0, 0, -1), (0, 0, 1))
e_all, e_odd, e_even = (len(b) - rank_int(comm_cols(Teh, b)) for b in (FAR, FARODD, FAREVEN))
sJ, tJ, rJ = -1, 1, -1
kry = -(sJ - tJ) * (sJ * tJ + sJ * rJ + tJ * rJ)
check((e_all, e_odd, e_even) == (4, 2, 2) and kry != 0,
      f"J = (-1, 1, -1), f_u = -ẑ, f_v = ẑ: far {e_all}, of odd height {e_odd}, of even height {e_even}, at a seat that "
      f"is not blind (the Krylov factor -(J_du - J_dv)·e2 = {kry})")

say()
say("   the witness's pins (TrianglePalindromeLocusWitnessTests, rows off the plane), read here exactly:")
PINS = [((2, 2, 3), (1, 0, 0), (1, 2, 1), (1, -2, -1), 1, 1), ((2, 2, 3), (0, 0, 0), (1, 2, 1), (2, 1, -1), 1, 1),
        ((1, 1, 3), (2, 1, 0), (0, 0, 1), (0, 0, -1), 1, 1), ((2, 2, 0), (2, 1, 0), (5, 5, 5), (7, 1, -5), 1, 1),
        ((2, 2, 0), (2, 1, 0), (5, 5, 5), (7, 1, 5), 1, 0), ((2, 2, 3), (1, 0, 0), (1, 2, 1), (1, 2, 1), 2, 1),
        ((2, 2, 3), (1, 0, 0), (1, 2, 1), (2, 1, -1), 1, 0), ((2, -2, 3), (0, 0, 0), (1, 2, 1), (2, 1, 1), 1, 0),
        ((6, -3, 2), (0, 0, 0), (1, 0, 1), (2, 1, 2), 1, 0), ((6, -3, 2), (0, 0, 0), (1, 0, 2), (4, 0, 4), 1, 0),
        ((1, 3, 0), (0, 0, 0), (1, 2, 3), (2, 1, -1), 1, 0), ((2, 2, 3), (1, 0, 1), (1, 2, 1), (1, -2, -1), 1, 0),
        ((2, 2, 3), (1, 0, 1), (1, 2, 0), (1, -2, 0), 1, 0), ((1, 2, 3), (0, 0, 0), (0, 0, 1), (0, 0, 1), 4, 0),
        ((2, 2, 3), (0, 0, 0), (1, 2, 1), (1, 2, 1), 3, 2), ((2, 2, 3), (1, 0, 0), (1, 2, 1), (1, 2, 2), 1, 0)]
pin_ok = True
for (J, fd, fu, fv, n_, f_) in PINS:
    T = terms_of(J, fd, fu, fv)
    got = (32 - rank_int(comm_cols(T, NEAR)), 32 - rank_int(comm_cols(T, FAR)))
    pin_ok &= got == (n_, f_)
    say(f"     J = {J}, f_d = {fd}, f_u = {fu}, f_v = {fv}: (near, far) = {got}")
check(pin_ok, f"the {len(PINS)} pinned rows read as the witness's tests pin them")


def witness_grid():
    for J in [(1, 2, 3), (2, 2, 3), (2, -2, 3), (2, -2, 0), (3, 6, -2), (0, 2, 3), (1, 3, 0), (6, -3, 2), (2, 2, 0)]:
        for fd in [(0, 0, 0), (1, 0, 0), (2, 1, 0), (0, 0, 1), (1, 0, 2)]:
            nn = fd[0]**2 + fd[1]**2
            g = (1, 2, 1)
            if nn == 0:
                bu, bv = (1, 2, 1), (2, 1, -1)
            else:
                dg = g[0] * fd[0] + g[1] * fd[1]
                bu = (nn * g[0], nn * g[1], nn * g[2])
                bv = (2 * dg * fd[0] - nn * g[0], 2 * dg * fd[1] - nn * g[1], -nn * g[2])
            pairs = [(bu, bv), (bu, (bv[0] + 1, bv[1], bv[2])), (bu, (bv[0], bv[1], -bv[2])), (g, g), (g, (-1, -2, -1)),
                     ((1, 0, 1), (2, 1, 0)), ((1, 2, 1), (2, 1, 1)), ((1, 2, 3), (2, 1, -1)), ((0, 0, 1), (0, 0, -1)),
                     ((3, 1, 2), (-1, 2, 1)), ((1, 0, 2), (4, 0, 4))]
            for fu, fv in pairs:
                yield J, fd, fu, fv


wg = Counter()
for (J, fd, fu, fv) in witness_grid():
    near, far, exact = ends(J, fd, fu, fv)
    wg['rows'] += 1
    wg['palindromic'] += far > 0 and near == far
    k = rule(J, fd, fu, fv)
    wg['B'] += k == 'Boff'
    wg['S'] += (k is None and not fd[2] and (fd[0] or fd[1]) and J[0] == J[1] and fu == fv)
    wg['z on d'] += fd[2] != 0
    wg['disagree'] += (far > 0 and near == far) != (k is not None)
say(f"   the witness's grid off the plane, read here: {dict(wg)}")
check(wg['rows'] == 495 and wg['disagree'] == 0 and wg['palindromic'] == wg['B'],
      "the witness's grid: 495 rows, none disagreeing, the palindromic rows exactly those in B off the plane")

say()
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED")
import os
_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'triangle_z_part_gate.txt')
with open(_out, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('\n'.join(OUT) + '\n')
