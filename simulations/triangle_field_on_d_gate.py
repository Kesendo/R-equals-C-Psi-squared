"""Gate: the Heisenberg triangle with a planar field on the dephased seat (PROOF_ONE_DEPHASED_SEAT, section
"The triangle with a planar field on d").

Setting: sites (d, u, v) = (0, 1, 2), Heisenberg bonds s·σ_d·σ_u + t·σ_d·σ_v + r·σ_u·σ_v with real couplings, at most
one zero, the jump Z on d alone, planar fields f_d != 0, f_u, f_v (no z part). A turn about z keeps the bonds and the
jump and turns every planar field by one angle, so f_d = p·X_d with p != 0 below. Exact throughout (sympy over Q(i),
integer ranks), no threshold. Five parts.

1. The letter split. W anticommutes with Z_d exactly when W = X_d (x) alpha + Y_d (x) beta, and with Σ_P = s·P_u + t·P_v
   and H_uv = r·σ_u·σ_v + f_u·σ_u + f_v·σ_v,
      [H, W] = I_d (x) ([Σ_X, alpha] + [Σ_Y, beta]) + i·Z_d (x) ({Σ_X, beta} - {Σ_Y, alpha} + 2p·beta)
             + X_d (x) ([H_uv, alpha] - i{Σ_Z, beta}) + Y_d (x) ([H_uv, beta] + i{Σ_Z, alpha}),
   checked symbolically in all eight parameters on the 32 strings; the field on d enters the Z_d part alone. The turn
   about z is checked symbolically too.
2. The hop. The I_d and Z_d parts are the two hop blocks: with H_01 = p + 2(s·σ⁻_u + t·σ⁻_v), sigma⁻ = |1><0| (A in
   the code), they vanish exactly when M01 = alpha - i·beta satisfies H_01·H_01^dag·M01 = M01·H_01^dag·H_01 (and
   M10 = alpha + i·beta = H_01^-1 M01 H_01^dag), so the hop space is H_01 times the commutant of H_01^dag·H_01.
   det H_01 = p^4; the discriminant of H_01^dag·H_01's characteristic
   polynomial is 2^24·p^6·s^4·t^4·Q with Q a sum of terms each >= 0 and 2p^4(s^6 + t^6) > 0, so for real couplings
   and p != 0 the spectrum is simple exactly when s·t != 0, and the hop space is then four-dimensional; at s = 0 it is
   eight-dimensional (u a spectator of the hop). Bases for the strata s != +-t, s = t, s = -t, s = 0 are checked to lie
   in the hop space symbolically, each element with a string the others lack (so independent on the whole stratum),
   the bases at s = +-t as the limits W1' = +-W1/t^2 and W2' = t·(W1 + W2)/Δ, and the hop space's dimension is read
   exactly at rows of each.
3. The strata. The X_d and Y_d parts on each basis, symbolic; every row the hand elimination uses is pinned as the
   proof displays it, and the elimination's case logic is printed (it is argued in the proof, not gated here),
   and beside the hand steps the nullity of each small system is read at generic values and at the special points
   where a step's case turns.
   Result: far elements only in A (every field on f_d's line) for s != +-t; in A, B on the axis (s = t, f_v the mirror
   image of f_u in f_d's line) and S (s = t, f_u = f_v) for s = t; in A for s = -t; in A and D (s = 0, f_u ⊥ f_d,
   f_d + f_u + f_v = 0) for s = 0, and D' (u <-> v) for t = 0.
4. The elements, symbolic on their families: A's (f_d-hat·σ)^(x)3, B's 2·X_d (x) SWAP·X_uX_v, D's
   t·q·YYY + p·q·XYI + r·p·XXX (W^2 a positive multiple of 1, operator Schmidt ranks 2, 2, 3 across d|uv, u|dv,
   v|du), and S's X_d (x) P_s, singular, so S outside A does not pair by F158's Lemma 3 (its far end one-dimensional).
5. An exact census by a second route, the full commutator on the 64 Pauli strings with integer ranks (Bareiss, every
   division asserted exact), independent of the split: a grid in two directions of f_d and planted rows of every class
   with near misses, f_d in five directions, (near, far) tabulated per class against the rule
      far end nonzero <=> A, B on the axis, S, D or D';   pairs <=> A, B on the axis, D or D'.

Run: python simulations/triangle_field_on_d_gate.py   (about two and a half minutes; writes
simulations/results/triangle_field_on_d_gate.txt)
"""
import itertools
import random
import sys
from collections import Counter

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
# Pauli algebra on (u, v) and on (d, u, v)
I2 = sp.eye(2)
PX = sp.Matrix([[0, 1], [1, 0]])
PY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
PZ = sp.Matrix([[1, 0], [0, -1]])
P1 = {'I': I2, 'X': PX, 'Y': PY, 'Z': PZ}
LET = 'IXYZ'


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.kronecker_product(out, m)
    return out


STR2 = [a + b for a in LET for b in LET]
MAT2 = {x: kron(P1[x[0]], P1[x[1]]) for x in STR2}
STR3 = [''.join(x) for x in itertools.product(LET, repeat=3)]
MAT3 = {x: kron(P1[x[0]], P1[x[1]], P1[x[2]]) for x in STR3}


def two(a, b):
    return MAT2[a + b]


def coeffs(M, strs, mats, dim):
    return [sp.expand((mats[x] * M).trace() / dim) for x in strs]


def zero(M):
    return all(sp.expand(e) == 0 for e in M)


s, t, r, p = sp.symbols('s t r p', real=True)
ux, uy, vx, vy = sp.symbols('u_x u_y v_x v_y', real=True)


def H3(s_, t_, r_, fd, fu, fv):
    def op(l, site):
        ls = ['I', 'I', 'I']
        ls[site] = l
        return MAT3[''.join(ls)]
    H = sp.zeros(8)
    for (i, j, J) in ((0, 1, s_), (0, 2, t_), (1, 2, r_)):
        for l in 'XYZ':
            H += J * op(l, i) * op(l, j)
    for site, f in ((0, fd), (1, fu), (2, fv)):
        H += f[0] * op('X', site) + f[1] * op('Y', site)
    return H


def letter_parts(s_, t_, r_, p_, fu, fv, alpha, beta):
    SX = s_ * two('X', 'I') + t_ * two('I', 'X')
    SY = s_ * two('Y', 'I') + t_ * two('I', 'Y')
    SZ = s_ * two('Z', 'I') + t_ * two('I', 'Z')
    Huv = (r_ * (two('X', 'X') + two('Y', 'Y') + two('Z', 'Z'))
           + fu[0] * two('X', 'I') + fu[1] * two('Y', 'I') + fv[0] * two('I', 'X') + fv[1] * two('I', 'Y'))
    EI = SX * alpha - alpha * SX + SY * beta - beta * SY
    EZ = SX * beta + beta * SX - (SY * alpha + alpha * SY) + 2 * p_ * beta
    EX = Huv * alpha - alpha * Huv - sp.I * (SZ * beta + beta * SZ)
    EY = Huv * beta - beta * Huv + sp.I * (SZ * alpha + alpha * SZ)
    return EI, EZ, EX, EY


# ---------------------------------------------------------------------------------------------------------------
say("1. The letter split, symbolic in s, t, r, p, f_u, f_v")
H = H3(s, t, r, (p, 0), (ux, uy), (vx, vy))
ok = True
for x in STR2:
    for which in 'XY':
        alpha = MAT2[x] if which == 'X' else sp.zeros(4)
        beta = MAT2[x] if which == 'Y' else sp.zeros(4)
        W = kron(PX, alpha) + kron(PY, beta)
        EI, EZ, EX, EY = letter_parts(s, t, r, p, (ux, uy), (vx, vy), alpha, beta)
        rhs = kron(I2, EI) + sp.I * kron(PZ, EZ) + kron(PX, EX) + kron(PY, EY)
        if not zero(H * W - W * H - rhs):
            ok = False
check(ok, "[H, X_d(x)alpha + Y_d(x)beta] = I_d(x)E_I + i Z_d(x)E_Z + X_d(x)E_X + Y_d(x)E_Y on all 32 strings")
far_strings = [x for x in STR3 if x[0] in 'XY']
check(all(zero(MAT3['ZII'] * MAT3[x] + MAT3[x] * MAT3['ZII']) for x in far_strings)
      and len(far_strings) == 32, "the 32 strings with X or Y on d are exactly the strings anticommuting with Z_d")
# the turn about z: U = (c - i sn Z)^(x)3, c^2 + sn^2 = 1, carries every planar field by the angle theta
mm = sp.symbols('m', real=True)
c, sn = (1 - mm**2) / (1 + mm**2), 2 * mm / (1 + mm**2)       # cos and sin of theta/2, c^2 + sn^2 = 1
U1 = c * I2 - sp.I * sn * PZ
U = kron(U1, U1, U1)
cth, sth = c**2 - sn**2, 2 * c * sn


def rot(f):
    return (cth * f[0] - sth * f[1], sth * f[0] + cth * f[1])


lhs = U * H * U.H
rhs = H3(s, t, r, rot((p, 0)), rot((ux, uy)), rot((vx, vy)))
check(all(sp.simplify(e) == 0 for e in (lhs - rhs)), "a turn about z keeps the bonds and turns every planar field by one angle")
check(all(sp.simplify(e) == 0 for e in (U * MAT3['ZII'] * U.H - MAT3['ZII'])),
      "the turn keeps the jump Z_d, so both ends map to themselves")

# ---------------------------------------------------------------------------------------------------------------
say()
say("2. The hop: M01 = H_01·Ω with Ω in the commutant of H_01^dag·H_01; when that commutant is four-dimensional")
sm = sp.Matrix([[0, 0], [1, 0]])
A = p * sp.eye(4) + 2 * (s * kron(sm, I2) + t * kron(I2, sm))
check(sp.expand(A.det() - p**4) == 0, "det H_01 = p^4, so the hop block is invertible exactly when p != 0")
Hd = H3(s, t, r, (p, 0), (ux, uy), (vx, vy))
check(zero(Hd[0:4, 4:8] - A), "H_01 = p + 2(s σ⁻_u + t σ⁻_v) is H's block from d = 1 to d = 0 (basis order d, u, v; |0> = up)")
AdA = A.H * A
check(coeffs(AdA, STR2, MAT2, 4) == [sp.expand(x) for x in
      [p**2 + 2 * s**2 + 2 * t**2, 2 * p * t, 0, 2 * t**2, 2 * p * s, 2 * s * t, 0, 0, 0, 0, 2 * s * t, 0, 2 * s**2, 0, 0, 0]],
      "H_01^dag·H_01 = p^2 + 2(s^2 + t^2) + 2p(s X_u + t X_v) + 2s^2 Z_u + 2t^2 Z_v + 2st(X_uX_v + Y_uY_v)")
lam = sp.symbols('lam')
disc = sp.discriminant(sp.Poly((AdA - lam * sp.eye(4)).det(), lam).as_expr(), lam)
Q = (p**6 * (s**2 - t**2)**2 + 2 * p**4 * (s**6 + t**6)
     + p**2 * ((s**4 - t**4)**2 + 10 * s**2 * t**2 * (s**2 - t**2)**2 + 13 * s**4 * t**4)
     + 4 * s**2 * t**2 * (s**2 + t**2)**3)
check(sp.expand(disc - 2**24 * p**6 * s**4 * t**4 * Q) == 0,
      "disc = 2^24 p^6 s^4 t^4 Q, Q = p^6(s^2-t^2)^2 + 2p^4(s^6+t^6) + p^2[(s^4-t^4)^2 + 10s^2t^2(s^2-t^2)^2 + 13s^4t^4] + 4s^2t^2(s^2+t^2)^3")
say("   (each term of Q is >= 0 for real s, t, p, and 2p^4(s^6 + t^6) > 0 when p != 0 and (s, t) != 0: the spectrum of")
say("    the Hermitian H_01^dag·H_01 is simple, its commutant four-dimensional, exactly when p·s·t != 0; read, not gated)")
M_v = sp.Matrix([[p**2 + 4 * t**2, 2 * t * p], [2 * t * p, p**2]])
check(zero(AdA.subs(s, 0) - kron(I2, M_v)) and sp.expand(sp.discriminant((M_v - lam * sp.eye(2)).det(), lam) - 16 * t**2 * (p**2 + t**2)) == 0,
      "at s = 0, H_01^dag·H_01 = 1_u (x) M_v with disc(M_v) = 16t^2(p^2 + t^2) > 0: two simple eigenvalues on v")

# the hop equations in letters
def hop_system(s_, t_, p_):
    SX = s_ * two('X', 'I') + t_ * two('I', 'X')
    SY = s_ * two('Y', 'I') + t_ * two('I', 'Y')
    cols = []
    for which in 'ab':
        for x in STR2:
            a = MAT2[x] if which == 'a' else sp.zeros(4)
            b = MAT2[x] if which == 'b' else sp.zeros(4)
            EI = SX * a - a * SX + SY * b - b * SY
            EZ = SX * b + b * SX - (SY * a + a * SY) + 2 * p_ * b
            cols.append(coeffs(EI, STR2, MAT2, 4) + coeffs(EZ, STR2, MAT2, 4))
    return sp.Matrix(cols).T


def in_hop(s_, t_, p_, alpha, beta):
    SX = s_ * two('X', 'I') + t_ * two('I', 'X')
    SY = s_ * two('Y', 'I') + t_ * two('I', 'Y')
    EI = SX * alpha - alpha * SX + SY * beta - beta * SY
    EZ = SX * beta + beta * SX - (SY * alpha + alpha * SY) + 2 * p_ * beta
    return zero(EI) and zero(EZ)


# the closed form of the hop space: M01 = A·(A^dag A)^k, k = 0..3 (A is H_01), split into alpha and beta
def split(d0, d1):
    return (d0 + d1) / 2, (d1 - d0) / (2 * sp.I)      # d0 = M01 = alpha - i beta, d1 = M10 = alpha + i beta


okc = True
for k in range(4):
    d0 = A * AdA**k
    d1 = A.inv() * d0 * A.H
    al, be = split(d0, d1)
    okc = okc and in_hop(s, t, p, al.applyfunc(sp.simplify), be.applyfunc(sp.simplify))
check(okc, "M01 = H_01·(H_01^dag H_01)^k (k = 0..3) with M10 = H_01^-1 M01 H_01^dag solve both hop parts, symbolically")
check(len(hop_system(s, t, p).nullspace()) == 4, "the hop parts have a four-dimensional solution space at generic s, t, p")

D = s**2 - t**2
BASES = {
    'generic': (s, t, {
        'W0': (two('X', 'X'), sp.zeros(4)),
        'WH': (p * sp.eye(4) + s * two('X', 'I') + t * two('I', 'X'), s * two('Y', 'I') + t * two('I', 'Y')),
        'W1': (p * D / t * two('X', 'I') - s * t * sp.eye(4) + s**2 * two('Y', 'Y') + s * t * two('Z', 'Z'), D * two('X', 'Y')),
        'W2': (p * D / s * two('I', 'X') + s * t * sp.eye(4) - t**2 * two('Y', 'Y') - s * t * two('Z', 'Z'), D * two('Y', 'X')),
    }),
    's=t': (t, t, {
        'W0': (two('X', 'X'), sp.zeros(4)),
        "W1'": (-sp.eye(4) + two('Y', 'Y') + two('Z', 'Z'), sp.zeros(4)),
        'WH': (p * sp.eye(4) + t * two('X', 'I') + t * two('I', 'X'), t * two('Y', 'I') + t * two('I', 'Y')),
        "W2'": (p * two('X', 'I') + p * two('I', 'X') + t * two('Y', 'Y'), t * two('X', 'Y') + t * two('Y', 'X')),
    }),
    's=-t': (-t, t, {
        'W0': (two('X', 'X'), sp.zeros(4)),
        "W1'": (-sp.eye(4) - two('Y', 'Y') + two('Z', 'Z'), sp.zeros(4)),
        '-WH': (-p * sp.eye(4) + t * two('X', 'I') - t * two('I', 'X'), t * two('Y', 'I') - t * two('I', 'Y')),
        "W2'": (p * two('X', 'I') - p * two('I', 'X') + t * two('Y', 'Y'), t * two('X', 'Y') + t * two('Y', 'X')),
    }),
    's=0': (0, t, dict(
        [(f'e{k}', (two(O, 'X'), sp.zeros(4))) for k, O in enumerate('IXYZ')]
        + [(f'e{4 + k}', (p * two(O, 'I'), t * two(O, 'Y'))) for k, O in enumerate('IXYZ')])),
}
for name, (sv, tv, basis) in BASES.items():
    inh = all(in_hop(sv, tv, p, a, b) for (a, b) in basis.values())
    vecs = sp.Matrix([coeffs(a, STR2, MAT2, 4) + coeffs(b, STR2, MAT2, 4) for (a, b) in basis.values()])
    check(inh, f"stratum {name}: its {len(basis)} basis elements solve the hop parts symbolically")
    rk = vecs.rank()
    check(rk == len(basis), f"stratum {name}: they are independent (rank {rk}) at generic values")
    priv = []
    for i in range(len(basis)):
        own = [j for j in range(32) if vecs[i, j] != 0 and all(vecs[k, j] == 0 for k in range(len(basis)) if k != i)]
        priv.append(own)
    names_ = list(basis)
    labels = ['a.' + x for x in STR2] + ['b.' + x for x in STR2]
    nz = {'generic': [s, t, s - t, s + t], 's=t': [t], 's=-t': [t], 's=0': [t, p]}[name]
    def nonzero_there(e):
        fs = [f for f, _ in sp.factor_list(e)[1]]
        return all(any(sp.simplify(f - g) == 0 or sp.simplify(f + g) == 0 for g in nz + [p]) for f in fs)
    okp = all(any(nonzero_there(vecs[i, j]) for j in priv[i]) for i in range(len(basis)))
    check(okp, f"stratum {name}: each element has a string the others lack, with a coefficient nonzero on the whole stratum "
               + "(" + ", ".join(f"{names_[i]}: {labels[priv[i][0]]}" for i in range(len(basis))) + ")")

# the bases at s = +-t as limits of the generic one
W1a, W1b = BASES['generic'][2]['W1']
W2a, W2b = BASES['generic'][2]['W2']
Ea = (p / t) * two('X', 'I') + (p / s) * two('I', 'X') + two('Y', 'Y')
Eb = two('X', 'Y') + two('Y', 'X')
check(zero((W1a + W2a - D * Ea).applyfunc(sp.simplify)) and zero((W1b + W2b - D * Eb).applyfunc(sp.simplify)),
      "W1 + W2 = Δ·(X_d (x) [(p/t)X_u + (p/s)X_v + Y_uY_v] + Y_d (x) (X_uY_v + Y_uX_v)), no Δ left in the bracket")
okl = True
for name, sv, sign in (('s=t', t, 1), ('s=-t', -t, -1)):
    b = BASES[name][2]
    okl &= zero((W1a.subs(s, sv) - sign * t**2 * b["W1'"][0]).applyfunc(sp.simplify)) and zero(W1b.subs(s, sv))
    okl &= zero((t * Ea.subs(s, sv) - b["W2'"][0]).applyfunc(sp.simplify)) and zero((t * Eb - b["W2'"][1]).applyfunc(sp.simplify))
check(okl, "at s = t: W1 = t^2·W1', and at s = -t: W1 = -t^2·W1'; at both, W2' = t·(W1 + W2)/Δ")

# exact dimension of the hop space at rows of each stratum (a nonzero p, s, t; rational, signs mixed)
rows_hop = {'generic': [(1, 2, 1), (-3, 2, 5), (2, -1, -1), (1, 3, 7)],
            's=t': [(2, 2, 1), (-1, -1, 3), (3, 3, -2)],
            's=-t': [(2, -2, 1), (-1, 1, 3)],
            's=0': [(0, 2, 1), (0, -3, 5)]}
for name, rws in rows_hop.items():
    dims = [32 - hop_system(sp.Integer(a), sp.Integer(b), sp.Integer(c_)).rank() for (a, b, c_) in rws]
    want = 8 if name == 's=0' else 4
    check(all(d_ == want for d_ in dims), f"stratum {name}: hop space of dimension {want} at {rws} (read {dims})")
    if name != 's=0':
        sv, tv, basis = BASES[name]
        okb = True
        for (a, b, c_) in rws:
            subs = {s: a, t: b, p: c_}
            vecs = sp.Matrix([[sp.expand(x.subs(subs)) for x in coeffs(al, STR2, MAT2, 4) + coeffs(be, STR2, MAT2, 4)]
                              for (al, be) in basis.values()])
            okb = okb and vecs.rank() == 4
        check(okb, f"stratum {name}: the basis is independent at those rows, so it spans the hop space there")

# ---------------------------------------------------------------------------------------------------------------
say()
say("3. The strata: the X_d and Y_d parts on each basis, and the hand elimination pinned to its rows")


def table(name):
    sv, tv, basis = BASES[name]
    rows = {}
    for k, (key, (a, b)) in enumerate(basis.items()):
        _, _, EX, EY = letter_parts(sv, tv, r, p, (ux, uy), (vx, vy), a, b)
        for lab, E in (('X_d', EX), ('Y_d', EY)):
            for x, cf in zip(STR2, coeffs(E, STR2, MAT2, 4)):
                if cf != 0:
                    rows.setdefault(f"{lab}.{x}", {})[key] = sp.factor(cf / (2 * sp.I))
    return rows


def show(rows):
    for key in sorted(rows):
        say(f"     {key}: " + ", ".join(f"{k}: {v}" for k, v in rows[key].items()))


def row_is(rows, key, want):
    got = rows.get(key, {})
    return set(got) == set(want) and all(sp.expand(got[k] - want[k]) == 0 for k in want)


# generic stratum
Tg = table('generic')
say("   stratum s != +-t, s·t != 0 (unknowns c0, cH, c1, c2 on W0 = XXX, WH = H's flip part, W1, W2):")
show(Tg)
pins = {
    'Y_d.ZY': {'W1': -uy * D}, 'Y_d.YZ': {'W2': -vy * D},
    'X_d.IZ': {'WH': -t * vy, 'W2': -p * vy * D / s}, 'X_d.ZI': {'WH': -s * uy, 'W1': -p * uy * D / t},
    'X_d.XZ': {'W0': -vy, 'W1': s * t * uy, 'W2': -s * t * uy}, 'X_d.ZX': {'W0': -uy, 'W1': s * t * vy, 'W2': -s * t * vy},
}
check(all(row_is(Tg, k, v) for k, v in pins.items()), "the six rows the elimination uses read as the proof displays them")
say("   hand elimination: Y_d.ZY gives u_y·c1 = 0, Y_d.YZ gives v_y·c2 = 0 (s^2 != t^2);")
say("     u_y, v_y != 0: c1 = c2 = 0, X_d.IZ gives cH = 0, X_d.XZ gives c0 = 0;")
say("     u_y != 0 = v_y: c1 = 0, X_d.ZI gives cH = 0, X_d.ZX gives c0 = 0, X_d.XZ gives c2 = 0;")
say("     u_y = 0 != v_y: c2 = 0, X_d.IZ gives cH = 0, X_d.XZ gives c0 = 0, X_d.ZX gives c1 = 0;")
say("     u_y = v_y = 0 is A (every field on f_d's line), with XXX among its far elements.")


def nullity(rows, keys, subs):
    M = sp.Matrix([[sp.expand(sp.S(rows[rk].get(k, 0)).subs(subs)) for k in keys] for rk in rows])
    return len(keys) - M.rank()


kg = ['W0', 'WH', 'W1', 'W2']
check(nullity(Tg, kg, {}) == 0 and nullity(Tg, kg, {uy: 0}) == 0 and nullity(Tg, kg, {vy: 0}) == 0
      and nullity(Tg, kg, {uy: 0, vy: 0}) >= 1,
      "beside it, the generic nullity of the 4-column system: 0, 0, 0 off A, at least 1 on A")

# s = t
Tt = table('s=t')
say("   stratum s = t != 0 (unknowns c0, c1, cH, c2 on W0 = XXX, W1' = X_d(x)(YY + ZZ - 1), WH, W2'):")
show(Tt)
pins_t = {
    'Y_d.IZ': {'WH': t * (p + vx)}, 'Y_d.ZI': {'WH': t * (p + ux)},
    'Y_d.XZ': {'WH': t**2, "W2'": t * (p + vx)}, 'Y_d.ZX': {'WH': t**2, "W2'": t * (p + ux)},
    'Y_d.YZ': {"W2'": -t * vy}, 'Y_d.ZY': {"W2'": -t * uy},
    'X_d.XZ': {'W0': -vy, "W1'": uy}, 'X_d.ZX': {'W0': -uy, "W1'": vy},
    'X_d.YZ': {"W1'": vx - ux, 'WH': -t**2, "W2'": t * vx}, 'X_d.ZY': {"W1'": ux - vx, 'WH': -t**2, "W2'": t * ux},
    'X_d.IZ': {'WH': -t * vy, "W2'": -p * vy}, 'X_d.ZI': {'WH': -t * uy, "W2'": -p * uy},
}
check(all(row_is(Tt, k, v) for k, v in pins_t.items()) and set(Tt) == set(pins_t),
      "the twelve rows read as the proof displays them, and there are no others")
say("   hand elimination: cH != 0 needs u_x = v_x = -p (Y_d.IZ, Y_d.ZI), and then Y_d.XZ reads t^2·cH = 0: cH = 0;")
say("     c2 != 0 then needs u_y = v_y = 0 (Y_d.YZ, Y_d.ZY) and u_x = v_x = -p (Y_d.XZ, Y_d.ZX), and X_d.YZ reads -t·p·c2 = 0: c2 = 0;")
say("     left: (v_x - u_x)·c1 = 0 and [[-v_y, u_y], [-u_y, v_y]]·(c0, c1) = 0, determinant u_y^2 - v_y^2:")
say("     c1 = 0 needs u_y = v_y = 0 (A); c1 != 0 needs u_x = v_x and u_y = v_y (S, c0 = c1) or u_y = -v_y (B on the axis, c0 = -c1).")
e0a, e1a = BASES['s=t'][2]['W0'][0], BASES['s=t'][2]["W1'"][0]
Ps = (sp.eye(4) - (two('X', 'X') + two('Y', 'Y') + two('Z', 'Z'))) / 4
SWAP = (sp.eye(4) + two('X', 'X') + two('Y', 'Y') + two('Z', 'Z')) / 2
check(zero(e0a + e1a + 4 * Ps), "c0 = c1: X_d (x) (XX + YY + ZZ - 1) = -4·X_d (x) P_s, the singlet element of S")
check(zero(e0a - e1a - 2 * SWAP * two('X', 'X')), "c0 = -c1: X_d (x) (1 + XX - YY - ZZ) = 2·X_d (x) SWAP·X_uX_v, B's element on the axis x")
kt = list(BASES['s=t'][2])
check(nullity(Tt, kt, {}) == 0 and nullity(Tt, kt, {vx: ux, vy: uy}) == 1 and nullity(Tt, kt, {vx: ux, vy: -uy}) == 1
      and nullity(Tt, kt, {uy: 0, vy: 0}) >= 1 and nullity(Tt, kt, {vx: ux + 1, vy: uy}) == 0,
      "beside it, generic nullities: 0 off the classes, 1 on S, 1 on B's axis, >= 1 on A, 0 at f_v = f_u + x")
check(nullity(Tt, kt, {ux: -p, vx: -p}) == 0 and nullity(Tt, kt, {ux: -p, vx: -p, uy: 0, vy: 0}) >= 1
      and nullity(Tt, kt, {ux: -p, vx: -p, vy: -uy}) == 1 and nullity(Tt, kt, {ux: -p, vx: -p, vy: uy}) == 1,
      "at the turning point u_x = v_x = -p: 0 off the classes, 1 on S and on B's axis, >= 1 on A")

# s = -t
Tm = table('s=-t')
say("   stratum s = -t != 0:")
show(Tm)
pins_m = {
    'Y_d.IZ': {"W1'": -2 * t, '-WH': -t * (p + vx)}, 'Y_d.ZI': {"W1'": 2 * t, '-WH': t * (p + ux)},
    'Y_d.XZ': {'-WH': -t * (2 * r - t), "W2'": t * (p + vx)}, 'Y_d.ZX': {'-WH': t * (2 * r + t), "W2'": t * (p + ux)},
    'Y_d.YZ': {"W2'": -t * vy}, 'Y_d.ZY': {"W2'": -t * uy},
    'X_d.XZ': {'W0': -vy, "W1'": uy}, 'X_d.ZX': {'W0': -uy, "W1'": vy},
}
check(all(row_is(Tm, k, v) for k, v in pins_m.items()), "the rows the elimination uses read as the proof displays them")
say("   hand elimination: c2 != 0 needs u_y = v_y = 0 (Y_d.YZ, Y_d.ZY), class A; with c2 = 0, Y_d.XZ and Y_d.ZX read")
say("     (2r - t)·cH = 0 = (2r + t)·cH, not both coefficients zero at t != 0: cH = 0; Y_d.IZ then gives c1 = 0, and")
say("     X_d.XZ, X_d.ZX give v_y·c0 = u_y·c0 = 0: outside A no far element.")
km = list(BASES['s=-t'][2])
check(nullity(Tm, km, {}) == 0 and nullity(Tm, km, {vx: ux, vy: -uy}) == 0 and nullity(Tm, km, {vx: ux, vy: uy}) == 0
      and nullity(Tm, km, {uy: 0, vy: 0}) >= 1,
      "beside it, generic nullities: 0 off A (B's and S's field patterns included), >= 1 on A")
check(all(nullity(Tm, km, {**sub, r: rr}) == 0 for rr in (t / 2, -t / 2)
          for sub in ({}, {vx: ux, vy: -uy}, {vx: ux, vy: uy}, {ux: -p, vx: -p})),
      "at the turning points 2r = +-t: 0 off A, with B's and S's field patterns and u_x = v_x = -p among them")

# s = 0
T0 = table('s=0')
say("   stratum s = 0, t, r != 0 (unknowns c0..c7; e_k = X_d O_u X_v and e_(4+k) = p X_d O_u + t Y_d O_u Y_v, O = I, X, Y, Z):")
show(T0)
pins_0 = {
    'Y_d.IX': {'e7': -r * t}, 'X_d.IY': {'e3': r}, 'Y_d.ZI': {'e5': -r * t}, 'Y_d.ZX': {'e4': -r * t},
    'X_d.ZY': {'e0': r, 'e5': -p * r}, 'X_d.IZ': {'e0': -vy, 'e2': -r},
    'X_d.ZI': {'e2': r, 'e5': -p * uy, 'e6': p * ux}, 'Y_d.YZ': {'e6': t * (p + vx)}, 'Y_d.ZY': {'e5': -t * uy, 'e6': t * ux},
    'X_d.XZ': {'e1': -vy, 'e6': -p * r}, 'X_d.ZX': {'e1': -uy, 'e2': ux, 'e6': p * r},
}
check(all(row_is(T0, k, v) for k, v in pins_0.items()), "the rows the elimination uses read as the proof displays them")
say("   hand elimination: Y_d.IX, X_d.IY, Y_d.ZI, Y_d.ZX give c7 = c3 = c5 = c4 = 0; X_d.ZY then c0 = 0, X_d.IZ c2 = 0;")
say("     c6 != 0 needs u_x = 0 (X_d.ZI) and v_x = -p (Y_d.YZ); X_d.XZ and X_d.ZX read v_y·c1 = -p·r·c6 and")
say("     u_y·c1 = p·r·c6, so v_y = -u_y != 0 and c1 = p·r·c6/u_y: family D (f_u ⊥ f_d, f_d + f_u + f_v = 0);")
say("     c6 = 0 leaves u_y·c1 = v_y·c1 = 0 and c1 on XXX: class A.")
q = sp.symbols('q', real=True)
k0 = [f'e{k}' for k in range(8)]
check(nullity(T0, k0, {}) == 0 and nullity(T0, k0, {ux: 0, vx: -p, vy: -uy}) == 1
      and nullity(T0, k0, {ux: 0, vx: -p, vy: -uy + 1}) == 0 and nullity(T0, k0, {uy: 0, vy: 0}) >= 1,
      "beside it, generic nullities: 0 off the classes, 1 on D, 0 one unit off D, >= 1 on A")
check(nullity(T0, k0, {vx: -p, vy: -uy}) == 0 and nullity(T0, k0, {ux: 0, vy: -uy}) == 0
      and nullity(T0, k0, {ux: 0, vx: -p}) == 0,
      "with one of D's three conditions left free (u_x, v_x or v_y generic): 0")

# ---------------------------------------------------------------------------------------------------------------
say()
say("4. The elements on their families, symbolic")
Zd = MAT3['ZII']


def far_ok(W, Hm):
    return zero(Hm * W - W * Hm) and zero(Zd * W + W * Zd)


a_, b_ = sp.symbols('a b', real=True)
HA = H3(s, t, r, (p, 0), (a_, 0), (b_, 0))
check(far_ok(MAT3['XXX'], HA), "A: XXX commutes with H when every field lies along x, and anticommutes with Z_d")
HB = H3(s, s, r, (p, 0), (ux, uy), (ux, -uy))
WB = kron(PX, 2 * SWAP * two('X', 'X'))
check(far_ok(WB, HB) and zero(WB * WB - 4 * sp.eye(8)),
      "B on the axis: 2·X_d (x) SWAP·X_uX_v is far at s = t, f_v = (u_x, -u_y), and squares to 4·1")
HD = H3(0, t, r, (p, 0), (0, q), (-p, -q))
WD = t * q * MAT3['YYY'] + p * q * MAT3['XYI'] + r * p * MAT3['XXX']
check(far_ok(WD, HD) and zero(WD * WD - (t**2 * q**2 + r**2 * p**2 + p**2 * q**2) * sp.eye(8)),
      "D: t·q·YYY + p·q·XYI + r·p·XXX is far at s = 0, f_u = (0, q), f_v = -(f_d + f_u), W^2 = (t^2q^2 + r^2p^2 + p^2q^2)·1")
cD = {x: sp.expand((MAT3[x] * WD).trace() / 8) for x in STR3}


def schmidt(cut):
    rows_, cols_ = LET, [a + b for a in LET for b in LET]
    def coef(a, rest):
        st = ['', '', '']
        st[cut] = a
        others = [k for k in range(3) if k != cut]
        st[others[0]], st[others[1]] = rest[0], rest[1]
        return cD[''.join(st)]
    return sp.Matrix([[coef(a, rest) for rest in cols_] for a in rows_]).rank()


check([schmidt(0), schmidt(1), schmidt(2)] == [2, 2, 3],
      "D's element has operator Schmidt ranks 2, 2, 3 across d|uv, u|dv, v|du at generic p, q, r, t")
HS = H3(s, s, r, (p, 0), (ux, uy), (ux, uy))
WS = kron(PX, Ps)
check(far_ok(WS, HS) and WS.rank() == 2, "S: X_d (x) P_s is far at s = t, f_u = f_v, any f_d along x, and has rank 2 of 8")
say("   On S outside A the far end is spanned by that singular element (stratum s = t, c0 = c1), so it holds no")
say("   invertible element and the row does not pair (F158's Lemma 3); the census reads its near end as 2.")
check(nullity(Tt, kt, {vx: ux, vy: uy}) == 1, "S outside A: the far end is one-dimensional (generic u_x, u_y)")

# ---------------------------------------------------------------------------------------------------------------
say()
say("5. The census: the full commutator on Pauli strings, integer ranks, against the rule")


def _mul1(a, b):
    if a == 0:
        return 0, b
    if b == 0:
        return 0, a
    if a == b:
        return 0, 0
    return (1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - a - b


S3 = list(itertools.product(range(4), repeat=3))
IDX = {x: k for k, x in enumerate(S3)}
NEAR = [x for x in S3 if x[0] in (0, 3)]
FAR = [x for x in S3 if x[0] in (1, 2)]


def comm_cols(terms, basis):
    cols = []
    for b in basis:
        col = [0] * 64
        for (a, h) in terms:
            ph, out = 0, []
            for x, y in zip(a, b):
                qq, cc = _mul1(x, y)
                ph += qq
                out.append(cc)
            if ph % 2 == 1:                       # anticommuting: [P_a, P_b]/(2i) = i^(ph - 1) P_c
                col[IDX[tuple(out)]] += (1 if (ph - 1) % 4 == 0 else -1) * h
        cols.append(col)
    return cols


def rank_int(cols):
    M = [list(x) for x in cols]
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


def terms_of(s_, t_, r_, fd, fu, fv):
    T = []
    for (i, j, J) in ((0, 1, s_), (0, 2, t_), (1, 2, r_)):
        if J:
            for l in (1, 2, 3):
                a = [0, 0, 0]
                a[i] = a[j] = l
                T.append((tuple(a), J))
    for site, f in ((0, fd), (1, fu), (2, fv)):
        for l, comp in ((1, f[0]), (2, f[1])):
            if comp:
                a = [0, 0, 0]
                a[site] = l
                T.append((tuple(a), comp))
    return T


def ends(s_, t_, r_, fd, fu, fv):
    T = terms_of(s_, t_, r_, fd, fu, fv)
    return 32 - rank_int(comm_cols(T, NEAR)), 32 - rank_int(comm_cols(T, FAR))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def classes(s_, t_, r_, fd, fu, fv):
    k = ''
    if cross(fd, fu) == 0 and cross(fd, fv) == 0:
        k += 'A'
    nn = dot(fd, fd)
    if s_ == t_ and (2 * dot(fu, fd) * fd[0] - nn * fu[0], 2 * dot(fu, fd) * fd[1] - nn * fu[1]) == (nn * fv[0], nn * fv[1]):
        k += 'B'
    if s_ == t_ and fu == fv:
        k += 'S'
    if s_ == 0 and t_ and r_ and dot(fd, fu) == 0 and fv == (-(fd[0] + fu[0]), -(fd[1] + fu[1])):
        k += 'D'
    if t_ == 0 and s_ and r_ and dot(fd, fv) == 0 and fu == (-(fd[0] + fv[0]), -(fd[1] + fv[1])):
        k += "D'"
    return k or '-'


tab = {}
bad = []
seen = set()


def read(label, s_, t_, r_, fd, fu, fv):
    if [s_, t_, r_].count(0) > 1:
        return
    near, far = ends(s_, t_, r_, fd, fu, fv)
    seen.add((s_, t_, r_, fd, fu, fv))
    k = classes(s_, t_, r_, fd, fu, fv)
    tab.setdefault((label, k), Counter())[(near, far)] += 1
    if (far > 0) != (k != '-') or (near == far) != (k not in ('-', 'S')):
        bad.append((label, (s_, t_, r_), fd, fu, fv, near, far, k))


Jset = [-1, 1, 2, 0]
Fset = [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 1), (2, 1), (1, -1), (-1, 0), (0, -1), (1, 2), (-1, -1), (2, -1)]
for (s_, t_, r_) in itertools.product(Jset, repeat=3):
    for fd in ((1, 0), (2, 1)):
        for fu in Fset:
            for fv in Fset:
                read('grid', s_, t_, r_, fd, fu, fv)


def perp(v):
    return (-v[1], v[0])


random.seed(7)
JJ = [-3, -2, -1, 1, 2, 3, 5]
for fd in [(1, 0), (2, 1), (-1, 3), (3, -2), (0, 1)]:
    for _ in range(40):
        t_, r_ = random.choice(JJ), random.choice(JJ)
        k = random.choice([-3, -2, -1, 1, 2, 3])
        fu = (k * perp(fd)[0], k * perp(fd)[1])
        fv = (-(fd[0] + fu[0]), -(fd[1] + fu[1]))
        read('D', 0, t_, r_, fd, fu, fv)
        read("D'", t_, 0, r_, fd, fv, fu)
        read('D, f_v off by x', 0, t_, r_, fd, fu, (fv[0] + 1, fv[1]))
        read('D, f_u off perp', 0, t_, r_, fd, (fu[0] + fd[0], fu[1] + fd[1]), (-(2 * fd[0] + fu[0]), -(2 * fd[1] + fu[1])))
        read('D, s = 1', 1, t_, r_, fd, fu, fv)
        s_ = random.choice(JJ)
        g = (random.randint(-3, 3), random.randint(-3, 3))
        nn = dot(fd, fd)
        fuB = (nn * g[0], nn * g[1])
        fvB = (2 * dot(g, fd) * fd[0] - nn * g[0], 2 * dot(g, fd) * fd[1] - nn * g[1])
        read('B', s_, s_, r_, fd, fuB, fvB)
        read('B, t = s + 1', s_, s_ + 1, r_, fd, fuB, fvB)
        read('B, f_d turned', s_, s_, r_, perp(fd), fuB, fvB)
        read('S', s_, s_, r_, fd, g, g)
        read('S, f_v off by x', s_, s_, r_, fd, g, (g[0] + 1, g[1]))
        a1, b1 = random.randint(-3, 3), random.randint(-3, 3)
        s3, t3 = random.choice(JJ), random.choice(JJ)
        read('A', s3, t3, r_, fd, (a1 * fd[0], a1 * fd[1]), (b1 * fd[0], b1 * fd[1]))
        read('A, f_u off the line', s3, t3, r_, fd, (a1 * fd[0] - fd[1], a1 * fd[1] + fd[0]), (b1 * fd[0], b1 * fd[1]))
        read('s = -t', s3, -s3, r_, fd, g, (b1, a1))

n_rows = sum(sum(cn.values()) for cn in tab.values())
say(f"   {n_rows} reads of {len(seen)} distinct rows; (near, far) by planted label and class (A, B on the axis, S, D, D', '-' for none):")
for (lab, k) in sorted(tab):
    say(f"     {lab:20s} {k:5s} {dict(sorted(tab[(lab, k)].items()))}")
check(not bad, f"every row: far end nonzero exactly in A, B, S, D, D'; pairs exactly in A, B, D, D' ({len(bad)} disagree)")
for b in bad[:10]:
    say(f"     disagrees: {b}")
alone = {k for (_, k) in tab}
check({'A', 'B', 'S', 'D', "D'", '-'} <= alone, "each class is met alone, outside the others")
check(all(nf == (2, 1) for (lab, k), cn in tab.items() if k == 'S' for nf in cn),
      "every row of S outside A, B and D reads near 2, far 1")
check(all(nf == (1, 1) for (lab, k), cn in tab.items() if k in ('B', 'D', "D'") for nf in cn),
      "every row of B, D or D' alone reads near 1, far 1")

say()
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED")
import os
_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'triangle_field_on_d_gate.txt')
with open(_out, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('\n'.join(OUT) + '\n')
raise SystemExit(1 if FAILS else 0)
