"""One nucleus, two bases: the rate relations of one spin-3/2 sample, and one published data set read against them
(experiments/ONE_NUCLEUS_TWO_BASES.md).

The generator is the one of experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md: jumps ZZ at 3J0/2, ZX and ZY at 3J1/2,
XI and YI at 3J2/2, each g(P rho P - rho), splitting q*ZZ, in Souza et al.'s labelling (|3/2>, |1/2>, |-1/2>,
|-3/2> = |00>, |01>, |10>, |11>). Built twice: from the strings, and from the tensors T_{2,m} with the Lindblad form
of the page's Conventions; compared entry by entry, exactly.

  1. at zero splitting: the eight distinct non-zero rates and the five relations among them, written out, with their
     rank (exact);
  1b. the ties as counts: XI, YI, ZX, ZY are the two bits' Majoranas, ZZ their parity; with free jump rates an even
     product decays at twice the rates of the Majoranas it holds, an odd one at twice the parity's plus those it
     lacks (read off the generator's diagonal); the first and fifth relation on class means; a sixth jump XX as the
     control; the turn about the field on the Majoranas; the counted letter cube's points at two rates; with a
     splitting, each odd string and its complement a closed pair, all four centred at 2g(ZZ) plus the four
     Majoranas' rates (c2 at the sample's rates), as the sodium page's map S(X) = ZZ X^T gives for free rates;
  2. with a splitting: the populations, the central transition and the triple quantum keep their rates; a satellite
     cell and a double-quantum cell have the diagonal element -c2 = -3(J0 + J1 + J2), and past the bright pairs'
     exceptional point both eigenvalues of each pair have real part -c2 (exact), so satellites and double quanta
     decay at one rate; at q = 0 the sum of the two satellite cells decays at 3(J0 + J1), the sum of the two
     double-quantum cells (I+^2, rank 2) at 3(J0 + 2J1 + J2) and their difference at 3(J0 + J2) (the controls);
     no choice of the levels' signs gives the paper's pairing; from I_z, under the generator at q = 0 and the pulses, nothing of rank 2 is reached (for every J, by
     block decoupling, with the dynamic frequency shift too; at q = 1 and at J0 = J1 = J2, read at points);
     the tensors and weights the page's table names;
  3. the paper's three-rate inversion (satellite sum at C(J0 + J1), as at zero splitting; double-quantum sum at
     C(J0 + J2), which at zero splitting is the rate of their difference), fed rates from past the exceptional
     point (exact);
  4. R. Auccaise et al., J. Magn. Reson. 192, 17 (2008), Tables 1 and 2: J1 + J2 from the coherences against the
     populations (arithmetic on printed numbers), and the printed a1 against a2 with the error model stated: the
     printed errors in quadrature, and a bound that holds however the one-digit numbers "4 +- 1" are read.
"""
import sys
import sympy as sp

CLAIMS = []


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {label} {detail}")


J0, J1, J2, q, t = sp.symbols("J0 J1 J2 q t", positive=True)
I2 = sp.eye(2)
PX = sp.Matrix([[0, 1], [1, 0]])
PY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
PZ = sp.diag(1, -1)
LET = {"I": I2, "X": PX, "Y": PY, "Z": PZ}
LABELS = [a + b for a in "IXYZ" for b in "IXYZ"]
S = {s: sp.kronecker_product(LET[s[0]], LET[s[1]]) for s in LABELS}
JUMPS = {"ZZ": 3 * J0 / 2, "ZX": 3 * J1 / 2, "ZY": 3 * J1 / 2, "XI": 3 * J2 / 2, "YI": 3 * J2 / 2}
c1 = 3 * (J1 + J2)
c2 = 3 * (J0 + J1 + J2)


def vec(m):
    return m.reshape(16, 1)  # row-major


def supop(f):
    cols = []
    for k in range(16):
        e = sp.zeros(16, 1)
        e[k] = 1
        cols.append(vec(f(e.reshape(4, 4))))
    return sp.Matrix.hstack(*cols)


def gen_strings(qq, jumps=JUMPS):
    H = qq * S["ZZ"]
    return supop(lambda r: -sp.I * (H * r - r * H)
                 + sum((c * (S[g] * r * S[g] - r) for g, c in jumps.items()), sp.zeros(4)))


def gen_tensors(qq):
    ms = [sp.Rational(3, 2) - k for k in range(4)]
    Iz = sp.diag(*ms)
    Ip = sp.zeros(4)
    for k in range(1, 4):
        Ip[k - 1, k] = sp.sqrt(sp.Rational(15, 4) - ms[k] * (ms[k] + 1))
    Im = Ip.T
    T = {0: (3 * Iz ** 2 - sp.Rational(15, 4) * sp.eye(4)) / sp.sqrt(6),
         1: -(Iz * Ip + Ip * Iz) / 2, -1: (Iz * Im + Im * Iz) / 2,
         2: Ip ** 2 / 2, -2: Im ** 2 / 2}
    Jm = {0: J0, 1: J1, -1: J1, 2: J2, -2: J2}
    wq = qq / (sp.sqrt(6) / 2)  # q*ZZ = omega_Q * T20, with sqrt(6)*T20 = 3*ZZ
    H = wq * T[0]

    def f(r):
        out = -sp.I * (H * r - r * H)
        for m, A in T.items():
            out += Jm[m] * (A * r * A.H - (A.H * A * r + r * A.H * A) / 2)
        return out
    return supop(f)


def cell(a, b):
    m = sp.zeros(4)
    m[a, b] = 1
    return m


def anticommutes(a, b):
    return (S[a] * S[b] + S[b] * S[a]) == sp.zeros(4)


print("0. the generator, built twice")
Ls, Lt = gen_strings(q), gen_tensors(q)
claim("from the strings = from the tensors T_{2,m}, entry by entry, q symbolic",
      sp.simplify(Ls - Lt) == sp.zeros(16))

print("1. zero splitting: rates and relations")
L0 = gen_strings(0)
rate = {}
for s in LABELS:
    out = L0 * vec(S[s])
    lam = sp.simplify((vec(S[s]).H * out)[0] / 4)
    rate[s] = sp.expand(-lam)
    assert sp.simplify(out + rate[s] * vec(S[s])) == sp.zeros(16, 1), s
by_count = {s: sp.expand(2 * sum(c for g, c in JUMPS.items() if anticommutes(s, g))) for s in LABELS}
claim("every string an eigenoperator, its rate twice the summed rate of the jumps it anticommutes with",
      all(sp.simplify(rate[s] - by_count[s]) == 0 for s in LABELS))
distinct = sorted({r for r in rate.values() if r != 0}, key=str)
claim("eight distinct non-zero rates", len(distinct) == 8, str(distinct))
r = rate
relations = {
    "mean: central and triple quantum = mean of the two longitudinal rates, r(XX) = (r(IZ) + r(ZI))/2":
        r["XX"] - (r["IZ"] + r["ZI"]) / 2,
    "sum: quadrupolar order = sum of the two longitudinal rates, r(ZZ) = r(IZ) + r(ZI)":
        r["ZZ"] - (r["IZ"] + r["ZI"]),
    "satellites: r(ZX) - r(IX) = r(ZI)": r["ZX"] - r["IX"] - r["ZI"],
    "double quanta: r(XI) - r(XZ) = r(IZ)": r["XI"] - r["XZ"] - r["IZ"],
    "bright against dark: r(IX) - r(XZ) = (r(IZ) - r(ZI))/2": r["IX"] - r["XZ"] - (r["IZ"] - r["ZI"]) / 2,
}
for label, e in relations.items():
    claim(label, sp.simplify(e) == 0)
# independence: the five relations as vectors over the eight rates
R8 = [r["IZ"], r["ZI"], r["ZZ"], r["XX"], r["IX"], r["XZ"], r["ZX"], r["XI"]]
x = sp.symbols("x0:8")
lin = [x[3] - (x[0] + x[1]) / 2, x[2] - x[0] - x[1], x[6] - x[4] - x[1], x[7] - x[5] - x[0],
       x[4] - x[5] - (x[0] - x[1]) / 2]
M = sp.Matrix([[sp.diff(e, xi) for xi in x] for e in lin])
G = sp.Matrix([[sp.diff(rr, Jk) for Jk in (J0, J1, J2)] for rr in R8])
claim("the five relations independent (rank 5), the eight rates of rank 3, 5 + 3 = 8",
      M.rank() == 5 and G.rank() == 3 and (M * G) == sp.zeros(5, 3))
# the exact condition for the mean relation, rates g of the five jumps free
gs = sp.symbols("gZZ gZX gZY gXI gYI", positive=True)
free = dict(zip(["ZZ", "ZX", "ZY", "XI", "YI"], gs))
rf = {s: sp.expand(2 * sum(c for g_, c in free.items() if anticommutes(s, g_))) for s in LABELS}
cond = (free["ZX"] - free["ZY"]) - (free["XI"] - free["YI"])
claim("with five free rates, r(XX) - mean = +cond and r(YY) - mean = -cond, cond = g(ZX) - g(ZY) - (g(XI) - g(YI))",
      sp.simplify(rf["XX"] - (rf["IZ"] + rf["ZI"]) / 2 - cond) == 0
      and sp.simplify(rf["YY"] - (rf["IZ"] + rf["ZI"]) / 2 + cond) == 0)
cond2 = (free["ZX"] - free["ZY"]) + (free["XI"] - free["YI"])
claim("on XY and YX the mean relation reads -+cond2, cond2 = g(ZX) - g(ZY) + (g(XI) - g(YI)): the whole class "
      "{XX, YY, XY, YX} at the mean needs g(ZX) = g(ZY) and g(XI) = g(YI)",
      sp.simplify(rf["XY"] - (rf["IZ"] + rf["ZI"]) / 2 + cond2) == 0
      and sp.simplify(rf["YX"] - (rf["IZ"] + rf["ZI"]) / 2 - cond2) == 0
      and sp.solve([cond, cond2], [free["ZY"], free["YI"]], dict=True) == [{free["ZY"]: free["ZX"], free["YI"]: free["XI"]}])
claim("with five free rates the sum, satellite and double-quantum relations hold identically (they test which jumps "
      "act), and the fifth, read on IX and XZ, is -cond, so it needs the same condition as the first",
      sp.simplify(rf["ZZ"] - rf["IZ"] - rf["ZI"]) == 0
      and sp.simplify(rf["ZX"] - rf["IX"] - rf["ZI"]) == 0 and sp.simplify(rf["ZY"] - rf["IY"] - rf["ZI"]) == 0
      and sp.simplify(rf["XI"] - rf["XZ"] - rf["IZ"]) == 0 and sp.simplify(rf["YI"] - rf["YZ"] - rf["IZ"]) == 0
      and sp.simplify(rf["IX"] - rf["XZ"] - (rf["IZ"] - rf["ZI"]) / 2 + cond) == 0)
# control: unequal rates inside the J1 pair break the mean relation
bad = dict(JUMPS)
bad["ZY"] = 3 * J1
Lb = gen_strings(0, bad)
rb = {s: sp.expand(-(vec(S[s]).H * Lb * vec(S[s]))[0] / 4) for s in ("XX", "IZ", "ZI")}
claim("control: ZY at a rate other than ZX's, XI and YI equal, breaks the mean relation",
      sp.simplify(rb["XX"] - (rb["IZ"] + rb["ZI"]) / 2) != 0)

print("1b. the ties as counts of Majoranas")
# the four jumps XI, YI (site 1) and ZX, ZY (site 2) are the Jordan-Wigner Majoranas of the two bits, ZZ their parity
MAJ = ["XI", "YI", "ZX", "ZY"]
Lfree = gen_strings(0, free)
rdiag = {s: sp.expand(-(vec(S[s]).H * Lfree * vec(S[s]))[0] / 4) for s in LABELS}


def as_product(s):
    for k in range(5):
        for sub in __import__("itertools").combinations(MAJ, k):
            M = sp.eye(4)
            for m in sub:
                M = M * S[m]
            if any(M == c * S[s] for c in (1, -1, sp.I, -sp.I)):
                return sub


subsets = {s: as_product(s) for s in LABELS}
count = {s: 2 * sum(free[m] for m in sub) if len(sub) % 2 == 0
         else 2 * (free["ZZ"] + sum(free[m] for m in MAJ if m not in sub)) for s, sub in subsets.items()}
claim("with five free rates, read off the generator's diagonal: every string is one product of the four Majoranas, "
      "an even product decays at twice the summed rate of the Majoranas it holds, an odd one at twice the sum of the "
      "parity's rate and the rates of those it lacks", len(set(subsets.values())) == 16
      and all(sp.simplify(rdiag[s] - count[s]) == 0 for s in LABELS))
mean = lambda ss: sum(rdiag[s] for s in ss) / len(ss)
claim("with five free rates the first and the fifth relation hold on class means: the mean of {XX, XY, YX, YY} is "
      "(r(IZ) + r(ZI))/2, and the mean of {IX, IY} less that of {XZ, YZ} is (r(IZ) - r(ZI))/2",
      sp.simplify(mean(["XX", "XY", "YX", "YY"]) - (rdiag["IZ"] + rdiag["ZI"]) / 2) == 0
      and sp.simplify(mean(["IX", "IY"]) - mean(["XZ", "YZ"]) - (rdiag["IZ"] - rdiag["ZI"]) / 2) == 0)
claim("each of the four Majoranas lights exactly two of {XX, XY, YX, YY}, the parity none",
      all(sum(anticommutes(s, m) for s in ("XX", "XY", "YX", "YY")) == 2 for m in MAJ)
      and not any(anticommutes(s, "ZZ") for s in ("XX", "XY", "YX", "YY")))
# control: a sixth jump that is no Majorana and not the parity, XX, breaks the sum relation
six = dict(free)
six["XX"] = sp.Symbol("gXX", positive=True)
L6 = gen_strings(0, six)
r6 = {s: sp.expand(-(vec(S[s]).H * L6 * vec(S[s]))[0] / 4) for s in ("ZZ", "IZ", "ZI")}
claim("control: a sixth jump XX breaks r(ZZ) = r(IZ) + r(ZI)", sp.simplify(r6["ZZ"] - r6["IZ"] - r6["ZI"]) != 0)
# the turn about the field: I_z = ZI + IZ/2 turns each site's two Majoranas into each other, site 1 at twice the angle
ph = sp.Symbol("phi", real=True)
Uz = sp.diag(*[sp.exp(-sp.I * ph * v) for v in (sp.Rational(3, 2), sp.Rational(1, 2), -sp.Rational(1, 2), -sp.Rational(3, 2))])
turned = {m: sp.simplify(Uz * S[m] * Uz.H) for m in MAJ}
# a dissipator that commutes with every turn about the field has one rate for the two Majoranas of each bit
Dfree = lambda M: sum((free[j] * (S[j] * M * S[j] - M) for j in free), sp.zeros(4))
eqs = []
for s in LABELS:
    eqs += list(sp.expand(sp.simplify(Dfree(Uz * S[s] * Uz.H) - Uz * Dfree(S[s]) * Uz.H)))
coeffs = set()
for e in eqs:
    e = sp.expand(sp.expand_trig(sp.simplify(e.rewrite(sp.cos))))
    if e != 0:
        coeffs |= set(sp.Poly(e, sp.cos(ph), sp.sin(ph)).coeffs())
sol = sp.solve(list(coeffs), [free["YI"], free["ZY"]], dict=True)
# the conditions are exactly these two: the coefficient equations span the same space as the two differences
span_ok = sp.Matrix([[sp.diff(c_, g_) for g_ in gs] for c_ in coeffs]).rank() == 2
claim("a dissipator of the five free rates commutes with every turn about the field exactly when g(XI) = g(YI) and "
      "g(ZX) = g(ZY), g(ZZ) free", sol == [{free["YI"]: free["XI"], free["ZY"]: free["ZX"]}] and span_ok, str(sol))
claim("I_z = ZI + IZ/2, and the turn about the field by phi turns XI into cos(2 phi) XI + sin(2 phi) YI and ZX into "
      "cos(phi) ZX + sin(phi) ZY (site 1 at twice the angle of site 2)",
      sp.diag(sp.Rational(3, 2), sp.Rational(1, 2), -sp.Rational(1, 2), -sp.Rational(3, 2)) == S["ZI"] + S["IZ"] / 2
      and sp.simplify(turned["XI"] - (sp.cos(2 * ph) * S["XI"] + sp.sin(2 * ph) * S["YI"])) == sp.zeros(4)
      and sp.simplify(turned["ZX"] - (sp.cos(ph) * S["ZX"] + sp.sin(ph) * S["ZY"])) == sp.zeros(4))
# the letter cube read as a count, (k_Z, k_X, k_Y), cannot tell the sites apart: the points that hold strings at two
# rates are exactly the pairs the swap of the two bits exchanges
pt = lambda s: (sum(c in "XY" for c in s), sum(c in "YZ" for c in s), sum(c in "XZ" for c in s))
points = {}
for s in LABELS:
    points.setdefault(pt(s), []).append(s)
two_str = {p_: ss for p_, ss in points.items() if len(ss) == 2}
rfree_xy = (rdiag["XY"], rdiag["YX"])
claim("the sixteen strings sit on ten points, six of them hold two strings, each a swap pair; XY and YX share c1 at "
      "the sample's rates and part with free rates",
      len(points) == 10 and len(two_str) == 6 and all(ss[0] == ss[1][::-1] for ss in two_str.values())
      and sp.simplify(r["XY"] - r["YX"]) == 0 and sp.simplify(rfree_xy[0] - rfree_xy[1]) != 0)
split_pts = {p_: ss for p_, ss in points.items() if len({r[s] for s in ss}) > 1}
pairs_ok = all(len(ss) == 2 and ss[0] == ss[1][::-1] for ss in split_pts.values())
diffs = {tuple(ss): sp.factor(r[ss[1]] - r[ss[0]]) for ss in split_pts.values()}
print("    points of the counted cube that hold two rates: " + "; ".join(f"{a}/{b}: {d}" for (a, b), d in diffs.items()))
claim("on the counted letter cube exactly five points hold strings at two rates, each a pair the swap of the two bits "
      "exchanges; at the four odd pairs the two rates differ by c1 = 3(J1 + J2), at IZ/ZI by 6(J2 - J1)",
      len(split_pts) == 5 and pairs_ok
      and all(sp.simplify(abs(d) - c1) == 0 or sp.simplify(d ** 2 - 36 * (J2 - J1) ** 2) == 0 for d in diffs.values())
      and sum(1 for d in diffs.values() if sp.simplify(d ** 2 - c1 ** 2) == 0) == 4)
# with a splitting, free rates: the parity carries each odd string onto its complement, a closed pair with one centre
Lqf = gen_strings(q, free)
centre = 2 * free["ZZ"] + sum(free[m] for m in MAJ)
ok_pairs = True
for s in LABELS:
    if len(subsets[s]) % 2 == 0:
        continue
    comp = [k for k in LABELS if set(subsets[k]) == set(MAJ) - set(subsets[s])][0]
    bi = [vec(S[s]) / 2, vec(S[comp]) / 2]
    blk = sp.Matrix(2, 2, lambda a_, b_: sp.simplify((bi[a_].H * Lqf * bi[b_])[0]))
    out_of = [sp.simplify(Lqf * bi[j] - (bi[0] * blk[0, j] + bi[1] * blk[1, j])) for j in (0, 1)]
    ok_pairs &= (all(o == sp.zeros(16, 1) for o in out_of) and sp.simplify(-(blk[0, 0] + blk[1, 1]) / 2 - centre) == 0
                 and sp.simplify(blk[0, 1] * blk[1, 0] + 4 * q ** 2) == 0)
claim("with five free rates and a splitting, each odd string and its complement (the parity times it) form a closed "
      "pair coupled by 2q, every pair centred at 2g(ZZ) plus the four Majoranas' rates, which is c2 at the sample's rates",
      ok_pairs and sp.simplify(centre.subs(dict(zip(gs, JUMPS.values()))) - c2) == 0)
# the sodium page's one map S(X) = ZZ X^T, read with free rates: S L S^-1 = -L + 2 L_ZZ - 2 (sum of the Majoranas' rates)
Ssup = supop(lambda X: S["ZZ"] * X.T)
LZZ = supop(lambda X: free["ZZ"] * (S["ZZ"] * X * S["ZZ"] - X))
claim("with five free rates and a splitting, S(X) = ZZ X^T gives S L S^-1 = -L + 2 L_ZZ - 2(g(XI) + g(YI) + g(ZX) + g(ZY))",
      sp.simplify(Ssup * Lqf * Ssup.inv() - (-Lqf + 2 * LZZ - 2 * sum(free[m] for m in MAJ) * sp.eye(16))) == sp.zeros(16))
# the pairs the splitting mixes sum to 6(J0 + J1 + J2) at the sample's rates (the Key page's pair sums)
claim("at the sample's rates the four pairs the splitting mixes, {IX, ZY}, {IY, ZX}, {XI, YZ}, {YI, XZ}, each sum to "
      "6(J0 + J1 + J2)", all(sp.simplify(r[a_] + r[b_] - 2 * c2) == 0 for a_, b_ in (("IX", "ZY"), ("IY", "ZX"), ("XI", "YZ"), ("YI", "XZ"))))
d_a = sp.simplify((rdiag["ZY"] - rdiag["IX"]) / 2)
d_b = sp.simplify((rdiag["ZX"] - rdiag["IY"]) / 2)
claim("the satellite pairs' half-differences: g(XI) + g(YI) + g(ZX) - g(ZY) for IX and ZY, the same with ZX and ZY "
      "exchanged for IY and ZX; both 3J2 at the sample's rates",
      sp.simplify(d_a - (free["XI"] + free["YI"] + free["ZX"] - free["ZY"])) == 0
      and sp.simplify(d_b - (free["XI"] + free["YI"] + free["ZY"] - free["ZX"])) == 0
      and all(sp.simplify(d_.subs(dict(zip(gs, JUMPS.values()))) - 3 * J2) == 0 for d_ in (d_a, d_b)))

print("2. with a splitting")
dark = {"IZ": 6 * J1, "ZI": 6 * J2, "ZZ": 6 * (J1 + J2)}
ok = all(sp.simplify(Ls * vec(S[s]) + v * vec(S[s])) == sp.zeros(16, 1) for s, v in dark.items())
central = vec(S["XX"] + S["YY"])
tq = vec(S["XX"] - S["YY"])
ok = ok and all(sp.simplify(Ls * w + c1 * w) == sp.zeros(16, 1)
                for w in (central, tq, vec(S["XY"] - S["YX"]), vec(S["XY"] + S["YX"])))
claim("populations at 6J1, 6J2, 6(J1 + J2), central transition and triple quantum at c1, at every q", ok)
idx = lambda a, b: 4 * a + b
sat = [idx(0, 1), idx(2, 3)]
dq = [idx(0, 2), idx(1, 3)]
diag_ok = all(sp.simplify(Ls[k, k] + c2 - s * 2 * sp.I * q) == 0
              for k, s in zip(sat + dq, (-1, 1, -1, 1)))
claim("a satellite cell and a double-quantum cell: diagonal element -c2 -+ 2iq, c2 = 3(J0 + J1 + J2)", diag_ok)
claim("c2 is the mean of each cell's two string classes: (r(IX) + r(ZX))/2 = (r(XZ) + r(XI))/2 = c2",
      sp.simplify((r["IX"] + r["ZX"]) / 2 - c2) == 0 and sp.simplify((r["XZ"] + r["XI"]) / 2 - c2) == 0)
for name, ks, Jx in (("satellites", sat, J2), ("double quanta", dq, J1)):
    B = Ls.extract(ks, ks)
    # the pair is closed: nothing else feeds it
    others = [k for k in range(16) if k not in ks]
    closed = Ls.extract(ks, others) == sp.zeros(2, 14) and Ls.extract(others, ks) == sp.zeros(14, 2)
    ev = list(B.eigenvals().keys())
    want = [-c2 + sp.sqrt(9 * Jx ** 2 - 4 * q ** 2), -c2 - sp.sqrt(9 * Jx ** 2 - 4 * q ** 2)]
    match = all(any(sp.simplify(e - w) == 0 for e in ev) for w in want)
    claim(f"{name}: a closed pair with eigenvalues -c2 +- sqrt(9{Jx}^2 - 4q^2), real part -c2 once 2q > 3{Jx}",
          closed and match)
# control: at q = 0 the sum of the two satellite cells decays at 3(J0 + J1), not at c2
v = vec(cell(0, 1) + cell(2, 3))
claim("control: at q = 0 the sum of the two satellite cells decays at 3(J0 + J1), not at c2",
      sp.simplify(L0 * v + 3 * (J0 + J1) * v) == sp.zeros(16, 1) and sp.simplify(3 * (J0 + J1) - c2) != 0)
# the double quanta at q = 0: the sum of the two cells is sigma+ (x) I, the XI class, T_{2,+-2} (prop. to I+^2),
# at 3(J0 + 2J1 + J2); the difference is sigma+ (x) Z, the XZ class, T_{3,+-2}, at 3(J0 + J2)
vs, vd = vec(cell(0, 2) + cell(1, 3)), vec(cell(0, 2) - cell(1, 3))
claim("at q = 0 the sum of the two double-quantum cells decays at 3(J0 + 2J1 + J2), their difference at 3(J0 + J2)",
      sp.simplify(L0 * vs + 3 * (J0 + 2 * J1 + J2) * vs) == sp.zeros(16, 1)
      and sp.simplify(L0 * vd + 3 * (J0 + J2) * vd) == sp.zeros(16, 1))
ms_ = [sp.Rational(3, 2) - k for k in range(4)]
Ip = sp.zeros(4)
for k in range(1, 4):
    Ip[k - 1, k] = sp.sqrt(sp.Rational(15, 4) - ms_[k] * (ms_[k] + 1))
claim("I+^2 is 2*sqrt(3) times the sum of the two double-quantum cells, and at J0 = J1 = J2 it decays at the "
      "rank-2 rate of quadrupolar order ZZ, not at the rate 3(J0 + J2)",
      Ip ** 2 == 2 * sp.sqrt(3) * (cell(0, 2) + cell(1, 3))
      and sp.simplify((3 * (J0 + 2 * J1 + J2) - r["ZZ"]).subs({J0: J1, J2: J1})) == 0
      and sp.simplify((3 * (J0 + J2) - r["ZZ"]).subs({J0: J1, J2: J1})) != 0)
# the levels' phases: rho_ij -> e^{i(phi_i - phi_j)} rho_ij turns rho01 against rho23 and rho02 against rho13 by one
# and the same angle, so a sum becomes a difference in both pairs or in neither
import itertools
hit = []
for s_ in itertools.product((1, -1), repeat=4):
    vsat = vec(s_[0] * s_[1] * cell(0, 1) + s_[2] * s_[3] * cell(2, 3))
    vdq = vec(s_[0] * s_[2] * cell(0, 2) + s_[1] * s_[3] * cell(1, 3))
    if sp.simplify(L0 * vsat + 3 * (J0 + J1) * vsat) == sp.zeros(16, 1) and \
            sp.simplify(L0 * vdq + 3 * (J0 + J2) * vdq) == sp.zeros(16, 1):
        hit.append(s_)
# a cell combination c01 + e^{i t} c23 is an eigenvector of L0 only for e^{i t} = +-1, so signs are all the phases
th_ = sp.symbols("theta", real=True)
w_ = vec(cell(0, 1) + sp.exp(sp.I * th_) * cell(2, 3))
Lw = sp.simplify(L0 * w_)
ratio_ok = sp.simplify(Lw[idx(0, 1)] * w_[idx(2, 3)] - Lw[idx(2, 3)] * w_[idx(0, 1)])
kfac = sp.simplify(ratio_ok / (sp.exp(2 * sp.I * th_) - 1))
w2_ = vec(cell(0, 2) + sp.exp(sp.I * th_) * cell(1, 3))
Lw2 = sp.simplify(L0 * w2_)
ratio2 = sp.simplify(Lw2[idx(0, 2)] * w2_[idx(1, 3)] - Lw2[idx(1, 3)] * w2_[idx(0, 2)])
kfac2 = sp.simplify(ratio2 / (sp.exp(2 * sp.I * th_) - 1))
claim("over every choice of the levels' signs, no basis has the satellite sum at 3(J0 + J1) and the double-quantum sum "
      "at 3(J0 + J2) together (the paper's pairing); and a pair combination is an eigenvector only at relative phase "
      "+-1 (its residual is a nonzero constant times e^(2i theta) - 1)",
      hit == [] and all(sp.simplify(sp.diff(k_, th_)) == 0 and sp.simplify(k_) != 0 for k_ in (kfac, kfac2)))
# below the exceptional point the pair's rates move from their q = 0 values only at order q^2
ser = sp.series(sp.sqrt(9 * J2 ** 2 - 4 * q ** 2), q, 0, 4).removeO()
print(f"    (arithmetic on the eigenvalue formula: sqrt(9J2^2 - 4q^2) = {ser} + O(q^4), the rates move at order q^2)")
# reachability: from I_z, under the generator at q = 0 and the pulses (I_x, I_y, I_z), no rank-2 part is reached;
# at q != 0 it is (the control). Rank 2 is the span of the five Dirac strings ZZ, ZX, ZY, XI, YI.
Ipm = Ip
Ix_, Iy_, Iz_ = (Ipm + Ipm.T) / 2, (Ipm - Ipm.T) / (2 * sp.I), sp.diag(*ms_)
ad = lambda A: supop(lambda rr: A * rr - rr * A)
pt = {J0: 7, J1: 3, J2: 2}


def reach(Lq):
    gens = [Lq, ad(Ix_), ad(Iy_), ad(Iz_)]
    basis = sp.Matrix(vec(Iz_))
    frontier = [vec(Iz_)]
    while frontier:
        new = []
        for v_ in frontier:
            for g_ in gens:
                w = sp.expand(g_ * v_)
                cand = basis.row_join(w)
                if cand.rank() > basis.rank():
                    basis = cand
                    new.append(w)
        frontier = new
    return basis


rank2 = sp.Matrix.hstack(*[vec(S[s]) for s in ("ZZ", "ZX", "ZY", "XI", "YI")])
R0 = reach(gen_strings(0).subs(pt))
Rq = reach(gen_strings(1).subs(pt))
others = sp.Matrix.hstack(*[vec(S[s]) for s in LABELS if s not in ("ZZ", "ZX", "ZY", "XI", "YI")])
decoupled = all((rank2.H * G_ * others).applyfunc(sp.simplify) == sp.zeros(5, 11)
                for G_ in (gen_strings(0), ad(Ix_), ad(Iy_), ad(Iz_)))
K1s, K2s = sp.symbols("K1 K2", real=True)
Hds = -3 * K2s * S["ZI"] - 3 * K1s * S["IZ"]     # the sodium page's dynamic frequency shift, a diagonal Hamiltonian
Lds = gen_strings(0) + supop(lambda rr: -sp.I * (Hds * rr - rr * Hds))
decoupled_ds = (rank2.H * Lds * others).applyfunc(sp.simplify) == sp.zeros(5, 11)
claim("for every J, the generator at q = 0 (with the dynamic frequency shift or without it) and the pulses I_x, I_y, "
      "I_z map the eleven strings outside rank 2 into themselves, and I_z lies among them: nothing of rank 2 is "
      "reached from it", decoupled and decoupled_ds and (rank2.H * vec(Iz_)) == sp.zeros(5, 1))
in_R0 = lambda s: R0.row_join(vec(S[s])).rank() == R0.rank()
Req = reach(gen_strings(0).subs({J0: 2, J1: 2, J2: 2}))
in_Req = lambda s: Req.row_join(vec(S[s])).rank() == Req.rank()
claim(f"read at points: at J = (7, 3, 2) the span reached from I_z has dimension {R0.rank()} at q = 0, holding IX and "
      f"XZ (R2f and the rank-3 double quantum); at J0 = J1 = J2 it has dimension {Req.rank()}, rank 1 alone, without "
      f"them; at q = 1 it has dimension {Rq.rank()} and reaches rank 2 (the control)",
      in_R0("IX") and in_R0("XZ") and not in_R0("XI") and Req.rank() == 3 and not in_Req("IX")
      and (rank2.H * Rq).applyfunc(sp.simplify) != sp.zeros(5, Rq.shape[1]))
# the tensors the table names, and the weights
T21 = Iz_ * Ip + Ip * Iz_
claim("I_z I+ + I+ I_z = 2 sqrt(3) (cell 01 - cell 23), the ZX class, rank 2, at 3(J0 + J1 + 2J2) = r(ZX)",
      T21 == 2 * sp.sqrt(3) * (cell(0, 1) - cell(2, 3))
      and sp.simplify(L0 * vec(T21) + r["ZX"] * vec(T21)) == sp.zeros(16, 1))
Im_ = Ip.T
T32 = Im_ * Ip ** 3 - Ip ** 3 * Im_
claim("[I-, I+^3] (rank 3, order 2) is proportional to cell 02 - cell 13, the XZ class",
      sp.simplify(T32 - T32[0, 2] * (cell(0, 2) - cell(1, 3))) == sp.zeros(4) and T32[0, 2] != 0)
claim("cell 01 = (IX + i IY + ZX + i ZY)/4 and cell 02 = (XI + i YI + XZ + i YZ)/4: half on each of two classes",
      cell(0, 1) == (S["IX"] + sp.I * S["IY"] + S["ZX"] + sp.I * S["ZY"]) / 4
      and cell(0, 2) == (S["XI"] + sp.I * S["YI"] + S["XZ"] + sp.I * S["YZ"]) / 4)
coef = lambda M, s: sp.simplify((S[s].H * M).trace() / 4)
w_z = {s: coef(Iz_, s) for s in ("ZI", "IZ")}
nz = sum(c ** 2 for c in w_z.values())
w_x = {"sat": coef(Ix_, "IX") ** 2, "ct": coef(Ix_, "XX") ** 2 + coef(Ix_, "YY") ** 2}
claim("I_z = ZI + IZ/2 and I_x = (sqrt3/2) IX + (XX + YY)/2, so the longitudinal weights 4/5 : 1/5 and the transverse "
      "3/5 : 2/5",
      Iz_ == S["ZI"] + S["IZ"] / 2 and sp.simplify(Ix_ - (sp.sqrt(3) / 2 * S["IX"] + (S["XX"] + S["YY"]) / 2)) == sp.zeros(4)
      and w_z["ZI"] ** 2 / nz == sp.Rational(4, 5) and w_x["sat"] / (w_x["sat"] + w_x["ct"]) == sp.Rational(3, 5))
eq = {J0: J1, J2: J1}
ten = [s for s in LABELS if s not in ("II", "ZZ", "ZX", "ZY", "XI", "YI")]
claim("at J0 = J1 = J2 the ten strings of ranks 1 and 3 decay at 6J1 and the five of rank 2 at 12J1",
      all(sp.simplify(r[s].subs(eq) - 6 * J1) == 0 for s in ten)
      and all(sp.simplify(r[s].subs(eq) - 12 * J1) == 0 for s in ("ZZ", "ZX", "ZY", "XI", "YI")))

print("3. a zero-splitting inversion fed rates from past the exceptional point")
C = sp.symbols("C", positive=True)
a1, a2, a3 = C * (J0 + J1 + J2), C * (J0 + J1 + J2), C * (J1 + J2)  # rates past the point, trade's units C*J
j0, j1, j2 = sp.symbols("j0 j1 j2")
sol = sp.solve([C * (j0 + j1) - a1, C * (j0 + j2) - a2, C * (j1 + j2) - a3], [j0, j1, j2], dict=True)[0]
claim("the inversion a1 = C(J0 + J1), a2 = C(J0 + J2), a3 = C(J1 + J2) returns J1' = J2' = (J1 + J2)/2 and "
      "J0' = J0 + (J1 + J2)/2",
      sp.simplify(sol[j1] - (J1 + J2) / 2) == 0 and sp.simplify(sol[j2] - (J1 + J2) / 2) == 0
      and sp.simplify(sol[j0] - (J0 + (J1 + J2) / 2)) == 0)
print("    (the third equation alone gives J1' + J2' = a3/C, whatever a1 and a2 are: arithmetic, not a check)")

print("4. Auccaise et al. 2008, Table 1 (units 1e-9 s, errors as printed)")
coh = {"J1": (sp.Rational(4), sp.Rational(1)), "J2": (sp.Rational(34, 10), sp.Rational(5, 10))}
pop = {"J1": (sp.Rational(38, 10), sp.Rational(4, 10)), "J2": (sp.Rational(36, 10), sp.Rational(3, 10))}
s_coh = coh["J1"][0] + coh["J2"][0]
s_pop = pop["J1"][0] + pop["J2"][0]
e_coh = sp.sqrt(coh["J1"][1] ** 2 + coh["J2"][1] ** 2)
e_pop = sp.sqrt(pop["J1"][1] ** 2 + pop["J2"][1] ** 2)
diff = s_coh - s_pop
comb = sp.sqrt(e_coh ** 2 + e_pop ** 2)
print(f"    J1 + J2: coherences {s_coh} +- {float(e_coh):.2f}, populations {s_pop} +- {float(e_pop):.2f}, "
      f"difference {diff} +- {float(comb):.2f} (errors added in quadrature, as if independent)")
print("    (transcription, not a check: 7.4 in both columns; C cancels, a3/C against the mean of 2CJ1 and 2CJ2 "
      "over 2C)")
print("    (the J1 + J2 comparison is arithmetic on printed numbers: the difference is exactly 0)")
t2 = {"J1": [sp.Rational(x, 10) for x in (38, 31, 35, 31, 32, 38)],
      "J2": [sp.Rational(x, 10) for x in (36, 36, 37, 37, 38, 36)]}
m1, m2 = sum(t2["J1"]) / 6, sum(t2["J2"]) / 6
print(f"    Table 2, six states: J1 + J2 = " + ", ".join(str(a + b) for a, b in zip(t2["J1"], t2["J2"]))
      + f"; means {float(m1):.3f} + {float(m2):.3f}")
claim("Table 2's six-state means round to the printed 3.4 and 3.7", round(float(m1), 1) == 3.4 and round(float(m2), 1) == 3.7)
J0p = sp.Rational(14)
a1C, a2C = J0p + coh["J1"][0], J0p + coh["J2"][0]
print(f"    a1/C = J0 + J1 = {a1C}, a2/C = J0 + J2 = {a2C} (the printed densities are an exact inversion of the "
      f"three rates); the generator's zero-splitting assignment would put a2/C at J1 + J2 = {s_pop} above a1/C")
# error model: a1/C - a2/C = J1 - J2 of the coherences column; its error is between |1 - 0.5| and 1 + 0.5 under any
# correlation of the two printed errors (Cauchy-Schwarz), 1.12 in quadrature
e1, e2 = coh["J1"][1], coh["J2"][1]
sig_q, sig_lo, sig_hi = sp.sqrt(e1 ** 2 + e2 ** 2), abs(e1 - e2), e1 + e2
d = a1C - a2C
print(f"    a1/C - a2/C = {d}, error {float(sig_q):.2f} in quadrature, between {sig_lo} and {sig_hi} under any "
      f"correlation (1e-9 s; a1 - a2 about 7 s^-1)")
claim("a1/C - a2/C = J1 - J2 of the coherences column is within one quadrature error of zero at the printed values: "
      "consistent with one rate c2", abs(d) <= sig_q)
pred = s_pop                                     # the generator's zero-splitting assignment: a2/C - a1/C = J1 + J2
dev = pred - (a2C - a1C)
sig_pred = pop["J1"][1] ** 2 + pop["J2"][1] ** 2
sig_tot = sp.sqrt(sig_hi ** 2 + sig_pred)
print(f"    the generator's zero-splitting assignment puts a2/C - a1/C at J1 + J2 = {pred}; measured {a2C - a1C}; "
      f"deviation {dev} = {float(dev / sig_q):.1f} quadrature errors, {float(dev / sig_tot):.2f} at the worst "
      f"correlation with the prediction's own error added in quadrature")
sig_qt = sp.sqrt(sig_q ** 2 + sig_pred)
print(f"    at the printed values: {float(dev / sig_qt):.1f} quadrature errors, the prediction's own error included")
# however the one-digit coherence numbers "4 +- 1" are read (value 3.5 to 4.5, error 0.5 to 1.5) and the two-digit
# ones (+-0.05): the smallest deviation over the largest quadrature error
h = sp.Rational(1, 20)
dev_min = (pred - 2 * h) - ((coh["J2"][0] + h) - (coh["J1"][0] - sp.Rational(1, 2)))
sig_max = sp.sqrt(sp.Rational(3, 2) ** 2 + (e2 + h) ** 2 + (pop["J1"][1] + h) ** 2 + (pop["J2"][1] + h) ** 2)
print(f"    however the printed digits are read: at least {float(dev_min / sig_max):.1f} quadrature errors")
sig_lin = sp.sqrt((sp.Rational(3, 2) + e2 + h) ** 2 + (pop["J1"][1] + h) ** 2 + (pop["J2"][1] + h) ** 2)
print(f"    with the two coherence errors added linearly as well: at least {float(dev_min / sig_lin):.2f}")
claim("that model lies about six and a half quadrature errors off at the printed values, and more than three however "
      "the printed digits are read and however the two coherence errors combine",
      sp.Rational(13, 2) - sp.Rational(1, 10) < dev / sig_qt < sp.Rational(13, 2) + sp.Rational(1, 10)
      and dev_min / sig_lin > 3)
d12 = pop["J1"][0] - pop["J2"][0]
e12 = sp.sqrt(pop["J1"][1] ** 2 + pop["J2"][1] ** 2)
print(f"    populations J1 - J2 = {d12} +- {float(e12):.2f}")
claim("the paper's own form puts C(J1 - J2) between a1 and a2; the populations' J1 - J2 = 0.2 +- 0.5 is consistent "
      "with zero, so that form is not told apart from one rate here", abs(d12) <= e12)

J0true = J0p - s_pop / 2
Cnum = sp.Rational(12, 10) * 10 ** 10
print(f"    a Lorentzian spread of nu_Q with half-width {float(Cnum * s_pop / 2 * sp.Rational(1, 10 ** 9) / (2 * sp.pi)):.1f} "
      f"Hz (C (J1 + J2)/2 / 2 pi) adds to J0' as much as the split reading's (J1 + J2)/2")
print(f"    read past the exceptional point: J0 = J0' - (J1 + J2)/2 = {J0true} (= J0' - J1' under the split "
      f"reading, where J1' = (J1 + J2)/2), and J0 - J1 = {J0true - pop['J1'][0]} (populations J1)")

print(f"\n{sum(CLAIMS)}/{len(CLAIMS)} checks hold")
sys.exit(0 if all(CLAIMS) else 1)
