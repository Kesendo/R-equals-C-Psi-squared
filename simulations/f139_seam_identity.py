"""F139 gate: the seam identity, the F134 wall as a Chebyshev divisor.

THE IDENTITY (G5): build, a priori from the F133 letter system (no table input),
the one-variable polynomial

    Φ_k(y) = (y²−1)·φ₀ + y·φ₁ + φ₂,   φ_i(y) = Σ_t ψ_{i,t}·y^{15−t},
    ψ_{i,t} = R⁻[e_i(P)·e_t({C_M})],

where P = {e₁..e₅, 𝟙} is the 6-letter window, {C_M} the 15 cosine pairs of the
Core (the spin cube minus its poles, negation-closed), and R⁻ the signed W(C₅)
read-off functional G ↦ Σ_{ε,σ} sgn(ε)sgn(σ)·G(ε∘μ + 𝟙 − σδ) at μ = ν+ρ₅,
δ = (2,1,0,−1,−2), all in x = t′² units. Then

    Φ_k = S₁₀·Q_k + R_k   with   Q_k = (−1)^k·P_k   and   deg R_k ≤ 4+k,

where S₁₀(y) = sin(11θ)/sinθ (y = 2cosθ) is the minimal polynomial of the Niven
nodes 2cos(rπ/11) and P_k are the committed F134 column polynomials. Because
(z−z⁻¹)·S₁₀(y) = z¹¹−z⁻¹¹ is the committed wall factor and S₁₀·S_j has
Chebyshev-S coefficients symmetric about 10 on [5,15] while R_k feeds only
degrees ≤ 4+k (below the k-window's low edge), the identity IMPLIES the F134
two-row reflection law: the wall 11 is a Chebyshev divisor, the quotients are
the committed columns, the remainders are the beyond-wall shadow.

GATES:
  G1  the slice identity [t₁^{2μ₁}]X₆ = (−1)^d·Δ̂₅·t′^{−2𝟙}·e_d(Λ₃₆), d = 18−μ₁,
      sampled against the committed F133 read-off Xc at μ₁ ∈ {6, 9, 11, 13, 16}.
  G2  the P-skew ψ_{6−i,t} = −ψ_{i,t} (whole grid; forces ψ₃ ≡ 0).
  G3  the degree lemma deg_y Φ_k ≤ 15.
  G4  the a-priori table: n_(λ₁,k) = (−1)^{λ₁}·b_{λ₁+5} (Chebyshev-S coefficients
      of Φ_k) equals the committed table on ALL parity-live strip cells incl.
      zeros (27 cells).
  G5  the seam division: Q_k = (−1)^k·P_k exactly for k = 0..5, deg R_k ≤ 4+k,
      and the k = 0 remainder is identically ZERO (Φ₀ = S₁₀·(S₁−S₅)).
  G6  the fence: three-row ν = (k₁,1) (l = 1, the law holds there) has remainder
      degree ≤ 4; ν ∈ {(2,2),(3,2),(4,2)} (l = 2, the famous breaks) has
      remainder degree ≥ 8: the breaks are remainder overflow.
  G7  corruption control: bumping one ψ entry must break G4.

Grade: the chain is a priori (letters → ψ → Φ → division; the committed table
enters only as the gate) and every step is exact ℤ; the F127-wall epistemic
class. C# second implementation: WSymplecticClosedForm.AnalyzeSeamIdentity.

Run: python simulations/f139_seam_identity.py   (~2-3 min, exit 0; the time is
the F133 halves build for G1, the seam chain itself takes seconds)
"""
import sys, os, itertools
from collections import defaultdict

sys.path.insert(0, "simulations")

RHO5 = (5, 4, 3, 2, 1)
DELTA = (2, 1, 0, -1, -2)
CUBE = [M for M in itertools.product((1, -1), repeat=5)]
ONE = (1, 1, 1, 1, 1)
MONE = (-1, -1, -1, -1, -1)
CORE = [M for M in CUBE if M not in (ONE, MONE)]
PAIRS_M = [M for M in CORE if M[0] == 1]          # one representative per ± pair
P_LETTERS = [tuple(1 if i == v else 0 for i in range(5)) for v in range(5)] + [ONE]
assert len(PAIRS_M) == 15 and len(P_LETTERS) == 6

FAIL = 0
def gate(name, ok, detail=""):
    global FAIL
    print("%s  %s%s" % (name, "OK " if ok else "FAIL", (" | " + detail) if detail else ""))
    if not ok:
        FAIL += 1

# ---------------------------------------------------------------- e-tables ----
def e_table_poly(letterpolys, maxdeg):
    T = [defaultdict(int) for _ in range(maxdeg + 1)]
    T[0][(0, 0, 0, 0, 0)] = 1
    for LP in letterpolys:
        for j in range(maxdeg, 0, -1):
            for k, c in list(T[j - 1].items()):
                for w, cw in LP.items():
                    T[j][tuple(x + y for x, y in zip(k, w))] += c * cw
    return [dict((k, c) for k, c in d.items() if c) for d in T]

ECOS = e_table_poly([{M: 1, tuple(-x for x in M): 1} for M in PAIRS_M], 15)
EP = e_table_poly([{p: 1} for p in P_LETTERS], 6)

perms5 = list(itertools.permutations(range(5)))
def _psgn(p):
    s = 1
    for i in range(5):
        for j in range(i + 1, 5):
            if p[i] > p[j]:
                s = -s
    return s
WPAIRS = [(_psgn(p) * s, p, eps)
          for p in perms5
          for eps, s in ((e, (lambda t: t[0] * t[1] * t[2] * t[3] * t[4])(e))
                         for e in itertools.product((1, -1), repeat=5))]

def convolve(A, B):
    out = defaultdict(int)
    small, big = (A, B) if len(A) < len(B) else (B, A)
    for k1, c1 in small.items():
        for k2, c2 in big.items():
            out[tuple(x + y for x, y in zip(k1, k2))] += c1 * c2
    return {k: c for k, c in out.items() if c}

def psi(i, t, nu):
    G = convolve(EP[i], ECOS[t])
    mu = tuple(v + r for v, r in zip(nu, RHO5))
    tot = 0
    for s, p, eps in WPAIRS:
        arg = tuple(eps[q] * mu[q] + 1 - DELTA[p[q]] for q in range(5))
        tot += s * G.get(arg, 0)
    return tot

def Phi_poly(nu, psi_bump=None):
    out = [0] * 18
    for i, place in ((0, 'a'), (1, 'b'), (2, 'c')):
        for t in range(16):
            v = psi(i, t, nu)
            if psi_bump is not None and psi_bump == (i, t):
                v += 1
            if not v:
                continue
            p = 15 - t
            if i == 0:
                out[p + 2] += v
                out[p] -= v
            elif i == 1:
                out[p + 1] += v
            else:
                out[p] += v
    return out

# Chebyshev S-polynomials: S_m(2cosθ) = sin((m+1)θ)/sinθ
S = [[1], [0, 1]]
for m in range(2, 20):
    row = [0] + S[m - 1]
    row = [(row[p] if p < len(row) else 0) - (S[m - 2][p] if p < len(S[m - 2]) else 0)
           for p in range(max(len(row), len(S[m - 2])))]
    S.append(row)
S10 = S[10]

def to_S_basis(A):
    A = list(A)
    b = {}
    for m in range(len(A) - 1, -1, -1):
        c = A[m] if m < len(A) else 0
        if c:
            b[m] = c
            for p, sc in enumerate(S[m]):
                A[p] -= c * sc
    return b

def poly_divmod(A, B):
    A = list(A)
    q = [0] * (max(len(A) - len(B) + 1, 1))
    db = len(B) - 1
    while any(A):
        da = max(p for p, a in enumerate(A) if a)
        if da < db:
            break
        c = A[da] // B[db]
        q[da - db] = c
        for p, bc in enumerate(B):
            A[p + da - db] -= c * bc
    return q, A

def poly_deg(A):
    return max((p for p, a in enumerate(A) if a), default=-1)

# ------------------------------------------------------- G1: the slice identity
print("G1: building the F133 halves for the slice spot-gate...", flush=True)
from f133_w_closed_form import build_halves, pack, PACKB
P1, P2 = build_halves()
_xc_cache = {}
def Xc(nu6):
    nu6 = tuple(nu6)
    v = _xc_cache.get(nu6)
    if v is not None:
        return v
    t = 0
    tp = pack(list(nu6)) + PACKB
    g = P1.get
    for k2, c2 in P2:
        w = g(tp - k2)
        if w:
            t += c2 * w
    _xc_cache[nu6] = t
    return t

ELAM36 = e_table_poly([{p: 1} for p in P_LETTERS[:5]]
                      + [{M: 1} for M in CUBE if M != MONE], 12)
A_DELTA = defaultdict(int)
for p in perms5:
    A_DELTA[tuple(DELTA[p[i]] for i in range(5))] += _psgn(p)

import random
random.seed(139)
bad = 0
ntest = 0
for mu1 in (6, 9, 11, 13, 16):
    d = 18 - mu1
    obj = defaultdict(int)
    for k1, c1 in A_DELTA.items():
        for k2, c2 in ELAM36[d].items():
            obj[tuple(x + y - 1 for x, y in zip(k1, k2))] += c1 * c2
    obj = {k: ((-1) ** d) * c for k, c in obj.items() if c}
    for _ in range(25):
        w = tuple(random.choice((-1, 1)) * random.randint(0, 7) for _ in range(5))
        ntest += 1
        if obj.get(w, 0) != Xc((2 * mu1,) + tuple(2 * x for x in w)):
            bad += 1
gate("G1", bad == 0, "slice identity vs Xc: %d/%d sampled coefficients match" % (ntest - bad, ntest))

# ------------------------------------------------------------- G2: the P-skew
NUS_ONEROW = [(k, 0, 0, 0, 0) for k in range(6)]
skew_bad = []
for nu in NUS_ONEROW:
    for i in range(7):
        for t in range(16):
            if psi(6 - i, t, nu) != -psi(i, t, nu):
                skew_bad.append((nu[0], i, t))
gate("G2", not skew_bad, "psi-skew whole grid (6 nu x 7 i x 16 t)%s"
     % ("" if not skew_bad else " first fail %s" % (skew_bad[0],)))

# --------------------------------------------------- G3: degree lemma + build
PHIS = {k: Phi_poly((k, 0, 0, 0, 0)) for k in range(6)}
gate("G3", all(poly_deg(PHIS[k]) <= 15 for k in range(6)),
     "deg Phi_k = %s (all <= 15)" % [poly_deg(PHIS[k]) for k in range(6)])

# --------------------------------------------------- G4: the a-priori table
HERE = os.path.dirname(os.path.abspath(__file__))
committed = {}
with open(os.path.join(HERE, 'results', 'f133_w_closed_form', 'chiC_coeffs.txt')) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts, nl = line.split(';')
        lam = tuple(int(x) for x in parts.split(',')) if parts else ()
        if len(lam) <= 2:
            committed[(lam[0] if lam else 0, lam[1] if len(lam) > 1 else 0)] = int(nl)

BS = {k: to_S_basis(PHIS[k]) for k in range(6)}
mism = []
cells = 0
for lam1 in range(0, 11):
    for k in range(0, min(lam1, 5) + 1):
        if (lam1 + k) % 2:
            continue
        cells += 1
        v = ((-1) ** lam1) * BS[k].get(lam1 + 5, 0)
        if v != committed.get((lam1, k), 0):
            mism.append((lam1, k, v, committed.get((lam1, k), 0)))
gate("G4", not mism, "a-priori table vs committed: %d/%d strip cells (incl. zeros)%s"
     % (cells - len(mism), cells, "" if not mism else " first %s" % (mism[0],)))

# --------------------------------------------------- G5: the seam division
# The committed F134 column polynomials P_k in the Chebyshev-S basis
# (P_0 = -(C10+C6) = S1-S5, P_1 = (t^4-t^-4)^2 = S4-S2-2S0, P_2 = 3C2 = 3S1,
#  P_3 = -(C4+3) = -(S2+2S0), P_4 = C2 = S1, P_5 = 0; C_a in t-units, z = t^2).
PK_S = {0: {1: 1, 5: -1},
        1: {0: -2, 2: -1, 4: 1},
        2: {1: 3},
        3: {0: -2, 2: -1},
        4: {1: 1},
        5: {}}
g5_ok = True
detail = []
for k in range(6):
    q, r = poly_divmod(PHIS[k], S10)
    bq = {m: c for m, c in to_S_basis(q).items() if c}
    want = {m: ((-1) ** k) * c for m, c in PK_S[k].items() if c}
    dr = poly_deg(r)
    ok = (bq == want) and (dr <= 4 + k) and (k != 0 or dr == -1)
    g5_ok = g5_ok and ok
    detail.append("k=%d degR=%d" % (k, dr))
gate("G5", g5_ok, "Q_k = (-1)^k P_k and deg R_k <= 4+k; k=0 remainder zero | " + " ".join(detail))

# --------------------------------------------------- G6: the fence
fence_ok = True
fdetail = []
for nu, lo, hi in (((2, 1, 0, 0, 0), None, 4), ((3, 1, 0, 0, 0), None, 4),
                   ((2, 2, 0, 0, 0), 8, None), ((3, 2, 0, 0, 0), 8, None),
                   ((4, 2, 0, 0, 0), 8, None)):
    _, r = poly_divmod(Phi_poly(nu), S10)
    dr = poly_deg(r)
    ok = (hi is None or dr <= hi) and (lo is None or dr >= lo)
    fence_ok = fence_ok and ok
    fdetail.append("nu=%s degR=%d" % (str(nu[:2]), dr))
gate("G6", fence_ok, "l=1 remainders small, l=2 remainders overflow | " + " ".join(fdetail))

# --------------------------------------------------- G7: corruption control
# (bump a PARITY-LIVE entry: for k = 0, i = 0 the live t are even; t = 4 feeds
#  the odd-degree part y^13 - y^11 that the k = 0 row actually reads)
Phi_bumped = Phi_poly((0, 0, 0, 0, 0), psi_bump=(0, 4))
bb = to_S_basis(Phi_bumped)
broken = any(((-1) ** lam1) * bb.get(lam1 + 5, 0) != committed.get((lam1, 0), 0)
             for lam1 in range(0, 11, 2))
gate("G7", broken, "bumped psi_{0,4} breaks the k=0 table reproduction")

print()
if FAIL:
    print("F139 gate: %d FAILURES" % FAIL)
    sys.exit(1)
print("F139 gate: ALL GREEN")
