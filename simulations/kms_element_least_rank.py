"""The least operator Schmidt rank of a thermal palindrome element on two sites
(experiments/THE_CUBE_SITE_BY_SITE.md, finding 1).

Setting: two sites, hopping (J/2)(XX + YY) with J = 1, optional (Δ/2) ZZ, cooling |0><1| at γ↓ = (3, 2) and heating
|1><0| at γ↑ = (1, 1), c = Σ(γ↓ + γ↑) = 7. An element is an invertible superoperator E with E·L·E⁻¹ = −L − c. The
Z-corner bound (Z⊗N ↦ ρ_ss, MirrorLocalityBoundClaim) gives Schmidt-rank(E) ≥ Schmidt-rank(ρ_ss) = 4 across the cut.

  A  exact: the element spaces (dim 36, 28 with ZZ; 18 keeping the coherence order d = n_ket − n_bra of
     THE_ONE_SQUARE §10), the split of the Schmidt rank (for an element that keeps d the realigned matrix is block
     diagonal by the order shift of the first site, so the rank is the sum of the block ranks), the steady state's split
     (1, 2, 1), also with ZZ, so the bound reads 2 + 1 + 1; F137's S_XY reflects L onto −L† − c, so U = E·S_XY⁻¹
     solves L·U = U·L† and Γ = U·T solves Lᵀ = Γ⁻¹LΓ, both carrying the trace onto ρ_ss; for M_β at one temperature
     Γ = Γ_β; across the gradient, without ZZ, no weight by ρ_ss (left, right, or the symmetric ρ^½·X·ρ^½, the last
     numerically) is one, with or without the transpose; the solution spaces of keeping and reversing d span 28 of 36
     (20 of 28 with ZZ), and with ZZ no element reverses d (L on the order-1 sector shares no eigenvalue with its
     conjugate).
  B  exact over Q(i), at this one rate point: no invertible element that keeps d (with and without ZZ), and none that
     reverses d (without ZZ), has Schmidt rank 4. Rank 4 forces R_{±2} = 0, rank R_{±1} ≤ 1, rank R_0 ≤ 2 (for an
     element reversing d, on its map U = E·S_XY⁻¹, which keeps d); with one Rabinowitsch variable per order sector and
     charts only over the two scalings that are symmetries (the global scale and ω^d for any nonzero complex ω), every
     Gröbner basis is {1}. Controls: the exact product element at one temperature satisfies every equation with every
     sector invertible; at the gradient the sectors d ≥ 0 alone are not empty; the full system was not decided on a case with
     solutions (sympy did not finish one).
  C  numerical, gated as a law: an element keeping d with Schmidt rank 5 (block ranks 0, 1, 3, 1, 0), embedded to 80
     digits; its rank tails follow the number of digits kept (20, 40, 80), the law a true solution follows and a near
     miss does not; a numerical point, not an exact one.
  D  the reading, without an explanation: on the maps U that keep d (the elements that reverse d) the rank-4
     equations of the sectors d ≥ 0 alone have six solutions modulo the scalings, with those sectors invertible
     (exact); at the rank-5 element R_0 has the rank 3 of its population block (its fourth singular value follows the
     digits kept, a law), and the hop coherence and the double coherences sit at the values that block dictates (a
     consistency reading, numerical).

Runtime about 25 minutes (sympy Gröbner over Q(i)).
"""
from __future__ import annotations

import itertools
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).parent))
sys.stdout.reconfigure(encoding='utf-8')

from framework import cube_tables as ct  # noqa: E402

R = sp.Rational
X, Y, Z = ct.PAULI['X'], ct.PAULI['Y'], ct.PAULI['Z']
QI = sp.QQ_I
FAILURES = []


def check(name, ok, detail=''):
    print(f"  {'PASS' if ok else 'FAIL'}: {name}" + (f"  [{detail}]" if detail else ''))
    if not ok:
        FAILURES.append(name)


def generator(down, up, J, delta):
    """L on vec_c(rho) (column-major, index i + 4j for |i><j|), two sites, exact over Q(i)."""
    on = lambda l, A: sp.kronecker_product(*[A if k == l else sp.eye(2) for k in range(2)])
    H = R(J, 2) * (on(0, X) * on(1, X) + on(0, Y) * on(1, Y)) + R(delta, 2) * on(0, Z) * on(1, Z)
    Id = sp.eye(4)
    L = -sp.I * (sp.kronecker_product(Id, H) - sp.kronecker_product(H.T, Id))
    for l in range(2):
        for A, rate in ((ct.LOWER, down[l]), (ct.RAISE, up[l])):
            Jm = on(l, sp.Matrix(A))
            JdJ = Jm.H * Jm
            L += rate * (sp.kronecker_product(Jm.conjugate(), Jm)
                         - R(1, 2) * (sp.kronecker_product(Id, JdJ) + sp.kronecker_product(JdJ.T, Id)))
    return L.applyfunc(sp.expand)


IDX = [(k % 4, k // 4) for k in range(16)]                      # vec index -> (row i, column j) of |i><j|
POP = lambda i: bin(i).count('1')
CHARGE = [POP(i) - POP(j) for i, j in IDX]                       # coherence order d = n_ket - n_bra (THE_ONE_SQUARE §10)
UNK = [(a, b) for a in range(16) for b in range(16) if CHARGE[a] == CHARGE[b]]
BIT = lambda x, s: (x >> (1 - s)) & 1                            # site s of a two-site index, s = 0 the high bit
UNITS = [(0, 0), (0, 1), (1, 0), (1, 1)]                         # one-site units (i, j): n0, q = |0><1|, p = |1><0|, n1
UNIT_NAME = {(0, 0): 'n0', (1, 1): 'n1', (1, 0): 'p', (0, 1): 'q'}
LOC = lambda u: u[0] - u[1]
OPS = [(o, i) for o in UNITS for i in UNITS]                    # one-site superoperator entries (out, in)
SHIFT = lambda op: LOC(op[0]) - LOC(op[1])
OP_NAME = lambda op: UNIT_NAME[op[1]] + '->' + UNIT_NAME[op[0]]


def unit(k, s):
    i, j = IDX[k]
    return (BIT(i, s), BIT(j, s))


def element_space(L, c, frame):
    """Charge-preserving solutions over Q(i) as coefficient vectors on UNK (canonical rref nullspace basis).
    frame 'E': E·L + (L + c)·E = 0;  frame 'U': L·U − U·L† = 0."""
    pos = {u: q for q, u in enumerate(UNK)}
    Ld = L.H
    rows = []
    for a in range(16):
        for b in range(16):
            row = [0] * len(UNK)
            for k in range(16):
                if frame == 'E':
                    if (a, k) in pos and L[k, b] != 0:
                        row[pos[(a, k)]] += L[k, b]
                    if (k, b) in pos and L[a, k] != 0:
                        row[pos[(k, b)]] += L[a, k]
                else:
                    if (k, b) in pos and L[a, k] != 0:
                        row[pos[(k, b)]] += L[a, k]
                    if (a, k) in pos and Ld[k, b] != 0:
                        row[pos[(a, k)]] -= Ld[k, b]
            if frame == 'E' and (a, b) in pos:
                row[pos[(a, b)]] += c
            if any(x != 0 for x in row):
                rows.append(row)
    N = DomainMatrix([[QI.from_sympy(sp.sympify(x)) for x in rw] for rw in rows], (len(rows), len(UNK)), QI)
    N = N.nullspace().to_Matrix()
    return [[N[q, p] for p in range(len(UNK))] for q in range(N.rows)]


def full_space_dim(L, c):
    """dim of all solutions of E·L + (L + c)·E = 0 (no charge restriction): 256 − rank of the Sylvester operator."""
    S = sp.kronecker_product(sp.eye(16), L.T) + sp.kronecker_product(L + c * sp.eye(16), sp.eye(16))
    return 256 - DomainMatrix.from_Matrix(S).convert_to(QI).rank()


def as_matrix(basis, coeffs):
    E = sp.zeros(16, 16)
    for (a, b), col in zip(UNK, zip(*basis)):
        E[a, b] = sum(cf * x for cf, x in zip(coeffs, col) if x != 0)
    return E


def realigned_blocks(get):
    """Blocks of the realigned matrix by the charge shift D of site 1; get(a, b) returns the entry E[a, b]."""
    Rm = {}
    for a in range(16):
        for b in range(16):
            Rm[((unit(a, 0), unit(b, 0)), (unit(a, 1), unit(b, 1)))] = (a, b)
    blocks = {}
    for D in range(-2, 3):
        rows = [op for op in OPS if SHIFT(op) == D]
        cols = [op for op in OPS if SHIFT(op) == -D]
        blocks[D] = [[get(*Rm[(r, cl)]) if (r, cl) in Rm else 0 for cl in cols] for r in rows]
    return blocks


def sector(m):
    return [k for k in range(16) if CHARGE[k] == m]


def is_zero(M):
    return M.applyfunc(sp.expand).is_zero_matrix is True


def qi_rank(M):
    return DomainMatrix.from_Matrix(sp.Matrix(M).applyfunc(sp.expand)).convert_to(QI).rank()


def build_system(basis):
    """Rank-4 equations and sector determinants over Q(i), on the symbols t_0..t_{n-1}."""
    n = len(basis)
    t = sp.symbols(f't0:{n}')
    forms = {u: [basis[q][p] for q in range(n)] for p, u in enumerate(UNK)}
    lin = lambda fv: sp.Poly(sum(x * tq for x, tq in zip(fv, t) if x != 0), *t, domain=QI)
    zero = [0] * n
    blk = realigned_blocks(lambda a, b: forms.get((a, b), zero))
    supp = lambda fv: {q for q, x in enumerate(fv) if x != 0}

    def det_forms(M):
        if len(M) == 1:
            return lin(M[0][0])
        tot = sp.Poly(0, *t, domain=QI)
        for j in range(len(M)):
            if any(M[0][j]):
                tot += lin(M[0][j]) * det_forms([row[:j] + row[j + 1:] for row in M[1:]]) * (1 if j % 2 == 0 else -1)
        return tot

    eqs = [lin(fv) for D in (2, -2) for row in blk[D] for fv in row if any(fv)]
    for D, r in ((1, 2), (-1, 2), (0, 3)):
        M = blk[D]
        for ri in itertools.combinations(range(len(M)), r):
            for ci in itertools.combinations(range(len(M[0])), r):
                p = det_forms([[M[i][j] for j in ci] for i in ri])
                if not p.is_zero:
                    eqs.append(p)
    eqs = list({p.as_expr(): p for p in eqs}.values())
    secvars = {m: sorted({q for a in sector(m) for b in sector(m) for q in supp(forms[(a, b)])}) for m in range(-2, 3)}
    dets = {m: det_forms([[forms[(a, b)] for b in sector(m)] for a in sector(m)]) for m in range(-2, 3)}
    return t, eqs, secvars, dets, blk


def charts(group):
    for j in range(len(group)):
        yield {group[i]: 0 for i in range(j)} | {group[j]: 1}


def decide(t, eqs, secvars, dets, stages=('A', 'B')):
    """Gröbner over Q(i) per chart: stage A = equations in the sectors d = 0, 1, 2 only; stage B = all.
    Only the two symmetries are normalised: the global scale (a chart in d = 0) and the phase (a chart in d = 1)."""
    in_a = set(secvars[0]) | set(secvars[1]) | set(secvars[2])
    is_a = lambda e: {int(str(g)[1:]) for g in e.as_expr().free_symbols} <= in_a
    zs = sp.symbols('z_m2 z_m1 z_0 z_1 z_2')
    results = []
    groups = [[t[q] for q in secvars[m]] for m in (0, 1)]
    for c0, c1 in itertools.product(charts(groups[0]), charts(groups[1])):
        sub = c0 | c1
        ds = {m: sp.expand(dets[m].as_expr().subs(sub)) for m in range(-2, 3)}
        if any(d == 0 for d in ds.values()):
            continue
        za = [z * ds[m] - 1 for z, m in zip(zs, range(-2, 3)) if m >= 0]
        zb = [z * ds[m] - 1 for z, m in zip(zs, range(-2, 3)) if m < 0]
        ea = [e for e in (sp.expand(p.as_expr().subs(sub)) for p in eqs if is_a(p)) if e != 0]
        eb = [e for e in (sp.expand(p.as_expr().subs(sub)) for p in eqs if not is_a(p)) if e != 0]
        gens = sorted(set().union(*[e.free_symbols for e in ea + eb + za + zb]), key=str)
        gens_a = sorted(set().union(*[e.free_symbols for e in ea + za]), key=str)
        GA = sp.groebner(ea + za, *gens_a, order='grevlex', domain=QI)
        empty_a = list(GA.exprs) == [1]
        empty_b = empty_a or list(sp.groebner(list(GA.exprs) + eb + zb, *gens, order='grevlex', domain=QI).exprs) == [1]
        results.append(({str(k): v for k, v in sub.items()}, empty_a, empty_b, GA, gens))
    return results


def main():
    t_start = time.time()
    L0 = generator((3, 2), (1, 1), 1, 0)
    L1 = generator((3, 2), (1, 1), 1, 1)

    print("\nA  the element spaces and the local-charge split")
    check("all elements: dim 36 without ZZ, 28 with ZZ",
          (full_space_dim(L0, 7), full_space_dim(L1, 7)) == (36, 28))
    basis_e = element_space(L0, 7, 'E')
    check("elements keeping the coherence order d: dim 18", len(basis_e) == 18)
    rng = np.random.default_rng(1)
    coeffs = [sp.Rational(int(rng.integers(-9, 10)), int(rng.integers(1, 9))) + sp.I * int(rng.integers(-5, 6))
              for _ in range(len(basis_e))]
    E = as_matrix(basis_e, coeffs)
    check("a random element that keeps d is one: E·L + (L + c)·E = 0", is_zero(E * L0 + (L0 + 7 * sp.eye(16)) * E))
    blk = realigned_blocks(lambda a, b: E[a, b])
    full = sp.zeros(16, 16)
    for a in range(16):
        for b in range(16):
            r = OPS.index((unit(a, 0), unit(b, 0)))
            cl = OPS.index((unit(a, 1), unit(b, 1)))
            full[r, cl] = E[a, b]
    off = [full[r, cl] for r, op1 in enumerate(OPS) for cl, op2 in enumerate(OPS) if SHIFT(op1) != -SHIFT(op2)]
    check("its realigned matrix vanishes off the blocks pairing order shift δ on site 1 with −δ on site 2, "
          "so the Schmidt rank is the sum of the block ranks", all(sp.expand(v) == 0 for v in off),
          f"{len(off)} off-block entries")
    ns = DomainMatrix.from_Matrix(L0).convert_to(QI).nullspace().to_Matrix()
    check("the steady state is unique (kernel of L of dimension 1)", ns.rows == 1, f"{ns.rows}")
    rho = sp.Matrix(4, 4, lambda i, j: ns[0, i + 4 * j])
    rho = rho / rho.trace()
    split = {}
    for d in (-1, 0, 1):
        rows = [(i1, j1) for i1 in (0, 1) for j1 in (0, 1) if i1 - j1 == d]
        cols = [(i2, j2) for i2 in (0, 1) for j2 in (0, 1) if i2 - j2 == -d]
        split[d] = qi_rank([[rho[2 * i1 + i2, 2 * j1 + j2] for (i2, j2) in cols] for (i1, j1) in rows])
    check("the steady state splits (1, 2, 1) by the order of site 1: the bound is 2 + 1 + 1 on R_0, R_{±1}",
          split == {-1: 1, 0: 2, 1: 1}, f"{split}")
    XX, YY = sp.kronecker_product(X, X), sp.kronecker_product(Y, Y)
    Sxy = sp.kronecker_product(YY.T, XX)                          # vec(X⊗X·ρ·Y⊗Y) = (Y⊗Y)ᵀ ⊗ (X⊗X)·vec(ρ)
    check("S_XY (ρ ↦ X⊗X·ρ·Y⊗Y) carries L onto −L† − c (F137), the product reflection, without ZZ and with ZZ",
          all(is_zero(Sxy * L * Sxy.inv() + L.H + 7 * sp.eye(16)) for L in (L0, L1)))
    zero_sector = lambda basis: sum(1 for vec in basis if all(CHARGE[a] == 0 for (a, b), x in zip(UNK, vec) if x != 0))
    flat = lambda M: [M[a, b] for a in range(16) for b in range(16)]
    span = {}
    for name, L in (("without ZZ", L0), ("with ZZ", L1)):
        keep, rev = element_space(L, 7, 'E'), element_space(L, 7, 'U')
        unit_vec = lambda n, q: [1 if k == q else 0 for k in range(n)]
        vecs = [flat(as_matrix(keep, unit_vec(len(keep), q))) for q in range(len(keep))]
        vecs += [flat(as_matrix(rev, unit_vec(len(rev), q)) * Sxy) for q in range(len(rev))]
        span[name] = (len(keep), len(rev), zero_sector(keep), qi_rank(vecs))
    check("keeping and reversing d do not exhaust the elements: their solution spaces share the order-zero sector and span 28 of 36, "
          "with ZZ 20 of 28 (exact joint rank)", span == {"without ZZ": (18, 18, 8, 28), "with ZZ": (18, 10, 8, 20)},
          f"{span}")
    xs = sp.Symbol('x')
    s1 = sector(1)
    gcd_deg, touches = {}, {}
    for name, L in (("without ZZ", L0), ("with ZZ", L1)):
        p = sp.Poly(L.extract(s1, s1).charpoly(xs).as_expr().expand(), xs, extension=sp.I)
        q = sp.Poly(L.extract(s1, s1).conjugate().charpoly(xs).as_expr().expand(), xs, extension=sp.I)
        gcd_deg[name] = sp.gcd(p, q).degree()
        touches[name] = any(abs(CHARGE[a]) == 1 for vec in element_space(L, 7, 'U') for (a, b), x in zip(UNK, vec) if x != 0)
    check("with ZZ no element reverses d: L on the order-1 sector shares no eigenvalue with its conjugate, so the maps U "
          "that keep d vanish on the orders ±1 and every U·S_XY is singular; without ZZ that spectrum is conjugation-closed",
          gcd_deg == {"without ZZ": 4, "with ZZ": 0} and touches == {"without ZZ": True, "with ZZ": False},
          f"gcd degrees {gcd_deg}, U touches d = ±1 {touches}")
    U = E * Sxy.inv()
    ivec = sp.Matrix(16, 1, lambda k, _: 1 if IDX[k][0] == IDX[k][1] else 0)
    rvec = sp.Matrix(16, 1, lambda k, _: rho[IDX[k][0], IDX[k][1]])
    Ui = U * ivec
    check("so U = E·S_XY⁻¹ is a metric, L·U = U·L†, and it carries the trace onto the steady state",
          is_zero(L0 * U - U * L0.H) and qi_rank(sp.Matrix.hstack(Ui, rvec)) == 1)
    T = sp.zeros(16, 16)
    for i in range(4):
        for j in range(4):
            T[j + 4 * i, i + 4 * j] = 1                            # vec(ρᵀ)
    G = U * T
    check("L† = T·Lᵀ·T, so Γ = U·T is a weight of PROOF_ANTILINEAR_TRIANGLE's form, L·Γ = Γ·Lᵀ, with Γ(I) ∝ ρ_ss",
          is_zero(L0.H - T * L0.T * T) and is_zero(L0 * G - G * L0.T) and qi_rank(sp.Matrix.hstack(G * ivec, rvec)) == 1)
    sand = lambda A, B: sp.kronecker_product(B.T, A)                # ρ ↦ A·ρ·B
    Lone = generator((3, 3), (1, 1), 1, 0)
    gibbs = sp.diag(1, sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 9)) * sp.Rational(9, 16)
    for side, (A, B) in (("left", (gibbs, sp.eye(4))), ("right", (sp.eye(4), gibbs))):
        W = sand(A, B)
        check(f"one temperature: the {side} weight by the Gibbs product is a metric composed with the transpose and not without it "
              "(the split is free at one ratio; the symmetric split gives M_β)", is_zero(Lone * W * T - W * T * Lone.H) and not is_zero(Lone * W - W * Lone.H))
    for side, (A, B) in (("left", (rho, sp.eye(4))), ("right", (sp.eye(4), rho))):
        W = sand(A, B)
        check(f"gradient: the {side} weight by ρ_ss is no metric, with or without the transpose",
              not is_zero(L0 * W * T - W * T * L0.H) and not is_zero(L0 * W - W * L0.H))
    from scipy.linalg import sqrtm
    num = lambda M: np.array(M.evalf(30).tolist(), dtype=complex)

    Tn = np.array(T.tolist(), dtype=float)

    def symmetric_residual(L, state, with_t):
        r = sqrtm(num(state))
        G = np.kron(r.T, r)                                         # ρ^½·X·ρ^½, the proposition's own form
        G = G @ Tn if with_t else G
        Ln = num(L)
        return np.linalg.norm(Ln @ G - G @ Ln.T) / (np.linalg.norm(Ln) * np.linalg.norm(G))
    one = symmetric_residual(Lone, gibbs, False)
    grad = [symmetric_residual(L0, rho, w) for w in (False, True)]
    check("the proposition's symmetric weight ρ^½·X·ρ^½ (numerical): a weight at one temperature (control); across the "
          "gradient no weight, with the transpose or without it", one < 1e-13 and min(grad) > 1e-3,
          f"relative residual {one:.1e} at one ratio, {grad[0]:.2e} and {grad[1]:.2e} across the gradient")
    ns1 = DomainMatrix.from_Matrix(L1).convert_to(QI).nullspace().to_Matrix()
    rho1 = sp.Matrix(4, 4, lambda i, j: ns1[0, i + 4 * j])
    rho1 = rho1 / rho1.trace()
    split1 = {}
    for d in (-1, 0, 1):
        rows = [(i1, j1) for i1 in (0, 1) for j1 in (0, 1) if i1 - j1 == d]
        cols = [(i2, j2) for i2 in (0, 1) for j2 in (0, 1) if i2 - j2 == -d]
        split1[d] = qi_rank([[rho1[2 * i1 + i2, 2 * j1 + j2] for (i2, j2) in cols] for (i1, j1) in rows])
    check("with ZZ the steady state splits (1, 2, 1) as well, so the same pattern is forced",
          split1 == {-1: 1, 0: 2, 1: 1}, f"{split1}")

    print("\nB  Schmidt rank 4 excluded for elements that keep d and for those that reverse it (Gröbner over Q(i))")
    for name, L, frame in (("elements keeping d", L0, 'E'), ("elements keeping d, with ZZ", L1, 'E'),
                           ("maps U keeping d, i.e. elements U·S_XY reversing d", L0, 'U')):
        t0 = time.time()
        t, eqs, secvars, dets, _ = build_system(element_space(L, 7, frame))
        res = decide(t, eqs, secvars, dets)
        survivors = [r[0] for r in res if not r[1]]
        check(f"{name}: every chart empty for the full system ({len(res)} charts, {len(eqs)} equations)",
              len(res) > 0 and all(r[2] for r in res), f"{time.time() - t0:.0f} s")
        check(f"{name}: anti-vacuity, the sectors d ≥ 0 alone are not empty", len(survivors) >= 1, f"surviving chart {survivors}")
    t, eqs, secvars, dets, _ = build_system(element_space(generator((3, 3), (1, 1), 1, 0), 8, 'E'))
    site = {((0, 0), (1, 1)): -3, ((1, 1), (0, 0)): 1, ((1, 0), (1, 0)): 1, ((0, 1), (0, 1)): -3}
    Eb = sp.zeros(16, 16)
    for a in range(16):
        for b in range(16):
            Eb[a, b] = site.get((unit(a, 0), unit(b, 0)), 0) * site.get((unit(a, 1), unit(b, 1)), 0)
    Lb = generator((3, 3), (1, 1), 1, 0)
    check("control: at one temperature (3, 3 | 1, 1) the product map n0 ↦ n1, n1 ↦ −3·n0, p ↦ p, q ↦ −3·q on each site is an element",
          is_zero(Eb * Lb + (Lb + 8 * sp.eye(16)) * Eb) and Eb.det() != 0)
    basis_b = element_space(Lb, 8, 'E')
    Mx = sp.Matrix([[basis_b[q][p] for q in range(len(basis_b))] for p in range(len(UNK))])
    x, params = Mx.gauss_jordan_solve(sp.Matrix([Eb[a, b] for (a, b) in UNK]))
    sub = {t[q]: x[q] for q in range(len(basis_b))}
    check("control: it satisfies every rank-4 equation and every sector determinant is nonzero",
          params.shape[0] == 0 and all(sp.simplify(e.as_expr().subs(sub)) == 0 for e in eqs)
          and all(sp.simplify(d.as_expr().subs(sub)) != 0 for d in dets.values()), f"{len(eqs)} equations")

    print("\nC  an element keeping d with Schmidt rank 5, gated as a law")
    rank5(basis_e)

    print("\nD  the reading, on the metrics U that keep d: the half d ≥ 0 is rigid; R_0 at the rank-5 point")
    reading(L0, basis_e)

    print(f"\n{'ALL PASS' if not FAILURES else 'FAILURES: ' + ', '.join(FAILURES)}  ({time.time() - t_start:.0f} s)")
    return 0 if not FAILURES else 1


# Rank-5 point, 80 digits, coordinates on the frame-E basis of element_space(L0, 7, 'E') scaled by s_q = max(|Re|+|Im|)
# of basis vector q (scale and phase fixed by t_0 = t_7 = 1, the curve parameter by t_12 = 5 + 10i).
RANK5_POINT = [
    ('1.0', '+0.0'),
    ('1.2918301945742049446132930985034906313011985028298927830171861120471744305718407', '+0.010518073675159321140912530706369308653965637297866625617248351550407306489370833'),
    ('0.0027744308306794899624378476238016765126764109399318577416118750265123528947005208', '-0.34075569638706385892852531038946521438739183670868717371450136136518102316957577'),
    ('0.25389423496143464691818733747895159605747375143021231789229889851521496847558232', '-7.2717934137679169045298843383677308683933563890668171094425149383079577979804719'),
    ('-0.19174884112586206380701675316917055625146836013663366444804206609005711050441042', '-0.0031393814507482127254075380443817821461532984661042145539469385806263224317075495'),
    ('0.86975890574146579351054718413027640110930218602986797180716460644083646462073463', '+0.014240007706307565795013849350347187329343581492846692302903967765053498863104754'),
    ('0.051174622473901403145059919971510455254007935206425300963885427920349698091792565', '-1.2955166003594648969647434634957097060179035383898091692103169441516356669716471'),
    ('1.0', '+0.0'),
    ('0.082484990250714348396196691929956506716592842925270492501558898343620865802336637', '-4.0077739840107732271770505773238339323885227631920813998820634016886558090739771'),
    ('-0.23137784942458552562063382512929941224397896382387290048314386675837186986361962', '+0.000020192632895407363797457150920931833273087661429931392838742680304633041682392227'),
    ('-0.25575866077068218686149094475401388303520006286526407640814782257947059038018479', '+0.000022320376646281732642790406961795987821412169622983989973168762192333683973213968'),
    ('0.018338221130151233849192807537876844458849573086821793604839460637869648670053244', '-0.79638956587050529158160229484680008654705415813325997884096547626703083947871046'),
    ('5.0', '+10.0'),
    ('6.9993976127370043795080265065033403607837358644789017964734151732796213232809358', '+8.6702582861272962141727253382198373962409887471046432180955019429566803275449067'),
    ('0.0072538961105744671230503233086665595236929923011504335407697904477666129747061768', '-0.89092378636552375502445300586412774845002347667751767652063976559421794480198339'),
    ('0.0055880820477557516545533635323592059051498246855003970546487161213975639442680797', '-0.68632844206993897452542031546526945783290808072817533668589950515205058753808752'),
    ('5.7474143860420068243179627777538800673434755535765127590704668099299722560645807', '+7.1194079778676962016890052813907865105450584048125580995542815660727194589789139'),
    ('2.6284447992587344268439566490155968319119606090775610706322143358373762069591328', '+3.4085704874064245740530341661603142083075213823756997946731509100980536451693895'),
]


def point_elements(basis, digits_list=(20, 40, 80)):
    """The rank-5 point as a 16×16 mpmath matrix, the coordinates cut to each number of digits kept."""
    mp.mp.dps = 100
    scale = [max(abs(sp.re(v)) + abs(sp.im(v)) for v in vec if v != 0) for vec in basis]
    to_mp = lambda v: mp.mpc(str(sp.re(sp.N(v, 110))), str(sp.im(sp.N(v, 110))))
    cols = [[to_mp(x / scale[q]) if x != 0 else None for q, x in enumerate(col)] for col in zip(*basis)]
    out = {}
    for digits in digits_list:
        u = [mp.mpc(mp.mpf(re_[:digits + 3]), mp.mpf(im_[:digits + 3])) for re_, im_ in RANK5_POINT]
        E = mp.matrix(16, 16)
        for (a, b), col in zip(UNK, cols):
            E[a, b] = mp.fsum(u[q] * w for q, w in enumerate(col) if w is not None)
        out[digits] = E
    return out


def mp_sv(M):
    return sorted([abs(s) for s in mp.svd_c(mp.matrix(M), compute_uv=False)], reverse=True)


def rank5(basis):
    tails = {}
    for digits, E in point_elements(basis).items():
        blk = realigned_blocks(lambda a, b: E[a, b])
        sv = mp_sv
        s1, sm1, s0 = sv(blk[1]), sv(blk[-1]), sv(blk[0])
        top = max(s0[0], s1[0], sm1[0])
        tails[digits] = max(s1[1], sm1[1], s0[3], abs(blk[2][0][0]), abs(blk[-2][0][0])) / top
        if digits == 80:
            cond = np.linalg.cond(np.array(E.tolist(), dtype=complex))
            kept = min(float(s0[2] / top), float(s1[0] / top), float(sm1[0] / top))
    check("block ranks (0, 1, 3, 1, 0): the tails follow the digits kept (log10 tail within 2.5 of minus the digits)",
          all(abs(float(mp.log10(tails[d])) + d) < 2.5 for d in tails), ', '.join(f"{d} digits: {float(tails[d]):.1e}" for d in tails))
    check("the kept singular values stay O(1) and the element is invertible", kept > 1e-3 and cond < 100,
          f"cond {cond:.1f}, smallest kept {kept:.2e}")


def reading(L, basis_e):
    basis = element_space(L, 7, 'U')
    t, eqs, secvars, dets, _ = build_system(basis)
    res = decide(t, eqs, secvars, dets)
    survivors = [r for r in res if not r[1]]
    check("on the maps U that keep d, the half d ≥ 0 survives in exactly one chart", len(survivors) == 1,
          f"{survivors[0][0] if survivors else None}")
    GA = survivors[0][3]
    GL = GA.fglm('lex').exprs
    last = GL[-1]
    xv = list(last.free_symbols)[0]
    facs = sp.factor_list(last, xv, extension=sp.I)[1]
    shape = all(len(g.free_symbols - {xv}) == 1 and sp.degree(g, list(g.free_symbols - {xv})[0]) == 1 for g in GL[:-1])
    squarefree = sp.degree(sp.gcd(last, sp.diff(last, xv)), xv) == 0
    check("it is rigid: a lex basis in shape position over a squarefree eliminant of degree 6, three quadratic factors "
          "(exactly six solutions modulo the scalings)",
          shape and squarefree and sp.degree(last, xv) == 6 and sorted(sp.degree(f, xv) for f, _ in facs) == [2, 2, 2],
          f"{len(GL)} generators")
    names = [OP_NAME(op) for op in OPS if SHIFT(op) == 0]
    pops = [k for k, nm in enumerate(names) if 'p' not in nm and 'q' not in nm]
    a, b = names.index('p->p'), names.index('q->q')
    pop_sv = {}
    for digits, Em in point_elements(basis_e).items():
        blk0 = realigned_blocks(lambda r, s: Em[r, s])[0]
        s = mp_sv([[blk0[r][cl] for cl in pops] for r in pops])
        pop_sv[digits] = [x / s[0] for x in s]
    check("at the rank-5 element the population block of R_0 has rank 3: the third singular value stays O(10⁻⁴), "
          "the fourth follows the digits kept (log10 within 2.5 of minus the digits)",
          all(pop_sv[d][2] > 1e-6 and abs(float(mp.log10(pop_sv[d][3])) + d) < 2.5 for d in pop_sv),
          f"80 digits: {', '.join(mp.nstr(x, 3) for x in pop_sv[80])}; fourth at 20, 40: "
          f"{mp.nstr(pop_sv[20][3], 2)}, {mp.nstr(pop_sv[40][3], 2)}")
    E = np.array(point_elements(basis_e, (40,))[40].tolist(), dtype=complex)
    R0 = np.array(realigned_blocks(lambda r, s: E[r, s])[0], dtype=complex)
    Kp = np.linalg.pinv(R0[np.ix_(pops, pops)], rcond=1e-10)
    need = lambda r, cl: (R0[np.ix_([r], pops)] @ Kp @ R0[np.ix_(pops, [cl])])[0, 0]
    ratios = [R0[r, cl] / need(r, cl) for r, cl in ((a, b), (b, a), (a, a), (b, b))]
    check("consistency: so the hop coherence p⊗q, its mirror and the double coherences take the values that block dictates",
          all(abs(x - 1) < 1e-8 for x in ratios), f"ratios {[complex(np.round(x, 10)) for x in ratios]}")


if __name__ == '__main__':
    sys.exit(main())
