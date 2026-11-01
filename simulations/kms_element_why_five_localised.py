"""Why five, localised: where the fifth term of a thermal palindrome element is forced
(experiments/THE_CUBE_SITE_BY_SITE.md, finding 1; the second route to ThermalElementLeastRankWitness.Localisation).

Setting as in kms_element_least_rank.py: two sites, H = J·(σ⁺σ⁻ + σ⁻σ⁺) with J = 1 (no ZZ), cooling |0><1| at
γ↓ = (3, 2) and heating |1><0| at γ↑ = (1, 1), c = 7. The Z-corner bound reads 2 + 1 + 1 = 4 on the blocks of the first
site's order shift δ; no element of Schmidt rank 4 exists (exact, that script and the witness); one of rank 5 is found
numerically with block ranks (0, 1, 3, 1, 0). This script asks WHERE rank 4 is refused.

  A  the classical shadow, exact over ℚ: the 2-bit Markov chain with the same two baths and the hop as a symmetric
     rate k between 10 and 01 (an irreversible chain, the cycle 00 → 10 → 01 → 00 carrying γ↑₀γ↓₁ : γ↑₁γ↓₀ = 2 : 3).
     Its transposer G (L·G = G·Lᵀ, G invertible) has least operator Schmidt rank 2 across the two bits, the rank of
     its steady state; 1 at one temperature with the hop, 1 at the gradient without it (two cases). The bound is the
     identity route, Lᵀ·1 = 0 so G·1 ∝ π, the image of Γ(I) ∝ ρ_ss. So the classical cycle alone forces nothing
     beyond the state; that the quantum surplus 4 → 5 is carried by the coherent current the rate model lacks is a
     reading by analogy, the chain being no restriction of the quantum problem.
  B  the point controls, exact: at one temperature the product weight Γ_β = π^½(·)π^½ solves L·Γ = Γ·Lᵀ_comp
     (Lᵀ_comp = T·L†·T, the Γ frame, E keeping d) and U_β = Γ_β·T solves L·U = U·L† (the U frame, E reversing d),
     with block ranks (0, 0, 1, 0, 0) (the U-frame check follows from the Γ-frame one by T² = 1, so it is a check
     of the code paths, not an independent control); at the gradient neither that product nor the gradient's own
     product weight π₀^½ ⊗ π₁^½ (·) π₀^½ ⊗ π₁^½ solves the Γ frame.
  C  the sector chain, modulo one prime p ≡ 1 (mod 4) with i = √−1 mod p (a scout, the exact decision is the
     witness's replay over ℚ(i)): the intertwining equations decouple by the coherence order d, so a set S of
     sectors can be imposed alone; the rank-4 shape (R_±2 = 0, rank R_±1 ≤ 1, rank R_0 ≤ 2) is asked only of the
     minors whose entries all lie in sectors of S, with every sector of S invertible (Rabinowitsch). The witness keeps
     the larger set of the equations whose VARIABLES lie in S (a minor with out-of-S entries that cancel is kept
     there, dropped here) and replays over ℚ(i), so its EMPTY is the proof and this script's EMPTY implies it, while a
     point here is implied by a point there; the witness also decides the mirror half d ≤ 0 and the case with ZZ,
     which this script does not; the mirror X ↦ T·X̄·T (L keeps Hermiticity, so T·L·T = L̄) carries the solutions
     of S onto those of −S and swaps the blocks δ ↔ −δ, which the shape respects. In both frames, without ZZ:
       S = {0}            the quiet edge with the hop coherence |10><01| and its mirror:  solutions
       S = {0, +1, +2}    one one-lit sector with the quiet edge (the witness's stage A; its mirror
                          S = {0, −1, −2} is the witness's): solutions
       S = {0, +1, −1}    both one-lit sectors with the quiet edge:                       EMPTY
       S = all            the full system (exact in the witness):                         EMPTY
     The fifth term is forced by the pair d = +1 and d = −1 together with d = 0: the cells with one site lit,
     read as σ⁺ or σ⁻ (d is dual to the lit site's choice between X and Y, THE_ONE_SQUARE §10); the double
     coherences are not needed for the refusal.

Runtime: A and B seconds; C about an hour (sympy's Gröbner modulo p; the 16- and 18-parameter runs take 5 to 30
minutes each). `--quick` runs A, B and the S = {0} rows only.
"""
from __future__ import annotations

import sys
import time
from itertools import combinations

import sympy as sp

sys.stdout.reconfigure(encoding='utf-8')
FAILURES: list[str] = []
P = 1000033            # ≡ 1 (mod 4)


def check(name, ok, detail=''):
    print(f"  {'PASS' if ok else 'FAIL'}: {name}" + (f"  [{detail}]" if detail else ''))
    if not ok:
        FAILURES.append(name)


# ----- A: the classical shadow ------------------------------------------------------------------------------------

def classical_chain(gu0, gd0, gu1, gd1, k):
    """Column-stochastic generator on the states (a, b) ↦ 2a + b: site baths and the hop 10 ↔ 01 at rate k."""
    L = sp.zeros(4, 4)

    def add(i, j, r):
        L[j, i] += r
        L[i, i] -= r
    for b in (0, 1):
        add(0*2 + b, 1*2 + b, gu0); add(1*2 + b, 0*2 + b, gd0)
    for a in (0, 1):
        add(a*2 + 0, a*2 + 1, gu1); add(a*2 + 1, a*2 + 0, gd1)
    add(2, 1, k); add(1, 2, k)
    return L


def realign4(G):
    """G[(a', b'), (a, b)] ↦ R[(a', a), (b', b)]: the operator Schmidt rank across the bits is rank R."""
    R = sp.zeros(4, 4)
    for a in (0, 1):
        for b in (0, 1):
            for a2 in (0, 1):
                for b2 in (0, 1):
                    R[2*a2 + a, 2*b2 + b] = G[2*a2 + b2, 2*a + b]
    return R


def classical_least_rank(L):
    """The least Schmidt rank of an invertible transposer, exact over ℚ: minors + Rabinowitsch on det G."""
    syms = sp.symbols('g0:16')
    G = sp.Matrix(4, 4, syms)
    (vec,) = sp.linsolve(list(L*G - G*L.T), syms)
    free = sorted(set().union(*[v.free_symbols for v in vec]), key=str)
    G = sp.Matrix(4, 4, list(vec))
    R = realign4(G)
    detG = sp.expand(G.det())
    z = sp.Symbol('z')
    for r in (1, 2, 3):
        minors = [sp.expand(R.extract(list(rows), list(cols)).det())
                  for rows in combinations(range(4), r + 1) for cols in combinations(range(4), r + 1)]
        minors = [m for m in minors if m != 0]
        if not minors:
            return r, len(free)
        gb = sp.groebner(minors + [z*detG - 1], *free, z, order='grevlex', domain='QQ')
        if list(gb.exprs) != [1]:
            return r, len(free)
    return 4, len(free)


def part_a():
    print("A. the classical shadow (exact over ℚ)")
    for label, args, expected in (
        ('one temperature (3,3)/(1,1), hop 1', (1, 3, 1, 3, 1), 1),
        ('gradient (3,2)/(1,1), hop 1', (1, 3, 1, 2, 1), 2),
        ('gradient (3,2)/(1,1), hop 1/7', (1, 3, 1, 2, sp.Rational(1, 7)), 2),
        ('gradient (3,2)/(1,1), no hop', (1, 3, 1, 2, 0), 1),
    ):
        L = classical_chain(*args)
        pi = L.nullspace()[0]
        pi = pi / sum(pi)
        state_rank = sp.Matrix(2, 2, list(pi)).rank()
        rank, dim = classical_least_rank(L)
        check(f"{label}: least transposer rank {rank} = steady-state rank {state_rank} = {expected}",
              rank == expected and state_rank == expected, f"solution space dim {dim}, π = {list(pi)}")


# ----- the quantum generator ----------------------------------------------------------------------------------------

def kron(a, b):
    return sp.kronecker_product(a, b)


def superop_of(f, dim=4):
    cols = []
    for idx in range(dim*dim):
        a, b = divmod(idx, dim)
        Eab = sp.zeros(dim, dim); Eab[a, b] = 1
        out = f(Eab)
        cols.append(sp.Matrix(dim*dim, 1, [out[i // dim, i % dim] for i in range(dim*dim)]))
    return sp.Matrix.hstack(*cols)


def liouvillian(gd, gu, J=1):
    """Cell basis |a><b| ↦ 4a + b, site 0 the high bit; σ⁺ = |1><0|, cooling σ⁻ at γ↓, heating σ⁺ at γ↑."""
    sp_ = sp.Matrix([[0, 0], [1, 0]]); sm_ = sp.Matrix([[0, 1], [0, 0]]); I2 = sp.eye(2)
    S = [(kron(sp_, I2), kron(sm_, I2)), (kron(I2, sp_), kron(I2, sm_))]
    H = J*(S[0][0]*S[1][1] + S[0][1]*S[1][0])

    def D(rho):
        out = sp.zeros(4, 4)
        for l in (0, 1):
            for (Jop, rate) in ((S[l][1], gd[l]), (S[l][0], gu[l])):
                JJ = Jop.H*Jop
                out += rate*(Jop*rho*Jop.H - sp.Rational(1, 2)*(JJ*rho + rho*JJ))
        return out
    return superop_of(lambda r: -sp.I*(H*r - r*H) + D(r))


def transpose_superop():
    T = sp.zeros(16, 16)
    for a in range(4):
        for b in range(4):
            T[4*b + a, 4*a + b] = 1
    return T


TT = transpose_superop()


def adjoint_for(L, frame):
    """The right-hand generator: L† for the U frame, T·L†·T = Lᵀ_comp for the Γ frame."""
    return L.H if frame == 'U' else TT*L.H*TT


# realigned blocks by the first site's order shift δ
SITE_CELLS = [(0, 0), (0, 1), (1, 0), (1, 1)]          # (ket bit, bra bit)
MAPS = [(o, i) for o in SITE_CELLS for i in SITE_CELLS]


def dloc(c):
    return c[0] - c[1]


def realigned_blocks(X):
    blocks = {}
    for delta in (-2, -1, 0, 1, 2):
        rows = [m for m in MAPS if dloc(m[0]) - dloc(m[1]) == delta]
        cols = [m for m in MAPS if dloc(m[0]) - dloc(m[1]) == -delta]
        Rm = sp.zeros(len(rows), len(cols))
        for ri, (o0, i0) in enumerate(rows):
            for ci, (o1, i1) in enumerate(cols):
                out = 4*(2*o0[0] + o1[0]) + (2*o0[1] + o1[1])
                inn = 4*(2*i0[0] + i1[0]) + (2*i0[1] + i1[1])
                Rm[ri, ci] = X[out, inn]
        blocks[delta] = Rm
    return blocks


# ----- B: the point controls ----------------------------------------------------------------------------------------

def product_weight(gd, gu):
    W = kron(sp.diag(sp.sqrt(gd[0]), sp.sqrt(gu[0])), sp.diag(sp.sqrt(gd[1]), sp.sqrt(gu[1])))
    return superop_of(lambda r: W*r*W)


def part_b():
    print("B. the point controls (exact)")
    G = product_weight((3, 3), (1, 1))
    L = liouvillian((3, 3), (1, 1))
    zero = sp.zeros(16, 16)
    check("one temperature: Γ_β solves L·Γ = Γ·T L† T (the Γ frame)", (L*G - G*adjoint_for(L, 'Gamma')).applyfunc(sp.simplify) == zero)
    check("one temperature: U_β = Γ_β·T solves L·U = U·L† (the U frame)", (L*G*TT - G*TT*adjoint_for(L, 'U')).applyfunc(sp.simplify) == zero)
    check("Γ_β invertible", sp.simplify(G.det()) != 0)
    ranks = tuple(realigned_blocks(G)[d].rank() for d in (-2, -1, 0, 1, 2))
    check(f"Γ_β block ranks {ranks} = (0, 0, 1, 0, 0)", ranks == (0, 0, 1, 0, 0))
    L2 = liouvillian((3, 2), (1, 1))
    check("gradient: the one-temperature product solves neither frame",
          (L2*G - G*adjoint_for(L2, 'Gamma')).applyfunc(sp.simplify) != zero and
          (L2*G*TT - G*TT*adjoint_for(L2, 'U')).applyfunc(sp.simplify) != zero)
    G2 = product_weight((3, 2), (1, 1))
    check("gradient: its own product weight π₀^½ ⊗ π₁^½ does not solve the Γ frame either",
          (L2*G2 - G2*adjoint_for(L2, 'Gamma')).applyfunc(sp.simplify) != zero)


# ----- C: the sector chain modulo p ----------------------------------------------------------------------------------

def modp_expr(expr, gens, p, r):
    """A polynomial over ℚ(i) in gens ↦ its image over GF(p), i ↦ r with r² = −1."""
    poly = sp.Poly(sp.expand(expr), *gens, domain='EX')
    terms = []
    for mon, c in poly.terms():
        a, b = sp.sympify(c).as_real_imag()
        assert a.is_Rational and b.is_Rational, (a, b)   # exact coefficients only; a float here would be a defect
        v = (int(a.p)*pow(int(a.q), -1, p) + int(b.p)*pow(int(b.q), -1, p)*r) % p
        if v:
            terms.append(v*sp.Mul(*[g**e for g, e in zip(gens, mon)]))
    return sp.Add(*terms) if terms else sp.Integer(0)


def red(expr, gens, p):
    return sp.Poly(sp.expand(expr), *gens, modulus=p).as_expr()


def sector_chain(gd, gu, sectors, frame, p=P):
    """Impose the intertwining equations on the sectors in `sectors` only, ask the rank-4 shape of the minors whose
    entries all lie in those sectors, with those sectors invertible. Returns ('EMPTY' | 'solutions', dim, n_cond)."""
    assert sp.isprime(p) and p % 4 == 1
    r = int(sp.sqrt_mod(-1, p))
    L = liouvillian(gd, gu)
    Ladj = adjoint_for(L, frame)
    pc = lambda n: bin(n).count('1')
    dtot = {i: pc(i // 4) - pc(i % 4) for i in range(16)}
    by_d = {d: [i for i in range(16) if dtot[i] == d] for d in (-2, -1, 0, 1, 2)}
    X = sp.zeros(16, 16); syms = []
    for d, ids in by_d.items():
        for o in ids:
            for i in ids:
                s = sp.Symbol(f'x_{o}_{i}'); syms.append(s); X[o, i] = s
    eqs = []
    for d in sectors:
        ids = by_d[d]
        Ld = L.extract(ids, ids); Ad = Ladj.extract(ids, ids); Xd = X.extract(ids, ids)
        eqs += [sp.expand(e) for e in (Ld*Xd - Xd*Ad)]
    (vec,) = sp.linsolve(eqs, syms)
    free = sorted(set().union(*[v.free_symbols for v in vec]), key=str)
    invisible = set()
    for d, ids in by_d.items():
        if d not in sectors:
            invisible |= {X[o, i] for o in ids for i in ids}
    visible_free = [f for f in free if f not in invisible]
    Xs = X.subs(dict(zip(syms, vec)))
    Xp = Xs.applyfunc(lambda e: (e if e in invisible else modp_expr(e, free, p, r)) if e != 0 else sp.Integer(0))
    blocks = realigned_blocks(Xp)

    def minors(Mx, k):
        out = []
        for a in combinations(range(Mx.rows), k):
            for b in combinations(range(Mx.cols), k):
                sub = Mx.extract(list(a), list(b))
                if sub.free_symbols & invisible:
                    continue
                out.append(red(sub.det(method='berkowitz'), visible_free, p))
        return out
    conds = [b for b in (blocks[2][0, 0], blocks[-2][0, 0]) if not (sp.sympify(b).free_symbols & invisible)]
    conds += minors(blocks[1], 2) + minors(blocks[-1], 2) + minors(blocks[0], 3)
    conds = [c for c in conds if c != 0]
    zs = []
    for d in sectors:
        ids = by_d[d]
        det = red(Xp.extract(ids, ids).det(method='berkowitz'), visible_free, p)
        z = sp.Symbol(f'z{d}'); zs.append(z); conds.append(z*det - 1)
    gb = sp.groebner(conds, *visible_free, *zs, order='grevlex', modulus=p)
    return ('EMPTY' if list(gb.exprs) == [1] else 'solutions'), len(visible_free), len(conds)


def part_c(quick):
    print(f"C. the sector chain modulo p = {P} (both frames; a scout, the exact decision is the witness's)")
    rows = [((0,), 'solutions', 8)]
    if not quick:
        rows += [((0, 1, 2), 'solutions', 13), ((0, 1, -1), 'EMPTY', 16), ((-2, -1, 0, 1, 2), 'EMPTY', 18)]
    for frame in ('Gamma', 'U'):
        for sectors, expected, dim in rows:
            t = time.time()
            verdict, nvis, ncond = sector_chain((3, 2), (1, 1), list(sectors), frame)
            check(f"{frame} frame, sectors {list(sectors)}: {verdict} (expected {expected}), {nvis} visible parameters",
                  verdict == expected and nvis == dim, f"{ncond} conditions, {time.time()-t:.0f} s")
            sys.stdout.flush()


if __name__ == '__main__':
    quick = '--quick' in sys.argv
    t0 = time.time()
    part_a()
    part_b()
    part_c(quick)
    print(f"\n{'ALL PASS' if not FAILURES else 'FAILURES: ' + ', '.join(FAILURES)}  ({(time.time()-t0)/60:.1f} min)")
    sys.exit(1 if FAILURES else 0)
