"""Where the sharp subspaces of two coupled sodium nuclei grow, and where the kernel of the √3 coupling stops carrying them.

Reading for experiments/THE_SODIUM_CHAIN.md (the section on the full secular dipolar coupling). Two spin-3/2 nuclei with
identical rates (J0, J1, J2), one offset ω, an equal splitting q(ZZ₁ + ZZ₂) and the full secular dipolar coupling
d(3I_zI_z − I·I), the generator built by sodium_chain_mirrors.lindbladian. K_m is the largest L-invariant subspace
inside the cells with m bright nuclei, U_m the same inside the operators that the √3 part H₃ annihilates from both
sides. Generically dim K_1 = dim U_1 = 12 and dim K_2 = dim U_2 = 1 (the page's E₁, E₂).

What this script reads, over GF(p):
(1) The offset drops out: L preserves the coherence order and ω acts on order o as the scalar −iωo, which moves no
    invariant subspace. So the parameters are (J0, J1, J2, q, d).
(2) Per coherence order o and m = 1, 2, the set where dim K_m (dim U_m) exceeds its generic value is the common zero set
    of the maximal minors of the Krylov constraint matrix T = [Π̄·L^k on the m cells]. Along a random line in parameter
    space two random combinations of the maximal minors, det(R·T·C), are interpolated, their gcd taken and factored
    against the candidate components. A codimension-one component meets a generic line, so a gcd with no factor left
    over says that the list is the whole codimension-one part (for that line and that prime: evidence, not a proof).
(3) Random points on the intersections of the special hyperplanes, up to three at a time: dim K_m against dim U_m.

The plane q = d ≠ 0, where U_1 ≠ K_1, is gated exactly in SodiumChainWitness (gate 7): there the √3 part commutes with the
rest of H, and it and ½{Z₁ + Z₃, H₃} are eigenoperators of L in the one-bright cells, outside the kernel, a pair about
−s₁ that S does not keep.

Run: python simulations/sodium_chain_exceptional_set.py   (about seven minutes)
"""
from __future__ import annotations

import contextlib
import io
import random
import sys
from itertools import combinations

import numpy as np
import sympy as sp

sys.path.insert(0, __file__.rsplit('simulations', 1)[0] + 'simulations')
import sodium_chain_mirrors as scm  # noqa: E402

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

RESULTS = []


def check(name, ok):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)


# ---------------------------------------------------------------- the linear pieces of L

n, nq, dim = 2, 4, 16
with contextlib.redirect_stdout(io.StringIO()):
    RATIONAL, ROOT = scm.gate_dictionary()
ZERO = (0, 0, 0)
PIECES = {      # L is linear in every parameter; the pieces are taken at the verifier's J = 2 (jump rate 3), twice the
               # page's unit; the components below (J1 = 0, J2 = 0 and ratios of q to d) do not depend on that
    'J0': scm.lindbladian(n, {}, (2, 0, 0), q=(0,), w=(0,))[0],
    'J1': scm.lindbladian(n, {}, (0, 2, 0), q=(0,), w=(0,))[0],
    'J2': scm.lindbladian(n, {}, (0, 0, 2), q=(0,), w=(0,))[0],
    'q': scm.lindbladian(n, {}, ZERO, q=(1,), w=(0,))[0],          # H ∋ 4q·ZZ_k
    'w': scm.lindbladian(n, {}, ZERO, q=(0,), w=(1,))[0],
    'drat': scm.lindbladian(n, RATIONAL, ZERO, q=(0,), w=(0,))[0],  # H ∋ d·4(3IzIz − I·I), rational part
    'd3': scm.lindbladian(n, {}, ZERO, q=(0,), w=(0,), root=ROOT, r3=1.0)[0],
}
H3 = sum(c * scm.string(scm.place(nq, 0, w_)) for w_, c in ROOT.items())
W = sum(scm.bright(n))
LEVEL_2M = np.array([(3 - 2 * (a >> 2)) + (3 - 2 * (a & 3)) for a in range(dim)])
ORDER = np.array([(LEVEL_2M[a // dim] - LEVEL_2M[a % dim]) // 2 for a in range(dim * dim)])


class Field:
    def __init__(self, p):
        assert sp.isprime(p) and p % 12 == 1
        self.p = p
        self.r3 = sp.sqrt_mod(3, p)
        self.P = {k: scm.to_mod(v, p) for k, v in PIECES.items()}
        h = scm.to_mod(H3, p)
        self.left = np.kron(h, np.eye(dim, dtype=np.int64)) % p
        self.right = np.kron(np.eye(dim, dtype=np.int64), h.T) % p

    def L(self, J0, J1, J2, q, d, w=0):
        p, P = self.p, self.P
        return (J0 % p * P['J0'] + J1 % p * P['J1'] + J2 % p * P['J2'] + q % p * P['q'] + w % p * P['w']
                + d % p * P['drat'] + (d * self.r3 % p) * P['d3']) % p

    def hull_dim(self, Lp, m, withU, o=None):
        p = self.p
        idx = np.arange(dim * dim) if o is None else np.where(ORDER == o)[0]
        Lb = Lp[np.ix_(idx, idx)]
        inside = W[idx] == m
        I = np.eye(len(idx), dtype=np.int64)
        rows = [I[~inside]]
        if withU:
            rows += [self.left[:, idx], self.right[:, idx]]
        C, _ = scm.rref_mod(np.vstack(rows) % p, p)
        while True:
            Cn, _ = scm.rref_mod(np.vstack([C, C @ Lb % p]) % p, p)
            if Cn.shape[0] == C.shape[0]:
                return len(idx) - C.shape[0]
            C = Cn


# ---------------------------------------------------------------- (1) the offset drops out

F1 = Field(10009)
rnd = random.Random(2026)
pts = [[rnd.randrange(1, F1.p) for _ in range(5)] for _ in range(2)]
same = all(F1.hull_dim(F1.L(*pt, w=0), m, U) == F1.hull_dim(F1.L(*pt, w=rnd.randrange(1, F1.p)), m, U)
           for pt in pts for m in (1, 2) for U in (False, True))
check('(1) the offset does not move dim K_m or dim U_m (two random points, m = 1, 2)', same)
check('(1) generic dimensions K_1 = U_1 = 12, K_2 = U_2 = 1 at a random point',
      [F1.hull_dim(F1.L(*pts[0]), m, U) for m in (1, 2) for U in (False, True)] == [12, 12, 1, 1])


# ---------------------------------------------------------------- (2) the codimension-one components, per block

P2 = next(x for x in range(4000000, 4100000) if x % 12 == 1 and sp.isprime(x))   # small enough for int64 products
F2 = Field(P2)
t = sp.symbols('t')


def poly_mod(xs, ys, p):
    coef = list(ys)
    for j in range(1, len(xs)):
        for i in range(len(xs) - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) * pow((xs[i] - xs[i - j]) % p, p - 2, p) % p
    poly = sp.Poly(coef[-1], t, modulus=p)
    for i in range(len(xs) - 2, -1, -1):
        poly = poly * sp.Poly(t - xs[i], t, modulus=p) + sp.Poly(coef[i], t, modulus=p)
    return poly


def krylov_rows(F, pt, o, m, withU):
    p = F.p
    Lp = F.L(*pt)
    blk = np.where(ORDER == o)[0]
    Lb = Lp[np.ix_(blk, blk)]
    inside = W[blk] == m
    E, O = np.where(inside)[0], np.where(~inside)[0]
    X = np.eye(len(blk), dtype=np.int64)[:, E]
    rows = []
    if withU:
        rows += [F.left[:, blk][:, E], F.right[:, blk][:, E]]
    for _ in range(len(blk)):
        X = Lb @ X % p
        rows.append(X[O])
        if withU:
            full = np.zeros((dim * dim, X.shape[1]), np.int64)
            full[blk] = X
            rows += [F.left @ full % p, F.right @ full % p]
    return np.vstack(rows) % p, len(E), len(blk)


def exceptional_poly(F, o, m, withU, A, B, rnd):
    p = F.p
    T0, nE, nblk = krylov_rows(F, [(a + 5 * b) % p for a, b in zip(A, B)], o, m, withU)
    r = scm.rref_mod(T0, p)[0].shape[0]
    if r == 0:
        return None, nE - r
    polys = []
    for _ in range(2):
        R = np.array([[rnd.randrange(p) for _ in range(T0.shape[0])] for _ in range(r)], np.int64)
        C = np.array([[rnd.randrange(p) for _ in range(r)] for _ in range(nE)], np.int64)
        xs = list(range(1, r * nblk + 3))
        ys = [scm.det_mod(R @ (krylov_rows(F, [(a + x * b) % p for a, b in zip(A, B)], o, m, withU)[0] @ C % p) % p, p)
              for x in xs]
        polys.append(poly_mod(xs, ys, p))
    return sp.gcd(polys[0], polys[1]), nE - r


def factor_against(g, A, B, p):
    J0, J1, J2, q, d = [sp.Poly(a + b * t, t, modulus=p) for a, b in zip(A, B)]
    known = {'d': d, 'J1': J1, 'J2': J2, 'q - d': q - d, '2q - d': 2 * q - d, '2q - 3d': 2 * q - 3 * d,
             '28(q - d)² − 25d²': 28 * (q - d) ** 2 - 25 * d ** 2}
    found, rest = {}, g
    for name, f in known.items():
        while True:
            quo, rem = sp.div(rest, f)
            if not rem.is_zero:
                break
            rest = quo
            found[name] = found.get(name, 0) + 1
    return found, rest.degree()


rnd = random.Random(3)
A = [rnd.randrange(P2) for _ in range(5)]
B = [rnd.randrange(P2) for _ in range(5)]
components = {}
for withU in (False, True):
    for m in (1, 2):
        for o in range(0, 6):
            if not np.any((ORDER == o) & (W == m)):
                continue
            g, generic = exceptional_poly(F2, o, m, withU, A, B, rnd)
            label = f"{'U' if withU else 'K'}_{m} order {o}"
            if g is None:
                print(f'      {label}: the whole block is invariant generically ({generic})')
                components[(withU, m, o)] = {}
                continue
            found, left_over = factor_against(g, A, B, P2)
            found.pop('d', None)     # at d = 0 nothing moves brightness: K_m is the whole block, and V is everything
            components[(withU, m, o)] = found
            print(f'      {label}: generic {generic}, components {found or "none"}, unexplained degree {left_over}',
                  flush=True)
            check(f'(2) {label}: the gcd factors completely into the listed components', left_over == 0)

allK = {c for (U, m, o), f in components.items() if not U for c in f}
check('(2) the codimension-one components of K_1 and K_2 are J1 = 0, J2 = 0, q = d, 2q = d, 2q = 3d and '
      '28(q − d)² = 25d² (besides d = 0)',
      allK == {'J1', 'J2', 'q - d', '2q - d', '2q - 3d', '28(q - d)² − 25d²'})
diff = {key: (components[(False,) + key[1:]], f) for key, f in components.items()
        if key[0] and components[(False,) + key[1:]] != f}
check('(2) U has the same components with the same multiplicities in every block, except q = d, which only K_1 '
      'has (order 0)',
      list(diff) == [(True, 1, 0)] and diff[(True, 1, 0)][0] == {'q - d': 2} and diff[(True, 1, 0)][1] == {})


# ---------------------------------------------------------------- (3) the intersections

P3 = next(x for x in range(10009, 20000) if x % 12 == 1 and sp.isprime(x) and pow(7, (x - 1) // 2, x) == 1)
F3 = Field(P3)
s7 = sp.sqrt_mod(7, P3)
c_plus = (1 + 5 * pow(2 * s7, P3 - 2, P3)) % P3
c_minus = (1 - 5 * pow(2 * s7, P3 - 2, P3)) % P3
PLANES = {'J0': [1, 0, 0, 0, 0], 'J1': [0, 1, 0, 0, 0], 'J2': [0, 0, 1, 0, 0], 'q': [0, 0, 0, 1, 0],
          'q=d': [0, 0, 0, 1, -1], '2q=d': [0, 0, 0, 2, -1], '2q=3d': [0, 0, 0, 2, -3],
          'q=c+d': [0, 0, 0, 1, -c_plus], 'q=c-d': [0, 0, 0, 1, -c_minus],
          'J1=J2': [0, 1, -1, 0, 0], 'J0=J1': [1, -1, 0, 0, 0], 'J0=J2': [1, 0, -1, 0, 0], 'J1=-J2': [0, 1, 1, 0, 0]}
rnd = random.Random(11)
read, bad, deficits, on_plane, table = 0, [], {}, [], {}
for depth in (1, 2, 3):
    for combo in combinations(PLANES, depth):
        Bm = scm.null_mod(np.array([PLANES[c] for c in combo], np.int64) % P3, P3, 5)
        if Bm.shape[1] == 0:
            continue
        pt = [int(x) for x in Bm @ np.array([rnd.randrange(1, P3) for _ in range(Bm.shape[1])], np.int64) % P3]
        if pt[4] == 0:
            continue
        Lp = F3.L(*pt)
        dims = {(m, U): F3.hull_dim(Lp, m, U) for m in (1, 2) for U in (False, True)}
        table[combo] = (dims[(1, False)], dims[(1, True)], dims[(2, False)], dims[(2, True)])
        read += 1
        if 'q=d' in combo:
            on_plane.append(combo)
        if any(dims[(m, False)] != dims[(m, True)] for m in (1, 2)):
            bad.append(combo)
            deficits[combo] = tuple(dims[(m, False)] - dims[(m, True)] for m in (1, 2))
print('      (dim K_1, dim U_1, dim K_2, dim U_2) on each plane:')
for c in PLANES:
    if (c,) in table:
        print(f'        {c:7s} {table[(c,)]}')
for c in [('J1', 'q=d'), ('J1', 'J2', 'q=d')]:
    print(f'        {" & ".join(c):15s} {table[c]}')
EXPECTED = {('J0',): (12, 12, 1, 1), ('q',): (12, 12, 1, 1), ('J1',): (44, 44, 3, 3), ('J2',): (12, 12, 2, 2),
            ('2q=d',): (28, 28, 1, 1), ('2q=3d',): (28, 28, 1, 1), ('q=c+d',): (16, 16, 1, 1), ('q=c-d',): (16, 16, 1, 1),
            ('q=d',): (14, 12, 1, 1), ('J1', 'q=d'): (46, 44, 5, 5), ('J1', 'J2', 'q=d'): (76, 72, 6, 6)}
check('(3) the dimensions the page quotes: K_1 44, 28, 28, 16, 16 and K_2 3, 2 on the components, 14/12 on q = d, '
      '46/44 with J1 = 0 and 76/72 with J1 = J2 = 0 inside it',
      all(table.get(k) == v for k, v in EXPECTED.items()))
print(f'      {read} strata read; U ≠ K on {len(bad)}; the deficits (m = 1, m = 2) that occur: '
      f'{sorted(set(deficits.values()))}')
check('(3) U_m ≠ K_m only on strata inside q = d', bad and all('q=d' in c for c in bad))
check('(3) every stratum inside q = d has U_1 ≠ K_1, and U_2 = K_2 on every stratum read',
      set(on_plane) == set(bad) and all(v[0] > 0 and v[1] == 0 for v in deficits.values()))

print()
failed = [n_ for n_, ok in RESULTS if not ok]
print(f'{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed' + (f'; FAILED: {failed}' if failed else ''))
sys.exit(1 if failed else 0)
