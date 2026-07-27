"""Gate for F145 and F146: the seed triplet and the scalar count.

Self-contained: every object is rebuilt here from the mode index up, and nothing
is imported from the other gates or from the framework, so agreement with
`eta_ceiling_reduction.py` is evidence rather than a shared bug.  The energy
classification is additionally cross-checked against a float computation, since
it is the one place where an exact-arithmetic slip would be invisible.

Run:
    python simulations/scalar_count.py             # about a minute
    python simulations/scalar_count.py --deep      # + rungs 9 and 10, + M = 44

The statement being gated, in one breath.  On the open XY chain of N sites with
M := N + 1, the states that attain the disagreement floor of F144 on a diagonal
rung are the SO(3)-invariant multilinear couplings of one spin-1 per chiral pair,
so there are C(floor(N/2), l)*R_l of them with R_l the Riordan number, and they
are spanned by products of two connected blocks.

What it checks, block by block:

  W1  F145, the seed triplet.  On each chiral pair {a, abar} the three cell
      patterns u+ (pair in the bra), u0 (the F143 seed, one mode on both sides
      with a minus between the two choices) and u- (pair in the ket) form a
      spin-1 multiplet of the F142 spin ladder: S- lowers u+ to u0 to u- and
      annihilates u-, S+ mirrors it, and the two ladder coefficients multiply
      to 2, which is the spin-1 invariant.
  W2  the two connected blocks ARE the two invariants of SO(3) in that basis:
      the 2-block is the metric u0u0 - 2(u+u- + u-u+) and the 3-block is minus
      the volume form, both exactly, cell by cell, over the integers.
  W3  every block product is a maximizer: all three Corollary 7.3 conditions
      vanish on it exactly, over the integers.
  W4  F146, the count, squeezed from both sides so that it is exact and not
      extrapolated: rank_p(products) <= dim <= nullity_p(conditions), and the
      two meet at C(floor(N/2), l)*R_l.
  W5  the relations among the products are the CLASSICAL syzygies: evaluating
      the same partitions as products of dot and triple products of random
      vectors in R^3 gives the same rank and the same relation SUBSPACE.
  W6  step (A) of the proof: at a chiral-only rung the triplet restriction
      removes nothing, and at a resonant N it removes exactly the surplus while
      the contract value stays.
  W7  the three proved families (M prime, M = 2p, M = 2^a) carry no non-chiral
      coincidence at any rung searched, and the measured minimal rung j(M)
      elsewhere, including the forward prediction j(4p) = (p+1)/2 at M = 44
      under --deep.

Every rank over GF(p) is one-sided in a direction that is stated where it is
used: a rank mod p can only drop, so it bounds a dimension from below, and a
nullity mod p can only grow, so it bounds one from above.
"""
import argparse
import itertools
import math
import sys
import time
from fractions import Fraction

import numpy as np

FAILURES = []
CHECKS = 0

RIORDAN = {2: 1, 3: 1, 4: 3, 5: 6, 6: 15, 7: 36, 8: 91, 9: 232, 10: 603}


def check(label, ok, detail=""):
    global CHECKS
    CHECKS += 1
    if not ok:
        FAILURES.append(label)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok


# --------------------------------------------------------------------------
# modes, exact Slater energies, cells
# --------------------------------------------------------------------------

def cyclotomic(n):
    """Phi_n as a coefficient list, highest degree first, by exact division."""
    poly = [1] + [0] * (n - 1) + [-1]                 # x^n - 1
    for d in range(1, n):
        if n % d:
            continue
        div = cyclotomic(d)
        q, r = [0] * (len(poly) - len(div) + 1), list(poly)
        for i in range(len(q)):
            c = r[i] // div[0]
            q[i] = c
            if c:
                for j in range(len(div)):
                    r[i + j] -= c * div[j]
        poly = q
    return poly


def poly_rem(num, phi):
    num = list(num)
    dp = len(phi) - 1
    for i in range(len(num) - dp):
        c = num[i]
        if c:
            for j in range(len(phi)):
                num[i + j] -= c * phi[j]
    return tuple(num[len(num) - dp:])


def energy_classes(N, l):
    """l-subsets of the modes grouped by EXACT Slater energy.

    The energy of a subset is the sum of eps_k = 2J cos(pi k / M); written in the
    2M-th roots of unity it is an algebraic integer, and two subsets agree iff
    their coefficient vectors agree modulo the cyclotomic polynomial.  No
    tolerance is involved.
    """
    twoM = 2 * (N + 1)
    phi = cyclotomic(twoM)
    buckets = {}
    for S in itertools.combinations(range(1, N + 1), l):
        v = [0] * twoM
        for k in S:
            v[k] += 1
            v[twoM - k] += 1
        buckets.setdefault(poly_rem(v[::-1], phi), []).append(S)
    return list(buckets.values())


def sign_remove(S, m):
    return -1 if sum(1 for j in S if j < m) % 2 else 1


def w0_energy_classes_agree_with_floats(N, l):
    """Cross-check: the exact classes and the float energies must agree."""
    M = N + 1
    eps = {k: 2 * math.cos(math.pi * k / M) for k in range(1, N + 1)}
    groups = energy_classes(N, l)
    lab = {}
    for i, g in enumerate(groups):
        for S in g:
            lab[S] = i
    vals = {S: sum(eps[k] for k in S) for S in lab}
    items = sorted(vals, key=lambda S: vals[S])
    for A, B in zip(items, items[1:]):
        same_float = abs(vals[A] - vals[B]) < 1e-9
        if same_float != (lab[A] == lab[B]):
            return False
    return True


# --------------------------------------------------------------------------
# the triplet, and the two ladders
# --------------------------------------------------------------------------

def triplet(M, x):
    """The three cell patterns of the chiral pair of x, as partial cells."""
    a, ab = sorted((x, M - x))
    return {1: {((a, ab), ()): 1},
            0: {((a,), (a,)): 1, ((ab,), (ab,)): -1},
            -1: {((), (a, ab)): 1}}


def s_minus(N, v):
    """The F142 spin ladder lowering: move a mode from the bra to the ket at the
    CHIRAL partner index, which is the staggering in the mode basis."""
    M = N + 1
    out = {}
    for (A, B), c in v.items():
        sB = set(B)
        for b in range(1, N + 1):
            bb = M - b
            if bb not in A or b in sB:
                continue
            A2 = tuple(x for x in A if x != bb)
            B2 = tuple(sorted(B + (b,)))
            out[(A2, B2)] = out.get((A2, B2), 0) + \
                c * sign_remove(A, bb) * sign_remove(B2, b)
    return {k: x for k, x in out.items() if x}


def s_plus(N, v):
    M = N + 1
    out = {}
    for (A, B), c in v.items():
        sA = set(A)
        for b in B:
            bb = M - b
            if bb in sA:
                continue
            B2 = tuple(x for x in B if x != b)
            A2 = tuple(sorted(A + (bb,)))
            out[(A2, B2)] = out.get((A2, B2), 0) + \
                c * sign_remove(B, b) * sign_remove(A2, bb)
    return {k: x for k, x in out.items() if x}


def proportional(u, v):
    """(cu, cv) with cv*u == cu*v cell by cell, or None."""
    if not u or not v:
        return None
    r = None
    for k in set(u) | set(v):
        a, b = u.get(k, 0), v.get(k, 0)
        if (a == 0) != (b == 0):
            return None
        if a and b:
            if r is None:
                r = (a, b)
            elif a * r[1] != b * r[0]:
                return None
    return r


# --------------------------------------------------------------------------
# the blocks and their products
# --------------------------------------------------------------------------

def merge(A1, A2):
    arr = list(A1) + list(A2)
    sgn = 1
    for i in range(len(arr)):
        for j in range(len(arr) - 1 - i):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                sgn = -sgn
    return sgn, tuple(arr)


def merge_cells(parts):
    out = {((), ()): 1}
    for blk in parts:
        new = {}
        for (A1, B1), c1 in out.items():
            for (A2, B2), c2 in blk.items():
                if set(A1) & set(A2) or set(B1) & set(B2):
                    continue
                sA, A = merge(A1, A2)
                sB, B = merge(B1, B2)
                new[(A, B)] = new.get((A, B), 0) + c1 * c2 * sA * sB
        out = {k: v for k, v in new.items() if v}
    return out


def block2(M, a, b):
    """The 2-block, six terms with entries in {1, 2}: the closed form of
    ETA_CEILING_REDUCTION, restated here so this gate stands alone."""
    pp, qq, rr, ss = sorted([a, M - a, b, M - b])
    return {((pp, qq), (pp, qq)): 1, ((pp, rr), (pp, rr)): -1,
            ((qq, ss), (qq, ss)): -1, ((rr, ss), (rr, ss)): 1,
            ((pp, ss), (qq, rr)): -2, ((qq, rr), (pp, ss)): -2}


def block3(M, a, b, c):
    """The 3-block in closed form (F146): twelve terms, entries +-1.  One chiral
    pair full in the bra, another full in the ket, one shared mode of the third."""
    P = [tuple(sorted((x, M - x))) for x in sorted((a, b, c))]
    out = {}
    for i, j in itertools.permutations(range(3), 2):
        k = ({0, 1, 2} - {i, j}).pop()
        for t, m in enumerate(P[k]):
            A = tuple(sorted(P[i] + (m,)))
            B = tuple(sorted(P[j] + (m,)))
            out[(A, B)] = (1 if i < j else -1) * (1 if t == 0 else -1)
    return out


def metric(M, x, y):
    """u0(P)u0(Q) - 2[u+(P)u-(Q) + u-(P)u+(Q)], the SO(3) metric in this basis."""
    tx, ty = triplet(M, x), triplet(M, y)
    out = {}
    for (mx, my, co) in [(0, 0, 1), (1, -1, -2), (-1, 1, -2)]:
        for k, v in merge_cells([tx[mx], ty[my]]).items():
            out[k] = out.get(k, 0) + co * v
    return {k: v for k, v in out.items() if v}


def volume(M, x, y, z):
    """The antisymmetric invariant: the three m-values are +1, 0, -1 in some
    order, with the permutation sign."""
    T = [triplet(M, p) for p in (x, y, z)]
    out = {}
    for perm in itertools.permutations(range(3)):
        pl, sgn = list(perm), 1
        for i in range(3):
            for j in range(2 - i):
                if pl[j] > pl[j + 1]:
                    pl[j], pl[j + 1] = pl[j + 1], pl[j]
                    sgn = -sgn
        ms = {perm[0]: 1, perm[1]: 0, perm[2]: -1}
        for k, v in merge_cells([T[i][ms[i]] for i in range(3)]).items():
            out[k] = out.get(k, 0) + sgn * v
    return {k: v for k, v in out.items() if v}


def partitions_23(items):
    """Set partitions of `items` into blocks of size 2 and 3."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for k in (1, 2):
        for others in itertools.combinations(rest, k):
            blk = (first,) + others
            left = [x for x in rest if x not in others]
            for p in partitions_23(left):
                yield [blk] + p


def block_products(N, l, all_choices=False):
    """Every product of 2- and 3-blocks on l chiral pairs.

    With `all_choices` the products are built over EVERY l-subset of the chiral
    pairs, which is what a lower bound on the whole maximizer space needs; the
    factor C(floor(N/2), l) is then measured rather than assumed."""
    M = N + 1
    chir = [a for a in range(1, M // 2 + 1) if a < M - a]
    choices = (list(itertools.combinations(chir, l)) if all_choices
               else [tuple(chir[:l])])
    parts, vecs = [], []
    for choice in choices:
        for part in partitions_23(list(choice)):
            blocks = [block2(M, *blk) if len(blk) == 2 else block3(M, *blk)
                      for blk in part]
            parts.append(part)
            vecs.append(merge_cells(blocks))
    return parts, vecs


def satisfies_conditions(N, v):
    """All three Corollary 7.3 conditions, exactly, over the integers."""
    M = N + 1
    for (A, B) in v:
        if {M - b for b in B} & set(A):
            return False, "X(a,abar)"
    acc = {}
    for (A, B), val in v.items():
        for x in set(A) & set(B):
            key = (min(x, M - x), tuple(y for y in A if y != x),
                   tuple(y for y in B if y != x))
            acc[key] = acc.get(key, 0) + val * sign_remove(A, x) * sign_remove(B, x)
    if any(acc.values()):
        return False, "X(a,a)+X(abar,abar)"
    if any(s_minus(N, v).values()):
        return False, "S-"
    return True, "all"


# --------------------------------------------------------------------------
# the maximizer space: delete, pair, then rank the spin ladder
# --------------------------------------------------------------------------

class Signed:
    """Union-find with a sign and a zero class."""

    def __init__(self):
        self.par, self.sgn, self.zero = {}, {}, set()

    def find(self, x):
        if x not in self.par:
            self.par[x], self.sgn[x] = x, 1
            return x, 1
        s, root = 1, x
        while self.par[root] != root:
            s *= self.sgn[root]
            root = self.par[root]
        cur, cs = x, s
        while self.par[cur] != cur:
            nxt, ns = self.par[cur], self.sgn[cur]
            self.par[cur], self.sgn[cur] = root, cs
            cur, cs = nxt, cs * ns
        return root, s

    def union(self, x, y, s):
        rx, sx = self.find(x)
        ry, sy = self.find(y)
        if rx == ry:
            if sx != s * sy:
                self.kill(rx)
            return
        self.par[rx], self.sgn[rx] = ry, s * sy * sx
        if rx in self.zero:
            self.zero.discard(rx)
            self.zero.add(ry)

    def kill(self, x):
        self.zero.add(self.find(x)[0])

    def is_zero(self, x):
        return self.find(x)[0] in self.zero


def is_triplet_cell(N, A, B):
    M = N + 1
    sA, sB = set(A), set(B)
    seen = set()
    for x in list(sA | sB):
        p = tuple(sorted((x, M - x)))
        if p[0] == p[1]:
            return False
        if p in seen:
            continue
        seen.add(p)
        a, ab = p
        if (a in sA, ab in sA, a in sB, ab in sB) in [
                (True, True, False, False), (False, False, True, True),
                (True, False, True, False), (False, True, False, True)]:
            continue
        return False
    return True


def reduce_space(N, l, triplet_only=False):
    """Cells of V0 after the two chiral defect conditions, and the S- rows on
    what is left.  The deletion removes coordinates and the pairing glues them,
    so no elimination happens before the last step."""
    M = N + 1
    pairs = []
    for grp in energy_classes(N, l):
        comps = [{M - b for b in S} for S in grp]
        sets = [set(S) for S in grp]
        for A in grp:
            for iB, B in enumerate(grp):
                if comps[iB] & set(A):
                    continue
                if triplet_only and not is_triplet_cell(N, A, B):
                    continue
                pairs.append((A, B))
    ipair = {p: i for i, p in enumerate(pairs)}
    uf = Signed()
    for i, (A, B) in enumerate(pairs):
        for a in set(A) & set(B):
            ab = M - a
            A2 = tuple(x for x in A if x != a)
            B2 = tuple(x for x in B if x != a)
            s1 = sign_remove(A, a) * sign_remove(B, a)
            if ab in A2 or ab in B2:
                uf.kill(i)
                continue
            A3, B3 = tuple(sorted(A2 + (ab,))), tuple(sorted(B2 + (ab,)))
            j = ipair.get((A3, B3))
            if j is None:
                uf.kill(i)
                continue
            uf.union(i, j, -s1 * sign_remove(A3, ab) * sign_remove(B3, ab))
    comp, roots = [None] * len(pairs), {}
    for i in range(len(pairs)):
        if uf.is_zero(i):
            continue
        r, s = uf.find(i)
        if r not in roots:
            roots[r] = len(roots)
        comp[i] = (roots[r], s)
    rows = {}
    for i, (A, B) in enumerate(pairs):
        if comp[i] is None:
            continue
        k, sgn = comp[i]
        sB = set(B)
        for b in range(1, N + 1):
            bb = M - b
            if bb not in A or b in sB:
                continue
            A2 = tuple(x for x in A if x != bb)
            B2 = tuple(sorted(B + (b,)))
            row = rows.setdefault((A2, B2), {})
            row[k] = row.get(k, 0) + sign_remove(A, bb) * sign_remove(B2, b) * sgn
    return pairs, comp, len(roots), [r for r in rows.values() if any(r.values())]


def psi_rows(N, l, pairs, comp, ncomp):
    """The lowest-weight rows Psi = sum_a X(a,a) on the surviving variables.

    Condition (ii) implies Psi v = 0, so the lowest-weight condition needs no
    separate imposition; that is an argument, and this builds the rows so the
    claim can be MEASURED rather than asserted."""
    rows = {}
    for i, (A, B) in enumerate(pairs):
        if comp[i] is None:
            continue
        k, sgn = comp[i]
        for a in set(A) & set(B):
            key = (tuple(x for x in A if x != a), tuple(x for x in B if x != a))
            row = rows.setdefault(key, {})
            row[k] = row.get(k, 0) + sign_remove(A, a) * sign_remove(B, a) * sgn
    return [r for r in rows.values() if any(r.values())]


PS = 999983


def _prev_prime(n):
    while True:
        n -= 1
        if all(n % d for d in range(2, int(n ** 0.5) + 1)):
            return n


def rank_mod_p(rows, ncols, p=PS, target=None):
    """Rank over GF(p), kept in reduced echelon form as float64 so the dot
    products go through BLAS.  Entries stay below p and n*p^2 < 2^53, so the
    arithmetic is exact."""
    while ncols * (p - 1) ** 2 >= 2 ** 53:
        p = _prev_prime(p // 2)
    P = np.zeros((ncols, ncols), dtype=np.float64)
    piv = np.zeros(ncols, dtype=np.intp)
    n = 0
    for r in rows:
        v = np.zeros(ncols, dtype=np.float64)
        for c, val in r.items():
            v[c] = val % p
        if n:
            coef = v[piv[:n]]
            if coef.any():
                v -= coef @ P[:n]
                v %= p
        nz = np.nonzero(v)[0]
        if not nz.size:
            continue
        c = int(nz[0])
        v = (v * pow(int(v[c]), p - 2, p)) % p
        col = P[:n, c]
        hit = np.nonzero(col)[0]
        if hit.size:
            P[hit] = (P[hit] - col[hit, None] * v) % p
        P[n] = v
        piv[n] = c
        n += 1
        if target is not None and n >= target:
            return n            # reaching the target already bounds the nullity
    return n


def rank_mod_p_exact(rows, p=2147483647):
    """The same rank by sparse exact elimination at a prime near 2^31, kept as an
    independent path: it shares no code with the float64 one."""
    pivots = {}
    rank = 0
    for r in rows:
        row = {c: v % p for c, v in r.items() if v % p}
        while row:
            c = min(row)
            if c in pivots:
                f = row[c]
                for cc, vv in pivots[c].items():
                    nv = (row.get(cc, 0) - f * vv) % p
                    if nv:
                        row[cc] = nv
                    else:
                        row.pop(cc, None)
            else:
                inv = pow(row[c], p - 2, p)
                pivots[c] = {cc: (vv * inv) % p for cc, vv in row.items()}
                rank += 1
                break
    return rank


def maximizer_dimension(N, l, triplet_only=False, exact_too=False, target_dim=None,
                        with_psi=False):
    """`target_dim` stops the rank as soon as the nullity has come DOWN to it.
    Since a nullity mod p can only be too large, reaching the target already proves
    dim <= target_dim, which with the product rank from below is the whole squeeze;
    it saves most of the rows at the largest rung."""
    pairs, comp, ncomp, rows = reduce_space(N, l, triplet_only)
    if with_psi:
        rows = rows + psi_rows(N, l, pairs, comp, ncomp)
    tgt = None if target_dim is None else ncomp - target_dim
    d = ncomp - rank_mod_p(rows, ncomp, target=tgt)
    if exact_too:
        d2 = ncomp - rank_mod_p_exact(rows)
        if d != d2:
            return None, len(pairs), ncomp
    return d, len(pairs), ncomp


# --------------------------------------------------------------------------
# the classical side, and the coincidences
# --------------------------------------------------------------------------

def classical_matrix(l, nsamp=4000, seed=11):
    rng = np.random.default_rng(seed)
    parts = list(partitions_23(list(range(l))))
    V = rng.normal(size=(nsamp, l, 3))
    rows = []
    for part in parts:
        val = np.ones(nsamp)
        for blk in part:
            if len(blk) == 2:
                val *= np.einsum('ni,ni->n', V[:, blk[0]], V[:, blk[1]])
            else:
                i, j, k = sorted(blk)
                val *= np.linalg.det(V[:, [i, j, k]])
        rows.append(val)
    return np.array(rows), parts


def left_null(A, tol=1e-10):
    """Kernel of c -> sum_p c_p A_p.  The row scaling is undone afterwards, so
    the kernel returned is the kernel of the ORIGINAL rows; forgetting that
    compares two different subspaces."""
    n = np.maximum(np.linalg.norm(A, axis=1), 1e-300)
    An = A / n[:, None]
    G = An @ An.T
    w, U = np.linalg.eigh((G + G.T) / 2)
    keep = w <= tol * max(w.max(), 1e-300)
    return U[:, keep] / n[:, None], int((~keep).sum())


def span_dim(X, tol=1e-8):
    X = X / np.maximum(np.linalg.norm(X, axis=0, keepdims=True), 1e-300)
    sv = np.linalg.svd(X, compute_uv=False)
    return int((sv > tol * sv[0]).sum())


def nonchiral_coincidence(N, l):
    """A pair of equal-energy l-subsets differing by something other than whole
    chiral pairs, with the deleted cells excluded, or None."""
    M = N + 1
    for grp in energy_classes(N, l):
        if len(grp) == 1:
            continue
        sets = [set(S) for S in grp]
        for i, X in enumerate(sets):
            Xc = {M - x for x in X}
            for k in range(i + 1, len(sets)):
                Y = sets[k]
                if X & Y or Xc & Y:
                    continue
                if all((M - x) in X for x in X) and all((M - y) in Y for y in Y):
                    continue
                return sorted(X), sorted(Y)
    return None


def minimal_resonant_rung(N, jmax):
    for j in range(2, jmax + 1):
        hit = nonchiral_coincidence(N, j)
        if hit:
            return j, hit
    return None, None


# --------------------------------------------------------------------------
# the blocks of the gate
# --------------------------------------------------------------------------

def w1_seed_triplet():
    print("\nW1  F145: the three cell patterns of a chiral pair are a spin-1 "
          "multiplet of the F142 ladder")
    for N in (6, 8, 9, 10, 12, 15):
        M = N + 1
        bad, tot = 0, 0
        for x in [a for a in range(1, M // 2 + 1) if a < M - a]:
            t = triplet(M, x)
            tot += 1
            ok = (proportional(s_minus(N, t[1]), t[0])
                  and proportional(s_minus(N, t[0]), t[-1])
                  and not s_minus(N, t[-1])
                  and proportional(s_plus(N, t[-1]), t[0])
                  and proportional(s_plus(N, t[0]), t[1])
                  and not s_plus(N, t[1]))
            bad += not ok
        # the two lowering coefficients multiply to 2, the spin-1 invariant
        t = triplet(M, 1)
        c1 = proportional(s_minus(N, t[1]), t[0])
        c2 = proportional(s_minus(N, t[0]), t[-1])
        prod = abs(c1[0] * c2[0]) if c1 and c2 else 0
        check(f"N={N}: all {tot} chiral pairs give a spin-1 multiplet, "
              f"S- coefficients multiply to 2", bad == 0 and prod == 2)


def w2_blocks_are_the_invariants(deep):
    print("\nW2  the 2-block IS the SO(3) metric and the 3-block IS the volume, "
          "cell by cell")
    for N in (6, 7, 8, 9, 10, 12, 13, 15, 16) if deep else (6, 7, 8, 9, 10, 12):
        M = N + 1
        chir = [a for a in range(1, M // 2 + 1) if a < M - a]
        def same(u, v, sign):
            r = proportional(u, v)
            return r is not None and r[0] == sign * r[1]
        bad2 = sum(not same(block2(M, x, y), metric(M, x, y), 1)
                   for x, y in itertools.combinations(chir, 2))
        bad3 = sum(not same(block3(M, *t), volume(M, *t), -1)
                   for t in itertools.combinations(chir, 3))
        n2 = math.comb(len(chir), 2)
        n3 = math.comb(len(chir), 3)
        check(f"N={N}: {n2} 2-blocks equal the metric, {n3} 3-blocks equal minus "
              f"the volume", bad2 == 0 and bad3 == 0)


def w3_products_are_maximizers(deep):
    print("\nW3  every block product satisfies all three Corollary 7.3 conditions, "
          "over the integers")
    cases = [(8, 4), (10, 5), (12, 6), (15, 7)] + ([(16, 8)] if deep else [])
    for (N, l) in cases:
        parts, vecs = block_products(N, l)
        bad = [p for p, v in zip(parts, vecs) if not satisfies_conditions(N, v)[0]]
        check(f"N={N} l={l}: {len(parts)} products, all maximizers", not bad)


def w4_the_count(deep):
    print("\nW4  F146: the count, squeezed between the product rank from below and "
          "the condition nullity from above")
    cases = [(8, 3), (10, 4), (10, 5), (12, 6), (15, 7), (16, 8)]
    if deep:
        cases += [(18, 9), (21, 10)]
    for (N, l) in cases:
        p = N // 2
        contract = math.comb(p, l) * RIORDAN[l]
        upper, ncells, ncomp = maximizer_dimension(
            N, l, exact_too=(l <= 6), target_dim=(contract if l >= 9 else None))
        parts, vecs = block_products(N, l, all_choices=True)
        keys = sorted({k for v in vecs for k in v})
        ik = {k: i for i, k in enumerate(keys)}
        lower = rank_mod_p([{ik[k]: c for k, c in v.items()} for v in vecs], len(keys))
        if l <= 7:
            with_psi, _, _ = maximizer_dimension(N, l, with_psi=True)
            check(f"N={N} l={l}: imposing the lowest-weight rows Psi as well leaves "
                  f"the nullity at {with_psi}, so (ii) really implies Psi v = 0",
                  with_psi == upper)
        check(f"N={N} l={l}: {ncells} cells -> {ncomp} variables, "
              f"{len(parts)} products give {lower} <= dim <= {upper}, "
              f"and C(p,l)*R_l = {contract}",
              lower == upper == contract,
              f"[R_{l} = {RIORDAN[l]}, matchings would give "
              f"{math.prod(range(1, l, 2)) if l % 2 == 0 else 0}]" if l % 2 == 0 else "")


def w5_relations_are_classical(deep):
    print("\nW5  the relations among the products are the classical SO(3) syzygies")
    C2 = -4.0                      # fitted at l = 6, prediction at 7 and 8
    cases = [(10, 5), (12, 6), (15, 7)] + ([(16, 8)] if deep else [])
    for (N, l) in cases:
        parts, vecs = block_products(N, l)
        keys = sorted({k for v in vecs for k in v})
        ik = {k: i for i, k in enumerate(keys)}
        A = np.zeros((len(vecs), len(keys)))
        for r, v in enumerate(vecs):
            for k, c in v.items():
                A[r, ik[k]] = c
        B, parts2 = classical_matrix(l)
        t = np.array([sum(1 for b in p if len(b) == 3) for p in parts])
        Ko, ro = left_null(A)
        Kc, rc = left_null(B)
        Kcs = Kc * (C2 ** ((t - t.min()) // 2))[:, None]
        joint = span_dim(np.hstack([Ko, Kcs])) if Ko.shape[1] else 0
        check(f"N={N} l={l}: rank {ro} = {rc}, relations {Ko.shape[1]} = "
              f"{Kc.shape[1]}, joint span {joint}: the SAME relation space",
              ro == rc == RIORDAN[l] and joint == Ko.shape[1] == Kc.shape[1])


def w6_triplet_sector(deep):
    print("\nW6  step (A): at a chiral-only rung the triplet restriction removes "
          "nothing; at a resonant N it removes exactly the surplus")
    for (N, l) in [(8, 3), (10, 4), (12, 6), (15, 7)] + ([(16, 8)] if deep else []):
        d1, c1, _ = maximizer_dimension(N, l)
        d2, c2, _ = maximizer_dimension(N, l, triplet_only=True)
        check(f"N={N} l={l} (chiral-only rung): {c1} cells, restriction removes "
              f"{c1 - c2}, dim {d1} = {d2}", c1 == c2 and d1 == d2)
    for (N, l) in [(11, 2), (11, 4), (14, 4), (19, 3)] + ([(20, 3)] if deep else []):
        d1, c1, _ = maximizer_dimension(N, l)
        d2, c2, _ = maximizer_dimension(N, l, triplet_only=True)
        contract = math.comb(N // 2, l) * RIORDAN[l]
        check(f"N={N} l={l} (resonant): dim {d1} > triplet part {d2} = "
              f"C(p,l)*R_l = {contract}, surplus {d1 - d2} entirely non-triplet",
              d2 == contract and d1 > d2)


def zero_mode_cells(N, l):
    """Surviving cells whose difference contains the self-paired zero mode M/2.

    The proof excludes it by a weight parity (a chiral pair contributes 4, the zero
    mode 2, the sum has weight 4d), so this must be 0 at every even M."""
    M = N + 1
    if M % 2:
        return 0, 0
    z = M // 2
    total = zero = 0
    for grp in energy_classes(N, l):
        for A in grp:
            sA = set(A)
            for B in grp:
                if {M - b for b in B} & sA:
                    continue
                total += 1
                if z in (sA - set(B)) | (set(B) - sA):
                    zero += 1
    return total, zero


def w7_the_three_families(deep):
    print("\nW7  the proved families carry no non-chiral coincidence, and the "
          "minimal resonant rung elsewhere")
    for (N, tag, jmax) in [(6, "M=7 prime", 3), (10, "M=11 prime", 5),
                           (12, "M=13 prime", 6), (16, "M=17 prime", 8),
                           (18, "M=19 prime", 9),
                           (9, "M=10 = 2*5", 4), (13, "M=14 = 2*7", 6),
                           (21, "M=22 = 2*11", 10),
                           (7, "M=8 = 2^3", 3), (15, "M=16 = 2^4", 7),
                           (31, "M=32 = 2^5", 4)]:
        j, hit = minimal_resonant_rung(N, jmax)
        check(f"N={N} ({tag}): clean to rung {jmax}", j is None,
              "" if j is None else f"resonant at {j}: {hit}")
    for (N, expect) in [(11, 2), (14, 2), (17, 2), (19, 3), (20, 2), (26, 3),
                        (27, 4), (34, 4), (24, 5)]:
        j, hit = minimal_resonant_rung(N, expect)
        check(f"M={N + 1}: minimal resonant rung j = {expect}", j == expect,
              f"witness {hit[0]} vs {hit[1]}" if hit else "")
    for (N, l) in [(7, 3), (9, 4), (13, 4), (15, 4), (21, 3)]:
        total, zero = zero_mode_cells(N, l)
        check(f"M={N + 1} l={l}: the zero mode M/2 appears in none of the {total} "
              f"surviving differences, as the weight parity requires", zero == 0)
    if deep:
        # the forward prediction: j(4p) = (p+1)/2, so M = 44 is clean at rung 5
        # and breaks first at rung 6
        clean5 = nonchiral_coincidence(43, 5) is None
        hit6 = nonchiral_coincidence(43, 6)
        check("M=44 (p=11): clean at rung 5 and first break at rung 6, "
              "as (p+1)/2 predicts", clean5 and hit6 is not None,
              f"witness {hit6[0]} vs {hit6[1]}" if hit6 else "")


def w0_energy_classes(deep):
    print("\nW0  the exact energy classification agrees with the float energies")
    for (N, l) in [(8, 3), (10, 4), (12, 4), (14, 3)] + ([(16, 4)] if deep else []):
        check(f"N={N} l={l}: exact classes = float classes",
              w0_energy_classes_agree_with_floats(N, l))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deep", action="store_true",
                    help="the larger rungs and the M = 44 prediction")
    args = ap.parse_args()
    t0 = time.time()
    print("scalar-count gate (F145, F146)")
    w0_energy_classes(args.deep)
    w1_seed_triplet()
    w2_blocks_are_the_invariants(args.deep)
    w3_products_are_maximizers(args.deep)
    w4_the_count(args.deep)
    w5_relations_are_classical(args.deep)
    w6_triplet_sector(args.deep)
    w7_the_three_families(args.deep)
    print(f"\n{CHECKS} checks in {time.time() - t0:.0f}s")
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES:")
        for f in FAILURES:
            print("   " + f)
        sys.exit(1)
    print("all green")


if __name__ == "__main__":
    main()
