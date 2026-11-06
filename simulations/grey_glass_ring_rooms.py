"""The grey glass on the even XY ring past N = 8: the room at Hamming distance N/2 in the block (N/2, 2).

Companion gate to docs/THE_GREY_GLASS.md (the paragraphs "Every division into pairs, and the rings to
N = 12"), kept apart from simulations/grey_glass_gate.py because at N = 12 it ranks matrices of thousands of
rows and runs long. Must print
"grey glass ring rooms: ALL GREEN".

Spin language, bit l = site l, bonds J_l between the sites l and l + 1, the wrap bond J_w between N - 1 and 0.
K = H_M (x) 1 - 1 (x) H_2 is the Hamiltonian part on the block (M, 2), M = N/2 (L_H = -i K); the room is the
largest subspace K keeps inside the cells at Hamming distance M.

  R1 the state with one pair full, for every division into pairs: for each of the (N - 1)!! divisions of the
     sites into pairs, no hop stays inside the configurations with one pair full and the rest empty, and H keeps
     a state on them exactly when the division is kept by the couplings between pairs (positive couplings) and
     its pairs, joined by those bonds, form a bipartite graph; the state is then the staggered sum; on the
     uniform ring these are the N/2 bond axes and, at even M, the antipodal map (exact, exhaustive at N = 6, 8,
     10, 12).
  R1b couplings of both signs: on random rings at N = 6 and 8, over every division, H keeps a sum of full pairs
     exactly when a sign per pair, D, gives P h P = -D h D on every bond between pairs (P the exchange inside every pair,
     h the spin hopping matrix, J_w on the wrap): a two-excitation state is a symmetric matrix psi, H acts as h psi + psi h off the
     diagonal, and the sum of full pairs is psi = D P.
  R2 the room from above by momentum: on the uniform ring the simultaneous rotation commutes with K and with the
     plane; per eigenvalue, the room over GF(p) as nested kernels V_{i+1} = {v in V_i : K v in V_i}. The
     rational room's integer lattice is a direct summand that K keeps inside the plane; reduced modulo p it
     stays kept, inside the plane and of the same dimension, so the room over GF(p) bounds the rational room
     from above.
  R3 the room from above by free fermions: at a prime where the single-particle matrix (wrap entry -J_w, both
     sides holding an even number of particles) splits, the minors det U[A, S] are the many-body eigenvectors
     (checked against the spin Hamiltonian modulo p), K is diagonal on their products, and the room is the sum
     over the frequencies of dim E_omega minus the rank of the off-plane rows; a subset of rows only lowers that
     rank, so this too is an upper bound. On the uniform ring at N = 12 it meets R2: two routes, one number.
  R4 the room from below: the ladder images of (M, 0) (shown kept by the exact integer intertwining K lad = lad H_M), the
     bond-axis products and the antipodal products, every row exact, at distance M alone and kept; their rank
     modulo p bounds their rational rank from below.

The counts: on the uniform ring at N = 10 the room is 10, the five bond axes' products; at N = 12 it is 938 =
C(12, 6) + 6 * 2 + 2; on rings with J_l = J_{l+M} and no bond axis 926, with two perpendicular bond axes 930;
on a ring with one bond axis and no antipodal symmetry 926.
"""
import random
import sys
import time
from collections import defaultdict
from itertools import combinations, product
from math import comb

import flint
import numpy as np
import scipy.sparse as sps
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

FAILS = []


def check(label, ok, detail=""):
    print(("  PASS " if ok else "  FAIL ") + label + (f"   [{detail}]" if detail else ""), flush=True)
    if not ok:
        FAILS.append(label)


def subsets(n, p):
    return [A for A in range(2 ** n) if bin(A).count("1") == p]


def ring_sector(n, p, J, Jw):
    bs = subsets(n, p); ix = {b: i for i, b in enumerate(bs)}
    r, c, v = [], [], []
    for A in bs:
        for l in range(n - 1):
            if ((A >> l) & 1) != ((A >> (l + 1)) & 1):
                r.append(ix[A ^ (3 << l)]); c.append(ix[A]); v.append(J[l])
        if (A & 1) != ((A >> (n - 1)) & 1):
            r.append(ix[A ^ (1 | 1 << (n - 1))]); c.append(ix[A]); v.append(Jw)
    return sps.csr_matrix((v, (r, c)), shape=(len(bs), len(bs)), dtype=np.int64), bs


def turn(n, J, Jw):
    m = n // 2
    Hm, bm = ring_sector(n, m, J, Jw); H2, b2 = ring_sector(n, 2, J, Jw)
    K = (sps.kron(Hm, sps.identity(len(b2), dtype=np.int64)) - sps.kron(sps.identity(len(bm), dtype=np.int64), H2)).tocsr()
    return K, bm, b2, Hm, H2


def ring_from(n, spec):
    """None: uniform; a list of N couplings: J_0..J_{N-1}, the last the wrap; a list of M: repeated, J_l = J_{l+M}"""
    if spec is None:
        return [1] * (n - 1), 1
    Jall = list(spec) if len(spec) == n else list(spec) + list(spec)
    return Jall[:n - 1], Jall[n - 1]


def rank_q(A_):
    return DomainMatrix([[ZZ(int(x)) for x in row] for row in A_], A_.shape, ZZ).convert_to(sp.QQ).rank()


def kernel_int(Hd, idx, size, e=0):
    A = Hd[:, idx] - e * np.eye(Hd.shape[0], dtype=np.int64)[:, idx]
    out = []
    for v in DomainMatrix([[ZZ(int(x)) for x in row] for row in A], A.shape, ZZ).convert_to(sp.QQ).nullspace().to_Matrix().tolist():
        l_ = sp.ilcm(*[sp.Rational(x).q for x in v]); w = np.zeros(size, dtype=np.int64)
        w[idx] = [int(sp.Rational(x) * l_) for x in v]
        out.append(w)
    return out


def nmod(A_, p):
    A_ = np.asarray(A_)
    return flint.nmod_mat(A_.shape[0], A_.shape[1], [int(x) for x in (A_ % p).ravel()], p)


# ---------------- R1 ----------------
def divisions(sites):
    if not sites:
        yield []
        return
    a = sites[0]
    for i in range(1, len(sites)):
        rest = sites[1:i] + sites[i + 1:]
        for m_ in divisions(rest):
            yield [(a, sites[i])] + m_


def full_pair_states(n, J, Jw):
    """for every division into pairs: (kernel dimension of H on the one-full-pair configurations, predicted dimension)"""
    Jb = {}
    for l in range(n - 1):
        Jb[frozenset((l, l + 1))] = J[l]
    Jb[frozenset((n - 1, 0))] = Jw
    H2, b2 = ring_sector(n, 2, J, Jw); H2d = H2.toarray(); i2 = {b: i for i, b in enumerate(b2)}
    out = []
    for P_ in divisions(list(range(n))):
        partner = {}
        for a, b in P_:
            partner[a] = b; partner[b] = a
        idx = [i2[(1 << a) | (1 << b)] for a, b in P_]
        sub = H2d[:, idx]
        inside = not np.any(H2d[np.ix_(idx, idx)])            # no hop stays inside: kept = annihilated
        dim = len(idx) - (rank_q(sub) if np.any(sub) else 0)
        stag = True
        if dim == 1:                                             # the kernel vector: equal size, opposite across bonds
            v = kernel_int(H2d, idx, len(b2))[0][idx]
            stag = all(abs(x) == abs(v[0]) and x for x in v)
            pidx = {}
            for k_, (a, b) in enumerate(P_):
                pidx[a] = k_; pidx[b] = k_
            stag &= all(v[pidx[x]] == -v[pidx[y]] for e in Jb for x, y in [tuple(e)] if pidx[x] != pidx[y])
        # predicted: every bond between pairs maps to a bond of the same coupling, and the pair graph is bipartite
        kept = all(Jb.get(frozenset((partner[x], partner[y])), 0) == c
                   for e, c in Jb.items() for x, y in [tuple(e)] if partner[x] != y)
        pid = {}
        for k, (a, b) in enumerate(P_):
            pid[a] = k; pid[b] = k
        adj = defaultdict(set)
        for e in Jb:
            x, y = tuple(e)
            if pid[x] != pid[y]:
                adj[pid[x]].add(pid[y]); adj[pid[y]].add(pid[x])
        col, bip = {0: 0}, True; stack = [0]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in col:
                    col[w] = 1 - col[u]; stack.append(w)
                elif col[w] == col[u]:
                    bip = False
        connected = len(col) == len(P_)
        out.append((P_, dim if inside and stag else -1, int(kept and bip and connected)))
    return out


print("R1 the state with one pair full, for every division into pairs")
ok1, info1 = True, {}
for n, spec in ((6, None), (8, None), (10, None), (12, None), (8, [1, 2, 1, 3]), (12, [1, 2, 3, 4, 5, 6]),
                (12, [1, 2, 3, 2, 1, 5]), (12, [1, 2, 3, 4, 5, 7, 5, 4, 3, 2, 1, 6])):
    J, Jw = ring_from(n, spec)
    res = full_pair_states(n, J, Jw)
    ok1 &= all(d == pred for _, d, pred in res)
    found = [P_ for P_, d, _ in res if d]
    info1[(n, tuple(spec) if spec else "uniform")] = (len(res), len(found))
    if spec is None:
        m = n // 2
        axes = [{tuple(sorted((l, (2 * r + n - 1 - l) % n))) for l in range(n)} for r in range(m)]
        anti = [{(d, d + m) for d in range(m)}] if m % 2 == 0 else []
        ok1 &= sorted(map(sorted, map(set, found))) == sorted(map(sorted, axes + anti))
check("no hop stays inside the configurations with one pair full; H keeps a state on them exactly when the "
      "couplings between pairs map to equal couplings under the division's involution and the pairs form a "
      "bipartite graph, and it is the staggered sum; on the uniform ring exactly the N/2 bond axes and, at even M, "
      "the antipodal map (exact, every division)",
      ok1, f"(divisions, admitted) {info1}")



def staggering_exists(n, J, Jw, P_):
    """is there a sign per pair, D, with P h P = -D h D on every bond between pairs (h the spin hopping matrix)"""
    h = np.zeros((n, n), dtype=np.int64)
    for l in range(n - 1):
        h[l, l + 1] = h[l + 1, l] = J[l]
    h[0, n - 1] = h[n - 1, 0] = Jw
    iota, pid = {}, {}
    for k, (a, b) in enumerate(P_):
        iota[a], iota[b] = b, a; pid[a] = pid[b] = k
    for sg in product((1, -1), repeat=len(P_)):
        if all(h[iota[u], iota[v]] == -sg[pid[u]] * sg[pid[v]] * h[u, v]
               for u in range(n) for v in range(n) if pid[u] != pid[v]):
            return True
    return False


print("R1b couplings of both signs: the relation P h P = -D h D")
ok1b, seen1b, rng1b = True, [0, 0], random.Random(5)
for n in (6, 8):
    for trial in range(6):
        J = [rng1b.choice((1, 2, -1, -2)) for _ in range(n - 1)]; Jw = rng1b.choice((1, 2, -1, -2))
        if trial == 0:                                       # the reflection's off-axis bonds turned against their images
            J = [1, 2, 1, -2, -1] if n == 6 else [1, 2, 3, 1, -3, -2, -1]; Jw = 1
        H2, b2 = ring_sector(n, 2, J, Jw); H2d = H2.toarray(); i2 = {b_: i for i, b_ in enumerate(b2)}
        for P_ in divisions(list(range(n))):
            idx = [i2[(1 << a) | (1 << b_)] for a, b_ in P_]
            dim = len(idx) - rank_q(H2d[:, idx])
            st = staggering_exists(n, J, Jw, P_)
            ok1b &= dim == int(st) and not np.any(H2d[np.ix_(idx, idx)])
            seen1b[int(st)] += 1
check("with couplings of both signs (random rings at N = 6 and 8, every division): H keeps a sum of full pairs "
      "exactly when a staggering D of the pairs gives P h P = -D h D on every bond between pairs",
      ok1b and seen1b[1] > 0 and seen1b[0] > 0, f"(no staggering, staggering) {seen1b}")

# ---------------- R2 ----------------
def root_of_order(N, p):
    assert (p - 1) % N == 0, "the momentum route needs p = 1 mod N"
    for g in range(2, p):
        z = pow(g, (p - 1) // N, p)
        if all(pow(z, N // q, p) != 1 for q in sp.primefactors(N)):
            return z


def room_by_momentum(n, p):
    m = n // 2; J, Jw = ring_from(n, None)
    K, bm, b2, _, _ = turn(n, J, Jw)
    im = {b: i for i, b in enumerate(bm)}; i2 = {b: i for i, b in enumerate(b2)}
    nb = len(b2); dim = len(bm) * nb; mask = (1 << n) - 1
    rot = lambda A: ((A << 1) | (A >> (n - 1))) & mask
    cellrot = np.array([im[rot(A)] * nb + i2[rot(B)] for A in bm for B in b2], dtype=np.int64)
    ham = np.array([bin(A ^ B).count("1") for A in bm for B in b2])
    rep = -np.ones(dim, dtype=np.int64); tpow = np.zeros(dim, dtype=np.int64); size = {}
    for c in range(dim):
        if rep[c] >= 0:
            continue
        x, t = c, 0
        while rep[x] < 0:
            rep[x] = c; tpow[x] = t; x = cellrot[x]; t += 1
        size[c] = t
    z = root_of_order(n, p); Kc = K.tocsc(); total = 0
    for mom in range(n):
        orbs = [o for o, s in size.items() if (mom * s) % n == 0]
        oi = {o: i for i, o in enumerate(orbs)}
        Pidx = [o for o in orbs if ham[o] == m]; Qidx = [o for o in orbs if ham[o] != m]
        pi = {o: i for i, o in enumerate(Pidx)}; qi = {o: i for i, o in enumerate(Qidx)}
        nP, nQ = len(Pidx), len(Qidx)
        KPP = flint.nmod_mat(nP, nP, p); KQP = flint.nmod_mat(nQ, nP, p)
        acc = defaultdict(int)
        for o in Pidx:
            col = Kc.getcol(o)
            for cp, val in zip(col.indices, col.data):
                O2 = rep[cp]
                if O2 in oi:
                    acc[(O2, o)] = (acc[(O2, o)] + int(val) * pow(z, int(mom * tpow[cp]) % n, p)) % p
        for (O2, o), v in acc.items():
            if O2 in pi:
                KPP[pi[O2], pi[o]] = v
            else:
                KQP[qi[O2], pi[o]] = v
        W = flint.nmod_mat(nP, nP, [1 if i == j else 0 for i in range(nP) for j in range(nP)], p)
        C, r = None, nP
        while r:
            blocks = [KQP * W] if nQ else []
            if C is not None:
                blocks.append(C * (KPP * W))
            if not blocks:
                break
            ent, nr = [], 0
            for B in blocks:
                ent += [int(x) for x in B.entries()]; nr += B.nrows()
            X, nul = flint.nmod_mat(nr, r, ent, p).nullspace()
            if nul == r:
                break
            xe = X.entries()
            W = W * flint.nmod_mat(r, nul, [int(xe[i * r + j]) for i in range(r) for j in range(nul)], p); r = nul
            if not r:
                break
            A_, na = W.transpose().nullspace(); ae = A_.entries()
            C = flint.nmod_mat(na, nP, [int(ae[j * nP + i]) for i in range(na) for j in range(nP)], p) if na else None
        total += r
    return total


# ---------------- R3 ----------------
def split_prime(h, start):
    n = h.shape[0]; p = start
    while True:
        p = int(sp.nextprime(p))
        _, facs = nmod(h, p).charpoly().factor()
        if all(f.degree() == 1 for f, _ in facs):
            cols = []
            for f, _ in facs:
                lam = int(-f[0]) % p
                X, nul = nmod(h - lam * np.eye(n, dtype=np.int64), p).nullspace(); xe = X.entries()
                cols += [(lam, [int(xe[i * n + j]) for i in range(n)]) for j in range(nul)]
            if len(cols) == n:
                return p, cols


def minors(U, n, k, p):
    rows = [[l for l in range(n) if (A >> l) & 1] for A in subsets(n, k)]
    modes = list(combinations(range(n), k))
    D = np.zeros((len(rows), len(modes)), dtype=np.int64)
    for j, S in enumerate(modes):
        sub = U[:, list(S)]
        for i, R in enumerate(rows):
            D[i, j] = int(flint.nmod_mat(k, k, [int(x) for x in sub[R, :].ravel()], p).det())
    return D, modes


def room_by_fermions(n, spec, start=10 ** 9, slack=40, seed=1):
    m = n // 2; J, Jw = ring_from(n, spec)
    assert m % 2 == 0, "the wrap entry -J_w holds for an even number of particles on both sides"
    h = np.zeros((n, n), dtype=np.int64)
    for l in range(n - 1):
        h[l, l + 1] = h[l + 1, l] = J[l]
    h[0, n - 1] = h[n - 1, 0] = -Jw
    p, cols = split_prime(h, start)
    lam = np.array([c[0] for c in cols], dtype=np.int64); U = np.array([c[1] for c in cols], dtype=np.int64).T % p
    Dm, mm = minors(U, n, m, p); D2, m2 = minors(U, n, 2, p)
    Em = np.array([sum(lam[list(S)]) % p for S in mm]); E2 = np.array([sum(lam[list(S)]) % p for S in m2])
    Hm, bm = ring_sector(n, m, J, Jw); H2, b2 = ring_sector(n, 2, J, Jw)
    good = (not np.any(((Hm @ Dm) - Dm * Em[None, :]) % p) and not np.any(((H2 @ D2) - D2 * E2[None, :]) % p)
            and nmod(Dm, p).rank() == len(bm) and nmod(D2, p).rank() == len(b2))
    ham = np.array([[bin(A ^ B).count("1") for B in b2] for A in bm])
    Qa, Qb = np.nonzero(ham != m)
    groups = defaultdict(list)
    for s in range(len(mm)):
        for t in range(len(m2)):
            groups[(Em[s] - E2[t]) % p].append((s, t))
    rng = random.Random(seed); total = 0
    for st in groups.values():
        d = len(st); s_ = np.array([x for x, _ in st]); t_ = np.array([y for _, y in st])
        pick = np.arange(len(Qa)) if d + slack >= len(Qa) else np.array(rng.sample(range(len(Qa)), d + slack))
        g = d - nmod((Dm[Qa[pick]][:, s_] * D2[Qb[pick]][:, t_]) % p, p).rank()
        if g and len(pick) < len(Qa):
            g = d - nmod((Dm[Qa][:, s_] * D2[Qb][:, t_]) % p, p).rank()
        total += g
    return total, p, good


# ---------------- R4 ----------------
def jw(A, l):
    return -1 if bin(A & ((1 << l) - 1)).count("1") % 2 else 1


def phi_sp(n, p, q):
    bp, bq, bP, bQ = subsets(n, p), subsets(n, q), subsets(n, p + 1), subsets(n, q + 1)
    iP = {b: i for i, b in enumerate(bP)}; iQ = {b: i for i, b in enumerate(bQ)}
    r, c, v = [], [], []
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(n):
                if not (A >> l) & 1 and not (B >> l) & 1:
                    r.append(iP[A | 1 << l] * len(bQ) + iQ[B | 1 << l]); c.append(i * len(bq) + j); v.append(jw(A, l) * jw(B, l))
    return sps.csr_matrix((v, (r, c)), shape=(len(bP) * len(bQ), len(bp) * len(bq)), dtype=np.int64)


def splus_sp(n, p, q):
    bp, bq, bP, bQ = subsets(n, p), subsets(n, q), subsets(n, p + 1), subsets(n, q - 1)
    iP = {b: i for i, b in enumerate(bP)}; iQ = {b: i for i, b in enumerate(bQ)}
    r, c, v = [], [], []
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(n):
                if not (A >> l) & 1 and (B >> l) & 1:
                    r.append(iP[A | 1 << l] * len(bQ) + iQ[B ^ 1 << l]); c.append(i * len(bq) + j)
                    v.append(jw(A, l) * jw(B ^ 1 << l, l) * (-1) ** l)
    return sps.csr_matrix((v, (r, c)), shape=(len(bP) * len(bQ), len(bp) * len(bq)), dtype=np.int64)


def room_from_below(n, spec, p=1000000009):
    """(rank of all rows mod p, rank without the antipodal rows, ladder kept, every row at distance M and kept)"""
    m = n // 2; J, Jw = ring_from(n, spec); Jb = list(J) + [Jw]
    axes = [r for r in range(m) if all(Jb[b] == Jb[(2 * r - 2 - b) % n] for b in range(n))]
    anti = all(Jb[b] == Jb[(b + m) % n] for b in range(n))
    K, bm, b2, Hm, H2 = turn(n, J, Jw); nb = len(b2)
    ham = np.array([bin(A ^ B).count("1") for A in bm for B in b2])
    Hmd, H2d = Hm.toarray(), H2.toarray()
    lad = (phi_sp(n, m - 1, 1) @ splus_sp(n, m - 1, 1).T.tocsr()).tocsr()       # S- to (M - 1, 1), Phi to (M, 2)
    lad_kept = abs(K @ lad - lad @ Hm).sum() == 0                                 # exact integer intertwining
    ladder = [np.asarray(lad[:, j].toarray()).ravel() for j in range(lad.shape[1])] if lad_kept else []
    axis_rows = []
    for r in axes:
        mir = lambda l: (2 * r + n - 1 - l) % n
        pairs = {tuple(sorted((l, mir(l)))) for l in range(n)}
        one = [i for i, s in enumerate(bm) if all(((s >> a) & 1) != ((s >> b) & 1) for a, b in pairs)]
        sym = [i for i, s in enumerate(b2) if all(((s >> a) & 1) == ((s >> b) & 1) for a, b in pairs)]
        psis = [v for e in range(-2 * max(Jb), 2 * max(Jb) + 1) for v in kernel_int(Hmd, one, len(bm), e)]
        axis_rows += [np.outer(a, b).ravel() for a in psis for b in kernel_int(H2d, sym, nb, 0)]
    anti_rows, n_quiet = [], None
    if anti:
        im = {b: i for i, b in enumerate(bm)}; i2 = {b: i for i, b in enumerate(b2)}
        S = [im[sum(1 << (d + m * ((b >> d) & 1)) for d in range(m))] for b in range(1 << m)]
        S2 = [i2[(1 << d) | (1 << (d + m))] for d in range(m)]
        quiet = kernel_int(H2d, S2, nb); n_quiet = len(quiet)
        anti_rows = [np.outer(a, b).ravel() for a in kernel_int(Hmd, S, len(bm)) for b in quiet]
    prod = axis_rows + anti_rows
    ok = True
    if prod:
        Pr = np.array(prod); KP = (K @ Pr.T).T
        sup = np.nonzero(np.any(np.vstack([Pr, KP]), axis=0))[0]
        ok = rank_q(np.vstack([Pr, KP])[:, sup]) == rank_q(Pr[:, sup])
    allr = ladder + prod
    ok &= not np.any(np.array(allr)[:, ham != m]) if allr else True
    cols = np.nonzero(ham == m)[0]
    rk = lambda R: nmod(np.array(R)[:, cols], p).rank() if R else 0
    return rk(allr), rk(ladder + axis_rows), lad_kept, ok, axes, (anti, n_quiet)


print("R2 + R4 the uniform ring at N = 10: the bond axes alone")
t = time.time()
up10 = [room_by_momentum(10, p) for p in (1000000021, 1000000181)]
lo10, lo10_noanti, lad10, ok10, ax10, an10 = room_from_below(10, None)
check("N = 10, uniform: the room is 10 (from above modulo two primes by momentum, from below by the five bond "
      "axes' products, exact rows); the ladder images of (5, 0) are not kept, there is no sum of full antipodal pairs",
      up10 == [10, 10] and lo10 == 10 and not lad10 and ok10 and an10 == (True, 0) and len(ax10) == 5,
      f"up {up10}, low {lo10}, axes {ax10}, antipodal (kept, sums of full pairs) {an10}, {time.time() - t:.0f} s")

print("R2 + R3 + R4 the uniform ring at N = 12: two routes from above, one from below")
t = time.time()
up12m = room_by_momentum(12, 1000000009)
up12f, p12, good12 = room_by_fermions(12, None)
lo12, lo12_noanti, lad12, ok12, ax12, an12 = room_from_below(12, None)
check("N = 12, uniform: the room is 938 = C(12, 6) + 6 * 2 + 2, by momentum and by free fermions from above "
      "(different routes and primes) and by the ladder images, the six bond axes' products and the two antipodal "
      "products from below; without the antipodal rows 936",
      (up12m, up12f, lo12, lo12_noanti) == (938, 938, 938, 936) and good12 and lad12 and ok12
      and p12 != 1000000009 and len(ax12) == 6 and an12 == (True, 1),
      f"momentum {up12m} (p = 1000000009), fermions {up12f} (p = {p12}), low {lo12}, axes {ax12}, "
      f"antipodal (kept, sums of full pairs) {an12}, {time.time() - t:.0f} s")

print("R3 + R4 three rings at N = 12 that the rotation does not keep")
for spec, want, desc in (([1, 2, 3, 4, 5, 6], 926, "J_l = J_{l+M}, no bond axis"),
                         ([1, 2, 3, 2, 1, 5], 930, "J_l = J_{l+M}, two perpendicular bond axes"),
                         ([1, 2, 3, 4, 5, 7, 5, 4, 3, 2, 1, 6], 926, "one bond axis, no antipodal symmetry")):
    t = time.time()
    up, pf, good = room_by_fermions(12, spec)
    lo, lo_na, ladk, okr, axes, anti = room_from_below(12, spec)
    pred = comb(12, 6) + 2 * len(axes) + (2 if anti[0] else 0)
    ok_sym = {926: (len(axes), anti[0]) in ((0, True), (1, False)), 930: (len(axes), anti[0]) == (2, True)}[want]
    check(f"N = 12, {desc}: the room is {want} = C(12, 6) + 2 * {len(axes)} + {2 if anti[0] else 0}, from above by free "
          f"fermions modulo one prime, from below by exact rows",
          up == lo == pred == want and good and ladk and okr and ok_sym and (not anti[0] or anti[1] == 1),
          f"up {up} (p = {pf}), low {lo}, axes {axes}, antipodal {anti}, {time.time() - t:.0f} s")

print("grey glass ring rooms:", "ALL GREEN" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
