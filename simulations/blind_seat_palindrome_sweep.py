"""A measurement, not a proof: at N = 4 (Heisenberg bonds on K4 with signed integer couplings, some zero, the jump
Z on site 0 alone, letter fields of signed magnitudes on every site), every palindromic row that is not a colouring
has its dephased seat blind in the signed one-magnon Laplacian (a row with a site symmetry fixing
site 0 is blind for that reason alone, so those rows are counted apart) (blind = N - rank of the Krylov matrix of e_0,
exact over Q; F157). The ends are ranked mod two primes (the smaller nullity, an upper bound on the rational one).
Run: python simulations/blind_seat_palindrome_sweep.py SEED ROWS   (the stored run: seeds 1..8, 15000 rows each,
output in simulations/results/blind_seat_palindrome_sweep.txt)"""
import itertools, random, sys
import numpy as np
import sympy as sp

N = 4
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PRIMES = [2147483629, 2147483587]
LET = {'X': 1, 'Y': 2, 'Z': 3}

def mul1(p, q):
    if p == 0: return 0, q
    if q == 0: return 0, p
    if p == q: return 0, 0
    return (1 if (p, q) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - p - q

def mul(s, t):
    ph, out = 0, []
    for p, q in zip(s, t):
        k, c = mul1(p, q); ph += k; out.append(c)
    return ph % 4, tuple(out)

def anti(s, t): return sum(1 for p, q in zip(s, t) if p and q and p != q) % 2 == 1

def rank_mod(M, p):
    A = M.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        r += 1
        if r == rows: break
    return r

def ends(terms):
    out = []
    for far in (False, True):
        cols = [s for s in itertools.product(range(4), repeat=N) if (s[0] in (1, 2)) == far]
        idx, entries = {}, []
        for w, t in terms:
            for j, s in enumerate(cols):
                if anti(t, s):
                    k, u = mul(t, s)
                    if u not in idx: idx[u] = len(idx)
                    entries.append((idx[u], j, (1 if k == 1 else -1) * w))
        M = np.zeros((max(1, len(idx)), len(cols)), dtype=np.int64)
        for i, j, v in entries: M[i, j] += v
        out.append(min(len(cols) - rank_mod(M, p) for p in PRIMES))
    return out

def terms_of(J, word, hs):
    t = []
    for (p, q), w in zip(K4, J):
        for c in (1, 2, 3):
            s = [0] * N; s[p] = c; s[q] = c; t.append((w, tuple(s)))
    for l, (ch, w) in enumerate(zip(word, hs)):
        if ch != '.' and w:
            s = [0] * N; s[l] = LET[ch]; t.append((w, tuple(s)))
    return t

def connected(J):
    adj = {i: set() for i in range(N)}
    for (p, q), w in zip(K4, J):
        if w: adj[p].add(q); adj[q].add(p)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x] - seen: seen.add(y); st.append(y)
    return len(seen) == N

AXES = {'x': {'X': ('X', 1), 'Y': ('Y', -1), 'Z': ('Z', -1)}, 'y': {'X': ('X', -1), 'Y': ('Y', 1), 'Z': ('Z', -1)},
        'x+y': {'X': ('Y', 1), 'Y': ('X', 1), 'Z': ('Z', -1)}, 'x-y': {'X': ('Y', -1), 'Y': ('X', -1), 'Z': ('Z', -1)}}

def site_symmetry(J, word, hs):
    """An involutive permutation fixing site 0 and preserving the couplings, with a pi rotation about an axis
    perpendicular to z (x, y, x + y or x - y suffice for letter fields) carrying each field onto its image site."""
    w = {e: x for e, x in zip(K4, J)}
    f = [(None, 0) if ch == '.' or h == 0 else (ch, h) for ch, h in zip(word, hs)]
    for perm in itertools.permutations(range(N)):
        if perm[0] != 0 or any(perm[perm[i]] != i for i in range(N)): continue
        if any(w[tuple(sorted((perm[p], perm[q])))] != x for (p, q), x in w.items()): continue
        for R in AXES.values():
            if all(((None, 0) if f[l][0] is None else (R[f[l][0]][0], R[f[l][0]][1] * f[l][1])) == f[perm[l]]
                   for l in range(N)):
                return True
    return False

def blind(J, seat=0):
    L = sp.zeros(N, N)
    for (p, q), w in zip(K4, J):
        L[p, p] += w; L[q, q] += w; L[p, q] -= w; L[q, p] -= w
    v = sp.zeros(N, 1); v[seat] = 1
    return N - sp.Matrix.hstack(*[L ** k * v for k in range(N)]).rank()

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]))
    vals = [-3, -2, -1, 1, 2, 3, 6]
    stats = {'draws': 0, 'palindromic_not_coloured': 0, 'with_site_symmetry': 0, 'without_site_symmetry': 0,
             'seat_blind': 0, 'seat_not_blind': 0}
    for _ in range(int(sys.argv[2])):
        J = [rng.choice(vals + [0, 0]) for _ in K4]
        kind = rng.randrange(4)
        if kind == 1: p, q = rng.sample(range(6), 2); J[q] = J[p]
        elif kind == 2:
            for _ in range(2): p, q = rng.sample(range(6), 2); J[q] = J[p]
        elif kind == 3: p, q = rng.sample(range(6), 2); J[q] = -J[p]
        if not connected(J): continue
        word = ''.join(rng.choice('.XYZ') for _ in range(N))
        m = rng.choice([1, 2, 3]); hs = [rng.choice([m, -m, m, rng.choice([1, 2, 3, 5])]) for _ in range(N)]
        stats['draws'] += 1
        nr, fr = ends(terms_of(J, word, hs))
        if not (fr >= 1 and nr == fr): continue
        if any(all(ch == '.' or w == 0 or ch == p for ch, w in zip(word, hs)) for p in 'XY'): continue
        stats['palindromic_not_coloured'] += 1
        if site_symmetry(J, word, hs): stats['with_site_symmetry'] += 1
        else:
            stats['without_site_symmetry'] += 1
            print('no site symmetry:', J, word, hs, nr, fr, 'blind', blind(J), flush=True)
        if blind(J) >= 1: stats['seat_blind'] += 1
        else:
            stats['seat_not_blind'] += 1
            print('PALINDROME, NOT COLOURED, SEAT NOT BLIND:', J, word, hs, nr, fr, flush=True)
    print(stats, flush=True)
