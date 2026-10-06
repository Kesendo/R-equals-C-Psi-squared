"""Census: when is F120's rung bound m* <= 2k + 1 tight?

F120 (docs/ANALYTICAL_FORMULAS.md, PROOF_MOMENT_TOWER_PUMP_CHANNEL.md §4): with t_j(l) = Tr(Z_l H^j) and k the first j at
which some t_j(l) is nonzero, the F117 deg-1 class gives P_{2k+1,1} = (2k+1) C(2k,k) sum_l t_k^2 > 0, so m* <= 2k + 1, where
m* is the first odd m with p_m(gamma) = Tr((A + gamma Q)^m) not identically zero (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE).

The supertrace factorization Tr(Q A^{2u}) = (-1)^u sum_l sum_j (-1)^j C(2u, j) t_j t_{2u-j} pairs t_j with t_{2u-j}, one
index at most u, so P_{m,1} = m Tr(Q A^{m-1}) vanishes for every odd m < 2k + 1. Hence the bound is tight exactly when the
gamma^1 coefficient of p_{m*} is nonzero, that is, when the deg-1 class is among the classes firing at m*.

This script reads every H = sum of distinct non-identity Pauli strings with unit coefficients, up to site permutation, at
N = 3 (pairs and triples) and N = 4 (pairs), with the tower read to j = 2^N, which by Cayley-Hamilton (t_0 = 0) decides whether any rung fires (m* read up to 9 where
none does): k, the f87 girth ell, m*, the gamma-degrees firing at m*, whether
F H F = -H for F = X^N, and checks:
  C1  m* <= 2k + 1 wherever k exists (the bound);
  C2  tight <=> gamma^1 coefficient of p_{m*} nonzero (the equivalence above);
  C3  where F H F = -H, t_ell != 0 => m* = 2 ell + 1 (the girth dichotomy; outside that setting it can fail);
  C4  under F H F = -H every t_j with j even vanishes, so k is odd and 2k + 1 = 3 mod 4;
  C5  the gamma^3 coefficient of p_5 is P_{5,3} = 20 4^N [6 sum_{|S|=3} h_S^2 + (3N - 2) sum_l c_l^2] for every Hermitian H,
      h_S the coefficient of the weight-3 Z-string Z_S and c_l that of Z_l (on the census, every coefficient is 1; also on
      random integer coefficients). Derivation: the ten words of three Q and two A give P_{5,3} = 5 Tr(Q^3 A^2) +
      5 Tr(Q A Q^2 A); with Q = sum_l Z_l (x) Z_l and A = -i(H (x) 1 - 1 (x) H^T), every term carrying the bare trace of a
      product of an odd number of Z's vanishes, and each trace reduces to 2 sum_{l,m,n} Tr(Z_l Z_m Z_n H)^2;
      an ordered triple of distinct sites hits a weight-3 string six times, and the 3N - 2 triples with a repeat give Z_n;
  C6  at k >= 3 a weight-3 Z-string gives m* = 5, and at k = 3 the bound is slack exactly when one is present.
  C7  k = 1 exactly when a single-site Z term is present.
k <= 2 is tight for every H: k = 1 iff a single-site Z term is present iff p_3 = 6 4^N sum_l c_l^2 gamma != 0
(PROOF_F87_WINDOWED_MONOMIAL_CONVERSE, the cell-free m = 3 face), and k = 2 leaves p_3 = 0, so m* >= 5 = 2k + 1. At k = 3,
m* is 5 or 7, and p_5 = P_{5,3} gamma^3 (P_{5,1} = 0 below 2k + 1, P_{5,5} = Tr Q^5 = 0 since Q's diagonal is symmetric
about 0), and with c = 0 C5 makes P_{5,3} = 120 4^N sum_{|S|=3} h_S^2, so slack exactly when a weight-3 Z-string is present (C6).
The p_m are exact integer polynomials: Tr(M^m) at integer gamma modulo four primes near 2^20 (float64 matmul is exact
while d2 p^2 < 2^53), joined by CRT and interpolated over gamma = 0..m; the CRT range exceeds the a-priori bound
4^N (||A|| + N gamma)^m at every node used (asserted). Output: simulations/results/f120_tightness_census.txt.
Run from the repo root (about 50 minutes): python simulations/f120_tightness_census.py"""
import itertools, math, os, sys
from fractions import Fraction
import numpy as np

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)
FAILS = []
def check(name, ok, detail=""):
    say(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

P1 = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
      "Y": np.array([[0, -1j], [1j, 0]], dtype=complex), "Z": np.array([[1, 0], [0, -1]], dtype=complex)}
def string_op(s):
    out = np.array([[1]], dtype=complex)
    for c in s: out = np.kron(out, P1[c])
    return out
def ham(strings):
    return sum(string_op(s) for s in strings)
def zl(N, l): return string_op("".join("Z" if i == l else "I" for i in range(N)))

PRIMES = [1048573, 1048583, 1048601, 1048609]
def crt(res):
    M = 1
    for p in PRIMES: M *= p
    x = 0
    for r, p in zip(res, PRIMES):
        Mi = M // p; x = (x + r * Mi * pow(Mi, -1, p)) % M
    return x - M if x > M // 2 else x
def generators(H, N):
    d = 2 ** N; I = np.eye(d)
    A = -1j * np.kron(H, I) + 1j * np.kron(I, H.T)
    Q = sum(np.kron(zl(N, l), zl(N, l)) for l in range(N)).real
    Ar, Ai = np.round(A.real).astype(np.int64), np.round(A.imag).astype(np.int64)
    assert np.array_equal(Ar + 1j * Ai, A)
    return Ar, Ai, np.round(Q).astype(np.int64)
def traces(Ar, Ai, Q, g, mmax):
    """exact integer Tr(M^m), m = 1..mmax, at integer gamma g (real part; the imaginary part is checked zero)"""
    per = []
    for p in PRIMES:
        R = ((Ar + g * Q) % p).astype(np.float64); Im = (Ai % p).astype(np.float64)
        PR, PI = R.copy(), Im.copy(); row = [(int(round(np.trace(PR))) % p, int(round(np.trace(PI))) % p)]
        for m in range(2, mmax + 1):
            PR, PI = (PR @ R - PI @ Im) % p, (PR @ Im + PI @ R) % p
            row.append((int(round(np.trace(PR))) % p, int(round(np.trace(PI))) % p))
        per.append(row)
    out = []
    for m in range(mmax):
        re = crt([per[i][m][0] for i in range(len(PRIMES))]); im = crt([per[i][m][1] for i in range(len(PRIMES))])
        assert im == 0, "Tr(M^m) has an imaginary part"
        out.append(re)
    return out
def poly_coeffs(vals):
    """integer polynomial through (g, vals[g]), g = 0..n-1, low degree first (exact, Newton form)"""
    n = len(vals); coef = [Fraction(v) for v in vals]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / j
    # Newton basis prod (x - g) -> monomials
    poly = [Fraction(0)] * n
    basis = [Fraction(1)]
    for i in range(n):
        for d, b in enumerate(basis): poly[d] += coef[i] * b
        nb = [Fraction(0)] * (len(basis) + 1)
        for d, b in enumerate(basis): nb[d + 1] += b; nb[d] -= i * b
        basis = nb
    assert all(c.denominator == 1 for c in poly)
    return [int(c) for c in poly]
def moments(H, N, mmax):
    """exact p_m(gamma) coefficient lists for odd m <= mmax"""
    Ar, Ai, Q = generators(H, N)
    normA = 2 * np.abs(H).sum(axis=1).max()
    nodes = range(mmax + 1)
    bound = (4 ** N) * (normA + N * mmax) ** mmax
    M = 1
    for p in PRIMES: M *= p
    assert 2 * bound < M, "CRT range too small"
    tr = {g: traces(Ar, Ai, Q, g, mmax) for g in nodes}
    return {m: poly_coeffs([tr[g][m - 1] for g in range(m + 1)]) for m in range(1, mmax + 1, 2)}
def tower(H, N, jmax):
    out = {}
    Hp = np.eye(2 ** N, dtype=complex)
    for j in range(1, jmax + 1):
        Hp = Hp @ H
        out[j] = [int(round((np.trace(zl(N, l) @ Hp)).real)) for l in range(N)]
    return out
def girth(H, N):
    if np.any(np.abs(np.diag(H)) > 1e-12): return 1
    d = 2 ** N; adj = [[j for j in range(d) if j != i and abs(H[i, j]) > 1e-12] for i in range(d)]
    best = None
    for s in range(d):
        dist = {s: 0}; par = {s: -1}; frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in dist: dist[v] = dist[u] + 1; par[v] = u; nxt.append(v)
                    elif v != par[u] and dist[v] == dist[u]:
                        c = 2 * dist[u] + 1; best = c if best is None else min(best, c)
            frontier = nxt
    return 0 if best is None else best
def canonical(strings, N):
    return min(tuple(sorted("".join(s[p] for p in perm) for s in strings)) for perm in itertools.permutations(range(N)))

def census(N, size, jmax, mcap):
    labels = ["".join(t) for t in itertools.product("IXYZ", repeat=N) if set(t) != {"I"}]
    seen = set(); rows = []
    F = string_op("X" * N)
    for combo in itertools.combinations(labels, size):
        c = canonical(combo, N)
        if c in seen: continue
        seen.add(c)
        H = ham(c)
        tw = tower(H, N, jmax)
        k = next((j for j in range(1, jmax + 1) if any(tw[j])), None)
        ell = girth(H, N)
        mmax = 2 * k + 1 if k else mcap
        pm = moments(H, N, mmax)
        mstar = next((m for m in range(1, mmax + 1, 2) if any(pm[m])), None)
        degs = tuple(d for d, c_ in enumerate(pm[mstar]) if c_) if mstar else ()
        chiral = np.array_equal(F @ H @ F, -H)
        w3 = sum(1 for t in c if set(t) <= {"I", "Z"} and t.count("Z") == 3)
        w1 = sum(1 for t in c if set(t) <= {"I", "Z"} and t.count("Z") == 1)
        rows.append(dict(H=c, k=k, ell=ell, mstar=mstar, degs=degs, tw=tw, pm=pm, chiral=chiral, w3=w3, w1=w1))
    return rows

def report(rows, N, size):
    say(f"\n=== N = {N}, {size}-string Hamiltonians, {len(rows)} up to site permutation ===")
    c1 = all(r["mstar"] is not None and r["mstar"] <= 2 * r["k"] + 1 for r in rows if r["k"])
    check(f"C1 N={N} size {size}: m* <= 2k + 1 wherever k exists", c1)
    c2 = all(((r["mstar"] == 2 * r["k"] + 1) == (1 in r["degs"])) for r in rows if r["k"])
    check(f"C2 N={N} size {size}: tight <=> the gamma^1 coefficient of p_m* is nonzero", c2)
    c3rows = [r for r in rows if r["chiral"] and r["ell"] >= 1 and r["ell"] in r["tw"] and any(r["tw"][r["ell"]])]
    c3 = all(r["mstar"] == 2 * r["ell"] + 1 for r in c3rows)
    check(f"C3 N={N} size {size}: F H F = -H and t_ell != 0 => m* = 2 ell + 1", c3, f"{len(c3rows)} rows")
    c4 = all(all(not any(r["tw"][j]) for j in r["tw"] if j % 2 == 0) for r in rows if r["chiral"])
    check(f"C4 N={N} size {size}: F H F = -H => every even rung of the tower vanishes", c4)
    c5rows = [r for r in rows if 5 in r["pm"]]
    c5 = all(r["pm"][5][3] == 20 * 4 ** N * (6 * r["w3"] + (3 * N - 2) * r["w1"]) for r in c5rows)
    check(f"C5 N={N} size {size}: P_5,3 = 20 4^N [6 (weight-3 Z-strings) + (3N - 2) (single-site Z terms)]", c5, f"{len(c5rows)} rows")
    c6rows = [r for r in rows if r["k"] and r["k"] >= 3]
    c6 = all((r["w3"] > 0) <= (r["mstar"] == 5) and (r["k"] != 3 or (r["mstar"] == 5) == (r["w3"] > 0)) for r in c6rows)
    check(f"C6 N={N} size {size}: at k >= 3 a weight-3 Z-string gives m* = 5, and at k = 3 the bound is slack exactly when one is present",
          c6, f"{len(c6rows)} rows with k >= 3, {sum(1 for r in c6rows if r['k'] >= 4)} with k >= 4")
    c7 = all((r["k"] == 1) == (r["w1"] > 0) for r in rows)
    check(f"C7 N={N} size {size}: k = 1 exactly when a single-site Z term is present", c7)
    table = {}
    for r in rows:
        if not r["k"]: key = ("no k", r["ell"], r["mstar"], r["degs"], r["chiral"])
        else: key = ("tight" if r["mstar"] == 2 * r["k"] + 1 else "slack", r["k"], r["ell"], r["mstar"], r["degs"], r["chiral"])
        table.setdefault(key, []).append(r["H"])
    say("  verdict, k, ell, m*, gamma-degrees at m*, FHF = -H : count (first example)")
    for key in sorted(table, key=str):
        say(f"  {key}: {len(table[key])}  ({'+'.join(table[key][0])})")
    return table

if __name__ == "__main__":
    t3p = report(census(3, 2, 8, 9), 3, 2)
    t3t = report(census(3, 3, 8, 9), 3, 3)
    t4p = report(census(4, 2, 16, 9), 4, 2)
    # anchors from the registry
    def one(strings, N, jmax=8):
        H = ham(strings); tw = tower(H, N, jmax)
        k = next((j for j in range(1, jmax + 1) if any(tw[j])), None)
        pm = moments(H, N, 2 * k + 1)
        ms = next(m for m in range(1, 2 * k + 2, 2) if any(pm[m]))
        return k, ms, pm[ms]
    import random
    rnd = random.Random(20261006); c5r = []
    for N in (3, 4):
        labels = ["".join(t) for t in itertools.product("IXYZ", repeat=N) if set(t) != {"I"}]
        zs = [t for t in labels if set(t) <= {"I", "Z"}]
        for _ in range(12):
            terms = set(rnd.sample(labels, rnd.randint(1, 4))) | set(rnd.sample(zs, rnd.randint(1, 3)))
            coef = {t: rnd.choice([1, 2, -1, 3, -2]) for t in terms}
            H = sum(cc * string_op(t) for t, cc in coef.items())
            w3 = sum(cc * cc for t, cc in coef.items() if set(t) <= {"I", "Z"} and t.count("Z") == 3)
            w1 = sum(cc * cc for t, cc in coef.items() if set(t) <= {"I", "Z"} and t.count("Z") == 1)
            c5r.append(moments(H, N, 5)[5][3] == 20 * 4 ** N * (6 * w3 + (3 * N - 2) * w1))
    check(f"C5 random integer coefficients, N = 3 and 4: P_5,3 = 20 4^N [6 sum h_S^2 + (3N - 2) sum c_l^2]", all(c5r), f"{sum(c5r)}/{len(c5r)}")
    k, ms, co = one(["IIY", "ZZY", "ZZZ"], 3)
    check("anchor: Y2 + Z0Z1Y2 + Z0Z1Z2 has k = 3, m* = 5, p5 = 7680 gamma^3 (F120 (a), slack)", (k, ms, co) == (3, 5, [0, 0, 0, 7680, 0, 0]),
          f"k={k}, m*={ms}, {co}")
    say("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
    os.makedirs("simulations/results", exist_ok=True)
    open("simulations/results/f120_tightness_census.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
    sys.exit(1 if FAILS else 0)
