"""Gate: a universal lower bound on the inner real degeneracies d_real(k), and the chain's commutant counts.

DEGENERACY_PALINDROME counts the real Liouvillian eigenvalues at Re = -2 gamma k (canonical book) for Heisenberg bonds under
uniform Z-dephasing. A coherence |i><j| decays at -2 gamma Hamming(i, j), so every X of pure XY-weight k (supported on
coherences at Hamming distance k) with [H, X] = 0 is a real eigenvector at -2 gamma k. Write c_k for the dimension of that
weight-k commutant; then c_k <= geometric multiplicity <= algebraic multiplicity, so d_real(k) >= c_k, with equality
exactly when no real mode at -2 gamma k mixes weights.

The bound (PROOF_UNIFORM_LAW Lemma A1 read as a count): H = sum J_ij (X_i X_j + Y_i Y_j + Z_i Z_j) = sum J_ij (2 SWAP_ij - 1), so every operator invariant under all
site permutations commutes with H, whatever the graph and the weights. Summing the Pauli strings of one letter multiset
(n_I, n_X, n_Y, n_Z) gives such an operator; distinct multisets have disjoint string supports, so they are independent; and
with n_X + n_Y = k there are (k+1)(N-k+1) of them. Hence c_k >= (k+1)(N-k+1) for every graph and every set of weights.

G1 (exact) the symmetrised multiset operators (of pure XY-weight k by construction) commute with H on the chain, ring,
   star and complete graph with random integer weights, N = 3, 4, 5; control: a single non-symmetrised string does not.
G2 (exact) c_k per joint-popcount block: exact lower bound from the multiset operators, upper bound from ranks modulo two
   primes; where they meet c_k is exact. Chain, unit weights, N = 3..8 (k <= N/2 at N = 7, 8; c_k = c_{N-k} by Pi
   conjugation): c_k = (k+1)(N-k+1) for 3 <= k <= N-3.
G2b (exact) chain, N = 5..8: the multiset operators together with the exact rational weight-2 commutant vectors of the
   (1,1) and (N-1,N-1) blocks span 3N - 5 + 2 floor(N/2) (a rank modulo p of exact elements is a lower bound), meeting
   the modular upper bound, so c_2 = 3N - 5 + 2 floor(N/2) exactly.
G2c (exact) N = 4: c_2 = 13 on the chain and 23 on the ring by ranks over Q, one below the stored d_real(2) = 14 and 24.
G4 (measured) chain N = 8, block (4,4): two real eigenvalues on the center line against a weight-4 commutant of one there.
G5 (measured reading) ring, star, complete at N = 4..6: modular c_k equals the stored d_real except the ring at N = 4, k = 2.
G3 (reading) the stored d_real rows (chain N = 2..7 from DEGENERACY_PALINDROME; ring, star, complete N = 4..6 read from
   simulations/results/topology_degeneracy_comparison.txt) sit at or above the bound in every entry.
Output: simulations/results/dreal_lower_bound_gate.txt. Run from the repo root: python simulations/dreal_lower_bound_gate.py"""
import itertools, os, random, sys
from math import comb
import numpy as np

OUT = []; FAILS = []
def say(s): print(s); OUT.append(s)
def check(name, ok, detail=""):
    say(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

P1 = {"I": np.eye(2, dtype=np.int64), "X": np.array([[0, 1], [1, 0]], dtype=np.int64),
      "Z": np.array([[1, 0], [0, -1]], dtype=np.int64)}
def string(s):
    """Pauli string as an exact complex-integer matrix (Y = i X Z written as complex)"""
    out = np.array([[1]], dtype=complex)
    for c in s:
        m = np.array([[0, -1j], [1j, 0]]) if c == "Y" else P1[c].astype(complex)
        out = np.kron(out, m)
    return out
def graph(N, kind):
    if kind == "chain": return [(i, i + 1) for i in range(N - 1)]
    if kind == "ring": return [(i, (i + 1) % N) for i in range(N)]
    if kind == "star": return [(0, i) for i in range(1, N)]
    return list(itertools.combinations(range(N), 2))
def heisenberg(N, bonds, J):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (a, b), j in zip(bonds, J):
        for P in "XYZ":
            s = ["I"] * N; s[a] = P; s[b] = P; H += j * string("".join(s))
    return H

# ---------------------------------------------------------------- G1
rnd = random.Random(20261007)
for N in (3, 4, 5):
    ok = True
    for kind in ("chain", "ring", "star", "complete"):
        bonds = graph(N, kind); J = [rnd.choice([1, 2, 3, -1, -2]) for _ in bonds]
        H = heisenberg(N, bonds, J)
        for k in range(N + 1):
            ms = [(nI, nX, nY, N - k - nI) for nX in range(k + 1) for nY in [k - nX] for nI in range(N - k + 1)]
            for nI, nX, nY, nZ in ms:
                letters = "I" * nI + "X" * nX + "Y" * nY + "Z" * nZ
                strings = set(itertools.permutations(letters))
                S = sum(string("".join(t)) for t in strings)
                ok &= np.array_equal(H @ S - S @ H, np.zeros_like(H))
    check(f"G1 N={N}: every symmetrised multiset operator commutes with H on chain, ring, star, complete (random integer weights)", ok)
H4 = heisenberg(4, graph(4, "chain"), [1, 2, 3])
X0 = string("XIII")
check("G1 control: the single string XIII does not commute with a weighted chain", not np.array_equal(H4 @ X0 - X0 @ H4, np.zeros_like(H4)))

# ---------------------------------------------------------------- G2
PRIMES = [1000003, 998244353]
def rank_mod(M, p):
    M = M % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p); M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r: M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows: break
    return r
def ck_upper(N, bonds, J, k, p):
    """sum over joint-popcount blocks of dim{X in block, Hamming(i,j) = k : [H, X] = 0} read modulo p"""
    d = 2 ** N
    Hr = np.round(heisenberg(N, bonds, J).real).astype(np.int64)
    assert np.array_equal(Hr, heisenberg(N, bonds, J))
    pop = [bin(i).count("1") for i in range(d)]
    total = 0
    for pa in range(N + 1):
        A = [i for i in range(d) if pop[i] == pa]
        for pb in range(N + 1):
            B = [j for j in range(d) if pop[j] == pb]
            cols = [(i, j) for i in A for j in B if bin(i ^ j).count("1") == k]
            if not cols: continue
            rows = {(a, b): r for r, (a, b) in enumerate((a, b) for a in A for b in B)}
            M = np.zeros((len(rows), len(cols)), dtype=np.int64)
            for c, (i, j) in enumerate(cols):
                for a in A:                       # (H X)_{a j} gets H[a, i]
                    if Hr[a, i]: M[rows[(a, j)], c] += Hr[a, i]
                for b in B:                       # -(X H)_{i b} gets -H[j, b]
                    if Hr[j, b]: M[rows[(i, b)], c] -= Hr[j, b]
            total += len(cols) - rank_mod(M, p)
    return total
ok_mid = True; ok_two = True; rows_read = []
for N in range(3, 9):
    bonds = graph(N, "chain"); J = [1] * len(bonds)
    for k in range(N + 1):
        if N >= 7 and k not in (0, 1, 2, 3, N // 2): continue
        ups = [ck_upper(N, bonds, J, k, p) for p in PRIMES]
        lower = (k + 1) * (N - k + 1)
        assert ups[0] == ups[1], (N, k, ups)
        c = ups[0]
        rows_read.append((N, k, lower, c))
        if 3 <= k <= N - 3: ok_mid &= c == lower
        if k == 2 and N >= 5: ok_two &= c == 3 * N - 5 + 2 * (N // 2)
say("G2 chain c_k read (N, k, bound, c_k): " + " ".join(f"({N},{k},{lo},{c})" for N, k, lo, c in rows_read))
check("G2 chain, unit weights, N = 3..8: c_k equals the bound (k+1)(N-k+1) for 3 <= k <= N-3 (exact: the multiset lower bound, "
      "a theorem, meets the two-prime upper bound)", ok_mid)
check("G2 chain: c_2 = 3N - 5 + 2 floor(N/2) for N = 5..8 (the upper bound modulo two primes; G2b supplies the "
      "matching lower bound)", ok_two)


# ---------------------------------------------------------------- G2b: c_2 exact for the chain, N = 5..8
import sympy as sp
def block_matrix(N, Hr, pa, pb, k):
    d = 2 ** N; pop = [bin(i).count("1") for i in range(d)]
    A = [i for i in range(d) if pop[i] == pa]; B = [j for j in range(d) if pop[j] == pb]
    cols = [(i, j) for i in A for j in B if bin(i ^ j).count("1") == k]
    rows = {(a, b): r for r, (a, b) in enumerate((a, b) for a in A for b in B)}
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    for c, (i, j) in enumerate(cols):
        for a in A:
            if Hr[a, i]: M[rows[(a, j)], c] += Hr[a, i]
        for b in B:
            if Hr[j, b]: M[rows[(i, b)], c] -= Hr[j, b]
    return M, cols
def multiset_vector(N, nX, nY, nZ):
    """the symmetrised weight-2 multiset operator as a dict (i, j) -> Gaussian integer (k = nX + nY = 2)"""
    d = 2 ** N; out = {}
    def esym(vals, r):
        e = [1] + [0] * r
        for v in vals:
            for t in range(r, 0, -1): e[t] += e[t - 1] * v
        return e[r]
    for i in range(d):
        for a, b in itertools.combinations(range(N), 2):
            j = i ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
            def el(site, P):
                bi, bj = (i >> (N - 1 - site)) & 1, (j >> (N - 1 - site)) & 1
                return 1 if P == "X" else (-1j if bi == 0 else 1j)
            xy = 0
            for pa_, pb_ in set(itertools.permutations("X" * nX + "Y" * nY)):
                xy += el(a, pa_) * el(b, pb_)
            rest = [(-1) ** ((i >> (N - 1 - l)) & 1) for l in range(N) if l not in (a, b)]
            v = xy * esym(rest, nZ)
            if v: out[(i, j)] = out.get((i, j), 0) + v
    return out
def sqrtm1(p):
    for g in range(2, 400):
        x = pow(g, (p - 1) // 4, p)
        if x * x % p == p - 1: return x
def rank_vectors_mod(vecs, p):
    ip = sqrtm1(p); keys = sorted({k for v in vecs for k in v}); idx = {k: n for n, k in enumerate(keys)}
    M = np.zeros((len(vecs), len(keys)), dtype=np.int64)
    for r, v in enumerate(vecs):
        for k_, z in v.items():
            z = complex(z); M[r, idx[k_]] = (int(round(z.real)) + int(round(z.imag)) * ip) % p
    return rank_mod(M.T.copy(), p) if M.shape[0] > M.shape[1] else rank_mod(M, p)
ok2b = True; rows2b = []
for N in range(5, 9):
    bonds = graph(N, "chain"); Hc = heisenberg(N, bonds, [1] * len(bonds)); Hr = np.round(Hc.real).astype(np.int64)
    vecs = []
    for nX in range(3):
        nY = 2 - nX
        for nZ in range(N - 1):
            v = multiset_vector(N, nX, nY, nZ)
            Sm = np.zeros((2 ** N, 2 ** N), dtype=complex)
            for (i, j), z in v.items(): Sm[i, j] = z
            assert np.array_equal(Hc @ Sm - Sm @ Hc, np.zeros_like(Hc)), "multiset operator does not commute"
            vecs.append(v)
    corner = 0
    for pa in (1, N - 1):
        M, cols = block_matrix(N, Hr, pa, pa, 2)
        for ns in sp.Matrix(M).nullspace():
            den = sp.ilcm(*[sp.fraction(x)[1] for x in ns]) if any(x != 0 for x in ns) else 1
            v = {cols[c]: int(ns[c] * den) for c in range(len(cols)) if ns[c] != 0}
            Sm = np.zeros((2 ** N, 2 ** N), dtype=complex)
            for (i, j), z in v.items(): Sm[i, j] = z
            assert np.array_equal(Hc @ Sm - Sm @ Hc, np.zeros_like(Hc)), "corner vector does not commute"
            vecs.append(v); corner += 1
    lower = max(rank_vectors_mod(vecs, p) for p in (998244353, 469762049))  # a rank modulo p never exceeds the rank over Q(i)
    upper = ck_upper(N, bonds, [1] * len(bonds), 2, PRIMES[0])
    rows2b.append((N, lower, upper, corner))
    ok2b &= lower == upper == 3 * N - 5 + 2 * (N // 2)
say("G2b chain c_2 (N, lower from exact elements, upper modulo p, corner-commutant vectors): " + " ".join(map(str, rows2b)))
check("G2b chain, N = 5..8: the multiset operators and the exact (1,1) and (N-1,N-1) weight-2 commutant vectors span "
      "3N - 5 + 2 floor(N/2), meeting the upper bound: c_2 exact", ok2b)

# ---------------------------------------------------------------- G2c: the ring at N = 4 also has a weight-mixing real mode
def c2_exact(N, kind):
    Hr = np.round(heisenberg(N, graph(N, kind), [1] * len(graph(N, kind))).real).astype(np.int64); c = 0
    for pa in range(N + 1):
        for pb in range(N + 1):
            M, cols = block_matrix(N, Hr, pa, pb, 2)
            if cols: c += len(cols) - sp.Matrix(M).rank()
    return c
c2c, c2r = c2_exact(4, "chain"), c2_exact(4, "ring")
check("G2c N = 4: c_2 = 13 on the chain and 23 on the ring exactly (ranks over Q), one below the stored d_real(2) = 14 and 24",
      (c2c, c2r) == (13, 23), f"chain {c2c}, ring {c2r}")

# ---------------------------------------------------------------- G4: the chain's center block at N = 8 (measured)
N = 8; d = 2 ** N; H8 = heisenberg(N, graph(N, "chain"), [1] * (N - 1))
A = [i for i in range(d) if bin(i).count("1") == 4]; n = len(A); g_ = 0.05
Hs = H8[np.ix_(A, A)]; Iblk = np.eye(n)
L44 = -1j * (np.kron(Hs, Iblk) - np.kron(Iblk, Hs.T)) - 2 * g_ * np.diag([bin(a ^ b).count("1") for a in A for b in A])
ev = np.linalg.eigvals(L44); dist = np.sort(np.abs(ev + 2 * g_ * 4))
on_line = int(np.sum(dist < 1e-12)); gap = dist[on_line] if on_line < len(dist) else float("inf")
M44, _ = block_matrix(N, np.round(H8.real).astype(np.int64), 4, 4, 4)
c4_block = M44.shape[1] - rank_mod(M44.copy(), PRIMES[1])
check("G4 chain N = 8, block (4,4): two eigenvalues on the center line (within 1e-12 of -8 gamma, the next at least 1e-4 away; "
      "float, measured) against a weight-4 commutant of one there (modulo p): one weight-mixing real mode",
      on_line == 2 and gap > 1e-4 and c4_block == 1, f"on line {on_line}, gap {gap:.2e}, commutant {c4_block}")

# ---------------------------------------------------------------- G5: ring, star, complete at N = 4..6 (measured reading)
import re
_txt = open("simulations/results/topology_degeneracy_comparison.txt", encoding="utf-8").read()
rows5 = {kind: {} for kind in ("ring", "star", "complete")}
for _part in re.split(r"\n=+\nN = ", _txt)[1:]:
    _N = int(_part.split("\n", 1)[0])
    for kind in rows5:
        _m = re.search(rf"\n\s*{kind}\s*: d_real = \[([0-9, ]+)\]", _part)
        if _m and _N in (4, 5, 6): rows5[kind][_N] = [int(x) for x in _m.group(1).split(",")]
assert all(len(v) == 3 for v in rows5.values()), "could not read every topology row"
mismatch = []
for kind, byN in rows5.items():
    for N, row in byN.items():
        bonds = graph(N, kind)
        for k in range(N + 1):
            c = ck_upper(N, bonds, [1] * len(bonds), k, PRIMES[1])
            if c != row[k]: mismatch.append((kind, N, k, c, row[k]))
say(f"G5 modular c_k against the stored d_real on ring, star, complete, N = 4..6: differences (kind, N, k, c_k, d_real) = {mismatch}")
check("G5 the modular commutant counts equal the stored d_real everywhere except the ring at N = 4, k = 2 (23 against 24)",
      mismatch == [("ring", 4, 2, 23, 24)])

# ---------------------------------------------------------------- G3
stored = {"chain": {2: [3, 4, 3], 3: [4, 6, 6, 4], 4: [5, 8, 14, 8, 5], 5: [6, 10, 14, 14, 10, 6],
                    6: [7, 12, 19, 16, 19, 12, 7], 7: [8, 14, 22, 20, 20, 22, 14, 8]}}
import re
txt = open("simulations/results/topology_degeneracy_comparison.txt", encoding="utf-8").read()
for kind in ("ring", "star", "complete"): stored[kind] = {}
for part in re.split(r"\n=+\nN = ", txt)[1:]:
    N = int(part.split("\n", 1)[0])
    for kind in ("ring", "star", "complete"):
        m = re.search(rf"\n\s*{kind}\s*: d_real = \[([0-9, ]+)\]", part)
        if m and N in (4, 5, 6): stored[kind][N] = [int(x) for x in m.group(1).split(",")]
say("G3 rows read from the topology comparison: " + str({k: v for k, v in stored.items() if k != "chain"}))
assert all(len(stored[k]) == 3 for k in ("ring", "star", "complete")), "could not read every topology row" 
ok3 = all(v >= (k + 1) * (N - k + 1) for kind in stored for N, row in stored[kind].items() for k, v in enumerate(row))
eq = sorted({(kind, N, k) for kind in stored for N, row in stored[kind].items() for k, v in enumerate(row) if v == (k + 1) * (N - k + 1)})
check("G3 every stored d_real entry (chain N = 2..7 as in DEGENERACY_PALINDROME, ring/star/complete N = 4..6 from the topology comparison) sits at or above "
      "(k+1)(N-k+1)", ok3, f"{len(eq)} entries meet it")

say("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
os.makedirs("simulations/results", exist_ok=True)
open("simulations/results/dreal_lower_bound_gate.txt", "w", encoding="utf-8", newline="\n").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
