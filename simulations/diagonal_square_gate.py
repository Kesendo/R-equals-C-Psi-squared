# diagonal_square_gate.py
#
# Gate for docs/THE_ONE_SQUARE.md. On a square S x S with the label swap t
# and a label involution c ("the opposite"), the diagonal is Fix(t) and the
# anti-diagonal is Fix(t o (c x c)); the swap and the one-sided flip 1 x c generate
# the eight symmetries of the square, the swap and the two-sided flip c x c only
# four. Every square here has the four; where the repo owns a one-sided move it has
# all eight (cells, blocks, mode indices); the corner and the rate square have four.
#
# Checks (all must pass; prints "diagonal square gate: ALL GREEN"). Everything is
# exact: permutations of integer indices, Fractions, integer Pauli matrices.
#   S1 the cell square (|a><b|, N = 3, 4, 5): <R, D> (R: rho -> rho F, D: rho -> rho^T,
#      F = X^N) closes to exactly the eight named maps 1, D, R, FF (rho -> F rho F),
#      FF.D, FF.R, Pi_Z = R.D (rho -> rho^T F), Pi_Y (rho -> F rho^T); D fixes exactly
#      a = b, FF.D exactly b = not a, the other five non-identity maps fix no cell;
#      {1, FF, D, FF.D} keep the disagreement set a xor b (hence k and the dephasing)
#      and {R, FF.R, Pi_Z, Pi_Y} complement it (k to N - k), so the latter map the
#      diagonal onto the anti-diagonal
#   S1b each named map is its operator formula on every matrix unit (N = 2, 3)
#   S2 the block square (p, q) = (popcount a, popcount b), p the ket side (the cited
#      proofs write (bra, ket); every statement here is symmetric): all eight maps descend to
#      maps of (p, q), distinct, composing as the cell maps do; each named map goes to
#      its named image (D -> transpose, FF.D -> (N-q, N-p), FF -> half-turn, R and FF.R
#      -> the two folds, Pi_Z and Pi_Y -> the two quarter-turns), with the fixed sets
#      p = q, p + q = N, the centre (even N), and the two midlines (even N); on the
#      charges (p - q, p + q - N) the first half changes signs, the second exchanges them
#   S3 the site square of the frozen divisor (single-excitation cells (a, b), N = 3..8):
#      t fixes the N populations, tauQ: (a, b) -> (R b, R a) the N cells (a, R a),
#      R x R the centre at odd N only, t o tauQ = R x R; the cell's dephasing rate
#      -2 (gamma_a + gamma_b) (0 on populations) is exactly -4 gbar on every
#      off-centre anti-diagonal cell of a random rational locus profile; off the
#      locus some anti-diagonal cell misses -4 times the mean (the control); the
#      Hamiltonian part K(rho) = h rho - rho h changes sign under t and tauQ and
#      commutes with R x R and with the one-sided reflection (a, b) -> (a, R b);
#      the half-turn R x R and tauQ negate the recentred rates (rate + 4 gbar) on every cell but the
#      populations, which keep 4 gbar (the tax); R x R negates every cell at gbar = 0
#      and fails off the locus (control); the one-sided reflection neither keeps nor
#      negates them
#   S4 the rate square: on profiles x in Q^N, F71 (x -> R x) and R90 (x -> 2 avg - R x)
#      are commuting involutions with product x -> 2 avg - x, so <F71, R90> has four
#      elements and none of order four; F71 fixes the palindromic profiles and R90 the
#      profiles with x_l + x_R(l) = 2 avg, exercised on constructed members of each and
#      on random profiles; the product fixes only the uniform profile; the four are
#      distinct; control: an R90 built on the wrong mean does not fix the locus
#   S5 the other sense of the word (N = 2, 3): on Pauli strings D is the sign
#      (-1)^n_Y and FF.D the sign (-1)^n_Z, diagonal matrices; R and FF.R fix no string
#      up to sign
#   S6 controls, through the same classifiers: a one-sided flip that flips only one
#      bit of b breaks the S1 classification, and a cell map that is not built
#      from the square's moves does not descend to (p, q)
#   S7 the mode-index square of SeedRungGramClaim: G = 2 + [a = c] + [a + c = M] is
#      invariant under the swap and both one-sided flips (all eight), n = 2..8; the control drops the anti-diagonal term
#
# Runtime: under a second.
import sys
import random
from fractions import Fraction
from itertools import product

import numpy as np

random.seed(2026)
FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILURES.append(name)


def pc(x):
    return bin(x).count("1")


def named_maps(N):
    """The eight elements of <R, D> on cells (a, b), by their formulas on rho."""
    f = (1 << N) - 1
    return {
        "1": lambda a, b: (a, b),
        "D": lambda a, b: (b, a),                  # rho^T
        "R": lambda a, b: (a, b ^ f),              # rho F
        "FF": lambda a, b: (a ^ f, b ^ f),         # F rho F
        "FF.D": lambda a, b: (b ^ f, a ^ f),       # F rho^T F
        "FF.R": lambda a, b: (a ^ f, b),           # F rho
        "Pi_Z": lambda a, b: (b, a ^ f),           # rho^T F
        "Pi_Y": lambda a, b: (b ^ f, a),           # F rho^T
    }


def closure(N, gens):
    """Close the maps generated by gens (cell -> cell) as tuples over all cells."""
    cells = [(a, b) for a in range(1 << N) for b in range(1 << N)]
    idx = {c: i for i, c in enumerate(cells)}
    ident = tuple(range(len(cells)))
    gen_t = [tuple(idx[g(*c)] for c in cells) for g in gens]
    seen, frontier = {ident}, [ident]
    while frontier:
        nxt = []
        for p in frontier:
            for g in gen_t:
                q = tuple(g[p[i]] for i in range(len(cells)))
                if q not in seen:
                    seen.add(q)
                    nxt.append(q)
        frontier = nxt
    return cells, idx, seen


def classify_cells(N, cells, maps):
    """Return (ok, detail) for the S1 claims about a dict of named maps."""
    f = (1 << N) - 1
    diag = {c for c in cells if c[0] == c[1]}
    anti = {c for c in cells if c[1] == c[0] ^ f}
    keep, flip = {"1", "FF", "D", "FF.D"}, {"R", "FF.R", "Pi_Z", "Pi_Y"}
    ok, bad = True, []
    for name, m in maps.items():
        img = {c: m(*c) for c in cells}
        fixed = {c for c in cells if img[c] == c}
        k_keep = all(pc(img[c][0] ^ img[c][1]) == pc(c[0] ^ c[1]) for c in cells)
        k_flip = all(pc(img[c][0] ^ img[c][1]) == N - pc(c[0] ^ c[1]) for c in cells)
        want_fixed = {"1": set(cells), "D": diag, "FF.D": anti}.get(name, set())
        s_keep = all(img[c][0] ^ img[c][1] == c[0] ^ c[1] for c in cells)
        s_comp = all(img[c][0] ^ img[c][1] == c[0] ^ c[1] ^ f for c in cells)
        this = (fixed == want_fixed and (k_keep and s_keep if name in keep else k_flip and s_comp)
                and (name in keep or {img[c] for c in diag} == anti))
        if not this:
            bad.append(name)
        ok &= this
    return ok, bad


# ---------- S1: the cell square ----------
print("S1 the cell square |a><b|: <R, D> and its fixed cells")
for N in (3, 4, 5):
    maps = named_maps(N)
    cells, idx, group = closure(N, [maps["D"], maps["R"]])
    named = {tuple(idx[m(*c)] for c in cells) for m in maps.values()}
    ok_c, bad = classify_cells(N, cells, maps)
    check(f"N={N}: <R, D> is exactly the eight named maps; D fixes exactly a = b, FF.D "
          "exactly b = not a, the other five non-identity maps fix no cell; {1, FF, D, FF.D} "
          "keep the disagreement set a xor b (so k), {R, FF.R, Pi_Z, Pi_Y} complement it (k to "
          "N - k) and map the diagonal onto the anti-diagonal", len(group) == 8 and group == named and ok_c,
          f"{len(group)} elements; failing: {bad}")

# S1b: each named map is the operator formula its comment names, on matrix units
Fm = lambda N: np.fliplr(np.eye(1 << N, dtype=int))      # X^N permutes b -> not b
formulas = {
    "1": lambda r, F: r, "D": lambda r, F: r.T, "R": lambda r, F: r @ F,
    "FF": lambda r, F: F @ r @ F, "FF.D": lambda r, F: F @ r.T @ F, "FF.R": lambda r, F: F @ r,
    "Pi_Z": lambda r, F: r.T @ F, "Pi_Y": lambda r, F: F @ r.T,
}
for N in (2, 3):
    F = Fm(N)
    maps = named_maps(N)
    ok = True
    for name, form in formulas.items():
        for a in range(1 << N):
            for b in range(1 << N):
                E = np.zeros((1 << N, 1 << N), dtype=int)
                E[a, b] = 1
                a2, b2 = maps[name](a, b)
                W = np.zeros_like(E)
                W[a2, b2] = 1
                ok &= np.array_equal(form(E, F), W)
    check(f"N={N}: each of the eight cell maps is its operator formula (rho^T, rho F, F rho F, "
          "F rho^T F, F rho, rho^T F, F rho^T) on every matrix unit", ok)

# ---------- S2: the block square ----------
print()
print("S2 the block square (p, q): the eight maps descend, one by one")


def block_image(N, m):
    img = {}
    for a in range(1 << N):
        for b in range(1 << N):
            pq, out = (pc(a), pc(b)), tuple(pc(x) for x in m(a, b))
            if img.setdefault(pq, out) != out:
                return None
    return img


for N in (3, 4, 5):
    maps = named_maps(N)
    pts = [(p, q) for p in range(N + 1) for q in range(N + 1)]
    want = {
        "1": lambda p, q: (p, q), "D": lambda p, q: (q, p),
        "R": lambda p, q: (p, N - q), "FF": lambda p, q: (N - p, N - q),
        "FF.D": lambda p, q: (N - q, N - p), "FF.R": lambda p, q: (N - p, q),
        "Pi_Z": lambda p, q: (q, N - p), "Pi_Y": lambda p, q: (N - q, p),
    }
    imgs = {name: block_image(N, m) for name, m in maps.items()}
    ok = all(v is not None for v in imgs.values())
    ok &= all(all(imgs[n][pq] == want[n](*pq) for pq in pts) for n in maps) if ok else False
    ok &= len({tuple(sorted(v.items())) for v in imgs.values()}) == 8 if ok else False
    if ok:   # composition: the image of "g then h" is "image g then image h"
        for g in maps:
            for h in maps:
                gh = lambda a, b, g=g, h=h: maps[h](*maps[g](a, b))
                ok &= block_image(N, gh) == {pq: imgs[h][imgs[g][pq]] for pq in pts}
    # the two charges d = p - q (zero on the diagonal) and m = p + q - N (zero on the
    # anti-diagonal): the first half only changes their signs, the second exchanges them
    for n, f in want.items():
        acts = {((p - q, p + q - N), (f(p, q)[0] - f(p, q)[1], f(p, q)[0] + f(p, q)[1] - N)) for p, q in pts}
        # exact signed action on (d, m): (exchange?, sign on first, sign on second)
        law = {"1": (0, 1, 1), "D": (0, -1, 1), "FF": (0, -1, -1), "FF.D": (0, 1, -1),
               "R": (1, 1, 1), "FF.R": (1, -1, -1), "Pi_Z": (1, 1, -1), "Pi_Y": (1, -1, 1)}[n]
        ex, s1, s2 = law
        ok &= all(y == ((s1 * x[1], s2 * x[0]) if ex else (s1 * x[0], s2 * x[1])) for x, y in acts)
    fix = {n: [pq for pq in pts if want[n](*pq) == pq] for n in want}
    even = N % 2 == 0
    ok &= fix["D"] == [(p, p) for p in range(N + 1)]
    ok &= fix["FF.D"] == [(p, N - p) for p in range(N + 1)]
    ok &= fix["FF"] == ([(N // 2, N // 2)] if even else [])
    ok &= fix["R"] == ([(p, N // 2) for p in range(N + 1)] if even else [])
    ok &= fix["FF.R"] == ([(N // 2, q) for q in range(N + 1)] if even else [])
    ok &= fix["Pi_Z"] == fix["Pi_Y"] == ([(N // 2, N // 2)] if even else [])
    check(f"N={N}: all eight maps descend to (p, q), each to its named image, distinct, "
          "composing as the cell maps do; fixed sets p = q (D), p + q = N (FF.D), the "
          "centre (FF, Pi_Z, Pi_Y) and the two midlines (R, FF.R), all at even N only "
          "except the two diagonals; on the charges (p - q, p + q - N) the first half changes "
          "signs only and the second half exchanges them, each by its exact signed law", ok)

# ---------- S3: the site square of the frozen divisor ----------
print()
print("S3 the site square (a, b) of the single-excitation corner")
for N in range(3, 9):
    Rs = lambda a: N - 1 - a
    sq = [(a, b) for a in range(N) for b in range(N)]
    t = lambda c: (c[1], c[0])
    tq = lambda c: (Rs(c[1]), Rs(c[0]))
    hh = lambda c: (Rs(c[0]), Rs(c[1]))
    fix_t = [c for c in sq if t(c) == c]
    fix_tq = [c for c in sq if tq(c) == c]
    fix_h = [c for c in sq if hh(c) == c]
    ok = fix_t == [(a, a) for a in range(N)] and sorted(fix_tq) == sorted((a, Rs(a)) for a in range(N))
    ok &= fix_h == ([((N - 1) // 2, (N - 1) // 2)] if N % 2 else [])
    ok &= all(t(tq(c)) == hh(c) for c in sq)
    gbar = Fraction(random.randint(3, 9), random.randint(1, 3))
    g = [None] * N
    for a in range(N // 2):
        d = Fraction(random.randint(-9, 9), random.randint(1, 5))
        g[a], g[Rs(a)] = gbar + d, gbar - d
    if N % 2:
        g[N // 2] = gbar
    rate = lambda c, g: Fraction(0) if c[0] == c[1] else -2 * (g[c[0]] + g[c[1]])
    ok &= all(rate(c, g) == 0 for c in fix_t)
    ok &= all(rate(c, g) == -4 * gbar for c in fix_tq if c[0] != c[1])
    g_off = [x + (Fraction(1, 7) if i == 0 else 0) for i, x in enumerate(g)]
    m_off = sum(g_off) / N
    ok &= any(rate(c, g_off) != -4 * m_off for c in fix_tq if c[0] != c[1])
    # the one-sided reflection (a, b) -> (a, R b): it commutes with the corner's Hamiltonian
    # part rho -> h rho - rho h (h the single-excitation hopping with an R-invariant
    # diagonal), but not with the rates
    h = np.zeros((N, N), dtype=int)
    for a in range(N - 1):
        h[a, a + 1] = h[a + 1, a] = 2
    h[0, 0] = h[N - 1, N - 1] = -2
    PR = np.fliplr(np.eye(N, dtype=int))
    rho = np.array([[random.randint(-5, 5) for _ in range(N)] for _ in range(N)])
    K = lambda r: h @ r - r @ h
    ok &= np.array_equal(K(rho @ PR), K(rho) @ PR)
    ok &= np.array_equal(K(rho.T), -K(rho).T)                       # t turns K over
    ok &= np.array_equal(K(PR @ rho @ PR), PR @ K(rho) @ PR)           # R x R keeps K
    ok &= np.array_equal(K((PR @ rho @ PR).T), -(PR @ K(rho) @ PR).T)  # tauQ turns K over
    delta = [x - gbar for x in g]
    ok &= any(delta[a] + delta[b] != delta[a] + delta[Rs(b)] for a in range(N) for b in range(N) if a != b and a != Rs(b))
    # the corner's paying move is the half-turn R x R: on the locus it negates the
    # recentred rate (rate + 4 gbar) on every cell except the populations, which keep
    # 4 gbar (the frozen divisor's tax); at gbar = 0 that value is 0 and it negates all
    rec = lambda c, g, gb: rate(c, g) + 4 * gb
    ok &= all(rec(hh(c), g, gbar) == -rec(c, g, gbar) for c in sq if c[0] != c[1])
    ok &= all(rec(c, g, gbar) == 4 * gbar and rec(hh(c), g, gbar) == rec(c, g, gbar) for c in fix_t)
    # tauQ has the same rate parity (t fixes the populations and the rates are t-symmetric):
    # odd off the diagonal, even on the populations; that evenness is the proof's tax
    ok &= all(rec(tq(c), g, gbar) == -rec(c, g, gbar) for c in sq if c[0] != c[1])
    ok &= all(rec(tq(c), g, gbar) == rec(c, g, gbar) for c in fix_t)
    # control: off the locus the half-turn no longer negates the off-diagonal rates
    m_off2 = sum(g_off) / N
    ok &= any(rec(hh(c), g_off, m_off2) != -rec(c, g_off, m_off2) for c in sq if c[0] != c[1])
    # the one-sided reflection neither keeps nor negates the recentred rates
    one = lambda c: (c[0], Rs(c[1]))
    ok &= any(rec(one(c), g, gbar) != -rec(c, g, gbar) for c in sq if c[0] != c[1] and one(c)[0] != one(c)[1])
    g0 = [x - gbar for x in g]                      # the same offsets at mean zero
    ok &= all(rec(hh(c), g0, 0) == -rec(c, g0, 0) for c in sq)
    check(f"N={N}: t fixes the N populations, tauQ the N cells (a, R a), R x R the centre "
          "at odd N only, t o tauQ = R x R; on the locus every off-centre anti-diagonal "
          "cell has dephasing rate -4 gbar exactly; off the locus one misses -4 times the mean; the "
          "one-sided reflection (a, b) -> (a, R b) commutes with the Hamiltonian part and moves the rates; "
          "the half-turn R x R and tauQ negate the recentred rates except on the populations (4 gbar), "
          "R x R on all cells at gbar = 0, and not off the locus; the one-sided reflection neither "
          "keeps nor negates them", ok)

# ---------- S4: the rate square ----------
print()
print("S4 the rate square: F71 and R90 on profiles, a Klein four-group")
for N in range(3, 9):
    avg = lambda x: sum(x) / len(x)
    f71 = lambda x: x[::-1]
    r90 = lambda x: [2 * avg(x) - v for v in x[::-1]]
    half = lambda x: [2 * avg(x) - v for v in x]
    ok = True
    for _ in range(20):
        x = [Fraction(random.randint(-20, 20), random.randint(1, 6)) for _ in range(N)]
        ok &= f71(f71(x)) == x and r90(r90(x)) == x and half(half(x)) == x
        ok &= f71(r90(x)) == r90(f71(x)) == half(x)
        ok &= avg(r90(x)) == avg(x)
        pal = [x[min(l, N - 1 - l)] for l in range(N)]
        anti = [avg(x) + (x[l] - x[N - 1 - l]) / 2 for l in range(N)]
        ok &= f71(pal) == pal and r90(anti) == anti
        ok &= all(anti[l] + anti[N - 1 - l] == 2 * avg(anti) for l in range(N))
        ok &= (r90(x) == x) == all(x[l] + x[N - 1 - l] == 2 * avg(x) for l in range(N))
        ok &= (f71(x) == x) == all(x[l] == x[N - 1 - l] for l in range(N))
        ok &= (half(x) == x) == (len(set(x)) == 1)
    x = [Fraction(random.randint(-20, 20), random.randint(1, 6)) for _ in range(N)]
    while len({tuple(f71(x)), tuple(r90(x)), tuple(half(x)), tuple(x)}) < 4:
        x = [Fraction(random.randint(-20, 20), random.randint(1, 6)) for _ in range(N)]
    ok &= len({tuple(f71(x)), tuple(r90(x)), tuple(half(x)), tuple(x)}) == 4
    uni = [Fraction(3, 2)] * N
    ok &= half(uni) == uni
    # control: an R90 built on the mean of the first half only does not fix the locus
    bad_r90 = lambda x: [2 * avg(x[:N // 2]) - v for v in x[::-1]]
    loc = [Fraction(1), Fraction(-3)] + [Fraction(0)] * (N - 4) + [Fraction(3), Fraction(-1)] if N >= 4 else None
    if loc is not None:
        ok &= r90(loc) == loc and bad_r90(loc) != loc
    check(f"N={N}: F71, R90 and their product x -> 2 avg - x are commuting involutions "
          "(four elements, none of order four); F71 fixes the palindromic profiles, R90 "
          "those with x_l + x_R(l) = 2 avg; the product fixes only the uniform profile", ok)

# ---------- S5: the other sense of the word ----------
print()
print("S5 the other sense: the diagonal mirrors are diagonal matrices on Pauli strings")
P = {"I": np.array([[1, 0], [0, 1]], dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
     "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.array([[1, 0], [0, -1]], dtype=complex)}
for N in (2, 3):
    F = P["X"]
    for _ in range(N - 1):
        F = np.kron(F, P["X"])
    ok = True
    for s in product("IXYZ", repeat=N):
        M = np.array([[1]], dtype=complex)
        for ch in s:
            M = np.kron(M, P[ch])
        sign = lambda img: next((sg for sg in (1, -1) if np.array_equal(img, sg * M)), None)
        ok &= sign(M.T) == (-1) ** s.count("Y")
        ok &= sign(F @ M.T @ F) == (-1) ** s.count("Z")
        ok &= sign(M @ F) is None and sign(F @ M) is None
    check(f"N={N}: D = (-1)^n_Y and FF.D = (-1)^n_Z on every Pauli string; R and FF.R "
          "fix no string up to sign", ok)

# ---------- S6: controls ----------
print()
print("S6 controls, through the same classifiers")
N = 3
maps = named_maps(N)
cells, _, _ = closure(N, [maps["D"], maps["R"]])
bad_maps = dict(maps)
bad_maps["R"] = lambda a, b: (a, b ^ 1)            # flips one ket bit only
ok_c, bad = classify_cells(N, cells, bad_maps)
check("a one-sided flip of a single bit of b fails the S1 classification", not ok_c, f"failing: {bad}")
check("a cell map not built from the square's moves (a, b) -> (a, b xor 1) does not "
      "descend to (p, q)", block_image(N, lambda a, b: (a, b ^ 1)) is None)

# ---------- S7: the mode-index square of the seed-rung Gram ----------
print()
print("S7 the mode-index square: the scaled seed-rung Gram carries all eight")
for n in range(2, 9):
    M = n + 1
    modes = range(1, M)
    G = lambda a, c: 2 + (a == c) + (a + c == M)
    flip = lambda a: M - a
    ok = all(G(c, a) == G(a, c) and G(a, flip(c)) == G(a, c) and G(flip(a), c) == G(a, c)
             for a in modes for c in modes)
    ok &= all(G(a, a) == (4 if 2 * a == M else 3) for a in modes)
    ok &= all(G(a, M - a) == (4 if 2 * a == M else 3) for a in modes)
    G_bad = lambda a, c: 2 + (a == c)          # control: without the anti-diagonal term
    ok &= not all(G_bad(a, flip(c)) == G_bad(a, c) for a in modes for c in modes)
    check(f"n={n}: G = 2 + [a = c] + [a + c = M] (M = n + 1) is invariant under the swap and "
          "under each one-sided flip; 3 on both lines, 4 where they cross (odd n)", ok)

print()
if FAILURES:
    print(f"diagonal square gate: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("diagonal square gate: ALL GREEN")
