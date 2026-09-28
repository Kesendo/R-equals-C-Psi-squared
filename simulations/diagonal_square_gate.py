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
#   S8 the cube of the three one-site squares: each letter's coordinates (k_Z, k_X, k_Y),
#      k_P = 1 iff it anticommutes with P, put I, X, Y, Z on the even tetrahedron, and the
#      face k_P = 0 is the diagonal of P's one-site square; Q_P sigma = (N - 2 k_P) sigma
#      on every string (N = 2, 3) with k = (n_X+n_Y, n_Y+n_Z, n_X+n_Z), so Pauli dephasing
#      is linear in k; the strings fill the C(N+3, 3) even lattice points of the scaled
#      tetrahedron (N = 1..6); the letter skeletons of R, D, h_zx, t_yz generate all 24
#      letter permutations, R flipping (k_Z, k_Y), the letter swaps exchanging axes, D
#      moving no letter; one-coordinate flips leave the tetrahedron (a statement); a set
#      of letter rates pairs exactly when some g_P = 0, by the half-turn about that axis,
#      and the N-site half-turn pairs every rate to 2 N sigma (isotropic control fails);
#      the centre holds N!/((N/4)!)^4 strings when 4 | N and none otherwise (N = 1..8);
#      the copy cube's disagreement patterns are the same four points; the polarity cube
#      is k mod 2 plus the transpose's n_Y = (k_Z + k_X - k_Y)/2
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

# ---------- S8: the cube of the three squares ----------
print()
print("S8 the cube: the three one-site squares as the face pairs of one cube")
LET = {"I": np.array([[1, 0], [0, 1]], dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
       "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.array([[1, 0], [0, -1]], dtype=complex)}
AXES = ("Z", "X", "Y")                       # coordinate order (k_Z, k_X, k_Y)


def as_letter(M, scale=1):
    """M = scale * phase * letter, phase in {1, -1, i, -i}; returns (letter, phase) or None."""
    for L, A in LET.items():
        for sg in (1, -1, 1j, -1j):
            if np.array_equal(M, scale * sg * A):
                return L, sg
    return None


def string_matrix(s):
    M = np.array([[1]], dtype=complex)
    for ch in s:
        M = np.kron(M, LET[ch])
    return M


# (a) the letter coordinates, read off the matrices: k_P = 1 iff the letter anticommutes with P
coord = {L: tuple(int(np.array_equal(LET[L] @ LET[P], -(LET[P] @ LET[L]))) for P in AXES) for L in LET}
EVEN_TET = {(0, 0, 0), (1, 0, 1), (1, 1, 0), (0, 1, 1)}
check("the letters sit on the even tetrahedron: I = 000, X = 101, Y = 110, Z = 011 in (k_Z, k_X, k_Y), "
      "and k_Y = k_Z xor k_X",
      coord == {"I": (0, 0, 0), "X": (1, 0, 1), "Y": (1, 1, 0), "Z": (0, 1, 1)}
      and set(coord.values()) == EVEN_TET and all(k[2] == k[0] ^ k[1] for k in coord.values()))
# the face k_P = 0 is the diagonal of P's one-site square: the letters diagonal in P's eigenbasis
H2 = np.array([[1, 1], [1, -1]], dtype=complex)          # sqrt2 * Hadamard: Z-basis -> X-basis
S2 = np.array([[1, 1], [1j, -1j]], dtype=complex)        # sqrt2 * (columns = Y eigenvectors)
BASIS = {"Z": (np.eye(2, dtype=complex), 1), "X": (H2, 2), "Y": (S2, 2)}
ok = True
for P, (B, norm) in BASIS.items():
    for L in LET:
        M = B.conj().T @ LET[L] @ B                      # = norm * (letter in P's eigenbasis)
        diagonal = M[0, 1] == 0 and M[1, 0] == 0
        antidiag = M[0, 0] == 0 and M[1, 1] == 0
        ok &= (diagonal and not antidiag) if coord[L][AXES.index(P)] == 0 else (antidiag and not diagonal)
check("each face pair is one one-site square: k_P = 0 exactly for the letters diagonal in P's "
      "eigenbasis, k_P = 1 exactly for the anti-diagonal ones", ok)

# (b) the three dephasing diagonals Q_P(sigma) = sum_l P_l sigma P_l are the cube's axes:
#     Q_P sigma = (N - 2 k_P) sigma with (k_Z, k_X, k_Y) = (n_X + n_Y, n_Y + n_Z, n_X + n_Z),
#     and Pauli dephasing sum_P gamma_P (Q_P - N) is linear in the coordinates
for N in (2, 3):
    ok = True
    for s in product("IXYZ", repeat=N):
        S = string_matrix(s)
        nX, nY, nZ = s.count("X"), s.count("Y"), s.count("Z")
        k = (nX + nY, nY + nZ, nX + nZ)
        for i, P in enumerate(AXES):
            Q = sum(string_matrix(tuple(P if j == l else "I" for j in range(N))) @ S
                    @ string_matrix(tuple(P if j == l else "I" for j in range(N))) for l in range(N))
            ok &= np.array_equal(Q, (N - 2 * k[i]) * S)
    check(f"N={N}: Q_P sigma = (N - 2 k_P) sigma on every Pauli string for P = Z, X, Y, with "
          "(k_Z, k_X, k_Y) = (n_X+n_Y, n_Y+n_Z, n_X+n_Z), so sum_P g_P (Q_P - N) is "
          "-2(g_Z k_Z + g_X k_X + g_Y k_Y) on each string", ok)

# (c) the strings fill the even lattice points of the tetrahedron scaled by N
for N in range(1, 7):
    pts = {(s.count("X") + s.count("Y"), s.count("Y") + s.count("Z"), s.count("X") + s.count("Z"))
           for s in product("IXYZ", repeat=N)}
    tet = {(a, b, c) for a in range(N + 1) for b in range(N + 1) for c in range(N + 1)
           if (a + b + c) % 2 == 0 and a + b + c <= 2 * N and a <= b + c and b <= a + c and c <= a + b}
    binom = (N + 1) * (N + 2) * (N + 3) // 6
    check(f"N={N}: the strings reach exactly the {binom} = C(N+3, 3) even lattice points of the "
          "tetrahedron with vertices I^N, X^N, Y^N, Z^N", pts == tet and len(pts) == binom,
          f"{len(pts)} points")

# (d) the symmetry: the one-site letter skeletons of R, D, h_zx, t_yz generate all 24 permutations
#     of the letters; each acts on the cube as an axis permutation followed by an even translation
T2 = LET["Y"] + LET["Z"]                                 # sqrt2 * the Y<->Z transposition Clifford
skel = {
    "R": {L: as_letter(LET[L] @ LET["X"])[0] for L in LET},
    "D": {L: as_letter(LET[L].T)[0] for L in LET},
    "h_zx": {L: as_letter(H2 @ LET[L] @ H2.conj().T, 2)[0] for L in LET},
    "t_yz": {L: as_letter(T2 @ LET[L] @ T2.conj().T, 2)[0] for L in LET},
}
group = {tuple("IXYZ")}
frontier = [tuple("IXYZ")]
while frontier:
    nxt = []
    for g in frontier:
        for m in skel.values():
            h = tuple(m[x] for x in g)
            if h not in group:
                group.add(h)
                nxt.append(h)
    frontier = nxt
check("the skeletons of R, D, h_zx, t_yz generate all 24 letter permutations (T_d on the "
      "tetrahedron)", len(group) == 24, f"{len(group)} permutations")
check("R is the translation by (1, 0, 1): it pays the Z and Y dephasing and not the X; "
      "the letter swap P <-> Q swaps the axes k_P and k_Q (h_zx: k_Z <-> k_X, t_yz: k_Z <-> k_Y); "
      "D moves no letter",
      all(coord[skel["R"][L]] == tuple(x ^ y for x, y in zip(coord[L], (1, 0, 1))) for L in LET)
      and all(skel["D"][L] == L for L in LET)
      and all(coord[skel["h_zx"][L]] == (coord[L][1], coord[L][0], coord[L][2]) for L in LET)
      and all(coord[skel["t_yz"][L]] == (coord[L][2], coord[L][1], coord[L][0]) for L in LET))
check("(a statement about the four points, not a control) flipping one coordinate, or all three, "
      "leaves the tetrahedron",
      all(not {tuple(x ^ y for x, y in zip(k, v)) for k in EVEN_TET} <= EVEN_TET
          for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]))

# (e) palindromic letter rates = a half-turn about a silent letter's axis
MATCHINGS = [(("I", "X"), ("Y", "Z")), (("I", "Y"), ("X", "Z")), (("I", "Z"), ("X", "Y"))]


def rate(k, g):
    return 2 * sum(g[P] * k[i] for i, P in enumerate(AXES))


def half_turn(k, P):   # the half-turn about P's axis flips the other two coordinates
    return tuple(c if i == AXES.index(P) else 1 - c for i, c in enumerate(k))


ok = True
for trip in product(range(4), repeat=3):
    g = dict(zip(AXES, trip))
    got = {m for m in MATCHINGS
           if rate(coord[m[0][0]], g) + rate(coord[m[0][1]], g) == rate(coord[m[1][0]], g) + rate(coord[m[1][1]], g)}
    pred = set()
    for P in AXES:
        if g[P] == 0:
            pairs = {tuple(sorted((L, next(M for M in LET if coord[M] == half_turn(coord[L], P))))) for L in LET}
            pred |= {m for m in MATCHINGS if set(m) == pairs}
    ok &= got == pred
    sig = sum(g.values())
    for P in AXES:   # the half-turn about P pairs I with P; its two pair sums are 2(sig -+ g_P)
        partner = next(M for M in LET if coord[M] == half_turn(coord["I"], P))
        rest = [M for M in LET if M not in ("I", partner)]
        ok &= partner == P and rate(coord[P], g) == 2 * (sig - g[P]) \
            and rate(coord[rest[0]], g) + rate(coord[rest[1]], g) == 2 * (sig + g[P])
check("per site, the half-turn about P's axis pairs I with P and its pair sums are 2(sigma - g_P) "
      "and 2(sigma + g_P); a set of letter rates (g_X, g_Y, g_Z in {0..3}) pairs to a constant sum exactly "
      "when some g_P = 0, and every pairing that works is the half-turn about such a letter's axis", ok)
ok = True
for P in AXES:
    g = {Q: (0 if Q == P else Fraction(random.randint(1, 9), random.randint(1, 5))) for Q in AXES}
    sig = sum(g.values())
    for N in (2, 3):
        for s in product("IXYZ", repeat=N):
            k = tuple(sum(coord[ch][i] for ch in s) for i in range(3))
            kk = tuple(sum(half_turn(coord[ch], P)[i] for ch in s) for i in range(3))
            ok &= rate(k, g) + rate(kk, g) == 2 * N * sig
check("N = 2, 3: with g_P = 0, the site-wise half-turn about P's axis pairs every string's rate "
      "with its partner's to 2 N sigma (arithmetic of the linear form, a consequence of the per-site check)", ok)
g_dep = {Q: Fraction(1, 3) for Q in AXES}
check("(a sub-case of the sweep above, read out) at isotropic depolarizing (g_X = g_Y = g_Z) no pairing of the four letters works",
      not any(rate(coord[m[0][0]], g_dep) + rate(coord[m[0][1]], g_dep)
              == rate(coord[m[1][0]], g_dep) + rate(coord[m[1][1]], g_dep) for m in MATCHINGS))

# (f) the cube's centre (N/2, N/2, N/2): n_X = n_Y = n_Z = N/4, every Q_P = 0 there
for N in range(1, 9):
    cnt = sum(1 for s in product("IXYZ", repeat=N)
              if 2 * (s.count("X") + s.count("Y")) == N and 2 * (s.count("Y") + s.count("Z")) == N
              and 2 * (s.count("X") + s.count("Z")) == N)
    if N % 4 == 0:
        q = N // 4
        f = 1
        for i in range(1, N + 1):
            f *= i
        qf = 1
        for i in range(1, q + 1):
            qf *= i
        expect = f // qf ** 4
    else:
        expect = 0
    check(f"N={N}: strings at the cube's centre: {cnt} (N!/((N/4)!)^4 when 4 | N, else 0)", cnt == expect)

# (g) the copy cube: three copies a, b, c per site disagree in the pattern (a^b, b^c, c^a),
#     which takes exactly the four values of the even tetrahedron
copy = {(a ^ b, b ^ c, c ^ a) for a, b, c in product((0, 1), repeat=3)}
check("(a statement, not a control) the copy cube: the disagreement patterns (a xor b, b xor c, c xor a) of three bits are exactly "
      "the letters' four points; no triple disagrees three times", copy == EVEN_TET)

# (h) the polarity cube of the Pi factorization section 7 is this cube read mod 2, plus the transpose
N = 3
ok = True
for s in product("IXYZ", repeat=N):
    S = string_matrix(s)
    k = tuple(sum(coord[ch][i] for ch in s) for i in range(3))
    Zn, Xn = string_matrix("Z" * N), string_matrix("X" * N)
    bit_a = 0 if np.array_equal(Zn @ S @ Zn, S) else 1
    bit_b = 0 if np.array_equal(Xn @ S @ Xn, S) else 1
    y_par = 0 if np.array_equal(S.T, S) else 1
    ok &= (bit_a, bit_b) == (k[0] % 2, k[1] % 2) and y_par == ((k[0] + k[1] - k[2]) // 2) % 2
check("N=3: bit_a = k_Z mod 2, bit_b = k_X mod 2 (the characters of Ad Z^N, Ad X^N), and the "
      "transpose's y_par = n_Y = (k_Z + k_X - k_Y)/2 mod 2, on every string", ok)

print()
if FAILURES:
    print(f"diagonal square gate: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("diagonal square gate: ALL GREEN")
