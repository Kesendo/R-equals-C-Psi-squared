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
#   S9 the Hamiltonian on the cube: every two-site Pauli term sends a letter pair to 0 or to
#      2*phase*(Pa)(x)(Qb) and changes the weight w = (k_Z + k_X + k_Y)/2 by exactly one; a
#      bond P(x)P steps +-2 along one axis, never its own; on XX+YY the delta k_Z = 0 moves are
#      the k_X and k_Y steps; a one-site field keeps w and steps on a face diagonal; at N = 3
#      every entry of ad_H for XYZ chains is one axis step (control: an X field adds face
#      diagonals); the quarter-turn about Z exchanges k_X and k_Y and commutes with L for
#      Heisenberg, XY and XXZ under Z-dephasing (control: X field), fixes Z_1 and Z_1 Z_2 and
#      sends X_1 Y_2 to -Y_1 X_2, whose centroid leaves the plane (float control); the Hermitian part of L in
#      the string basis is -2 gamma diag(k_Z) exactly; the populations and the XOR strings are
#      the opposite edges I^N - Z^N and X^N - Y^N, exchanged by R's letter flip
#   S10 the cells and the strings, two bases of one pattern: at every disagreement pattern x the
#      products X^x Z^z and the cells |b xor x><b| are related by a Walsh transform,
#      X^x Z^z = sum_b (-1)^(z.b) |b xor x><b| (N = 1..4; control: the sign read on the ket fails);
#      the Pauli string is i^|x and z| X^x Z^z, and inside a pattern the two bases are mutually
#      unbiased (N = 1..3); one site is PROOF_CHROMATICITY's Lemma 2 (control: the other
#      orientation fails); the quarter-turn about Z multiplies the cell |a><b| by i^(p - q),
#      p = popcount a (N = 1..4; control: i^(q - p) fails), and its square, the conjugation by Z^N,
#      is (-1)^(p - q) on the cells and (-1)^k_Z on the strings; at N = 3, which processes are
#      diagonal on the strings, diagonal on the cells, and keep (p, q), against a fixed table whose
#      every column holds both answers; correlated Z noise with real covariances sends a string to
#      -4 C_kl Z_k Z_l sigma per pair of lit sites and does nothing else (control: counting pairs with
#      any lit site fails), and collective dephasing has largest off-diagonal entry 4 gamma; a string
#      spreads over several blocks (N = 3) and a cell over every point of its height plane,
#      (k+1)(N-k+1) of them (N = 1..4); amplitude damping's off-diagonal part is the jump |1><1| -> |0><0| on the cells and
#      the arrow I -> Z (+rate) on the strings, two different operators; F98's Krawtchouk expansion is
#      the transform at x = 0 with |S| = k_X = k_Y, and X^N Z_S is (-i)^|S| times the string with Y on S
#      and X elsewhere (control: the Krawtchouk number without its sign fails); the doubled site of
#      PROOF_FROZEN_BAND_SO4 at one site: two commuting su(2) (eta on the populations, spin on the
#      coherences), the letters as their x-eigenvectors, X rho X = 2(eta^x + S^x), Y rho Y =
#      2(eta^x - S^x), Z rho Z = 1 - 8 (S^z)^2, the quarter-turn about Z = exp(i pi S^z), R exchanging
#      the doublets, D = 1 on the eta doublet and 2 S^x on the spin one (these are the local, string-free
#      operators; the frozen-band proof's ladders carry Jordan-Wigner strings); at N = 3 correlated Z noise
#      is -8 sum_{k,l} a_kl S^z_k S^z_l, the k = l terms the local price and the k != l terms the exchange;
#      the quarter-turn's Z4 on the cells has the multiplicities sum_{d = r mod 4} C(2N, N + d), the multiset
#      of ONE_FOUR_THESIS's Majorana-degree count (N = 1..8), while the two gradings differ and agree mod 2
#
# Runtime: under a second.
import sys
import random
from fractions import Fraction
from itertools import product
from math import comb

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
print("S9 the Hamiltonian on the cube: how ad_H moves a string between lattice points")


def kvec(s):
    return tuple(sum(coord[ch][i] for ch in s) for i in range(3))


def weight(s):
    return sum(ch != "I" for ch in s)


# (a) one bond term P(x)Q on one letter pair a(x)b: the commutator is 0 or 2*phase*(Pa)(x)(Qb),
#     exactly; when it is not 0 the weight changes by exactly one, and for P = Q the step is
#     +-2 along exactly one axis, never along P's own
ok_form, ok_w, ok_axis, ok_nonzero = True, True, True, True
steps_by_term = {}
for P, Q in product("XYZ", repeat=2):
    T = string_matrix(P + Q)
    steps = set()
    for a, b in product("IXYZ", repeat=2):
        S = string_matrix(a + b)
        C = T @ S - S @ T
        if not C.any():
            continue
        pa, pb = as_letter(LET[P] @ LET[a]), as_letter(LET[Q] @ LET[b])
        target = pa[0] + pb[0]
        found = any(np.array_equal(C, 2 * sg * string_matrix(target)) for sg in (1, -1, 1j, -1j))
        ok_form &= found
        dk = tuple(x - y for x, y in zip(kvec(target), kvec(a + b)))
        steps.add(dk)
        ok_w &= abs(weight(target) - weight(a + b)) == 1
        if P == Q:
            nz = [i for i in range(3) if dk[i] != 0]
            ok_axis &= len(nz) == 1 and abs(dk[nz[0]]) == 2 and AXES[nz[0]] != P
    steps_by_term[P + Q] = steps
    ok_nonzero &= len(steps) > 0
check("every two-site Pauli term P(x)Q sends a letter pair to 0 or to 2*phase*(Pa)(x)(Qb), entry by "
      "entry (all 9 terms, all 16 pairs)", ok_form and ok_nonzero)
check("every nonzero move of a two-site term changes the weight w = (k_Z + k_X + k_Y)/2 by exactly one",
      ok_w)
check("a bond term P(x)P steps by +-2 along exactly one axis, and never along its own axis P",
      ok_axis and all({i for d in steps_by_term[P + P] for i in range(3) if d[i]} == {i for i in range(3) if AXES[i] != P}
                      for P in "XYZ"))
# the old "n_XY changes by 0 or +-2": on an XX+YY bond the 0 are the k_X and k_Y steps
xy = steps_by_term["XX"] | steps_by_term["YY"]
check("on an XX+YY bond the moves with delta k_Z = 0 are exactly the +-2 steps along k_X or k_Y, and the "
      "delta k_Z = +-2 moves are the k_Z-axis steps",
      {d for d in xy if d[0] == 0} == {(0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)}
      and {d for d in xy if d[0] != 0} == {(2, 0, 0), (-2, 0, 0)})
# a one-site field P keeps the weight and moves on the face diagonal inside the plane k_P = const
ok_f = True
for P in "XYZ":
    for a in "IXYZ":
        C = LET[P] @ LET[a] - LET[a] @ LET[P]
        if not C.any():
            continue
        t = as_letter(LET[P] @ LET[a])[0]
        dk = tuple(x - y for x, y in zip(coord[t], coord[a]))
        i = AXES.index(P)
        ok_f &= a != "I" and t != "I" and dk[i] == 0 and sorted(abs(x) for x in dk) == [0, 1, 1] and sum(dk) == 0
check("a one-site field P keeps the weight and steps along a face diagonal (+-1, -+1) in the plane k_P = const",
      ok_f)

# (b) the whole ad_H at N = 3 on the Pauli strings, integer matrices: every nonzero entry of an
#     XYZ-type chain (any integer bond weights, no field) joins two strings one +-2 axis step apart;
#     control: a transverse field adds weight-keeping face-diagonal entries
N = 3
strs = ["".join(t) for t in product("IXYZ", repeat=N)]
B0 = np.array([string_matrix(s).reshape(-1) for s in strs]).T      # entries 0, +-1, +-i
d = 2 ** N


def pauli_super(M):
    """a superoperator on vec (row-major) written in the string basis; exact since d is a power of 2"""
    return (B0.conj().T @ M @ B0) / d


def site_op(P, l):
    return string_matrix("".join(P if k == l else "I" for k in range(N)))


def xyz_chain(jx, jy, jz):
    return sum(c * site_op(P, i) @ site_op(P, i + 1) for i in range(N - 1) for c, P in ((jx, "X"), (jy, "Y"), (jz, "Z")))


Id = np.eye(d)


def adH(H):
    return pauli_super(np.kron(H, Id) - np.kron(Id, H.T))


def moves_of(A):
    return {tuple(x - y for x, y in zip(kvec(strs[i]), kvec(strs[j]))) for i, j in zip(*np.nonzero(A))}


ok_chain = True
for jx, jy, jz in [(1, 1, 1), (1, 1, 0), (1, 1, 2), (0, 0, 1), (3, -2, 5)]:
    mv = moves_of(adH(xyz_chain(jx, jy, jz)))
    ok_chain &= all(sorted(abs(x) for x in m) == [0, 0, 2] for m in mv) and len(mv) > 0
check("N=3: every nonzero entry of ad_H for an XYZ chain (Heisenberg, XY, XXZ at Delta = 2, Ising, "
      "(3, -2, 5)) joins strings one +-2 axis step apart", ok_chain)
mv_f = moves_of(adH(xyz_chain(0, 0, 1) + sum(site_op("X", l) for l in range(N))))
check("control: an X field on the Ising chain adds face-diagonal moves of weight change 0",
      any(sorted(abs(x) for x in m) == [0, 1, 1] and sum(m) == 0 for m in mv_f))

# (c) the mirror plane k_X = k_Y: the quarter-turn about Z, rho -> S rho S^dag with S = diag(1, i)^N,
#     is a signed permutation of the strings exchanging k_X and k_Y, and it commutes exactly with
#     L = -i ad_H + gamma (Q_Z - N) for Heisenberg, XY and XXZ; an X field breaks it (control)
S1 = np.array([[1, 0], [0, 1j]])
SN = np.array([[1]], dtype=complex)
for _ in range(N):
    SN = np.kron(SN, S1)
Rot = pauli_super(np.kron(SN, SN.conj()))
ok_perm = True
for j in range(len(strs)):
    col = Rot[:, j]
    nz = np.nonzero(col)[0]
    ok_perm &= len(nz) == 1 and abs(col[nz[0]]) == 1
    i = nz[0]
    kj, ki = kvec(strs[j]), kvec(strs[i])
    ok_perm &= ki == (kj[0], kj[2], kj[1])
check("N=3: the quarter-turn about Z is a signed permutation of the strings that exchanges k_X and k_Y",
      ok_perm)
g = Fraction(1, 4)
QZ = sum(np.kron(site_op("Z", l), site_op("Z", l).conj()) for l in range(N))
LD = pauli_super(float(g) * (QZ - N * np.kron(Id, Id)))


def L_of(H):
    return -1j * adH(H) + LD


ok_c = all(np.array_equal(Rot @ L_of(xyz_chain(*j)) - L_of(xyz_chain(*j)) @ Rot, np.zeros_like(Rot))
           for j in [(1, 1, 1), (1, 1, 0), (1, 1, 2)])
Hf = xyz_chain(0, 0, 1) + sum(site_op("X", l) for l in range(N))
ctrl = not np.array_equal(Rot @ L_of(Hf) - L_of(Hf) @ Rot, np.zeros_like(Rot))
check("N=3: it commutes with L entry by entry for Heisenberg, XY and XXZ under Z-dephasing "
      "(gamma = 1/4); control: with an X field on the Ising chain it does not", ok_c and ctrl)

# the quarter-turn fixes the populations Z_1 and Z_1 Z_2, so their whole trajectory stays fixed and
# symmetric about the plane; a string that merely lies on the plane need not be fixed: X_1 Y_2 goes
# to -Y_1 X_2, and its centroid leaves the plane (a float reading on the Heisenberg chain)
e = lambda st: np.eye(len(strs))[strs.index(st)].astype(complex)
fixed = all(np.array_equal(Rot @ e(st), e(st)) for st in ("ZII", "ZZI"))
xy_moves = np.array_equal(Rot @ e("XYI"), -e("YXI"))
check("N=3: the quarter-turn fixes Z_1 and Z_1 Z_2 exactly and sends X_1 Y_2 (on the plane too) to -Y_1 X_2",
      fixed and xy_moves)
Lh = L_of(xyz_chain(1, 1, 1))
w_, V_ = np.linalg.eig(Lh)
Vi = np.linalg.inv(V_)
Kx = np.array([kvec(s_)[1] for s_ in strs])
Ky = np.array([kvec(s_)[2] for s_ in strs])


def plane_gap(st):
    c0 = Vi @ e(st)
    worst = 0.0
    for t in np.linspace(0.0, 8.0, 81):
        v = V_ @ (np.exp(w_ * t) * c0)
        pw = np.abs(v) ** 2
        worst = max(worst, abs((pw @ Kx - pw @ Ky) / pw.sum()))
    return worst


gz, gxy = plane_gap("ZII"), plane_gap("XYI")
check("control (float reading, Heisenberg N=3, gamma = 1/4): the centroid of X_1 Y_2 leaves the plane "
      "k_X = k_Y; Z_1's stays on it by the exact invariance above (its float value read, not gated)",
      gxy > 0.1, f"max |<k_X> - <k_Y>| = {gxy:.3f} for X_1 Y_2, {gz:.1e} for Z_1")

# (d) the budget: in the string basis the Hermitian part of L is -2 gamma diag(k_Z) exactly, so
#     d/dt ||rho||^2 = -4 gamma <k_Z> ||rho||^2 and ln ||rho(t)||^2 = -4 gamma int <k_Z> dt
ok_b = True
for j in [(1, 1, 1), (1, 1, 0), (3, -2, 5)]:
    L = L_of(xyz_chain(*j))
    herm = (L + L.conj().T) / 2
    ok_b &= np.array_equal(herm, np.diag([-2 * float(g) * kvec(s)[0] for s in strs]).astype(complex))
Lf = L_of(Hf)
ok_b &= np.array_equal((Lf + Lf.conj().T) / 2, np.diag([-2 * float(g) * kvec(s)[0] for s in strs]).astype(complex))
check("N=3: the Hermitian part of L in the string basis is exactly -2 gamma diag(k_Z), with or without a "
      "field (the Hamiltonian part is skew)", ok_b)

# (e) the two ends of the palindrome are opposite edges: the populations {I, Z}^N are the strings
#     with k_Z = 0 and lie on the edge I^N - Z^N; the XOR strings {X, Y}^N are those with k_Z = N, on the
#     edge X^N - Y^N; R's letter flip (I <-> X, Y <-> Z) carries the one edge onto the other
ok_e = True
for M in range(1, 6):
    for s in product("IXYZ", repeat=M):
        k = kvec(s)
        ok_e &= (k[0] == 0) == all(ch in "IZ" for ch in s)
        ok_e &= (k[0] == M) == all(ch in "XY" for ch in s)
        if k[0] == 0:
            ok_e &= k[1] == k[2]                      # on the segment from (0,0,0) to (0,N,N)
        if k[0] == M:
            ok_e &= k[1] + k[2] == M                  # on the segment from (N,0,N) to (N,N,0)
        flip = "".join({"I": "X", "X": "I", "Y": "Z", "Z": "Y"}[ch] for ch in s)
        ok_e &= kvec(flip) == (M - k[0], k[1], M - k[2])
check("N=1..5: populations = strings with k_Z = 0, on the edge I^N - Z^N; XOR strings = k_Z = N, on the "
      "opposite edge X^N - Y^N; R's letter flip maps (k_Z, k_X, k_Y) to (N - k_Z, k_X, N - k_Y), one edge "
      "onto the other", ok_e)

# ---------- S10: the cells and the strings, two bases of one pattern ----------
print()
print("S10 the cells and the strings: two bases of one pattern")


def xz_operator(x, z):
    """X^x Z^z as a product site by site (the Hermitian string is i^(#Y) times this, Y = i X Z)."""
    M = np.array([[1]], dtype=complex)
    for xb, zb in zip(x, z):
        M = np.kron(M, np.linalg.matrix_power(LET["X"], xb) @ np.linalg.matrix_power(LET["Z"], zb))
    return M


def cell_matrix(a, b):
    """|a><b| for bit tuples a, b; site 0 is the leftmost tensor factor, the top bit of the index."""
    M = np.zeros((2 ** len(a), 2 ** len(a)), dtype=complex)
    M[int("".join(map(str, a)), 2), int("".join(map(str, b)), 2)] = 1
    return M


def dot(u, v):
    return sum(i * j for i, j in zip(u, v))


def flip(b, x):
    return tuple(i ^ j for i, j in zip(b, x))


# (a) the Walsh transform between the strings and the cells of one disagreement pattern
ok_w, ok_wc = True, True
for N_ in range(1, 5):
    bits = list(product((0, 1), repeat=N_))
    sign_on_ket_fails = False
    for x in bits:
        for z in bits:
            lhs = xz_operator(x, z)
            ok_w &= np.array_equal(lhs, sum((-1) ** dot(z, b) * cell_matrix(flip(b, x), b) for b in bits))
            sign_on_ket_fails |= not np.array_equal(
                lhs, sum((-1) ** dot(z, flip(b, x)) * cell_matrix(flip(b, x), b) for b in bits))
    ok_wc &= sign_on_ket_fails
check("N=1..4: X^x Z^z = sum_b (-1)^(z.b) |b xor x><b| for every disagreement pattern x and every z, so the "
      "strings and the cells with disagreement set x are two bases of one 2^N space", ok_w)
check("control: reading the sign on the ket, (-1)^(z.(b xor x)), fails at every N = 1..4", ok_wc)
# the Pauli string is i^|x and z| X^x Z^z (a Y is i X Z on its site), and inside a pattern every string meets every
# cell with modulus 1 and the cells of other patterns with 0: the two bases are mutually unbiased there
ok_h, ok_mub = True, True
for N_ in range(1, 4):
    for s in product("IXYZ", repeat=N_):
        x = tuple(int(ch in "XY") for ch in s)
        z = tuple(int(ch in "ZY") for ch in s)
        Ms = string_matrix(s)
        ok_h &= np.array_equal(Ms, (1j) ** (dot(x, z) % 4) * xz_operator(x, z))
        for a in product((0, 1), repeat=N_):
            for b in product((0, 1), repeat=N_):
                entry = Ms[int("".join(map(str, a)), 2), int("".join(map(str, b)), 2)]
                ok_mub &= (abs(entry) == 1) if flip(b, x) == a else (entry == 0)
check("N=1..3: every Pauli string is i^|x and z| X^x Z^z", ok_h)
check("N=1..3: inside its pattern x every string meets every cell |b xor x><b| with modulus 1 and meets no other cell: "
      "the two bases of one pattern are mutually unbiased", ok_mub)

# (b) one site: PROOF_CHROMATICITY's Lemma 2
lemma2 = {(0, 0): (LET["I"] + LET["Z"]) / 2, (1, 1): (LET["I"] - LET["Z"]) / 2,
          (0, 1): (LET["X"] + 1j * LET["Y"]) / 2, (1, 0): (LET["X"] - 1j * LET["Y"]) / 2}
check("one site (PROOF_CHROMATICITY Lemma 2): |0><0| = (I+Z)/2, |1><1| = (I-Z)/2, |0><1| = (X+iY)/2, "
      "|1><0| = (X-iY)/2: the dark letters span the populations, the lit letters the coherences",
      all(np.array_equal(cell_matrix((a,), (b,)), M) for (a, b), M in lemma2.items()))
check("control: the coherences the other way round, |0><1| = (X-iY)/2, fail",
      not np.array_equal(cell_matrix((0,), (1,)), (LET["X"] - 1j * LET["Y"]) / 2))

# (c) the quarter-turn about Z on the cells: |a><b| -> i^(p - q) |a><b|, p = popcount a (the ket, as in §3)
ok_q, ok_qc = True, True
for N_ in range(1, 5):
    S_ = np.array([[1]], dtype=complex)
    for _ in range(N_):
        S_ = np.kron(S_, S1)
    bits = list(product((0, 1), repeat=N_))
    other_orientation_fails = False
    for a in bits:
        for b in bits:
            C = cell_matrix(a, b)
            img = S_ @ C @ S_.conj().T
            ok_q &= np.array_equal(img, (1j) ** ((sum(a) - sum(b)) % 4) * C)
            other_orientation_fails |= not np.array_equal(img, (1j) ** ((sum(b) - sum(a)) % 4) * C)
    ok_qc &= other_orientation_fails
check("N=1..4: the quarter-turn about Z, S = diag(1, i) on every site, multiplies the cell |a><b| by i^(p - q), "
      "p = popcount a, q = popcount b: diagonal on the cells, a signed permutation of the strings (S9)", ok_q)
check("control: the other orientation, i^(q - p), fails at every N = 1..4", ok_qc)
# its square is the conjugation by Z^N (S^2 = Z): on every cell it is the sign (-1)^(p - q), on every string (-1)^k_Z,
# so d mod 2 is k_Z mod 2, which is bit_a
ok_sq = True
for N_ in range(1, 5):
    S_ = np.array([[1]], dtype=complex)
    for _ in range(N_):
        S_ = np.kron(S_, S1)
    ZN_ = string_matrix("Z" * N_)
    ok_sq &= np.array_equal(S_ @ S_, ZN_)
    for a in product((0, 1), repeat=N_):
        for b in product((0, 1), repeat=N_):
            C = cell_matrix(a, b)
            ok_sq &= np.array_equal(ZN_ @ C @ ZN_, (-1) ** (sum(a) - sum(b)) * C)
    if N_ <= 3:
        for s in product("IXYZ", repeat=N_):
            M_ = string_matrix(s)
            ok_sq &= np.array_equal(ZN_ @ M_ @ ZN_, (-1) ** kvec(s)[0] * M_)
check("N=1..4: (S^N)^2 = Z^N; the conjugation by Z^N is (-1)^(p - q) on every cell and (-1)^k_Z on every string "
      "(N <= 3): d mod 2 is k_Z mod 2, bit_a", ok_sq)

# (d) which processes are diagonal on which basis, N = 3 (S9's strs, pauli_super, site_op, Id; row-major vec,
#     rho -> A rho B is kron(A, B^T)); integer rates and couplings, exact comparisons
assert N == 3 and d == 8
cells3 = [(a, b) for a in range(d) for b in range(d)]       # row-major: index a * d + b


def conj_super(P):                                          # rho -> P rho P, P a Hermitian string
    return np.kron(P, P.conj())


def ham_super(h):                                           # rho -> -i [h, rho]
    return -1j * (np.kron(h, Id) - np.kron(Id, h.T))


COV = np.array([[0, 2, -1], [2, 0, 3], [-1, 3, 0]])         # covariances C_kl of correlated Z noise


def corr_super(C=COV):                                      # sum_{k != l} C_kl (Z_k rho Z_l - {Z_l Z_k, rho}/2)
    out = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        for l in range(N):
            if k != l:
                zk, zl = site_op("Z", k), site_op("Z", l)
                out += C[k, l] * (np.kron(zk, zl.T) - 0.5 * np.kron(zl @ zk, Id) - 0.5 * np.kron(Id, (zl @ zk).T))
    return out


def site_matrix(A, l):                                     # a one-site 2x2 matrix on site l of the N = 3 chain
    out = np.array([[1]], dtype=complex)
    for k in range(N):
        out = np.kron(out, A if k == l else np.eye(2))
    return out


one = np.kron(Id, Id)
SM0 = site_matrix(np.array([[0, 1], [0, 0]]), 0)           # sigma^- = |0><1| on site 0 (Z|0> = |0>, as in §9)
DAMP = 2 * (np.kron(SM0, SM0.conj()) - 0.5 * np.kron(SM0.conj().T @ SM0, Id) - 0.5 * np.kron(Id, (SM0.conj().T @ SM0).T))
X0, Z0 = site_op("X", 0), site_op("Z", 0)                 # the jump (X + Z)/sqrt2 at rate 2, written without sqrt2
TILT = (np.kron(X0, X0) + np.kron(X0, Z0) + np.kron(Z0, X0) + np.kron(Z0, Z0)) - 2 * one
PROCESSES = {
    "Z dephasing, rates 1, 2, 3 per site": sum(g_ * (conj_super(site_op("Z", l)) - one) for l, g_ in enumerate((1, 2, 3))),
    "X dephasing, rates 1, 2, 3": sum(g_ * (conj_super(site_op("X", l)) - one) for l, g_ in enumerate((1, 2, 3))),
    "Y dephasing, rates 2, 1, 1": sum(g_ * (conj_super(site_op("Y", l)) - one) for l, g_ in enumerate((2, 1, 1))),
    "correlated Z dephasing, the covariances": corr_super(),
    "Z field 2 Z_0 + Z_2": ham_super(2 * site_op("Z", 0) + site_op("Z", 2)),
    "ZZ bond on sites 0, 1": ham_super(site_op("Z", 0) @ site_op("Z", 1)),
    "XX + YY bond on sites 0, 1": ham_super(site_op("X", 0) @ site_op("X", 1) + site_op("Y", 0) @ site_op("Y", 1)),
    "XX bond alone on sites 1, 2": ham_super(site_op("X", 1) @ site_op("X", 2)),
    "X field on site 1": ham_super(site_op("X", 1)),
    "amplitude damping, rate 2 on site 0": DAMP,
    "dephasing along (X + Z)/sqrt2, rate 2 on site 0": TILT,
}
EXPECT = {   # (diagonal on the strings, diagonal on the cells, keeps (p, q)); each column holds both answers
    "Z dephasing, rates 1, 2, 3 per site": (True, True, True),
    "X dephasing, rates 1, 2, 3": (True, False, False),
    "Y dephasing, rates 2, 1, 1": (True, False, False),
    "correlated Z dephasing, the covariances": (False, True, True),
    "Z field 2 Z_0 + Z_2": (False, True, True),
    "ZZ bond on sites 0, 1": (False, True, True),
    "XX + YY bond on sites 0, 1": (False, False, True),
    "XX bond alone on sites 1, 2": (False, False, False),
    "X field on site 1": (False, False, False),
    "amplitude damping, rate 2 on site 0": (False, False, False),
    "dephasing along (X + Z)/sqrt2, rate 2 on site 0": (False, False, False),
}


def is_diagonal(M):
    return np.array_equal(M, np.diag(np.diag(M)))


def keeps_charges(M):
    rows, cols = np.nonzero(M)
    return all((pc(cells3[i][0]), pc(cells3[i][1])) == (pc(cells3[j][0]), pc(cells3[j][1])) for i, j in zip(rows, cols))


got = {name: (is_diagonal(pauli_super(M)), is_diagonal(M), keeps_charges(M)) for name, M in PROCESSES.items()}
# every column of EXPECT holds both answers, so a reader of either basis that answered the same for every process
# would turn some row red; the rows are the control
for name in PROCESSES:
    check(f"N=3, {name}: diagonal on the strings {got[name][0]}, on the cells {got[name][1]}, keeps (p, q) "
          f"{got[name][2]}", got[name] == EXPECT[name])

# (e) correlated Z noise on a string: -4 C_kl Z_k Z_l sigma for every pair k < l of lit sites, nothing else
ok_c, one_lit_fails = True, False
for s in strs:
    sig = string_matrix(s)
    v = (PROCESSES["correlated Z dephasing, the covariances"] @ sig.reshape(-1)).reshape(d, d)
    lit = [ch in "XY" for ch in s]
    pairs = [(k, l) for k in range(N) for l in range(k + 1, N)]
    want = sum((-4 * COV[k, l] * site_op("Z", k) @ site_op("Z", l) @ sig for k, l in pairs if lit[k] and lit[l]),
               np.zeros((d, d), dtype=complex))
    loose = sum((-4 * COV[k, l] * site_op("Z", k) @ site_op("Z", l) @ sig for k, l in pairs if lit[k] or lit[l]),
                np.zeros((d, d), dtype=complex))
    ok_c &= np.array_equal(v, want)
    one_lit_fails |= not np.array_equal(v, loose)
check("N=3: correlated Z noise sends a string to -4 C_kl Z_k Z_l sigma for each pair of its lit sites and to "
      "nothing else: on the strings a covariance moves letters at fixed k_Z, it is no height", ok_c)
check("control: counting every pair with a lit site, not only the pairs lit at both, fails", one_lit_fails)
# collective dephasing (one jump sum_k Z_k at rate gamma): every C_kl = gamma; at gamma = 1 the largest off-diagonal
# entry on the strings is 4 = 4 C_kl, so the absorption theorem's 0.2 at gamma = 0.05 is this 4 gamma
collective = pauli_super(corr_super(np.ones((N, N)) - np.eye(N)))
off = collective - np.diag(np.diag(collective))
check("N=3, collective dephasing at gamma = 1 (all C_kl = 1): the largest off-diagonal entry on the strings is 4, "
      "and every nonzero one is +-4", np.max(np.abs(off)) == 4 and set(np.abs(off[off != 0])) == {4.0})

# (f) neither labelling is a function of the other: a string spreads over several blocks (N = 3), and a single cell
#     spreads over every point of its height plane k_Z = k, all (k + 1)(N - k + 1) of them (N = 1..4)
blocks_of = {s: {(pc(a), pc(b)) for a in range(d) for b in range(d) if string_matrix(s)[a, b] != 0} for s in strs}
ok_plane = True
for N_ in range(1, 5):
    strs_ = ["".join(t) for t in product("IXYZ", repeat=N_)]
    mats_ = {s: string_matrix(s) for s in strs_}
    for a in range(2 ** N_):
        for b in range(2 ** N_):
            k_ = pc(a ^ b)
            pts = {kvec(s) for s in strs_ if mats_[s][a, b] != 0}
            plane = {kvec(s) for s in strs_ if kvec(s)[0] == k_}
            ok_plane &= pts == plane and len(plane) == (k_ + 1) * (N_ - k_ + 1)
check("N=3: the string ZII touches the four blocks (p, p), p = 0..3",
      blocks_of["ZII"] == {(0, 0), (1, 1), (2, 2), (3, 3)})
check("N=1..4: every cell |a><b| spreads over every point of its height plane k_Z = popcount(a xor b), all "
      "(k + 1)(N - k + 1) of them", ok_plane)

# (h) amplitude damping keeps d and splits differently in the two bases (N = 3, site 0, rate 2): on the cells its only
#     off-diagonal entries are the jump |1><1| -> |0><0| at site 0 (p and q down by one); on the strings they are the
#     arrow I -> Z at site 0; the two off-diagonal parts are different operators
cells_off = DAMP - np.diag(np.diag(DAMP))
ok_cell_jump = len(np.nonzero(cells_off)[0]) > 0
for i, j in zip(*np.nonzero(cells_off)):
    (a1, b1), (a0, b0) = cells3[i], cells3[j]
    top = 1 << (N - 1)                                      # site 0 is the top bit
    ok_cell_jump &= bool(a0 & top) and bool(b0 & top) and a1 == a0 ^ top and b1 == b0 ^ top
P_damp = pauli_super(DAMP)
strings_off = P_damp - np.diag(np.diag(P_damp))
ok_arrow = len(np.nonzero(strings_off)[0]) > 0
for i, j in zip(*np.nonzero(strings_off)):                  # the arrow I -> Z with coefficient +rate, sign included
    ok_arrow &= strs[j][0] == "I" and strs[i] == "Z" + strs[j][1:] and strings_off[i, j] == 2
same_operator = np.array_equal(cells_off, B0 @ strings_off @ B0.conj().T / d)
check("N=3, amplitude damping on site 0: on the cells the only off-diagonal entries are the jump |1><1| -> |0><0| "
      "there, lowering p and q by one", ok_cell_jump)
check("N=3, amplitude damping on site 0: on the strings the only off-diagonal entries are the arrow I -> Z there, "
      "coefficient +rate", ok_arrow)
check("N=3: the two off-diagonal parts are different operators (the jump is not the arrow)", not same_operator)

# (g) F98's Krawtchouk expansion is the transform at x = 0; at x = 1...1 the same expansion carries X^N
ok_k, ok_kc = True, True
for N_ in range(1, 5):
    bits = list(product((0, 1), repeat=N_))
    XN_ = string_matrix("X" * N_)
    unsigned_fails = False
    for n_ in range(N_ + 1):
        P = sum(cell_matrix(b, b) for b in bits if sum(b) == n_)
        acc, acc_unsigned = np.zeros_like(P), np.zeros_like(P)
        for S in bits:                                       # S as an indicator tuple
            k_ = sum(S)
            # F98's K_|S|(n; N) = sum over b with popcount n of (-1)^|S & b|: j of the |S| sites in S set, n - j outside
            K = sum((-1) ** j * comb(k_, j) * comb(N_ - k_, n_ - j) for j in range(min(k_, n_) + 1))
            K_unsigned = sum(comb(k_, j) * comb(N_ - k_, n_ - j) for j in range(min(k_, n_) + 1))
            zs = "".join("Z" if t else "I" for t in S)
            ys = "".join("Y" if t else "X" for t in S)
            ZS = string_matrix(zs)
            acc = acc + K * ZS
            acc_unsigned = acc_unsigned + K_unsigned * ZS
            ok_k &= kvec(zs) == (0, k_, k_) and kvec(ys) == (N_, k_, N_ - k_)
            ok_k &= np.array_equal(XN_ @ ZS, (-1j) ** (k_ % 4) * string_matrix(ys))    # X Z = -i Y on each site of S
        ok_k &= np.array_equal(2 ** N_ * P, acc)
        unsigned_fails |= not np.array_equal(2 ** N_ * P, acc_unsigned)
    ok_kc &= unsigned_fails
check("N=1..4: P_n = 2^-N sum_S K Z_S with K = sum_j (-1)^j C(s, j) C(N - s, n - j) at s = |S| (F88b's and "
      "PopcountCoherencePi2Odd's K_n(s; N), which is F98's K_|S|(n; N)) is the transform at x = 0, Z_S at "
      "(0, |S|, |S|), so |S| is k_X = k_Y; X^N Z_S is (-i)^|S| times the string with Y on S and X elsewhere, at "
      "(N, |S|, N - |S|)", ok_k)
check("control: the Krawtchouk number without its sign (-1)^j fails at every N = 1..4", ok_kc)

# (i) the doubled site of PROOF_FROZEN_BAND_SO4, the ket as one species and the bra as the other: on one site the cells
#     |0><0|, |1><1| are the eta doublet (empty, doubly occupied) and |1><0|, |0><1| the spin doublet (up, down); the maps
#     rho -> s+ rho s- and rho -> s+ rho s+ raise them (the two ladders Phi and S+ of that proof at one site); the cells are
#     the eigenbasis of both z-components and the letters that of both x-components
sp, sm, n1 = np.array([[0, 0], [1, 0]]), np.array([[0, 1], [0, 0]]), np.array([[0, 0], [0, 1]])
I2_ = np.eye(2)


def sup(A, B):                                              # rho -> A rho B on one site, row-major vec
    return np.kron(A, B.T)


def comm(A, B):
    return A @ B - B @ A


eta_p, eta_m = sup(sp, sm), sup(sm, sp)
spin_p, spin_m = sup(sp, sp), sup(sm, sm)
eta_z = 0.5 * (sup(n1, I2_) + sup(I2_, n1)) - 0.5 * np.eye(4)
spin_z = 0.5 * (sup(n1, I2_) - sup(I2_, n1))
eta_x, spin_x = (eta_p + eta_m) / 2, (spin_p + spin_m) / 2
ok_su2 = (np.array_equal(comm(eta_z, eta_p), eta_p) and np.array_equal(comm(eta_p, eta_m), 2 * eta_z)
          and np.array_equal(comm(spin_z, spin_p), spin_p) and np.array_equal(comm(spin_p, spin_m), 2 * spin_z)
          and all(not comm(A, B).any() for A in (eta_p, eta_m, eta_z) for B in (spin_p, spin_m, spin_z)))
check("one site: rho -> s+ rho s- and rho -> s+ rho s+, their adjoints, and the Cartans eta^z = (n_ket + n_bra - 1)/2 "
      "and S^z = (n_ket - n_bra)/2 form two commuting su(2), eta on |0><0|, |1><1| and spin on |1><0|, |0><1|", ok_su2)


def vec(M):
    return np.asarray(M, dtype=complex).reshape(-1)


letters_ok = (np.array_equal(eta_x @ vec(LET["I"]), 0.5 * vec(LET["I"]))
              and np.array_equal(eta_x @ vec(LET["Z"]), -0.5 * vec(LET["Z"]))
              and np.array_equal(spin_x @ vec(LET["X"]), 0.5 * vec(LET["X"]))
              and np.array_equal(spin_x @ vec(LET["Y"]), -0.5 * vec(LET["Y"]))
              and not (spin_x @ vec(LET["I"])).any() and not (spin_x @ vec(LET["Z"])).any()
              and not (eta_x @ vec(LET["X"])).any() and not (eta_x @ vec(LET["Y"])).any())
check("one site: I and Z are the eta^x = +1/2, -1/2 eigenvectors and X and Y the S^x = +1/2, -1/2 eigenvectors, each "
      "pair annihilated by the other su(2)'s x-component", letters_ok)
ok_axes = (np.array_equal(sup(LET["X"], LET["X"]), 2 * (eta_x + spin_x))
           and np.array_equal(sup(LET["Y"], LET["Y"]), 2 * (eta_x - spin_x))
           and np.array_equal(sup(LET["Z"], LET["Z"]), np.eye(4) - 8 * spin_z @ spin_z))
check("one site: X rho X = 2(eta^x + S^x), Y rho Y = 2(eta^x - S^x), Z rho Z = 1 - 8 (S^z)^2 (-1 on the spin doublet); "
      "the sums over the sites, Q_X, Q_Y and Q_Z = N - 2 (singly occupied sites), follow site by site", ok_axes)
two_s = [int(t) for t in 2 * np.real(np.diag(spin_z))]
check("one site: the quarter-turn about Z, rho -> S rho S^dagger, is exp(i pi S^z), a half-turn of the spin",
      np.array_equal(sup(S1, S1.conj().T), np.diag([(1j) ** (t % 4) for t in two_s])))
R1 = sup(I2_, LET["X"])                                     # rho -> rho X, its own inverse
D1 = np.zeros((4, 4))
for a_ in range(2):
    for b_ in range(2):
        D1[2 * b_ + a_, 2 * a_ + b_] = 1                    # rho -> rho^T
P_eta = np.diag([1.0, 0.0, 0.0, 1.0])                       # the cells |0><0|, |1><1|
check("one site: R exchanges the doublets, R eta^z R = S^z and R S^z R = eta^z, which is section 3's (d, m) -> (m, d); the "
      "transpose D is 1 on the eta doublet and 2 S^x on the spin doublet, (-1)^n_Y on the strings",
      np.array_equal(R1 @ eta_z @ R1, spin_z) and np.array_equal(R1 @ spin_z @ R1, eta_z)
      and np.array_equal(D1, P_eta + 2 * spin_x))


def spin_z_at(l):                                           # S^z of site l on the N = 3 chain, as a superoperator
    nl = site_matrix(n1, l)
    return 0.5 * (np.kron(nl, Id) - np.kron(Id, nl.T))


a_full = np.eye(N) + COV                                   # a = gamma I + C at gamma = 1
full_corr = sum(conj_super(site_op("Z", l)) - one for l in range(N)) + corr_super()
ising = -8 * sum(a_full[k, l] * spin_z_at(k) @ spin_z_at(l) for k in range(N) for l in range(N))
local_part = -8 * sum(spin_z_at(l) @ spin_z_at(l) for l in range(N))
check("N=3: correlated Z noise with a = gamma I + C is -8 sum_{k,l} a_kl S^z_k S^z_l, the shared field's form: the k = l "
      "terms are the local price -2 gamma on each lit site, the k != l terms the Ising exchange",
      np.array_equal(full_corr, ising) and np.array_equal(sum(conj_super(site_op("Z", l)) - one for l in range(N)),
                                                          local_part))

# (j) the quarter-turn's Z4 on the cells: multiplicities sum_{d = r mod 4} C(2N, N + d), the same multiset as the
#     Majorana-degree count 4^(N-1) + 2^(N-1) Re(i^(N-r)) of ONE_FOUR_THESIS, never Pi's uniform 4^(N-1); the two
#     gradings differ (at one site Z has d = 0 but Majorana degree 2) and agree mod 2
ok_mult = True
for N_ in range(1, 9):
    cells_mult = sorted(sum(comb(2 * N_, N_ + dd) for dd in range(-N_, N_ + 1) if (dd - r) % 4 == 0) for r in range(4))
    a4, b2 = 4 ** (N_ - 1), 2 ** (N_ - 1)
    thesis = sorted(a4 + b2 * [1, 0, -1, 0][(N_ - r) % 4] for r in range(4))      # Re(i^k), k = 0..3
    ok_mult &= cells_mult == thesis and cells_mult != [a4] * 4
    if N_ <= 4:                                                                  # counted from the cells themselves
        direct = sorted(sum(1 for a in range(2 ** N_) for b in range(2 ** N_) if (pc(a) - pc(b) - r) % 4 == 0)
                        for r in range(4))
        ok_mult &= direct == cells_mult
check("N=1..8: the quarter-turn's Z4 has the multiplicities sum_{d = r mod 4} C(2N, N + d) on the cells (counted "
      "directly for N <= 4), the multiset of the Majorana-degree count 4^(N-1) + 2^(N-1) Re(i^(N-r)), never Pi's "
      "uniform 4^(N-1)", ok_mult)
majorana_deg = {"I": 0, "X": 1, "Y": 1, "Z": 2}                                  # one site
check("one site: the two gradings differ, Z is fixed by the quarter-turn (d = 0) and has Majorana degree 2, and they "
      "agree mod 2 on every letter (k_Z = Majorana degree mod 2)",
      np.array_equal(S1 @ LET["Z"] @ S1.conj().T, LET["Z"]) and majorana_deg["Z"] % 4 != 0
      and all(coord[L][0] == majorana_deg[L] % 2 for L in "IXYZ"))

print()
if FAILURES:
    print(f"diagonal square gate: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("diagonal square gate: ALL GREEN")
