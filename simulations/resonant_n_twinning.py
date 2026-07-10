"""The resonance criterion and the twinning law at the resonant N (F89, 2026-07-10).

Self-asserting. Two independent halves that must agree.

THE COMBINATORIAL HALF (exact, no matrices, no floats). The (N1') theorem gives
nullity(K66) = 3*Z3 with Z3 = #{a<b<c : lam_a + lam_b + lam_c = 0} over the one-magnon
energies lam_k = 2 cos(k pi / n), n = N + 1. Writing lam_k = zeta^k + zeta^-k with zeta a
primitive 2n-th root of unity turns a vanishing triple into a conjugation-symmetric
VANISHING SUM OF SIX 2n-TH ROOTS OF UNITY. Conway & Jones, "Trigonometric diophantine
equations", Acta Arith. 30 (1976) 229-240, classify these: Theorem 6 lists the minimal
vanishing sums of weight <= 9 (the unique weight-6 one is -a-a^2+b+b^2+b^3+b^4, a, b
primitive of order 3, 5), and Theorem 7 the rational-angle cosine relations. Exactly three
families survive, and SINCE N IS ODD, n IS EVEN -- the counts below need that:

  TRIV   three antipodal pairs (R2+R2+R2):  {k, n-k, n/2}          count (N-1)/2
  ROT3   two rotated cube-root triples (R3+R3), iff 3 | n          count n/3 - 2
  PENT   the minimal weight-6 sum (R5:R3), iff 15 | n              count exactly 2

  (For n even, 15 | n iff 30 | n. For ODD n the ROT3 count is 2*floor((n-3)/6) instead;
   that case never arises here, but the formula below must not be quoted outside n even.)

  =>  Z3(N) = (N-1)/2 + [3|n]*(n/3 - 2) + 2*[15|n],  and  RESONANT <=> 3 | (N+1), N >= 11.

The trivial triples are exactly the SELF-MIRROR ones: if {a,b,c} = {n-a,n-b,n-c} then the
involution k -> n-k on a 3-set has a fixed point n/2 and swaps the other two, so the triple
is {a, n-a, n/2}; conversely those vanish by chirality lam_{n-k} = -lam_k. Hence the extras
come in MIRROR PAIRS, and their number is even.

THE MATRIX HALF (float64). T = diag((-1)^(a0+a1+b0)) satisfies T K T = -K, so T preserves
ker(K66) and realizes the mirror k -> n-k. With u_k(z) = sin(k(z+1) pi/n) one has u_{n-k}(z)
= (-1)^z u_k(z), so the Slater determinant obeys v_tau'(z) = -(-1)^(z1+z2+z3) v_tau(z) (the
minus is the reversal permutation of the three modes). A SELF-MIRROR triple therefore has
v_tau supported on ODD site-sum, where T acts by -1: its 3 dims are all class O, which IS
the baseline purity ker(K66) subset O. A mirror PAIR is swapped by T, so it splits 3E + 3O:

  dim_E ker(K66) = 3 * (#extras) / 2,   and   ker(K66) subset O  <=>  N not resonant.

On a mirror pair, ker(K66) is spanned by the gauged Slater determinants of the two mode
triples, which are orthogonal, and T swaps them sector by sector. So Heff := P_ker K62 K26
P_ker restricted to the pair is [[X, Y], [Y, X]] with E-block X+Y and O-block X-Y. Measured:
Y = 0 (to 1e-15). Hence E-block = O-block = X as MATRICES, and the twinning is an
orthogonality, not a coincidence:

  <K26 w_tau, K26 w_tau'> = 0   for mirror partners tau, tau'.

spec(X) = {high, low, 0} with (high, low) = (18/n, 6/n) on a ROT3 pair and (24/n, 8/n) on
the PENT pair. Both ratios are 3:1, both rational; the O-spectrum is otherwise irrational.
Each mirror pair therefore yields 2 coupled E-levels (each exactly twinned) and 1 uncoupled.

WHY IT MATTERS. F89_SEED_EXISTENCE_REDUCTION.md names an untwinned coupled class-E level at
a resonant N as the cheapest falsification of the cell law. This runs that test at the next
resonant N = 23 (and, with '29', at the first N where the PENT family appears). Both hold.

  python simulations/resonant_n_twinning.py        # combinatorics N=3..29 + matrices to N=23, ~15 s
  python simulations/resonant_n_twinning.py 29     # adds N=29 (dim 10962, ~90 s, ~2 GB)
"""
import sys
from itertools import combinations

import numpy as np
from sympy import Poly, cyclotomic_poly, symbols

_x = symbols('x')


# ---------------------------------------------------------------- combinatorics (exact)

def closed_form_z3(N):
    """Z3 by the classification. Valid for n = N+1 EVEN, i.e. N odd."""
    n = N + 1
    assert n % 2 == 0, "the ROT3 count n/3-2 assumes n even"
    triv = (N - 1) // 2
    rot3 = max(0, n // 3 - 2) if n % 3 == 0 else 0
    pent = 2 if n % 15 == 0 else 0
    return triv, rot3, pent


def z3_exact(N):
    """Count vanishing triples exactly in Z[x]/Phi_{2n}(x), classified. No floats."""
    n, m = N + 1, 2 * (N + 1)
    phi = Poly(cyclotomic_poly(m, _x), _x)
    d = phi.degree()
    pw = []
    for e in range(m):
        c = Poly(_x ** e, _x).rem(phi).all_coeffs()
        pw.append([0] * (d - len(c)) + [int(v) for v in c])

    def kind(a, b, c):
        if (n // 2) in (a, b, c) and any(u + v == n for u, v in combinations((a, b, c), 2)):
            return 0                                             # TRIV
        if n % 3 == 0:
            s, ex = 2 * n // 3, {e % m for k in (a, b, c) for e in (k, -k)}
            if any({e, (e + s) % m, (e + 2 * s) % m} <= ex for e in ex):
                return 1                                         # ROT3
        return 2                                                 # PENT

    counts, triples = [0, 0, 0], [[], [], []]
    for a, b, c in combinations(range(1, N + 1), 3):
        s = [0] * d
        for e in (a, m - a, b, m - b, c, m - c):
            v = pw[e]
            for i in range(d):
                s[i] += v[i]
        if not any(s):
            k = kind(a, b, c)
            counts[k] += 1
            triples[k].append((a, b, c))
    return counts, triples


# ---------------------------------------------------------------------- matrices (float)

def _hop(N, k):
    st = list(combinations(range(N), k))
    idx = {s: i for i, s in enumerate(st)}
    H = np.zeros((len(st), len(st)))
    for s in st:
        occ = set(s)
        for i in range(N - 1):
            j = i + 1
            if (i in occ) ^ (j in occ):
                src, dst = (i, j) if i in occ else (j, i)
                H[idx[tuple(sorted((occ - {src}) | {dst}))], idx[s]] += 1.0
    return H


def blocks(N):
    """K66, K26 and the bipartite sign on the -6 rung."""
    H2, H1 = _hop(N, 2), _hop(N, 1)
    K = -(np.kron(H2, np.eye(N)) - np.kron(np.eye(H2.shape[0]), H1))
    basis = [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]
    A = np.array([-2.0 * len(set(a) ^ set(b)) for a, b in basis])
    t = np.array([(-1.0) ** (a[0] + a[1] + b[0]) for a, b in basis])
    m2, m6 = np.isclose(A, -2.0), np.isclose(A, -6.0)
    assert np.allclose(K, K.T) and np.allclose(t[:, None] * K * t[None, :], -K), "T K T = -K"
    return K[np.ix_(m6, m6)], K[np.ix_(m2, m6)], t[m6]


def check_N(N, tol=1e-8):
    (triv, rot3, pent), triples = z3_exact(N)
    assert (triv, rot3, pent) == closed_form_z3(N), f"N={N}: classification != closed form"
    n, Z3 = N + 1, triv + rot3 + pent
    extras = rot3 + pent
    assert extras % 2 == 0, "extras must pair up under the mirror"

    K66, K26, t6 = blocks(N)
    w, V = np.linalg.eigh(K66)
    ker = V[:, np.abs(w) < tol]
    assert ker.shape[1] == 3 * Z3, f"N={N}: nullity(K66)={ker.shape[1]} != 3*Z3={3*Z3}"

    tw, tv = np.linalg.eigh(ker.T @ (t6[:, None] * ker))
    assert np.allclose(np.abs(tw), 1.0, atol=1e-7), "T|ker is not an involution"
    E, O = ker @ tv[:, tw > 0], ker @ tv[:, tw < 0]
    assert E.shape[1] == 3 * extras // 2, f"N={N}: dim_E={E.shape[1]} != 3*extras/2"

    def levels(W):
        if W.shape[1] == 0:
            return np.array([])
        M = K26 @ W
        return np.sort(np.linalg.eigvalsh(M.T @ M))[::-1]

    lE, lO = levels(E), levels(O)
    if E.shape[1] and O.shape[1]:
        off = np.linalg.norm(E.T @ (K26.T @ K26) @ O)
        assert off < 1e-8, f"N={N}: Heff not class-diagonal ({off:.2e})"

    coup = lE[lE > 1e-7]
    assert len(coup) == extras, f"N={N}: {len(coup)} coupled E-levels != {extras} extras"

    # the kill test: every coupled class-E level has an exact class-O twin
    worst = 0.0
    for v in coup:
        g = float(np.min(np.abs(lO - v)))
        worst = max(worst, g)
        assert g < 1e-9 * max(1.0, abs(v)), f"N={N}: UNTWINNED coupled E-level {v}"

    # the closed form for the coupled levels: (18/n, 6/n) per ROT3 pair, (24/n, 8/n) per PENT pair
    want = sorted([18 / n] * (rot3 // 2) + [6 / n] * (rot3 // 2)
                  + [24 / n] * (pent // 2) + [8 / n] * (pent // 2), reverse=True)
    assert np.allclose(coup, want, atol=1e-9), f"N={N}: levels {coup} != closed form {want}"

    print(f"  N={N:>3} n={n:>3} | Z3={Z3:>3} = {triv}T+{rot3}R+{pent}P | "
          f"nullity(K66)={3*Z3:>3} = {E.shape[1]}E+{O.shape[1]}O | "
          f"coupled {len(coup):>2} | worst twin gap {worst:.1e} | "
          f"{'RESONANT' if extras else ''}")
    return triples


def check_Y_is_zero(N, triples, tol=1e-9):
    """The mechanism: <K26 w_tau, K26 w_tau'> = 0 on every mirror pair."""
    n = N + 1
    K66, K26, _ = blocks(N)
    H3 = _hop(N, 3)
    basis = [(a, b) for a in combinations(range(N), 2) for b in combinations(range(N), 1)]
    b6 = [(a, b[0]) for (a, b) in basis if b[0] not in a]
    tri3 = list(combinations(range(N), 3))
    ti = {t: i for i, t in enumerate(tri3)}
    U = np.array([[np.sin(k * (j + 1) * np.pi / n) for j in range(N)] for k in range(1, N + 1)])

    def slater(tau):
        v = np.array([np.linalg.det(np.array([[U[k - 1, z] for z in zz] for k in tau]))
                      for zz in tri3])
        return v / np.linalg.norm(v)

    def lift(tau, s):
        v, w = slater(tau), np.zeros(len(b6))
        for i, (a, b) in enumerate(b6):
            z = tuple(sorted(a + (b,)))
            if list(z).index(b) == s:
                w[i] = ((-1.0) ** b) * v[ti[z]]
        return w / np.linalg.norm(w)

    # a self-mirror (TRIV) triple has its Slater det supported on ODD site-sum, where T = -1:
    # that is why its 3 kernel dims are class O, and the baseline purity ker(K66) subset O.
    even = np.array([i for i, z in enumerate(tri3) if sum(z) % 2 == 0])
    for tau in triples[0]:
        assert np.linalg.norm(slater(tau)[even]) < 1e-9, f"N={N}: self-mirror {tau} has even support"

    # the lifts really are kernel vectors of K66, and all 3*Z3 of them span it (Piece 3's gauge
    # identity, machine-checked here rather than assumed)
    allw = np.array([lift(t, s) for t in sum(triples, []) for s in range(3)]).T
    assert np.linalg.norm(K66 @ allw) < 1e-8, f"N={N}: lifts are not in ker(K66)"
    assert np.linalg.matrix_rank(allw, tol=1e-8) == 3 * len(sum(triples, [])), \
        f"N={N}: lifts do not span ker(K66)"

    extras = triples[1] + triples[2]
    seen, worst_Y, worst_ov = set(), 0.0, 0.0
    for tau in extras:
        taup = tuple(sorted(n - k for k in tau))
        if tau in seen:
            continue
        seen |= {tau, taup}
        assert np.linalg.norm(H3 @ slater(tau)) < 1e-9, f"slater{tau} not in ker(H3)"
        W = np.array([lift(tau, s) for s in range(3)]).T
        Wp = np.array([lift(taup, s) for s in range(3)]).T
        worst_ov = max(worst_ov, float(np.linalg.norm(W.T @ Wp)))
        Y = (K26 @ W).T @ (K26 @ Wp)
        worst_Y = max(worst_Y, float(np.linalg.norm(Y)))
        X = (K26 @ W).T @ (K26 @ W)
        eE = np.sort(np.linalg.eigvalsh(X + Y))[::-1]
        eO = np.sort(np.linalg.eigvalsh(X - Y))[::-1]
        assert np.allclose(eE, eO, atol=1e-9), f"N={N} pair {tau}: E/O blocks not isospectral"
    assert worst_Y < tol, f"N={N}: cross block Y != 0 ({worst_Y:.2e})"
    print(f"  N={N:>3} | mirror pairs {len(seen)//2} | max <W,W'> {worst_ov:.1e} | "
          f"max ||Y|| {worst_Y:.1e}  -> twinning IS an orthogonality")


def main():
    extra = [int(v) for v in sys.argv[1:]]

    print("The resonance criterion (exact, Z[x]/Phi_2n): TRIV + ROT3 + PENT")
    for N in range(3, 30, 2):
        (t, r, p), _ = z3_exact(N)
        assert (t, r, p) == closed_form_z3(N), f"N={N}: closed form fails"
    print("  N = 3..29: classification == closed form Z3 = (N-1)/2 + [3|n](n/3-2) + 2[15|n]  OK")
    res = [N for N in range(3, 30, 2) if sum(closed_form_z3(N)[1:]) > 0]
    assert res == [11, 17, 23, 29], f"resonant set {res} != [11, 17, 23, 29]"
    print(f"  resonant N <= 29: {res}   (= 3 | N+1 and N >= 11; next after 17 is 23)")

    print("\nThe matrix half: nullity, class split, coupled levels, the twin gap")
    tri = {}
    for N in [5, 7, 9, 11, 13, 17, 23] + [n for n in extra if n not in (5, 7, 9, 11, 13, 17, 23)]:
        tri[N] = check_N(N)

    print("\nThe mechanism: the cross block Y on each mirror pair")
    for N in [11, 17, 23] + [n for n in extra if n >= 29]:
        check_Y_is_zero(N, tri[N])

    print("\nAll assertions passed. Twinning HOLDS at every resonant N tested; "
          "the cell law survives its cheapest kill.")


if __name__ == "__main__":
    main()
