"""F50 at exceptional couplings: where the count d_real(Re = -2γ) = 2N (+ δ_G) is exceeded.

H = J·Σ_bonds (XX + YY + ZZ) (Pauli convention, J = 1 here), uniform Z-dephasing γ, L acting on
|x⟩⟨y| with the dephasing diagonal -2γ·Hamming(x, y). L is block diagonal in the joint popcount
(p, q) = (popcount x, popcount y).

The theorem this gate checks (docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, "The count at exceptional
couplings"):
  (a) By the Absorption Theorem, Re λ = -2γ⟨n_XY⟩ for every eigenvector, so a real eigenvalue at
      -2γ needs ⟨n_XY⟩ = 1. In a block (p, q) every coherence has Hamming distance ≥ |p - q| and of
      the parity of p - q. Odd p - q: weight ≥ 1 everywhere, so the mode is pure weight 1 (the
      commutant F50 counts). Even p ≠ q: weight ≥ 2, impossible. Only the diagonal blocks (p, p),
      1 ≤ p ≤ N-1, can hold a real -2γ mode that is not pure weight 1.
  (b) In (p, p) the Hamming distances are even, never 1, so P_p(γ) = det(B_pp(γ) + 2γ·I) has
      leading coefficient Π (2 - 2·Hamming) ≠ 0: a nonzero polynomial, finitely many roots. The
      exceptional set E(N, G) is the set of its positive real roots over all p.
  (c) Off E the converse holds: every real eigenvalue at -2γ is a weight-1 commutant mode.

G1  the leading coefficient of every P_p is the product predicted in (b) (exact).
G2  the off-diagonal even blocks (p, q), p ≠ q, have no positive root of det(B + 2γ) (exact; a
    consistency check of (a), whose argument needs no computation).
G3  E(N, G) for N = 2, 3, 4 on chain, ring, star, complete, pinned to exact minimal polynomials.
G4  at N = 2, 3 the multiplicity of -2γ at each exceptional point, exactly (rank over the
    algebraic numbers): algebraic and geometric, and the weight content of the extra kernel.
G5  (measured) at N = 4 the eigensolver count at -2γ at each exceptional point and at a generic
    point, window 1e-6 on |λ + 2γ| (an EP2 splits by ~sqrt(eps)·‖L‖ ≈ 1e-7).
G6  the chain's smallest exceptional γ/J is 1/Q*_gap(N), the Heisenberg spectral-gap threshold that
    simulations/absorption_ladder_regimes.py bisects to six decimals (0.500000, 0.800243,
    1.342243, and 1.819350 under --n5): agreement within the printed precision, 5e-7 in Q.
Run with --n5 to add the N = 5 chain (about twenty minutes).
"""
import sys
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ_I

g = sp.symbols('g')
R = ZZ_I[g]
FAIL = []
sys.stdout.reconfigure(encoding="utf-8")


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def graph(N, kind):
    if kind == 'chain': return [(i, i + 1) for i in range(N - 1)]
    if kind == 'ring': return [(i, (i + 1) % N) for i in range(N)]
    if kind == 'star': return [(0, i) for i in range(1, N)]
    if kind == 'complete': return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(kind)


def ham_matrix(N, bonds):
    """XX + YY + ZZ = 2·SWAP - I per bond, integer entries."""
    d = 2 ** N
    H = [[0] * d for _ in range(d)]
    bit = lambda x, l: (x >> (N - 1 - l)) & 1
    for a, b in bonds:
        for x in range(d):
            H[x][x] -= 1
            y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b)) if bit(x, a) != bit(x, b) else x
            H[y][x] += 2
    return H


def block(N, H, p, q):
    d = 2 ** N
    pc = lambda x: bin(x).count('1')
    idx = [(x, y) for x in range(d) for y in range(d) if pc(x) == p and pc(y) == q]
    pos = {e: i for i, e in enumerate(idx)}
    return idx, pos


def shifted_block_poly(N, H, p, q):
    """det(B_pq(g) + 2g·I) in ZZ_I[g], and the Hamming distances on the block's diagonal."""
    d = 2 ** N
    idx, pos = block(N, H, p, q)
    n = len(idx)
    rows = [[R.zero] * n for _ in range(n)]
    I = R.convert(sp.I); G = R.convert(g)
    hams = []
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if H[xp][x]: rows[pos[(xp, y)]][c] += -I * R.convert(H[xp][x])
        for yp in range(d):
            if H[y][yp]: rows[pos[(x, yp)]][c] += I * R.convert(H[y][yp])
        h = bin(x ^ y).count('1'); hams.append(h)
        rows[c][c] += G * R.convert(2 - 2 * h)
    return n, sp.expand(R.to_sympy(DomainMatrix(rows, (n, n), R).det())), hams


def positive_roots(P):
    cs = sp.Poly(P, g).all_coeffs()
    A = sp.Poly([sp.re(c) for c in cs], g); B = sp.Poly([sp.im(c) for c in cs], g)
    Q = A if B.is_zero else (B if A.is_zero else sp.gcd(A, B))
    if Q.degree() <= 0: return {}
    out = {}
    for r in Q.real_roots():
        if r > 0: out[r] = out.get(r, 0) + 1
    return out


def exact_block(N, H, p, gval):
    """B_pp at an exact algebraic g as a sympy Matrix (for exact ranks at N <= 3)."""
    d = 2 ** N
    idx, pos = block(N, H, p, p)
    n = len(idx)
    M = sp.zeros(n, n)
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if H[xp][x]: M[pos[(xp, y)], c] += -sp.I * H[xp][x]
        for yp in range(d):
            if H[y][yp]: M[pos[(x, yp)], c] += sp.I * H[y][yp]
        M[c, c] += -2 * gval * bin(x ^ y).count('1')
    return M, idx


def full_L(N, bonds, gam):
    X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]); I2 = np.eye(2)
    def op(l, P):
        r = np.eye(1)
        for j in range(N): r = np.kron(r, P if j == l else I2)
        return r
    d = 2 ** N
    H = sum(op(a, P) @ op(b, P) for a, b in bonds for P in (X, Y, Z))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = op(l, Z); L += gam * (np.kron(Zl, Zl.T) - np.eye(d * d))
    return L


# E(N, G) pinned: (N, kind) -> {p: [(minimal polynomial in x = γ/J, root index among positive roots ascending, multiplicity in g)]}
x = sp.symbols('x')
PINNED = {
    (2, 'chain'):    {1: [(x - 2, 1)]},
    (3, 'chain'):    {1: [(x**4 + x**2 - 4, 1), (x**2 - 3, 1)], 2: [(x**4 + x**2 - 4, 1), (x**2 - 3, 1)]},
    (3, 'complete'): {1: [(x**2 - 3, 2)], 2: [(x**2 - 3, 2)]},
    (4, 'chain'):    {1: [(x**4 + 4*x**2 - 4, 1), (x**4 + 4*x**2 - 16, 1), (x**4 - 4*x**2 - 4, 1)],
                      2: [(9*x**12 + 132*x**10 + 68*x**8 - 1696*x**6 - 2240*x**4 + 1280*x**2 + 256, 1),
                          (3*x**8 + 36*x**6 - 84*x**4 - 656*x**2 + 384, 1),
                          (3*x**6 + 4*x**4 - 20*x**2 - 16, 1),
                          (9*x**12 + 132*x**10 + 68*x**8 - 1696*x**6 - 2240*x**4 + 1280*x**2 + 256, 1),
                          (3*x**8 + 36*x**6 - 84*x**4 - 656*x**2 + 384, 1)],
                      3: [(x**4 + 4*x**2 - 4, 1), (x**4 + 4*x**2 - 16, 1), (x**4 - 4*x**2 - 4, 1)]},
    (4, 'ring'):     {1: [(x - 2, 2)], 2: [(3*x**4 + 28*x**2 - 64, 1), (3*x**2 - 16, 1)], 3: [(x - 2, 2)]},
    (4, 'star'):     {1: [(x**4 + 6*x**2 - 3, 2)],
                      2: [(3*x**2 - 4, 2), (x**2 - 3, 2), (3*x**2 - 16, 1)],
                      3: [(x**4 + 6*x**2 - 3, 2)]},
    (4, 'complete'): {2: [(3*x**2 - 16, 3)]},
}

E = {}
for (N, kind), pinned in PINNED.items():
    H = ham_matrix(N, graph(N, kind))
    found = {}
    for p in range(N + 1):
        for q in range(N + 1):
            if (p - q) % 2 or (p, q) > (q, p): continue
            n, P, hams = shifted_block_poly(N, H, p, q)
            if p == q:
                lead = sp.Poly(P, g).LC()
                pred = sp.prod([2 - 2 * h for h in hams])
                check(f"G1 N={N} {kind} block ({p},{p}): leading coefficient = Π(2 - 2·Hamming) = {pred}, nonzero (a construction check)",
                      pred != 0 and sp.simplify(lead - pred) == 0, f"got {lead}")
                r = positive_roots(P)
                if r: found[p] = r
            else:
                r = positive_roots(P)
                check(f"G2 N={N} {kind} block ({p},{q}): no positive root of det(B + 2γ)", not r, f"{r}")
    # compare with pinned minimal polynomials
    ok = set(found) == set(pinned)
    detail = []
    for p in pinned:
        if p not in found: ok = False; continue
        got = sorted(found[p].items(), key=lambda t: float(t[0]))
        # each found root must be a root of exactly the pinned polynomial at the same position
        if len(got) != len(pinned[p]): ok = False; detail.append(f"p={p} count {len(got)} vs {len(pinned[p])}"); continue
        for (r, m), (mp, mw) in zip(got, pinned[p]):
            if m != mw or sp.Poly(sp.minimal_polynomial(r, x), x).monic() != sp.Poly(mp, x).monic(): ok = False; detail.append(f"p={p} root {float(r):.10f}")
    E[(N, kind)] = found
    if kind == 'chain':
        per = {p: len(found.get(p, {})) for p in range(1, N)}
        check(f"G3b N={N} chain: block (p,p) has C(N,p) - 1 distinct exceptional points: {per}",
              all(per[p] == sp.binomial(N, p) - 1 for p in per))
    pts = sorted({r for p in found for r in found[p]}, key=float)
    check(f"G3 N={N} {kind}: E = {{{', '.join(f'{float(r):.10f}' for r in pts)}}} (γ/J), pinned minimal polynomials",
          ok, '; '.join(detail))

# G4: exact multiplicities at N = 2, 3, ranks over the field Q(i, γ/J)
def to_field(M, gv):
    K = sp.QQ.algebraic_field(sp.I) if gv.is_Rational else sp.QQ.algebraic_field(sp.I, gv)
    rows = [[K.from_sympy(M[i, j]) for j in range(M.shape[1])] for i in range(M.shape[0])]
    return DomainMatrix(rows, M.shape, K), K

def generalized_kernel_dim(A):
    n = A.shape[0]
    P = A
    prev = -1
    for k in range(1, n + 1):
        dim = n - P.rank()
        if dim == prev: return dim
        prev = dim; P = P * A
    return prev

EXPECT_G4 = {  # (N, kind, γ/J) -> (extra algebraic, extra geometric) summed over the diagonal blocks
    (2, 'chain', sp.Integer(2)): (2, 1),
    (3, 'chain', sp.sqrt(3)): (2, 2),
    (3, 'chain', sp.sqrt((sp.sqrt(17) - 1) / 2)): (2, 2),
    (3, 'complete', sp.sqrt(3)): (4, 4),
}
for (N, kind, gv), (ea, eg) in EXPECT_G4.items():
    H = ham_matrix(N, graph(N, kind))
    alg = geo = 0
    weights_ok = True
    for p in range(1, N):
        M, idx = exact_block(N, H, p, gv)
        A, K = to_field(M + 2 * gv * sp.eye(M.shape[0]), gv)
        geo_p = A.shape[0] - A.rank()
        alg_p = generalized_kernel_dim(A)
        geo += geo_p; alg += alg_p
        if geo_p:
            # a vector is in the kernel iff it is a combination of these rows; every kernel vector must mix
            # Hamming 0 with Hamming >= 2, so test that no kernel vector is supported on one weight alone:
            # the kernel meets the span of the Hamming-0 coordinates (and of the >= 2 ones) only in 0
            Ns = A.nullspace().to_Matrix()
            for keep in (lambda h: h == 0, lambda h: h >= 2):
                cols = [k for k, (xx, yy) in enumerate(idx) if not keep(bin(xx ^ yy).count('1'))]
                sub = Ns[:, cols]
                if sub.rank() < Ns.shape[0]: weights_ok = False
    check(f"G4 N={N} {kind} γ/J = {gv}: extra multiplicity at -2γ algebraic {alg}, geometric {geo}; every extra kernel vector mixes Hamming 0 with Hamming ≥ 2",
          (alg, geo) == (ea, eg) and weights_ok, f"expected ({ea}, {eg})")

# G5: measured eigensolver counts at N = 4 (algebraic: eigenvalues within 1e-6 of -2γ; geometric: nullity of
# L + 2γ, singular values below 1e-8·σ_max). The pinned values are the tabulated ones.
G5_COUNTS = {'chain': [9, 9, 10, 9, 10, 9, 9, 10], 'ring': [9, 16, 9], 'star': [12, 10, 10, 9], 'complete': [11]}
G5_NULL = {'chain': [9, 9, 10, 9, 10, 9, 9, 10], 'ring': [9, 12, 9], 'star': [12, 10, 10, 9], 'complete': [11]}
for kind in ('chain', 'ring', 'star', 'complete'):
    bonds = graph(4, kind)
    base = len(np.where(np.abs(np.linalg.eigvals(full_L(4, bonds, 1.3)) + 2.6) < 1e-6)[0])
    pts4 = sorted({r for p in E[(4, kind)] for r in E[(4, kind)][p]}, key=float)
    counts = []
    for r in pts4:
        gam = float(r)
        ev = np.linalg.eigvals(full_L(4, bonds, gam))
        counts.append((round(gam, 6), int(np.sum(np.abs(ev + 2 * gam) < 1e-6))))
    geos = []
    for r in pts4:
        gam = float(r); A = full_L(4, bonds, gam) + 2 * gam * np.eye(256)
        s = np.linalg.svd(A, compute_uv=False)
        geos.append((round(gam, 6), int(np.sum(s < 1e-8 * s[0]))))
    print(f"G5 (measured) N=4 {kind}: generic γ/J = 1.3 count {base}; at E: {counts}; nullity of L + 2γ (singular values below 1e-8·σ_max): {geos}")
    check(f"G5 N=4 {kind}: generic count 8 = 2N, and the counts and nullities at E are the tabulated ones",
          base == 8 and [c for _, c in counts] == G5_COUNTS[kind] and [c for _, c in geos] == G5_NULL[kind])

# G6: the chain's smallest exceptional γ/J against Q*_gap, two routes. (i) The spectrum: just below γ* the gap is
# 2γ, just above it is smaller and set by a real mode (relative offset 1e-6; the eigenvalue error is ~1e-14·‖L‖,
# the crossing moves the slow mode by ~1e-6·γ). (ii) The six decimals absorption_ladder_regimes.py bisects.
def gap_sides(N, gstar, eps=1e-6):
    out = []
    for gam in (gstar * (1 - eps), gstar * (1 + eps)):
        ev = np.linalg.eigvals(full_L(N, graph(N, 'chain'), gam))
        nz = ev[np.abs(ev) > 1e-9]
        k = np.argmin(np.abs(nz.real))
        out.append((gam, -nz[k].real, abs(nz[k].imag)))
    (g1, gap1, _), (g2, gap2, im2) = out
    return abs(gap1 - 2 * g1) < 1e-9 and 2 * g2 - gap2 > 1e-8 and im2 < 1e-9, (gap1 - 2 * g1, 2 * g2 - gap2, im2)

QGAP = {2: 0.500000, 3: 0.800243, 4: 1.342243}
for N, q in QGAP.items():
    xmin = min((r for p in E[(N, 'chain')] for r in E[(N, 'chain')][p]), key=float)
    Qe = 1 / float(xmin)
    ok, d = gap_sides(N, float(xmin))
    check(f"G6 N={N} chain: at γ*(1 ∓ 1e-6), γ* = {float(xmin):.10f} the smallest point of E, the gap is 2γ below and a real mode below 2γ above (offsets {d[0]:.1e}, {d[1]:.1e}, Im {d[2]:.1e}); J/γ* = {Qe:.7f} matches Q*_gap({N}) = {q:.6f} to the printed precision",
          ok and abs(Qe - q) <= 5e-7)

if '--n5' in sys.argv:
    # the N = 5 chain (about twenty minutes): thirteen exceptional points, the smallest at 1/Q*_gap(5)
    H = ham_matrix(5, graph(5, 'chain'))
    pts = set(); per5 = {}; min2 = None
    for p in range(1, 5):
        n, P, hams = shifted_block_poly(5, H, p, p)
        r = positive_roots(P); per5[p] = len(r)
        if p == 2: min2 = min(r, key=float)
        pts |= set(r)
    pts = sorted(pts, key=float)
    print(f"N=5 chain: E = {{{', '.join(f'{float(r):.10f}' for r in pts)}}}")
    check("G3 N=5 chain: thirteen exceptional points", len(pts) == 13)
    check(f"G3b N=5 chain: per block C(5,p) - 1 = 4, 9, 9, 4 points, {per5}", per5 == {1: 4, 2: 9, 3: 9, 4: 4})
    deg = sp.Poly(sp.minimal_polynomial(pts[0], x), x).degree()
    check(f"G3c N=5 chain: the smallest point is algebraic of degree {deg} (48 expected) and comes from the (2, 2) block", deg == 48 and min2 == pts[0])
    Qe = 1 / float(pts[0])
    ok, d = gap_sides(5, float(pts[0]))
    check(f"G6 N=5 chain: gap 2γ below γ*, a real mode below 2γ above (offsets {d[0]:.1e}, {d[1]:.1e}, Im {d[2]:.1e}); J/γ* = {Qe:.7f} matches Q*_gap(5) = 1.819350 to the printed precision",
          ok and abs(Qe - 1.819350) <= 5e-7)

print()
print("FAILED: " + ", ".join(FAIL) if FAIL else "ALL GATES PASS")
sys.exit(1 if FAIL else 0)
