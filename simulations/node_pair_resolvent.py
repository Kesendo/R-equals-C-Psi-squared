"""The reduced resolvent between two nodes, and what one detuned bond does to a blind mode.

A seat watched by Z-dephasing is blind to the modes with a node there (F157). Such a mode
carries no light and does not decay. Detune one bond: the mode rotates, leans onto the seat,
and picks up light of order epsilon^2. This file gates the theorem behind that, the criterion
it implies, and the structure of the resulting bond profile.

Setting: open chain of N = 2m+1 sites, uniform hopping J, bond b (between sites b and b+1)
scaled to J(1+eps), Z-dephasing at rate gamma on the CENTRE seat c = m only, single-excitation
sector and within it the (1,1) joint-popcount block. At eps = 0,
psi_k(x) = sqrt(2/(N+1)) sin(pi k (x+1)/(N+1)) and E_k = 2J cos(pi k/(N+1)) (F2b); the modes
blind to the centre are exactly the even k, and there are m of them.

WHAT IS GATED, AND WHAT IS ONLY READ.

  GATED, and able to fail:
    G1   Theorem 1: <x|R_k|y> = 0 whenever psi_k has a node at BOTH x and y.  Exact, by
         minimal polynomial over Q, at the centre seat and off it.
    G1b  anti-vacuity: the same sum on node/non-node pairs, where it must NOT vanish.
    G1c  the zero mode's stronger vanishing: <x|R_k0|y> = 0 for every SAME-PARITY pair,
         node or not, forced by the chiral pairing E_k = -E_{N+1-k}.
    G2   c_k = (1/2)*F65(2k) on the bonds whose profile level is 1/2.  Exact, per mode.
         M1 substitutes psi_k(0)^2 and must go red.
    G3   the profile sum_k c_k(b) = (1/2) min(b+1, m-1-b), exact, every bond; and the
         LEVEL STRUCTURE, that the whole c-vector is a function of the level alone.
         M2 shifts the formula by one and must go red.
    G4   the one-seat light additivity w(u v^dagger) = a + b - 2ab, exactly, on rationals.
         M3 drops the overlap term.  NOTE: this is an algebraic identity about the charged
         cell set, not a physical measurement; it cannot fail for a physical reason.
    G5   the peripheral dimension is m^2 + 1, by an EXACT rank over GF(p) at two primes with
         i a square root of -1, no eigensolver.  M4 moves the charged-cell set to another
         seat and requires the dimension to CHANGE (a mutation of the object, not of the
         prediction).
    G6   the one-node criterion, by the EXACT route: F157's blind count of the DETUNED chain
         by Krylov rank over Q at two rational eps, against the node prediction, at the
         CENTRE seat and every non-incident bond.  Replaces a float threshold with the
         repo's own instrument.
    G7   control, and a fence on this file's own scope: a bond INCIDENT on the watched seat
         is struck away with the seat, so both principal blocks are eps-free and the blind
         count is unchanged at EVERY eps, not merely to second order.  Owned and called
         trivial by PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA; gated here so this file cannot
         re-sell it.
    G9   Corollary B at ALL orders: a node at an end of the moved bond makes the left
         block's characteristic polynomial vanish at E_k IDENTICALLY in eps, so the mode
         keeps its energy and its node at every eps, not merely to second order.  The
         control requires the no-node cells NOT to be identically zero.
    G8   the rate statement, correctly scoped: 2 gamma (c_i + c_j) is the DIAGONAL of the
         second-order effective operator in the dyad basis.  It is the rate on every
         omega != 0 peripheral block (gated against the generator), and it is NOT on the
         omega = 0 block, which always carries vec(I) at rate exactly 0 because z^2 = I.
         The omega = 0 mismatch at N = 9 is gated as a POSITIVE fact, so a version of this
         file that "fixed" the mismatch would go red.

  WHAT THE MUTATIONS ACTUALLY MUTATE.  M4 mutates the OBJECT (the charged-cell set moves
  to another seat and the dimension collapses to 1).  M1, M2 and M3 mutate the PREDICTION
  against a fixed measurement; each shows the comparison can return the other answer, which
  is weaker, and M3 is entailed by G4 (given w = a+b-2ab, w != a+b is exactly ab != 0).
  G3's level-structure cells that compare a bond with its whole-chain mirror image cannot
  fail, the reflection commuting with h; at N = 5 every class is such a pair, so that N
  carries no level-structure content.

  READ, not gated:
    R1   the cubic germ's coefficient kappa_N, read rather than assumed.
    R2   the two regimes in (gamma, eps, J, N): the node law's coefficient is J-free and
         gamma-free; beyond it the rate saturates, and the saturated value and the crossing
         are read with their J and N dependence.

Run:  python simulations/node_pair_resolvent.py
"""

from fractions import Fraction
import sympy as sp
import numpy as np

EPSM = np.finfo(float).eps
NOISE_FACTOR = 100.0     # a rate below 100*eps_machine*||A|| is unresolved, not zero
PRIMES = (2013265921, 1004535809)          # both 1 mod 4, so sqrt(-1) exists
_T = sp.Symbol("t")


# ------------------------------------------------------------------ exact comparisons

def exact_zero(expr):
    """True iff expr is exactly 0, certified by its minimal polynomial over Q.

    sympy.simplify does not close on the nested radicals from N = 13 up; the minimal
    polynomial is an exact certificate and is fast.
    """
    return sp.minimal_polynomial(expr, _T) == _T


def exact_rational(expr):
    """expr as a Fraction if it is rational, else None.  Exact, no float."""
    mp = sp.Poly(sp.minimal_polynomial(expr, _T), _T)
    if mp.degree() != 1:
        return None
    a, b = mp.all_coeffs()
    return Fraction(int(-b), int(a))


# ------------------------------------------------------------------------- the system

def hopping(N, b, eps, J=2.0):
    """Open chain, hopping J, bond b scaled to J(1+eps)."""
    H = np.zeros((N, N))
    t = [J] * (N - 1)
    t[b] = J * (1.0 + eps)
    for j in range(N - 1):
        H[j, j + 1] = H[j + 1, j] = t[j]
    return H


def A_generator(N, c, b, eps, gamma, J=2.0):
    """Row-stack vec of L_A(X) = -i[h,X] + gamma (z X z - X), z = I - 2|c><c|."""
    h = hopping(N, b, eps, J)
    I = np.eye(N)
    z = np.eye(N)
    z[c, c] = -1.0
    return (-1j * (np.kron(h, I) - np.kron(I, h.T))
            + gamma * (np.kron(z, z.T) - np.eye(N * N)))


def sym_modes(N):
    """Exact sine modes and energies (hopping 2, so E_k = 4 cos(pi k/(N+1)))."""
    s = sp.pi / (N + 1)
    amp = sp.sqrt(sp.Rational(2, N + 1))
    return (lambda k, x: amp * sp.sin(s * k * (x + 1))), (lambda k: 4 * sp.cos(s * k))


def blind_modes(N, seat):
    """Mode indices with a node at `seat`: exactly the k with (N+1) | k(seat+1)."""
    return [k for k in range(1, N + 1) if (k * (seat + 1)) % (N + 1) == 0]


def amplitude_exact(N, seat, b, k):
    """<seat| R_k V_b |psi_k>, the first-order seat amplitude of a blind mode."""
    psi, E = sym_modes(N)
    blind = set(blind_modes(N, seat))
    total = sum(psi(j, seat) * 2 * (psi(j, b) * psi(k, b + 1) + psi(j, b + 1) * psi(k, b))
                / (E(k) - E(j))
                for j in range(1, N + 1) if j not in blind)
    return sp.simplify(total)


def resolvent_entry(N, seat, k, x, y):
    """<x| R_k |y>, the reduced resolvent at E_k, exact."""
    psi, E = sym_modes(N)
    return sp.simplify(sum(psi(j, x) * psi(j, y) / (E(k) - E(j))
                           for j in range(1, N + 1) if j != k))


def level(N, b):
    """The predicted profile level of bond b, as an exact Fraction."""
    m = (N - 1) // 2
    bb = b if b < m else (N - 2 - b)     # mirror the RIGHT half onto the left
    return Fraction(1, 2) * min(bb + 1, m - 1 - bb)


# ---------------------------------------------------------------- exact rank machinery

def sqrt_minus_one(p):
    for a in range(2, 200):
        r = pow(a, (p - 1) // 4, p)
        if (r * r) % p == p - 1:
            return r
    raise RuntimeError("no square root of -1")


def rref_rank(rows, ncols, p):
    mat = [r[:] for r in rows]
    rank, pivots = 0, []
    for col in range(ncols):
        piv = next((r for r in range(rank, len(mat)) if mat[r][col] % p), None)
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        inv = pow(mat[rank][col], p - 2, p)
        mat[rank] = [(v * inv) % p for v in mat[rank]]
        for r in range(len(mat)):
            if r != rank and mat[r][col] % p:
                f = mat[r][col]
                mat[r] = [(a - f * bb) % p for a, bb in zip(mat[r], mat[rank])]
        pivots.append(col)
        rank += 1
        if rank == len(mat):
            break
    return rank, pivots, mat[:rank]


def nullspace(rows, ncols, p):
    rank, pivots, red = rref_rank(rows, ncols, p)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for f in free:
        vec = [0] * ncols
        vec[f] = 1
        for i, pc in enumerate(pivots):
            vec[pc] = (-red[i][f]) % p
        basis.append(vec)
    return basis


def A_generator_modp(N, c, p, i_unit, gamma=1):
    d = N * N
    h = [[0] * N for _ in range(N)]
    for j in range(N - 1):
        h[j][j + 1] = h[j + 1][j] = 2
    z = [1] * N
    z[c] = -1
    A = [[0] * d for _ in range(d)]
    for a in range(N):
        for bb in range(N):
            row = a * N + bb
            for x in range(N):
                if h[a][x]:
                    A[row][x * N + bb] = (A[row][x * N + bb] - i_unit * h[a][x]) % p
            for y in range(N):
                if h[y][bb]:
                    A[row][a * N + y] = (A[row][a * N + y] + i_unit * h[y][bb]) % p
            A[row][row] = (A[row][row] + gamma * (z[a] * z[bb] - 1)) % p
    return A


def peripheral_dimension_exact(N, c, p, charged_seat=None):
    """dim of the largest A-invariant subspace inside ker(charged-cell projector).

    A mode has Re lambda = 0 iff its charged weight vanishes (the Hermitian part of A is
    -2 gamma times that projector), so this subspace IS the peripheral space.  `charged_seat`
    defaults to c; passing a different seat is the mutation M4, which breaks the object.
    """
    i_unit = sqrt_minus_one(p)
    A = A_generator_modp(N, c, p, i_unit)
    d = N * N
    cs = c if charged_seat is None else charged_seat
    charged = [x for x in range(d) if (x // N == cs) != (x % N == cs)]
    ann = [[1 if x == q else 0 for x in range(d)] for q in charged]
    basis = nullspace(ann, d, p)
    while True:
        if not basis:
            return 0
        ann = nullspace([[bv[x] for x in range(d)] for bv in basis], d, p)
        AB_cols = []
        for bv in basis:
            col = [0] * d
            for x in range(d):
                bx = bv[x]
                if bx:
                    for yq in range(d):
                        if A[yq][x]:
                            col[yq] = (col[yq] + A[yq][x] * bx) % p
            AB_cols.append(col)
        AB = [[sum(w[yq] * col[yq] for yq in range(d)) % p for col in AB_cols] for w in ann]
        keep = (nullspace(AB, len(basis), p) if AB
                else [[1 if t == s else 0 for t in range(len(basis))]
                      for s in range(len(basis))])
        if len(keep) == len(basis):
            return len(basis)
        basis = [[sum(kv[s] * basis[s][x] for s in range(len(basis))) % p for x in range(d)]
                 for kv in keep]


def blind_count_exact(N, seat, b, eps):
    """F157's blind count of the DETUNED chain, by Krylov rank over Q.  Exact, no float.

    blind(seat) = N - rank of the Krylov matrix [e_j, h e_j, h^2 e_j, ...].
    """
    t = [sp.Integer(2)] * (N - 1)
    t[b] = 2 * (1 + eps)
    h = sp.zeros(N, N)
    for j in range(N - 1):
        h[j, j + 1] = h[j + 1, j] = t[j]
    v = sp.zeros(N, 1)
    v[seat] = 1
    cols = []
    for _ in range(N):
        cols.append(v)
        v = h * v
    return N - sp.Matrix.hstack(*cols).rank()


# ------------------------------------------------------------------------------- gates

RESULTS = []


def check(name, ok, detail):
    """detail REPORTS the measurement; it never restates the intent."""
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def gate_G1():
    print("\nG1  Theorem 1: <x|R_k|y> = 0 whenever psi_k has a node at both x and y")
    for N in (7, 9, 11, 13):
        pairs, zeros = 0, 0
        for k in range(1, N + 1):
            nodes = [x for x in range(N) if (k * (x + 1)) % (N + 1) == 0]
            for x in nodes:
                for y in nodes:
                    pairs += 1
                    zeros += exact_zero(resolvent_entry(N, 0, k, x, y))
        check(f"G1 N={N}", pairs > 0 and zeros == pairs,
              f"{zeros}/{pairs} node-node pairs certified 0 by minimal polynomial")
    print("G1b anti-vacuity: one node, one NON-node, where the same sum must not vanish")
    for N in (7, 9, 11):
        nz, tot = 0, 0
        for k in range(1, N + 1):
            nodes = [x for x in range(N) if (k * (x + 1)) % (N + 1) == 0]
            non = [x for x in range(N) if (k * (x + 1)) % (N + 1) != 0]
            for x in nodes:
                for y in non:
                    tot += 1
                    nz += not exact_zero(resolvent_entry(N, 0, k, x, y))
        check(f"G1b N={N}", tot > 0 and nz == tot,
              f"{nz}/{tot} node/non-node pairs are NONZERO")
    print("G1c the zero mode vanishes on every SAME-PARITY pair, node or not")
    for N in (7, 11, 15):
        k0 = (N + 1) // 2
        same, diff_nz, tot_same, tot_diff = 0, 0, 0, 0
        for x in range(N):
            for y in range(N):
                v = resolvent_entry(N, 0, k0, x, y)
                if (x - y) % 2 == 0:
                    tot_same += 1
                    same += exact_zero(v)
                else:
                    tot_diff += 1
                    diff_nz += not exact_zero(v)
        check(f"G1c N={N} k0={k0}", same == tot_same and diff_nz > 0,
              f"{same}/{tot_same} same-parity entries are 0; "
              f"{diff_nz}/{tot_diff} opposite-parity entries are NONZERO")


def gate_G2():
    print("\nG2  c_k = (1/2)*F65(2k) = (2/(N+1)) sin^2(2 pi k/(N+1)) on the level-1/2 bonds")
    for N in (7, 9, 11):
        c = (N - 1) // 2
        psi, _ = sym_modes(N)
        ks = blind_modes(N, c)
        half = [b for b in range(N - 1) if level(N, b) == Fraction(1, 2)]
        hit = all(exact_zero(amplitude_exact(N, c, b, k) ** 2 - psi(k, 1) ** 2)
                  for b in half for k in ks)
        mut = sum(1 for k in ks
                  if not exact_zero(amplitude_exact(N, c, 0, k) ** 2 - psi(k, 0) ** 2))
        check(f"G2 N={N}", hit and len(half) > 1,
              f"{len(half)} level-1/2 bonds {half}, all {len(ks)} modes equal "
              f"(2/(N+1))sin^2(2 pi k/(N+1)) exactly")
        check(f"M1 N={N} (psi_k(0) instead)", mut > 0,
              f"{mut}/{len(ks)} modes disagree with psi_k(0)^2")


def gate_G3():
    print("\nG3  the profile and the level structure, exactly, every bond")
    for N in (5, 7, 9, 11):
        c = (N - 1) // 2
        ks = blind_modes(N, c)
        # the individual c_k are algebraic and generally irrational (at N = 9 they are
        # (5 +- sqrt 5)/40); only their SUM is rational, so they are carried symbolically
        # and compared by exact_zero, while the sum is read as an exact Fraction.
        vecs = {b: tuple(sp.simplify(amplitude_exact(N, c, b, k) ** 2) for k in ks)
                for b in range(N - 1)}
        sums = {b: exact_rational(sum(v)) for b, v in vecs.items()}
        ok = all(sums[b] == level(N, b) for b in range(N - 1))
        m = (N - 1) // 2
        shifted = {b: Fraction(1, 2) * min((b if b < m else N - 2 - b) + 1,
                                           m - (b if b < m else N - 2 - b))
                   for b in range(N - 1)}
        mut = sum(1 for b in range(N - 1) if sums[b] != shifted[b])
        check(f"G3 N={N} profile", ok and None not in sums.values(),
              f"{[str(sums[b]) for b in range(N - 1)]}")
        check(f"M2 N={N} (shifted formula)", mut > 0,
              f"{mut}/{N - 1} bonds disagree with the shifted formula")
        check(f"G3 N={N} not constant", len(set(sums.values())) >= 2,
              f"{len(set(sums.values()))} distinct levels in the profile")
        # the level structure: the whole vector is a function of the level alone
        by_level = {}
        for b, v in vecs.items():
            by_level.setdefault(level(N, b), []).append(v)
        same = all(all(exact_zero(x - y) for v in vs for x, y in zip(v, vs[0]))
                   for vs in by_level.values())
        check(f"G3 N={N} level structure", same,
              f"{len(by_level)} levels, each with one c-vector across "
              f"{[sum(1 for b in range(N - 1) if level(N, b) == L) for L in sorted(by_level)]}"
              f" bonds")


def gate_G4():
    print("\nG4  one-seat light additivity w(u v^dagger) = a + b - 2ab (an ALGEBRAIC identity")
    print("    about the charged cell set: it cannot fail for a physical reason)")
    N, c = 7, 3
    rng = np.random.default_rng(20260912)
    bad = mut = overlapping = drawn = 0
    for _ in range(200):
        u = [Fraction(int(v), 7) for v in rng.integers(-9, 10, N)]
        v = [Fraction(int(x), 5) for x in rng.integers(-9, 10, N)]
        nu, nv = sum(x * x for x in u), sum(x * x for x in v)
        if nu == 0 or nv == 0:
            continue
        drawn += 1
        a, bb = u[c] * u[c] / nu, v[c] * v[c] / nv
        w = sum(u[i] * u[i] * v[j] * v[j] for i in range(N) for j in range(N)
                if (i == c) != (j == c)) / (nu * nv)
        bad += w != a + bb - 2 * a * bb
        if a * bb != 0:
            overlapping += 1
            mut += w != a + bb
    check("G4 rational pairs", bad == 0, f"{bad} violations in {drawn} drawn pairs")
    check("M3 (a+b, no overlap term)", mut == overlapping and overlapping > 0,
          f"{mut}/{overlapping} overlapping pairs disagree with a+b")


def gate_G5():
    print("\nG5  peripheral dimension = m^2 + 1, exact rank over GF(p), two primes")
    for N in (3, 5, 7, 9):
        c = m = (N - 1) // 2
        dims = [peripheral_dimension_exact(N, c, p) for p in PRIMES]
        check(f"G5 N={N}", all(d == m * m + 1 for d in dims),
              f"dims {dims} at two primes against m^2+1 = {m * m + 1}")
        if N >= 5:
            other = 0 if c != 0 else 1
            mut = [peripheral_dimension_exact(N, c, p, charged_seat=other) for p in PRIMES]
            check(f"M4 N={N} (charged cells moved to seat {other})",
                  all(d != m * m + 1 for d in mut),
                  f"the object changes: dims {mut} instead of {m * m + 1}")


def gate_G6():
    print("\nG6  the one-node criterion, by F157's own instrument in EXACT arithmetic:")
    print("    blind(seat) of the DETUNED chain = N - rank Krylov(e_seat) over Q.")
    print("    A mode blind at eps = 0 stays blind at eps != 0 exactly when it has a node")
    print("    at an end of the moved bond, so the surviving count is the node count.")
    agree = mismatch = 0
    nontrivial = 0
    for N in (5, 7, 9, 11):
        c = (N - 1) // 2
        for b in range(N - 1):
            if b in (c - 1, c):
                continue                      # G7's trivial cells, excluded as the source does
            pred = sum(1 for k in blind_modes(N, c)
                       if (k * (b + 1)) % (N + 1) == 0 or (k * (b + 2)) % (N + 1) == 0)
            for eps in (sp.Rational(7, 5), sp.Rational(-3, 11)):
                got = blind_count_exact(N, c, b, eps)
                if got == pred:
                    agree += 1
                else:
                    mismatch += 1
            nontrivial += pred != len(blind_modes(N, c))
    check("G6 exact blind count vs the node prediction", mismatch == 0 and agree > 0,
          f"{agree}/{agree + mismatch} (N, bond, eps) cells agree over N = 5..11 at two "
          f"rational eps; {nontrivial} of the bonds predict a count BELOW the unperturbed "
          f"one, so the criterion is not the trivial always-blind answer")


def gate_G7():
    print("\nG7  control and fence: a bond INCIDENT on the seat is struck away with it, so")
    print("    the blind count is unchanged at EVERY eps, not merely to second order.")
    print("    Owned by PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA, which calls these cells")
    print("    trivial and excludes them; gated here so this file cannot re-sell them.")
    for N in (7, 9, 11):
        c = m = (N - 1) // 2
        inc = [blind_count_exact(N, c, c - 1, e)
               for e in (sp.Rational(1, 3), sp.Rational(7, 5), sp.Integer(9))]
        far = blind_count_exact(N, c, 0, sp.Rational(1, 3))
        check(f"G7 N={N}", all(v == m for v in inc) and far != m,
              f"incident bond keeps {inc} = m = {m} at eps = 1/3, 7/5 and 9; the end bond "
              f"drops to {far}")


def gate_G8():
    print("\nG8  the rate statement, correctly scoped.  2 gamma (c_i + c_j) is the DIAGONAL")
    print("    of the second-order effective operator in the dyad basis.  It is the rate on")
    print("    every omega != 0 peripheral block, and NOT on the omega = 0 block.")
    gamma, eps = 0.3, 1e-4
    print("  G8a  vec(I) is an exact steady state at every eps, because z^2 = I")
    for N, e in ((7, 0.3), (9, 0.13), (11, 0.77)):
        c = (N - 1) // 2
        r = np.linalg.norm(A_generator(N, c, 0, e, 0.7) @ np.eye(N).reshape(-1))
        check(f"G8a N={N} eps={e}", r == 0.0, f"||A vec(I)|| = {r:.1e} (exactly 0.0 required)")
    print("  G8b  omega != 0 blocks: the diagonal IS the rate")
    print("  G8c  omega = 0 block: the diagonal is NOT the rate, and the mismatch is gated")
    for N in (7, 9, 11):
        c = (N - 1) // 2
        s = np.pi / (N + 1)
        cs = {k: (2.0 / (N + 1)) * np.sin(2 * s * k) ** 2 for k in blind_modes(N, c)}
        Es = {k: 4 * np.cos(s * k) for k in cs}
        A0 = A_generator(N, c, 0, 0.0, gamma)
        normA = np.linalg.norm(A0, 2)
        w0 = np.linalg.eigvals(A0)
        per = [k for k in range(N * N) if abs(w0[k].real) < 1e-10 * max(1.0, gamma)]
        w = np.linalg.eigvals(A_generator(N, c, 0, eps, gamma))
        used, meas = set(), []   # normA is used by G8c's resolution floor below
        for k in per:
            for idx in np.argsort(np.abs(w - w0[k])):
                if idx not in used:
                    used.add(idx)
                    meas.append((w0[k].imag, -w[idx].real / (gamma * eps * eps)))
                    break
        # off-diagonal dyads, grouped by frequency
        offdiag = {}
        for i in cs:
            for j in cs:
                if i == j:
                    continue
                offdiag.setdefault(round(Es[i] - Es[j], 6), []).append(2 * (cs[i] + cs[j]))
        ok_nonzero, checked = True, 0
        for omega, pred in offdiag.items():
            if abs(omega) < 1e-9:
                continue
            got = sorted(r for o, r in meas if abs(o - omega) < 1e-6)
            if len(got) != len(pred):
                ok_nonzero = False
                continue
            checked += len(got)
            for a, b in zip(got, sorted(pred)):
                if abs(a - b) > 3e-3 * max(1.0, b):
                    ok_nonzero = False
        check(f"G8b N={N}", ok_nonzero and checked > 0,
              f"{checked} rates on omega != 0 blocks match 2 gamma (c_i+c_j) "
              f"within the eps-germ")
        # The omega = 0 diagonal multiset has m+1 entries: 4 c_k for each diagonal dyad,
        # and I_E's own entry, which is 4*(sum_k c_k)/(m+1) (its first-order change is
        # -sum_i (d_i delta_i^dag + h.c.), so its charged weight is 2 sum_k c_k over
        # ||I_E||^2 = m+1).  Comparing against the m dyad entries alone could not fail,
        # because the lists have different lengths by construction.
        m = (N - 1) // 2
        zero_meas = sorted(round(r, 4) for o, r in meas if abs(o) < 1e-9)
        zero_diag = sorted([round(4 * cs[k], 4) for k in cs]
                           + [round(4 * sum(cs.values()) / (m + 1), 4)])
        differ = (len(zero_meas) != len(zero_diag)
                  or any(abs(a - b) > 1e-3 for a, b in zip(zero_meas, zero_diag)))
        # vec(I) has rate exactly 0; that exactness is G8a's, tested on the generator.
        # Here the rate comes from an eigensolver, so the honest statement is that the
        # block carries a rate below what the read can resolve: NOISE_FACTOR * eps_machine
        # * ||A||, divided by the gamma eps^2 the column is normalised by.
        floor = NOISE_FACTOR * EPSM * normA / (gamma * eps * eps)
        smallest = min(abs(r) for o, r in meas if abs(o) < 1e-9)
        check(f"G8c N={N}", differ and smallest < floor,
              f"omega=0 measured {zero_meas} vs the full diagonal multiset {zero_diag} "
              f"(m dyads + I_E, same length): they differ; the smallest measured rate is "
              f"{smallest:.1e} against a resolution floor of {floor:.1e}")


def gate_G9():
    """The criterion at ALL orders, by the route the parent named and did not take."""
    print("\nG9  Corollary B at ALL orders, not merely first.  Expand the left block's")
    print("    characteristic polynomial along the moved bond:")
    print("      chi_L = P*Q - (J(1+eps))^2 * P'*Q',   P = chi(0..b), Q = chi(b+1..c-1),")
    print("      P' = chi(0..b-1), Q' = chi(b+2..c-1).")
    print("    A node at b splits the left block there, so by (J3) E_k is a root of both")
    print("    pieces: P'(E_k) = Q(E_k) = 0, and BOTH terms vanish.  A node at b+1 gives")
    print("    P(E_k) = Q'(E_k) = 0, same conclusion.  So chi_L(E_k; eps) is identically")
    print("    zero in eps and the mode keeps its energy AND its node at every eps.")
    lam, e = sp.Symbol("lam"), sp.Symbol("e")

    def chi_left(N, c, b, eps):
        t = [sp.Integer(2)] * (N - 1)
        t[b] = 2 * (1 + eps)
        M = sp.zeros(c, c)
        for j in range(c - 1):
            M[j, j + 1] = M[j + 1, j] = t[j]
        return sp.expand((lam * sp.eye(c) - M).det())

    hit = miss = 0
    control_nonzero = control_tot = 0
    for N in (5, 7, 9, 11, 13):
        c = (N - 1) // 2
        s = sp.pi / (N + 1)
        for b in range(c - 1):                       # left half, non-incident
            P = chi_left(N, c, b, e)
            for k in range(2, N + 1, 2):
                Ek = 4 * sp.cos(s * k)
                val = sp.simplify(sp.expand(P.subs(lam, Ek)))
                node = ((k * (b + 1)) % (N + 1) == 0) or ((k * (b + 2)) % (N + 1) == 0)
                if node:
                    if sp.simplify(val) == 0:
                        hit += 1
                    else:
                        miss += 1
                else:
                    control_tot += 1
                    control_nonzero += sp.simplify(val) != 0
    check("G9 all-orders criterion", miss == 0 and hit > 0,
          f"{hit}/{hit + miss} node-carrying (N, bond, mode) cells give chi_L(E_k; eps) "
          f"identically 0 in eps over N = 5..13")
    check("G9 control (no node at either end)", control_nonzero == control_tot
          and control_tot > 0,
          f"{control_nonzero}/{control_tot} cells WITHOUT a node are not identically zero, "
          f"so the criterion is not a property of chi_L alone")


# -------------------------------------------------------------------------------- reads

def read_R1():
    print("\nR1  READ: the cubic germ.  slowest/(2 gamma eps^2 min(c_i+c_j)) = 1 + kappa_N eps.")
    print("    kappa_7 = +1 is the published N=7 germ; the others are read, not assumed.")
    gamma = 0.3
    for N in (5, 7, 9, 11, 13):
        c = (N - 1) // 2
        s = np.pi / (N + 1)
        cs = [(2.0 / (N + 1)) * np.sin(2 * s * k) ** 2 for k in range(2, N + 1, 2)]
        pred = min(p_ for p_ in sorted({round(x + y, 12) for x in cs for y in cs})
                   if p_ > 1e-12)
        A0 = A_generator(N, c, 0, 0.0, gamma)
        normA = np.linalg.norm(A0, 2)
        w0 = np.linalg.eigvals(A0)
        per = [k for k in range(N * N) if abs(w0[k].real) < 1e-10 * max(1.0, gamma)]
        out = []
        for eps in (1e-2, 1e-3, 1e-4):
            w = np.linalg.eigvals(A_generator(N, c, 0, eps, gamma))
            used, rates = set(), []
            for k in per:
                for idx in np.argsort(np.abs(w - w0[k])):
                    if idx not in used:
                        used.add(idx)
                        rates.append(-w[idx].real)
                        break
            pos = sorted(r for r in rates if r > NOISE_FACTOR * EPSM * normA)
            ratio = pos[0] / (2 * gamma * pred * eps * eps)
            out.append((eps, ratio, (ratio - 1.0) / eps))
        print(f"    N={N:3d}  2*min(c_i+c_j)={2 * pred:.8f}")
        for eps, ratio, kappa in out:
            print(f"          eps={eps:<8g} ratio={ratio:.8f}   kappa_N = {kappa:+.4f}")
        if N == 5:
            print("          (N=5's three reads disagree by an order of magnitude: the rate")
            print("           is ~1e-9 at eps=1e-4, so this row reads kappa_5 = 0 plus noise)")


def read_R2():
    print("\nR2  READ: two regimes, and what each depends on.")
    N, c = 7, 3
    print("\n    (a) the node law's coefficient, over gamma and over J")
    for J in (0.5, 1.0, 2.0, 5.0):
        row = []
        for gamma in (0.03, 0.3, 3.0):
            A0 = A_generator(N, c, 0, 0.0, gamma, J)
            normA = np.linalg.norm(A0, 2)
            w0 = np.linalg.eigvals(A0)
            per = [k for k in range(N * N) if abs(w0[k].real) < 1e-10 * max(1.0, gamma)]
            w = np.linalg.eigvals(A_generator(N, c, 0, 1e-4, gamma, J))
            used, rates = set(), []
            for k in per:
                for idx in np.argsort(np.abs(w - w0[k])):
                    if idx not in used:
                        used.add(idx)
                        rates.append(-w[idx].real)
                        break
            pos = sorted(r for r in rates if r > NOISE_FACTOR * EPSM * normA)
            row.append(pos[0] / (gamma * 1e-8))
        print(f"      J={J:<5} slowest/(gamma eps^2) at gamma = 0.03, 0.3, 3: "
              f"{['%.6f' % v for v in row]}")
    print("      -> the coefficient is free of BOTH gamma and J.")
    print("\n    (b) the saturated rate, deep beyond the window (eps = 0.3)")
    for J in (0.5, 1.0, 2.0, 4.0):
        vals = []
        for gamma in (300.0, 1000.0):
            A0 = A_generator(N, c, 0, 0.0, gamma, J)
            normA = np.linalg.norm(A0, 2)
            w0 = np.linalg.eigvals(A0)
            per = [k for k in range(N * N) if abs(w0[k].real) < 1e-10 * max(1.0, gamma)]
            w = np.linalg.eigvals(A_generator(N, c, 0, 0.3, gamma, J))
            used, rates = set(), []
            for k in per:
                for idx in np.argsort(np.abs(w - w0[k])):
                    if idx not in used:
                        used.add(idx)
                        rates.append(-w[idx].real)
                        break
            pos = sorted(r for r in rates if r > NOISE_FACTOR * EPSM * normA)
            vals.append(pos[0] * gamma / (J * J))
        print(f"      J={J:<5} rate*gamma/J^2 at gamma = 300, 1000: {['%.6f' % v for v in vals]}")
    print("      -> at N = 7 this is a genuine plateau: A(7) = 1/8, flat in eps.  The next")
    print("         block prints the same quantity at other N, and it is NOT a plateau")
    print("         everywhere: sweeping eps at gamma = 1000 gives a flat value at N = 5, 7")
    print("         and (up to eps ~ 0.3) at N = 11, and a value that falls by 3.5x and 9x")
    print("         across the same sweep at N = 9 and 13.  So the numbers below are reads")
    print("         at eps = 0.3, not constants, except at N = 5 and 7.")
    for Nn in (5, 7, 9, 11):
        cc = (Nn - 1) // 2
        gamma, J = 1000.0, 2.0
        A0 = A_generator(Nn, cc, 0, 0.0, gamma, J)
        normA = np.linalg.norm(A0, 2)
        w0 = np.linalg.eigvals(A0)
        per = [k for k in range(Nn * Nn) if abs(w0[k].real) < 1e-10 * gamma]
        w = np.linalg.eigvals(A_generator(Nn, cc, 0, 0.3, gamma, J))
        used, rates = set(), []
        for k in per:
            for idx in np.argsort(np.abs(w - w0[k])):
                if idx not in used:
                    used.add(idx)
                    rates.append(-w[idx].real)
                    break
        pos = sorted(r for r in rates if r > NOISE_FACTOR * EPSM * normA)
        print(f"      N={Nn}: A(N) = rate*gamma/J^2 = {pos[0] * gamma / (J * J):.6f}")
    print("\n    (c) the two laws cross where they agree: (1/2) gamma eps^2 = A(N) J^2/gamma,")
    print("        i.e. gamma*eps = J*sqrt(2 A(N)), which at N=7 is gamma*eps = J/2.")
    print("    TWO FENCES ON (b) AND (c).  First, the branch this read follows is the one")
    print("    continued from eps = 0, and in the saturated regime it is NOT the block's")
    print("    slowest: at gamma=300, eps=0.1, J=2 the tracked branch gives rate*gamma =")
    print("    0.500 while the block's slowest gives 0.4525.  Second, the saturated law is")
    print("    read here, not derived, and it is not identified with any owned closed form.")
    print("    The Absorption Theorem's Zeno asymptote is for UNIFORM dephasing, and its J is")
    print("    the PAULI coupling whose hopping element is 2J, so a comparison has to cross a")
    print("    convention the repo has been caught by before; it is not made here.")


# --------------------------------------------------------------------------------- main

if __name__ == "__main__":
    print("The reduced resolvent between two nodes: gates and reads")
    for g in (gate_G1, gate_G2, gate_G3, gate_G4, gate_G5, gate_G6, gate_G7, gate_G8, gate_G9):
        g()
    read_R1()
    read_R2()
    npass = sum(1 for _, ok, _ in RESULTS if ok)
    print(f"\n{npass}/{len(RESULTS)} gates pass.")
    if npass != len(RESULTS):
        print("FAILURES:")
        for name, ok, detail in RESULTS:
            if not ok:
                print(f"  {name}: {detail}")
        raise SystemExit(1)
