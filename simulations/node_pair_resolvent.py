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
    G2   c_k = psi_k(1)^2 = (2/(N+1))*sin^2(2*pi*k/(N+1)) on the bonds whose
         profile level is 1/2.  Exact, per mode, directly against the sine eigenvector.
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
    G6   the one-node sufficient condition, by the EXACT route: F157's blind count of the
         DETUNED chain by Krylov rank over Q at two rational eps, against the node lower
         bound, at the CENTRE seat and every non-incident bond.  The sampled generic-count
         equality is a reading, not a converse.  The epsilon=-2 sign-gauge return is an
         exact control where extra blind modes survive.
    G7   control, and a fence on this file's own scope: a bond INCIDENT on the watched seat
         is struck away with the seat, so both principal blocks are eps-free and the blind
         count is unchanged at EVERY eps, not merely to second order.  Owned and called
         trivial by PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA; gated here so this file cannot
         re-sell it.
    G9   Corollary B at ALL orders on the zero-free knob domain: a node at an end of the
         moved bond makes the left block's characteristic polynomial vanish at E_k
         IDENTICALLY in eps, so the mode keeps its energy and its node for eps != -1, not
         merely to second order.  On the UNIFORM CENTRE-WATCHED family the exact determinant
         factorization closes the pointwise converse for r=1+eps outside r=0,+1,-1.  The
          control requires the no-node cells NOT to be identically zero, using an exact
          coefficient-wise polynomial oracle with a seven-root zero identity as its own
          fail-open control.  G9 fences the three exceptional ratios, and M8 rejects a
          linear-in-r mutation of the factorization.
    G8   the rate statement, correctly scoped: 2 gamma (c_i + c_j) is the DIAGONAL of the
         second-order effective operator in the dyad basis.  It is the rate on every
         omega != 0 peripheral block (gated against the generator), and it is NOT on the
         omega = 0 block, which always carries vec(I) at rate exactly 0 because z^2 = I.
         The omega = 0 mismatch at N = 9 is gated as a POSITIVE fact, so a version of this
         file that "fixed" the mismatch would go red.  M5 rejects a truncated acquisition;
         M6a-c reject full-length acquisitions containing NaN or either infinity.  G8b also
          validates the whole acquisition before matching.  Its three-acquisition germ gate
          spans two decades and requires the normalized error to shrink linearly in epsilon,
          after an explicit eigensolver-rounding reserve.  M7a-c replace a nonzero-frequency
          branch with NaN or either infinity, and M9 scales every prediction by 1.003; all
          must be rejected through that same multi-decade path.

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


def exact_polynomial_nonzero(expr, symbol):
    """True iff a polynomial in `symbol` has a coefficient certified nonzero.

    `simplify(expr) != 0` is fail-open for unsimplified algebraic zeroes.  A polynomial is
    identically zero exactly when every algebraic coefficient is certified zero.
    """
    coefficients = sp.Poly(sp.expand(expr), symbol).all_coeffs()
    return any(not exact_zero(coefficient) for coefficient in coefficients)


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


def complete_zero_frequency_mismatch(measured, predicted, expected_count, tolerance=1e-3):
    """Whether two complete equal-size zero-frequency multisets disagree.

    Missing or non-finite branches are failed acquisition, not evidence of disagreement.
    """
    return (len(measured) == expected_count
            and len(predicted) == expected_count
            and all(np.isfinite(value) for values in (measured, predicted) for value in values)
            and any(abs(a - b) > tolerance for a, b in zip(measured, predicted)))


def nonzero_frequency_normalized_error(measured, predicted_by_frequency,
                                       frequency_tolerance=1e-6):
    """Maximum normalized rate error after a complete, finite frequency match."""
    expected_count = sum(len(rates) for omega, rates in predicted_by_frequency.items()
                         if np.isfinite(omega) and abs(omega) >= frequency_tolerance)
    acquired_count = sum(1 for omega, _ in measured
                         if np.isfinite(omega) and abs(omega) >= frequency_tolerance)
    if (not all(np.isfinite(omega) and np.isfinite(rate) for omega, rate in measured)
            or not all(np.isfinite(omega) and np.isfinite(rate)
                       for omega, rates in predicted_by_frequency.items() for rate in rates)):
        return None, acquired_count, expected_count

    predicted = {omega: rates for omega, rates in predicted_by_frequency.items()
                 if abs(omega) >= frequency_tolerance}
    acquired = [(omega, rate) for omega, rate in measured
                if abs(omega) >= frequency_tolerance]
    if len(acquired) != expected_count:
        return None, len(acquired), expected_count

    errors = []
    for omega, rates in predicted.items():
        got = sorted(rate for measured_omega, rate in acquired
                     if abs(measured_omega - omega) < frequency_tolerance)
        if len(got) != len(rates):
            return None, len(acquired), expected_count
        for actual, expected in zip(got, sorted(rates)):
            errors.append(abs(actual - expected) / max(1.0, abs(expected)))
    if len(errors) != expected_count:
        return None, len(acquired), expected_count
    return max(errors, default=0.0), len(acquired), expected_count


def complete_nonzero_frequency_match(measured, predicted_by_frequency,
                                     frequency_tolerance=1e-6, rate_tolerance=3e-3):
    """Whether the complete nonzero-frequency acquisition matches its grouped prediction."""
    error, acquired_count, expected_count = nonzero_frequency_normalized_error(
        measured, predicted_by_frequency, frequency_tolerance)
    return error is not None and error <= rate_tolerance, acquired_count, expected_count


def complete_nonzero_frequency_germ(acquisitions, predicted_by_frequency, gamma,
                                    linear_error_budget=2.0, shrink_slack=1.25):
    """Gate complete finite rate acquisitions at three or more points spanning two decades.

    Each tuple is (epsilon, measured, ||A_0||).  In normalized rate units, each acquisition
    gets the explicit error budget `linear_error_budget*|epsilon| + rounding_reserve`, where
    `rounding_reserve = NOISE_FACTOR*eps_machine*||A_0||/(gamma*epsilon^2)`.  After that
    reserve is removed, adjacent maximum errors must shrink linearly with epsilon; the named
    `shrink_slack` allows for sampling a finite interval rather than the asymptotic limit.
    """
    if (not np.isfinite(gamma) or gamma <= 0
            or any(not np.isfinite(eps) or eps == 0 or not np.isfinite(norm_a) or norm_a < 0
                   for eps, _, norm_a in acquisitions)):
        return False, []

    rows = []
    for eps, measured, norm_a in sorted(acquisitions, key=lambda item: -abs(item[0])):
        rounding_reserve = NOISE_FACTOR * EPSM * norm_a / (gamma * eps * eps)
        budget = linear_error_budget * abs(eps) + rounding_reserve
        matched, acquired_count, expected_count = complete_nonzero_frequency_match(
            measured, predicted_by_frequency, rate_tolerance=budget)
        error, _, _ = nonzero_frequency_normalized_error(measured, predicted_by_frequency)
        rows.append({"eps": abs(eps), "matched": matched, "acquired": acquired_count,
                     "expected": expected_count, "error": error, "budget": budget,
                     "rounding_reserve": rounding_reserve})

    distinct_epsilons = {row["eps"] for row in rows}
    spans_two_decades = (len(rows) >= 3 and len(distinct_epsilons) == len(rows)
                         and rows[0]["eps"] / rows[-1]["eps"] >= 100.0 * (1.0 - 1e-12))
    if (not spans_two_decades
            or any(not row["matched"] or row["error"] is None for row in rows)):
        return False, rows

    for large, small in zip(rows, rows[1:]):
        epsilon_ratio = small["eps"] / large["eps"]
        large_resolved = max(0.0, large["error"] - large["rounding_reserve"])
        small_resolved = max(0.0, small["error"] - small["rounding_reserve"])
        if small_resolved > shrink_slack * epsilon_ratio * large_resolved:
            return False, rows
    return True, rows


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
    print("\nG2  c_k = psi_k(1)^2 = (2/(N+1)) sin^2(2 pi k/(N+1)) on the level-1/2 bonds")
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
    print("\nG6  the one-node sufficient condition, by F157's own instrument in EXACT arithmetic:")
    print("    blind(seat) of the DETUNED chain = N - rank Krylov(e_seat) over Q.")
    print("    A node at an end of the moved bond guarantees survival for eps != -1.")
    print("    Equality with that lower bound at the two sampled generic knobs is a reading,")
    print("    not the converse: eps = -2 is a sign-gauge return with extra blind modes.")
    lower_bound_holds = lower_bound_fails = 0
    sampled_equal = sampled_unequal = 0
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
                if got >= pred:
                    lower_bound_holds += 1
                else:
                    lower_bound_fails += 1
                if got == pred:
                    sampled_equal += 1
                else:
                    sampled_unequal += 1
            nontrivial += pred != len(blind_modes(N, c))
    check("G6 sufficient node lower bound", lower_bound_fails == 0 and lower_bound_holds > 0,
          f"{lower_bound_holds}/{lower_bound_holds + lower_bound_fails} (N, bond, eps) cells "
          f"meet blind count >= node count over N = 5..11 at eps = 7/5 and -3/11")
    check("G6 sampled generic-count reading", sampled_unequal == 0 and sampled_equal > 0,
          f"{sampled_equal}/{sampled_equal + sampled_unequal} sampled cells equal the node lower "
          f"bound; {nontrivial} bonds put that bound below the unperturbed count")

    sign_gauge_got = blind_count_exact(5, 2, 0, sp.Integer(-2))
    sign_gauge_pred = sum(1 for k in blind_modes(5, 2)
                          if (k * 1) % 6 == 0 or (k * 2) % 6 == 0)
    check("G6 control epsilon=-2 sign-gauge return",
          sign_gauge_pred == 0 and sign_gauge_got == 2,
          f"N=5, seat=2, bond=0: node lower bound {sign_gauge_pred}, exact blind count "
          f"{sign_gauge_got}; changing one bond from +2 to -2 is gauge-equivalent")


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
    gamma = 0.3
    germ_epsilons = (1e-2, 1e-3, 1e-4)
    print("  G8a  float implementation control for the algebraic identity A[I] = 0")
    for N, e in ((7, 0.3), (9, 0.13), (11, 0.77)):
        c = (N - 1) // 2
        r = np.linalg.norm(A_generator(N, c, 0, e, 0.7) @ np.eye(N).reshape(-1))
        check(f"G8a N={N} eps={e}", r == 0.0,
              f"float64 construction gives ||A vec(I)|| = {r:.1e} (bit-zero required)")
    print("  G8b multi-decade germ: omega != 0 blocks have the predicted rate")
    print("       in three acquisitions spanning two decades; normalized error must shrink")
    print("       linearly in epsilon after the explicit eigensolver-rounding reserve")
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
        acquisitions = []
        for eps in germ_epsilons:
            w = np.linalg.eigvals(A_generator(N, c, 0, eps, gamma))
            used, meas = set(), []
            for k in per:
                for idx in np.argsort(np.abs(w - w0[k])):
                    if idx not in used:
                        used.add(idx)
                        meas.append((w0[k].imag, -w[idx].real / (gamma * eps * eps)))
                        break
            acquisitions.append((eps, meas, normA))
        # off-diagonal dyads, grouped by frequency
        offdiag = {}
        for i in cs:
            for j in cs:
                if i == j:
                    continue
                offdiag.setdefault(round(Es[i] - Es[j], 6), []).append(2 * (cs[i] + cs[j]))
        ok_nonzero, germ_rows = complete_nonzero_frequency_germ(
            acquisitions, offdiag, gamma)
        germ_detail = "; ".join(
            f"eps={row['eps']:.0e}: {row['acquired']}/{row['expected']}, "
            f"error={row['error'] if row['error'] is None else format(row['error'], '.3e')}, "
            f"budget={row['budget']:.3e}, reserve={row['rounding_reserve']:.3e}"
            for row in germ_rows)
        check(f"G8b N={N} three-point germ",
              ok_nonzero and germ_rows[0]["expected"] > 0, germ_detail)
        if N == 7:
            middle_eps, middle_meas, _ = acquisitions[1]
            nonzero_index = next(index for index, (omega, _) in enumerate(middle_meas)
                                 if abs(omega) >= 1e-6)
            for mutation_name, label, nonfinite in (
                    ("M7a G8b rejects full-length NaN", "NaN", np.nan),
                    ("M7b G8b rejects full-length +Inf", "+Inf", np.inf),
                    ("M7c G8b rejects full-length -Inf", "-Inf", -np.inf)):
                mutated = list(middle_meas)
                omega, _ = mutated[nonzero_index]
                mutated[nonzero_index] = (omega, nonfinite)
                mutated_acquisitions = list(acquisitions)
                mutated_acquisitions[1] = (middle_eps, mutated, normA)
                accepted, mutation_rows = complete_nonzero_frequency_germ(
                    mutated_acquisitions, offdiag, gamma)
                mutation_row = mutation_rows[1]
                check(mutation_name, not accepted,
                      f"replaced one nonzero-frequency rate at omega={omega:.6f}; "
                      f"acquisition remained {mutation_row['acquired']}/"
                      f"{mutation_row['expected']}; multi-decade path accepted={accepted}")
            scaled_prediction = {omega: [1.003 * rate for rate in rates]
                                 for omega, rates in offdiag.items()}
            mutation_accepted, mutation_rows = complete_nonzero_frequency_germ(
                acquisitions, scaled_prediction, gamma)
            check("M9 G8b constant 1.003 coefficient", not mutation_accepted,
                  f"same three-acquisition path accepted={mutation_accepted}; errors="
                  f"{[None if row['error'] is None else float(row['error']) for row in mutation_rows]}")
        eps, meas, _ = acquisitions[-1]  # G8c uses the smallest-epsilon acquisition below
        # The omega = 0 diagonal multiset has m+1 entries: 4 c_k for each diagonal dyad,
        # and I_E's own entry, which is 4*(sum_k c_k)/(m+1) (its first-order change is
        # -sum_i (d_i delta_i^dag + h.c.), so its charged weight is 2 sum_k c_k over
        # ||I_E||^2 = m+1).  Comparing against the m dyad entries alone could not fail,
        # because the lists have different lengths by construction.
        m = (N - 1) // 2
        zero_meas = sorted(round(r, 4) for o, r in meas if abs(o) < 1e-9)
        zero_diag = sorted([round(4 * cs[k], 4) for k in cs]
                           + [round(4 * sum(cs.values()) / (m + 1), 4)])
        expected_zero = m + 1
        differ = complete_zero_frequency_mismatch(zero_meas, zero_diag, expected_zero)
        # Algebra gives vec(I) exact rate 0; G8a only checks the float construction.
        # Here the rate comes from an eigensolver, so the honest statement is that the
        # block carries a rate below what the read can resolve: NOISE_FACTOR * eps_machine
        # * ||A||, divided by the gamma eps^2 the column is normalised by.
        floor = NOISE_FACTOR * EPSM * normA / (gamma * eps * eps)
        smallest = min(abs(r) for o, r in meas if abs(o) < 1e-9)
        check(f"G8c N={N}", differ and smallest < floor,
              f"omega=0 measured {zero_meas} vs the full diagonal multiset {zero_diag} "
              f"(m dyads + I_E, {len(zero_meas)}/{expected_zero} branches acquired): they "
              f"differ; the smallest measured rate is "
              f"{smallest:.1e} against a resolution floor of {floor:.1e}")
        truncated = [min(zero_meas, key=abs)]
        check(f"M5 G8c N={N} (retain only stationary branch)",
              not complete_zero_frequency_mismatch(truncated, zero_diag, expected_zero),
              f"truncated acquisition has {len(truncated)}/{expected_zero} branches and is "
              f"rejected before spectral comparison")

    finite_measured = [0.0, 0.25, 0.75]
    finite_predicted = [0.0, 0.5, 0.5]
    finite_control = complete_zero_frequency_mismatch(finite_measured, finite_predicted, 3)
    for label, nonfinite in (("NaN", np.nan), ("+Inf", np.inf), ("-Inf", -np.inf)):
        measured_mutation = [0.0, nonfinite, 0.75]
        predicted_mutation = [0.0, nonfinite, 0.5]
        measured_accepted = complete_zero_frequency_mismatch(measured_mutation, finite_predicted, 3)
        predicted_accepted = complete_zero_frequency_mismatch(finite_measured, predicted_mutation, 3)
        check(f"M6 G8c rejects full-length {label}",
              finite_control and not measured_accepted and not predicted_accepted,
              f"finite control accepted={finite_control}; 3/3 measured mutation accepted="
              f"{measured_accepted}; 3/3 predicted mutation accepted={predicted_accepted}")


def gate_G9():
    """The all-orders sufficient direction and the uniform-centre pointwise converse."""
    print("\nG9  Corollary B at ALL orders for eps != -1, not merely first.  Expand the left block's")
    print("    characteristic polynomial along the moved bond:")
    print("      chi_L = P*Q - (J(1+eps))^2 * P'*Q',   P = chi(0..b), Q = chi(b+1..c-1),")
    print("      P' = chi(0..b-1), Q' = chi(b+2..c-1).")
    print("    A node at b splits the left block there, so by (J3) E_k is a root of both")
    print("    pieces: P'(E_k) = Q(E_k) = 0, and BOTH terms vanish.  A node at b+1 gives")
    print("    P(E_k) = Q'(E_k) = 0, same conclusion.  So chi_L(E_k; eps) is identically")
    print("    zero in eps; while the chain stays zero-free, the mode keeps its energy AND node.")
    lam, e, r = sp.symbols("lam e r")

    seven_root_zero = sum(sp.cos(2 * sp.pi * k / 7) for k in range(1, 7)) + 1
    seven_root_nonzero = sum(sp.cos(2 * sp.pi * k / 7) for k in range(1, 6)) + 1
    zero_classified_nonzero = exact_polynomial_nonzero((1 + e) * seven_root_zero, e)
    nonzero_classified_nonzero = exact_polynomial_nonzero((1 + e) * seven_root_nonzero, e)
    check("G9-control exact polynomial oracle (seven-root zero identity)",
          not zero_classified_nonzero and nonzero_classified_nonzero,
          f"six-root sum + 1 classified nonzero={zero_classified_nonzero}; the genuine "
          f"five-root nonzero trigonometric sequence classified nonzero="
          f"{nonzero_classified_nonzero}")

    def path_poly(n):
        """Characteristic polynomial of the uniform n-site path with hopping 2."""
        p0, p1 = sp.Integer(1), lam
        if n == 0:
            return p0
        if n == 1:
            return p1
        for _ in range(2, n + 1):
            p0, p1 = p1, sp.expand(lam * p1 - 4 * p0)
        return p1

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
                    if not exact_polynomial_nonzero(val, e):
                        hit += 1
                    else:
                        miss += 1
                else:
                    control_tot += 1
                    control_nonzero += exact_polynomial_nonzero(val, e)
    check("G9 all-orders sufficient identity", miss == 0 and hit > 0,
          f"{hit}/{hit + miss} node-carrying (N, bond, mode) cells give chi_L(E_k; eps) "
          f"identically 0 in eps over N = 5..13")
    check("G9 control (no node at either end)", control_nonzero == control_tot
          and control_tot > 0,
          f"{control_nonzero}/{control_tot} cells WITHOUT a node are not identically zero, "
          f"so the criterion is not a property of chi_L alone")

    factorized = mutated = total = 0
    for m in range(2, 7):
        N, c = 2 * m + 1, m
        for b in range(m - 1):
            direct = sp.expand(chi_left(N, c, b, r - 1))
            correction = 4 * (r ** 2 - 1) * path_poly(b) * path_poly(m - b - 2)
            formula = sp.expand(path_poly(m) - correction)
            wrong = sp.expand(path_poly(m)
                              - 4 * (r - 1) * path_poly(b) * path_poly(m - b - 2))
            total += 1
            factorized += sp.expand(direct - formula) == 0
            mutated += sp.expand(direct - wrong) != 0
    check("G9b uniform-centre determinant factorization", factorized == total and total > 0,
          f"{factorized}/{total} left nonincident (m,b) cells exactly satisfy "
          f"chi_L=P_m-4(r^2-1)P_b P_(m-b-2)")
    check("M8 G9b linear-in-r correction", mutated == total,
          f"{mutated}/{total} cells reject replacing r^2-1 by r-1")

    iff_hit = iff_total = node_positive = node_zero = 0
    for N in (5, 7, 9, 11, 13):
        c = (N - 1) // 2
        for b in range(c - 1):
            node_count = sum(1 for k in blind_modes(N, c)
                             if (k * (b + 1)) % (N + 1) == 0
                             or (k * (b + 2)) % (N + 1) == 0)
            node_positive += node_count > 0
            node_zero += node_count == 0
            for ratio in (sp.Rational(2, 3), sp.Rational(3, 2), sp.Integer(-2)):
                got = blind_count_exact(N, c, b, ratio - 1)
                iff_total += 1
                iff_hit += got == node_count
    check("G9c uniform-centre pointwise iff outside r=0,+1,-1",
          iff_hit == iff_total and node_positive > 0 and node_zero > 0,
          f"{iff_hit}/{iff_total} exact Krylov counts equal the endpoint-node count at "
          f"r=2/3,3/2,-2; bond populations include {node_positive} positive and "
          f"{node_zero} zero-node cells")

    N, c, b = 7, 3, 0
    exceptional = {ratio: blind_count_exact(N, c, b, sp.Integer(ratio) - 1)
                   for ratio in (-1, 0, 1)}
    check("G9d exceptional-set fence r=0,+1,-1",
          exceptional[-1] == c and exceptional[1] == c and exceptional[0] < c,
          f"N=7, seat=3, bond=0 exact blind counts are {exceptional}: r=+1 is the "
          f"unperturbed chain, r=-1 is its sign gauge, and r=0 cuts the zero-free chain")


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
