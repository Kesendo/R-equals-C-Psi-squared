"""The odd orders of an F129 collision gap on F160's road, and where they stop.

STATEMENT (derived 2026-09-02, proved in docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md):

  The open chain of N sites gets its wrap bond switched on at strength u (F160's road,
  n = N + 1).  A level moves as E_k(u) = E_k + d_1 u + d_2 u^2 + ..., and each of the
  FIVE COMPUTED coefficients is, in closed form, a signed combination of neighbouring
  evaluations of the comb under an integer MULTIPLIER: a difference of two through third
  order, three rungs at fifth.  The shape past the fifth is open.  Writing eta_k = (-1)^(k+1) for the chain reflection's sign on mode k
  (F75; carried in PROOF_K_PARTNERSHIP, derived in PROOF_C1_MIRROR_SYMMETRY) and

      M_j(tau)   = sum_{k in tau} cos(j k pi / n)
      X_2j(tau)  = sum_{k in tau} eta_k cos(2 j k pi / n)  =  - M_{n+2j}(tau),

  the gap of a collision pair (tau, sigma), Delta of any of these being tau minus sigma:

      c_1 = (2/n) (dX_0 - dX_2)                                   [F160 Theorem G]
      c_2 = ((n-3)/n^2) (dM_1 - dM_3) = -((n-3)/n^2) dM_3         [dM_1 = 0 IS the collision]
      c_3 = (2(n-2)(n-4)/(3 n^3)) (dX_2 - dX_4)
      c_5 = F [ (5/4)(n-1) dX_2 + (3n^2/2 - 8n + 8) dX_4 - (3/4)(2n-3)(n-3) dX_6 ],
            F = 4(n-6)(n-2)/(15 n^5)

  X_2j is the comb read under the multiplier k -> k(n+2j), and that map is a Galois
  automorphism of Q(zeta_2n) exactly when gcd(n+2j, 2n) = 1, which at ODD n is
  gcd(j, n) = 1.  So at odd n every collision has dX_2 = dX_4 = 0 by a theorem, hence
  c_3 = 0 for EVERY collision pair, standing or separating.  At even n no rung is a
  Galois image of dM_1 at all, and the twelve non-mirror standing pairs of the census are
  carried instead by the Conway-Jones ROT3 shape (F89_SEED_EXISTENCE_REDUCTION, Seed.cs's
  TripleFamily.Rot3): if a PARITY-UNIFORM triple's +-label set mod n is a union of two
  cosets of the order-3 subgroup, then X_2j = 0 is FORCED for every j with 3 not dividing
  j.  Parity-uniformity is free at even n (6|n makes every coset parity-homogeneous) and
  load-bearing at odd n, where (1,2,4) at n = 9 is ROT3, mixed, and has X_2 nonzero.
  The shape gives c_3 = 0; c_1 = (4/n)(o_tau - o_sigma) needs the two triples to share a
  parity class as well, and 58 collision pairs of this census are both-ROT3 with c_1 != 0.

  Where it stops is j = 3 in both cases, for two different reasons.  At odd n: F129 fires
  only at 3|n or 10|n and 10|n forces n even, so EVERY odd firing modulus has 3|n (one line
  of arithmetic on the condition, not a measurement) and j = 3 is the first multiplier that
  is not an automorphism.  At even n: 3|j is what collapses the ROT3 coset.  X_6 is the
  first surviving rung and it CAN enter at FIFTH order.  Whether it does is dX_6 != 0, an
  exact decision per pair rather than a theorem about all of them; it holds on all 212 and
  fails on exactly the 11 Theta-mirror pairs.  That accounts for the measured "the odd part
  starts at u^5" of experiments/THE_COMB_ON_THE_ROAD.md.

  WORD FENCE.  "Multiplier ladder" here is the sequence of evaluations M_{n+2j} of ONE
  comb and nothing else: not F142's spin ladder (eta_ladder_breakinput.py), not
  Pi2DyadicLadderClaim's dyadic ladder, not ClockHandLadderClaim's.  "Rung" is one j, and
  NOT the glossary's rung, which is a rate (2 gamma, 2N gamma) or a Q threshold, nor the
  arc sideways_spin_ladder's unit of checking.  "Road" is F160's u axis, as on
  THE_COMB_ON_THE_ROAD.

GATES (exit 0 iff all pass).  One is not exact: L4, an error-model law on an eigensolver,
where no exact route to an eigenvalue exists.  Everything else, L13 included, is decided
exactly in Z[zeta_2n] or symbolically in n.
  [L1]  the characteristic polynomial against the actual matrix, exact, N = 3..9
  [L2]  d_1..d_5 in closed form, symbolic in n, at eta = +1 and eta = -1 separately
  [L2b] the rank-one lemma: V is 2 eta a a^T inside each reflection sector, so d_m carries
        eta^m and the even orders carry no sign at all
  [L3]  the same five as multiplier differences (the cos(2j theta) basis)
  [L4]  the series against 60-digit eigenvalues: residual/u^6 CONVERGES across two decades
  [L5]  the census, from the committed exact layer: 2558 / 2335 / 223
  [L6]  odd n: dX_2 = dX_4 = 0 on every collision pair, hence c_3 = 0 on all 627
  [L7]  odd n, general rung: dX_2j = 0 on every collision whenever gcd(j, n) = 1
  [L8]  NEGATIVE CONTROL: at odd n, dX_2 = 0 EXACTLY on the collisions and nowhere else
  [L9]  even n: the 23 standing, the 11 Theta-mirror, the 12 as parity-matched ROT3, and
        parity-uniformity forced at even n but not at odd n
  [L10] the ROT3 rung lemma, forced direction, plus the BREAK-INPUT at 3|j
  [L11] where it stops: dX_6 = 0 on exactly the 11 mirror pairs, so != 0 on the other 212
  [L12] the second order: c_2 = -((n-3)/n^2) dM_3, the 60, and n = 20's twenty forced
  [L13] the SHARP criterion: c_2 = 0 exactly when SOME minimal tiling of the pair's vanishing
        sum is entirely 3-free (a tiling is not unique, so "no piece is" would be no property
        of the sum); EXACT, by enumerating the vanishing subsets in Z[zeta_2n]; and inside
        this census the 60 are F129 family C, which is not the only 3-free family
  [L14] what the sign of u is: u = +1 and u = -1 are the ring's two parity combs

Run: python simulations/collision_gap_odd_orders.py           (~1 min, measured quiet)
     python simulations/collision_gap_odd_orders.py --fast    (skips L4 and L8's big n; it
                                                               buys little, prefer the full run)
"""
import itertools
import random
import sys

import numpy as np
import sympy as sp
import mpmath as mp

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import f129_level_collision_law as f129            # noqa: E402  the exact ring
import comb_road_f129 as road                      # noqa: E402  the census and the solver

FAST = "--fast" in sys.argv
FIRING = [9, 12, 15, 18, 20, 21, 24, 27, 30]
FAIL = []


def gate(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- the exact ring, cached
_BASIS = {}


def basis(n):
    """The reduction of zeta_2n^e for e = 0..2n-1, from the COMMITTED primitive."""
    if n not in _BASIS:
        _BASIS[n] = np.array([f129.root_sum_vec(n, [e]) for e in range(2 * n)], dtype=np.int64)
    return _BASIS[n]


def twoM(n, tri, j):
    """2 * M_j(tri) in Z[zeta_2n], as an integer vector."""
    b = basis(n)
    return sum(b[(j * k) % (2 * n)] + b[(-j * k) % (2 * n)] for k in tri)


def twoX(n, tri, j):
    """2 * X_2j(tri) in Z[zeta_2n]; eta_k = (-1)^(k+1)."""
    b = basis(n)
    out = np.zeros(b.shape[1], dtype=np.int64)
    for k in tri:
        sgn = 1 if k % 2 else -1
        out += sgn * (b[(2 * j * k) % (2 * n)] + b[(-2 * j * k) % (2 * n)])
    return out


def zero(v):
    return not np.any(v)


# ================================================================== L1
print("L1  the characteristic polynomial against the actual matrix, exact")
xs, us = sp.symbols("x u")
ok = True
for N in range(3, 10):
    H = sp.zeros(N, N)
    for j in range(N - 1):
        H[j, j + 1] = H[j + 1, j] = 1
    H[0, N - 1] = H[N - 1, 0] = us
    ok = ok and sp.expand((xs * sp.eye(N) - H).det()
                          - (sp.chebyshevu(N, xs / 2)
                             - us**2 * sp.chebyshevu(N - 2, xs / 2) - 2 * us)) == 0
gate("det(xI - H(u)) = U_N(x/2) - u^2 U_{N-2}(x/2) - 2u, N = 3..9", ok)

# ================================================================== L2, L3
print()
print("L2  d_1..d_5 in closed form, symbolic in n")
th = sp.Symbol("theta")
nn = sp.Symbol("n", integer=True, positive=True)
eta = sp.Symbol("eta")
s, c = sp.sin(th), sp.cos(th)


def at_root(e):
    e = sp.expand_trig(sp.expand(e))
    return sp.simplify(sp.expand(e.subs({sp.sin(nn * th): 0, sp.cos(nn * th): -eta})))


def d_dx(e):
    return sp.diff(e, th) / (-2 * sp.sin(th))


A_th = sp.sin(nn * th) / s
B_th = sp.expand_trig(sp.sin((nn - 2) * th)) / s
Ad, e = [at_root(A_th)], A_th
for _ in range(5):
    e = d_dx(e)
    Ad.append(at_root(e))
Bd, e = [at_root(B_th)], B_th
for _ in range(3):
    e = d_dx(e)
    Bd.append(at_root(e))

ORDER = 5
dsym = sp.symbols("d1:6")
delta = sum(dsym[i] * us**(i + 1) for i in range(ORDER))
eqn = sp.expand(sum(Ad[m] * delta**m / sp.factorial(m) for m in range(1, ORDER + 1))
                - us**2 * sum(Bd[m] * delta**m / sp.factorial(m) for m in range(4)) - 2 * us)
sol = {}
for m in range(1, ORDER + 1):
    sol[dsym[m - 1]] = sp.simplify(sp.solve(sp.expand(eqn.coeff(us, m).subs(sol)), dsym[m - 1])[0])
D = [sol[dsym[i]] for i in range(ORDER)]


def both_signs(a, b):
    """a == b at eta = +1 and at eta = -1, exactly.  eta is a SIGN, not a symbol."""
    return all(sp.simplify(sp.expand_trig(sp.expand(a.subs(eta, v) - b.subs(eta, v)))) == 0
               for v in (1, -1))


F5 = sp.Rational(4, 15) * (nn - 6) * (nn - 2) / nn**5
CLAIM = [
    eta * 4 * s**2 / nn,
    4 * (nn - 3) * s**2 * c / nn**2,
    sp.Rational(4, 3) * (nn - 2) * (nn - 4) / nn**3 * eta * s**2 * (4 * c**2 - 1),
    None,
    eta * (F5 * sp.Rational(5, 4) * (nn - 1) * sp.cos(2 * th)
           + F5 * (sp.Rational(3, 2) * nn**2 - 8 * nn + 8) * sp.cos(4 * th)
           - F5 * sp.Rational(3, 4) * (2 * nn - 3) * (nn - 3) * sp.cos(6 * th)),
]
gate("d_1 = eta (4/n) sin^2       [= F160 Theorem G, the road's velocity]",
     both_signs(D[0], CLAIM[0]))
gate("d_2 = (4(n-3)/n^2) sin^2 cos                              [eta-free]",
     both_signs(D[1], CLAIM[1]))
gate("d_3 = (4(n-2)(n-4)/(3n^3)) eta sin^2 (4cos^2 - 1)",
     both_signs(D[2], CLAIM[2]))
gate("d_5 = eta F [ (5/4)(n-1) cos2th + (3n^2/2-8n+8) cos4th - (3/4)(2n-3)(n-3) cos6th ]",
     both_signs(D[4], CLAIM[4]))

print()
print("L2b  why only the ODD orders carry the sign: V is rank one inside each parity sector")
ok_v, ok_rank = True, True
for N in range(4, 10):
    m_ = N + 1

    def psi(k, j):
        return sp.sqrt(sp.Rational(2, 1) / m_) * sp.sin(sp.pi * k * (j + 1) / m_)

    for k in range(1, N + 1):
        for l in range(1, N + 1):
            vkl = sp.simplify(psi(k, 0) * psi(l, N - 1) + psi(k, N - 1) * psi(l, 0))
            claim = sp.simplify(psi(k, 0) * psi(l, 0)
                                * ((-1)**(k + 1) + (-1)**(l + 1)))
            ok_v = ok_v and sp.simplify(vkl - claim) == 0
            if (k + l) % 2:                       # opposite reflection parity
                ok_rank = ok_rank and sp.simplify(vkl) == 0
gate("<psi_k|V|psi_l> = psi_k(0) psi_l(0) (eta_k + eta_l), N = 4..9, exact", ok_v)
gate("hence V vanishes between opposite-parity modes: it is 2 eta a a^T in each sector",
     ok_rank)
gate("so the series runs in the scale 2 eta u, and d_m = eta^m times an eta-free function",
     all(sp.simplify(sp.expand_trig(sp.expand(
         D[i].subs(eta, -1) - (-1)**(i + 1) * D[i].subs(eta, 1)))) == 0
         for i in range(ORDER)))

# and the same conclusion for EVERY m, not just the five computed: put x = 2cos(phi),
# phi = theta + eps, into the transcendental form of the characteristic equation.  With
# sin(n theta) = 0 and cos(n theta) = -eta it becomes an equation in which eta appears
# ONLY through s = eta*u, so eps = eps(s) and hence delta = delta(s) to all orders.
eps = sp.Symbol("epsilon")
phi = th + eps
raw = sp.sin(nn * phi) - us**2 * sp.sin((nn - 2) * phi) - 2 * us * sp.sin(phi)
raw = sp.expand_trig(sp.expand(raw)).subs({sp.sin(nn * th): 0, sp.cos(nn * th): -eta})
target = sp.sin(nn * eps) - (eta * us)**2 * sp.sin(nn * eps - 2 * phi) + 2 * (eta * us) * sp.sin(phi)
# NOTE.  Lemma B also quotes the sector resolvent R_k = (n-3)cos(theta)/(2n).  That is NOT
# an independent claim and gets no gate: rank-one second-order perturbation theory gives
# d_2 = (2 eta)^2 a_k^2 R_k, so R_k is d_2 (gated at L2) divided by 4 a_k^2.
gate("and for EVERY m: in the angle equation eta enters only through v = eta*u",
     all(sp.simplify(sp.expand_trig(sp.expand(
         (-v * raw).subs(eta, v) - target.subs(eta, v)))) == 0 for v in (1, -1)),
     "sin(n eps) - v^2 sin(n eps - 2 phi) + 2 v sin phi = 0, no eta in it")

# The parity split (odd orders on even multipliers, even orders on odd ones) is a theorem
# at EVERY m, and it comes out of the angle equation itself rather than out of the comb.  Write
# E(theta, eps, v) = sin(n eps) - v^2 sin(n eps - 2 phi) + 2 v sin phi with phi = theta + eps.
# Two substitutions send E to -E, hence fix its zero set:
#     (theta, eps, v) -> (pi - theta, -eps, -v)   and   (theta, eps) -> (-theta, -eps).
# By uniqueness of the analytic branch eps(v; theta) with eps(0) = 0, the first gives
# eps(-v; pi - theta) = -eps(v; theta), hence delta(-v; pi - theta) = -delta(v; theta), hence
# D_m(pi - theta) = (-1)^(m+1) D_m(theta) as functions of a FREE theta; the second gives
# D_m(-theta) = D_m(theta), so D_m is a cosine polynomial in the first place.  A term
# cos(r theta) in D_m therefore needs r = m+1 mod 2.  No comb points, no chiral K, no
# case split on the parity of N, and no bound on the degree of D_m.
_eps, _v = sp.symbols("epsilon v")
_phi = th + _eps
_E = sp.sin(nn * _eps) - _v**2 * sp.sin(nn * _eps - 2 * _phi) + 2 * _v * sp.sin(_phi)
gate("the angle equation is ODD under (th, eps, v) -> (pi - th, -eps, -v)",
     sp.simplify(sp.expand_trig(sp.expand(
         _E.subs({th: sp.pi - th, _eps: -_eps, _v: -_v}) + _E))) == 0,
     "so D_m(pi - th) = (-1)^(m+1) D_m(th) for EVERY m, not only the five computed")
gate("and ODD under (th, eps) -> (-th, -eps), so D_m is even in th and a cosine polynomial",
     sp.simplify(sp.expand_trig(sp.expand(
         _E.subs({th: -th, _eps: -_eps}) + _E))) == 0)
gate("the five computed D_m obey it, as they must: D_m(pi - th) = (-1)^(m+1) D_m(th)",
     all(sp.simplify(sp.expand_trig(sp.expand(
         (D[i].subs(eta, 1).subs(th, sp.pi - th)
          - (-1)**(i + 2) * D[i].subs(eta, 1)))))== 0 for i in range(ORDER)))

print()
print("L3  the same coefficients as multiplier differences")
gate("d_1 = (2/n) eta (1 - cos 2th)          -> c_1 ~ dX_0 - dX_2",
     both_signs(D[0], sp.Rational(2, 1) / nn * eta * (1 - sp.cos(2 * th))))
gate("d_2 = ((n-3)/n^2) (cos th - cos 3th)   -> c_2 ~ dM_1 - dM_3",
     both_signs(D[1], (nn - 3) / nn**2 * (sp.cos(th) - sp.cos(3 * th))))
gate("d_3 = (2(n-2)(n-4)/(3n^3)) eta (cos 2th - cos 4th)  -> c_3 ~ dX_2 - dX_4",
     both_signs(D[2], 2 * (nn - 2) * (nn - 4) / (3 * nn**3) * eta
                * (sp.cos(2 * th) - sp.cos(4 * th))))
gate("d_5's X_6 coefficient vanishes for no n = 9..60 (all four of its factors are "
     "positive there; the no-X_0 half is pinned by L2's form check)",
     all(sp.simplify((sp.Rational(3, 4) * (2 * nn - 3) * (nn - 3) * F5).subs(nn, v)) != 0
         for v in range(9, 61)))
gate("the four d_5 coefficients sum to zero  [d_5 -> 0 as theta -> 0]",
     sp.simplify(F5 * (sp.Rational(5, 4) * (nn - 1) + sp.Rational(3, 2) * nn**2 - 8 * nn + 8
                       - sp.Rational(3, 4) * (2 * nn - 3) * (nn - 3))) == 0)

# ================================================================== L4
print()
print("L4  the series against 60-digit eigenvalues of the real matrix")
if FAST:
    print("       (skipped under --fast)")
else:
    mp.mp.dps = 60
    for (N, k) in ((14, 4), (14, 9), (19, 7)):
        m = N + 1
        et = 1 if k % 2 else -1
        co = [mp.mpf(str(sp.N(di.subs({nn: m, eta: et}).subs(th, sp.pi * k / m), 50))) for di in D]
        # The error model, not a window: the first omitted term is d_6 u^6, so
        # residual / u^6 must CONVERGE to |d_6| as u falls.  Gate the convergence.
        q = []
        for e in (3, 4, 5, 6):
            uu = mp.mpf(10)**(-e)
            lv = road.road_levels_mp(N, uu)[k]
            ser = 2 * mp.cos(mp.pi * k / m) + sum(co[i] * uu**(i + 1) for i in range(ORDER))
            q.append(abs(lv - ser) / uu**6)
        # The model is q(u) = d_6 + c*u + O(u^2), so successive differences fall by exactly
        # one decade and the RATIO of differences itself approaches 1/10 from its own O(u).
        # Four values of u give two ratios, so the law can be seen converging rather than
        # merely landing inside a window once.
        r1 = (q[2] - q[1]) / (q[1] - q[0])
        r2 = (q[3] - q[2]) / (q[2] - q[1])
        d6 = (q[3] - r2 * q[2]) / (1 - r2)            # Richardson limit
        gate(f"N={N} k={k}: q = residual/u^6 has error O(u), and the difference ratio "
             f"CONVERGES to 1/10",
             abs(r2 - mp.mpf("0.1")) < abs(r1 - mp.mpf("0.1")) < mp.mpf("0.002"),
             f"{mp.nstr(r1, 8)} -> {mp.nstr(r2, 8)}, q -> d_6 = {mp.nstr(d6, 8)}")

# ================================================================== L5
print()
print("L5  the census, from the committed exact layer")
PAIRS = {n: road.exact_collision_pairs(n) for n in FIRING}
INFO = {}
for n in FIRING:
    rows = []
    for tau, sig in PAIRS[n]:
        dX = {j: twoX(n, tau, j) - twoX(n, sig, j) for j in (0, 1, 2, 3)}
        rows.append(dict(
            tau=tau, sig=sig,
            c1zero=zero(dX[0] - dX[1]), c3zero=zero(dX[1] - dX[2]),
            c2zero=zero(twoM(n, tau, 3) - twoM(n, sig, 3)),
            X0=zero(dX[0]), X2=zero(dX[1]), X4=zero(dX[2]), X6=zero(dX[3]),
            mirror=tuple(sorted(n - k for k in tau)) == tuple(sorted(sig))))
    INFO[n] = rows
tot = sum(len(v) for v in INFO.values())
stand = [r for n in FIRING for r in INFO[n] if r["c1zero"]]
gate("2558 collision pairs, 2335 separate at first order, 223 stand",
     (tot, tot - len(stand), len(stand)) == (2558, 2335, 223),
     f"{tot} / {tot - len(stand)} / {len(stand)}")
gate("REGRESSION PIN (not a second route): dX_0 - dX_2 retyped here equals the committed "
     "n_times_D_vec on every pair; both are 2(dX_0 - dX_2) through the same primitive",
     all(r["c1zero"] == (not np.any(road.n_times_D_vec(n, r["tau"], r["sig"])))
         for n in FIRING for r in INFO[n]))

# ================================================================== L6, L7
print()
print("L6  odd n: the Galois kill at the second and third rung")
odd = [r for n in FIRING if n % 2 for r in INFO[n]]
gate("dX_2 = 0 and dX_4 = 0 on every collision pair at odd n",
     all(r["X2"] and r["X4"] for r in odd), f"{len(odd)} pairs")
# c3zero = zero(dX_2 - dX_4) follows from the row above over exact integers; what is
# independent here is the POPULATION, so that is what this row asserts.
gate("and that population is all 627 odd-n collision pairs, so c_3 = 0 on every one",
     len(odd) == 627 and all(r["c3zero"] for r in odd), f"{len(odd)} pairs")
# c1zero == X0 follows from dX_2 = 0 above and would be f(x) == f(x); what can fail is
# the INTEGER form of dX_0 against its cyclotomic vector.
gate("dX_0 equals its integer form 2(o_tau - o_sigma) exactly, on every odd-n pair",
     all(zero(twoX(n, r["tau"], 0) - twoX(n, r["sig"], 0)
              - 4 * (sum(k % 2 for k in r["tau"]) - sum(k % 2 for k in r["sig"]))
              * basis(n)[0])          # twoX carries the factor 2, so the difference is 4*(...)
         for n in FIRING if n % 2 for r in INFO[n]),
     "so c_1 = (4/n)(o_tau - o_sigma), and c_1 = 0 is the equal-odd-count condition")

print()
print("L7  odd n, the general rung: dX_2j = 0 whenever gcd(j, n) = 1")
ok = True
for n in [x for x in FIRING if x % 2]:
    for j in range(1, 9):
        if sp.gcd(j, n) != 1:
            continue
        ok = ok and all(zero(twoX(n, t, j) - twoX(n, s, j)) for t, s in PAIRS[n])
gate("every collision pair, every coprime rung j <= 8, at n = 9, 15, 21, 27", ok)
gate("the criterion is gcd(n + 2j, 2n) = 1, and at odd n that is gcd(j, n) = 1",
     all(sp.gcd(n + 2 * j, 2 * n) == sp.gcd(j, n)
         for n in range(9, 60, 2) for j in range(1, 9)))

# ================================================================== L8
print()
print("L8  the kill fires on collisions and nowhere else.  At odd n this is the CONVERSE")
print("     of Theorem D, proved two paragraphs up, so it is a regression check on the")
print("     implementation rather than a control that could distinguish a false law.")
rng = random.Random(20260902)
for n in ([9, 15] if FAST else [9, 15, 21, 27]):
    tris = [t for t in itertools.combinations(range(1, n), 3) if f129.is_clean(n, t)]
    coll = {(a, b) for a, b in PAIRS[n]} | {(b, a) for a, b in PAIRS[n]}
    allp = list(itertools.combinations(tris, 2))
    sample = allp if len(allp) <= 40000 else rng.sample(allp, 40000)
    bad = [(a, b) for a, b in sample
           if zero(twoX(n, a, 1) - twoX(n, b, 1)) != ((a, b) in coll)]
    how = "all" if len(allp) <= 40000 else "sampled from %d" % len(allp)
    gate(f"n={n}: dX_2 = 0 exactly on the collisions ({len(sample)} pairs, {how})",
         not bad, f"{len(bad)} exceptions" if bad else "")

# ================================================================== L9
print()
print("L9  even n: the 23 that stand, the 11 that are Theta-mirror, the 12 that are ROT3")
ev = [r for n in FIRING if n % 2 == 0 for r in INFO[n] if r["c1zero"]]
gate("23 standing pairs at even n, 11 of them Theta-mirror",
     (len(ev), sum(r["mirror"] for r in ev)) == (23, 11),
     f"{len(ev)} / {sum(r['mirror'] for r in ev)}")
gate("every one has dX_0 = dX_2 = dX_4 = 0",
     all(r["X0"] and r["X2"] and r["X4"] for r in ev))
gate("the 11 Theta-mirror pairs have dX_2j = 0 at every rung j <= 3",
     all(r["X6"] for r in ev if r["mirror"]))


def rot3_cosets(n, tri):
    """The +-label set of tri mod n as a union of two cosets of the order-3 subgroup, or None."""
    if n % 3:
        return None
    step = n // 3
    pm = sorted({e % n for k in tri for e in (k, -k)})
    if len(pm) != 6:
        return None
    cs = sorted({tuple(sorted((a + i * step) % n for i in range(3))) for a in pm})
    return cs if len(cs) == 2 and all(set(x) <= set(pm) for x in cs) else None


twelve = [(n, r) for n in FIRING if n % 2 == 0
          for r in INFO[n] if r["c1zero"] and not r["mirror"]]
gate("the other 12 are non-mirror, and BOTH triples of each are doubled-label ROT3",
     len(twelve) == 12 and all(rot3_cosets(n, t) is not None
                               for n, r in twelve for t in (r["tau"], r["sig"])),
     f"{len(twelve)} pairs")
gate("and both triples of each share a parity class, which is what c_1 = 0 needs",
     all(len({k % 2 for k in r["tau"]}) == 1 and len({k % 2 for k in r["sig"]}) == 1
         and r["tau"][0] % 2 == r["sig"][0] % 2 for n, r in twelve))
_r3 = [(n, t) for n in FIRING if n % 3 == 0
       for t in itertools.combinations(range(1, n), 3)
       if f129.is_clean(n, t) and rot3_cosets(n, t)]
_r3even = [(n, t) for n, t in _r3 if n % 2 == 0]
gate("at EVEN n parity-uniformity is not a second hypothesis: 6|n forces it on every ROT3 triple",
     len(_r3even) == 80 and all(len({k % 2 for k in tt}) == 1 for _, tt in _r3even),
     f"{len(_r3even)} triples at n = 12, 18, 24, 30")
gate("at ODD n it IS a hypothesis: (1,2,4) at n = 9 is clean, ROT3, mixed, and X_2 != 0",
     rot3_cosets(9, (1, 2, 4)) is not None and f129.is_clean(9, (1, 2, 4))
     and len({k % 2 for k in (1, 2, 4)}) == 2 and not zero(twoX(9, (1, 2, 4), 1)))
_both = [(n, r) for n in FIRING if n % 3 == 0 for r in INFO[n]
         if rot3_cosets(n, r["tau"]) and rot3_cosets(n, r["sig"])
         and len({k % 2 for k in r["tau"]}) == 1 and len({k % 2 for k in r["sig"]}) == 1
         and not r["c1zero"]]        # parity-uniform: the hypothesis they refute
_odd = sum(1 for n, _ in _both if n % 2)
_lo = min(_both, key=lambda x: (x[0], x[1]["tau"], x[1]["sig"]))
gate("the ROT3 shape alone does NOT give c_1 = 0: 58 collision pairs are both-ROT3 with c_1 != 0",
     len(_both) == 58,
     f"smallest modulus carrying one: n = {_lo[0]}, {_lo[1]['tau']} ~ {_lo[1]['sig']}, "
     f"c_1 = {sp.Rational(4 * (sum(k % 2 for k in _lo[1]['tau']) - sum(k % 2 for k in _lo[1]['sig'])), _lo[0])}")
gate("and they span BOTH parities of n, so the shape fails on either side",
     (_odd, len(_both) - _odd) == (30, 28), f"{_odd} at odd n, {len(_both) - _odd} at even n")

# ================================================================== L10
print()
print("L10  the ROT3 rung lemma: forced for 3 not dividing j, a coincidence at 3|j")
forced, collapse, broke_zero, broke_nonzero = True, True, False, False
for n in [x for x in FIRING if x % 3 == 0]:
    tris = [t for t in itertools.combinations(range(1, n), 3)
            if f129.is_clean(n, t) and rot3_cosets(n, t)]
    for t in tris:
        et = 1 if t[0] % 2 else -1
        if not all(k % 2 == t[0] % 2 for k in t):
            continue                                  # the lemma is for parity-uniform triples
        a = rot3_cosets(n, t)[0][0]
        for j in range(1, 10):
            v = twoX(n, t, j)
            if j % 3:
                forced = forced and zero(v)
            else:
                b = basis(n)
                pred = et * 3 * (b[(2 * j * a) % (2 * n)] + b[(-2 * j * a) % (2 * n)])
                collapse = collapse and zero(v - pred)
                broke_zero = broke_zero or zero(v)
                broke_nonzero = broke_nonzero or not zero(v)
gate("X_2j = 0 for EVERY parity-uniform ROT3 triple at every 3|n <= 30 and every 3 not | j",
     forced)
gate("at 3|j the coset collapses: X_2j = 3 eta cos(2 pi j a / n), a the coset label", collapse)
_ex24 = (1, 7, 9)
gate("BREAK-INPUT: at 3|j the value is NOT forced, both outcomes occur in the same census",
     broke_zero and broke_nonzero
     and rot3_cosets(24, _ex24) is not None
     and not zero(twoX(24, _ex24, 6)) is False and zero(twoX(24, _ex24, 6))
     and not zero(twoX(24, _ex24, 3)),
     "pinned on n=24 (1,7,9), a=1: zero at j=6 (cos(pi/2)), nonzero at j=3 (cos(pi/4))")

# ================================================================== L11
print()
print("L11  where the ladder stops")
# That every odd firing modulus has 3|n is one line of arithmetic on F129's condition
# (10|n forces n even), NOT a measurement.  An earlier version of this block gated it by
# reading the same condition it was testing, which is f(x) == f(x) and cannot fail.
gate("among the 223 standing pairs, dX_6 = 0 holds exactly for the 11 Theta-mirror ones",
     all(r["X6"] == r["mirror"] for r in stand),
     f"{sum(1 for r in stand if r['X6'])} with dX_6 = 0, "
     f"{sum(1 for r in stand if r['mirror'])} Theta-mirror")
gate("so the other 212 all have dX_6 != 0: the first surviving rung is j = 3, and c_5 != 0",
     sum(1 for r in stand if not r["X6"]) == 212,
     f"{sum(1 for r in stand if not r['X6'])} pairs")

# ================================================================== L12
print()
print("L12  the second order")
c2z = {n: sum(r["c2zero"] for r in INFO[n]) for n in FIRING}
gate("60 collision pairs have c_2 = 0, matching R4's independently measured count",
     sum(c2z.values()) == 60, ", ".join(f"n={k}: {v}" for k, v in c2z.items() if v))
gate("at n = 20 all 20 pairs have c_2 = 0: gcd(3, 20) = 1, so k -> 3k is an automorphism",
     c2z[20] == len(PAIRS[20]) == 20)
gate("no standing pair has c_2 = 0: every one leaves at second order",
     not any(r["c2zero"] for r in stand))

# ================================================================== L13
print()
print("L13  the sharp criterion: a multiplier must be invertible on every PIECE's RATIO-ORDER")
print("       (decided here by enumerating the vanishing subsets exactly in Z[zeta_2n]; the")
print("        committed two-prime decomposition is kept only as a cross-check)")
import f129_family_inventory as inv                # noqa: E402  the committed decomposition

def three_free_cover(n, exps):
    """Is there an exact cover of exps by MINIMAL vanishing pieces, every ratio-order
    coprime to 3?  Exact: a subset vanishes iff its sum is 0 in Z[zeta_2n], decided on
    the same integer vectors the rest of this file uses."""
    b_, m_, k = basis(n), 2 * n, len(exps)
    rows = np.array([b_[e % m_] for e in exps], dtype=np.int64)
    sums = np.zeros((1 << k, rows.shape[1]), dtype=np.int64)
    for i in range(k):                       # subset sums by doubling
        lo = 1 << i
        sums[lo:2 * lo] = sums[:lo] + rows[i]
    van = [msk for msk in range(1, 1 << k) if not sums[msk].any()]
    minimal = []
    for msk in sorted(van, key=lambda x: bin(x).count("1")):
        if not any((s & msk) == s for s in minimal):
            minimal.append(msk)
    ok3 = []
    for msk in minimal:
        idx = [i for i in range(k) if msk >> i & 1]
        g = m_
        for i in idx:
            g = np.gcd(g, (exps[i] - exps[idx[0]]) % m_)
        if (m_ // int(g)) % 3:
            ok3.append(msk)
    full = (1 << k) - 1
    seen = {}

    def cover(rem):
        if rem == 0:
            return True
        if rem in seen:
            return seen[rem]
        lo = rem & -rem
        out = any(cover(rem & ~s) for s in ok3 if (s & lo) and (s & rem) == s)
        seen[rem] = out
        return out

    return cover(full)


three_free_matches, fams, exact_matches = True, {}, True
for n in FIRING:
    mps = (inv._ModP(n, 0), inv._ModP(n, 1))
    for r in INFO[n]:
        exps, dd = inv.twelve_roots(n, r["tau"], r["sig"])
        parts = inv.piece_decomposition(n, exps, mps)
        three_free = all(o % 3 for _, o in parts)
        three_free_matches = three_free_matches and (three_free == r["c2zero"])
        exact_matches = exact_matches and (three_free_cover(n, exps) == r["c2zero"])
        if r["c2zero"]:
            fam = next((k for k, v in inv.FAMILIES.items() if (v[1], v[2]) == (dd, parts)), "?")
            fams[(n, fam)] = fams.get((n, fam), 0) + 1
gate("c_2 = 0 EXACTLY when SOME minimal tiling is entirely 3-free, on all 2558 pairs",
     exact_matches,
     "decided by enumerating the vanishing subsets in Z[zeta_2n], no primes and no greedy choice")
gate("and the committed greedy tiling agrees with that on every pair",
     three_free_matches, "so the inventory's two-prime decomposition costs nothing here")

# BREAK-INPUT for the WORDING: a minimal tiling is not unique, so "no piece has a
# 3-divisible ratio-order" is not a property of the sum.  (15,23,25) ~ (17,19,29) at n = 30
# tiles BOTH as family C's 3-free ((2,2),(5,5),(5,5)) and as ((6,30),(6,30)), which is not.
# The right form, and the one the sufficient direction actually needs, is SOME tiling is 3-free.
_n, _tau, _sig = 30, (15, 23, 25), (17, 19, 29)
_ex, _d = inv.twelve_roots(_n, _tau, _sig)
_m = 2 * _n
_tt = sp.Symbol("t")
_phi = sp.Poly(sp.cyclotomic_poly(_m, _tt), _tt)
_red = [sp.rem(sp.Poly(_tt**(e % _m), _tt), _phi) for e in range(_m)]
def _van(sub):
    acc = sp.Poly(0, _tt)
    for e in sub:
        acc = acc + _red[e % _m]
    return acc.is_zero
def _ord(sub):
    g = _m
    for e in sub:
        g = sp.gcd(g, (e - sub[0]) % _m)
    return _m // int(g)
_all = list(range(len(_ex)))
_v = [s for r in range(2, len(_ex) + 1) for s in itertools.combinations(_all, r)
      if _van([_ex[i] for i in s])]
_min = [s for s in _v if not any(set(o) < set(s) for o in _v)]
def _tile(rem, cur):
    if not rem:
        yield list(cur)
        return
    lo = min(rem)
    for s in _min:
        if lo in s and set(s) <= rem:
            yield from _tile(rem - set(s), cur + [s])
_sigs = {tuple(sorted((len(s), _ord([_ex[i] for i in s])) for s in T))
         for T in _tile(set(_all), [])}
gate("BREAK-INPUT on the wording: a minimal tiling is NOT unique, so the criterion must say "
     "SOME tiling is 3-free and not that no piece is",
     _sigs == {((2, 2), (5, 5), (5, 5)), ((6, 30), (6, 30))},
     f"(15,23,25) ~ (17,19,29) at n = 30 tiles two ways: {sorted(_sigs)}")
gate("the 60 are one F129 family, C (zero mode + two R5 pieces), at its own door 10|n",
     set(f for _, f in fams) == {"C"} and sum(fams.values()) == 60, str(fams))
gate("and their count is that family's committed closed form 2(n-10)",
     all(cnt == 2 * (n - 10) for (n, _), cnt in fams.items()))
_free = sorted(k for k, v in inv.FAMILIES.items() if all(o % 3 for _, o in v[2]))
gate("but C is not the only 3-free family in the committed table: L is too, door 70|n",
     _free == ["C", "L"])
_pred = lambda n: (2 * (n - 10) if n % 10 == 0 else 0) + (20 if n % 70 == 0 else 0)
gate("so the count reads 2(n-10)*[10|n] + 20*[70|n], each term behind its own family's door",
     all(_pred(n) == sum(fams.get((n, f), 0) for f in ("C", "L")) for n in FIRING),
     f"census: {[(n, _pred(n)) for n in FIRING if _pred(n)]}; "
     f"past it: n=70 -> {_pred(70)}, n=210 -> {_pred(210)}, and n=9 -> {_pred(9)} not -2")

# ================================================================== L14
print()
print("L14  what the SIGN of u is: the two fermion-parity combs of the ring")
kk = sp.Symbol("k")
NN_ = sp.Symbol("N", integer=True, positive=True)
G = sp.sin((NN_ + 1) * kk) - sp.Symbol("uu")**2 * sp.sin((NN_ - 1) * kk) \
    - 2 * sp.Symbol("uu") * sp.sin(kk)
gate("u = +1: G = 2 sin k (cos Nk - 1), zeros 2*pi*m/N, the PERIODIC ring comb",
     sp.simplify(sp.expand_trig(G.subs(sp.Symbol("uu"), 1)
                                - 2 * sp.sin(kk) * (sp.cos(NN_ * kk) - 1))) == 0)
gate("u = -1: G = 2 sin k (cos Nk + 1), zeros (2m+1)*pi/N, the ANTI-PERIODIC comb",
     sp.simplify(sp.expand_trig(G.subs(sp.Symbol("uu"), -1)
                                - 2 * sp.sin(kk) * (sp.cos(NN_ * kk) + 1))) == 0)
# the two combs are named in the labels above, so they get asserted and not merely printed:
# G must vanish at every comb point of the open interval, and the counts must be right.
ok_p = ok_a = True
cnt_p = cnt_a = 0
for Nc in range(4, 13):
    for m_ in range(0, Nc):
        kp = sp.Rational(2 * m_, Nc) * sp.pi
        if 0 < kp < sp.pi:
            cnt_p += 1
            ok_p = ok_p and sp.simplify(G.subs({sp.Symbol("uu"): 1, NN_: Nc, kk: kp})) == 0
        ka = sp.Rational(2 * m_ + 1, Nc) * sp.pi   # m_ from 0: k = pi/N is interior
        if 0 < ka < sp.pi:
            cnt_a += 1
            ok_a = ok_a and sp.simplify(G.subs({sp.Symbol("uu"): -1, NN_: Nc, kk: ka})) == 0
gate("and G really vanishes at every interior point of those two combs, N = 4..12",
     ok_p and ok_a, f"{cnt_p} periodic points and {cnt_a} anti-periodic ones, all exact")
print("       NOTE (a reading, not a gate): flipping one bond of a loop flips the product of")
print("       the hoppings around it and no gauge undoes that, so u and -u sit half a flux")
print("       quantum apart.  The two combs above are PROOF_RING_GAP_DOMINANCE's two")
print("       Jordan-Wigner parity sectors, and the ODD part of a level in u is therefore its")
print("       response to WHICH sector, the even part its response to how strongly the loop")
print("       is closed.  That reading is what makes c_3 = 0 say something, for a pair that")
print("       also STANDS at first order: such a collision cannot tell the ring's two parity")
print("       sectors apart below fifth order.  A pair with c_1 != 0 tells them apart at once.")

print()
print(f"{'ALL PASS' if not FAIL else str(len(FAIL)) + ' FAILED: ' + '; '.join(FAIL)}")
sys.exit(1 if FAIL else 0)
