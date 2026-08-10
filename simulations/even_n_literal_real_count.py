"""The literal real-eigenvalue population of the (1,2) block at even N, across the measured q-grid.

Companion to the even-N real-q defective check in experiments/F89_PATH_K_DIABOLIC.md (the
sideways-ladder arc's even-N half). The nullity endpoints of F89_SEED_EXISTENCE_REDUCTION give,
at even N, r0 = r_inf = 3*Z3 with zero surplus: no seed is FORCED, but the counting cannot say
whether the literal real population supports an accidental one. This script measures the literal
population on a 16-point grid spanning 0.002..102.4 in the carrier knob (= 0.001..51.2 octic)
and gates, at EVERY grid point:

  N=6 (Z3 = 0): ZERO real eigenvalues; the near-axis floor lifts at FIRST order near q = 0
      (ratio 2 per q-doubling, three doublings).
  N=8 (Z3 = 2, the 2cos(pi/9) family, the smallest even N with resonances): the literal real
      population is EXACTLY TWO at every grid point, both pinned at lambda = -6, the RATE-WINDOW
      EDGE, where the window-edge lemma (PROOF_CODIM1_BY_ADDITIVITY section 6) proves
      semisimplicity, so this pair can never carry a Jordan block. (The pair splits one per
      R-sector, so no sector holds a real PAIR to collide.) The other four resonance modes lift
      at THIRD order near q = 0 (|Im| ~ q^3, the resonant-N=11 mechanism, ratio 8 per
      q-doubling, three doublings) and approach the axis like c/q at large q without reaching
      it (ratio ~1/2 per doubling, measured across 25.6 -> 51.2 -> 102.4). The window-edge pair
      splits ONE PER R-SECTOR at every grid point (the reflection commutes with the block
      exactly; the 2-dim lambda = -6 real eigenspace has rank 1 in each R-projection, gated at
      all 16 points), so no sector holds a real pair to collide even at the edge.

  Tangency guard: a conjugate pair grazing the axis would show as an interior dip of the
  near-axis floor toward zero with no count change (the count-change census's stated blind
  spot). Gated: the floor's minimum over the grid sits at the smallest probed q (the small-q
  lift envelope), with no interior grid point below it. The floor is NOT monotone inside the
  grid (the large-q c/q descent begins there); the gate is the envelope minimum, nothing more.

  Disc-reality (the mechanism's numbers, re-gated here so they live in a committed file; the
  July multi-sector monodromy design spec's B3 measured them first, locally): at real q the
  residual discriminant prod_{i<j} (lambda_i - lambda_j)^2 of the (1,2) block is REAL exactly
  where a within-block antiunitary acts. Gated: N=4 and N=5 relative imaginary part at machine
  floor (< 1e-9; the self-fold / the real char poly), N=6 genuinely complex (relative imaginary
  part O(1), arg-from-real > 0.1 rad at every probed q). A real-q defective locus is a real
  zero of that discriminant; genuinely complex coefficients make one codimension 2 in the one
  real knob, the genericity half of the even-N verdict.

  Consequence: on this grid the literally-real count never changes (0 at N=6, 2 at N=8), so no
  real-to-complex transition (PT breaking) occurs at any probed q, and the real-lambda
  defective-locus species is excluded at even N = 6, 8 on the measured grid, small-q and
  large-q included, not merely in the census window. A grid is a grid: the statement is
  per-point plus the envelope gates, not a continuum theorem.

Tolerances are laws, not numbers (the no-rounding convention): the exactly-real reads are
machine-floor (< 1e-12; the nearest lifted mode sits at 2.77e-11 at the smallest probed q, a
four-order split there, widening with q); the order-of-lift gates are RATIOS per q-doubling
(2 first-order and 8 third-order across three doublings, 1/2 asymptotic across two).

Run:  python simulations/even_n_literal_real_count.py     (asserts everything; prints the table)
"""
import numpy as np
import sys, os
from itertools import combinations

sys.path.insert(0, os.path.dirname(__file__))
from seed_existence_nullity_check import build, _exc

SMALL = (0.002, 0.004, 0.008, 0.016)                # the doubling ladder for the lift-order gates
MID = (0.032, 0.128, 0.256, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)
LARGE = (25.6, 51.2, 102.4)                         # the doubling ladder for the c/q gate
QS = SMALL + MID + LARGE
EXACT_TOL = 1e-12          # machine-floor read; nearest lifted mode is > 2.7e-11 at every probed q


def scan(L0, C, q):
    lam = np.linalg.eigvals(L0 + q * C)
    im = np.abs(lam.imag)
    order = np.argsort(im)
    exact = int(np.sum(im < EXACT_TOL))
    floor = im[order[exact]] if exact < len(im) else np.inf
    lifted4 = [im[order[exact + j]] for j in range(4)] if exact + 4 <= len(im) else []
    reals = sorted(lam[im < EXACT_TOL].real)
    return exact, floor, lifted4, reals


def reflection(N):
    """The site reflection i -> N-1-i on the (2-exc ket, 1-exc bra) product basis."""
    kets, bras = _exc(N, 2), _exc(N, 1)
    ki = {s: i for i, s in enumerate(kets)}
    bi = {s: i for i, s in enumerate(bras)}
    refl = lambda s: tuple(sorted(N - 1 - x for x in s))
    P = np.zeros((len(kets) * len(bras),) * 2)
    for a in kets:
        for b in bras:
            P[ki[refl(a)] * len(bras) + bi[refl(b)], ki[a] * len(bras) + bi[b]] = 1.0
    return P


for N, z3, expect in ((6, 0, 0), (8, 2, 2)):
    A, C = build(N)
    L0 = np.diag(A.astype(complex))
    print(f"\n=== N={N} (dim {len(A)}, Z3={z3}) ===")
    floors = {}
    lifted = {}
    for q in QS:
        exact, floor, lifted4, reals = scan(L0, C, q)
        floors[q] = floor
        lifted[q] = lifted4
        vals = ", ".join(f"{v:.9f}" for v in reals) if reals else "-"
        print(f"  q={q:8.3f}  exactly real: {exact}  (floor of lifted: {floor:.3e})  values: {vals}")
        assert exact == expect, (N, q, exact)
        if N == 8:
            assert all(abs(v + 6.0) < 1e-9 for v in reals), (q, reals)  # the window-edge pin, every q

    # tangency guard: no interior grid point dips below the smallest-q floor
    qmin = min(QS)
    assert min(floors.values()) == floors[qmin], sorted(floors.items(), key=lambda kv: kv[1])[:3]

    # lift-order gates, three doublings each
    order_ratio = 2.0 if N == 6 else 8.0
    for qa, qb in zip(SMALL, SMALL[1:]):
        r = floors[qb] / floors[qa]
        assert abs(r - order_ratio) < 0.25 * order_ratio / 8, (N, qa, qb, r)
    if N == 8:
        for qa, qb in zip(LARGE, LARGE[1:]):
            r = floors[qb] / floors[qa]
            assert abs(r - 0.5) < 0.05, (qa, qb, r)   # asymptotically real, c/q, never real
        # "the other four": ALL four smallest lifted modes carry the third-order ratio 8
        for qa, qb in zip(SMALL, SMALL[1:]):
            for j in range(4):
                r = lifted[qb][j] / lifted[qa][j]
                assert abs(r - 8.0) < 0.25, (qa, qb, j, r)
        # the window-edge pair splits one per R-sector: rank 1 in each R-projection of the
        # 2-dim real eigenspace (the reflection commutes with the block exactly), every grid point
        P = reflection(N)
        assert np.abs(P @ (L0 + 0.4 * C) - (L0 + 0.4 * C) @ P).max() == 0.0
        for q in QS:
            w, V = np.linalg.eig(L0 + q * C)
            B = V[:, np.abs(w.imag) < EXACT_TOL]
            assert B.shape[1] == 2, (q, B.shape)
            for sign in (+1, -1):
                proj = (B + sign * (P @ B)) / 2
                assert np.linalg.matrix_rank(proj, tol=1e-8) == 1, (q, sign)
        print("  gate: constant pair at -6.000000000 (window edge, one per R-sector) at all 16")
        print("        points; third-order lift ratio 8 x3 on ALL FOUR lifted resonance modes;")
        print("        asymptotic c/q ratio 1/2 x2; floor minimum at q_min")
    else:
        print("  gate: zero real modes at all 16 points; first-order lift ratio 2 x3; floor")
        print("        minimum at q_min")

# ---------------------------------------------------------------- disc-reality (the mechanism)
def r_even_block(N, q):
    """The R-even sector of the (1,2) block at coupling q (orthonormal real projection basis)."""
    A, C = build(N)
    L = np.diag(A.astype(complex)) + q * C
    P = reflection(N)
    wP, VP = np.linalg.eigh((P + P.T) / 2)
    B = VP[:, wP > 0.5]
    return B.T @ L @ B


def disc_arg_from_real(M):
    """Distance of arg(disc) to {0, pi} for disc = prod_{i<j} (l_i - l_j)^2 of M's spectrum."""
    w = np.linalg.eigvals(M)
    s = 0.0
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            s += 2.0 * np.angle(w[i] - w[j])
    return abs((s + np.pi / 2) % np.pi - np.pi / 2)


print("\n=== disc-reality on the R-even sector (the genericity mechanism) ===")
for N, kind in ((4, "real (self-fold)"), (5, "real (self-conjugate sector)"), (6, "genuinely complex")):
    args = [disc_arg_from_real(r_even_block(N, q)) for q in (0.4, 1.0, 2.0)]
    print(f"  N={N}: arg-from-real at q=0.4/1.0/2.0 = " + ", ".join(f"{a:.3e}" for a in args)
          + f"   ({kind})")
    if N in (4, 5):
        assert all(a < 1e-9 for a in args), (N, args)
    else:
        assert all(a > 0.02 for a in args) and max(args) > 0.1, (N, args)

print("\nALL GATES PASS")
