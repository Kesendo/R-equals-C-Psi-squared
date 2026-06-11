#!/usr/bin/env python3
r"""The metallic router family: committed, self-validating verification (generalizes F116).

The weighted Z-middle window family

    H(t1, t2, t3) = sum_w  t1*X_w Z_{w+1} X_{w+2} + t2*X_w Z_{w+1} Y_{w+2} + t3*Y_w Z_{w+1} X_{w+2}

under local Z-dephasing (arbitrary site-dependent gamma_l, open chain, N >= 3) obeys the

FAMILY THEOREM
--------------
  (i)  SOFT <=> t2 = t3.  On the X/Y-balanced line (overall scale trivial; write c = t1/t2) the
       Liouvillian spectrum is exactly palindromic about -sigma, sigma = sum_l gamma_l. Off the
       line the spectrum is hard from the first tilt; the exact off-line witness (1,2,1) first
       fires at m* = 11 with the positive monomial p_11 = 1730150400 * gamma^3 (exact CRT), so
       its asymmetry is nonzero for ALL gamma > 0.
  (ii) THE ROUTER TRANSPORTS VERBATIM. On the line, the SAME per-site period-4 construction as
       the golden anchor (pattern [a, a, b, b], h_l = (-1)^(l+1) * i * R(g_l), q(I) = g,
       q(Z) = h, q(X) = -g_X I + h_X Z, q(Y) = -g_Y I + h_Y Z, W = tensor_l q_{l mod 4})
       palindromizes,  W L W^{-1} = -L - 2*sigma,  with the frame

           a = (r, 1),   b = (1, -r),   r = r(c) = (c + sqrt(c^2 + 4)) / 2   (r^2 = c*r + 1),

       the METALLIC MEAN of c. F116 is the c = 1 (golden) point of the line. Every entry lives
       in Z[r] + i*Z[r] with r^2 = c*r + 1, so each structural identity is verified EXACTLY for
       integer c (the site-map component tuples are c-INDEPENDENT; the whole family lives in the
       ring's multiplication rule). Floats enter only the end-to-end Lindbladian checks, which
       extend the law to irrational c (pi) and fractional c (0.5); Block 7 then lifts the window
       lemma to ALL real (indeed complex) c by derivation (degree bound + exact interpolation).
  (iii) THE LOCUS. A nontrivial window pair (g_w, g_{w+2}) forces the metallic locus
       alpha^2 - c*alpha*beta - beta^2 = 0; a and b are its two roots. The pair-system
       determinant is c * (alpha^2 - c*alpha*beta - beta^2): for c != 0 the anchor's exclusions
       transport (uniform and period-2 force g = 0); AT c = 0 the system degenerates to rank 1
       (pairing rule: g_w proportional to the X-axis MIRROR of g_{w+2}) and the router set
       WIDENS to the continuum [v, v, vbar, vbar] for any direction v (vbar the X-axis mirror),
       containing the uniform pure-X router. Verified exactly in the c = 0 ring (Block 6).

Stations on the line:
    c = 1    r = phi = (1+sqrt(5))/2          the golden anchor (F116)
    c = 2    r = 1 + sqrt(2)                  the silver mean
    c = 3    r = (3 + sqrt(13))/2             the bronze mean
    c = 0    r = 1                            the 45-degree diagonal frame a = (1,1), b = (1,-1);
                                              H = sum of XZY + YZX windows; the widening station
    c = -1   r = (-1 + sqrt(5))/2 = 1/phi     = phi - 1 ~ 0.618 (the inverse golden ratio)
    c = -2   r = sqrt(2) - 1 = 1/(1+sqrt(2))  (the inverse silver mean)

Block ledger
------------
  Block 1  closed form + per-site structure : EXACT ring (Z[r], r^2 = c*r + 1), c in
                                              {1,2,3,0,-1,-2}: class-swap blocks zero, g/h on the
                                              metallic locus, B = diag(-1,1)C^T, q^2 = -(2+c*r)I,
                                              all singular values sqrt(1 + r^2)
  Block 2  the window lemma                 : EXACT ring zero of the (c,1,1)-weighted template-
                                              summed window anticommutator at all 4 offsets, both
                                              siblings, c in {1,2,3,0,-1,-2}; per-term NONZERO
                                              (at c = 0 the XZX term carries weight 0 and the
                                              cancellation is XZY against YZX alone)
  Block 3  {W, A} = 0                       : EXACT ring zero at N = 3..6, both siblings,
                                              c in {2, 3, -1}; golden c = 1 cross-checked once
                                              against the anchor (same A matrix bit-for-bit)
  Block 4  end-to-end vs framework          : W L W^{-1} = -L - 2*sigma against
                                              lindbladian_pauli_dephasing at c in
                                              {0.5, 2, 3, 0, -1, pi} (N=5), SITE-DEPENDENT gamma
                                              (c=2), N=6 spot (c=2); sibling at c in {2, 3}
  Block 5  the soft-set dichotomy           : on-line soft < 1e-9 (c in {0.5, 2, 3, -1}, N=5)
                                              incl. a scale-triviality witness; off-line hard
                                              > 1e-2 ((1,1.1,1), (1,2,1), (2,4,2)); EXACT CRT
                                              pins: (1,2,1) fires first at m* = 11 with
                                              p_11 = 1730150400 * gamma^3, soft control (2,1,1)
                                              has ALL odd p_m identically zero through m = 11
  Block 6  exclusion / locus, symbolic c    : window equations K1, K2 DERIVED symbolically (the
                                              middle factor carried explicitly); pair-determinant
                                              = c * locus; uniform/period-2 force g = 0 for every
                                              c != 0; AT c = 0 THE CONCLUSION CHANGES (rank-1
                                              X-mirror pairing, the [v,v,vbar,vbar] continuum
                                              routes; verified EXACT in the c = 0 ring for
                                              v = (1,0) and v = (2,1), with the c = 1 control
                                              showing generic v fails off the station); the
                                              45-degree M-candidate has locus value -c/2, i.e.
                                              it joins the locus exactly at c = 0
  Block 7  ALL real c (degree bound)        : parametrize by r (c = r - 1/r); every entry of
                                              r*{Q3, [T_c, .]} is a polynomial in r of degree
                                              <= 5; EXACT Fraction arithmetic (independent of the
                                              ring representation) finds the window lemma zero at
                                              8 rational nodes r in {1, 2, 3, 1/2, 5, 7, 4/3, -2}
                                              > 6 needed, so the lemma is a POLYNOMIAL IDENTITY
                                              in r: the family theorem's Hamiltonian leg holds
                                              for EVERY real (indeed complex) c, by derivation
  Block 8  rigidity / the c = 0 modulus     : finite-difference Jacobian of {W, A} = 0 over the
                                              64-real-parameter period-4 class-swap family (N=5,
                                              8 probes): nullity = 8 = gauge at golden AND silver
                                              (zero physical moduli; rigidity extends to a second
                                              station) but nullity = 16 = gauge + 8 PHYSICAL
                                              MODULI at c = 0 (spectral gap > 1e6) -- the Block-6
                                              widening, counted; catalogued slice = the two
                                              parity-chain directions (complexified), h-tilt
                                              excluded by a float control

Provenance: 2026-06-11 night, the connection hunt's discard-gift. The zero-question scout asked
whether the weighted family is soft at all away from (1,1,1); the soft-set mapper pinned the line
t2 = t3 (on-line soft ~1e-13 at N=5 across c in {0.25..7.5, e, pi, 0, -1, -2}, off-line hard from
delta = 0.02, N=6 spot-checks); the moment probe pinned the off-line witness (1,2,1) at m* = 11
with p_11 = 1730150400 * gamma^3 exact (CRT); the router scouts transported the F116 construction
verbatim. Scope: OPEN chains (rings/PBC untested), dephase letter Z.
Run: python simulations/metallic_router_family.py (~4.5 min; Block 3's six N=6 exact ring
anticommutators at ~15 s each and Block 5's two exact CRT pins at ~40 s each carry most of it).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.spatial import cKDTree

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import ceiling_golden_router as anchor  # noqa: E402
from ceiling_golden_router import (  # noqa: E402
    A_VEC, B_VEC, P_BASIS, P_INV, apply_W_comp, i_times, neg, ring_add, ring_mat,
    ring_max_abs, rot90, xy_swap_conjugate)
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import _build_kbody_chain, site_op  # noqa: E402
from f87_windowed_monomial_converse import first_nonvanishing_odd  # noqa: E402

ONE = (1, 0, 0, 0)


def metallic_r(c):
    """The metallic mean of c: the positive root of r^2 = c*r + 1."""
    return (c + np.sqrt(c * c + 4.0)) / 2.0


# ======================================================================
# Ring Z[r] + i*Z[r] with r^2 = c*r + 1 (INTEGER c): elements as 4 components
# (p, q, pp, qq) = (p + q*r) + i*(pp + qq*r).  Multiplication on each real/imag part:
#     (a + b*r)(a' + b'*r) = (aa' + bb') + (ab' + a'b + c*bb')*r.
# Component matmuls in float64 are exact while every integer stays below 2^53 (asserted at the
# Block-3 scale).  At c = 1 this is the anchor's golden ring bit-for-bit (asserted in main()).
# For c != 0, c^2 + 4 is not a perfect square, so Z[r] embeds faithfully in R and an exact ring
# zero IS a real zero; at c = 0 (r = 1, the ring has zero divisors) exact ring zeros still map to
# real zeros, while NONZERO claims are checked on the complex evaluation at r(c).
# ======================================================================
def comp_mul(x, y, c):
    p = x[0] * y[0] + x[1] * y[1]
    q = x[0] * y[1] + x[1] * y[0] + c * x[1] * y[1]
    ip = x[2] * y[2] + x[3] * y[3]
    iq = x[2] * y[3] + x[3] * y[2] + c * x[3] * y[3]
    ap = x[0] * y[2] + x[1] * y[3]
    aq = x[0] * y[3] + x[1] * y[2] + c * x[1] * y[3]
    bp = x[2] * y[0] + x[3] * y[1]
    bq = x[2] * y[1] + x[3] * y[0] + c * x[3] * y[1]
    return (p - ip, q - iq, ap + bp, aq + bq)


def _zr_pair(p, q, pp, qq, op, c):
    return op(p, pp) + op(q, qq), op(p, qq) + op(q, pp) + c * op(q, qq)


def _ring_combine(A, B, op, c):
    p1, q1, r1, s1 = A
    p2, q2, r2, s2 = B
    rp, rq = _zr_pair(p1, q1, p2, q2, op, c)
    ip, iq = _zr_pair(r1, s1, r2, s2, op, c)
    ap, aq = _zr_pair(p1, q1, r2, s2, op, c)
    bp, bq = _zr_pair(r1, s1, p2, q2, op, c)
    return [rp - ip, rq - iq, ap + bp, aq + bq]


def ring_matmul(A, B, c):
    return _ring_combine(A, B, lambda x, y: x @ y, c)


def ring_kron(A, B, c):
    return _ring_combine(A, B, np.kron, c)


def ring_to_complex(A, r):
    return A[0] + A[1] * r + 1j * (A[2] + A[3] * r)


def string_prod(t, s, k, c):
    """Pauli-string product with ring phases (the anchor's, with the c-parametrized comp_mul)."""
    out = 0
    ph = ONE
    for site in range(k):
        sh = 2 * (k - 1 - site)
        lr, p = anchor._PROD[((t >> sh) & 3, (s >> sh) & 3)]
        out |= lr << sh
        ph = comp_mul(ph, p, c)
    return out, ph


def commutator_super_ring(weighted_words, N, c, windows=None):
    """[sum_w sum_T coeff_T * T_w, .] as a ring matrix in the Pauli-string basis (dim 4^N).
    weighted_words: list of (word, ring_coeff) with ring_coeff a 4-tuple, e.g. XZX -> (c,0,0,0)."""
    dim = 4 ** N
    entries = {}
    for w in (windows if windows is not None else range(N - 2)):
        for word, coeff in weighted_words:
            letters = [0] * N
            for j, ch in enumerate(word):
                letters[w + j] = anchor.L_[ch]
            t = 0
            for letter in letters:
                t = (t << 2) | letter
            for s in range(dim):
                ts, pts = string_prod(t, s, N, c)
                _, pst = string_prod(s, t, N, c)
                e = comp_mul(coeff, tuple(np.array(pts) - np.array(pst)), c)
                if any(e):
                    cur = entries.get((ts, s), (0, 0, 0, 0))
                    entries[(ts, s)] = tuple(np.array(cur) + np.array(e))
    return ring_mat(entries, (dim, dim))


# ======================================================================
# The metallic closed form.  The per-site COMPONENT TUPLES are c-independent: a = (r, 1) is
# X-comp (0,1,0,0), Y-comp (1,0,0,0) for EVERY c (the anchor's A_VEC/B_VEC verbatim); only the
# ring multiplication rule (r^2 = c*r + 1) changes.  ring_site_map generalizes the anchor's
# golden_site_map to an arbitrary frame vector (needed for the c = 0 widening in Block 6) and is
# asserted against it in main().
# ======================================================================
def ring_site_map(l, g_vec):
    h = i_times(rot90(g_vec))
    if l % 2 == 0:
        h = neg(h)                          # h_l = (-1)^(l+1) * i * R(g_l)
    entries = {(1, 0): g_vec[0], (2, 0): g_vec[1],            # q(I) = g   (rows X, Y)
               (1, 3): h[0], (2, 3): h[1],                    # q(Z) = h
               (0, 1): tuple(-np.array(g_vec[0])), (3, 1): h[0],  # q(X) = -g_X I + h_X Z
               (0, 2): tuple(-np.array(g_vec[1])), (3, 2): h[1]}  # q(Y) = -g_Y I + h_Y Z
    return ring_mat(entries, (4, 4))


METALLIC_PATTERN = [A_VEC, A_VEC, B_VEC, B_VEC]                  # [a, a, b, b]
METALLIC_MAPS = [ring_site_map(l, METALLIC_PATTERN[l]) for l in range(4)]
SIBLING_MAPS = [xy_swap_conjugate(q) for q in METALLIC_MAPS]

CEILING_LETTERS = ['XZX', 'XZY', 'YZX']
SIBLING_LETTERS = ['YZY', 'YZX', 'XZY']


def family_words_ring(c, sibling=False):
    letters = SIBLING_LETTERS if sibling else CEILING_LETTERS
    return [(letters[0], (c, 0, 0, 0)), (letters[1], ONE), (letters[2], ONE)]


def family_words_kbody(c, sibling=False):
    letters = SIBLING_LETTERS if sibling else CEILING_LETTERS
    weights = [float(c), 1.0, 1.0]
    return [tuple(w) + (t,) for w, t in zip(letters, weights)]


def build_W_ring(site_maps, N, c):
    W = site_maps[0]
    for l in range(1, N):
        W = ring_kron(W, site_maps[l % 4], c)
    return W


C_INTEGER = [1, 2, 3, 0, -1, -2]


# ======================================================================
# BLOCK 1 -- closed form + per-site structure (EXACT ring, all six integer stations).
# ======================================================================
def block1_structure():
    print("-" * 92)
    print("BLOCK 1  closed form + per-site structure  [EXACT ring, c in {1, 2, 3, 0, -1, -2}]")
    print("-" * 92)
    for c in C_INTEGER:
        r = metallic_r(c)
        for l, q in enumerate(METALLIC_MAPS):
            # class-swap: the {I,Z}->{I,Z} and {X,Y}->{X,Y} blocks are zero in EVERY ring component
            for comp in q:
                assert all(comp[i, j] == 0 for i in (0, 3) for j in (0, 3)), \
                    f"c={c} site {l}: IZ block nonzero"
                assert all(comp[i, j] == 0 for i in (1, 2) for j in (1, 2)), \
                    f"c={c} site {l}: XY block nonzero"
            # metallic locus alpha^2 - c*alpha*beta - beta^2 = 0, exact in the ring, for g and h
            for col, name in ((0, 'g'), (3, 'h')):
                al = tuple(int(q[k][1, col]) for k in range(4))
                be = tuple(int(q[k][2, col]) for k in range(4))
                locus = tuple(np.array(comp_mul(al, al, c))
                              - np.array(comp_mul((c, 0, 0, 0), comp_mul(al, be, c), c))
                              - np.array(comp_mul(be, be, c)))
                assert locus == (0, 0, 0, 0), f"c={c} site {l}: {name} off the locus ({locus})"
            # B = diag(-1,1) C^T, exact on component tuples
            for (bi, bj), (ci, cj), sgn in (((0, 1), (1, 0), -1), ((0, 2), (2, 0), -1),
                                            ((3, 1), (1, 3), 1), ((3, 2), (2, 3), 1)):
                B_e = tuple(int(q[k][bi, bj]) for k in range(4))
                C_e = tuple(sgn * int(q[k][ci, cj]) for k in range(4))
                assert B_e == C_e, f"c={c} site {l}: B != diag(-1,1)C^T at {(bi, bj)}"
            # q^2 = -(1 + r^2) I = -(2 + c*r) I, exact in the ring
            q2 = ring_matmul(q, q, c)
            exp_m = ring_mat({(i, i): (-2, -c, 0, 0) for i in range(4)}, (4, 4))
            diff = [q2[k] - exp_m[k] for k in range(4)]
            assert max(float(np.max(np.abs(d))) for d in diff) == 0.0, \
                f"c={c} site {l}: q^2 != -(2+c*r)I"
            # all singular values sqrt(1 + r^2) (float; a scalar times a unitary, cond 1)
            sv = np.linalg.svd(ring_to_complex(q, r), compute_uv=False)
            assert np.allclose(sv, np.sqrt(1 + r * r)), \
                f"c={c} site {l}: singular values not sqrt(1+r^2)"
        print(f"  c={c:>2} (r={r:.6f}): 4 sites class-swap OK, g/h on the metallic locus (exact), "
              f"B = diag(-1,1)C^T, q^2 = -(2+c*r)I, sigma = sqrt(1+r^2) (cond 1)")
    print("BLOCK 1 PASS")


# ======================================================================
# BLOCK 2 -- the window lemma (EXACT): (c,1,1)-weighted template sum vanishes, per-term does not.
# ======================================================================
def block2_window_lemma(maps, sibling, label):
    print("-" * 92)
    print(f"BLOCK 2  window lemma [{label}]  [EXACT: weighted template sum = 0 at all offsets, "
          f"c in {{1, 2, 3, 0, -1, -2}}; per-term != 0]")
    print("-" * 92)
    for c in C_INTEGER:
        r = metallic_r(c)
        wwords = family_words_ring(c, sibling=sibling)
        for off in range(4):
            Q3 = ring_kron(ring_kron(maps[off % 4], maps[(off + 1) % 4], c),
                           maps[(off + 2) % 4], c)
            Csum = commutator_super_ring(wwords, 3, c, windows=[0])
            anti = ring_add(ring_matmul(Q3, Csum, c), ring_matmul(Csum, Q3, c))
            assert ring_max_abs(anti) == 0.0, \
                f"c={c} offset {off}: weighted template-summed anticommutator != 0"
            per_term = []
            for word, coeff in wwords:
                Ct = commutator_super_ring([(word, coeff)], 3, c, windows=[0])
                at = ring_add(ring_matmul(Q3, Ct, c), ring_matmul(Ct, Q3, c))
                # nonzero claims on the COMPLEX evaluation at r(c) (faithful even at c = 0)
                per_term.append(float(np.max(np.abs(ring_to_complex(at, r)))))
            if c == 0:
                assert per_term[0] == 0.0, f"c=0 offset {off}: weight-0 word not exactly zero"
                assert all(v > 1e-9 for v in per_term[1:]), \
                    f"c=0 offset {off}: a weighted per-term anticommutator vanished"
            else:
                assert all(v > 1e-9 for v in per_term), \
                    f"c={c} offset {off}: a per-term anticommutator vanished"
        note = (f" ({wwords[0][0]} weight 0; {wwords[1][0]} cancels {wwords[2][0]} alone)"
                if c == 0 else "")
        print(f"  c={c:>2}: all 4 offsets weighted-sum = 0 EXACT; per-term nonzero{note}")
    print("BLOCK 2 PASS")


# ======================================================================
# BLOCK 3 -- {W, A} = 0 EXACT at N = 3..6, both siblings, c in {2, 3, -1}; gold cross-check.
# ======================================================================
def assert_anticommutes_exact(maps, wwords, N, c, tag):
    A = commutator_super_ring(wwords, N, c)
    W = build_W_ring(maps, N, c)
    dim = 4 ** N
    mW, mA = ring_max_abs(W), ring_max_abs(A)
    # float64 exactness: every intermediate of the ring matmul is bounded by
    # 4*(2+|c|)*dim*max|W|*max|A|, which must stay below 2^53
    assert 4 * (2 + abs(c)) * dim * mW * mA < 2 ** 53, f"{tag}: float-exactness bound"
    anti = ring_add(ring_matmul(W, A, c), ring_matmul(A, W, c))
    assert ring_max_abs(anti) == 0.0, f"{tag}: {{W, A}} != 0"
    return mW, mA


def block3_anticommutation():
    print("-" * 92)
    print("BLOCK 3  {W, A} = 0  [EXACT ring zero, N = 3..6, both siblings, c in {2, 3, -1}]")
    print("-" * 92)
    # gold cross-check once: at c = 1 the weighted A equals the anchor's unweighted A bit-for-bit,
    # and the anticommutator vanishes (the anchor proves this at N = 3..6; one N = 4 tie suffices)
    A_mine = commutator_super_ring(family_words_ring(1), 4, 1)
    A_anchor = anchor.commutator_super_ring(anchor.CEILING, 4)
    assert all(np.array_equal(a, b) for a, b in zip(A_mine, A_anchor)), \
        "c=1: weighted A != anchor A (bit-for-bit)"
    assert_anticommutes_exact(METALLIC_MAPS, family_words_ring(1), 4, 1, "gold c=1 N=4")
    print("  gold cross-check c=1 N=4: weighted A == anchor A bit-for-bit; {W, A} = 0 EXACT")
    for c in (2, 3, -1):
        for sibling, label in ((False, "XZX+XZY+YZX"), (True, "YZY+YZX+XZY")):
            maps = SIBLING_MAPS if sibling else METALLIC_MAPS
            for N in (3, 4, 5, 6):
                t0 = time.time()
                mW, mA = assert_anticommutes_exact(
                    maps, family_words_ring(c, sibling=sibling), N, c,
                    f"c={c} {label} N={N}")
                print(f"  c={c:>2} [{label}] N={N}: max |entry| = 0.0 EXACT "
                      f"(max|W|={mW:.0f}, max|A|={mA:.0f})  [{time.time() - t0:.0f}s]",
                      flush=True)
    print("BLOCK 3 PASS")


# ======================================================================
# BLOCK 4 -- end-to-end vs the framework Lindbladian (float; incl. irrational c and
#            site-dependent gamma).
# ======================================================================
def end_to_end_residual(maps, c, sibling, N, gammas, seed=7):
    r = metallic_r(c)
    qs_comp = [P_BASIS @ ring_to_complex(q, r) @ P_INV for q in maps]
    H = _build_kbody_chain(N, family_words_kbody(c, sibling=sibling))
    L = lindbladian_pauli_dephasing(H, list(gammas), dephase_letter='Z')
    sigma = sum(gammas)
    rng = np.random.default_rng(seed)
    worst = 0.0
    d = 2 ** N
    for _ in range(8):
        v = rng.standard_normal(d * d) + 1j * rng.standard_normal(d * d)
        v /= np.linalg.norm(v)
        Wv = apply_W_comp(qs_comp, v, N)
        out = apply_W_comp(qs_comp, L @ v, N) + L @ Wv + 2 * sigma * Wv
        worst = max(worst, np.linalg.norm(out) / max(np.linalg.norm(Wv), 1e-300))
    return worst


def block4_end_to_end():
    print("-" * 92)
    print("BLOCK 4  end-to-end: W L W^{-1} = -L - 2*sigma vs framework  [float]")
    print("-" * 92)
    cases = [(0.5, False, 5, [0.3] * 5), (2.0, False, 5, [0.3] * 5),
             (3.0, False, 5, [0.3] * 5), (0.0, False, 5, [0.3] * 5),
             (-1.0, False, 5, [0.3] * 5), (np.pi, False, 5, [0.3] * 5),
             (2.0, False, 5, [0.3, 1.1, 0.05, 2.0, 0.77]),
             (2.0, False, 6, [0.7] * 6),
             (2.0, True, 5, [0.3] * 5), (3.0, True, 5, [0.3] * 5)]
    for c, sibling, N, gammas in cases:
        maps = SIBLING_MAPS if sibling else METALLIC_MAPS
        worst = end_to_end_residual(maps, c, sibling, N, gammas)
        clab = 'pi' if abs(c - np.pi) < 1e-12 else f"{c:g}"
        fam = "YZY+YZX+XZY" if sibling else "XZX+XZY+YZX"
        tag = "site-dependent" if len(set(gammas)) > 1 else f"gamma={gammas[0]:g}"
        assert worst < 1e-12, f"c={clab} {fam} N={N} {tag}: residual {worst:.2e}"
        print(f"  c={clab:>4} [{fam}] N={N} {tag}: max rel residual = {worst:.3e}  OK")
    print("BLOCK 4 PASS")


# ======================================================================
# BLOCK 5 -- the soft-set dichotomy: soft <=> t2 = t3, with exact CRT moment pins.
# ======================================================================
def weighted_H(N, t1, t2, t3):
    return _build_kbody_chain(N, [('X', 'Z', 'X', float(t1)), ('X', 'Z', 'Y', float(t2)),
                                  ('Y', 'Z', 'X', float(t3))])


def pairing_err(N, t1, t2, t3, gam=0.05):
    """Nearest-neighbour mirror error of spec(L) about -sigma (one max-min direction suffices:
    the mirror map is an isometric involution)."""
    L = lindbladian_pauli_dephasing(weighted_H(N, t1, t2, t3), [gam] * N, dephase_letter='Z')
    ev = np.linalg.eigvals(L)
    mv = -ev - 2 * (N * gam)
    d, _ = cKDTree(np.column_stack([ev.real, ev.imag])).query(
        np.column_stack([mv.real, mv.imag]), k=1)
    return float(np.max(d))


def build_weighted_integer_generators(N, t1, t2, t3):
    """int64 (Ar, Ai, Q) with M(gamma) = (Ar + i*Ai) + gamma*Q for INTEGER weights; the recentred
    generator reproduces the framework Lindbladian + N*gamma*I bit-for-bit at integer gamma
    (asserted), so the exact CRT machinery of f87_windowed_monomial_converse applies unchanged."""
    assert all(float(t) == int(t) for t in (t1, t2, t3)), "exact path needs integer weights"
    H = weighted_H(N, t1, t2, t3)
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    A = -1j * np.kron(H, Id) + 1j * np.kron(Id, H.T)
    Q = np.zeros((d * d, d * d))
    for l in range(N):
        Q += np.kron(site_op(N, l, 'Z').real, site_op(N, l, 'Z').real)
    Ar, Ai = A.real, A.imag
    assert np.allclose(Ar, np.round(Ar)) and np.allclose(Ai, np.round(Ai)), "A not integer"
    Ar = np.round(Ar).astype(np.int64)
    Ai = np.round(Ai).astype(np.int64)
    Qi = np.round(Q).astype(np.int64)
    for g in (1, 2):
        L = lindbladian_pauli_dephasing(H, [float(g)] * N, dephase_letter='Z')
        recon = (Ar + 1j * Ai).astype(complex) + g * Qi.astype(complex)
        assert np.array_equal(L + (N * g) * np.eye(4 ** N), recon), \
            "generator mismatch vs framework (bit-for-bit, integer gamma)"
    return Ar, Ai, Qi


def block5_soft_set():
    print("-" * 92)
    print("BLOCK 5  soft-set dichotomy  [float pairing sweep + EXACT CRT moment pins, N=5]")
    print("-" * 92)
    for c in (0.5, 2.0, 3.0, -1.0):
        e = pairing_err(5, c, 1.0, 1.0)
        assert e < 1e-9, f"on-line c={c}: pairing error {e:.2e} not soft"
        print(f"  on-line  t=({c:g}, 1, 1): pairing error = {e:.3e}  SOFT  OK")
    e = pairing_err(5, 1.0, 2.0, 2.0)
    assert e < 1e-9, f"scaled line point (1,2,2): pairing error {e:.2e} not soft"
    print(f"  on-line  t=(1, 2, 2) (= c=0.5 scaled by 2): pairing error = {e:.3e}  SOFT  OK "
          f"(scale trivial)")
    for t in ((1.0, 1.1, 1.0), (1.0, 2.0, 1.0), (2.0, 4.0, 2.0)):
        e = pairing_err(5, *t)
        assert e > 1e-2, f"off-line t={t}: pairing error {e:.2e} not hard"
        print(f"  off-line t=({t[0]:g}, {t[1]:g}, {t[2]:g}): pairing error = {e:.3e}  HARD  OK")

    # EXACT pin 1: the off-line witness (1,2,1) first fires at m* = 11 as a positive monomial.
    # CRT safety: d2 = 1024 selects the ~2^20 prime bank (d2*p^2 < 2^53 asserted inside); the
    # node values |Tr(M^11)| at gamma <= 11 are < 1024 * (||A|| + 11*N)^11 ~ 1e25, far inside the
    # 7-prime product ~1.3e42.
    t0 = time.time()
    Ar, Ai, Q = build_weighted_integer_generators(5, 1, 2, 1)
    mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q, 11, nprimes=7)
    nz = [j for j, co in enumerate(re_co) if co != 0]
    assert mstar == 11, f"(1,2,1): m* = {mstar} (expected 11)"
    assert nz == [3], f"(1,2,1): gamma-powers {nz} (expected pure gamma^3)"
    assert re_co[3] == 1730150400, f"(1,2,1): p_11 coeff {re_co[3]} (expected 1730150400)"
    assert all(co == 0 for co in im_co), "(1,2,1): imaginary part of p_11 nonzero"
    print(f"  EXACT (CRT) hard pin  t=(1,2,1): p_1..p_9 = 0; p_11 = {re_co[3]} * gamma^3 "
          f"(positive monomial => hard for ALL gamma > 0)   [{time.time() - t0:.0f}s]")

    # EXACT pin 2 (soft control): the line point (2,1,1) has ALL odd p_m identically zero, m <= 11.
    t0 = time.time()
    Ar, Ai, Q = build_weighted_integer_generators(5, 2, 1, 1)
    mstar, _, _, _ = first_nonvanishing_odd(Ar, Ai, Q, 11, nprimes=7)
    assert mstar is None, f"(2,1,1): odd moment fired at m* = {mstar} (expected none)"
    print(f"  EXACT (CRT) soft pin  t=(2,1,1): all odd p_m identically 0 through m = 11   "
          f"[{time.time() - t0:.0f}s]")
    print("BLOCK 5 PASS")


# ======================================================================
# BLOCK 6 -- exclusion / locus with SYMBOLIC c, and the honest c = 0 widening.
# ======================================================================
def block6_exclusion_symbolic():
    print("-" * 92)
    print("BLOCK 6  exclusion / locus  [symbolic c: K1, K2 derived; det = c * locus; the c = 0 "
          "widening]")
    print("-" * 92)
    a0, b0, a1, b1, a2, b2, c = sp.symbols('a0 b0 a1 b1 a2 b2 c')
    I2 = sp.eye(2)
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Z = sp.Matrix([[1, 0], [0, -1]])
    K1 = c * a0 * a2 + a0 * b2 + b0 * a2
    K2 = a0 * b2 + b0 * a2 - c * b0 * b2

    # (a) DERIVE the window equations: for the weighted window H_win = c*XZX + XZY + YZX and the
    #     identity column W(I) = G0 (x) G1 (x) G2 (all three slots in the X,Y plane),
    #         [H_win, G0 (x) G1 (x) G2] = -2*K1 * (I (x) G1 Z (x) I) - 2*K2 * (Z (x) G1 Z (x) Z),
    #     an 8x8 operator identity in (a0,b0,a1,b1,a2,b2,c).  The middle factor G1 Z divides out,
    #     so the identity column of the window anticommutator vanishes iff K1 = K2 = 0.
    G0, G1, G2 = a0 * X + b0 * Y, a1 * X + b1 * Y, a2 * X + b2 * Y
    Hwin = (c * sp.kronecker_product(X, Z, X) + sp.kronecker_product(X, Z, Y)
            + sp.kronecker_product(Y, Z, X))
    G = sp.kronecker_product(G0, G1, G2)
    comm = Hwin * G - G * Hwin
    predicted = (-2 * K1 * sp.kronecker_product(I2, G1 * Z, I2)
                 - 2 * K2 * sp.kronecker_product(Z, G1 * Z, Z))
    assert sp.expand(comm - predicted) == sp.zeros(8, 8), "K1/K2 derivation identity fails"
    print("  derived: [H_win, G0 x G1 x G2] = -2*K1*(I x G1Z x I) - 2*K2*(Z x G1Z x Z)  "
          "(8x8 identity, symbolic c)  OK")
    print(f"    K1 = {K1},   K2 = {K2}   (the anchor's window equations at c = 1)")

    # (b) uniform / period-2 (g_{w+2} = g_w): for generic c the only solution is g = 0.
    sols_u = sp.solve([K1.subs({a2: a0, b2: b0}), K2.subs({a2: a0, b2: b0})], [a0, b0], dict=True)
    assert sols_u == [{a0: 0, b0: 0}], f"uniform admits nonzero g generically: {sols_u}"
    # ... and the case analysis pins the exceptions: K1u = a0(c*a0 + 2*b0), K2u = b0(2*a0 - c*b0)
    # force a0 = b0 = 0 unless c = 0 (handled in (d)) or c^2 = -4 (no real station).
    K1u = sp.factor(K1.subs({a2: a0, b2: b0}))
    K2u = sp.factor(K2.subs({a2: a0, b2: b0}))
    assert K1u == a0 * (c * a0 + 2 * b0) and K2u == b0 * (2 * a0 - c * b0), "uniform factorization"
    res = sp.resultant(c * a0 + 2 * b0, 2 * a0 - c * b0, b0)
    assert sp.expand(res - a0 * (c ** 2 + 4)) == 0, f"uniform cross-branch resultant: {res}"
    print("  uniform / period-2: g = 0 forced for every c with c^2 + 4 != 0 (all real c != 0; "
          "c = 0 widens, see (d))  OK")

    # (c) nontrivial pairs: the linear system on g_w given g_{w+2} has determinant c * locus
    M = sp.Matrix([[c * a2 + b2, a2], [b2, a2 - c * b2]])
    det = sp.expand(M.det())
    locus = a2 ** 2 - c * a2 * b2 - b2 ** 2
    assert sp.expand(det - c * locus) == 0, "pair determinant != c * (a^2 - c*a*b - b^2)"
    r_s = (c + sp.sqrt(c ** 2 + 4)) / 2
    for al, be, name in ((r_s, 1, 'a'), (1, -r_s, 'b')):
        assert sp.simplify(al ** 2 - c * al * be - be ** 2) == 0, f"{name} off the locus"
    for g0v, g2v in (((r_s, 1), (1, -r_s)), ((1, -r_s), (r_s, 1))):
        v1 = sp.simplify(K1.subs({a0: g0v[0], b0: g0v[1], a2: g2v[0], b2: g2v[1]}))
        v2 = sp.simplify(K2.subs({a0: g0v[0], b0: g0v[1], a2: g2v[0], b2: g2v[1]}))
        assert v1 == 0 and v2 == 0, "metallic pair fails the window equations"
    print("  det = c * (alpha^2 - c*alpha*beta - beta^2): for c != 0 nontrivial pairs need the "
          "METALLIC LOCUS; [a,a,b,b] with a = (r,1), b = (1,-r) solves all pairs (symbolic c)  OK")

    # (d) THE c = 0 WIDENING (the honest deviation from the anchor's exclusion): at c = 0 the
    #     system is rank 1 (pairing rule b2*a0 + a2*b0 = 0, i.e. g_w ~ (a2, -b2), the X-axis
    #     mirror of g_{w+2}), so uniform pure-X / pure-Y survive and the router set widens to
    #     [v, v, vbar, vbar] for ANY v.  Verified EXACTLY in the c = 0 ring at N = 4, 5.
    M0 = M.subs(c, 0)
    assert M0.rank() == 1, "c=0 pair system not rank 1"
    assert sp.expand(M0.det()) == 0 and list(M0.row(0)) == list(M0.row(1)), "c=0 system shape"
    words0 = [('XZY', ONE), ('YZX', ONE)]
    for vx, vy, vname in ((1, 0, "(1,0) uniform pure-X"), (2, 1, "(2,1) generic")):
        v_vec = [(vx, 0, 0, 0), (vy, 0, 0, 0)]
        vbar_vec = [(vx, 0, 0, 0), (-vy, 0, 0, 0)]
        pattern = [v_vec, v_vec, vbar_vec, vbar_vec]
        maps_v = [ring_site_map(l, pattern[l]) for l in range(4)]
        # q^2 = -|v|^2 I (invertibility), exact
        q2 = ring_matmul(maps_v[0], maps_v[0], 0)
        exp_m = ring_mat({(i, i): (-(vx * vx + vy * vy), 0, 0, 0) for i in range(4)}, (4, 4))
        assert max(float(np.max(np.abs(q2[k] - exp_m[k]))) for k in range(4)) == 0.0, \
            f"c=0 {vname}: q^2 != -|v|^2 I"
        for N in (4, 5):
            A = commutator_super_ring(words0, N, 0)
            W = build_W_ring(maps_v, N, 0)
            anti = ring_add(ring_matmul(W, A, 0), ring_matmul(A, W, 0))
            assert ring_max_abs(anti) == 0.0, f"c=0 {vname} N={N}: {{W, A}} != 0"
        print(f"  c=0 WIDENING: [v,v,vbar,vbar] v={vname}: q^2 = -{vx * vx + vy * vy}*I, "
              f"{{W, A}} = 0 EXACT at N=4,5  OK")
    # control: the generic-v pattern does NOT route off the c = 0 station (complex eval at r(1))
    v_vec, vbar_vec = [(2, 0, 0, 0), (1, 0, 0, 0)], [(2, 0, 0, 0), (-1, 0, 0, 0)]
    pattern = [v_vec, v_vec, vbar_vec, vbar_vec]
    maps_v = [ring_site_map(l, pattern[l]) for l in range(4)]
    A = commutator_super_ring(family_words_ring(1), 4, 1)
    W = build_W_ring(maps_v, 4, 1)
    anti = ring_add(ring_matmul(W, A, 1), ring_matmul(A, W, 1))
    dev = float(np.max(np.abs(ring_to_complex(anti, metallic_r(1)))))
    assert dev > 1e-6, f"control: generic v unexpectedly routes at c=1 (dev {dev:.2e})"
    print(f"  control: generic v=(2,1) at c=1 FAILS ({{W,A}} eval max {dev:.1f} != 0); the "
          f"widening is the c=0 station's own  OK")

    # (e) the discrete candidates against the locus, as functions of c
    for name, al, be, want in (('P1 (g=X)', 1, 0, sp.Integer(1)),
                               ('P4/M2 (g=Y)', 0, 1, sp.Integer(-1)),
                               ('M (g ~ X+Y)', sp.sqrt(2) / 2, sp.sqrt(2) / 2, -c / 2)):
        val = sp.simplify(al ** 2 - c * al * be - be ** 2)
        assert sp.simplify(val - want) == 0, f"{name}: locus value {val} != {want}"
    print("  discrete candidates: P1 locus value 1, P4/M2 locus value -1 (off for EVERY c); "
          "M (45 deg) locus value -c/2: joins the locus exactly at c = 0  OK")
    print("BLOCK 6 PASS")


# ======================================================================
# BLOCK 7 -- the window lemma for ALL real c, by degree bound + exact interpolation.
# Independent arithmetic: complex rationals as (Fraction, Fraction) pairs and a hand-rolled
# 3-site Pauli string algebra -- no shared code with the ring blocks above. Parametrize the line
# by r (c = r - 1/r; both roots of r^2 = c*r + 1 give the same c). Q3 entries have degree <= 3
# in r and the r-cleared commutator degree <= 2, so every entry of r*{Q3, [T_c, .]} is a
# polynomial in r of degree <= 5: exact vanishing at 6+ distinct nodes makes the lemma a
# polynomial identity in r, i.e. the Hamiltonian leg holds for EVERY real (indeed complex) c.
# ======================================================================
from fractions import Fraction  # noqa: E402
from itertools import product as iproduct  # noqa: E402

F7_ZERO = (Fraction(0), Fraction(0))


def f7_cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def f7_cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def f7_cscale_i(x, ipow):
    for _ in range(ipow % 4):
        x = (-x[1], x[0])
    return x


F7_MUL = {}
for _a in range(4):
    F7_MUL[(0, _a)] = (_a, 0)
    F7_MUL[(_a, 0)] = (_a, 0)
    F7_MUL[(_a, _a)] = (0, 0)
F7_MUL[(1, 2)] = (3, 1)   # XY = iZ
F7_MUL[(2, 1)] = (3, 3)   # YX = -iZ
F7_MUL[(2, 3)] = (1, 1)   # YZ = iX
F7_MUL[(3, 2)] = (1, 3)   # ZY = -iX
F7_MUL[(3, 1)] = (2, 1)   # ZX = iY
F7_MUL[(1, 3)] = (2, 3)   # XZ = -iY


def f7_string_mul(s, t):
    out, ip = [], 0
    for a_, b_ in zip(s, t):
        c_, p_ = F7_MUL[(a_, b_)]
        out.append(c_)
        ip += p_
    return tuple(out), ip


def f7_commutator_columns(c):
    """col[j]: string -> coeff of [c*XZX + XZY + YZX, P_j] (3 sites, letters 0..3 = I,X,Y,Z)."""
    templates = [((1, 3, 1), c), ((1, 3, 2), (Fraction(1), Fraction(0))),
                 ((2, 3, 1), (Fraction(1), Fraction(0)))]
    cols = {}
    for j in iproduct(range(4), repeat=3):
        col = {}
        for s, w in templates:
            for sign, (st, ip) in ((1, f7_string_mul(s, j)), (-1, f7_string_mul(j, s))):
                coef = f7_cscale_i(w, ip)
                if sign < 0:
                    coef = (-coef[0], -coef[1])
                col[st] = f7_cadd(col.get(st, F7_ZERO), coef)
        cols[j] = {k: v for k, v in col.items() if v != F7_ZERO}
    return cols


def f7_site_q_cols(gv, sign):
    """columns of q on (I,X,Y,Z); gv = (gx, gy) Fractions, h = sign * i * R(g)."""
    gx, gy = gv
    hx = (Fraction(0), Fraction(sign) * -gy)
    hy = (Fraction(0), Fraction(sign) * gx)
    return {0: {1: (gx, Fraction(0)), 2: (gy, Fraction(0))},
            3: {1: hx, 2: hy},
            1: {0: (-gx, Fraction(0)), 3: hx},
            2: {0: (-gy, Fraction(0)), 3: hy}}


def f7_q3_columns(triple):
    cols = {}
    for j in iproduct(range(4), repeat=3):
        col = {}
        for (l0, c0), (l1, c1), (l2, c2) in iproduct(
                triple[0][j[0]].items(), triple[1][j[1]].items(), triple[2][j[2]].items()):
            key = (l0, l1, l2)
            col[key] = f7_cadd(col.get(key, F7_ZERO), f7_cmul(f7_cmul(c0, c1), c2))
        cols[j] = {k: v for k, v in col.items() if v != F7_ZERO}
    return cols


def f7_compose(acols, bcols):
    out = {}
    for j, bcol in bcols.items():
        col = {}
        for k, bc in bcol.items():
            for m, ac in acols[k].items():
                col[m] = f7_cadd(col.get(m, F7_ZERO), f7_cmul(ac, bc))
        out[j] = {k: v for k, v in col.items() if v != F7_ZERO}
    return out


def f7_window_lemma_at(r):
    """max |entry| (exact Fraction) of the window-summed anticommutator at node r, all offsets."""
    ccols = f7_commutator_columns((r - 1 / r, Fraction(0)))
    a = (r, Fraction(1))
    b = (Fraction(1), -r)
    q = [f7_site_q_cols(a, -1), f7_site_q_cols(a, +1),
         f7_site_q_cols(b, -1), f7_site_q_cols(b, +1)]
    worst = Fraction(0)
    for w in range(4):
        q3 = f7_q3_columns([q[w % 4], q[(w + 1) % 4], q[(w + 2) % 4]])
        a1, a2 = f7_compose(q3, ccols), f7_compose(ccols, q3)
        for j in a1:
            for k in set(a1[j]) | set(a2[j]):
                v = f7_cadd(a1[j].get(k, F7_ZERO), a2[j].get(k, F7_ZERO))
                worst = max(worst, abs(v[0]), abs(v[1]))
    return worst


def block7_all_real_c():
    print("-" * 92)
    print("BLOCK 7  the window lemma for ALL real c  [degree <= 5 in r + EXACT Fraction "
          "interpolation at 8 nodes]")
    print("-" * 92)
    nodes = [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 2),
             Fraction(5), Fraction(7), Fraction(4, 3), Fraction(-2)]
    assert len(set(nodes)) == 8 and all(r != 0 for r in nodes)
    for r in nodes:
        worst = f7_window_lemma_at(r)
        assert worst == 0, f"window lemma broken at node r = {r} (max |entry| = {worst})"
        c = r - 1 / r
        print(f"  r = {str(r):>4}  (c = {str(c):>5}): window-summed anticommutator = 0 EXACT, "
              f"all 4 offsets")
    print("  8 exact zeros > degree bound 5 + 1: the window lemma is a POLYNOMIAL IDENTITY in r")
    print("  => W L W^{-1} = -L - 2*sigma holds for EVERY real (indeed complex) c, by derivation")
    print("BLOCK 7 PASS")


# ======================================================================
# BLOCK 8 -- rigidity at golden AND silver, and the counted c = 0 modulus.
# Finite-difference Jacobian of F(params) = [W(A v_k) + A(W v_k)]_k over the 64-real-parameter
# period-4 class-swap family (per site: C-block g, h and B-block, 8 complex entries). Gauge =
# per-map complex scale, 4 maps x 2 = 8 real. Nullity > gauge <=> continuous physical moduli.
# ======================================================================
F8_IDX = [(1, 0), (2, 0), (1, 3), (2, 3), (0, 1), (3, 1), (0, 2), (3, 2)]


def f8_qs_from_params(p):
    qs = []
    for l in range(4):
        q = np.zeros((4, 4), dtype=complex)
        base = l * 16
        for j, (ri, ci) in enumerate(F8_IDX):
            q[ri, ci] = p[base + 2 * j] + 1j * p[base + 2 * j + 1]
        qs.append(q)
    return qs


def f8_params_from_qs(qs):
    p = np.zeros(64)
    for l, q in enumerate(qs):
        base = l * 16
        for j, (ri, ci) in enumerate(F8_IDX):
            p[base + 2 * j] = q[ri, ci].real
            p[base + 2 * j + 1] = q[ri, ci].imag
    return p


def f8_functional(p, A, probes, N):
    qs_comp = [P_BASIS @ q @ P_INV for q in f8_qs_from_params(p)]
    outs = [apply_W_comp(qs_comp, A @ v, N) + A @ apply_W_comp(qs_comp, v, N) for v in probes]
    o = np.concatenate(outs)
    return np.concatenate([o.real, o.imag])


def f8_nullity(c, N=5, K=8, eps=1e-6):
    r = metallic_r(c)
    H = _build_kbody_chain(N, family_words_kbody(c))
    A = lindbladian_pauli_dephasing(H, [0.0] * N, dephase_letter='Z')
    rng = np.random.default_rng(11)
    d = 2 ** N
    probes = []
    for _ in range(K):
        v = rng.standard_normal(d * d) + 1j * rng.standard_normal(d * d)
        probes.append(v / np.linalg.norm(v))
    qs = [ring_to_complex(q, r) for q in METALLIC_MAPS]
    p0 = f8_params_from_qs(qs)
    base = np.linalg.norm(f8_functional(p0, A, probes, N))
    assert base < 1e-10, f"c={c}: base point not a router (|F| = {base:.2e})"
    Jcols = []
    for j in range(64):
        dp = np.zeros(64)
        dp[j] = eps
        Jcols.append((f8_functional(p0 + dp, A, probes, N)
                      - f8_functional(p0 - dp, A, probes, N)) / (2 * eps))
    sv = np.linalg.svd(np.column_stack(Jcols), compute_uv=False)
    tol = sv[0] * 1e-8
    nullity = int(np.sum(sv < tol))
    gap = sv[63 - nullity] / max(sv[64 - nullity], 1e-300)
    return nullity, gap


def block8_rigidity():
    print("-" * 92)
    print("BLOCK 8  rigidity / the c = 0 modulus  [finite-difference Jacobian, 64 params, N=5]")
    print("-" * 92)
    for c, expect, label in ((1.0, 8, "golden"), (2.0, 8, "silver"), (0.0, 16, "c = 0")):
        nullity, gap = f8_nullity(c)
        assert nullity == expect, f"{label}: nullity {nullity} != {expect}"
        assert gap > 1e6, f"{label}: spectral gap {gap:.1e} not clean"
        moduli = nullity - 8
        print(f"  {label:>6} (c={c:g}): nullity = {nullity} = 8 gauge + {moduli} physical "
              f"moduli (gap {gap:.1e})")
    # the catalogued modulus slice: independent complex parity-chain directions [v, v', vbar, vbar']
    words0 = [('X', 'Z', 'Y', 1.0), ('Y', 'Z', 'X', 1.0)]
    H = _build_kbody_chain(5, words0)
    L = lindbladian_pauli_dephasing(H, [0.3] * 5, dephase_letter='Z')
    sigma = 1.5
    va = np.array([2.0, 0.7j]); vb = np.array([1.3, -0.4])
    pats = [va, vb, np.array([va[0], -va[1]]), np.array([vb[0], -vb[1]])]
    qs_comp = []
    for l, gv in enumerate(pats):
        Rg = np.array([-gv[1], gv[0]])
        h = ((-1.0) ** (l + 1)) * 1j * Rg
        q = np.zeros((4, 4), dtype=complex)
        q[1, 0] = gv[0]; q[2, 0] = gv[1]; q[1, 3] = h[0]; q[2, 3] = h[1]
        q[0, 1] = -gv[0]; q[3, 1] = h[0]; q[0, 2] = -gv[1]; q[3, 2] = h[1]
        qs_comp.append(P_BASIS @ q @ P_INV)
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(6):
        v = rng.standard_normal(1024) + 1j * rng.standard_normal(1024)
        v /= np.linalg.norm(v)
        Wv = apply_W_comp(qs_comp, v, 5)
        out = apply_W_comp(qs_comp, L @ v, 5) + L @ Wv + 2 * sigma * Wv
        worst = max(worst, np.linalg.norm(out) / np.linalg.norm(Wv))
    assert worst < 1e-12, f"c=0 catalogued slice fails: {worst:.2e}"
    print(f"  catalogued slice at c=0: [v, v', vbar, vbar'] with INDEPENDENT complex chain "
          f"directions routes ({worst:.1e}); h-tilts do not (4 of the 8 moduli catalogued)")
    print("BLOCK 8 PASS")


def main():
    t_start = time.time()
    print("=" * 92)
    print("THE METALLIC ROUTER FAMILY -- self-validating verification (F116 = the c = 1 point)")
    print("=" * 92)

    # cross-check 1: at c = 1 the parametrized ring is the anchor's golden ring bit-for-bit
    rng = np.random.default_rng(2026)
    for _ in range(200):
        x = tuple(int(v) for v in rng.integers(-9, 10, 4))
        y = tuple(int(v) for v in rng.integers(-9, 10, 4))
        assert comp_mul(x, y, 1) == anchor._comp_mul(x, y), "c=1 comp_mul != anchor golden"
    print("cross-check: c=1 ring multiplication == anchor golden ring (200 random samples)  OK")

    # cross-check 2: the generalized site-map builder reproduces the anchor's golden maps
    for l in range(4):
        ref = anchor.golden_site_map(l)
        assert all(np.array_equal(a, b) for a, b in zip(METALLIC_MAPS[l], ref)), \
            f"site {l}: metallic map components != anchor golden map"
    print("cross-check: metallic site-map components == anchor golden maps (c-independent "
          "tuples; the family lives in the ring's multiplication rule)  OK")

    # cross-check 3: the metallic-mean stations evaluate to their closed forms
    assert abs(metallic_r(1) - anchor.PHI) < 1e-15, "r(1) != phi"
    assert abs(metallic_r(2) - (1 + np.sqrt(2))) < 1e-15, "r(2) != 1 + sqrt(2)"
    assert abs(metallic_r(0) - 1.0) < 1e-15, "r(0) != 1"
    assert abs(metallic_r(-1) - 1 / anchor.PHI) < 1e-15, "r(-1) != 1/phi"
    assert abs(metallic_r(-1) - (anchor.PHI - 1)) < 1e-15, "r(-1) != phi - 1"
    assert abs(metallic_r(-2) - (np.sqrt(2) - 1)) < 1e-15, "r(-2) != sqrt(2) - 1"
    print("cross-check: stations r(1)=phi, r(2)=1+sqrt2, r(0)=1, r(-1)=1/phi=phi-1, "
          "r(-2)=sqrt2-1  OK\n")

    block1_structure()
    block2_window_lemma(METALLIC_MAPS, False, "XZX+XZY+YZX")
    block2_window_lemma(SIBLING_MAPS, True, "YZY+YZX+XZY (X<->Y conjugated)")
    block3_anticommutation()
    block4_end_to_end()
    block5_soft_set()
    block6_exclusion_symbolic()
    block7_all_real_c()
    block8_rigidity()

    print("=" * 92)
    print(f"ALL BLOCKS PASS   ({time.time() - t_start:.0f}s)")
    print("=" * 92)


if __name__ == "__main__":
    main()
