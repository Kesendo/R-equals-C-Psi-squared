"""The closed form of the FULL cross form, exact over Q(i), term by term.

THEOREM (found 2026-07-14, lead-b scout; proved here):

    cross_form(x) = -(1/8) * Sum_{L in {+-1}^6, L1=+1} (Prod L) * T(L o x) * cot((L.x)/2)

equivalently: cross_form = -4 * (the totally-odd-in-each-variable part of T(x)*cot(s)),
s = (sum a + sum b)/2, with T the F127 core function (which since this morning has its
own closed form, simulations/f127_closed_form.py).

PROOF SHAPE: pure distributivity + a 288-fold elementary identity.  cross_form's literal
expansion (committed cross_triple_orthogonality.cross_form: prefactor cot(mu/2) times the
four-cotangent kernel Xh) is a flat list of 288 elementary terms, each

    (1/8) (-1)^(i+j) musign ssign c_xi c_up cot(mu/2) cot(arg/2),

where arg = mu + sxi*xi + sup*up is a FULL +-1 form L.x of all six angles (half a sheet
form).  The RHS flattens to 32 sheets x 9 blocks = 288 terms

    -(Prod L)/8 (-1)^(i+j) alpha_i(L o a) alpha_j(L o b) cot((L_i a_i + L_j b_j)/2) cot((L.x)/2).

Gates:
  [G1] BIJECTION: the (canonical sheet, block) tag hits each of the 288 cells exactly
       once on each side (exact integer combinatorics; this is the S2 incidence).
  [G2] TERM IDENTITY: for every cell, F-term == RHS-term as rational functions over
       Q(i) in z-coordinates (z_k = e^{i a_k}, w_k = e^{i b_k}; all arguments are
       integer angle combinations, so cot(./2) = i(M+1)/(M-1) with M a monomial).
  [G3] NUMERIC PIN: the flat F expansion reproduces the committed cross_form at random
       points (the expansion is distributivity, but pin it from below anyway), and the
       assembled RHS equals cross_form.

Consequences (stated, not gated here): with f127_closed_form.py's bracket substituted
for T, cross_form is fully explicit; the S3/S5 single-angle oddness of cross_form is
manifest (an odd-projection is totally odd); F127 (cross_form == 0 on V) becomes ONE
explicit 32-term flip-sum identity, no atoms, residues, or transport needed to state it.

Run: python simulations/f127_global_representation.py   (~1 min, exits 0 iff all gates pass)
"""
import itertools
import math
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, "simulations")
from cross_triple_orthogonality import cross_form

I = sp.I
Z = sp.symbols("z1 z2 z3")
W = sp.symbols("w1 w2 w3")
ZW = list(Z) + list(W)


def mono(vec):
    """The monomial e^{i * vec.(a,b)} in z-coordinates."""
    m = sp.Integer(1)
    for c, v in zip(vec, ZW):
        m *= v ** c
    return m


def cot_half(vec):
    """cot(vec.(a,b)/2) = i (M+1)/(M-1), M = mono(vec)."""
    m = mono(vec)
    return I * (m + 1) / (m - 1)


def sin_full(k, sign=1):
    """sin(sign * angle_k), angle order (a1,a2,a3,b1,b2,b3)."""
    v = ZW[k]
    return sign * (v - 1 / v) / (2 * I)


def canonical(vec):
    return tuple(vec) if vec[0] == 1 else tuple(-c for c in vec)


def F_term_list():
    """The 288 elementary terms of cross_form, symbolic, tagged (sheet, block)."""
    terms = []
    for i in range(3):
        ja, la = [t for t in range(3) if t != i]
        for j in range(3):
            jb, lb = [t for t in range(3) if t != j]
            for musign in (+1, -1):
                muvec = [0] * 6
                muvec[i] = 1
                muvec[3 + j] = musign
                for xi_name in ("psum", "pdif"):
                    # c_xi: alpha = sin a_l - sin a_j;  beta = -(sin a_l + sin a_j)
                    if xi_name == "psum":
                        cxi = sin_full(la) - sin_full(ja)
                        xivec = [0] * 6
                        xivec[ja] = xivec[la] = 1
                    else:
                        cxi = -(sin_full(la) + sin_full(ja))
                        xivec = [0] * 6
                        xivec[ja], xivec[la] = -1, 1
                    for up_name in ("psum", "pdif"):
                        if up_name == "psum":
                            cup = sin_full(3 + lb) - sin_full(3 + jb)
                            upvec = [0] * 6
                            upvec[3 + jb] = upvec[3 + lb] = 1
                        else:
                            cup = -(sin_full(3 + lb) + sin_full(3 + jb))
                            upvec = [0] * 6
                            upvec[3 + jb], upvec[3 + lb] = -1, 1
                        for sxi, sup, ssign in ((1, -1, 1), (-1, 1, 1),
                                                (1, 1, -1), (-1, -1, -1)):
                            vec = [muvec[k] + sxi * xivec[k] + sup * upvec[k]
                                   for k in range(6)]
                            L = canonical(vec)
                            coeff = sp.Rational(1, 8) * (-1) ** (i + j) * musign * ssign
                            term = coeff * cxi * cup * cot_half(muvec) * cot_half(vec)
                            terms.append((L, (i, j), term))
    return terms


def RHS_term_list():
    """The 32 x 9 terms of -(1/8) sum_L (Prod L) T(L o x) cot((L.x)/2), symbolic."""
    terms = []
    for L in [(1,) + t for t in itertools.product((1, -1), repeat=5)]:
        prodL = 1
        for c in L:
            prodL *= c
        for i in range(3):
            ja, la = [t for t in range(3) if t != i]
            for j in range(3):
                jb, lb = [t for t in range(3) if t != j]
                # alpha_i(L o a) = sin(L_la a_la) - sin(L_ja a_ja); same on b
                ai = sin_full(la, L[la]) - sin_full(ja, L[ja])
                bj = sin_full(3 + lb, L[3 + lb]) - sin_full(3 + jb, L[3 + jb])
                cvec = [0] * 6
                cvec[i], cvec[3 + j] = L[i], L[3 + j]
                term = (-sp.Rational(1, 8) * prodL * (-1) ** (i + j)
                        * ai * bj * cot_half(cvec) * cot_half(list(L)))
                terms.append((L, (i, j), term))
    return terms


def gate1_bijection(ft, rt):
    fkeys = sorted((L, blk) for L, blk, _ in ft)
    rkeys = sorted((L, blk) for L, blk, _ in rt)
    assert len(ft) == len(rt) == 288
    assert len(set(fkeys)) == 288, "F tags collide"
    assert len(set(rkeys)) == 288, "RHS tags collide"
    assert fkeys == rkeys, "tag sets differ"
    print("[G1] bijection: 288 = 32 sheets x 9 blocks, each cell once per side  PASS")


def gate2_terms(ft, rt):
    fd = {(L, blk): t for L, blk, t in ft}
    rd = {(L, blk): t for L, blk, t in rt}
    n = 0
    for key in fd:
        diff = sp.cancel(sp.together(fd[key] - rd[key]))
        assert diff == 0, f"term identity FAILED at {key}"
        n += 1
    print(f"[G2] term identity: {n}/288 cells equal exactly over Q(i)  PASS")


def gate3_numeric(ft, rt):
    rng = np.random.default_rng(9)
    lam_f = sp.lambdify(ZW, sum(t for _, _, t in ft), "numpy")
    lam_r = sp.lambdify(ZW, sum(t for _, _, t in rt), "numpy")
    worst_f = worst_r = 0.0
    for _ in range(40):
        ang = rng.uniform(0.3, 2 * np.pi - 0.3, 6)
        zw = np.exp(1j * ang)
        ref = cross_form(tuple(ang[:3]), tuple(ang[3:]))
        worst_f = max(worst_f, abs(complex(lam_f(*zw)) - ref))
        worst_r = max(worst_r, abs(complex(lam_r(*zw)) - ref))
    assert worst_f < 1e-8 and worst_r < 1e-8, (worst_f, worst_r)
    print(f"[G3] numeric pin vs committed cross_form, 40 pts: "
          f"flat-F dev {worst_f:.2e}, RHS dev {worst_r:.2e}  PASS")


if __name__ == "__main__":
    ft = F_term_list()
    rt = RHS_term_list()
    gate1_bijection(ft, rt)
    gate3_numeric(ft, rt)
    gate2_terms(ft, rt)
    print("ALL GATES PASS: cross_form = -(1/8) sum_L (Prod L) T(L o x) cot((L.x)/2), "
          "exact over Q(i), term by term.")
