#!/usr/bin/env python3
r"""F87 windowed-converse MONOMIAL THEOREM: committed, self-validating verification.

This script is the computational anchor of the windowed-converse monomial theorem (Phase B). It
consolidates the trustworthy parts of two WIP scouts (`_f87_monomial_law.py`, the structural law,
and `_f87_rwinding_verify.py`, the two-reflection R-winding proof) into ONE committed file whose own
`assert`s are the test: every block raises on failure, prints a single PASS line on success, and the
process exits 0 only if all blocks pass.

The object
----------
Recenter the Liouvillian at the palindrome center, on the d^2 coherence space (d = 2^N, basis
|i><j| column-stacked, so kron(.,I) acts on the BRA index i and kron(I,.) on the KET index j):

    M(gamma) = A + gamma*Q ,
    A = -i[H,.] = A_L + A_R ,   A_L = -i(H (x) I)  (bra/left hop),  A_R = +i(I (x) H^T) (ket hop),
    Q = sum_l Z_l (x) Z_l ,     Q diagonal, Q_{ij} = N - 2*popcount(i^j) in [-N, N].

H is the windowed diagonal-cell Pauli pair (uniform J = 1). For a flux pair (odd #Y, an X..Y bond)
H is COMPLEX, so A is a general Gaussian-integer matrix; we carry both Re A and Im A (the corrected
generator). M is built so that  M(gamma) = lindbladian_pauli_dephasing(H, [gamma]*N, 'Z') + N*gamma*I
exactly, i.e. the framework Lindblad builder PROOF_F103 uses, recentred by sigma = N*gamma. So the
generator is correct by construction.

Two involutions on coherence space (F = X^{(x)N}, F^2 = I):
    Fcal = F (x) F :   Fcal A Fcal = -A ,                           Fcal Q Fcal = +Q
    R    = I (x) F :   R A_L R = +A_L ,  R A_R R = -A_R ,           R Q R    = -Q

The theorem (verified here)
---------------------------
Soft (spectrum symmetric about 0) <=> H's basis-state hopping graph is bipartite (no odd cycle). For
a HARD (non-bipartite, or diagonal-lift) pair with effective odd-cycle order ell:
  - the first nonvanishing odd power sum p_{m*}(gamma) = Tr(M^m*) is a strictly POSITIVE MONOMIAL
    c*gamma^deg, deg in {1, 3};
  - m* = 2*ell + deg (so the cycle length reaches the spectrum through the ORDER, not the degree);
  - deg = 1 only for a single-site-Z diagonal lift, deg = 3 for every cycle pair and multi-Z lift.
A positive monomial has no positive root, so the asymmetry is nonzero for ALL gamma > 0 => HARD.

Block ledger (rigor labels)
----------------------------
  Block 1  structural law              : verified-exhaustively (exact CRT) at N=4,5; the ell=5
                                         instance (N=6) is RIGOROUS-GENERAL by Blocks 2+4 and is
                                         re-derived numerically only under --heavy.
  Block 2  two-reflection sign table   : RIGOROUS-GENERAL (the signs are operator identities),
                                         verified EXACTLY incl. a complex-H flux pair.
  Block 3  all-odd parity              : RIGOROUS-GENERAL (per-word, from Block 2), verified by an
                                         exhaustive word census at the first nonvanishing moment.
  Block 4  threshold #A >= 2*ell       : RIGOROUS-GENERAL (unsigned odd-girth path-existence),
                                         verified incl. the flux pair (signed (H^3)_ii = 0 but
                                         unsigned (|H|^3)_ii > 0).
  Block 5  soft control                : verified-exhaustively (all odd p_m == 0 exactly).
  Block 6  deg-1 positivity closed form: RIGOROUS-GENERAL closed form, verified against exact p_3.
  Block 7  R-deg & R-sign cell-wide    : verified-exhaustively over all 50 hard pairs of the N=4
                                         k=3 Z cell; the +N-Perron skew over the 16 pure cycles.

Run
---
    python simulations/f87_windowed_monomial_converse.py            # fast (N<=5), all blocks
    python simulations/f87_windowed_monomial_converse.py --heavy    # + N=6 ell=5 exact certificate
                                                                     #   (~50 min, one full p_13)
Exact integer/CRT arithmetic carries the power sums (true zeros are exact rational 0, cleanly
separated from any nonzero); float64 eigenvalues are used only as an independent cross-check.
"""
from __future__ import annotations

import os
import sys
from collections import Counter, deque
from fractions import Fraction
from functools import reduce
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import _build_kbody_chain, site_op  # noqa: E402

import sympy as sp  # noqa: E402


# ======================================================================
# Exact generators  M(gamma) = A + gamma*Q  =  (A_L + A_R) + gamma*Q.
#   A = -i[H,.] is Gaussian-integer (Re + i*Im); Q is integer diagonal.
#   Built exactly the way the framework Lindblad builder builds -i[H,.], so A + gamma*Q reproduces
#   lindbladian_pauli_dephasing(H, [gamma]*N, 'Z') + N*gamma*I bit-for-bit (asserted below).
# ======================================================================
def build_H(N, pair):
    return _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])


def build_AL_AR_Q(N, pair):
    """Complex A_L, A_R, Q (and H) with A_L + A_R = -i[H,.] and Q = sum_l Z_l (x) Z_l."""
    H = build_H(N, pair)
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    A_L = -1j * np.kron(H, Id)          # -i (H (x) I)   acts on the BRA index
    A_R = +1j * np.kron(Id, H.T)        # +i (I (x) H^T) acts on the KET index
    Q = np.zeros((d * d, d * d), dtype=complex)
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        Q = Q + np.kron(Zl, Zl.conj())
    return A_L, A_R, Q.real.astype(complex), H


def build_integer_generators(N, pair):
    """int64 (Ar, Ai, Q) with M(gamma) = (Ar + i*Ai) + gamma*Q.  All three are integer (asserted)."""
    A_L, A_R, Q, H = build_AL_AR_Q(N, pair)
    A = A_L + A_R
    Ar, Ai = A.real, A.imag
    assert np.allclose(Ar, np.round(Ar)), "Re A not integer"
    assert np.allclose(Ai, np.round(Ai)), "Im A not integer"
    assert np.allclose(Q.real, np.round(Q.real)), "Q not integer"
    return (np.round(Ar).astype(np.int64),
            np.round(Ai).astype(np.int64),
            np.round(Q.real).astype(np.int64))


def assert_generator_matches_framework(N, pair):
    """M(gamma) == framework Lindbladian + N*gamma*I at several gamma (correct by construction)."""
    Ar, Ai, Q = build_integer_generators(N, pair)
    H = build_H(N, pair)
    for g in (0.3, 1.0, 2.7):
        L = lindbladian_pauli_dephasing(H, [g] * N, dephase_letter='Z')
        recon = (Ar + 1j * Ai).astype(complex) + g * Q.astype(complex)
        assert np.allclose(L + (N * g) * np.eye(4 ** N), recon, atol=1e-9), \
            f"generator mismatch vs framework at N={N} gamma={g}"


# ======================================================================
# Exact Tr(M^m) at integer gamma via CRT over mod-prime matmuls (the scouts' exact arithmetic).
#   The float64 modular matmul (the only BLAS-fast route) is EXACT while d2*p^2 < 2^53; the prime
#   bank is auto-sized to that bound. An int64 split-contraction path is the redundancy cross-check.
# ======================================================================
_CRT_PRIMES_LARGE = [1048573, 1048583, 1048601, 1048609, 1048631, 1048633,
                     1048661, 1048681, 1048703, 1048783, 1048789, 1048799]
_CRT_PRIMES_SMALL = [262139, 262147, 262151, 262153, 262187, 262193, 262217, 262231,
                     262237, 262253, 262261, 262271, 262273, 262277, 262279, 262289]


def _select_primes(d2, nprimes):
    bank = _CRT_PRIMES_LARGE if d2 <= 2 ** 11 else _CRT_PRIMES_SMALL
    assert d2 * bank[-1] ** 2 < 2 ** 53, f"prime bank unsafe for d2={d2}"
    return bank[:nprimes]


def _traces_mod_p(Ar, Ai, Q, g, mmax, p):
    """Exact Tr(M^m) mod p, m=1..mmax, via float64 BLAS matmul (exact while d2*p^2 < 2^53)."""
    PRr = (Ar % p).astype(np.float64)
    QR = (Q % p).astype(np.float64)
    R = (PRr + float(int(g) % p) * QR) % p
    Im = (Ai % p).astype(np.float64)
    PR, PI = R.copy() % p, Im.copy() % p
    res = {1: (int(round(np.trace(PR))) % p, int(round(np.trace(PI))) % p)}
    for m in range(2, mmax + 1):
        nR = (PR @ R - PI @ Im) % p
        nI = (PR @ Im + PI @ R) % p
        PR, PI = nR, nI
        res[m] = (int(round(np.trace(PR))) % p, int(round(np.trace(PI))) % p)
    return res


def _traces_mod_p_int(Ar, Ai, Q, g, mmax, p):
    """Exact Tr(M^m) mod p via int64 matmul, contraction split in two halves (overflow-safe to
    d2 = 4096). Slower (not BLAS-backed); used only as a redundancy cross-check."""
    pp = int(p)
    R = (((Ar % pp) + (int(g) % pp) * (Q % pp)) % pp).astype(np.int64)
    Im = (Ai % pp).astype(np.int64)
    d2 = Ar.shape[0]
    half = d2 // 2

    def mul(XR, XI, YR, YI):
        nR = (XR[:, :half] @ YR[:half] + XR[:, half:] @ YR[half:]) % pp
        nR = (nR - (XI[:, :half] @ YI[:half] + XI[:, half:] @ YI[half:])) % pp
        nI = (XR[:, :half] @ YI[:half] + XR[:, half:] @ YI[half:]) % pp
        nI = (nI + (XI[:, :half] @ YR[:half] + XI[:, half:] @ YR[half:])) % pp
        return nR, nI

    PR, PI = R.copy(), Im.copy()
    res = {1: (int(np.trace(PR)) % pp, int(np.trace(PI)) % pp)}
    for m in range(2, mmax + 1):
        PR, PI = mul(PR, PI, R, Im)
        res[m] = (int(np.trace(PR)) % pp, int(np.trace(PI)) % pp)
    return res


def _crt(residues, mods):
    M = reduce(lambda a, b: a * b, mods)
    x = 0
    for r, mi in zip(residues, mods):
        Mi = M // mi
        x = (x + r * Mi * pow(Mi, -1, mi)) % M
    return x - M if x > M // 2 else x


def traces_up_to_m_crt(Ar, Ai, Q, g, mmax, nprimes=10, use_int=False):
    """Exact dict m -> (re_int, im_int) of Tr(M^m), m=1..mmax, at integer gamma g via CRT."""
    d2 = Ar.shape[0]
    ps = _select_primes(d2, nprimes)
    fn = _traces_mod_p_int if use_int else _traces_mod_p
    per_p = [fn(Ar, Ai, Q, g, mmax, p) for p in ps]
    out = {}
    for m in range(1, mmax + 1):
        re = _crt([per_p[i][m][0] % ps[i] for i in range(len(ps))], ps)
        im = _crt([per_p[i][m][1] % ps[i] for i in range(len(ps))], ps)
        out[m] = (re, im)
    return out


def pm_polynomials_exact(Ar, Ai, Q, odd_ms, nprimes=10, use_int=False):
    """EXACT polynomials p_m(gamma) = Tr(M^m) for each odd m, by Lagrange interpolation over integer
    nodes. Returns dict m -> (P_real_lowfirst, P_imag_lowfirst) of sympy Rationals, and the symbol."""
    gamma = sp.symbols('gamma')
    mmax = max(odd_ms)
    nodes = list(range(0, mmax + 1))
    node_traces = {g: traces_up_to_m_crt(Ar, Ai, Q, g, mmax, nprimes=nprimes, use_int=use_int)
                   for g in nodes}

    def lagrange(vals):
        expr = sp.Integer(0)
        for i, yi in enumerate(vals):
            term = sp.Rational(yi.numerator, yi.denominator)
            for j in range(len(nodes)):
                if j == i:
                    continue
                term *= (gamma - nodes[j]) / sp.Integer(nodes[i] - nodes[j])
            expr += term
        poly = sp.Poly(sp.expand(expr), gamma)
        return [poly.coeff_monomial(gamma ** d) for d in range(mmax + 1)]

    out = {}
    for m in odd_ms:
        re_vals = [Fraction(node_traces[g][m][0]) for g in nodes]
        im_vals = [Fraction(node_traces[g][m][1]) for g in nodes]
        out[m] = (lagrange(re_vals), lagrange(im_vals))
    return out, gamma


def first_nonvanishing_odd(Ar, Ai, Q, max_m, nprimes=10, use_int=False):
    """Return (m*, lowfirst real coeff list, lowfirst imag coeff list, polys_dict) where m* is the
    first odd m with a nonzero real coefficient; m* = None if all odd p_m (m<=max_m) vanish."""
    odd_ms = list(range(1, max_m + 1, 2))
    polys, _ = pm_polynomials_exact(Ar, Ai, Q, odd_ms, nprimes=nprimes, use_int=use_int)
    for m in odd_ms:
        re_co, _im = polys[m]
        if any(c != 0 for c in re_co):
            return m, polys[m][0], polys[m][1], polys
    return None, None, None, polys


# ======================================================================
# Graph / obstruction helpers (the effective odd-cycle order ell).
# ======================================================================
def klein_index(t):
    return ((t.count('X') + t.count('Y')) % 2, (t.count('Y') + t.count('Z')) % 2)


def y_parity(t):
    return t.count('Y') % 2


def has_nonzero_diagonal(H, tol=1e-9):
    return float(np.max(np.abs(np.diag(H)))) > tol


def adjacency(H, tol=1e-9):
    Adj = np.abs(H).copy()
    np.fill_diagonal(Adj, 0.0)
    return Adj


def shortest_odd_cycle_graph(H, tol=1e-9):
    """Shortest odd cycle of H's basis-state hopping graph (BFS layers). 0 if bipartite."""
    d = H.shape[0]
    A = (np.abs(H) > tol)
    np.fill_diagonal(A, False)
    adj = [np.where(A[u])[0] for u in range(d)]
    best = 10 ** 9
    for s in range(d):
        dist = {s: 0}
        par = {s: -1}
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    par[v] = u
                    q.append(v)
                elif v != par[u] and (dist[u] + dist[v] + 1) % 2 == 1:
                    best = min(best, dist[u] + dist[v] + 1)
    return 0 if best == 10 ** 9 else best


def effective_ell(N, pair):
    """ell for the monomial law: 1 if H has a nonzero diagonal (a {I,Z} diagonal lift, no hopping
    cycle needed); else the hopping graph's shortest odd cycle (0 = bipartite = soft)."""
    H = build_H(N, pair)
    if has_nonzero_diagonal(H):
        return 1, 'diagonal-lift'
    g = shortest_odd_cycle_graph(H)
    return g, ('bipartite/soft' if g == 0 else f'{g}-cycle')


def F_op(N):
    Fx = np.array([[0, 1], [1, 0]], dtype=complex)
    out = Fx
    for _ in range(N - 1):
        out = np.kron(out, Fx)
    return out


# ======================================================================
# Exhaustive per-word exact trace census (letters over {'AL','AR','Q'}) for the parity block.
# ======================================================================
def census_word_traces(mat_map, mstar, d2, nprimes=4):
    """Exact Gaussian-integer Tr over EVERY length-mstar word in the letters of mat_map (letter ->
    (Re,Im) int64 arrays). The float64 modular matmul is exact while d2*p^2 < 2^53; the per-prime
    mod-reduced matrices are built ONCE and reused across all len(letters)^mstar words. Yields
    (word, re_int, im_int)."""
    bank = _CRT_PRIMES_LARGE if d2 <= 2 ** 11 else _CRT_PRIMES_SMALL
    ps = bank[:nprimes]
    letters = sorted(mat_map)
    reduced = [{L: ((mat_map[L][0] % p).astype(np.float64),
                    (mat_map[L][1] % p).astype(np.float64)) for L in letters} for p in ps]
    for word in product(letters, repeat=mstar):
        rr, ri = [], []
        for p, red in zip(ps, reduced):
            d2_ = red[word[0]][0].shape[0]
            Pr = np.eye(d2_, dtype=np.float64)
            Pi = np.zeros((d2_, d2_), dtype=np.float64)
            for L in word:
                Rr, Ri = red[L]
                nr = (Pr @ Rr - Pi @ Ri) % p
                ni = (Pr @ Ri + Pi @ Rr) % p
                Pr, Pi = nr, ni
            rr.append(int(round(np.trace(Pr))) % p)
            ri.append(int(round(np.trace(Pi))) % p)
        yield word, _crt(rr, ps), _crt(ri, ps)


def split_int(A_L, A_R, Q):
    """int64 (Re,Im) letter map for the {A_L, A_R, Q} word census."""
    ALr, ALi = np.round(A_L.real).astype(np.int64), np.round(A_L.imag).astype(np.int64)
    ARr, ARi = np.round(A_R.real).astype(np.int64), np.round(A_R.imag).astype(np.int64)
    Qr = np.round(Q.real).astype(np.int64)
    z = np.zeros_like(Qr)
    return {'AL': (ALr, ALi), 'AR': (ARr, ARi), 'Q': (Qr, z)}


# ======================================================================
# Q restricted to ker(A) (the omega=0 commutant block); +N Perron / -N reflection readings.
# ======================================================================
def Q_on_kerA_block(H, N, tol=1e-7):
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < tol]
    n = len(modes)
    Qb = np.zeros((n, n), dtype=complex)
    for i, (a, b) in enumerate(modes):
        for j, (ap, bp) in enumerate(modes):
            Qb[i, j] = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
    return Qb, modes


# ======================================================================
# Independent numpy float64 eigenvalue cross-check of an exact power-sum polynomial.
# ======================================================================
def numpy_power_sum(N, pair, m, gammas=(0.5, 1.0, 2.0)):
    H = build_H(N, pair)
    out = []
    for g in gammas:
        L = lindbladian_pauli_dephasing(H, [g] * N, dephase_letter='Z')
        ev = np.linalg.eigvals(L + (N * g) * np.eye(4 ** N))
        out.append((g, complex(np.sum(ev ** m))))
    return out


# ======================================================================
# Canonical pairs (the diagonal-cell, dephase-letter Z, uniform J = 1).
# ======================================================================
DIAG_LIFT = [('I', 'I', 'Z'), ('I', 'Z', 'I')]           # ell=1, deg=1
K3_EVEN = [('X', 'X', 'Z'), ('X', 'Z', 'X')]             # ell=3, deg=3, real H
FLUX = [('I', 'X', 'Y'), ('X', 'I', 'Y')]                # ell=3, deg=3, COMPLEX H (odd #Y)
MULTIZ = [('X', 'X', 'Z'), ('Z', 'Z', 'Z')]              # ell=1 diag-lift, deg=3, m*=5
SOFT = [('X', 'X', 'Z'), ('Z', 'X', 'X')]                # bipartite (soft control)
ELL5_N6 = [('I', 'X', 'Z', 'X'), ('X', 'I', 'Z', 'X')]   # ell=5 (N=6, k=4), heavy


def _fmt(c):
    return str(sp.nsimplify(c))


# ======================================================================
# BLOCK 1 -- the structural law (documented values; m* = 2*ell + deg; deg in {1,3}).
# ======================================================================
def block1_structural_law(heavy=False):
    print("-" * 92)
    print("BLOCK 1  structural law  [verified-exhaustively (exact CRT) N=4,5; ell=5 is heavy/opt-in]")
    print("-" * 92)
    # (a) K3 triangle: first nonvanishing odd moment is the positive monomial p_9. max_m=9 exactly
    #     recovers the degree-9 polynomial p_9 from 10 integer nodes (a tighter cap is cheaper at N=5).
    for N, want in ((4, 2064384), (5, 16515072)):
        Ar, Ai, Q = build_integer_generators(N, K3_EVEN)
        mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q, 9)
        nz = [j for j, c in enumerate(re_co) if c != 0]
        assert all(c == 0 for c in im_co), f"K3 N={N}: imag part nonzero"
        assert mstar == 9, f"K3 N={N}: m*={mstar} (expected 9)"
        assert nz == [3], f"K3 N={N}: gamma-powers {nz} (expected pure gamma^3)"
        assert re_co[3] == want, f"K3 N={N}: p_9 coeff {re_co[3]} (expected {want})"
        print(f"  K3 (XXZ+XZX) N={N}: p_9 = {re_co[3]} * gamma^3   (m*=9=2*3+3, deg=3)  OK")

    # (b) single-site-Z diagonal lift: p_3 = 9216 * gamma (N=4); deg=1.
    Ar, Ai, Q = build_integer_generators(4, DIAG_LIFT)
    mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q, 9)
    nz = [j for j, c in enumerate(re_co) if c != 0]
    assert all(c == 0 for c in im_co), "diag-lift N=4: imag part nonzero"
    assert mstar == 3 and nz == [1] and re_co[1] == 9216, \
        f"diag-lift N=4: m*={mstar}, powers {nz}, coeff {re_co[1] if len(re_co) > 1 else None}"
    print(f"  single-Z lift (IIZ+IZI) N=4: p_3 = {re_co[1]} * gamma^1   (m*=3=2*1+1, deg=1)  OK")

    # (c) m* = 2*ell + deg confirmed across the tractable cases (ell=1 deg=1, ell=3 deg=3) and the
    #     multi-site-Z lift (ell=1 but deg=3, m*=5).  ell determined independently from the graph.
    for label, N, pair, exp_ell, exp_deg, exp_mstar in (
            ("ell=1 deg=1 (IIZ+IZI)", 4, DIAG_LIFT, 1, 1, 3),
            ("ell=3 deg=3 (XXZ+XZX)", 4, K3_EVEN, 3, 3, 9),
            ("ell=3 deg=3 flux (IXY+XIY)", 4, FLUX, 3, 3, 9),
            ("ell=1 deg=3 multi-Z (XXZ+ZZZ)", 4, MULTIZ, 1, 3, 5)):
        ell, _kind = effective_ell(N, pair)
        assert ell == exp_ell, f"{label}: effective ell {ell} != {exp_ell}"
        Ar, Ai, Q = build_integer_generators(N, pair)
        # max_m=9 exactly recovers each p_m up to degree 9 (the largest m* in this set is 9).
        mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q, 9)
        nz = [j for j, c in enumerate(re_co) if c != 0]
        deg = nz[0]
        assert len(nz) == 1, f"{label}: first moment not a monomial (powers {nz})"
        assert deg == exp_deg, f"{label}: deg {deg} != {exp_deg}"
        assert mstar == exp_mstar, f"{label}: m* {mstar} != {exp_mstar}"
        assert mstar == 2 * ell + deg, f"{label}: m* {mstar} != 2*ell+deg {2*ell+deg}"
        assert re_co[deg] > 0, f"{label}: leading coeff not positive"
        print(f"  {label}: ell={ell} deg={deg} m*={mstar}=2*{ell}+{deg}; "
              f"coeff {re_co[deg]} > 0  OK")

    # (d) the ell=5 instance.  RIGOROUS-GENERAL by Blocks 2+4: a pure cycle has deg=3 (Block 2's
    #     R-deg kills the deg-1 class), and the threshold (Block 4) gives #A=2*ell, so m*=2*ell+deg
    #     = 13 != 3*ell = 15.  The smallest realiser is N=6 (k=4); its exact p_13 (one ~50-min run)
    #     is asserted only under --heavy.  The symbolic prediction is asserted unconditionally.
    ell5, deg5 = 5, 3
    assert 2 * ell5 + deg5 == 13 and 3 * ell5 == 15, "ell=5 arithmetic"
    print(f"  ell=5 (symbolic): m* = 2*5+3 = 13 (deg=3 != ell=5); task conj m*=3*ell=15 FAILS  OK")
    if heavy:
        print("  [--heavy] N=6 k=4 ell=5 exact certificate (IXZX+XIZX); this takes ~50 min ...",
              flush=True)
        Ar, Ai, Q = build_integer_generators(6, ELL5_N6)
        # confirm low odd moments vanish and p_13 is the positive monomial 50381979648 * gamma^3
        odd_ms = list(range(1, 14, 2))
        polys, _ = pm_polynomials_exact(Ar, Ai, Q, odd_ms, nprimes=7, use_int=False)
        for m in (1, 3, 5, 7, 9, 11):
            assert all(c == 0 for c in polys[m][0]), f"ell=5 N=6: p_{m} not exactly 0"
        re13 = polys[13][0]
        nz = [j for j, c in enumerate(re13) if c != 0]
        assert all(c == 0 for c in polys[13][1]), "ell=5 N=6: p_13 imag nonzero"
        assert nz == [3], f"ell=5 N=6: p_13 powers {nz} (expected pure gamma^3)"
        assert re13[3] == 50381979648, f"ell=5 N=6: p_13 coeff {re13[3]} (expected 50381979648)"
        print(f"  [--heavy] N=6 ell=5: p_1..p_11 = 0; p_13 = {re13[3]} * gamma^3 (m*=13)  OK")
    print("BLOCK 1 PASS")


# ======================================================================
# BLOCK 2 -- the two-reflection sign table (incl. a complex-H flux pair).  RIGOROUS-GENERAL.
# ======================================================================
def block2_reflection_signs():
    print("-" * 92)
    print("BLOCK 2  two-reflection sign table  [RIGOROUS-GENERAL operator identities; exact check]")
    print("-" * 92)
    cases = [("real-H K3", 4, K3_EVEN), ("complex-H flux IXY+XIY", 4, FLUX),
             ("real-H K3 N=5", 5, K3_EVEN)]
    for label, N, pair in cases:
        A_L, A_R, Q, H = build_AL_AR_Q(N, pair)
        F = F_op(N)
        Fcal = np.kron(F, F)
        R = np.kron(np.eye(2 ** N, dtype=complex), F)
        assert np.allclose(Fcal @ Fcal, np.eye(Fcal.shape[0])), f"{label}: Fcal^2 != I"
        assert np.allclose(R @ R, np.eye(R.shape[0])), f"{label}: R^2 != I"
        is_complex = not np.allclose(H.imag, 0)
        # the six signs (F, R are involutions so F^{-1}=F, R^{-1}=R)
        assert np.allclose(Fcal @ (A_L + A_R) @ Fcal, -(A_L + A_R)), f"{label}: Fcal A Fcal != -A"
        assert np.allclose(Fcal @ Q @ Fcal, Q), f"{label}: Fcal Q Fcal != +Q"
        assert np.allclose(R @ A_L @ R, A_L), f"{label}: R A_L R != +A_L"
        assert np.allclose(R @ A_R @ R, -A_R), f"{label}: R A_R R != -A_R"
        assert np.allclose(R @ Q @ R, -Q), f"{label}: R Q R != -Q"
        # the driving lemma:  F H F = -H  and  F H^T F = (F H F)^T = -H^T
        assert np.allclose(F @ H @ F, -H), f"{label}: F H F != -H"
        assert np.allclose(F @ H.T @ F, -H.T), f"{label}: F H^T F != -H^T"
        print(f"  [{label}] complex-H={is_complex}: Fcal A Fcal=-A, Fcal Q Fcal=+Q; "
              f"R A_L R=+A_L, R A_R R=-A_R, R Q R=-Q  OK")
    print("BLOCK 2 PASS")


# ======================================================================
# BLOCK 3 -- all-odd parity at the first nonvanishing odd moment.
#
# RIGOROUS-GENERAL (the proof): trace-conjugation invariance + Block 2's signs act PER WORD. For a
# length-m word W in {A_L, A_R, Q} with counts (a, b, c) = (#A_L, #A_R, #Q):
#     Fcal W Fcal = (-1)^(a+b) W  =>  Tr(W) = (-1)^(a+b) Tr(W)  =>  Tr(W) = 0 unless a+b even,
#     R    W R    = (-1)^(b+c) W  =>  Tr(W) = (-1)^(b+c) Tr(W)  =>  Tr(W) = 0 unless b+c even.
# For ODD m = a+b+c: a+b even forces c odd; b+c even with c odd forces b odd, hence a odd. So every
# surviving word at any odd power sum has (a, b, c) ALL ODD.  We assert this implication over every
# count-triple (no matmuls), then confirm it by an EXHAUSTIVE word census at a cheap scale, and tie
# the larger K3/flux m*=9 to Block 1's exact p_9 (a pure #Q=3 monomial => only the all-odd class).
# ======================================================================
def block3_all_odd_parity():
    print("-" * 92)
    print("BLOCK 3  all-odd parity at m*  [RIGOROUS-GENERAL count-arithmetic; exhaustive census check]")
    print("-" * 92)

    # (a) the count-arithmetic implication (the proof), over EVERY triple summing to a small odd m.
    for m in (1, 3, 5, 7, 9, 11, 13):
        for a in range(m + 1):
            for b in range(m - a + 1):
                c = m - a - b
                survives = ((a + b) % 2 == 0) and ((b + c) % 2 == 0)   # Block 2's two sign-rules
                all_odd = (a % 2 == 1) and (b % 2 == 1) and (c % 2 == 1)
                assert survives == all_odd, f"parity arithmetic fails at (a,b,c)=({a},{b},{c}) m={m}"
    print("  count-arithmetic: Block-2 signs force every surviving odd-m word to all-ODD counts  OK")

    # (b) EXHAUSTIVE word census at a cheap scale (multi-Z lift, m*=5: 3^5 = 243 words). Every
    #     nonzero word's counts are all ODD, and the census sum equals the exact p_5 (61440).
    A_L, A_R, Q, H = build_AL_AR_Q(4, MULTIZ)
    three = split_int(A_L, A_R, Q)
    d2 = A_L.shape[0]
    buckets = Counter()
    total_re = 0
    for word, re, im in census_word_traces(three, 5, d2):
        if re != 0 or im != 0:
            cnts = (word.count('AL'), word.count('AR'), word.count('Q'))
            buckets[cnts] += 1
            total_re += re
            assert all(x % 2 == 1 for x in cnts), f"multi-Z census: even-count nonzero word {word}"
    assert total_re == 61440, f"multi-Z census sum {total_re} != exact p_5 61440"
    print(f"  exhaustive census (multi-Z m*=5, 3^5 words): nonzero buckets {sorted(buckets)}; "
          f"sum = {total_re} == exact p_5; all counts ODD  OK")

    # (c) the K3 (real-H) and flux (complex-H) witnesses at m*=9: Block 1 already proved p_9 is a
    #     PURE gamma^3 monomial, i.e. the ONLY surviving class is #Q=3 (=> #A=6); with the (a)
    #     arithmetic that class is the all-odd (#A_L,#A_R,#Q) class.  This covers the complex-H pair
    #     without the 3^9 enumeration.
    for label, pair, want in (("real-H K3", K3_EVEN, 2064384), ("complex-H flux", FLUX, 589824)):
        Ar, Ai, Q3 = build_integer_generators(4, pair)
        mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q3, 9)
        nz = [j for j, c in enumerate(re_co) if c != 0]
        assert mstar == 9 and nz == [3] and re_co[3] == want and all(c == 0 for c in im_co), \
            f"{label}: p_9 not the pure gamma^3 monomial {want} (got m*={mstar}, powers {nz})"
        print(f"  [{label}] m*=9: only #Q=3 (=> #A=6, all-odd) class survives; "
              f"p_9 = {want} * gamma^3  OK")
    print("BLOCK 3 PASS")


# ======================================================================
# BLOCK 4 -- the threshold #A >= 2*ell via the UNSIGNED odd-girth (path existence), incl. flux.
# ======================================================================
def block4_threshold():
    print("-" * 92)
    print("BLOCK 4  threshold #A >= 2*ell  [RIGOROUS-GENERAL: unsigned odd-girth path-existence]")
    print("-" * 92)
    # The bra index returns to i via #A_L hops of H; a closed odd walk of length k < ell does NOT
    # exist (k < odd-girth), so #A_L >= ell (and #A_R >= ell), hence #A >= 2*ell.  Existence is the
    # UNSIGNED statement (|H|^k)_ii > 0.  For the FLUX pair the signed amplitude (H^3)_ii cancels to
    # 0 yet the unsigned 3-walk (|H|^3)_ii > 0 EXISTS, so the threshold still holds.
    for label, N, pair, ell, complex_expected in (
            ("real-H K3", 4, K3_EVEN, 3, False),
            ("complex-H flux", 4, FLUX, 3, True)):
        H = build_H(N, pair)
        assert (not np.allclose(H.imag, 0)) == complex_expected, f"{label}: complex-H mismatch"
        assert shortest_odd_cycle_graph(H) == ell, f"{label}: graph odd-girth != {ell}"
        Adj = adjacency(H)
        Ak = np.eye(H.shape[0])
        Hk = np.eye(H.shape[0], dtype=complex)
        unsigned, signed = {}, {}
        for k in range(1, ell + 1):
            Ak = Ak @ Adj
            Hk = Hk @ H
            unsigned[k] = float(np.max(np.diag(Ak)))
            signed[k] = float(np.max(np.abs(np.diag(Hk))))
        # no odd closed walk shorter than ell
        for k in range(1, ell, 2):
            assert unsigned[k] < 1e-9, f"{label}: unsigned closed {k}-walk exists (< ell)"
        assert unsigned[ell] > 1e-9, f"{label}: no unsigned closed ell-walk at k=ell"
        print(f"  [{label}] (|H|^k)_ii odd k<{ell}: 0; (|H|^{ell})_ii = {unsigned[ell]:.0f} > 0  "
              f"[signed (H^{ell})_ii = {signed[ell]:.1f}]  => #A_L,#A_R >= {ell}, #A >= {2*ell}  OK")
    # the flux pair must be the cancellation witness: signed (H^3)_ii == 0 while unsigned > 0
    Hf = build_H(4, FLUX)
    Adjf = adjacency(Hf)
    sgn3 = float(np.max(np.abs(np.diag(np.linalg.matrix_power(Hf, 3)))))
    uns3 = float(np.max(np.diag(np.linalg.matrix_power(Adjf, 3))))
    assert sgn3 < 1e-9 < uns3, f"flux witness: signed (H^3)_ii={sgn3}, unsigned={uns3}"
    print(f"  flux witness: signed (H^3)_ii = {sgn3:.1f} (cancels) but unsigned (|H|^3)_ii "
          f"= {uns3:.0f} > 0  OK")
    print("BLOCK 4 PASS")


# ======================================================================
# BLOCK 5 -- soft control: all odd power sums are exactly 0 (no odd closed walk).
# ======================================================================
def block5_soft_control():
    print("-" * 92)
    print("BLOCK 5  soft control  [verified-exhaustively: all odd p_m == 0 exactly]")
    print("-" * 92)
    H = build_H(4, SOFT)
    assert fw.classify_pauli_pair(fw.ChainSystem(N=4), [tuple(SOFT[0]), tuple(SOFT[1])],
                                  dephase_letter='Z') == 'soft', "SOFT control not classified soft"
    assert shortest_odd_cycle_graph(H) == 0, "SOFT control hopping graph is not bipartite"
    # no odd closed walk at any odd length (unsigned)
    Adj = adjacency(H)
    Ak = np.eye(H.shape[0])
    for k in range(1, 12):
        Ak = Ak @ Adj
        if k % 2 == 1:
            assert float(np.max(np.diag(Ak))) < 1e-9, f"SOFT: unsigned odd {k}-walk exists"
    Ar, Ai, Q = build_integer_generators(4, SOFT)
    odd_ms = list(range(1, 12, 2))
    polys, _ = pm_polynomials_exact(Ar, Ai, Q, odd_ms, nprimes=10)
    for m in odd_ms:
        assert all(c == 0 for c in polys[m][0]) and all(c == 0 for c in polys[m][1]), \
            f"SOFT: p_{m} not exactly 0"
    print(f"  (XXZ+ZXX) bipartite; all odd p_m == 0 exactly through m={odd_ms[-1]}  OK")
    print("BLOCK 5 PASS")


# ======================================================================
# BLOCK 6 -- deg-1 positivity closed form  P_{3,1} = 6 * sum_x deg_A(x)(w(x) - N/2).
# ======================================================================
def block6_deg1_closed_form():
    print("-" * 92)
    print("BLOCK 6  deg-1 closed form  [RIGOROUS-GENERAL closed form; checked against exact p_3]")
    print("-" * 92)
    # P_{3,1} = 3 Tr(A^2 Q) = -3 sum_x deg_A(x) Q_x = 6 sum_x deg_A(x)(w(x)-N/2), deg_A(x)=sum_y|A_xy|^2
    for N, want in ((4, 9216), (5, 61440)):
        A_L, A_R, Q, H = build_AL_AR_Q(N, DIAG_LIFT)
        A = A_L + A_R
        qdiag = np.round(Q.real).astype(np.int64).diagonal()
        degA = np.sum(np.abs(A) ** 2, axis=1)
        w = (N - qdiag) // 2
        closed = int(round(6 * np.sum(degA * (w - N / 2))))
        # exact P_{3,1} from the polynomial (gamma^1 coefficient of p_3)
        Ar, Ai, Qi = build_integer_generators(N, DIAG_LIFT)
        polys, _ = pm_polynomials_exact(Ar, Ai, Qi, [3], nprimes=10)
        P31 = int(polys[3][0][1])
        assert P31 == want, f"deg-1 N={N}: exact P_31={P31} (expected {want})"
        assert closed == want, f"deg-1 N={N}: closed form={closed} (expected {want})"
        assert P31 > 0, f"deg-1 N={N}: P_31 not positive"
        print(f"  N={N}: exact P_{{3,1}} = {P31} == 6*sum deg_A(x)(w(x)-N/2) = {closed} > 0  OK")
    print("BLOCK 6 PASS")


# ======================================================================
# BLOCK 7 -- R-deg & R-sign cell-wide over the N=4 k=3 Z diagonal cell.
# ======================================================================
def block7_cell_wide():
    print("-" * 92)
    print("BLOCK 7  cell-wide R-deg & R-sign  [verified-exhaustively: all 50 hard pairs, N=4 k=3]")
    print("-" * 92)
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [''.join(t) for t in product('IXYZ', repeat=k)
             if not all(L == 'I' for L in t) and klein_index(''.join(t)) == (0, 1)]
    n_hard = 0
    n_pos_monomial = 0
    pure_cycle = 0
    perron_pure = 0
    failures = []
    for t1, t2 in combinations_with_replacement(terms, 2):
        if y_parity(t1) != y_parity(t2):
            continue
        if fw.classify_pauli_pair(chain, [tuple(t1), tuple(t2)], dephase_letter='Z') != 'hard':
            continue
        n_hard += 1
        Ar, Ai, Q = build_integer_generators(N, [tuple(t1), tuple(t2)])
        # every hard pair in this cell breaks by m*<=9 (cycles 9, lifts 3 or 5); max_m=11 finds the
        # first nonvanishing odd moment with margin and recovers its (degree <=9) polynomial exactly.
        # nprimes=6: the node value |Tr(M^11)| at N=4, gamma<=11 is < d2*44^11 ~ 4e20, well inside
        # the product of 6 large primes (~1.05e6 each, product ~1.3e36), so the CRT is exact.
        mstar, re_co, im_co, _ = first_nonvanishing_odd(Ar, Ai, Q, 11, nprimes=6)
        nz = [j for j, c in enumerate(re_co) if c != 0] if re_co is not None else []
        if mstar is not None and len(nz) == 1 and re_co[nz[0]] > 0 and all(c == 0 for c in im_co):
            n_pos_monomial += 1
        else:
            failures.append((t1, t2, mstar, nz))
        # the +N-Perron skew is the §7.5 reading for the PURE-CYCLE pairs (no diagonal in H)
        H = build_H(N, [tuple(t1), tuple(t2)])
        if not has_nonzero_diagonal(H):
            pure_cycle += 1
            ev = np.linalg.eigvals(Q_on_kerA_block(H, N)[0]).real
            if np.any(np.abs(ev - N) < 1e-6) and not np.any(np.abs(ev + N) < 1e-6):
                perron_pure += 1
    assert n_hard == 50, f"expected 50 hard pairs, got {n_hard}"
    assert n_pos_monomial == 50, f"positive-monomial only {n_pos_monomial}/50; fails {failures}"
    assert pure_cycle == 16, f"expected 16 pure-cycle pairs, got {pure_cycle}"
    assert perron_pure == 16, f"+N Perron skew only {perron_pure}/16 pure cycles"
    print(f"  hard pairs {n_hard}: first nonvanishing odd p_m is a positive monomial "
          f"{n_pos_monomial}/{n_hard}  OK")
    print(f"  pure-cycle pairs {pure_cycle}: +N present & -N absent on Q|ker(A) "
          f"{perron_pure}/{pure_cycle}  OK")
    print("BLOCK 7 PASS")


def main():
    heavy = ("--heavy" in sys.argv) or (os.environ.get("RCPSI_HEAVY") == "1")
    print("=" * 92)
    print("F87 WINDOWED-CONVERSE MONOMIAL THEOREM -- self-validating verification")
    print(f"  heavy (N=6 ell=5 exact certificate): {heavy}")
    print("=" * 92)

    # generator sanity: the recentred Liouvillian equals the framework Lindbladian + sigma.
    for N, pair in ((4, K3_EVEN), (4, FLUX), (5, K3_EVEN)):
        assert_generator_matches_framework(N, pair)
    print("generator check: M(gamma) == framework Lindbladian + N*gamma*I  OK\n")

    # second-method cross-check (independent of the exact CRT route): float64 eigenvalues of the
    # framework-built (L + sigma) reproduce the exact K3 p_9 polynomial at several gamma.
    poly = sum(c * sp.symbols('gamma') ** j for j, c in
               enumerate(first_nonvanishing_odd(*build_integer_generators(4, K3_EVEN), 9)[1]))
    g = sp.symbols('gamma')
    for gv, pv in numpy_power_sum(4, K3_EVEN, 9):
        exact = float(poly.subs(g, gv))
        rel = abs(exact - pv.real) / max(1.0, abs(exact))
        assert rel < 1e-6, f"numpy cross-check of K3 p_9 failed at gamma={gv} (rel {rel:.1e})"
    print("cross-check: numpy float64 eigenvalue p_9 matches exact CRT polynomial (K3)  OK")

    # arithmetic redundancy cross-check: the float64-BLAS modular matmul and the int64 split-
    # contraction matmul must give bit-identical exact traces (two independent exact routes).
    Ar, Ai, Q = build_integer_generators(4, K3_EVEN)
    tr_f = traces_up_to_m_crt(Ar, Ai, Q, 2, 9, nprimes=4, use_int=False)
    tr_i = traces_up_to_m_crt(Ar, Ai, Q, 2, 9, nprimes=4, use_int=True)
    assert tr_f == tr_i, f"float64-CRT vs int64-CRT trace mismatch: {tr_f} vs {tr_i}"
    print("cross-check: float64-CRT and int64-CRT exact traces agree bit-for-bit (K3)  OK\n")

    block1_structural_law(heavy=heavy)
    block2_reflection_signs()
    block3_all_odd_parity()
    block4_threshold()
    block5_soft_control()
    block6_deg1_closed_form()
    block7_cell_wide()

    print("=" * 92)
    print("ALL BLOCKS PASS")
    print("=" * 92)


if __name__ == "__main__":
    main()
