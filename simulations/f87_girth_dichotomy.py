#!/usr/bin/env python3
r"""F87 girth dichotomy: the deg-1 class in closed form; R-deg dissolves (committed verifier).

The single-Q coefficients of the recentred odd power sums factorize through the Z_l-weighted
moments t_j = Tr(Z_l H^j) (the supertrace identity), even-j moments die on the F-chirality,
the unsigned odd-girth kills t_j for j < ell, and at the first candidate moment the binomial
sum collapses to a SUM OF SQUARES:

    P_{2*ell+1, 1} = (2*ell+1) * C(2*ell, ell) * Sum_l t_ell^{(l)} ** 2  >= 0.

The GIRTH LADDER follows, replacing the former R-deg residual:
  - t_ell != 0 somewhere  =>  m* = 2*ell + 1, deg = 1, and p_{m*} = P_{2ell+1,1} * gamma is a
    POSITIVE monomial outright (threshold kills every odd moment below, and the #Q = 3 class at
    m* needs #A = 2*ell - 2 < 2*ell): hard at every gamma > 0, NO residual.
  - t_ell == 0 everywhere  =>  the deg-1 class is dead at 2*ell+1 AND 2*ell+3 (proven); if the
    gamma^3 class fires, p_{2*ell+3} = P_{2ell+3,3} * gamma^3 exactly (monomial PROVEN at that
    rung); if it does not, the ladder continues to higher odd rungs (k=4 witness IIXY+ZXZY:
    m* = 11 = 2*ell+5, p_11 = 86507520 * gamma^5, positive). The open residual (R-sign, ladder
    form) is that the first surviving class is single and positive.
The single-site-Z lift is the ell = 1 face: t_1 = 2^N * c_l, and the closed form reproduces
P_{3,1} = 6 * 4^N * Sum c_l^2 bit-exactly. The k=3 taxonomy ("deg = 1 only for a single-site-Z
lift", deg in {1,3}) was a k=3-cell fact: at k = 4 exactly 20 pure-cycle pairs carry t_3 != 0
and fire deg 1 at m* = 7 (first representative IXXZ+XIXZ, coefficient 573440), and the gamma^5
rung appears (IIXY+ZXZY).

Block ledger
------------
  Block 1  supertrace factorization      : RIGOROUS-GENERAL identity, checked EXACTLY against
                                           the direct superoperator at N=4 (5 canonical pairs
                                           + IXXZ+XIXZ, k = 1..3)
  Block 2  girth kill + closed form      : RIGOROUS-GENERAL; P_{2ell+1,1} formula asserted
                                           against the exact CRT gamma^1 coefficients for all
                                           four branch representatives (incl. the ell=1 face
                                           == the P31 = 6*4^N*sum c_l^2 identity)
  Block 3  k=4 dichotomy battery         : census of all 1056 same-y-par pairs (branch counts
                                           pinned), then exact verification of ALL 20 deg-1
                                           pure-cycle pairs (m* = 7, single power, coefficient
                                           == the closed form) + stratified samples verifying
                                           the ladder law + the pinned gamma^5 rung (IIXY+ZXZY)
  Block 4  k=3 regression                : the dichotomy re-derives Block 7 of the anchor: all
                                           16 pure cycles have t_3 == 0 (deg-3 branch), the 10
                                           diagonal lifts split by t_1 exactly as the cell law

Run: python simulations/f87_girth_dichotomy.py   (~20-30 min; Block 3's exact verifications at
N=5 dominate). Anchor machinery: f87_windowed_monomial_converse.py. 2026-06-10.
"""
from __future__ import annotations

import sys
from collections import Counter
from itertools import combinations_with_replacement, product as iproduct
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import site_op  # noqa: E402
import f87_windowed_monomial_converse as anchor  # noqa: E402


# ======================================================================
# Moments and the closed form.
# ======================================================================
def t_moment(H, N, j):
    """t_j^{(l)} = Tr(Z_l H^j) for every site l (exact integers for integer H)."""
    Hp = np.linalg.matrix_power(H, j)
    return [complex(np.trace(site_op(N, l, 'Z') @ Hp)) for l in range(N)]


def sum_sq(ts):
    return sum(t.real ** 2 + t.imag ** 2 for t in ts)


def closed_form_p2lp1(ell, ts):
    """P_{2*ell+1, 1} = (2*ell+1) * C(2*ell, ell) * Sum_l t_ell^2."""
    return (2 * ell + 1) * comb(2 * ell, ell) * sum_sq(ts)


# ======================================================================
# BLOCK 1 -- the supertrace factorization (exact vs direct superoperator).
# ======================================================================
def trQA2k_direct(H, N, k):
    d = 2 ** N
    A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    Q = sum(np.kron(site_op(N, l, 'Z'), site_op(N, l, 'Z').conj()) for l in range(N))
    return complex(np.trace(Q @ np.linalg.matrix_power(A, 2 * k)))


def trQA2k_formula(H, N, k):
    t = {j: t_moment(H, N, j) for j in range(2 * k + 1)}
    s = 0
    for l in range(N):
        for j in range(2 * k + 1):
            s += (-1) ** j * comb(2 * k, j) * t[2 * k - j][l] * t[j][l]
    return (-1) ** k * s


def block1_factorization():
    print("-" * 92)
    print("BLOCK 1  supertrace factorization  [RIGOROUS-GENERAL; exact vs direct superop, N=4]")
    print("-" * 92)
    cases = dict(anchor.__dict__[n] for n in ()) if False else {
        'K3': anchor.K3_EVEN, 'FLUX': anchor.FLUX, 'MULTIZ': anchor.MULTIZ,
        'DIAGLIFT': anchor.DIAG_LIFT, 'SOFT': anchor.SOFT,
        'IXXZ+XIXZ': [('I', 'X', 'X', 'Z'), ('X', 'I', 'X', 'Z')]}
    N = 4
    for name, pair in cases.items():
        if len(pair[0]) > N:
            continue
        H = anchor.build_H(N, pair)
        for k in (1, 2, 3):
            diff = abs(trQA2k_formula(H, N, k) - trQA2k_direct(H, N, k))
            assert diff == 0.0, f"{name} k={k}: factorization mismatch ({diff})"
        # even-j moments vanish (F-chirality)
        for j in (2, 4):
            assert max(abs(t) for t in t_moment(H, N, j)) == 0.0, f"{name}: t_{j} != 0"
        print(f"  {name:10s}: Tr(Q A^(2k)) == binomial bilinear in t_j (k=1,2,3, exact); "
              f"even-j moments = 0  OK")
    print("BLOCK 1 PASS")


# ======================================================================
# BLOCK 2 -- girth kill + the sum-of-squares closed form at every branch representative.
# ======================================================================
def gamma1_coeff_of_pm(pair, N, m, nprimes=8):
    """The exact gamma^1 coefficient of p_m via the anchor's CRT polynomial machinery."""
    Ar, Ai, Q = anchor.build_integer_generators(N, pair)
    polys, _ = anchor.pm_polynomials_exact(Ar, Ai, Q, [m], nprimes=nprimes)
    re_co, im_co = polys[m]
    assert all(c == 0 for c in im_co), f"p_{m} imaginary part nonzero"
    return int(re_co[1])


def block2_closed_form():
    print("-" * 92)
    print("BLOCK 2  girth kill + P_(2ell+1,1) = (2ell+1) C(2ell,ell) sum t_ell^2  [exact]")
    print("-" * 92)
    reps = [
        ("DIAGLIFT (ell=1, t_1 != 0)", anchor.DIAG_LIFT, 4, 1),
        ("MULTIZ (ell=1, t_1 == 0)", anchor.MULTIZ, 4, 1),
        ("K3 (ell=3, t_3 == 0)", anchor.K3_EVEN, 4, 3),
        ("IXXZ+XIXZ (ell=3, t_3 != 0)", [('I', 'X', 'X', 'Z'), ('X', 'I', 'X', 'Z')], 5, 3),
    ]
    for name, pair, N, ell in reps:
        H = anchor.build_H(N, pair)
        got_ell, _ = anchor.effective_ell(N, pair)
        assert got_ell == ell, f"{name}: ell {got_ell} != {ell}"
        # girth kill below ell (pure off-diagonal cases)
        if not anchor.has_nonzero_diagonal(H):
            for j in range(1, ell):
                assert max(abs(t) for t in t_moment(H, N, j)) == 0.0, \
                    f"{name}: t_{j} != 0 below the girth"
        ts = t_moment(H, N, ell)
        pred = closed_form_p2lp1(ell, ts)
        got = gamma1_coeff_of_pm(pair, N, 2 * ell + 1)
        assert got == round(pred), f"{name}: P_(2ell+1,1) {got} != closed form {pred}"
        print(f"  {name:30s}: sum t_ell^2 = {sum_sq(ts):.0f}; "
              f"P_({2*ell+1},1) = {got} == (2ell+1)C(2ell,ell)*sum  OK")
        if name.startswith("DIAGLIFT"):
            # the ell=1 face IS the P31 identity: 3*C(2,1)*sum t_1^2 = 6*sum(2^N c_l)^2
            assert got == 9216, "ell=1 face does not reproduce the P31 = 9216 identity"
            print(f"    (the ell=1 face reproduces P_(3,1) = 6*4^N*sum c_l^2 = 9216 bit-exact)")
    print("BLOCK 2 PASS")


# ======================================================================
# BLOCK 3 -- the k=4 dichotomy battery.
# ======================================================================
def k4_pairs():
    cells = [t for t in iproduct('IXYZ', repeat=4)
             if not all(c == 'I' for c in t)
             and (''.join(t).count('X') + ''.join(t).count('Y')) % 2 == 0
             and (''.join(t).count('Y') + ''.join(t).count('Z')) % 2 == 1]
    ypar = lambda t: ''.join(t).count('Y') % 2
    for c1, c2 in combinations_with_replacement(cells, 2):
        if ypar(c1) == ypar(c2):
            yield [tuple(c1), tuple(c2)]


def block3_k4_battery():
    print("-" * 92)
    print("BLOCK 3  k=4 dichotomy battery  [census of 1056 pairs + exact verification]")
    print("-" * 92)
    N = 5
    census = Counter()
    deg1_cycles = []
    samples = {('cycle', 'deg3'): [], ('lift', 'deg1'): [], ('lift', 'deg3'): []}
    for pair in k4_pairs():
        ell, kind = anchor.effective_ell(N, pair)
        if ell == 0:
            census['soft'] += 1
            continue
        H = anchor.build_H(N, pair)
        ts = t_moment(H, N, ell)
        nz = sum_sq(ts) > 1e-9
        diag = 'lift' if kind == 'diagonal-lift' else 'cycle'
        census[(diag, ell, 'deg1' if nz else 'deg3')] += 1
        if diag == 'cycle' and nz:
            deg1_cycles.append((pair, ell, sum_sq(ts)))
        else:
            key = (diag, 'deg1' if nz else 'deg3')
            if key in samples and len(samples[key]) < 4:
                samples[key].append((pair, ell, sum_sq(ts)))
    expected = {'soft': 636, ('cycle', 3, 'deg1'): 20, ('cycle', 3, 'deg3'): 172,
                ('lift', 1, 'deg1'): 122, ('lift', 1, 'deg3'): 106}
    assert dict(census) == expected, f"census {dict(census)} != {expected}"
    print(f"  census: {dict(census)}  (pinned)  OK")

    # exact verification of ALL 20 deg-1 pure cycles: m* = 7, single power gamma^1,
    # coefficient == the closed form, POSITIVE
    print(f"  verifying all {len(deg1_cycles)} deg-1 pure-cycle pairs exactly (m*=7, "
          f"coefficient == closed form):", flush=True)
    for pair, ell, ss in deg1_cycles:
        Ar, Ai, Q = anchor.build_integer_generators(N, pair)
        mstar, re_co, im_co, _ = anchor.first_nonvanishing_odd(Ar, Ai, Q, 9, nprimes=8)
        nz = [j for j, c in enumerate(re_co) if c != 0]
        pred = round(closed_form_p2lp1(ell, t_moment(anchor.build_H(N, pair), N, ell)))
        label = '+'.join(''.join(t) for t in pair)
        assert mstar == 2 * ell + 1, f"{label}: m* {mstar} != {2*ell+1}"
        assert nz == [1], f"{label}: powers {nz} != [1]"
        assert re_co[1] == pred and pred > 0, \
            f"{label}: coefficient {re_co[1]} != closed form {pred}"
        print(f"    {label}: m*=7, p_7 = {re_co[1]}*gamma (== closed form, > 0)  OK", flush=True)

    # stratified samples of the other branches: the LADDER law. On the t_ell = 0 branch the
    # gamma^1 class is dead at 2*ell+1 and 2*ell+3 (proven), but the gamma^3 class need not
    # fire: the first survivor can sit higher (k=4 witness IIXY+ZXZY: m* = 11 = 2*ell+5,
    # gamma^5, coefficient 86507520 > 0). The verified law: m* = 2*ell + deg with deg odd, a
    # SINGLE surviving class, positive coefficient.
    for (diag, branch), rows in samples.items():
        for pair, ell, ss in rows:
            Ar, Ai, Q = anchor.build_integer_generators(N, pair)
            mstar, re_co, im_co, _ = anchor.first_nonvanishing_odd(Ar, Ai, Q, 13, nprimes=6)
            nz = [j for j, c in enumerate(re_co) if c != 0]
            label = '+'.join(''.join(t) for t in pair)
            assert mstar is not None and len(nz) == 1, f"{label}: not a monomial ({mstar},{nz})"
            deg = nz[0]
            if branch == 'deg1':
                assert deg == 1 and mstar == 2 * ell + 1, f"{label}: ({mstar},{nz})"
            else:
                assert deg % 2 == 1 and deg >= 3 and mstar == 2 * ell + deg, \
                    f"{label}: ({mstar},{nz}) off the ladder"
            assert re_co[deg] > 0, f"{label}: coefficient {re_co[deg]} not positive"
            print(f"    {label} [{diag}/{branch}]: m*={mstar} = 2*{ell}+{deg}, single power "
                  f"gamma^{deg}, coeff {re_co[deg]} > 0  OK", flush=True)

    # the gamma^5 rung, pinned: the k=4 pair where the gamma^3 class ALSO dies
    pair5 = [('I', 'I', 'X', 'Y'), ('Z', 'X', 'Z', 'Y')]
    Ar, Ai, Q = anchor.build_integer_generators(N, pair5)
    mstar, re_co, im_co, _ = anchor.first_nonvanishing_odd(Ar, Ai, Q, 13, nprimes=6)
    nz = [j for j, c in enumerate(re_co) if c != 0]
    assert mstar == 11 and nz == [5] and re_co[5] == 86507520, \
        f"IIXY+ZXZY: ({mstar},{nz},{re_co[nz[0]] if nz else None})"
    print(f"    IIXY+ZXZY [cycle, gamma^5 rung]: m*=11 = 2*3+5, p_11 = 86507520*gamma^5 > 0 "
          f"(P_(2ell+3,3) = 0 here: the ladder continues past the gamma^3 rung)  OK")
    print("BLOCK 3 PASS")


# ======================================================================
# BLOCK 4 -- k=3 regression: the dichotomy re-derives the anchor's cell law.
# ======================================================================
def block4_k3_regression():
    print("-" * 92)
    print("BLOCK 4  k=3 regression  [all 16 pure cycles t_3 == 0; lifts split by t_1]")
    print("-" * 92)
    N = 4
    import framework as fw
    chain = fw.ChainSystem(N=N)
    terms = [''.join(t) for t in iproduct('IXYZ', repeat=3)
             if not all(c == 'I' for c in t) and anchor.klein_index(''.join(t)) == (0, 1)]
    n_cycle = n_lift1 = n_lift3 = 0
    for t1, t2 in combinations_with_replacement(terms, 2):
        if anchor.y_parity(t1) != anchor.y_parity(t2):
            continue
        if fw.classify_pauli_pair(chain, [tuple(t1), tuple(t2)], dephase_letter='Z') != 'hard':
            continue
        pair = [tuple(t1), tuple(t2)]
        H = anchor.build_H(N, pair)
        ell, kind = anchor.effective_ell(N, pair)
        ts = t_moment(H, N, ell)
        if kind == 'diagonal-lift':
            if sum_sq(ts) > 1e-9:
                n_lift1 += 1
            else:
                n_lift3 += 1
        else:
            assert sum_sq(ts) == 0.0, f"{t1}+{t2}: k=3 pure cycle with t_3 != 0?!"
            n_cycle += 1
    assert n_cycle == 16, f"pure cycles {n_cycle} != 16"
    assert n_lift1 + n_lift3 == 34, f"lifts {n_lift1}+{n_lift3} != 34"
    print(f"  16/16 pure cycles have t_3 == 0 exactly (deg-3 branch: the k=3 cell had no "
          f"deg-1 cycles, which is why the old taxonomy held there)")
    print(f"  lifts: {n_lift1} with t_1 != 0 (deg 1 at m*=3), {n_lift3} with t_1 == 0 "
          f"(deg 3 at m*=5)  OK")
    print("BLOCK 4 PASS")


# ======================================================================
# BLOCK 5 -- the GF(2)[x] production channels: B monomial <=> the deg-1 channel can open.
# ======================================================================
def _poly(letters, which):
    m = 0
    for i, c in enumerate(letters):
        if (which == 'a' and c in 'XY') or (which == 'b' and c in 'YZ'):
            m |= 1 << i
    return m


def _deg(p):
    return p.bit_length() - 1


def _gf2_mod(a, b):
    while b and _deg(a) >= _deg(b):
        a ^= b << (_deg(a) - _deg(b))
    return a


def _gf2_gcd(a, b):
    while b:
        a, b = b, _gf2_mod(a, b)
    return a


def _gf2_div(a, b):
    q = 0
    while a and _deg(a) >= _deg(b):
        sh = _deg(a) - _deg(b)
        q |= 1 << sh
        a ^= b << sh
    return q


def _gf2_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def _is_monomial(p):
    return p != 0 and (p & (p - 1)) == 0


def syzygy_B(pair):
    """B = (p2/g)q1 + (p1/g)q2 over GF(2)[x] (None for degenerate diagonal templates)."""
    p1, q1 = _poly(pair[0], 'a'), _poly(pair[0], 'b')
    p2, q2 = _poly(pair[1], 'a'), _poly(pair[1], 'b')
    if p1 == 0 or p2 == 0:
        return None
    g = _gf2_gcd(p1, p2)
    return _gf2_mul(_gf2_div(p2, g), q1) ^ _gf2_mul(_gf2_div(p1, g), q2)


def block5_production_channels():
    print("-" * 92)
    print("BLOCK 5  GF(2)[x] production channels  [B monomial <=> the deg-1 channel can open]")
    print("-" * 92)
    # every k=4 deg-1 pure cycle (t_ell != 0) must have B monomial (production realized);
    # K3 (deg-3 branch) must have B non-monomial (production impossible at the parity level);
    # the flux pair has B monomial yet t_3 = 0: production open but sign-cancelled
    # (y_par = 1: word reversal is sign-reversing; the Cartier-Foata kill).
    N = 5
    n_deg1 = 0
    for pair in k4_pairs():
        ell, kind = anchor.effective_ell(N, pair)
        if ell == 0 or kind == 'diagonal-lift':
            continue
        H = anchor.build_H(N, pair)
        if sum_sq(t_moment(H, N, ell)) > 1e-9:
            B = syzygy_B(pair)
            assert B is not None and _is_monomial(B), \
                f"deg-1 cycle {pair}: B = {bin(B) if B else None} not monomial"
            n_deg1 += 1
    assert n_deg1 == 20, f"deg-1 cycles with monomial B: {n_deg1} != 20"
    print(f"  all 20 k=4 deg-1 pure cycles have B monomial (the channel is open and fires)  OK")
    B_k3 = syzygy_B([('X', 'X', 'Z'), ('X', 'Z', 'X')])
    assert B_k3 is not None and not _is_monomial(B_k3), "K3: B unexpectedly monomial"
    B_flux = syzygy_B([('I', 'X', 'Y'), ('X', 'I', 'Y')])
    assert B_flux is not None and _is_monomial(B_flux), "flux: B unexpectedly non-monomial"
    Hf = anchor.build_H(4, anchor.FLUX)
    assert sum_sq(t_moment(Hf, 4, 3)) == 0.0, "flux: t_3 != 0?!"
    print(f"  K3: B = {bin(B_k3)} non-monomial (no production, t_3 = 0 by parity);  "
          f"flux: B = {bin(B_flux)} monomial yet t_3 = 0 (production sign-cancelled, "
          f"y_par = 1 reversal)  OK")
    print("BLOCK 5 PASS")


def main():
    print("=" * 92)
    print("F87 GIRTH DICHOTOMY -- the deg-1 class in closed form (R-deg dissolves)")
    print("=" * 92)
    block1_factorization()
    block2_closed_form()
    block3_k4_battery()
    block4_k3_regression()
    block5_production_channels()
    print("=" * 92)
    print("ALL BLOCKS PASS")
    print("=" * 92)


if __name__ == "__main__":
    main()
