#!/usr/bin/env python3
r"""F87 residual-frontier census: three certified results at the R-deg / R-sign frontier.

Consolidated from the WIP scout `_f87_rwinding_followup.py` (since removed), with all shared
machinery imported from the committed anchor `f87_windowed_monomial_converse.py` (the windowed
converse monomial theorem: M(gamma) = A + gamma*Q recentred Liouvillian, two-reflection parity,
m* = 2*ell + deg). This script is self-validating: every part raises on failure, prints a PASS
line on success, and the process exits 0 only if all parts pass. It certifies the three results
the multi-agent review (2026-06-09) flagged as the residual frontier worth banking; result (3)
is the exploration data that feeds the still-open R-sign residual.

(1) FLUX threshold, letter-robust form. The naive single-leg return ((-iH)^k)_{ii} is 0 for ALL
    odd k for the flux pair IXY+XIY (signed cancellation around the unsigned 3-cycle, the §7.10
    XX+YY effect), yet the flux pair is hard with m*=9, deg=3. So the single-leg signed proxy is
    the WRONG necessary object once Q is interleaved: Q is diagonal but does not commute with A_L
    (unless H is diagonal), so the bra-return amplitude is Q-weighted in a ket-dependent way; the
    unweighted (H^k)_{ii} is only a sufficient-condition proxy (exact for K3, under-reads flux).
    The CORRECT threshold statement is on the TOTAL #A: a surviving word needs #A_L, #A_R >= ell
    (unsigned odd-girth path-existence, per word), hence #A >= 2*ell. Certified here by the full
    exact 3^9 word census at m*=9 (letters {A_L, A_R, Q}, exact CRT Gaussian-integer traces):
    every nonzero word has (#A_L, #A_R, #Q) all ODD and #A_L, #A_R >= 3; the only surviving #Q
    CLASS-SUM is #Q=3 (=> minimal surviving total-#A = 6 = 2*ell); and the surviving
    (#A_L, #A_R) split of the #Q=3 class is exactly (3,3), with class-sum equal to the anchor's
    exact monomial coefficient P_{9,3} (589824 flux, 2064384 K3 control).

(2) The 10 Perron mismatches. Cell-wide over the 50 hard pairs of the N=4 k=3 Z diagonal cell,
    the §7.5 reading "+N present and -N absent on Q|ker(A)" holds 40/50. The 10 failures are
    EXACTLY the diagonal-lift pairs (nonzero diagonal in H): their break is the lifted diagonal
    (the deg=1 / multi-Z closed form), NOT a missing -N reflection mode, so -N may survive there.
    All 16 pure-cycle pairs satisfy the skew. This is the result that justified the committed
    anchor's Block 7 restricting the Perron check to the 16 pure cycles: the two positivity
    mechanisms (Perron skew for pure cycles, diagonal lift for the rest) are disjoint and both
    hold.

(3) deg-3 block-skew table. For the 16 pure-cycle hard pairs: P_{m*,3} (the exact leading
    monomial coefficient, all positive) against the block odd power sums Tr(Qb^3), Tr(Qb^5) of
    Q restricted to ker(A). The top-skew Tr(Qb^3) > 0 holds 16/16, consistent with the +N Perron
    skew; the full table is printed as exploration data for the open R-sign residual (is
    P_{m*,3} expressible via the +N-skewed block?).

Run
---
    python simulations/f87_residual_frontier_census.py

Part (1) is the heavy step: 2 pairs x 3^9 = 19683 exact-CRT word traces each. The CRT prime
count is auto-sized from the rigorous trace bound |Tr(W)| <= d^2 * prod ||letter||_inf, so the
census arithmetic is exact by construction.
"""
from __future__ import annotations

import sys
import time
from collections import Counter
from itertools import combinations_with_replacement, product
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
import f87_windowed_monomial_converse as anchor  # noqa: E402


# ======================================================================
# PART 1 -- flux threshold: the letter-robust total-#A form, by full exact word census at m*.
# ======================================================================
def _census_nprimes(three, d2, mstar):
    """Minimal CRT prime count covering the rigorous per-word trace bound
    |Tr(W)| <= d2 * prod_i ||letter_i||_inf <= d2 * cmax^mstar (induced inf-norm = max abs row
    sum). The float64 modular matmul itself is exact while d2 * p^2 < 2^53 (anchor invariant)."""
    cmax = max(float(np.max(np.sum(np.abs(Re.astype(np.float64) + 1j * Im.astype(np.float64)),
                                   axis=1)))
               for (Re, Im) in three.values())
    bound = int(d2 * cmax ** mstar) + 1
    bank = anchor._CRT_PRIMES_LARGE if d2 <= 2 ** 11 else anchor._CRT_PRIMES_SMALL
    M = 1
    for k, p in enumerate(bank, start=1):
        M *= p
        if M // 2 > bound:
            return max(k, 2), bound
    raise AssertionError(f"prime bank too small for bound {bound}")


def census_threshold(label, N, pair, mstar, ell, want_coeff):
    """Full exact 3^mstar word census in {A_L, A_R, Q} at the first nonvanishing odd moment.
    Certifies (per word) all-odd counts + #A_L,#A_R >= ell, and (per class) that the only
    surviving #Q class-sum is #Q = mstar - 2*ell with (#A_L,#A_R) split exactly (ell,ell)."""
    A_L, A_R, Q, H = anchor.build_AL_AR_Q(N, pair)
    three = anchor.split_int(A_L, A_R, Q)
    d2 = A_L.shape[0]
    nprimes, bound = _census_nprimes(three, d2, mstar)
    is_complex = not np.allclose(H.imag, 0)
    print(f"\n  [{label}]  N={N}  pair={''.join(pair[0])}+{''.join(pair[1])}  ell={ell}  "
          f"m*={mstar}  complex-H={is_complex}")
    print(f"    exact census: 3^{mstar} = {3 ** mstar} words, nprimes={nprimes} "
          f"(CRT bound {bound})", flush=True)

    qclass_re, qclass_im = Counter(), Counter()
    split_re, split_im = Counter(), Counter()
    nonzero_by_q = Counter()
    t0 = time.time()
    n_seen = 0
    for word, re, im in anchor.census_word_traces(three, mstar, d2, nprimes=nprimes):
        n_seen += 1
        if n_seen % 4000 == 0:
            print(f"    ... {n_seen}/{3 ** mstar} words ({time.time() - t0:.0f} s)", flush=True)
        if re == 0 and im == 0:
            continue
        a, b, c = word.count('AL'), word.count('AR'), word.count('Q')
        # per-word parity (Block 3 of the anchor, RIGOROUS-GENERAL): all counts ODD
        assert a % 2 == 1 and b % 2 == 1 and c % 2 == 1, \
            f"{label}: nonzero word with non-all-odd counts ({a},{b},{c}): {word}"
        # per-word threshold (Block 4, unsigned odd-girth path-existence): both legs >= ell
        assert a >= ell and b >= ell, \
            f"{label}: nonzero word below threshold (#A_L,#A_R)=({a},{b}) < ell={ell}: {word}"
        qclass_re[c] += re
        qclass_im[c] += im
        nonzero_by_q[c] += 1
        if c == mstar - 2 * ell:
            split_re[(a, b)] += re
            split_im[(a, b)] += im
    assert n_seen == 3 ** mstar, f"{label}: census enumerated {n_seen} != 3^{mstar} words"
    print(f"    census done in {time.time() - t0:.0f} s; nonzero words by #Q: "
          f"{dict(sorted(nonzero_by_q.items()))}")

    # class-sum ledger: imaginary parts vanish exactly (Tr(M^m) is real)
    assert all(v == 0 for v in qclass_im.values()), f"{label}: imaginary #Q class-sum nonzero"
    assert all(v == 0 for v in split_im.values()), f"{label}: imaginary split class-sum nonzero"
    surviving = sorted(c for c, v in qclass_re.items() if v != 0)
    print(f"    surviving #Q class-sums at m*: {surviving} "
          f"(=> #A = {[mstar - c for c in surviving]})")
    # the only surviving class-sum is #Q = deg = mstar - 2*ell  =>  minimal total-#A = 2*ell
    deg = mstar - 2 * ell
    assert surviving == [deg], f"{label}: surviving #Q classes {surviving} != [{deg}]"
    min_total_A = mstar - max(surviving)
    assert min_total_A == 2 * ell, f"{label}: minimal surviving total-#A {min_total_A} != 2*ell"
    # the (#A_L,#A_R) split of the surviving class is exactly (ell,ell)
    nonzero_split = {k: v for k, v in split_re.items() if v != 0}
    print(f"    nonzero (#A_L,#A_R) class-sums in the #Q={deg} (#A={2 * ell}) class: "
          f"{nonzero_split}")
    assert set(nonzero_split) == {(ell, ell)}, \
        f"{label}: split classes {sorted(nonzero_split)} != [({ell},{ell})]"
    assert nonzero_split[(ell, ell)] == want_coeff, \
        f"{label}: ({ell},{ell}) class-sum {nonzero_split[(ell, ell)]} != P_{{m*,{deg}}} {want_coeff}"

    # cross-anchor pin: the census class-sum equals the exact monomial coefficient P_{m*,deg}
    Ar, Ai, Qi = anchor.build_integer_generators(N, pair)
    ms, re_co, im_co, _ = anchor.first_nonvanishing_odd(Ar, Ai, Qi, mstar)
    nz = [j for j, cc in enumerate(re_co) if cc != 0]
    assert ms == mstar and nz == [deg] and all(cc == 0 for cc in im_co), \
        f"{label}: anchor polynomial m*={ms}, powers {nz} (expected m*={mstar}, [{deg}])"
    assert re_co[deg] == want_coeff, f"{label}: P_{{m*,{deg}}} {re_co[deg]} != {want_coeff}"
    assert qclass_re[deg] == re_co[deg], \
        f"{label}: census class-sum {qclass_re[deg]} != exact P_{{m*,{deg}}} {re_co[deg]}"
    print(f"    >>> minimal total-#A = {2 * ell} = 2*ell; split = ({ell},{ell}) only; "
          f"class-sum = exact P_{{{mstar},{deg}}} = {want_coeff}  OK")


def part1_flux_threshold():
    print("=" * 92)
    print("PART 1  FLUX threshold: min total-#A = 2*ell, (#A_L,#A_R) split = (ell,ell)  "
          "[exact 3^9 census]")
    print("=" * 92)
    # the flux witness first: the signed single-leg diagonal cancels ((H^3)_ii = 0) while the
    # unsigned 3-walk exists, so only the total-#A form of the threshold is letter-robust
    Hf = anchor.build_H(4, anchor.FLUX)
    sgn3 = float(np.max(np.abs(np.diag(np.linalg.matrix_power(Hf, 3)))))
    uns3 = float(np.max(np.diag(np.linalg.matrix_power(anchor.adjacency(Hf), 3))))
    assert sgn3 < 1e-9 < uns3, f"flux witness: signed (H^3)_ii={sgn3}, unsigned={uns3}"
    print(f"  flux single-leg witness: signed (H^3)_ii = {sgn3:.1f} (cancels) but unsigned "
          f"(|H|^3)_ii = {uns3:.0f} > 0")
    census_threshold("flux IXY+XIY", 4, anchor.FLUX, mstar=9, ell=3, want_coeff=589824)
    census_threshold("K3 XXZ+XZX (control)", 4, anchor.K3_EVEN, mstar=9, ell=3,
                     want_coeff=2064384)
    print("PART 1 PASS")


# ======================================================================
# PART 2 -- the 10 Perron mismatches are exactly the diagonal-lift pairs.
# ======================================================================
def _hard_pairs(chain, N, k):
    terms = [''.join(t) for t in product('IXYZ', repeat=k)
             if not all(L == 'I' for L in t) and anchor.klein_index(''.join(t)) == (0, 1)]
    for t1, t2 in combinations_with_replacement(terms, 2):
        if anchor.y_parity(t1) != anchor.y_parity(t2):
            continue
        if fw.classify_pauli_pair(chain, [tuple(t1), tuple(t2)], dephase_letter='Z') != 'hard':
            continue
        yield t1, t2


def part2_perron_mismatches(N=4, k=3):
    print("\n" + "=" * 92)
    print("PART 2  Perron reading cell-wide: the mismatches of '+N present & -N absent' "
          "on Q|ker(A)")
    print("=" * 92)
    chain = fw.ChainSystem(N=N)
    mismatches, matches = [], []
    for t1, t2 in _hard_pairs(chain, N, k):
        H = anchor.build_H(N, [tuple(t1), tuple(t2)])
        diag = anchor.has_nonzero_diagonal(H)
        Qb, _ = anchor.Q_on_kerA_block(H, N)
        ev = np.linalg.eigvals(Qb).real
        has_pN = bool(np.any(np.abs(ev - N) < 1e-6))
        has_mN = bool(np.any(np.abs(ev + N) < 1e-6))
        rec = (t1, t2, diag, has_pN, has_mN)
        (matches if (has_pN and not has_mN) else mismatches).append(rec)
    n_hard = len(matches) + len(mismatches)
    n_diag_mis = sum(1 for r in mismatches if r[2])
    n_cycle_match = sum(1 for r in matches if not r[2])
    print(f"  hard pairs: {n_hard}")
    print(f"  '+N present & -N absent' (the §7.5 skew): {len(matches)}")
    print(f"  mismatches (skew reading fails): {len(mismatches)}")
    print(f"    of the mismatches, DIAGONAL-LIFT (nonzero diag in H): "
          f"{n_diag_mis}/{len(mismatches)}")
    print(f"    of the skew-matches, NON-diagonal (pure cycle): "
          f"{n_cycle_match}/{len(matches)}")
    print(f"  the 10 mismatches (t1,t2,diag,+N,-N):")
    for r in mismatches:
        print(f"     {r[0]}+{r[1]}  diag={r[2]}  +N={r[3]}  -N={r[4]}")
    # the pinned numbers (recorded by the 2026-06-09 scout run)
    assert n_hard == 50, f"expected 50 hard pairs, got {n_hard}"
    assert len(matches) == 40, f"Perron skew holds {len(matches)}/50 (expected 40)"
    assert len(mismatches) == 10, f"expected 10 mismatches, got {len(mismatches)}"
    assert n_diag_mis == 10, \
        f"only {n_diag_mis}/10 mismatches are diagonal-lift (claim: ALL of them)"
    assert n_cycle_match == 16, \
        f"pure-cycle skew-matches {n_cycle_match} != 16 (claim: all 16 pure cycles match)"
    # so: the skew reading is the PURE-CYCLE positivity mechanism (16/16); the diagonal-lift
    # pairs break via the lifted diagonal (deg=1 / multi-Z closed form), and -N may survive there
    print("  >>> mismatches are exactly diagonal-lift pairs; all 16 pure cycles carry the skew")
    print("PART 2 PASS")


# ======================================================================
# PART 3 -- deg-3 block-skew table over the 16 pure-cycle hard pairs (open R-sign data).
# ======================================================================
def part3_deg3_block_skew(N=4, k=3):
    print("\n" + "=" * 92)
    print("PART 3  deg=3 pure-cycle hard pairs: P_(m*,3) sign vs Q|ker(A) block "
          "odd-power-sum skew")
    print("=" * 92)
    chain = fw.ChainSystem(N=N)
    rows = []
    for t1, t2 in _hard_pairs(chain, N, k):
        pair = [tuple(t1), tuple(t2)]
        H = anchor.build_H(N, pair)
        if anchor.has_nonzero_diagonal(H):
            continue  # pure-cycle only
        Ar, Ai, Qm = anchor.build_integer_generators(N, pair)
        # every pure cycle in this cell sits at m* = 9 (Block 7 of the anchor); max_m=11 finds
        # the first nonvanishing odd moment with margin, nprimes=6 per the anchor's CRT bound
        mstar, re_co, im_co, _ = anchor.first_nonvanishing_odd(Ar, Ai, Qm, 11, nprimes=6)
        nz = [j for j, c in enumerate(re_co) if c != 0]
        assert len(nz) == 1, f"{t1}+{t2}: first moment not a monomial (powers {nz})"
        deg, coeff = nz[0], re_co[nz[0]]
        Qb, _ = anchor.Q_on_kerA_block(H, N)
        ev = np.linalg.eigvals(Qb).real
        s3 = float(np.sum(ev ** 3))
        s5 = float(np.sum(ev ** 5))
        rows.append((t1, t2, mstar, deg, coeff, s3, s5))
    n = len(rows)
    n_deg3 = sum(1 for r in rows if r[3] == 3)
    n_pos = sum(1 for r in rows if r[4] > 0)
    n_s3pos = sum(1 for r in rows if r[5] > 1e-9)
    print(f"  pure-cycle hard pairs: {n}  (deg=3: {n_deg3})")
    print(f"  P_(m*,deg) > 0: {n_pos}/{n}")
    print(f"  block odd-power-sum Tr(Qb^3) > 0 (top-skew): {n_s3pos}/{n}")
    print(f"  {'t1':>6} {'t2':>6} {'m*':>3} {'deg':>3} {'P_(m*,deg)':>12} {'Tr(Qb^3)':>10} "
          f"{'Tr(Qb^5)':>10}")
    for (t1, t2, mstar, deg, coeff, s3, s5) in sorted(rows, key=lambda r: (r[2], r[0])):
        print(f"  {t1:>6} {t2:>6} {mstar:>3} {deg:>3} {str(coeff):>12} {s3:>10.2f} {s5:>10.2f}")
    # the pinned numbers (recorded by the 2026-06-09 scout run)
    assert n == 16, f"expected 16 pure-cycle hard pairs, got {n}"
    assert n_deg3 == 16, f"deg=3 only {n_deg3}/16"
    assert all(r[2] == 9 for r in rows), "a pure cycle sits off m*=9 (ell=3 cell)"
    assert n_pos == 16, f"P_(m*,3) > 0 only {n_pos}/16"
    assert n_s3pos == 16, f"Tr(Qb^3) > 0 only {n_s3pos}/16"
    print("  >>> all 16: positive monomial at m*=9 AND positive block top-skew Tr(Qb^3) "
          "(R-sign still open: no closed form ties them yet)")
    print("PART 3 PASS")


def main():
    t0 = time.time()
    print("=" * 92)
    print("F87 RESIDUAL-FRONTIER CENSUS -- self-validating (consolidated from "
          "_f87_rwinding_followup.py)")
    print("=" * 92)
    part1_flux_threshold()
    part2_perron_mismatches(N=4, k=3)
    part3_deg3_block_skew(N=4, k=3)
    print("\n" + "=" * 92)
    print(f"ALL PARTS PASS  ({time.time() - t0:.0f} s)")
    print("=" * 92)


if __name__ == "__main__":
    main()
