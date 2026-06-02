#!/usr/bin/env python3
"""Klein-routing of the palindrome's hidden symmetry Q across two-term bilinears.

THE_OTHER_SIDE (Q6/Q7) asks: when the canonical palindrome operator Π (P1)
fails to pair the Liouvillian spectrum, is there a different symmetry Q that
still closes it, and does the Klein-Vierergruppe index of the two bond
bilinears predict which of the three fates a Hamiltonian meets?

The three fates of a Z-dephased XY-type chain H = Σ_bonds (t1 + t2):
  truly : the canonical Π (P1) already pairs the spectrum (M = 0).
  soft  : Π fails (M ≠ 0) but the palindrome still holds; a hidden Q ≠ Π
          closes the spectrum-pairing. These are the benign ones.
  hard  : the palindrome is BROKEN; no Q closes it.

This probe is the data foundation: it classifies all C(9,2) = 36 unordered
two-term bond-bilinear Hamiltonians at N ∈ {3, 4, 5} under Z-dephasing, records
each term's Klein index (bit_a, bit_b), and an independent parity-break test
(does H commute with the global X-flip X^⊗N). The N=3 layer is self-validated
against the verbatim March hard-set (14 combos) and the parity counts (26, 12).

Conventions: work in LETTER strings ('X', 'Y', 'Z') throughout; never raw
indices, since standalone scripts use I,X,Y,Z = 0,1,2,3 while the framework
uses Klein bits (Y → (1,1), Z → (0,1)). All framework calls go through the
letter-tuple API.

Reuses the framework verbatim:
  fw.classify_pauli_pair(chain, [t1, t2])  → trichotomy class
  fw.klein_index('XY')                      → (bit_a, bit_b)
  fw.site_op(N, site, 'X')                  → single-site operator
  fw.ChainSystem(N, topology='chain')       → the chain (default γ₀ = 0.05)
"""
from __future__ import annotations

import sys
from itertools import combinations

import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
# The 9 traceless bond bilinears in canonical H_LABELS order (A,B for A,B in XYZ).
H_LABELS = [a + b for a in 'XYZ' for b in 'XYZ']
assert H_LABELS == ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']

# The C(9,2) = 36 unordered pairs, each a canonical 'A+B' label with the two
# bilinears in H_LABELS order (so 'XX+XY', never 'XY+XX').
COMBOS = list(combinations(H_LABELS, 2))
assert len(COMBOS) == 36

N_VALUES = [3, 4, 5]

PARITY_TOL = 1e-9

# The verbatim March N=3 hard list (14 combos): palindrome BROKEN.
MARCH_HARD_N3 = {
    'XX+XY', 'XX+YX', 'XY+XZ', 'XY+YY', 'XY+YZ', 'XY+ZX', 'XY+ZY',
    'XZ+YX', 'XZ+ZY', 'YX+YY', 'YX+YZ', 'YX+ZX', 'YX+ZY', 'YZ+ZX',
}


def combo_label(bil_a, bil_b):
    """Canonical 'A+B' label: the two bilinears ordered by H_LABELS index."""
    ia, ib = H_LABELS.index(bil_a), H_LABELS.index(bil_b)
    if ia <= ib:
        return f"{bil_a}+{bil_b}"
    return f"{bil_b}+{bil_a}"


def letters(bilinear):
    """'XY' → ('X', 'Y')."""
    return (bilinear[0], bilinear[1])


# ---------------------------------------------------------------------------
# Independent parity-break test (does NOT go through classify_pauli_pair)
# ---------------------------------------------------------------------------
def parity_break(N, bil_a, bil_b):
    """True iff H = Σ_bonds (t1 + t2) fails to commute with X^⊗N.

    H is built directly on the open chain bonds (k, k+1) via fw.site_op, the
    independent route, so this is a genuine cross-check on the framework's
    own bilinear builder inside classify_pauli_pair.

    X^⊗N is the global spin-flip; ‖[H, X^N]‖ > tol signals that the X-parity
    Z₂ symmetry is broken by this term pair.
    """
    d = 2 ** N
    bonds = [(k, k + 1) for k in range(N - 1)]
    H = np.zeros((d, d), dtype=complex)
    for (a_letter, b_letter) in (letters(bil_a), letters(bil_b)):
        for (i, j) in bonds:
            H = H + fw.site_op(N, i, a_letter) @ fw.site_op(N, j, b_letter)

    Xn = fw.site_op(N, 0, 'X')
    for site in range(1, N):
        Xn = Xn @ fw.site_op(N, site, 'X')

    comm = H @ Xn - Xn @ H
    return float(np.linalg.norm(comm)) > PARITY_TOL


# ---------------------------------------------------------------------------
# Per-N classification table
# ---------------------------------------------------------------------------
def classify_all(N):
    """Classify all 36 combos at chain length N. Returns list of row dicts."""
    chain = fw.ChainSystem(N, topology='chain')
    rows = []
    for bil_a, bil_b in COMBOS:
        label = combo_label(bil_a, bil_b)
        ka = fw.klein_index(bil_a)
        kb = fw.klein_index(bil_b)
        pbreak = parity_break(N, bil_a, bil_b)
        cls = fw.classify_pauli_pair(chain, [letters(bil_a), letters(bil_b)])
        rows.append({
            'label': label,
            'klein_a': ka,
            'klein_b': kb,
            'parity_break': pbreak,
            'class': cls,
        })
    return rows


def summarize(rows):
    """Return (class_counts, parity_break_count, hard_set)."""
    class_counts = {'truly': 0, 'soft': 0, 'hard': 0}
    pbreak_count = 0
    hard_set = set()
    for r in rows:
        class_counts[r['class']] += 1
        if r['parity_break']:
            pbreak_count += 1
        if r['class'] == 'hard':
            hard_set.add(r['label'])
    return class_counts, pbreak_count, hard_set


def print_table(N, rows):
    kl = fw.KLEIN_LABELS
    print(f"  {'combo':>8s}  {'klein(t1)':>11s}  {'klein(t2)':>11s}  "
          f"{'parity_break':>12s}  {'class':>6s}")
    print(f"  {'-' * 8}  {'-' * 11}  {'-' * 11}  {'-' * 12}  {'-' * 6}")
    for r in rows:
        ka, kb = r['klein_a'], r['klein_b']
        ka_str = f"{ka} {kl[ka]}"
        kb_str = f"{kb} {kl[kb]}"
        print(f"  {r['label']:>8s}  {ka_str:>11s}  {kb_str:>11s}  "
              f"{str(r['parity_break']):>12s}  {r['class']:>6s}")


def print_summary(N, rows):
    class_counts, pbreak_count, hard_set = summarize(rows)
    benign_breakers = sum(
        1 for r in rows if r['parity_break'] and r['class'] != 'hard'
    )
    hard_break_consistent = all(
        r['parity_break'] for r in rows if r['class'] == 'hard'
    )
    print()
    print(f"  N={N} counts: truly={class_counts['truly']}, "
          f"soft={class_counts['soft']}, hard={class_counts['hard']}")
    print(f"  N={N} parity_break: {pbreak_count} / 36  "
          f"(benign breakers [break ∧ ¬hard]: {benign_breakers})")
    print(f"  N={N} every hard combo breaks parity: {hard_break_consistent}")
    print(f"  N={N} hard set ({len(hard_set)}):")
    for label in sorted(hard_set):
        print(f"      {label}")
    return class_counts, pbreak_count, hard_set, benign_breakers


# ---------------------------------------------------------------------------
# N=3 validation (harness self-check against known March data)
# ---------------------------------------------------------------------------
def validate_n3(rows):
    class_counts, pbreak_count, hard_set = summarize(rows)
    benign_breakers = sum(
        1 for r in rows if r['parity_break'] and r['class'] != 'hard'
    )

    # 1. Exact hard set matches the verbatim March list of 14.
    assert hard_set == MARCH_HARD_N3, (
        f"N=3 hard set mismatch.\n"
        f"  missing (expected, not found): {sorted(MARCH_HARD_N3 - hard_set)}\n"
        f"  extra (found, not expected):   {sorted(hard_set - MARCH_HARD_N3)}"
    )
    assert len(hard_set) == 14, f"N=3 hard count {len(hard_set)} != 14"

    # 2. Exactly 26 combos break X-parity.
    assert pbreak_count == 26, f"N=3 parity_break count {pbreak_count} != 26"

    # 3. Exactly 12 benign parity-breakers (break ∧ not hard).
    assert benign_breakers == 12, (
        f"N=3 benign parity-breakers {benign_breakers} != 12"
    )

    # 4. Strict containment: every hard combo breaks parity.
    for r in rows:
        if r['class'] == 'hard':
            assert r['parity_break'], (
                f"N=3 hard combo {r['label']} does NOT break parity "
                f"(containment violated)"
            )

    print()
    print("  N=3 ASSERTIONS PASSED:")
    print("    [1] hard set == verbatim March list (14 combos)")
    print("    [2] count(parity_break) == 26")
    print("    [3] count(parity_break ∧ class != hard) == 12")
    print("    [4] every hard combo has parity_break == True (containment)")


# ===========================================================================
# THE ROUTING RULE (R2b): predict the fate from letters alone, no eigensolve
# ===========================================================================
# Z-dephasing splits each site into a DC bus {I, Z} (immune) and an AC bus
# {X, Y} (damped). The palindrome operator Q routes the AC bus, and per-site it
# comes in two crossover families: P1 routes the X-channel, P4 routes the
# Y-channel. A single uniform Q = M ⊗ M ⊗ ... must serve BOTH bond terms at
# once, so a Hamiltonian is palindromic iff its two terms share a valid uniform
# Q-family; failing that, it escapes via an alternating (site-parity) or a
# non-local (entangled) Q; failing THAT, the palindrome is broken (hard).
#
# fam(t) = the set of uniform Q-families that close term t alone:
#   fam(XX) = fam(YY) = {P1, P4}   both crossovers route a same-letter bilinear
#   fam(ZZ) = {all}                ZZ is dephasing-aligned, every Q fixes it
#   fam(YZ) = fam(ZY) = {P1}       only the X-router survives the Y·Z mix
#   fam(XZ) = fam(ZX) = {P4}       only the Y-router survives the X·Z mix
#   fam(XY) = fam(YX) = {}         no single uniform Q routes both X and Y
# This reproduces classify_pauli_pair bit-exactly across N = 3, 4, 5.
ALL_Q = frozenset({'P1', 'P4', 'M2'})
FAM = {
    'XX': frozenset({'P1', 'P4'}),
    'YY': frozenset({'P1', 'P4'}),
    'ZZ': ALL_Q,
    'YZ': frozenset({'P1'}),
    'ZY': frozenset({'P1'}),
    'XZ': frozenset({'P4'}),
    'ZX': frozenset({'P4'}),
    'XY': frozenset(),
    'YX': frozenset(),
}

# The two non-local escapes: same-site X&Y collision sitting over a shared dark
# Z port, rescued by an entangled (non-product) Q.
NON_LOCAL_PAIRS = {frozenset({'XZ', 'YZ'}), frozenset({'ZX', 'ZY'})}

MOTHER = {'XX', 'YY', 'ZZ'}


def route(bil_a, bil_b):
    """Predict the fate of H = Σ_bonds (bil_a + bil_b) from the letters alone.

    Implements the R2b routing rule (verified 36/36 against classify_pauli_pair
    at N = 3, 4, 5). No Liouvillian is built; the verdict comes from each term's
    valid uniform-Q-family set FAM[·].

    Returns a dict {fate, family, reason}:
      fate   ∈ {'truly', 'soft', 'hard'}
      family ∈ {'P1', 'uniform', 'alternating', 'non_local', None}
      reason : short human-readable string.
    """
    fa, fb = FAM[bil_a], FAM[bil_b]

    # 1. Shared uniform Q-family: the spectrum pairs under one product Q.
    if fa & fb:
        if bil_a in MOTHER and bil_b in MOTHER:
            return {'fate': 'truly', 'family': 'P1',
                    'reason': 'both Mother (XX/YY/ZZ); canonical Π = P1 pairs spectrum'}
        return {'fate': 'soft', 'family': 'uniform',
                'reason': f'shared uniform Q-family {{{", ".join(sorted(fa & fb))}}}'}

    # 2. Alternating (site-parity) escape: an XY/YX term plus an XY/YX or ZZ
    #    partner closes via P1 ⊗ M2 ⊗ P1 ⊗ ... (alternating per-site crossover).
    if ((bil_a in {'XY', 'YX'} or bil_b in {'XY', 'YX'})
            and bil_a in {'XY', 'YX', 'ZZ'} and bil_b in {'XY', 'YX', 'ZZ'}):
        return {'fate': 'soft', 'family': 'alternating',
                'reason': 'XY/YX with site-parity (alternating) Q'}

    # 3. Non-local (entangled) escape: same-site X&Y over a shared dark Z port.
    if frozenset({bil_a, bil_b}) in NON_LOCAL_PAIRS:
        return {'fate': 'soft', 'family': 'non_local',
                'reason': 'same-site X&Y collision over shared dark Z; entangled Q'}

    # 4. No Q closes it: the palindrome is broken.
    return {'fate': 'hard', 'family': None,
            'reason': 'no shared uniform Q, no alternating/non-local escape'}


# ===========================================================================
# THE SECOND CONDITION (R2c): WHY a hard combo is hard (structural reason)
# ===========================================================================
# Beyond breaking X-parity, a combo is hard iff it carries an irreducible,
# unroutable same-qubit X/Y demand. Read through the Klein cells of the two
# terms (M, F_a, C, F_b), the discriminator is a cell-pair rule with one
# refinement per split cell. Reproduces the 14 hard / 19 soft / 3 truly
# bit-exactly at N = 3, 4, 5.
#   M   = (0,0) = {XX, YY, ZZ}   Mother
#   F_a = (0,1) = {XY, YX}       Y-Father  (Π²-odd)
#   C   = (1,0) = {YZ, ZY}       Child     (dephase-aligned)
#   F_b = (1,1) = {XZ, ZX}       Z-Father  (Π²-odd)
KLEIN_CELL = {(0, 0): 'M', (0, 1): 'F_a', (1, 0): 'C', (1, 1): 'F_b'}


def _klein_index(bilinear):
    """(bit_a, bit_b) of a 2-letter bilinear: per-letter Klein bits summed mod 2."""
    a0, b0 = fw.LABEL_TO_INDEX[bilinear[0]]
    a1, b1 = fw.LABEL_TO_INDEX[bilinear[1]]
    return ((a0 + a1) % 2, (b0 + b1) % 2)


def _cell_pair(bil_a, bil_b):
    """The unordered Klein-cell pair of the two terms, e.g. ('F_a', 'M')."""
    return tuple(sorted([KLEIN_CELL[_klein_index(bil_a)],
                         KLEIN_CELL[_klein_index(bil_b)]]))


def _site_letters(bil_a, bil_b):
    """Letters landing on (site0, site1), one from each bond term."""
    return (bil_a[0], bil_b[0]), (bil_a[1], bil_b[1])


def _same_site_xy(bil_a, bil_b):
    """True iff one single qubit is asked to carry BOTH an X and a Y (the
    diplexer conflict: that site would need to be both P1 and P4 at once)."""
    s0, s1 = _site_letters(bil_a, bil_b)
    return ({'X', 'Y'} <= set(s0)) or ({'X', 'Y'} <= set(s1))


def q7_reason(bil_a, bil_b):
    """The structural reason a combo is (or is not) hard: the second condition.

    Returns a dict {hard, cell_pair, reason}:
      hard      : bool, True iff an irreducible unroutable same-qubit X/Y demand.
      cell_pair : the unordered Klein-cell pair, e.g. ('C', 'F_b').
      reason    : short string naming the conflict type.

    The hard-by-reason set equals the 14 hard combos bit-exactly at N = 3, 4, 5,
    with the lit-vs-dark ({F_a, M}) and same-site/different-site ({C, F_b})
    discriminations resolved.
    """
    k = _cell_pair(bil_a, bil_b)

    # Always-hard cells: two Fathers, or Child + Y-Father.
    if k == ('F_a', 'F_b'):
        return {'hard': True, 'cell_pair': k,
                'reason': 'F_a+F_b (Y-Father × Z-Father): both X and Y demanded, unroutable'}
    if k == ('C', 'F_a'):
        return {'hard': True, 'cell_pair': k,
                'reason': 'C+F_a (Child × Y-Father): irreducible same-qubit X/Y demand'}

    # Split cell A: Y-Father (XY/YX) + Mother (XX/YY/ZZ).
    #   lit Mother  (X or Y letter collides on one site)        -> hard
    #   dark Mother (ZZ; Z is dephasing-aligned, no X/Y demand) -> soft
    if k == ('F_a', 'M'):
        if _same_site_xy(bil_a, bil_b):
            return {'hard': True, 'cell_pair': k,
                    'reason': 'F_a+M lit Mother: X and Y forced onto one qubit'}
        return {'hard': False, 'cell_pair': k,
                'reason': 'F_a+M dark Mother (ZZ): dephasing-aligned, no X/Y conflict'}

    # Split cell B: Z-Father (XZ/ZX) + Child (YZ/ZY).
    #   X,Y on the SAME site over a shared dark Z -> entangled Q rescues -> soft
    #   X,Y on DIFFERENT sites                    -> no Q exists         -> hard
    if k == ('C', 'F_b'):
        if frozenset({bil_a, bil_b}) in NON_LOCAL_PAIRS:
            return {'hard': False, 'cell_pair': k,
                    'reason': 'C+F_b same-site X&Y over shared dark Z: entangled (non-local) rescue'}
        return {'hard': True, 'cell_pair': k,
                'reason': 'C+F_b lit X,Y on different sites: no Q closes it'}

    # All other cells (M+M, M+C, M+F_b, C+C, F_a+F_a, ...): never hard.
    return {'hard': False, 'cell_pair': k,
            'reason': 'no irreducible same-qubit X/Y conflict'}


# ---------------------------------------------------------------------------
# Routing-rule validation phase (R2b + R2c): no eigensolve, letters only
# ---------------------------------------------------------------------------
# Family ground truth (R2a, bit-exact Q construction). The 22 palindromic
# combos partitioned by the family of their hidden Q; the 14 hard have none.
ROUTE_FAMILY_GT = {
    'P1': {'XX+YY', 'XX+ZZ', 'YY+ZZ'},
    'uniform': {
        'XX+XZ', 'XX+YZ', 'XX+ZX', 'XX+ZY', 'XZ+YY', 'XZ+ZX', 'XZ+ZZ',
        'YY+YZ', 'YY+ZX', 'YY+ZY', 'YZ+ZY', 'YZ+ZZ', 'ZX+ZZ', 'ZY+ZZ',
    },
    'alternating': {'XY+YX', 'XY+ZZ', 'YX+ZZ'},
    'non_local': {'XZ+YZ', 'ZX+ZY'},
}


def validate_routing(N_values):
    """Assert route() == classify_pauli_pair fate for all 36 combos at every N
    (0 mismatches), and that route()'s family matches the R2a ground truth for
    the 22 palindromic combos. Returns the per-N mismatch counts (all 0)."""
    print()
    print("=" * 78)
    print("ROUTING RULE (R2b): route() vs classify_pauli_pair, letters only")
    print("=" * 78)

    mismatch_counts = {}
    for N in N_values:
        chain = fw.ChainSystem(N, topology='chain')
        mismatches = 0
        for bil_a, bil_b in COMBOS:
            r = route(bil_a, bil_b)
            cls = fw.classify_pauli_pair(chain, [letters(bil_a), letters(bil_b)])
            if r['fate'] != cls:
                mismatches += 1
                print(f"    N={N} MISMATCH {combo_label(bil_a, bil_b):>8s}: "
                      f"route={r['fate']} classify={cls}")
        mismatch_counts[N] = mismatches
        print(f"  N={N}: route vs classify_pauli_pair mismatches = {mismatches} / 36")

    # ASSERT: 0 mismatches at every N.
    for N, m in mismatch_counts.items():
        assert m == 0, f"N={N} route disagrees with classify_pauli_pair on {m} combo(s)"

    # ASSERT: family of the 22 palindromic combos matches the R2a ground truth.
    family_by_label = {}
    for bil_a, bil_b in COMBOS:
        family_by_label[combo_label(bil_a, bil_b)] = route(bil_a, bil_b)['family']
    for fam_name, labels in ROUTE_FAMILY_GT.items():
        for label in labels:
            assert family_by_label[label] == fam_name, (
                f"family mismatch {label}: route says {family_by_label[label]!r}, "
                f"ground truth {fam_name!r}"
            )

    # ASSERT: exactly the 14 hard combos carry family == None.
    none_family = {lab for lab, f in family_by_label.items() if f is None}
    assert none_family == MARCH_HARD_N3, (
        f"family==None set {sorted(none_family)} != 14 hard"
    )

    # Tally from route().
    fates = [route(a, b)['fate'] for a, b in COMBOS]
    fams = [route(a, b)['family'] for a, b in COMBOS]
    tally_fate = {f: fates.count(f) for f in ('truly', 'soft', 'hard')}
    tally_fam = {f: fams.count(f) for f in ('P1', 'uniform', 'alternating', 'non_local')}

    print()
    print("  ROUTING ASSERTIONS PASSED:")
    print("    [R2b] route fate == classify_pauli_pair for all 36 combos, all N "
          "(0 mismatches)")
    print("    [R2a] route family matches ground truth for all 22 palindromic combos")
    print(f"    fate tally:   {tally_fate}")
    print(f"    family tally: {tally_fam}  (P1=truly, rest=soft; hard has no family)")
    return mismatch_counts, tally_fate, tally_fam


def validate_q7_reason(N_values):
    """Assert q7_reason()'s hard-by-reason set equals the 14 hard combos, with
    the lit-vs-dark and same-site/different-site discriminations correct, and
    that it agrees with classify_pauli_pair on every combo at every N."""
    print()
    print("=" * 78)
    print("SECOND CONDITION (R2c): q7_reason() structural reason a combo is hard")
    print("=" * 78)

    # Hard-by-reason set (N-independent: pure letter/Klein-cell logic).
    hard_by_reason = set()
    for bil_a, bil_b in COMBOS:
        if q7_reason(bil_a, bil_b)['hard']:
            hard_by_reason.add(combo_label(bil_a, bil_b))

    print(f"  hard-by-reason set ({len(hard_by_reason)}):")
    for label in sorted(hard_by_reason):
        # Re-derive the term order for the reason string.
        a, b = label.split('+')
        info = q7_reason(a, b)
        print(f"      {label:>8s}  cell={str(info['cell_pair']):14s}  {info['reason']}")

    # ASSERT: hard-by-reason == the 14 hard combos, bit-exact.
    assert hard_by_reason == MARCH_HARD_N3, (
        f"q7_reason hard set mismatch.\n"
        f"  missing: {sorted(MARCH_HARD_N3 - hard_by_reason)}\n"
        f"  extra:   {sorted(hard_by_reason - MARCH_HARD_N3)}"
    )
    assert len(hard_by_reason) == 14, f"hard-by-reason count {len(hard_by_reason)} != 14"

    # ASSERT: the named discriminations resolve correctly.
    #   lit-vs-dark on cell {F_a, M}:
    assert q7_reason('XY', 'YY')['hard'] is True, "XY+YY (lit Y) should be hard"
    assert q7_reason('XY', 'ZZ')['hard'] is False, "XY+ZZ (dark Z) should be soft"
    assert q7_reason('XX', 'XY')['hard'] is True, "XX+XY (lit X) should be hard"
    #   same-site vs different-site on cell {C, F_b}:
    assert q7_reason('XZ', 'ZY')['hard'] is True, "XZ+ZY (X,Y diff sites) should be hard"
    assert q7_reason('YZ', 'ZX')['hard'] is True, "YZ+ZX (X,Y diff sites) should be hard"
    assert q7_reason('XZ', 'YZ')['hard'] is False, "XZ+YZ (X,Y same site) should be soft"
    assert q7_reason('ZX', 'ZY')['hard'] is False, "ZX+ZY (X,Y same site) should be soft"

    # ASSERT: hard-by-reason agrees with classify_pauli_pair at every N.
    for N in N_values:
        chain = fw.ChainSystem(N, topology='chain')
        for bil_a, bil_b in COMBOS:
            is_hard = q7_reason(bil_a, bil_b)['hard']
            cls_hard = (fw.classify_pauli_pair(
                chain, [letters(bil_a), letters(bil_b)]) == 'hard')
            assert is_hard == cls_hard, (
                f"N={N} q7_reason disagrees on {combo_label(bil_a, bil_b)}: "
                f"reason_hard={is_hard} classify_hard={cls_hard}"
            )

    print()
    print("  SECOND-CONDITION ASSERTIONS PASSED:")
    print("    [R2c] hard-by-reason set == 14 hard combos (bit-exact)")
    print("    lit (XY+YY) hard vs dark (XY+ZZ) soft  [cell {F_a, M}]")
    print("    diff-site (XZ+ZY, YZ+ZX) hard vs same-site (XZ+YZ, ZX+ZY) soft  "
          "[cell {C, F_b}]")
    print(f"    hard-by-reason agrees with classify_pauli_pair at all N ∈ {N_values}")
    return hard_by_reason


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Klein-routing of the palindrome's hidden symmetry Q (THE_OTHER_SIDE Q6/Q7)")
    print(f"  36 two-term bond-bilinear Hamiltonians, Z-dephasing, γ₀ = 0.05")
    print(f"  N ∈ {N_VALUES}; bilinears (H_LABELS order): {H_LABELS}")
    print(f"  class: truly = canonical Π pairs spectrum;  soft = hidden Q ≠ Π pairs;  "
          f"hard = palindrome broken")
    print()

    results = {}
    for N in N_VALUES:
        print("=" * 78)
        print(f"N = {N}")
        print("=" * 78)
        rows = classify_all(N)
        print_table(N, rows)
        summary = print_summary(N, rows)
        results[N] = (rows, summary)

        if N == 3:
            validate_n3(rows)

    # ----- Cross-N comparison (N=3 is the validated reference) -----
    print()
    print("=" * 78)
    print("CROSS-N COMPARISON (is N=3 a special case?)")
    print("=" * 78)
    ref_counts, ref_pbreak, ref_hard, ref_benign = results[3][1]
    print(f"  Reference N=3: counts truly/soft/hard = "
          f"{ref_counts['truly']}/{ref_counts['soft']}/{ref_counts['hard']}, "
          f"parity_break = {ref_pbreak}, benign = {ref_benign}")
    print()
    any_diff = False
    for N in N_VALUES:
        if N == 3:
            continue
        counts, pbreak, hard, benign = results[N][1]
        same_hard = (hard == ref_hard)
        same_counts = (counts == ref_counts)
        same_pbreak = (pbreak == ref_pbreak)
        differs = not (same_hard and same_counts and same_pbreak)
        any_diff = any_diff or differs
        tag = "DIFFERS from N=3" if differs else "matches N=3"
        print(f"  N={N}: counts truly/soft/hard = "
              f"{counts['truly']}/{counts['soft']}/{counts['hard']}, "
              f"parity_break = {pbreak}, benign = {benign}  -> {tag}")
        if not same_hard:
            added = sorted(hard - ref_hard)
            removed = sorted(ref_hard - hard)
            if added:
                print(f"      hard ADDED vs N=3 ({len(added)}): {added}")
            if removed:
                print(f"      hard REMOVED vs N=3 ({len(removed)}): {removed}")

    print()
    if any_diff:
        print("  FINDING: the trichotomy / parity structure is NOT N-invariant.")
        print("  N=3 is a special case; the Klein-routing of Q changes with chain")
        print("  length. See the per-N hard sets above for the exact migration.")
    else:
        print("  FINDING: the trichotomy / parity structure is identical across")
        print(f"  N ∈ {N_VALUES}. The Klein routing of Q is N-invariant in this range;")
        print("  N=3 is representative, not special.")

    # ----- Routing rule (R2b) + second condition (R2c): letters only -----
    validate_routing(N_VALUES)
    validate_q7_reason(N_VALUES)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
