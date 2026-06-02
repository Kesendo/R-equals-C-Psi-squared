"""Klein-routing classifier for two-term bond-bilinear palindrome fate (THE_OTHER_SIDE Q6/Q7).

A Z-dephased XY-type chain H = Σ_bonds (t1 + t2) meets one of three fates:
  truly : the canonical palindrome operator Π (P1) already pairs the spectrum (M = 0).
  soft  : Π fails (M ≠ 0) but the palindrome still holds; a hidden symmetry Q ≠ Π
          closes the spectrum-pairing. The benign ones.
  hard  : the palindrome is BROKEN; no Q closes it.

The fate is decided from the two bilinears' letters alone, no eigensolve. Z-dephasing
splits each site into a DC bus {I, Z} (immune) and an AC bus {X, Y} (damped). The
palindrome operator Q routes the AC bus, and per-site it comes in two crossover families:
P1 routes the X-channel, P4 routes the Y-channel. A single uniform product Q must serve
BOTH bond terms at once, so the chain is palindromic iff its two terms share a valid
uniform Q-family; failing that, it escapes via an alternating (site-parity) or a non-local
(entangled) Q; failing THAT, the palindrome is broken (hard).

This routing reproduces fw.classify_pauli_pair bit-exactly across N = 3, 4, 5
(0 mismatches over all 36 combos), validated in the probe
simulations/q6_klein_routing_two_term.py. The routing is N-independent; the optional
N here is used only for the parity-break cross-check ‖[H, X^⊗N]‖ > 0.

Public API:
  classify_two_term_palindrome(term1, term2, N=3, gamma_0=0.05)
"""
from __future__ import annotations

import numpy as np

from ..pauli import LABEL_TO_INDEX, site_op
from ..symmetry import klein_index, KLEIN_LABELS

_PARITY_TOL = 1e-9

# ---------------------------------------------------------------------------
# Routing rule (R2b): valid uniform-Q-family per term
# ---------------------------------------------------------------------------
# fam(t) = the set of uniform Q-families that close term t alone:
#   fam(XX) = fam(YY) = {P1, P4}   both crossovers route a same-letter bilinear
#   fam(ZZ) = {all}                ZZ is dephasing-aligned, every Q fixes it
#   fam(YZ) = fam(ZY) = {P1}       only the X-router survives the Y·Z mix
#   fam(XZ) = fam(ZX) = {P4}       only the Y-router survives the X·Z mix
#   fam(XY) = fam(YX) = {}         no single uniform Q routes both X and Y
_ALL_Q = frozenset({'P1', 'P4', 'M2'})
_FAM = {
    'XX': frozenset({'P1', 'P4'}),
    'YY': frozenset({'P1', 'P4'}),
    'ZZ': _ALL_Q,
    'YZ': frozenset({'P1'}),
    'ZY': frozenset({'P1'}),
    'XZ': frozenset({'P4'}),
    'ZX': frozenset({'P4'}),
    'XY': frozenset(),
    'YX': frozenset(),
}

# The two non-local escapes: same-site X&Y collision over a shared dark Z port,
# rescued by an entangled (non-product) Q.
_NON_LOCAL_PAIRS = {frozenset({'XZ', 'YZ'}), frozenset({'ZX', 'ZY'})}

_MOTHER = {'XX', 'YY', 'ZZ'}


def _route(bil_a, bil_b):
    """Predict the fate of H = Σ_bonds (bil_a + bil_b) from the letters alone.

    Returns {fate, family, reason} where fate ∈ {'truly','soft','hard'} and
    family ∈ {'P1','uniform','alternating','non_local', None}.
    """
    fa, fb = _FAM[bil_a], _FAM[bil_b]

    # 1. Shared uniform Q-family: the spectrum pairs under one product Q.
    if fa & fb:
        if bil_a in _MOTHER and bil_b in _MOTHER:
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
    if frozenset({bil_a, bil_b}) in _NON_LOCAL_PAIRS:
        return {'fate': 'soft', 'family': 'non_local',
                'reason': 'same-site X&Y collision over shared dark Z; entangled Q'}

    # 4. No Q closes it: the palindrome is broken.
    return {'fate': 'hard', 'family': None,
            'reason': 'no shared uniform Q, no alternating/non-local escape'}


# ---------------------------------------------------------------------------
# Second condition (R2c): WHY a hard combo is hard (structural Klein-cell reason)
# ---------------------------------------------------------------------------
# Beyond breaking X-parity, a combo is hard iff it carries an irreducible,
# unroutable same-qubit X/Y demand. Read through the Klein cells of the two
# terms, the discriminator is a cell-pair rule with one refinement per split cell.
#   M   = (0,0) = {XX, YY, ZZ}   Mother
#   F_a = (0,1) = {XY, YX}       Y-Father  (Π²-odd)
#   C   = (1,0) = {YZ, ZY}       Child     (dephase-aligned)
#   F_b = (1,1) = {XZ, ZX}       Z-Father  (Π²-odd)


def _cell_pair(bil_a, bil_b):
    """The unordered Klein-cell pair of the two terms, e.g. ('F_a', 'M')."""
    return tuple(sorted([KLEIN_LABELS[klein_index(bil_a)],
                         KLEIN_LABELS[klein_index(bil_b)]]))


def _site_letters(bil_a, bil_b):
    """Letters landing on (site0, site1), one from each bond term."""
    return (bil_a[0], bil_b[0]), (bil_a[1], bil_b[1])


def _same_site_xy(bil_a, bil_b):
    """True iff one single qubit is asked to carry BOTH an X and a Y."""
    s0, s1 = _site_letters(bil_a, bil_b)
    return ({'X', 'Y'} <= set(s0)) or ({'X', 'Y'} <= set(s1))


def _q7_reason(bil_a, bil_b):
    """The structural reason a combo is (or is not) hard, via Klein cells.

    Returns {hard, cell_pair, reason}. The hard-by-reason set equals the 14 hard
    combos bit-exactly at N = 3, 4, 5, with the lit-vs-dark ({F_a, M}) and
    same-site/different-site ({C, F_b}) discriminations resolved.
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
    if k == ('F_a', 'M'):
        if _same_site_xy(bil_a, bil_b):
            return {'hard': True, 'cell_pair': k,
                    'reason': 'F_a+M lit Mother: X and Y forced onto one qubit'}
        return {'hard': False, 'cell_pair': k,
                'reason': 'F_a+M dark Mother (ZZ): dephasing-aligned, no X/Y conflict'}

    # Split cell B: Z-Father (XZ/ZX) + Child (YZ/ZY).
    if k == ('C', 'F_b'):
        if frozenset({bil_a, bil_b}) in _NON_LOCAL_PAIRS:
            return {'hard': False, 'cell_pair': k,
                    'reason': 'C+F_b same-site X&Y over shared dark Z: entangled (non-local) rescue'}
        return {'hard': True, 'cell_pair': k,
                'reason': 'C+F_b lit X,Y on different sites: no Q closes it'}

    # All other cells (M+M, M+C, M+F_b, C+C, F_a+F_a, ...): never hard.
    return {'hard': False, 'cell_pair': k,
            'reason': 'no irreducible same-qubit X/Y conflict'}


# ---------------------------------------------------------------------------
# Parity-break cross-check (independent of the routing rule)
# ---------------------------------------------------------------------------

def _parity_break(N, bil_a, bil_b):
    """True iff H = Σ_bonds (bil_a + bil_b) fails to commute with X^⊗N.

    H is built directly on the open chain bonds (k, k+1) via site_op; X^⊗N is the
    global spin-flip. ‖[H, X^⊗N]‖ > tol signals the X-parity Z₂ symmetry is broken.
    """
    d = 2 ** N
    bonds = [(k, k + 1) for k in range(N - 1)]
    H = np.zeros((d, d), dtype=complex)
    for (a_letter, b_letter) in ((bil_a[0], bil_a[1]), (bil_b[0], bil_b[1])):
        for (i, j) in bonds:
            H = H + site_op(N, i, a_letter) @ site_op(N, j, b_letter)

    Xn = site_op(N, 0, 'X')
    for site in range(1, N):
        Xn = Xn @ site_op(N, site, 'X')

    comm = H @ Xn - Xn @ H
    return float(np.linalg.norm(comm)) > _PARITY_TOL


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _normalize_term(term):
    """Coerce a 2-char string ('XY') or 2-tuple (('X','Y')) into a canonical 'XY' string.

    Validates that each letter is a recognized Pauli ('I','X','Y','Z') and that the
    term has exactly two letters.
    """
    if isinstance(term, str):
        letters = term
    else:
        letters = ''.join(term)
    if len(letters) != 2:
        raise ValueError(f"term must be two letters; got {term!r}")
    for ch in letters:
        if ch not in LABEL_TO_INDEX:
            raise ValueError(f"unknown Pauli letter {ch!r} in term {term!r}")
    return letters


def classify_two_term_palindrome(term1, term2, N=3, gamma_0=0.05):
    """Classify the palindrome fate of H = Σ_bonds (term1 + term2) under Z-dephasing.

    Decides the fate (truly / soft / hard) and the Klein-routing of the hidden
    symmetry Q from the two bilinears' letters alone, no eigensolve. Reproduces
    fw.classify_pauli_pair bit-exactly across N = 3, 4, 5 (0 mismatches over all
    36 two-term combos). The routing is N-independent; N here drives only the
    parity-break cross-check.

    Args:
        term1, term2: each a 2-char Pauli string ('XY') or 2-tuple (('X','Y')).
        N: chain length, used only for the parity-break test (default 3).
        gamma_0: per-site dephasing rate. Carried for signature parity with the
            cockpit; the routing fate does not depend on it (default 0.05).

    Returns:
        dict with keys:
          term1, term2  : the canonical 'XY'-string forms.
          klein1, klein2: each term's Klein index (bit_a, bit_b) from klein_index.
          parity_break  : bool, True iff ‖[H, X^⊗N]‖ > 1e-9 (global X-parity broken).
          fate          : 'truly' | 'soft' | 'hard'.
          q_family      : 'P1' (truly) | 'uniform' | 'alternating' | 'non_local'
                          (soft) | None (hard).
          reason        : structural explanation. For hard combos this is the
                          q7 Klein-cell reading; for soft/truly it is the rule clause.

    Raises:
        ValueError: if a term is not exactly two recognized Pauli letters.
    """
    t1 = _normalize_term(term1)
    t2 = _normalize_term(term2)

    r = _route(t1, t2)
    fate = r['fate']
    q_family = r['family']

    if fate == 'hard':
        # Prefer the structural Klein-cell reason for hard combos.
        reason = _q7_reason(t1, t2)['reason']
    else:
        reason = r['reason']

    return {
        'term1': t1,
        'term2': t2,
        'klein1': klein_index(t1),
        'klein2': klein_index(t2),
        'parity_break': _parity_break(N, t1, t2),
        'fate': fate,
        'q_family': q_family,
        'reason': reason,
    }
