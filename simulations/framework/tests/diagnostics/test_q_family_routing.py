"""Tests for the two-term palindrome Klein-routing classifier.

The authoritative logic is the validated probe simulations/q6_klein_routing_two_term.py
(route() + q7_reason(), bit-exact 0 mismatches vs fw.classify_pauli_pair at N=3,4,5).
The key validation here: classify_two_term_palindrome's fate reproduces
fw.classify_pauli_pair for ALL 36 two-term bond-bilinear combos at N=3.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw

# The 9 traceless bond bilinears in canonical H_LABELS order.
H_LABELS = [a + b for a in 'XYZ' for b in 'XYZ']
COMBOS = list(combinations(H_LABELS, 2))


def test_all_36_fate_matches_classify_pauli_pair_N3():
    """fate reproduces fw.classify_pauli_pair for all 36 combos at N=3.

    This is THE validation: 0 mismatches against the bit-exact eigensolve route.
    """
    chain = fw.ChainSystem(3, topology='chain')
    mismatches = []
    for t1, t2 in COMBOS:
        result = fw.classify_two_term_palindrome(t1, t2, N=3)
        expected = fw.classify_pauli_pair(chain, [t1, t2])
        if result['fate'] != expected:
            mismatches.append((t1, t2, result['fate'], expected))
    assert mismatches == [], f"fate mismatches vs classify_pauli_pair: {mismatches}"


def test_family_truly_xx_yy():
    """XX+YY: both Mother, canonical Π pairs spectrum → truly / 'P1'."""
    r = fw.classify_two_term_palindrome('XX', 'YY')
    assert r['fate'] == 'truly'
    assert r['q_family'] == 'P1'


def test_family_soft_uniform_xx_xz():
    """XX+XZ: shared uniform Q-family {P4} → soft / 'uniform'."""
    r = fw.classify_two_term_palindrome('XX', 'XZ')
    assert r['fate'] == 'soft'
    assert r['q_family'] == 'uniform'


def test_family_soft_alternating_xy_yx():
    """XY+YX: site-parity (alternating) escape → soft / 'alternating'."""
    r = fw.classify_two_term_palindrome('XY', 'YX')
    assert r['fate'] == 'soft'
    assert r['q_family'] == 'alternating'


def test_family_soft_continuous_xz_yz():
    """XZ+YZ: same-site X&Y over shared dark Z → soft, closed by a local continuous map.

    Formerly labelled 'non_local'; the continuous per-site rotation M (M²=−I) is a genuine
    product Π, so the family is 'continuous' and the mirror is local. See
    test_crossover_product_pi for the bit-exact verification.
    """
    r = fw.classify_two_term_palindrome('XZ', 'YZ')
    assert r['fate'] == 'soft'
    assert r['q_family'] == 'continuous'


def test_family_hard_xx_xy():
    """XX+XY: no shared Q, no escape → hard / None."""
    r = fw.classify_two_term_palindrome('XX', 'XY')
    assert r['fate'] == 'hard'
    assert r['q_family'] is None


def test_parity_break_true_for_xx_xy():
    """XX+XY breaks the global X-parity (‖[H, X^⊗N]‖ > 0)."""
    r = fw.classify_two_term_palindrome('XX', 'XY', N=3)
    assert r['parity_break'] is True


def test_parity_break_false_for_xx_yy():
    """XX+YY commutes with the global X-flip → no parity break."""
    r = fw.classify_two_term_palindrome('XX', 'YY', N=3)
    assert r['parity_break'] is False


def test_input_accepts_string_and_tuple():
    """term arguments accepted as 2-char strings or 2-tuples, identically."""
    r_str = fw.classify_two_term_palindrome('XY', 'YX')
    r_tup = fw.classify_two_term_palindrome(('X', 'Y'), ('Y', 'X'))
    assert r_str['fate'] == r_tup['fate']
    assert r_str['q_family'] == r_tup['q_family']
    assert r_tup['term1'] == 'XY'
    assert r_tup['term2'] == 'YX'


def test_result_dict_keys():
    """The returned dict carries the full documented key set."""
    r = fw.classify_two_term_palindrome('YZ', 'ZY', N=3)
    assert set(r.keys()) == {
        'term1', 'term2', 'klein1', 'klein2',
        'parity_break', 'fate', 'q_family', 'reason',
    }
    assert r['klein1'] == fw.klein_index('YZ')
    assert r['klein2'] == fw.klein_index('ZY')


def test_hard_reason_is_structural_q7():
    """For a hard combo, reason is the q7 Klein-cell explanation, not the rule clause."""
    r = fw.classify_two_term_palindrome('XY', 'XZ')
    assert r['fate'] == 'hard'
    assert r['q_family'] is None
    # q7_reason for F_a+F_b names the unroutable both-X-and-Y demand.
    assert 'X' in r['reason'] and 'Y' in r['reason']
