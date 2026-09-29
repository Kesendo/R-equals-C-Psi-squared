"""F158 §(f5) odd-word pre-filter: can this generator be palindromic at all?

The question it answers: given H and Pauli-string jumps, is there a cheap certificate that NO
palindrome about -sigma exists, at any positive rates, before any Liouvillian is built?

F158 (docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md §(f5)): an invertible U with [H, U] = 0 and
U A_l U^-1 = -A_l makes every word in {H, A_1 .. A_m} with an odd number of A-letters traceless,
since conjugation by U fixes H and flips the sign of such a word. A nonzero trace of one such word
therefore rules the palindrome out. A finite budget of words that are all traceless proves nothing;
over every word the condition would also be sufficient (traces of all words determine a finite-
dimensional *-representation up to unitary equivalence, so (H, A) and (H, -A) would be unitarily
equivalent), but this function only ever looks at a finite budget.

Everything here is exact: H is a dict {Pauli string: rational coefficient}, the jumps are Pauli
strings, products are taken in the Pauli algebra with Gaussian-rational coefficients, and a trace is
2^N times the coefficient of the identity string. No matrices, no floats, no tolerance.

For a qubit chain with Z dephasing on site l the shortest word already decides a detuning:
Tr(H Z_l) = 2^N delta_l, with delta_l the coefficient of the Pauli string Z_l in H (a spin-convention
field h/2 reads 2^(N-1) h), so any longitudinal field on a dephased site rules the full palindrome out
(a uniform one included; the palindrome of the decay RATES can survive, which this does not see).
"""
from __future__ import annotations

import itertools
from fractions import Fraction

from ..pauli import pauli_product

__all__ = ["odd_word_obstruction"]

def _times_i_power(c, e):
    """multiply a Gaussian rational c = (re, im) by i**e"""
    re, im = c
    for _ in range(e % 4):
        re, im = -im, re
    return re, im


def _mul(A, B):
    """product of two operators given as {string: (re, im)}"""
    out = {}
    for s, (ar, ai) in A.items():
        for t, (br, bi) in B.items():
            e, u = pauli_product(s, t)
            c = _times_i_power((ar * br - ai * bi, ar * bi + ai * br), e)
            r0, i0 = out.get(u, (Fraction(0), Fraction(0)))
            out[u] = (r0 + c[0], i0 + c[1])
    return {u: c for u, c in out.items() if c != (0, 0)}


def _trace(A, n):
    re, im = A.get('I' * n, (Fraction(0), Fraction(0)))
    return (re * 2 ** n, im * 2 ** n)


def odd_word_obstruction(h_terms, jumps, max_power=4, max_jumps=3):
    """The shortest word with an odd number of jump letters whose trace is nonzero, or None.

    Args:
        h_terms: {Pauli string: rational coefficient} (int, Fraction, or anything Fraction accepts).
        jumps: list of Pauli strings, the active jumps (each a Hermitian involution).
        max_power: largest total power of H in a word.
        max_jumps: largest (odd) number of jump letters, 1 or 3.

    Returns:
        None when every word in the budget is traceless (no conclusion), or a dict with the word
        (as a readable string), its trace (re, im) as Fractions, and 'palindrome': False, meaning
        F158 rules the palindrome out for these jumps at every positive rates.
    """
    if not jumps:
        raise ValueError("F158 needs at least one jump")
    n = len(jumps[0])
    H = {s: (Fraction(c), Fraction(0)) for s, c in h_terms.items() if Fraction(c) != 0}
    one = {'I' * n: (Fraction(1), Fraction(0))}
    powers = [one]
    for _ in range(max_power):
        powers.append(_mul(powers[-1], H))
    J = [{a: (Fraction(1), Fraction(0))} for a in jumps]
    words = []
    for k in range(max_power + 1):                        # one jump letter: Tr(H^k A), cyclic
        for i in range(len(jumps)):
            words.append((k, f"H^{k}·A{i}" if k else f"A{i}", [k], [i]))
    if max_jumps >= 3:
        for total in range(max_power + 1):
            for a, b in itertools.product(range(total + 1), repeat=2):
                c = total - a - b
                if c < 0:
                    continue
                for trip in itertools.product(range(len(jumps)), repeat=3):
                    words.append((total, f"H^{a}·A{trip[0]}·H^{b}·A{trip[1]}·H^{c}·A{trip[2]}", [a, b, c], list(trip)))
    words.sort(key=lambda w: (w[0], len(w[3])))
    for _, name, hp, js in words:
        W = one
        for p, j in zip(hp, js):
            W = _mul(_mul(W, powers[p]), J[j])
        tr = _trace(W, n)
        if tr != (0, 0):
            return {'word': name, 'trace': tr, 'palindrome': False}
    return None
