#!/usr/bin/env python3
"""Does the two-term palindrome routing rotate with the dephasing axis?

The Klein-routing of the hidden symmetry Q (q6_klein_routing_two_term.py,
THE_OTHER_SIDE Q6/Q7) was worked out under Z-dephasing: the fate (truly / soft /
hard) of H = Σ_bonds (t1 + t2) is read from the two bilinears' letters, with the
"light" (damped) axis being {X, Y} and the dark axis {I, Z}.

Z-, X-, and Y-dephasing are related by single-qubit Cliffords (the Hadamard
X↔Z; the axis-swap Y↔Z), so the whole story should ROTATE: the X-dephasing fate
table should be the Z table with the bilinears relabelled X↔Z (the dark Mother
moving ZZ → XX, the lit channels rotating), and the Y table the Z table relabelled
Y↔Z. If it rotates cleanly the routing is not a Z-accident but SU(2)-covariant;
the framework's dissipator-resonance law already asserts this at the trichotomy
level (f77_trichotomy docstring), and this probe makes it explicit for the k=2
routing and pins the exact permutation, then checks that our Z-only rule
(fw.classify_two_term_palindrome) reproduces the X/Y fates once the combo is
relabelled.

Self-validating: at N=3 and N=4 the fate counts stay {3,19,14} under every
dephasing letter, the X/Y tables equal the relabelled Z table with zero
mismatches, and the relabelled Z-rule reproduces them.

Conventions: letter strings throughout; no raw indices.
"""
from __future__ import annotations

import sys
from itertools import combinations

sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


H_LABELS = [a + b for a in 'XYZ' for b in 'XYZ']
COMBOS = list(combinations(H_LABELS, 2))
assert len(COMBOS) == 36
N_VALUES = [3, 4, 5]
DEPHASE_LETTERS = ['Z', 'X', 'Y']

# The Clifford axis-swap that maps Z-dephasing onto each other dephasing letter:
# Hadamard X↔Z takes Z-deph → X-deph; the Y↔Z swap takes Z-deph → Y-deph.
RELABEL = {'Z': ('Z', 'Z'), 'X': ('X', 'Z'), 'Y': ('Y', 'Z')}


def letters(bil):
    return (bil[0], bil[1])


def combo_label(bil_a, bil_b):
    """Canonical 'A+B' label, the two bilinears ordered by H_LABELS index."""
    ia, ib = H_LABELS.index(bil_a), H_LABELS.index(bil_b)
    return f"{bil_a}+{bil_b}" if ia <= ib else f"{bil_b}+{bil_a}"


def swap_letter(ch, a, b):
    return b if ch == a else (a if ch == b else ch)


def relabel_bilinear(bil, a, b):
    return swap_letter(bil[0], a, b) + swap_letter(bil[1], a, b)


def relabel_combo(bil_a, bil_b, swap):
    """Relabel both bilinears under the letter swap, return the canonical label."""
    a, b = swap
    return combo_label(relabel_bilinear(bil_a, a, b), relabel_bilinear(bil_b, a, b))


def classify_under(N, dephase_letter):
    """Fate of all 36 combos under the given dephasing letter. {label: fate}."""
    chain = fw.ChainSystem(N, topology='chain')
    out = {}
    for bil_a, bil_b in COMBOS:
        fate = fw.classify_pauli_pair(
            chain, [letters(bil_a), letters(bil_b)], dephase_letter=dephase_letter)
        out[combo_label(bil_a, bil_b)] = fate
    return out


def counts(table):
    c = {'truly': 0, 'soft': 0, 'hard': 0}
    for f in table.values():
        c[f] += 1
    return c


def main():
    print("Does the two-term palindrome routing rotate with the dephasing axis?")
    print(f"  36 two-term combos, N ∈ {N_VALUES}, dephasing ∈ {DEPHASE_LETTERS}")
    print(f"  Hypothesis: X-deph table = Z table relabelled X↔Z; Y-deph = Z relabelled Y↔Z")
    print()

    for N in N_VALUES:
        print("=" * 74)
        print(f"N = {N}")
        print("=" * 74)
        tables = {L: classify_under(N, L) for L in DEPHASE_LETTERS}

        # 1. Fate counts invariant under the dephasing letter (SU(2) covariance).
        for L in DEPHASE_LETTERS:
            c = counts(tables[L])
            print(f"  {L}-dephasing counts: truly={c['truly']}, soft={c['soft']}, hard={c['hard']}")
            assert c == {'truly': 3, 'soft': 19, 'hard': 14}, \
                f"{L}-deph counts not {{3,19,14}} at N={N}: {c}"

        # 2. Each non-Z table equals the Z table relabelled by the axis swap.
        for L in ('X', 'Y'):
            swap = RELABEL[L]
            mism = []
            for bil_a, bil_b in COMBOS:
                label = combo_label(bil_a, bil_b)
                rotated = relabel_combo(bil_a, bil_b, swap)
                if tables[L][label] != tables['Z'][rotated]:
                    mism.append((label, rotated, tables[L][label], tables['Z'][rotated]))
            print(f"  {L}-deph vs Z relabelled {swap[0]}↔{swap[1]}: "
                  f"{len(mism)} mismatches / 36")
            for m in mism:
                print(f"      {m[0]}: {L}-fate {m[2]}  vs  Z[{m[1]}]={m[3]}")
            assert not mism, f"{L}-deph does not rotate cleanly from Z at N={N}"

        # 3. Our Z-only rule (the primitive), applied to the relabelled combo,
        #    reproduces the X/Y fate: the RULE itself is the rotated rule.
        for L in ('X', 'Y'):
            a, b = RELABEL[L]
            rule_mism = 0
            for bil_a, bil_b in COMBOS:
                label = combo_label(bil_a, bil_b)
                t1 = relabel_bilinear(bil_a, a, b)
                t2 = relabel_bilinear(bil_b, a, b)
                rule_fate = fw.classify_two_term_palindrome(t1, t2, N=N)['fate']
                if rule_fate != tables[L][label]:
                    rule_mism += 1
            print(f"  Z-rule relabelled {a}↔{b} reproduces {L}-deph fate: "
                  f"{36 - rule_mism}/36")
            assert rule_mism == 0, \
                f"Z-rule does not rotate to {L}-deph at N={N} ({rule_mism} off)"
        print()

    print("=" * 74)
    print("FINDING: the two-term routing is SU(2)-covariant. The X- and Y-dephasing")
    print("fate tables are the Z table relabelled by the Hadamard X↔Z and the swap")
    print("Y↔Z; the counts {3,19,14} are dephasing-invariant; and the Z-only rule")
    print("reproduces every X/Y fate once the combo is relabelled. The hidden-Q")
    print("routing rides on whichever axis the dephasing makes 'light', not on Z.")
    print("=" * 74)
    print("Done.")


if __name__ == "__main__":
    main()
