#!/usr/bin/env python3
"""V-Effect 14/19/3 — combinatorial derivation attempt.

The original V-Effect at N=3 (V_EFFECT_PALINDROME, 2026-03-26): 36 possible
two-bilinear Pauli-pair Hamiltonians, all palindromic at N=2, but 14 of
36 break palindromic spectrum at N=3. Framework refinement:
  3 truly (operator-equation holds)
  19 soft (spectrum-paired, operator-broken)
  14 hard (spectrum-broken)

Where do these specific counts come from?

This script enumerates all 36 unordered cross-pairs (T1 ≠ T2) of bilinears
T_k ∈ {X, Y, Z} × {X, Y, Z} (= 9 options), excludes self-pairs giving
C(9,2) = 36, and for each pair computes:

  Framework verdict:  truly / soft / hard (palindrome_residual + spectrum)
  Structural features:
    - bit_a parity sum (#X+#Y in (T1, T2))
    - bit_b parity sum (#Y+#Z)
    - both-parity-even content: #{T_k ∈ {XX, YY, ZZ}}
    - bond-flipped: T2 = reverse(T1)?
    - same-letter cross: T1 has letter X and T2 has letter Y at same site?

Cross-tabulate verdict against feature combos. Look for: does some
specific feature predict the truly/soft/hard verdict exactly?
"""
import math
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def classify(H, N, gamma):
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    if float(np.linalg.norm(M)) < 1e-10:
        return 'truly'
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    sigma_g = N * gamma
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_g
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        used[i] = True
        if best_j != i:
            used[best_j] = True
        max_err = max(max_err, float(dists[best_j]))
    return 'soft' if max_err < 1e-6 else 'hard'


def bit_a_b_of(letter):
    return {'I': (0, 0), 'X': (1, 0), 'Y': (1, 1), 'Z': (0, 1)}[letter]


def features_of_pair(t1, t2):
    """Structural features of the bilinear pair (T1, T2)."""
    out = {}
    # bit_a, bit_b for each letter
    a1 = [bit_a_b_of(c)[0] for c in t1]
    b1 = [bit_a_b_of(c)[1] for c in t1]
    a2 = [bit_a_b_of(c)[0] for c in t2]
    b2 = [bit_a_b_of(c)[1] for c in t2]

    out['bit_a_sum'] = (sum(a1) + sum(a2)) % 2
    out['bit_b_sum'] = (sum(b1) + sum(b2)) % 2
    out['both_parity_even_count'] = (
        (1 if (sum(a1) % 2 == 0 and sum(b1) % 2 == 0) else 0)
        + (1 if (sum(a2) % 2 == 0 and sum(b2) % 2 == 0) else 0)
    )
    out['t1_both_even'] = (sum(a1) % 2 == 0 and sum(b1) % 2 == 0)
    out['t2_both_even'] = (sum(a2) % 2 == 0 and sum(b2) % 2 == 0)
    out['bond_flipped'] = (t1[::-1] == t2)
    # Same-position match: do T1 and T2 share letter at same site?
    out['same_at_site_0'] = (t1[0] == t2[0])
    out['same_at_site_1'] = (t1[1] == t2[1])
    # Are T1 and T2 in the same "Heisenberg row" {XX, YY, ZZ}?
    heisenberg = {('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}
    out['t1_in_heisenberg'] = (t1 in heisenberg)
    out['t2_in_heisenberg'] = (t2 in heisenberg)
    return out


def main():
    N = 3
    GAMMA, J = 0.1, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    # The 36 V-Effect pairs: cross-pairs (T1 ≠ T2) of T_k ∈ {X,Y,Z}^2.
    # = C(9, 2) = 36
    pauli_xyz = ['X', 'Y', 'Z']
    bilinears = [(a, b) for a in pauli_xyz for b in pauli_xyz]  # 9 of them
    pairs = []
    for t1, t2 in combinations(bilinears, 2):
        # Sort so the unordered pair is canonical
        pair = tuple(sorted([t1, t2]))
        pairs.append(pair)
    pairs = list(set(pairs))
    pairs.sort()
    assert len(pairs) == 36, f"Expected 36 pairs, got {len(pairs)}"

    print(f"V-Effect 36-enum combinatorial derivation at N={N}")
    print(f"  γ_dephasing = {GAMMA}")
    print(f"  Total pairs: {len(pairs)}")
    print()

    rows = []
    for t1, t2 in pairs:
        terms = [(t1[0], t1[1], J), (t2[0], t2[1], J)]
        H = fw._build_bilinear(N, bonds, terms)
        cat = classify(H, N, GAMMA)
        feats = features_of_pair(t1, t2)
        label = f"{t1[0]}{t1[1]}+{t2[0]}{t2[1]}"
        rows.append({'label': label, 'category': cat, **feats,
                     't1': t1, 't2': t2})

    counts = Counter(r['category'] for r in rows)
    print(f"  Category counts: {dict(counts)}")
    print(f"  (Original V-Effect at N=3: 14 hard / 22 not-hard;")
    print(f"   framework refinement: 14 hard / 19 soft / 3 truly)")
    print()

    # ===== Tabulate by both_parity_even_count =====
    print(f"=" * 80)
    print(f"By both-parity-even count (terms in {{XX, YY, ZZ}}):")
    print(f"=" * 80)
    print()
    print(f"  {'count':>5s}  {'truly':>5s}  {'soft':>5s}  {'hard':>5s}  {'total':>5s}  members")
    for bp_count in [0, 1, 2]:
        sub = [r for r in rows if r['both_parity_even_count'] == bp_count]
        n_tr = sum(1 for r in sub if r['category'] == 'truly')
        n_so = sum(1 for r in sub if r['category'] == 'soft')
        n_hd = sum(1 for r in sub if r['category'] == 'hard')
        members = ', '.join(r['label'] for r in sub[:6])
        if len(sub) > 6:
            members += f" ... ({len(sub) - 6} more)"
        print(f"  {bp_count:>5d}  {n_tr:>5d}  {n_so:>5d}  {n_hd:>5d}  "
              f"{len(sub):>5d}  {members}")
    print()

    # ===== Bond-flipped vs not =====
    print(f"=" * 80)
    print(f"By bond-flipped pair (T2 = reverse(T1)):")
    print(f"=" * 80)
    print()
    for bf in [True, False]:
        sub = [r for r in rows if r['bond_flipped'] == bf]
        n_tr = sum(1 for r in sub if r['category'] == 'truly')
        n_so = sum(1 for r in sub if r['category'] == 'soft')
        n_hd = sum(1 for r in sub if r['category'] == 'hard')
        flag = "yes" if bf else "no"
        print(f"  bond-flipped={flag}: truly={n_tr}, soft={n_so}, "
              f"hard={n_hd}, total={len(sub)}")
        if bf:
            print(f"    members: {[r['label'] for r in sub]}")
    print()

    # ===== Truly count derivation =====
    print(f"=" * 80)
    print(f"TRULY count derivation:")
    print(f"=" * 80)
    print()
    truly_rows = [r for r in rows if r['category'] == 'truly']
    print(f"  All truly cases ({len(truly_rows)}):")
    for r in truly_rows:
        print(f"    {r['label']}: t1={r['t1']}, t2={r['t2']}, "
              f"both_parity_even_count={r['both_parity_even_count']}")
    print()
    print(f"  Hypothesis: truly iff both T1 AND T2 are in {{XX, YY, ZZ}}")
    truly_predicted = [r for r in rows if r['t1_in_heisenberg']
                       and r['t2_in_heisenberg']]
    actual_truly_set = {r['label'] for r in truly_rows}
    predicted_truly_set = {r['label'] for r in truly_predicted}
    if actual_truly_set == predicted_truly_set:
        print(f"  ✓ EXACT MATCH: 3 cross-pairs of {{XX, YY, ZZ}} ↔ 3 truly")
    else:
        print(f"  ✗ Mismatch.")
        print(f"    Actual: {sorted(actual_truly_set)}")
        print(f"    Predicted: {sorted(predicted_truly_set)}")
    print()

    # ===== Soft vs Hard among non-truly =====
    print(f"=" * 80)
    print(f"SOFT vs HARD among 33 non-truly:")
    print(f"=" * 80)
    print()

    non_truly = [r for r in rows if r['category'] != 'truly']
    print(f"  By bit_a_sum:")
    for ba in [0, 1]:
        sub = [r for r in non_truly if r['bit_a_sum'] == ba]
        n_so = sum(1 for r in sub if r['category'] == 'soft')
        n_hd = sum(1 for r in sub if r['category'] == 'hard')
        print(f"    bit_a_sum={ba}: soft={n_so}, hard={n_hd}, total={len(sub)}")
    print()

    print(f"  By bit_b_sum:")
    for bb in [0, 1]:
        sub = [r for r in non_truly if r['bit_b_sum'] == bb]
        n_so = sum(1 for r in sub if r['category'] == 'soft')
        n_hd = sum(1 for r in sub if r['category'] == 'hard')
        print(f"    bit_b_sum={bb}: soft={n_so}, hard={n_hd}, total={len(sub)}")
    print()

    print(f"  By both-parity-even count (in non-truly):")
    for bp in [0, 1]:
        sub = [r for r in non_truly if r['both_parity_even_count'] == bp]
        n_so = sum(1 for r in sub if r['category'] == 'soft')
        n_hd = sum(1 for r in sub if r['category'] == 'hard')
        print(f"    both_parity_even={bp}: soft={n_so}, hard={n_hd}, "
              f"total={len(sub)}")
    print()

    # ===== Joint feature ===== (bit_a_sum, bit_b_sum)
    print(f"  By (bit_a_sum, bit_b_sum):")
    print(f"  {'(a,b)':>6s}  {'soft':>5s}  {'hard':>5s}  {'total':>5s}")
    for ba in [0, 1]:
        for bb in [0, 1]:
            sub = [r for r in non_truly
                   if r['bit_a_sum'] == ba and r['bit_b_sum'] == bb]
            n_so = sum(1 for r in sub if r['category'] == 'soft')
            n_hd = sum(1 for r in sub if r['category'] == 'hard')
            print(f"  ({ba},{bb})  {n_so:>5d}  {n_hd:>5d}  {len(sub):>5d}")


if __name__ == "__main__":
    main()
