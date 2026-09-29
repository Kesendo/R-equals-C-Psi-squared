#!/usr/bin/env python3
"""F158 §(f5)'s odd-word pre-filter on the N = 3 census (docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md).

An invertible U with [H, U] = 0 and U A U^-1 = -A for every jump makes every word in {H, A_1 .. A_m}
with an odd number of A-letters traceless. `fw.odd_word_obstruction` looks for the shortest such word
with a nonzero trace, exactly (Pauli algebra over the Gaussian rationals); a nonzero trace rules the
palindrome out at every positive rate, and needs no Liouvillian and no rank.

The grid is the colouring page's N = 3 census (experiments/THE_PALINDROME_AS_A_COLOURING.md): graphs
P3, K3 and a bond with an isolated site, bonds ZZ, XX + YY and XX + YY + ZZ (weight 100), at most one
dephasing axis per site and at least one jump, every field pattern (fields 30, 22, 41 on sites 0, 1,
2), 36,288 rows. The palindrome verdict to compare with comes from F158's two kernels as ranks modulo
two primes (simulations/anticommuting_sum_census.py).

Checks (all must pass; prints "ALL CHECKS PASS"):
  (i)   the pre-filter never fires on a row the ranks call palindromic (the theorem says it cannot;
        a firing there would contradict F158 or expose a rank error);
  (ii)  the rows it certifies broken, pinned as measured: with words of one jump letter and total H
        power at most 4, and with up to three jump letters;
  (iii) the shortest firing word by power of H, pinned as measured;
  (iv)  a control that can fail: on the cascade bond of the colouring page with a field eps X on the
        first site the first firing word is H^2 A with trace 8 eps, and without that field none fires.

Run:  python simulations/f158_odd_word_prefilter.py
   >  simulations/results/f158_odd_word_prefilter.txt     (runtime about 1 minute on 22 cores)
"""
import collections
import itertools
import sys
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import anticommuting_sum_census as census   # noqa: E402  (its module body defines functions only)
import framework as fw                      # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def work(args):
    gname, edges, bset, deph = args
    jumps = [''.join(deph[l][0] if k == l else 'I' for k in range(3)) for l in range(3) if deph[l]]
    out = collections.Counter()
    for fields in itertools.product('IXYZ', repeat=3):
        terms = census.row_terms(3, edges, bset, fields)
        pal = census.nullity(terms, census.dark_of(deph)) == census.nullity(terms, census.lit_of(deph))
        h = {}
        for t, v in terms:
            h[t] = h.get(t, 0) + v
        o1 = fw.odd_word_obstruction(h, jumps, max_power=4, max_jumps=1)
        o3 = o1 or fw.odd_word_obstruction(h, jumps, max_power=4, max_jumps=3)
        key = 'palindromic' if pal else 'broken'
        out[key] += 1
        out[key + ', certified by one jump letter'] += o1 is not None
        out[key + ', certified by up to three'] += o3 is not None
        if o1:
            out['shortest word ' + o1['word'].split('·')[0]] += 1
    return out


if __name__ == "__main__":
    jobs = [(g, e, b, d) for g, e in census.GRAPHS[3].items() for b in census.BONDSETS
            for d in itertools.product([(), ('X',), ('Y',), ('Z',)], repeat=3) if any(d)]
    tot = collections.Counter()
    with Pool(22) as pool:
        for c in pool.imap_unordered(work, jobs, chunksize=2):
            tot.update(c)
    for k in sorted(tot):
        print(f"    {k}: {tot[k]}")
    check("(i) the pre-filter never fires on a row the ranks call palindromic",
          tot['palindromic'] == 9609 and tot['palindromic, certified by one jump letter'] == 0
          and tot['palindromic, certified by up to three'] == 0, f"{tot['palindromic']} palindromic rows")
    check("(ii) broken rows certified exactly: 24,671 of 26,679 by one jump letter, 25,831 by up to three",
          (tot['broken'], tot['broken, certified by one jump letter'], tot['broken, certified by up to three'])
          == (26679, 24671, 25831))
    check("(iii) shortest firing word: H^1 A on 17,091 rows, H^2 A on 6,292, H^3 A on 1,288",
          (tot['shortest word H^1'], tot['shortest word H^2'], tot['shortest word H^3']) == (17091, 6292, 1288))
    eps = Fraction(1, 3)
    h = {'XX': 1, 'YY': 1, 'YI': Fraction(3, 10), 'IZ': Fraction(11, 50), 'XI': eps}
    o = fw.odd_word_obstruction(h, ['IX'])
    none = fw.odd_word_obstruction({k: v for k, v in h.items() if k != 'XI'}, ['IX'])
    check("(iv) control: the cascade bond with a field eps X fires at H^2 A with trace 8 eps; without it, nothing",
          o is not None and o['word'].startswith('H^2') and o['trace'] == (8 * eps, 0) and none is None)
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL CHECKS PASS")
