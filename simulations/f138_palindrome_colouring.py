#!/usr/bin/env python3
"""F138's clauses as a colouring of the graph by lit letters (experiments/THE_PALINDROME_AS_A_COLOURING.md).

F158's far kernel is W = {W : [H, W] = 0, A_l W A_l = -W for every jump}; a Pauli string F in W
carries the palindrome (F158's sufficiency with an invertible element). For H a sum of distinct
Pauli strings with nonzero coefficients and single-letter Pauli jumps, F in W is a COLOURING: at
each site one letter (its colour) that every jump there lights, i.e. that anticommutes with every
dephasing letter of the site (any of I, X, Y, Z on an undephased site); a bond P(x)P is satisfied
when both or neither of its sites anticommute with P; a field P at a site is satisfied when the
colour there is I or P. Distinct terms T give distinct products T F, so [H, F] = 0 exactly when
every term commutes with F. (Not F103 section 7.12's 2-colouring, which picks per-site letters
that ANTIcommute with H.)

Stages (all must pass; prints "ALL STAGES PASS"):
  A  EXACT. The colouring rule, read through fw.cube_step (None exactly when two strings commute),
     agrees with the matrix commutator for every two-site term P(x)P and every one-site field against
     every letter pair, and a sum of distinct strings commutes with F exactly when every term does
     (checked on the Heisenberg bond against all 16 two-site F).
  B  F138's own domain, at most one dephasing axis per site, N = 3 on P3, K3 and a bond plus an
     isolated site, every nonempty bond letter set of {XX, YY, ZZ}, all 4^3 dephasing x 4^3 field
     patterns, fields at the committed magnitudes (0.30, 0.22, 0.41), positive signs. Ground truth:
     F158's count dim ker L = dim ker(L + 2 sigma), read by singular values (at least one jump; with
     none the spectrum of -i ad_H pairs about 0). Gated: (i) a colouring never exists without a
     palindrome; (ii) the rows where the palindrome holds and F138's clauses fail reproduce the
     counts F138's registry entry records (22 / 104 / 0 at a two-letter bond on P3 / bond+iso / K3,
     776 / 520 / 732 at one letter), an anchor this script did not produce; (iii) the number of those
     a colouring explains, pinned as measured (714 / 468 / 714 at one letter, 0 at two); (iv) the
     singular values separate by at least six decades.
  C  The rows F138's clauses accept where the palindrome fails, over dephasing patterns with up to
     three axes per site, the clauses read fairly: the two-axis ceiling is lifted only in a component
     that has bonds and whose bonds carry one letter, and three axes on one site fail clause 1 (F138:
     "The depolarizing channel is clause 1 failing"; MIRROR_SYMMETRY_PROOF: three axes on a single
     site break the Ising bond too). Tallied: every such row carries a two-axis site in a bonded
     component of one bond letter and three axes (a consistency check: at one axis per site the
     clauses imply a colouring); none has a colouring (F158's corollary, a check of the float verdicts); the
     singular values separate by at least six decades.
  D  F138's sentence "a three-axis component must then be field-free": on P3 with the ZZ bond and one
     axis per site (X, Z, Y), a field along the colour of its site pairs (an X field on the Z-dephased
     middle site), while a uniform X field and any Z field do not; the colouring predicts all three.
  E  EXACT. Beyond the colouring, one row worked out: P3 with H = a (XX + YY) on bond (0,1) and
     b (XX + YY) on bond (1,2), jumps X on site 0, Z on site 1, Y on site 2. U = a YYZ + b ZXX commutes
     with H, anticommutes with every jump and squares to (a^2 + b^2) I, entry by entry in integer
     arithmetic at (a, b) = (1, 1), (2, 3), (3, -1), (0, 1), (1, 0); for ab != 0 no colouring exists,
     at a = 0 or b = 0 U is itself a string; at a = b = 1 the far kernel is one-dimensional (the
     singular-value count), so U spans it. The equal-weight U is F138's; the weighted form was
     handed over by a second session (Codex, 2026-09-29) and is checked here.

Run:  python simulations/f138_palindrome_colouring.py
   >  simulations/results/f138_palindrome_colouring.txt     (runtime about 20 minutes)
"""
import collections
import itertools
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


LET = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], complex),
       'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.array([[1, 0], [0, -1]], complex)}
N = 3
d = 2 ** N
GRAPHS = {'P3': [(0, 1), (1, 2)], 'K3': [(0, 1), (1, 2), (0, 2)], 'bond+iso': [(0, 1)]}
BONDSETS = [s for r in (1, 2, 3) for s in itertools.combinations('XYZ', r)]
MAG = [0.30, 0.22, 0.41]
Id = np.eye(d)


def op(letters):
    M = np.array([[1]], complex)
    for c in letters:
        M = np.kron(M, LET[c])
    return M


def placed(letters_by_site):
    return ''.join(letters_by_site.get(k, 'I') for k in range(N))


def terms_of(edges, bset, fields):
    """the placed Pauli strings of H with their coefficients"""
    out = [(placed({a: P, b: P}), 1.0) for (a, b) in edges for P in bset]
    out += [(placed({l: P}), MAG[l]) for l, P in enumerate(fields) if P != 'I']
    return out


def H_of(terms):
    return sum(c * op(t) for t, c in terms) if terms else np.zeros((d, d), complex)


def colouring(terms, deph):
    """a string F, one colour per site, in W: every term commutes with F (fw.cube_step is None) and
    every colour anticommutes with the dephasing letters of its site"""
    allowed = [[P for P in 'XYZ' if P not in deph[l]] if deph[l] else list('IXYZ') for l in range(N)]
    for F in itertools.product(*allowed):
        F = ''.join(F)
        if all(fw.cube_step(t, F) is None for t, _ in terms):
            return F
    return None


def nullities(H, jumps):
    adH = np.kron(H, Id) - np.kron(Id, H.T)
    out = []
    for sgn in (-1, +1):
        M = np.vstack([adH] + [np.kron(A, A.conj()) + sgn * np.eye(d * d) for A in jumps])
        s = np.linalg.svd(M, compute_uv=False)
        out.append(s)
    return out


def components(edges):
    parent = list(range(N))

    def f(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in edges:
        parent[f(a)] = f(b)
    comp = collections.defaultdict(list)
    for l in range(N):
        comp[f(l)].append(l)
    return list(comp.values())


def f138_clauses(edges, bset, deph, fields):
    """F138 as its registry entry states it, per component carrying dephasing: (1) at most two
    distinct axes, a ceiling lifted only in a component that has bonds, all carrying one letter, and
    never for three axes on one site (the depolarizing channel is clause 1 failing); (2) the field
    has one common axis in the component, orthogonal to every dephasing axis present."""
    for comp in components(edges):
        axes = set().union(*[set(deph[l]) for l in comp])
        if not axes:
            continue
        has_bond = any(a in comp and b in comp for a, b in edges)
        lifted = has_bond and len(bset) == 1
        if any(len(deph[l]) == 3 for l in comp):
            return False
        if len(axes) > 2 and not lifted:
            return False
        fl = {fields[l] for l in comp if fields[l] != 'I'}
        if len(fl) > 1 or (fl & axes) or (len(axes) == 3 and fl):
            return False
    return True


def stage_a():
    print("## Stage A: the colouring rule against the matrix commutator")
    ok = True
    for P in 'XYZ':
        for a, b in itertools.product('IXYZ', repeat=2):
            T, S = P + P, a + b
            comm = np.kron(LET[P], LET[P]) @ np.kron(LET[a], LET[b]) - np.kron(LET[a], LET[b]) @ np.kron(LET[P], LET[P])
            rule = (sum(x != 'I' and y != 'I' and x != y for x, y in zip(T, S)) % 2 == 0)
            ok &= (fw.cube_step(T, S) is None) == (not comm.any()) == rule
        for a in 'IXYZ':
            comm = LET[P] @ LET[a] - LET[a] @ LET[P]
            ok &= (fw.cube_step(P, a) is None) == (not comm.any()) == (a in ('I', P))
    check("bond P(x)P and field P against every letter pair: fw.cube_step, the matrix commutator and "
          "the colouring rule agree (exact)", ok)
    heis = [('XX', 1), ('YY', 1), ('ZZ', 1)]
    Hm = sum(np.kron(LET[t[0]], LET[t[1]]) for t, _ in heis)
    ok = all((not (Hm @ np.kron(LET[F[0]], LET[F[1]]) - np.kron(LET[F[0]], LET[F[1]]) @ Hm).any())
             == all(fw.cube_step(t, F) is None for t, _ in heis)
             for F in (a + b for a, b in itertools.product('IXYZ', repeat=2)))
    check("the Heisenberg bond commutes with a two-site string exactly when each of its terms does "
          "(all 16 strings)", ok)


def stage_b():
    print()
    print("## Stage B: F138's domain, one axis per site; ground truth F158's count")
    t0 = time.time()
    single = [(), ('X',), ('Y',), ('Z',)]
    tab = collections.defaultdict(collections.Counter)
    kept_min, dropped_max = np.inf, 0.0
    string_without_palindrome = 0
    for g, edges in GRAPHS.items():
        for bset in BONDSETS:
            for deph in itertools.product(single, repeat=N):
                for fields in itertools.product('IXYZ', repeat=N):
                    terms = terms_of(edges, bset, fields)
                    F = colouring(terms, deph)
                    if not any(deph):
                        pal = True
                    else:
                        jumps = [op(placed({l: deph[l][0]})) for l in range(N) if deph[l]]
                        sN, sW = nullities(H_of(terms), jumps)
                        for s in (sN, sW):
                            kept_min = min(kept_min, s[s >= 1e-9].min())
                            if np.any(s < 1e-9):
                                dropped_max = max(dropped_max, s[s < 1e-9].max())
                        pal = int(np.sum(sN < 1e-9)) == int(np.sum(sW < 1e-9))
                    if F is not None and not pal:
                        string_without_palindrome += 1
                    if pal and not f138_clauses(edges, bset, deph, fields):
                        key = (g, len(bset))
                        tab[key]['exceptions'] += 1
                        tab[key]['coloured'] += F is not None
    for g in GRAPHS:
        for nl in (1, 2, 3):
            e = tab[(g, nl)]
            print(f"    {g:9s} {nl} bond letter(s): {e['exceptions']:4d} exceptions summed over the "
                  f"{sum(1 for b in BONDSETS if len(b) == nl)} bond set(s), {e['coloured']:4d} of them coloured")
    check("(i) a colouring never exists without the palindrome", string_without_palindrome == 0,
          f"{string_without_palindrome} rows")
    # F138's registry counts are per bond set; the three one-letter sets and the three two-letter
    # sets are related by relabelling letters and give the same count each
    anchor = {('P3', 1): 776, ('K3', 1): 732, ('bond+iso', 1): 520,
              ('P3', 2): 22, ('K3', 2): 0, ('bond+iso', 2): 104, ('P3', 3): 0, ('K3', 3): 0, ('bond+iso', 3): 0}
    per_set = {k: tab[k]['exceptions'] // (3 if k[1] < 3 else 1) for k in anchor}
    exact_multiple = all(tab[k]['exceptions'] == per_set[k] * (3 if k[1] < 3 else 1) for k in anchor)
    check("(ii) the exceptions per bond set reproduce F138's registry counts (776 / 732 / 520 at one "
          "letter, 22 / 0 / 104 at two, 0 at three)", exact_multiple and per_set == anchor,
          str(per_set))
    col = {k: tab[k]['coloured'] // (3 if k[1] < 3 else 1) for k in anchor}
    pinned = {('P3', 1): 714, ('K3', 1): 714, ('bond+iso', 1): 468,
              ('P3', 2): 0, ('K3', 2): 0, ('bond+iso', 2): 0, ('P3', 3): 0, ('K3', 3): 0, ('bond+iso', 3): 0}
    check("(iii) exceptions a colouring explains, per bond set, as measured on the first run: 714 / 714 / "
          "468 at one letter, 0 at two (a pin, so that a change is seen)", col == pinned, str(col))
    dec = np.log10(kept_min / max(dropped_max, 1e-300))
    check("(iv) singular values separate by at least six decades", dec >= 6,
          f"smallest kept {kept_min:.2e}, largest dropped {dropped_max:.2e}, {dec:.1f} decades  ({time.time() - t0:.0f} s)")


def stage_c():
    print()
    print("## Stage C: the rows F138's clauses accept where the palindrome fails")
    t0 = time.time()
    choices = [(), ('X',), ('Y',), ('Z',), ('X', 'Z'), ('X', 'Y'), ('Y', 'Z'), ('X', 'Y', 'Z')]
    bad = two_axis_one_letter = coloured = 0
    kept_min, dropped_max = np.inf, 0.0
    for g, edges in GRAPHS.items():
        for bset in BONDSETS:
            for deph in itertools.product(choices, repeat=N):
                if not any(deph):
                    continue
                for fields in itertools.product('IXYZ', repeat=N):
                    if not f138_clauses(edges, bset, deph, fields):
                        continue
                    terms = terms_of(edges, bset, fields)
                    jumps = [op(placed({l: P})) for l in range(N) for P in deph[l]]
                    sN, sW = nullities(H_of(terms), jumps)
                    for sv in (sN, sW):
                        kept_min = min(kept_min, sv[sv >= 1e-9].min())
                        if np.any(sv < 1e-9):
                            dropped_max = max(dropped_max, sv[sv < 1e-9].max())
                    if int(np.sum(sN < 1e-9)) != int(np.sum(sW < 1e-9)):
                        bad += 1
                        two_axis_one_letter += len(bset) == 1 and any(
                            any(a in comp and b in comp for a, b in edges)
                            and any(len(deph[l]) == 2 for l in comp)
                            and len(set().union(*[set(deph[l]) for l in comp])) == 3
                            for comp in components(edges))
                        coloured += colouring(terms, deph) is not None
    check("every row the fairly read clauses accept and the palindrome refuses has a two-axis site in a bonded "
          "component of one bond letter and three axes (consistency), and none has a colouring (F158's corollary)",
          bad > 0 and two_axis_one_letter == bad and coloured == 0,
          f"{bad} rows over the three graphs and seven bond sets, {two_axis_one_letter} of that kind, "
          f"{coloured} coloured  ({time.time() - t0:.0f} s)")
    dec = np.log10(kept_min / max(dropped_max, 1e-300))
    check("stage C singular values separate by at least six decades", dec >= 6,
          f"smallest kept {kept_min:.2e}, largest dropped {dropped_max:.2e}, {dec:.1f} decades")


def stage_d():
    print()
    print("## Stage D: 'a three-axis component must then be field-free', read on P3 with ZZ, axes X / Z / Y")
    edges, bset, deph = GRAPHS['P3'], ('Z',), (('X',), ('Z',), ('Y',))
    jumps = [op(placed({l: deph[l][0]})) for l in range(N)]
    rows = {'X field on the middle site': ('I', 'X', 'I'), 'uniform X field': ('X', 'X', 'X'),
            'Z field on the first site': ('Z', 'I', 'I'), 'no field': ('I', 'I', 'I')}
    ok = True
    for name, fields in rows.items():
        terms = terms_of(edges, bset, fields)
        sN, sW = nullities(H_of(terms), jumps)
        pal = int(np.sum(sN < 1e-9)) == int(np.sum(sW < 1e-9))
        F = colouring(terms, deph)
        ok &= pal == (F is not None)
        print(f"    {name:28s} palindrome {pal!s:5s} colouring {F}")
    check("the colouring decides these rows, and a field along the colour of its site pairs in a "
          "three-axis component", ok)


def stage_e():
    print()
    print("## Stage E: beyond the colouring, U = a YYZ + b ZXX on P3 with XX+YY bonds and jumps X, Z, Y")
    deph = (('X',), ('Z',), ('Y',))
    jumps = [op('XII'), op('IZI'), op('IIY')]
    ok = True
    for a, b in ((1, 1), (2, 3), (3, -1), (0, 1), (1, 0)):
        H = a * (op('XXI') + op('YYI')) + b * (op('IXX') + op('IYY'))
        U = a * op('YYZ') + b * op('ZXX')
        ok &= not (H @ U - U @ H).any()
        ok &= np.array_equal(U @ U, (a * a + b * b) * np.eye(d))
        ok &= all(not (A @ U + U @ A).any() for A in jumps)
        terms = ([('XXI', a), ('YYI', a)] if a else []) + ([('IXX', b), ('IYY', b)] if b else [])
        F = colouring(terms, deph)
        ok &= (F is None) == (a * b != 0)
        print(f"    (a, b) = ({a:2d}, {b:2d}): [H, U] = 0, U^2 = {a * a + b * b} I, anticommutes with the jumps; "
              f"colouring {F}")
    check("U = a YYZ + b ZXX carries the palindrome at every tested weight, exactly, and no single string "
          "does when ab != 0", ok)
    H = op('XXI') + op('YYI') + op('IXX') + op('IYY')
    sN, sW = nullities(H, jumps)
    nN, nW = int(np.sum(sN < 1e-9)), int(np.sum(sW < 1e-9))
    check("at a = b = 1 the far kernel is one-dimensional and the two counts agree", nW == 1 and nN == nW,
          f"dim N = {nN}, dim W = {nW}")


if __name__ == "__main__":
    stage_a()
    stage_b()
    stage_c()
    stage_d()
    stage_e()
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL STAGES PASS")
