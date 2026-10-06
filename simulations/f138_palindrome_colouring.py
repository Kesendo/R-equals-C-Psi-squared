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
  F  With commuting jumps, every palindrome has a lit string commuting with H in a frame the
     dissipator cannot see. For commuting single-letter
     jumps, S a Hermitian unitary in W, F any lit string and P = (1 + A)/2 for any one jump A, the
     unitary V = P + (1 - P) F S commutes with every jump and V^dagger F V = S, so F commutes with
     V H V^dagger. (i) EXACT (sympy): on stage E's row at (a, b) = (3, 4) and on a conditioned row
     (P3, ZZ bonds, jumps X on sites 0 and 2, fields 11/50 Z on site 1 and 41/100 Y on site 2, where
     U = Z (x) (P+ (x) (hY + Z) + P- (x) (hY - Z)) / sqrt(1 + h^2) takes its colour on site 2 from the
     sector of Z on the undephased site 1), for several F. (ii) the rows beyond the colouring of stage
     B, one bond set per letter count: a Hermitian unitary U is built from W by the polar and sign
     route (float), V as above with A the first jump and F the first lit string; gated: unitarity,
     [V, A] = 0 and [V H V^dagger, F] = 0 sit at least six decades below the control V = 1 (no lit
     string commutes with H on these rows). (iii) under uniform Z-dephasing V is diagonal, so the palindrome is a diagonal
     phase gauge to an X^N-symmetric Hamiltonian: on the 36 two-term bilinear pairs on the open chain
     at N = 3 (experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md: 22 palindromic, 14 hard;
     hypotheses/THE_OTHER_SIDE.md: 26 break the parity), every palindromic pair gets a diagonal V
     with [V H V^dagger, X^N] = 0, six decades below the parity breaking of H as written. (iv) without
     a jump that commutes with all the others the frame can fail: jumps X and Z on site 0 of three, the other sites undephased,
     an H whose W is spanned by Y (x) diag(1, 1, 1, -1); a frame image of a lit string has a second
     factor with eigenvalue split (2, 2) or (4, 0), so no frame exists although the palindrome holds.
  G  EXACT (sympy, sqrt 5). The golden router (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md, F116) is a
     colouring with colours on the lit circle: on the open chain with the sliding windows of
     XZX + XZY + YZX under Z-dephasing, G = (x)_l g_l with g_l in [a, a, b, b], a = phi X + Y,
     b = X - phi Y, commutes with H, anticommutes with every Z_l and squares to (1 + r^2)^N I
     (r = phi), at N = 3, 4, 5, so it lies in F158's far kernel; the same for the silver member
     (c = 2, r = 1 + sqrt 2) of the metallic family c XZX + XZY + YZX at N = 4; controls: the
     patterns [a, b, a, b], [a, a, a, a] and the 45 degree colour, and on the silver chain the golden
     colours and [a, b, a, b], do not commute with H. No
     colouring by the letters X and Y exists for c != 0 at N = 3..6, while at c = 0 exactly the four
     period-2 strings XX.., XYXY.., YXYX.., YY.. colour the chain.
  H  EXACT (Fraction arithmetic on Pauli strings, no float). Every row beyond the colouring of stage B
     (one bond set per letter count, 258 rows) has an element of W built from anticommuting sums: real
     combinations of pairwise anticommuting lit strings, tensored over the components of H and, where
     an undephased site carries a letter every term of H leaves alone (I or that letter), summed over
     its two sectors. [T, P] / 2i is a real rational combination of strings, so the kernel on each
     maximal anticommuting set of lit strings is an exact rational nullspace. Counted: 22 / 104 at two
     letters and 62 / 18 / 52 at one, by graph, and 74 single sums, 156 products, 28 conditioned by kind.
     Control: the same recursion with sets of size one (single strings, inside components and sectors
     too) finds nothing on any of these rows, and on every row of the same grid whose palindrome breaks
     the full recursion finds nothing; and every element found is rebuilt as a dense sympy matrix and
     checked against H and the jumps by matrix products, a route that does not use the string phase
     table (it would catch a wrong sign there, which the counts alone do not). (ii) the defect
     cascade on one bond, symbolic in the fields: H = XX + YY + h0 Y (x) I + h1 I (x) Z, jump X on the
     second site, C = ZZ - h0 IY - h0 h1 YZ commutes with H, anticommutes with the jump and squares
     to (1 + h0^2 + h0^2 h1^2) I, for all real h0, h1.
  I  EXACT (sympy ranks over the rationals). F138's clause 2 names a field axis; the axis is a letter,
     and a direction that is not a letter needs more (a rotation about the dephasing axis that turns it
     into a letter and keeps the bonds in the class, as for rotation-invariant bonds). Two sites, Z dephasing on both, fields 30 (X + Y) and 22 (X + Y): with bonds
     100 (XX + ZZ) the far kernel is empty while ker L is not (dimensions 1 and 0, no palindrome by
     F158); with the Heisenberg bond 100 (XX + YY + ZZ) the two dimensions agree (1 and 1).
  J  EXACT (sympy, symbolic). What the cascade predicts. (i) The carrier moves with the fields: the old C
     at fields (h, k) fails at (h + d, k) with [H, C] = 2i d XZ, while C(h + d, k) lies in W. (ii) A field
     e X on the first site is an obstruction: F158 (f5) makes every word with an odd number of jump letters
     traceless, and Tr(H^2 A) = 8e; the commutator on the lit strings has determinant e^2 on eight output
     strings at every h, k, so W = 0 for e != 0 (and W is 1-dimensional at e = 0). (iii) On stage E's row the adjoint dynamics L^dag(O) = i[H, O] + sum_l gamma_l (A_l O A_l - O)
     closes on P = YYZ, Q = ZXX, R = YIX: L^dag P = -2 sigma P - 2b R, L^dag Q = -2 sigma Q + 2a R,
     L^dag R = 2b P - 2a Q - 2(gamma_0 + gamma_2) R, so only a P + b Q decays with the single rate 2 sigma.

Run:  python simulations/f138_palindrome_colouring.py
   >  simulations/results/f138_palindrome_colouring.txt     (runtime about 30 minutes)
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


def op2(letters):
    return np.kron(LET[letters[0]], LET[letters[1]])


def lit_strings(deph, n=N):
    allowed = [[P for P in 'XYZ' if P not in deph[l]] if deph[l] else list('IXYZ') for l in range(n)]
    return [''.join(p) for p in itertools.product(*allowed)]


def hermitian_unitary_in_W(H, deph):
    """a Hermitian unitary element of W (float), from any invertible element: polar part W0, then
    K = e^{it} W0 + h.c. is Hermitian and invertible for generic t, and sign(K) lies in W"""
    S = lit_strings(deph)
    Ms = [op(s) for s in S]
    M = np.array([(H @ m - m @ H).ravel() for m in Ms]).T
    _, sv, vh = np.linalg.svd(M)
    basis = [sum(c * m for c, m in zip(v.conj(), Ms)) for v in vh[int(np.sum(sv > 1e-9)):]]
    rng = np.random.default_rng(7)
    U0 = sum(rng.normal() * b for b in basis)
    u, s, w = np.linalg.svd(U0)
    if s.min() < 1e-9:
        return None
    W0 = u @ w
    for t in rng.uniform(0, np.pi, size=20):
        K = np.exp(1j * t) * W0 + np.exp(-1j * t) * W0.conj().T
        ev, Q = np.linalg.eigh(K)
        if np.abs(ev).min() > 1e-6:
            return Q @ np.diag(np.sign(ev)) @ Q.conj().T
    return None


def frame(U, F, A):
    P = (np.eye(len(U)) + A) / 2
    return P + (np.eye(len(U)) - P) @ F @ U


def stage_f():
    import sympy as sp
    print()
    print("## Stage F: in a frame the dissipator cannot see, a lit string commutes with H")
    SL = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
          'Z': sp.diag(1, -1)}

    def sop(s):
        M = sp.Matrix([[1]])
        for c in s:
            M = sp.kronecker_product(M, SL[c])
        return M
    E, Z0 = sp.eye(d), sp.zeros(d)
    h1, h2 = sp.Rational(11, 50), sp.Rational(41, 100)
    Pp, Pm = sp.diag(1, 0), sp.diag(0, 1)
    rows = {
        "stage E's row, (a, b) = (3, 4)": (3 * (sop('XXI') + sop('YYI')) + 4 * (sop('IXX') + sop('IYY')),
                                           (3 * sop('YYZ') + 4 * sop('ZXX')) / 5,
                                           [sop('XII'), sop('IZI'), sop('IIY')], ['YYZ', 'ZXX', 'YXX', 'ZYX']),
        "conditioned row": (sop('ZZI') + sop('IZZ') + h1 * sop('IZI') + h2 * sop('IIY'),
                            (sp.kronecker_product(sp.kronecker_product(SL['Z'], Pp), h2 * SL['Y'] + SL['Z'])
                             + sp.kronecker_product(sp.kronecker_product(SL['Z'], Pm), h2 * SL['Y'] - SL['Z']))
                            / sp.sqrt(1 + h2 ** 2),
                            [sop('XII'), sop('IIX')], ['ZIZ', 'YXY', 'ZXY', 'YIZ']),
    }
    ok = True
    for name, (H, U, jumps, Fs) in rows.items():
        base = all(sp.simplify(x) == Z0 for x in (H * U - U * H, U * U - E, U - U.H)) and all(
            sp.simplify(A * U + U * A) == Z0 for A in jumps)
        for A in jumps:
            P = (E + A) / 2
            for f in Fs:
                F, V = sop(f), P + (E - P) * sop(f) * U
                ok &= base and all(sp.simplify(x) == Z0 for x in (
                    V * V.H - E, V.H * F * V - U, V * H * V.H * F - F * V * H * V.H)) and all(
                    sp.simplify(V * B - B * V) == Z0 for B in jumps)
        print(f"    {name}: S Hermitian unitary in W, frame checked for every jump as A and F in {Fs}")
    check("(i) V = P + (1 - P) F S is unitary, commutes with every jump, V^dagger F V = S and "
          "[V H V^dagger, F] = 0 (exact)", ok)
    no_string = colouring(terms_of(GRAPHS['P3'], ('Z',), ('I', 'Z', 'Y')), (('X',), (), ('X',))) is None
    check("(i'') the conditioned row has no colouring: no single lit string commutes with its H", no_string)
    G = sop('XZY')
    V = (2 * E + sp.I * G) / sp.sqrt(5)   # e^{i theta G / 2} with cos(theta) = 3/5
    _, U, jumps, _ = rows["stage E's row, (a, b) = (3, 4)"]
    check("(i') on stage E's row one frame is the rotation about the jump string XZY: V = e^{i theta XZY/2}, "
          "cos theta = a/r, commutes with every jump and V^dagger YYZ V = (a YYZ + b ZXX)/r (exact)",
          sp.simplify(V * V.H - E) == Z0 and all(sp.simplify(V * A - A * V) == Z0 for A in jumps)
          and sp.simplify(V.H * sop('YYZ') * V - U) == Z0)

    single = [(), ('X',), ('Y',), ('Z',)]
    worst, control, rows_done = 0.0, np.inf, collections.Counter()
    for gname, edges in GRAPHS.items():
        for bset in (('Z',), ('X', 'Y')):
            for deph in itertools.product(single, repeat=N):
                if not any(deph):
                    continue
                for fields in itertools.product('IXYZ', repeat=N):
                    terms = terms_of(edges, bset, fields)
                    if f138_clauses(edges, bset, deph, fields) or colouring(terms, deph) is not None:
                        continue
                    H = H_of(terms)
                    jumps = [op(placed({l: deph[l][0]})) for l in range(N) if deph[l]]
                    sN, sW = nullities(H, jumps)
                    if int(np.sum(sN < 1e-9)) != int(np.sum(sW < 1e-9)):
                        continue
                    U = hermitian_unitary_in_W(H, deph)
                    if U is None:
                        worst = np.inf
                        continue
                    f = lit_strings(deph)[0]
                    F = op(f)
                    V = frame(U, F, jumps[0])
                    Hv = V @ H @ V.conj().T
                    worst = max(worst, np.abs(V @ V.conj().T - Id).max(),
                                max(np.abs(V @ A - A @ V).max() for A in jumps), np.abs(Hv @ F - F @ Hv).max())
                    control = min(control, np.abs(H @ F - F @ H).max())
                    rows_done[(gname, len(bset))] += 1
    print(f"    rows per graph and bond-letter count: {dict(rows_done)}")
    expect = {('P3', 1): 62, ('K3', 1): 18, ('bond+iso', 1): 52, ('P3', 2): 22, ('bond+iso', 2): 104}
    dec = np.log10(control / worst)
    check("(ii) every row beyond the colouring (62 / 18 / 52 at one letter, 22 / 104 at two; the ZZ and XX + YY "
          "bond sets, A the first jump, F the first lit string) has a lit string commuting with H in a jump-commuting frame: the "
          "residuals sit at least six decades below the control V = 1, where by the row's definition F does not "
          "commute with H", dict(rows_done) == expect and dec >= 6,
          f"largest residual {worst:.1e}, smallest ||[H, F]|| {control:.2f}, {dec:.1f} decades")

    XN = op('XXX')
    Zj = [op(placed({l: 'Z'})) for l in range(N)]
    pal = hard = parity_breaking = 0
    worst, offdiag, control = 0.0, 0.0, np.inf
    kept_min, dropped_max = np.inf, 0.0
    for t1, t2 in itertools.combinations([a + b for a in 'XYZ' for b in 'XYZ'], 2):
        H = sum(op(placed({i: t[0], i + 1: t[1]})) for i in range(N - 1) for t in (t1, t2))
        breaks = np.abs(H @ XN - XN @ H).max()   # integer matrices: exact
        parity_breaking += bool((H @ XN - XN @ H).any())
        sN, sW = nullities(H, Zj)
        for sv in (sN, sW):
            kept_min = min(kept_min, sv[sv >= 1e-9].min())
            if np.any(sv < 1e-9):
                dropped_max = max(dropped_max, sv[sv < 1e-9].max())
        if int(np.sum(sN < 1e-9)) != int(np.sum(sW < 1e-9)):
            hard += 1
            continue
        pal += 1
        if breaks > 0:
            control = min(control, breaks)
        U = hermitian_unitary_in_W(H, (('Z',),) * N)
        if U is None:
            worst = np.inf
            continue
        V = frame(U, XN, Zj[0])
        Hv = V @ H @ V.conj().T
        offdiag = max(offdiag, np.abs(V - np.diag(np.diag(V))).max())
        worst = max(worst, np.abs(Hv @ XN - XN @ Hv).max(), np.abs(V @ V.conj().T - Id).max())
    dec = np.log10(control / max(worst, offdiag))
    dec_sv = np.log10(kept_min / max(dropped_max, 1e-300))
    check("(iii) uniform Z-dephasing, the 36 two-term pairs: 22 palindromic and 14 hard (the Klein routing's "
          "counts), 26 break the parity, and every palindromic pair is a diagonal phase gauge to "
          "[H', X^N] = 0: residual and off-diagonal part of V at least six decades below the parity breaking "
          "of the palindromic parity-breakers as written; the palindrome verdicts' singular values separate by six decades",
          (pal, hard, parity_breaking) == (22, 14, 26) and dec >= 6 and dec_sv >= 6,
          f"{pal} palindromic, {hard} hard, {parity_breaking} parity-breaking; largest residual {worst:.1e}, "
          f"largest off-diagonal entry of V {offdiag:.1e}, smallest ||[H, X^N]|| among the palindromic "
          f"parity-breakers {control:.2f}, {dec:.1f} decades; singular values {dec_sv:.1f} decades apart")

    # (iv) without a jump that commutes with all the others the frame can fail
    X, Z, Y = LET['X'], LET['Z'], LET['Y']
    e = np.eye(4)
    B1 = np.outer(e[0], e[3]) + np.outer(e[3], e[0])
    B2 = np.outer(e[1] + e[2], e[3]) + np.outer(e[3], e[1] + e[2])
    H = np.kron(np.eye(2), np.diag([1, 2, 3, 4])) + np.kron(X, B1) + np.kron(Z, B2)
    jumps = [np.kron(X, np.eye(4)), np.kron(Z, np.eye(4))]
    w = np.diag([1, 1, 1, -1])
    Uw = np.kron(Y, w)
    exact = (not (H @ Uw - Uw @ H).any()) and all(not (A @ Uw + Uw @ A).any() for A in jumps)
    sN, sW = nullities(H, jumps)
    nN, nW = int(np.sum(sN < 1e-9)), int(np.sum(sW < 1e-9))
    sv = np.concatenate([sN, sW])
    dec_sv = np.log10(sv[sv >= 1e-9].min() / max(sv[sv < 1e-9].max(), 1e-300))
    # the unitaries commuting with both jumps: the commutant of {X, Z} on site 0 is 1 (x) M_4, exactly
    units = [np.kron(np.eye(2), np.outer(e[i], e[j])) for i in range(4) for j in range(4)]
    others = [np.kron(P_, np.outer(e[i], e[j])) for P_ in (X, Y, Z) for i in range(4) for j in range(4)]
    comm = lambda M: np.concatenate([(M @ A - A @ M).ravel() for A in jumps])
    commutant_ok = all(not comm(M).any() for M in units) and np.linalg.matrix_rank(
        np.array([comm(M) for M in others]).T) == len(others)
    split_w = int(np.sum(np.diag(w) > 0))
    spectra = {tuple(sorted(np.round(np.linalg.eigvalsh(op2(f)), 9))) for f in itertools.product('IXYZ', repeat=2)}
    check("(iv) jumps X and Z on site 0, sites 1 and 2 undephased: the palindrome holds (dim N = dim W = 1), W is "
          "spanned by Y (x) diag(1, 1, 1, -1) (exact), while every frame image of a lit string is Y (x) v^dagger f v "
          "with f a two-site Pauli string (the unitaries commuting with both jumps are 1 (x) v: the commutant is "
          "1 (x) M_4, the 16 units commute exactly and the 48 X/Y/Z (x) units are independent under the commutator), "
          "whose eigenvalue split is (2, 2) or (4, 0), while diag(1, 1, 1, -1) splits (3, 1): no frame",
          exact and nN == nW == 1 and dec_sv >= 6 and commutant_ok and split_w == 3
          and all(sum(x > 0 for x in sp_) in (2, 4) for sp_ in spectra),
          f"dim N = {nN}, dim W = {nW}, singular values {dec_sv:.1f} decades apart; split of w ({split_w}, "
          f"{4 - split_w}); splits of the two-site strings {sorted({int(sum(x > 0 for x in s_)) for s_ in spectra})}")


def stage_g():
    import sympy as sp
    print()
    print("## Stage G: the golden router's G is a colouring by letters on the lit circle")
    SL = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
          'Z': sp.diag(1, -1)}

    def skron(ms):
        M = sp.Matrix([[1]])
        for m in ms:
            M = sp.kronecker_product(M, m)
        return M

    def window_H(n, c):
        H = sp.zeros(2 ** n)
        for w in range(n - 2):
            for t, coef in (('XZX', c), ('XZY', 1), ('YZX', 1)):
                H += coef * skron([SL['I']] * w + [SL[x] for x in t] + [SL['I']] * (n - w - 3))
        return H

    ok = True
    for n, c, r in ((3, 1, (1 + sp.sqrt(5)) / 2), (4, 1, (1 + sp.sqrt(5)) / 2), (5, 1, (1 + sp.sqrt(5)) / 2),
                    (4, 2, 1 + sp.sqrt(2))):
        a, b = r * SL['X'] + SL['Y'], SL['X'] - r * SL['Y']
        G = skron([[a, a, b, b][l % 4] for l in range(n)])
        H = window_H(n, c)
        Zs = [skron([SL['I']] * l + [SL['Z']] + [SL['I']] * (n - l - 1)) for l in range(n)]
        z0 = sp.zeros(2 ** n)
        this = (sp.simplify(H * G - G * H) == z0 and all(Zl * G + G * Zl == z0 for Zl in Zs)
                and sp.simplify(G * G - (1 + r ** 2) ** n * sp.eye(2 ** n)) == z0)
        ok &= this
        print(f"    N = {n}, c = {c}, r = {r}: [H, G] = 0, {{Z_l, G}} = 0, G^2 = (1 + r^2)^N I: {this}")
    phi = (1 + sp.sqrt(5)) / 2
    H4 = window_H(4, 1)
    controls = []
    H4s = window_H(4, 2)
    for Hc, r, pat in ((H4, phi, 'abab'), (H4, phi, 'aaaa'), (H4, 1, 'aabb'), (H4s, phi, 'aabb'),
                       (H4s, 1 + sp.sqrt(2), 'abab')):
        a, b = r * SL['X'] + SL['Y'], SL['X'] - r * SL['Y']
        Gc = skron([{'a': a, 'b': b}[pat[l]] for l in range(4)])
        controls.append(any(sp.simplify(x).equals(0) is False for x in (Hc * Gc - Gc * Hc)))
    check("(i) G = (x) [a, a, b, b] lies in the far kernel of the golden (N = 3, 4, 5) and silver (N = 4) window "
          "chains (exact); the controls [a, b, a, b], [a, a, a, a] and the 45 degree colour on the golden chain, and the "
          "golden colours and the pattern [a, b, a, b] on the silver chain, fail at N = 4 (an entry of [H, G] is exactly nonzero)",
          ok and all(controls))

    def commutes(s, t):
        return sum(x != 'I' and y != 'I' and x != y for x, y in zip(s, t)) % 2 == 0
    none_found, at_zero = True, True
    for n in (3, 4, 5, 6):
        for tpl in (('XZX', 'XZY', 'YZX'), ('XZY', 'YZX')):
            terms = [('I' * w + t + 'I' * (n - w - 3)) for w in range(n - 2) for t in tpl]
            sols = {''.join(F) for F in itertools.product('XY', repeat=n) if all(commutes(''.join(F), t) for t in terms)}
            if len(tpl) == 3:
                none_found &= not sols
            else:
                at_zero &= sols == {('XY' * n)[:n], ('YX' * n)[:n], 'X' * n, 'Y' * n}
    check("(ii) for c != 0 no colouring by the letters X and Y: every one of the 2^N strings fails a window "
          "template, N = 3..6; at c = 0 exactly the four period-2 strings colour the chain", none_found and at_zero)


_MUL = {('I', 'I'): (1, 'I'), ('I', 'X'): (1, 'X'), ('I', 'Y'): (1, 'Y'), ('I', 'Z'): (1, 'Z'),
        ('X', 'I'): (1, 'X'), ('X', 'X'): (1, 'I'), ('X', 'Y'): (1j, 'Z'), ('X', 'Z'): (-1j, 'Y'),
        ('Y', 'I'): (1, 'Y'), ('Y', 'X'): (-1j, 'Z'), ('Y', 'Y'): (1, 'I'), ('Y', 'Z'): (1j, 'X'),
        ('Z', 'I'): (1, 'Z'), ('Z', 'X'): (1j, 'Y'), ('Z', 'Y'): (-1j, 'X'), ('Z', 'Z'): (1, 'I')}


def _smul(s, t):
    ph, out = 1, []
    for a, b in zip(s, t):
        p, c = _MUL[(a, b)]
        ph *= p
        out.append(c)
    return ph, ''.join(out)


def _anti(s, t):
    return sum(x != 'I' and y != 'I' and x != y for x, y in zip(s, t)) % 2 == 1


def _rank_q(rows):
    """rank of a list of dict-rows over the rationals (Fraction entries)"""
    from fractions import Fraction
    rows = [dict(r) for r in rows if any(v != 0 for v in r.values())]
    rank, pivots = 0, []
    for r in rows:
        for (k, pr) in pivots:
            if r.get(k, 0) != 0:
                f = r[k] / pr[k]
                for kk, vv in pr.items():
                    r[kk] = r.get(kk, Fraction(0)) - f * vv
        nz = [k for k, v in r.items() if v != 0]
        if nz:
            pivots.append((nz[0], r))
            rank += 1
    return rank


def _clifford_kernel(terms, deph, max_size=None):
    """does W contain a real combination of pairwise anticommuting lit strings? terms: {string: Fraction}"""
    n = len(deph)
    allowed = [[P for P in 'XYZ' if P not in deph[l]] if deph[l] else list('IXYZ') for l in range(n)]
    S = [''.join(p) for p in itertools.product(*allowed)]
    adj = {s: {t for t in S if t != s and _anti(s, t)} for s in S}
    cliques = []

    def bk(R, P, X):
        if not P and not X:
            cliques.append(sorted(R))
            return
        for v in list(P):
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}
    bk(set(), set(S), set())
    if max_size is not None:
        cliques = [[s] for s in S]
    for C in cliques:
        # column for each string P in C: ad_H(P) / 2i as {output string: rational}; kernel iff rank < |C|
        cols = []
        for P in C:
            col = {}
            for T, c in terms.items():
                if _anti(T, P):
                    ph, Q = _smul(T, P)
                    col[Q] = col.get(Q, 0) + c * (1 if ph == 1j else -1)
            cols.append(col)
        if _rank_q(cols) < len(C):
            import sympy as sp
            rows = sorted({q for col in cols for q in col if col[q] != 0})
            if not rows:
                return {C[0]: 1}
            M = sp.Matrix([[sp.Rational(col.get(q, 0).numerator, col.get(q, 0).denominator)
                            if col.get(q, 0) != 0 else 0 for col in cols] for q in rows])
            v = M.nullspace()[0]
            return {P: v[i] for i, P in enumerate(C) if v[i] != 0}
    return None


def _clifford_explained(terms, deph, max_size=None):
    out = _clifford_element(terms, deph, max_size)
    return None if out is None else out[0]


def _clifford_element(terms, deph, max_size=None):
    """(kind, element) with the element a dict {string: rational}, built from anticommuting sums"""
    from fractions import Fraction
    n = len(deph)
    terms = {t: c for t, c in terms.items() if c != 0 and set(t) != {'I'}}
    parent = list(range(n))

    def f(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for t in terms:
        sites = [i for i, x in enumerate(t) if x != 'I']
        for a in sites[1:]:
            parent[f(a)] = f(sites[0])
    comps = collections.defaultdict(list)
    for i in range(n):
        comps[f(i)].append(i)
    comps = list(comps.values())
    if len(comps) > 1:
        elem = {'I' * n: 1}
        for c in comps:
            sub = {''.join(t[i] for i in c): v for t, v in terms.items() if all(t[i] == 'I' for i in range(n) if i not in c)}
            part = _clifford_element(sub, tuple(deph[i] for i in c), max_size)
            if part is None:
                return None
            new = {}
            for s0, c0 in elem.items():
                for s1, c1 in part[1].items():
                    lst = list(s0)
                    for k, i in enumerate(c):
                        lst[i] = s1[k]
                    new[''.join(lst)] = c0 * c1
            elem = new
        return 'product', elem
    kern = _clifford_kernel(terms, deph, max_size)
    if kern:
        return 'sum', kern
    for u in range(n):
        if deph[u]:
            continue
        for P in 'XYZ':
            if all(t[u] in ('I', P) for t in terms):
                subs = []
                for s in (1, -1):
                    sub = {}
                    for t, c in terms.items():
                        key = t[:u] + t[u + 1:]
                        sub[key] = sub.get(key, Fraction(0)) + c * (s if t[u] == P else 1)
                    subs.append(_clifford_element(sub, deph[:u] + deph[u + 1:], max_size))
                if None not in subs:
                    import sympy as sp
                    elem = {}
                    for sgn, part in zip((1, -1), subs):
                        for t, c in part[1].items():
                            for letter, w in (('I', sp.Rational(1, 2)), (P, sp.Rational(sgn, 2))):
                                key = t[:u] + letter + t[u:]
                                elem[key] = elem.get(key, 0) + w * c
                    return 'conditioned', {k: v for k, v in elem.items() if v != 0}
    return None


def stage_h():
    import sympy as sp
    from fractions import Fraction
    print()
    print("## Stage H: every row beyond the colouring is built from anticommuting sums (exact)")
    mag = [Fraction(30, 100), Fraction(22, 100), Fraction(41, 100)]
    single = [(), ('X',), ('Y',), ('Z',)]
    tally = collections.Counter()
    control_hits = broken_rows = broken_hits = verified = 0
    _SL = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
           'Z': sp.diag(1, -1)}

    def _sym_op(t):
        M = sp.Matrix([[1]])
        for ch in t:
            M = sp.kronecker_product(M, _SL[ch])
        return M
    for gname, edges in GRAPHS.items():
        for bset in (('Z',), ('X', 'Y')):
            for deph in itertools.product(single, repeat=N):
                if not any(deph):
                    continue
                for fields in itertools.product('IXYZ', repeat=N):
                    fterms = terms_of(edges, bset, fields)
                    if f138_clauses(edges, bset, deph, fields) or colouring(fterms, deph) is not None:
                        continue
                    H = H_of(fterms)
                    jumps = [op(placed({l: deph[l][0]})) for l in range(N) if deph[l]]
                    sN, sW = nullities(H, jumps)
                    terms = {placed({a: P, b: P}): Fraction(1) for (a, b) in edges for P in bset}
                    for l, P in enumerate(fields):
                        if P != 'I':
                            terms[placed({l: P})] = mag[l]
                    if int(np.sum(sN < 1e-9)) != int(np.sum(sW < 1e-9)):
                        # control that can fail: a row whose palindrome breaks must have no such element
                        broken_rows += 1
                        broken_hits += _clifford_explained(terms, deph) is not None
                        continue
                    found = _clifford_element(terms, deph)
                    kind = None if found is None else found[0]
                    tally[(gname, len(bset), kind)] += 1
                    if found is not None:
                        # independent route: dense sympy matrices, no string phase table
                        Gm = sum((c * _sym_op(t) for t, c in found[1].items()), sp.zeros(d))
                        Hm = sum((sp.Rational(c.numerator, c.denominator) * _sym_op(t) for t, c in terms.items()),
                                 sp.zeros(d))
                        Am = [_sym_op(placed({l: deph[l][0]})) for l in range(N) if deph[l]]
                        Z0 = sp.zeros(d)
                        verified += bool(sp.expand(Hm * Gm - Gm * Hm) == Z0
                                         and all(sp.expand(A * Gm + Gm * A) == Z0 for A in Am)
                                         and sp.simplify(Gm.det()) != 0)
                    # control: the same recursion with single strings only
                    control_hits += _clifford_explained(terms, deph, max_size=1) is not None
    for k in sorted(tally, key=str):
        print(f"    {k[0]:9s} {k[1]} bond letter(s), kind {k[2]}: {tally[k]}")
    per = collections.Counter()
    for (gname, nl, kind), v in tally.items():
        per[(gname, nl)] += v
    expect = {('P3', 1): 62, ('K3', 1): 18, ('bond+iso', 1): 52, ('P3', 2): 22, ('bond+iso', 2): 104}
    none = sum(v for (gname, nl, kind), v in tally.items() if kind is None)
    by_kind = collections.Counter()
    for (gname, nl, kind), v in tally.items():
        by_kind[kind] += v
    check("(i) every one of the 258 rows beyond the colouring has an element of W built from anticommuting sums, "
          "74 single sums, 156 products, 28 conditioned "
          "(exact rational kernels, kinds labelled by the first the recursion finds); controls: on every row of the "
          "same grid whose palindrome breaks the recursion finds nothing, and with single strings only it finds "
          "nothing on the 258; every element rebuilt densely and verified exactly", dict(per) == expect and none == 0 and control_hits == 0 and broken_rows > 0
          and broken_hits == 0 and dict(by_kind) == {'sum': 74, 'product': 156, 'conditioned': 28}
          and verified == 258,
          f"kinds {dict((k, v) for k, v in tally.items())}; rows without: {none}; single-string hits: {control_hits}; "
          f"broken rows {broken_rows}, hits there {broken_hits}; elements rebuilt as dense sympy matrices and "
          f"verified ([H, G] = 0, {{A, G}} = 0, det G != 0): {verified}")
    h0, h1 = sp.symbols('h0 h1', real=True)
    SX, SY, SZ, SI = (sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1), sp.eye(2))
    k = sp.kronecker_product
    Hb = k(SX, SX) + k(SY, SY) + h0 * k(SY, SI) + h1 * k(SI, SZ)
    Cb = k(SZ, SZ) - h0 * k(SI, SY) - h0 * h1 * k(SY, SZ)
    ok = (sp.simplify(Hb * Cb - Cb * Hb) == sp.zeros(4) and sp.simplify(k(SI, SX) * Cb + Cb * k(SI, SX)) == sp.zeros(4)
          and sp.simplify(Cb * Cb - (1 + h0 ** 2 + h0 ** 2 * h1 ** 2) * sp.eye(4)) == sp.zeros(4))
    Cbad = k(SZ, SZ) - h0 * k(SI, SY) + h0 * h1 * k(SY, SZ)
    ok &= sp.simplify(Hb * Cbad - Cbad * Hb) != sp.zeros(4)
    check("(ii) the defect cascade on one bond, symbolic: C = ZZ - h0 IY - h0 h1 YZ commutes with "
          "XX + YY + h0 Y0 + h1 Z1, anticommutes with the jump X1, C^2 = (1 + h0^2 + h0^2 h1^2) I; the sign-flipped "
          "last coefficient fails", ok)


def stage_i():
    import sympy as sp
    print()
    print("## Stage I: the field axis of clause 2 is a letter; another direction needs the bonds to follow it")
    SL = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
          'Z': sp.diag(1, -1)}

    def sop(t):
        return sp.kronecker_product(SL[t[0]], SL[t[1]])

    def dims(H):
        """dim ker L and dim ker(L + 2 sigma) by F158's Lemma 1: ker ad_H on the strings commuting with
        both Z jumps ({I, Z} per site) and on the lit ones ({X, Y} per site), exact rational ranks"""
        out = []
        for letters in ('IZ', 'XY'):
            cols = [sop(a + b) for a in letters for b in letters]
            M = sp.Matrix.hstack(*[(H * C - C * H).reshape(16, 1) for C in cols])
            out.append(len(cols) - M.rank())
        return out
    fields = 30 * (sop('XI') + sop('YI')) + 22 * (sop('IX') + sop('IY'))
    reduced = dims(100 * (sop('XX') + sop('ZZ')) + fields)
    heis = dims(100 * (sop('XX') + sop('YY') + sop('ZZ')) + fields)
    check("bonds XX + ZZ under Z dephasing, fields along X + Y: dim ker L = 1, dim ker(L + 2 sigma) = 0, no "
          "palindrome; the Heisenberg bond with the same fields: 1 and 1", reduced == [1, 0] and heis == [1, 1],
          f"XX + ZZ {reduced}, Heisenberg {heis}")


def stage_j():
    import sympy as sp
    print()
    print("## Stage J: what the cascade predicts (the carrier moves, an obstruction, one decay rate)")
    SL = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
          'Z': sp.diag(1, -1)}

    def sop(t):
        M = sp.Matrix([[1]])
        for ch in t:
            M = sp.kronecker_product(M, SL[ch])
        return M
    h, k, dd, e = sp.symbols('h k d e', real=True)
    Hb = lambda hh, kk, ee=0: sop('XX') + sop('YY') + hh * sop('YI') + kk * sop('IZ') + ee * sop('XI')
    C = lambda hh, kk: sop('ZZ') - hh * sop('IY') - hh * kk * sop('YZ')
    A = sop('IX')
    z4 = sp.zeros(4)
    moves = (sp.expand(Hb(h + dd, k) * C(h, k) - C(h, k) * Hb(h + dd, k) - 2 * sp.I * dd * sop('XZ')) == z4
             and sp.expand(Hb(h + dd, k) * C(h + dd, k) - C(h + dd, k) * Hb(h + dd, k)) == z4)
    check("(i) the carrier moves with the fields: [H(h + d, k), C(h, k)] = 2i d XZ, and C(h + d, k) commutes with "
          "H(h + d, k)", moves)
    tr = sp.expand((Hb(h, k, e) ** 2 * A).trace())

    def dimW(H):
        cols = [sop(a + b) for a in 'IXYZ' for b in 'YZ']      # the strings anticommuting with IX
        M = sp.Matrix.hstack(*[(H * c - c * H).reshape(16, 1) for c in cols])
        return len(cols) - M.rank()
    dims = (dimW(Hb(sp.Rational(3, 10), sp.Rational(11, 50), 0)), dimW(Hb(sp.Rational(3, 10), sp.Rational(11, 50),
                                                                         sp.Rational(1, 3))))
    cols = [a + b for a in 'IXYZ' for b in 'YZ']
    outs = ['IZ', 'IY', 'XI', 'XY', 'YZ', 'ZY', 'XZ', 'IX']
    Hs = Hb(h, k, e)
    M8 = sp.Matrix([[sp.expand(((Hs * sop(c) - sop(c) * Hs) * sop(o)).trace() / (8 * sp.I)) for c in cols]
                    for o in outs])
    det = sp.factor(M8.det())
    check("(ii) a field e X on the first site: Tr(H^2 A) = 8e; the commutator on the eight lit strings has "
          "determinant e^2 on eight output strings at every h, k, so W = 0 for e != 0; W is 1-dimensional at e = 0",
          tr == 8 * e and det == e ** 2 and dims == (1, 0), f"Tr = {tr}, det = {det}, dims {dims}")
    a, b, g0, g1, g2 = sp.symbols('a b gamma0 gamma1 gamma2', positive=True)
    H = a * (sop('XXI') + sop('YYI')) + b * (sop('IXX') + sop('IYY'))
    jumps = [(g0, sop('XII')), (g1, sop('IZI')), (g2, sop('IIY'))]
    Ldag = lambda O: sp.I * (H * O - O * H) + sum((g * (J * O * J - O) for g, J in jumps), sp.zeros(8))
    P, Q, R = sop('YYZ'), sop('ZXX'), sop('YIX')
    sig = g0 + g1 + g2
    z8 = sp.zeros(8)
    closes = (sp.expand(Ldag(P) - (-2 * sig * P - 2 * b * R)) == z8 and sp.expand(Ldag(Q) - (-2 * sig * Q + 2 * a * R)) == z8
              and sp.expand(Ldag(R) - (2 * b * P - 2 * a * Q - 2 * (g0 + g2) * R)) == z8)
    U, V = a * P + b * Q, b * P - a * Q
    single = sp.expand(Ldag(U) + 2 * sig * U) == z8 and sp.expand(Ldag(V) + 2 * sig * V) != z8
    check("(iii) the adjoint dynamics closes on YYZ, ZXX, YIX as stated; a YYZ + b ZXX decays with the single rate "
          "2 sigma and b YYZ - a ZXX does not", closes and single)


if __name__ == "__main__":
    stage_a()
    stage_b()
    stage_c()
    stage_d()
    stage_e()
    stage_f()
    stage_g()
    stage_h()
    stage_i()
    stage_j()
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL STAGES PASS")
