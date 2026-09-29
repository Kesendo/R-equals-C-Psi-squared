#!/usr/bin/env python3
"""The anticommuting-sum census at N = 3 and N = 4 (experiments/THE_PALINDROME_AS_A_COLOURING.md,
section "What a sum is").

For Pauli jumps (at most one dephasing axis per site) F158's two kernels are kernels of ad_H on
string spans: ker L on the strings that commute with every jump, ker(L + 2 sigma) (the far kernel W)
on the lit strings, those that anticommute with every jump. So the palindrome, dim ker L =
dim ker(L + 2 sigma), is two ranks of integer matrices: (ad_H P) / 2i is +- a string times the term's
coefficient, with bonds weighted 100 and fields 30, 22, 41, 17 (the colouring page's census tuple,
scaled to integers, plus a fourth site). The ranks are taken modulo two primes near 2^31, which must
agree (a bad reduction can only raise a nullity).

Every palindromic row either has a colouring (a single lit string commuting with H) or, beyond the
colouring, an element of W built from anticommuting sums: real combinations of pairwise
anticommuting lit strings, tensored over the components of H and summed over the two sectors of a
letter that an undephased site keeps, recursively. Each such element is built exactly (rational
nullspace on the anticommuting set, sector weights 1/2, scaled to integers) and verified exactly:
[H, G] = 0 and {A, G} = 0 in integer matrix arithmetic (the products are bounded entrywise by
|H| |G| and |A| |G|, checked to stay below 2^50, where float64 holds every partial sum exactly), and G invertible because its rank modulo a prime
p = 1 mod 4, with i sent to a square root of -1, is full, which forces det G != 0.

Stages (all must pass; prints "ALL STAGES PASS"):
  A  N = 3, graphs P3, K3, a bond with an isolated site; the rows beyond the colouring reproduce the
     colouring page's stage H: 62 / 18 / 52 at ZZ, 22 / 0 / 104 at XX + YY, 0 at XX + YY + ZZ, by kind
     52 + 10, 18, 52 and 22, 104; every one verified exactly.
  B  N = 4, graphs P4, S4 (star), C4 (ring), K4, bonds ZZ, XX + YY, XX + YY + ZZ, every dephasing
     pattern with at most one axis per site and at least one jump, every field pattern (783,360 rows):
     the counts pinned as measured, no palindromic row without a colouring or an anticommuting-sum
     element, and every such element verified exactly. The coloured counts do not depend on the
     graph (13215 at ZZ, 3795 at XX + YY and at XX + YY + ZZ): on a connected graph every bond passes
     the same class (one letter) or the same colour (two or three letters) along.
  C  Controls. The same recursion finds nothing on the rows whose palindrome breaks (all of them at
     N = 3; at N = 4 those at a fixed set of field patterns, every pattern with at most two fields,
     three equal fields on any three sites, the uniform and the cyclic mixed fields on all four,
     whose coverage is checked), the palindromic rows per graph and bond set are pinned,
     and with single strings only (so only sector sums of single strings could still fire) it finds
     nothing beyond the colouring; on every row, no row with a colouring, which certifies the
     palindrome exactly, is called broken by the ranks; every graph and bond set is present once with
     all its rows. Before the census the verifier rejects three fixed false objects (one failing each
     duty) and accepts the SWAP certificate of F138's coincident-magnitude family, on which the
     building rule finds nothing: that palindrome lies outside the rule's grammar.
  The run records its revision and the sha256 of this script, and every element beyond the colouring
  is written out, row and rational coefficients, to results/anticommuting_sum_certificates.json (the
  clique search iterates in sorted order, so the file does not depend on the hash seed).

Run:  python simulations/anticommuting_sum_census.py
   >  simulations/results/anticommuting_sum_census.txt     (runtime about 25 minutes on 22 cores)
"""
import collections
import itertools
import sys
import time
from pathlib import Path
from fractions import Fraction
from multiprocessing import Pool

import numpy as np

PRIMES = (2147483629, 2147483587)          # both near 2^31; the first is 1 mod 4
P_GAUSS = PRIMES[0]
_NR = next(a for a in range(2, 100) if pow(a, (P_GAUSS - 1) // 2, P_GAUSS) == P_GAUSS - 1)
SQRT_MINUS_ONE = pow(_NR, (P_GAUSS - 1) // 4, P_GAUSS)   # a^((p-1)/4) for a non-residue a
assert SQRT_MINUS_ONE * SQRT_MINUS_ONE % P_GAUSS == P_GAUSS - 1
LET = 'IXYZ'
MAG = [30, 22, 41, 17]
MUL = {('I', 'I'): (1, 'I'), ('I', 'X'): (1, 'X'), ('I', 'Y'): (1, 'Y'), ('I', 'Z'): (1, 'Z'),
       ('X', 'I'): (1, 'X'), ('X', 'X'): (1, 'I'), ('X', 'Y'): (1j, 'Z'), ('X', 'Z'): (-1j, 'Y'),
       ('Y', 'I'): (1, 'Y'), ('Y', 'X'): (-1j, 'Z'), ('Y', 'Y'): (1, 'I'), ('Y', 'Z'): (1j, 'X'),
       ('Z', 'I'): (1, 'Z'), ('Z', 'X'): (1j, 'Y'), ('Z', 'Y'): (-1j, 'X'), ('Z', 'Z'): (1, 'I')}
PAULI = {'I': np.array([[1, 0], [0, 1]], complex), 'X': np.array([[0, 1], [1, 0]], complex),
         'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.array([[1, 0], [0, -1]], complex)}
GRAPHS = {3: {'P3': [(0, 1), (1, 2)], 'K3': [(0, 1), (1, 2), (0, 2)], 'bond+iso': [(0, 1)]},
          4: {'P4': [(0, 1), (1, 2), (2, 3)], 'S4': [(0, 1), (0, 2), (0, 3)],
              'C4': [(0, 1), (1, 2), (2, 3), (3, 0)], 'K4': [(a, b) for a in range(4) for b in range(a + 1, 4)]}}
BONDSETS = (('Z',), ('X', 'Y'), ('X', 'Y', 'Z'))

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def anti(s, t):
    return sum(a != 'I' and b != 'I' and a != b for a, b in zip(s, t)) % 2 == 1


def smul(s, t):
    ph, out = 1, []
    for a, b in zip(s, t):
        p, c = MUL[(a, b)]
        ph *= p
        out.append(c)
    return ph, ''.join(out)


def rank_mod(M, p):
    M = M.copy() % p
    r, (m, n) = 0, M.shape
    for c in range(n):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0:
            continue
        i = r + piv[0]
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        others = np.nonzero(M[:, c])[0]
        others = others[others != r]
        if len(others):
            M[others] = (M[others] - M[others, c].reshape(-1, 1) * M[r]) % p
        r += 1
        if r == m:
            break
    return r


def ad_matrix(terms, cols):
    """(ad_H P) / 2i for P in cols as an integer matrix; rows are the output strings"""
    rowidx, entries = {}, []
    for j, P in enumerate(cols):
        for T, c in terms:
            if anti(T, P):
                ph, Q = smul(T, P)
                entries.append((rowidx.setdefault(Q, len(rowidx)), j, c if ph == 1j else -c))
    M = np.zeros((max(len(rowidx), 1), len(cols)), dtype=np.int64)
    for i, j, v in entries:
        M[i, j] += v
    return M


def nullity(terms, cols):
    if not cols:
        return 0
    M = ad_matrix(terms, cols)
    rs = {rank_mod(M, p) for p in PRIMES}
    if len(rs) != 1:
        raise RuntimeError('the two primes disagree on a rank')
    return len(cols) - rs.pop()


def strings(allowed):
    return [''.join(p) for p in itertools.product(*allowed)]


def lit_of(deph):
    return strings([[P for P in 'XYZ' if P not in d] if d else list(LET) for d in deph])


def dark_of(deph):
    return strings([['I'] + list(d) if d else list(LET) for d in deph])


def exact_kernel_vector(terms, cols):
    """a nonzero rational vector in the kernel of ad_H on span(cols)"""
    import sympy as sp
    M = ad_matrix(terms, cols)
    if not M.any():
        return {cols[0]: Fraction(1)}
    v = sp.Matrix(M.tolist()).nullspace()[0]
    return {P: Fraction(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for P, x in zip(cols, v) if x != 0}


def clique_element(terms, deph, max_size=None):
    lit = lit_of(deph)
    if max_size == 1:
        for s in lit:
            if nullity(terms, [s]) == 1:
                return {s: Fraction(1)}
        return None
    adj = {s: {t for t in lit if t != s and anti(s, t)} for s in lit}
    found = []

    def bk(R, P, X):
        if found:
            return
        if not P and not X:
            if nullity(terms, sorted(R)) > 0:
                found.append(sorted(R))
            return
        u = max(sorted(P | X), key=lambda v: len(adj[v] & P))     # sorted: the result must not
        for v in sorted(P - adj[u]):                                # depend on the hash seed
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}
    bk(set(), set(lit), set())
    return exact_kernel_vector(terms, found[0]) if found else None


def element(terms, deph, max_size=None):
    """(kind, {string: Fraction}) built from anticommuting sums, or None"""
    n = len(deph)
    terms = [(t, c) for t, c in terms if c != 0 and set(t) != {'I'}]
    parent = list(range(n))

    def f(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for t, _ in terms:
        sites = [i for i, x in enumerate(t) if x != 'I']
        for a in sites[1:]:
            parent[f(a)] = f(sites[0])
    comps = collections.defaultdict(list)
    for i in range(n):
        comps[f(i)].append(i)
    comps = list(comps.values())
    if len(comps) > 1:
        elem = {'I' * n: Fraction(1)}
        for c in comps:
            sub = [(''.join(t[i] for i in c), v) for t, v in terms if all(t[i] == 'I' for i in range(n) if i not in c)]
            part = element(sub, tuple(deph[i] for i in c), max_size)
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
    e = clique_element(terms, deph, max_size)
    if e is not None:
        return 'sum', e
    for u in range(n):
        if deph[u]:
            continue
        for P in 'XYZ':
            if all(t[u] in ('I', P) for t, _ in terms):
                parts = []
                for s in (1, -1):
                    acc = collections.defaultdict(int)
                    for t, c in terms:
                        acc[t[:u] + t[u + 1:]] += c * (s if t[u] == P else 1)
                    part = element(list(acc.items()), deph[:u] + deph[u + 1:], max_size)
                    if part is None:
                        break
                    parts.append(part)
                if len(parts) == 2:
                    elem = collections.defaultdict(Fraction)
                    for s, part in zip((1, -1), parts):
                        for t, c in part[1].items():
                            elem[t[:u] + 'I' + t[u:]] += c / 2
                            elem[t[:u] + P + t[u:]] += s * c / 2
                    return 'conditioned', {k: v for k, v in elem.items() if v != 0}
    return None


def dense(s):
    M = np.array([[1]], complex)
    for ch in s:
        M = np.kron(M, PAULI[ch])
    return M


def verify(terms, deph, elem):
    """exact: [H, G] = 0, {A, G} = 0 in integer arithmetic, G invertible by full rank mod p"""
    den = 1
    for c in elem.values():
        den = den * c.denominator // np.gcd(den, c.denominator)
    G = sum(int(c * den) * dense(s) for s, c in elem.items())
    H = sum(c * dense(t) for t, c in terms)
    jumps = [dense(''.join(d[0] if d and k == l else 'I' for k, d in enumerate(deph))) for l in range(len(deph)) if deph[l]]
    aH, aG = np.abs(H), np.abs(G)
    bounds = [aH @ aG, aG @ aH] + [np.abs(A) @ aG for A in jumps] + [aG @ np.abs(A) for A in jumps]
    if max(M.max() for M in bounds) >= 2 ** 50:
        raise RuntimeError('partial sums could exceed exact float64 range')
    ok = not (H @ G - G @ H).any() and all(not (A @ G + G @ A).any() for A in jumps)
    s = SQRT_MINUS_ONE
    Gp = (np.rint(G.real).astype(np.int64) % P_GAUSS + s * (np.rint(G.imag).astype(np.int64) % P_GAUSS)) % P_GAUSS
    return ok and rank_mod(Gp, P_GAUSS) == G.shape[0]


def placed(n, d):
    return ''.join(d.get(k, 'I') for k in range(n))


def row_terms(n, edges, bset, fields):
    return ([(placed(n, {a: P, b: P}), 100) for (a, b) in edges for P in bset]
            + [(placed(n, {l: P}), MAG[l]) for l, P in enumerate(fields) if P != 'I'])


def control_fields(n):
    """the field patterns at which broken rows are searched for an element (None: all of them): every
    pattern with at most two fields (any letters), three equal fields on any three sites, the uniform and
    the cyclic mixed fields on all sites"""
    if n == 3:
        return None
    pats = set()
    for f in itertools.product(LET, repeat=n):
        support = [l for l in range(n) if f[l] != 'I']
        if len(support) <= 2 or (len(support) == 3 and len({f[l] for l in support}) == 1):
            pats.add(''.join(f))
    for P in 'XYZ':
        pats.add(P * n)
    for k in range(3):
        pats.add(''.join('XYZ'[(k + l) % 3] for l in range(n)))
    return frozenset(pats)


def work(args):
    n, gname, edges, bset, deph, controls = args
    out = collections.Counter()
    certs = []
    for fi, fields in enumerate(itertools.product(LET, repeat=n)):
        terms = row_terms(n, edges, bset, fields)
        lit = lit_of(deph)
        coloured = any(all(not anti(T, s) for T, _ in terms) for s in lit)
        if nullity(terms, dark_of(deph)) != nullity(terms, lit):
            out['broken'] += 1
            out['broken but coloured'] += coloured      # a colouring certifies the palindrome exactly
            if controls is None or ''.join(fields) in controls:
                out['broken checked'] += 1
                out['broken with element'] += element(terms, deph) is not None
            continue
        out['palindromic'] += 1
        if coloured:
            out['coloured'] += 1
            continue
        found = element(terms, deph)
        out['beyond ' + ('none' if found is None else found[0])] += 1
        if found is not None:
            out['verified'] += verify(terms, deph, found[1])
            certs.append({'n': n, 'graph': gname, 'bonds': ''.join(bset), 'dephasing': [''.join(d) for d in deph],
                          'fields': ''.join(fields), 'kind': found[0],
                          'element': {k: str(v) for k, v in sorted(found[1].items())}})
        out['single-string hits'] += element(terms, deph, max_size=1) is not None
    return (gname, ''.join(bset)), out, certs


def census(n, workers, certs):
    controls = control_fields(n)
    jobs = [(n, g, e, b, d, controls) for g, e in GRAPHS[n].items() for b in BONDSETS
            for d in itertools.product([(), ('X',), ('Y',), ('Z',)], repeat=n) if any(d)]
    tot = collections.defaultdict(collections.Counter)
    t0 = time.time()
    with Pool(workers) as pool:
        for key, c, cs in pool.imap_unordered(work, jobs, chunksize=4):
            tot[key].update(c)
            certs.extend(cs)
    for key in sorted(tot):
        print(f"    {key[0]:9s} {key[1]:4s} {dict(sorted(tot[key].items()))}")
    print(f"    ({time.time() - t0:.0f} s)")
    return tot


def beyond(c):
    return {k[len('beyond '):]: v for k, v in c.items() if k.startswith('beyond ')}


def summary_checks(tot, expect, label, per_group, palindromic, all_broken_checked):
    keys_ok = set(tot) == set(expect)
    pal = {k: c['palindromic'] for k, c in tot.items()}
    check(f"{label}: the palindromic rows per graph and bond set, as measured (a pin, so that a change is seen)",
          pal == palindromic, str(dict(sorted(pal.items()))))
    if all_broken_checked:
        check(f"{label}: every broken row searched for an element",
              all(c.get('broken checked', 0) == c['broken'] for c in tot.values()))
    rows = {k: c['palindromic'] + c['broken'] for k, c in tot.items()}
    check(f"{label}: every graph and bond set present exactly once, each with all its rows ({per_group})",
          keys_ok and all(v == per_group for v in rows.values()),
          f"groups {len(tot)} of {len(expect)}; row counts {sorted(set(rows.values()))}")
    ok_kinds = keys_ok and all(beyond(tot[k]) == v for k, v in expect.items())
    none = sum(c.get('beyond none', 0) for c in tot.values())
    nb = sum(sum(beyond(c).values()) for c in tot.values())
    ver = sum(c.get('verified', 0) for c in tot.values())
    check(f"{label}: the rows beyond the colouring by graph, bond set and kind, as pinned; none without an "
          f"element", ok_kinds and none == 0, f"{nb} rows beyond, {none} without")
    check(f"{label}: every element verified exactly ([H, G] = 0, {{A, G}} = 0 in integers, full rank mod p)",
          ver == nb, f"{ver} of {nb}")
    hits = sum(c.get('single-string hits', 0) for c in tot.values())
    bchk = sum(c.get('broken checked', 0) for c in tot.values())
    bhit = sum(c.get('broken with element', 0) for c in tot.values())
    bc = sum(c.get('broken but coloured', 0) for c in tot.values())
    check(f"{label}: no row with a colouring (an exact certificate of the palindrome) is called broken by the "
          f"ranks, over every row", bc == 0, f"{bc} rows")
    check(f"{label}: controls: with single strings only (sector sums of them) nothing beyond the colouring; "
          f"the recursion finds nothing on the broken rows checked", hits == 0 and bchk > 0 and bhit == 0,
          f"single-string hits {hits}; broken rows checked {bchk}, with an element {bhit}")


def fixed_objects():
    """the verifier against fixed objects: three false ones, each failing one duty, and the SWAP
    certificate of F138's coincident-magnitude family, true but outside the building rule"""
    F = Fraction
    bad = [([('Z', 1)], (('X',),), {'Y': F(1)}, "[H, G] = 0 fails (H = Z, A = X, G = Y)"),
           ([('Z', 1)], (('X',),), {'I': F(1)}, "{A, G} = 0 fails (H = Z, A = X, G = I)"),
           ([('IZ', 1)], (('Z',), ()), {'XI': F(1), 'XZ': F(1)}, "invertibility fails (H = IZ, A = ZI, G = XI + XZ)")]
    for terms, deph, elem, what in bad:
        check(f"the verifier rejects a false object: {what}", not verify(terms, deph, elem))
    H = [(t, 100) for t in ('XXI', 'YYI', 'ZZI', 'IXX', 'IYY', 'IZZ')] + [('XII', 30), ('IIX', -30)]
    deph = ((), ('X',), ())
    U = {'ZZZ': F(1, 2), 'YZY': F(-1, 2), 'XZX': F(-1, 2), 'IZI': F(1, 2)}
    check("the SWAP certificate (SWAP_02 Z Z Z, F138's coincident-magnitude family) is verified, and the building "
          "rule finds nothing there: a palindrome outside its grammar", verify(H, deph, U) and element(H, deph) is None)


if __name__ == "__main__":
    import json
    import platform
    import subprocess
    workers = 22
    import hashlib
    import sympy
    here = Path(__file__).resolve()
    try:
        r1 = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=here.parent)
        r2 = subprocess.run(['git', 'status', '--porcelain', '--', here.name], capture_output=True, text=True,
                            cwd=here.parent)
        rev = r1.stdout.strip() if r1.returncode == 0 else 'unknown'
        dirty = r2.returncode != 0 or bool(r2.stdout.strip())
    except OSError:
        rev, dirty = 'unknown', True
    sha = hashlib.sha256(here.read_bytes()).hexdigest()
    print(f"run: {time.strftime('%Y-%m-%d %H:%M:%S')}, revision {rev}{' (this script modified)' if dirty else ''}, "
          f"script sha256 {sha}, python {platform.python_version()}, numpy {np.__version__}, sympy {sympy.__version__}, "
          f"{workers} workers")
    print()
    print("## Fixed objects")
    fixed_objects()
    certs = []
    print()
    print("## Stage A: N = 3, the colouring page's census")
    t3 = census(3, workers, certs)
    exp3 = {('P3', 'Z'): {'sum': 52, 'conditioned': 10}, ('K3', 'Z'): {'conditioned': 18},
            ('bond+iso', 'Z'): {'product': 52}, ('P3', 'XY'): {'sum': 22}, ('K3', 'XY'): {},
            ('bond+iso', 'XY'): {'product': 104}, ('P3', 'XYZ'): {}, ('K3', 'XYZ'): {}, ('bond+iso', 'XYZ'): {}}
    pal3 = {('P3', 'Z'): 1385, ('P3', 'XY'): 625, ('P3', 'XYZ'): 603, ('K3', 'Z'): 1341, ('K3', 'XY'): 603,
            ('K3', 'XYZ'): 603, ('bond+iso', 'Z'): 1795, ('bond+iso', 'XY'): 1379, ('bond+iso', 'XYZ'): 1275}
    summary_checks(t3, exp3, "N = 3", 63 * 64, pal3, True)
    print()
    print("## Stage B + C: N = 4, chain, star, ring, complete graph")
    t4 = census(4, workers, certs)
    cf = control_fields(4)
    check("N = 4: the broken-row control patterns carry every letter on every site, every pair of letters on every "
          "pair of sites, and some fields on three and on four sites",
          all(any(p[a] == P and p[b] == Q for p in cf) for a in range(4) for b in range(4) if a != b
              for P in 'XYZ' for Q in 'XYZ') and any(sum(ch != 'I' for ch in p) == 3 for p in cf)
          and any('I' not in p for p in cf), f"{len(cf)} of 256 patterns")
    exp4 = {('P4', 'Z'): {'sum': 516, 'conditioned': 160}, ('S4', 'Z'): {'sum': 618, 'conditioned': 190},
            ('C4', 'Z'): {'conditioned': 240}, ('K4', 'Z'): {'conditioned': 56},
            ('P4', 'XY'): {'sum': 40}, ('S4', 'XY'): {}, ('C4', 'XY'): {}, ('K4', 'XY'): {},
            ('P4', 'XYZ'): {}, ('S4', 'XYZ'): {}, ('C4', 'XYZ'): {}, ('K4', 'XYZ'): {}}
    pal4 = {('P4', 'Z'): 13891, ('P4', 'XY'): 3835, ('P4', 'XYZ'): 3795, ('S4', 'Z'): 14023, ('S4', 'XY'): 3795,
            ('S4', 'XYZ'): 3795, ('C4', 'Z'): 13455, ('C4', 'XY'): 3795, ('C4', 'XYZ'): 3795, ('K4', 'Z'): 13271,
            ('K4', 'XY'): 3795, ('K4', 'XYZ'): 3795}
    summary_checks(t4, exp4, "N = 4", 255 * 256, pal4, False)
    col = {k: c['coloured'] for k, c in t4.items()}
    check("N = 4: the coloured counts do not depend on the graph (13215 at ZZ, 3795 at XX + YY and at XX + YY + ZZ)",
          all(v == {'Z': 13215, 'XY': 3795, 'XYZ': 3795}[k[1]] for k, v in col.items()), str(sorted(set(col.values()))))
    rows = sum(c['palindromic'] + c['broken'] for c in t4.values())
    check("N = 4: every row counted (4 graphs x 3 bond sets x 255 dephasing patterns x 256 field patterns)",
          rows == 4 * 3 * 255 * 256, str(rows))
    out = Path(__file__).parent / 'results' / 'anticommuting_sum_certificates.json'
    out.write_bytes(json.dumps(sorted(certs, key=lambda c: (c['n'], c['graph'], c['bonds'], c['dephasing'],
                                                             c['fields'])), indent=0).encode('utf-8'))
    check("every element beyond the colouring written out with its row (n, graph, bonds, dephasing, fields, kind, "
          "rational coefficients)", len(certs) == 258 + 1820, f"{len(certs)} certificates -> {out.name}")
    print()
    if FAIL:
        print(f"{len(FAIL)} FAILURE(S): {FAIL}")
        sys.exit(1)
    print("ALL STAGES PASS")
