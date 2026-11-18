"""Gate: the block lemma for far elements on a sheet (PROOF_ONE_DEPHASED_SEAT, section on blind dephased seats,
bullet "The block from one excitation to N - 1" under the paragraph "The sheet and blindness, at N = 3 and N = 4").

Setting: Heisenberg bonds sum J_ij sigma_i.sigma_j on a connected graph (integer couplings, signs free), the jump Z on
seat d = 0, planar integer fields on the sites != d, none on d. A far element W ([H, W] = 0, {Z_d, W} = 0) homogeneous
under Z^(x)N commutes with the bonds and with the fields separately. Four parts, exact or certified, no threshold.

1. The block lemma, N = 3, 4, 5, no fields involved. The blocks T from the one-excitation sector to the
   (N - 1)-excitation sector with B_(N-1)·T = T·B_1 and Z_d·T = -T·Z_d (B the bond part, restricted from the full
   2^N x 2^N matrix built from Pauli matrices, not from the Laplacian identity) are solved exactly (sympy nullspace).
   Predicted: the solution space is one-dimensional, spanned by the block of Xi = X^(x)N, exactly when seat d is not
   blind (blind read independently, as N minus the rank of d's Krylov matrix of the signed Laplacian); at blind rows
   it is larger, 1 + dim commutant of the Laplacian on the Krylov complement, which the gate also computes.
2. N = 3, even height, bonds and jump alone. The far elements of even height of the bonds alone (no field: the
   field equations play no part in the hand step) are counted exactly (framework f158_end_count._kernel on the
   strings with X or Y at d and an even number of letters X or Y). Predicted: zero at every row where d is not
   blind, on distinct rows. Read beside it, not predicted by the step: the counts at blind rows of the swap kind
   (J_du = J_dv, e2 != 0) and of the reciprocal kind (e2 = 0), with and without fields, the nonzero ones exactly.
3. N = 3, odd height, at a seat that is not blind with the fields on one line (class A). The step's two outputs,
   T = tau·Xi on the one-excitation sector and tau·omega_j = chi·conj(omega_j) (omega_j = f_x + i f_y,
   chi = <top|W|vac>; the proof writes Xi and chi because Phi and c are spent in its section), are read off the exact far element of odd height: non-vacuous checks of the block lemma and
   of the vacuum equation, at rows where they have a solution. Controls: the fields not on one line at the same
   couplings give no far element of odd height (the step's conclusion).
4. Measured, N = 4 and N = 5: random connected rows (each edge of K_N present with probability 0.6, couplings from
   {-3, ..., 3} \\ {0}), fields drawn from a list of integer planar vectors, (0, 0) included. In the cell "d not blind,
   fields not on one line" the number of rows with a far element homogeneous under Z^(x)N, in either parity class, is
   counted on distinct rows by a rank modulo the prime p = 10^9 + 9 (p = 1 mod 4, i mapped to a square root of -1): a kernel modulo p
   is at least the rational kernel, so a zero there is an exact zero. Predicted: zero. Not proven at N >= 4. Control:
   the drawn rows with d not blind and the fields on one line all carry one (the colouring of class A).

Run: python simulations/blind_seat_block_lemma_gate.py   (about half a minute; writes results/blind_seat_block_lemma_gate.txt)
"""
import itertools
import random
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import pauli_product
from framework.diagnostics.f158_end_count import _kernel

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")   # the console's default code page cannot print every character below

OUT, FAILS = [], []


def say(s=""):
    print(s)
    OUT.append(s)


def check(cond, what):
    say(("  PASS  " if cond else "  FAIL  ") + what)
    if not cond:
        FAILS.append(what)


# ---- strings and the exact kernels -----------------------------------------------------------------------------------
P = 1000000009
I_P = next(pow(g, (P - 1) // 4, P) for g in range(2, 100) if pow(g, (P - 1) // 2, P) == P - 1)
assert I_P * I_P % P == P - 1
IPOW = [1, I_P, P - 1, P - I_P]


def anticomm(s, t):
    return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1


def height(s):
    return sum(1 for a in s if a in 'XY')


def h_terms(n, bonds, fields):
    t = {}
    for (i, j), J in bonds.items():
        for L in 'XYZ':
            s = ['I'] * n; s[i] = L; s[j] = L
            t[''.join(s)] = (J, 0)
    for j, (fx, fy) in fields.items():
        for L, c in (('X', fx), ('Y', fy)):
            if c:
                s = ['I'] * n; s[j] = L
                t[''.join(s)] = (c, 0)
    return t


def lit_span(n, parity, d=0):
    return [''.join(w) for w in itertools.product('IXYZ', repeat=n) if w[d] in 'XY' and height(w) % 2 == parity]


def rank_mod_p(M):
    M = M.copy() % P
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), P - 2, P)) % P
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
        if r == rows:
            break
    return r


def kernel_dim_mod_p(terms, span):
    rowidx, entries = {}, {}
    for col, s in enumerate(span):
        for Pp, (re, im) in terms.items():
            if not anticomm(Pp, s):
                continue
            e, u = pauli_product(Pp, s)
            r = rowidx.setdefault(u, len(rowidx))
            entries[(r, col)] = (entries.get((r, col), 0) + 2 * (re + im * I_P) * IPOW[e]) % P
    if not rowidx:
        return len(span)
    M = np.zeros((len(rowidx), len(span)), dtype=np.int64)
    for (r, c), v in entries.items():
        M[r, c] = v
    return len(span) - rank_mod_p(M)


def exact_far(n, bonds, fields, parity):
    ht = {s: sp.Integer(re) + sp.I * im for s, (re, im) in h_terms(n, bonds, fields).items()}
    span = lit_span(n, parity)
    return span, _kernel(ht, span)


# ---- graphs, blindness -----------------------------------------------------------------------------------------------
def laplacian(n, bonds):
    L = sp.zeros(n, n)
    for (i, j), J in bonds.items():
        L[i, j] -= J; L[j, i] -= J; L[i, i] += J; L[j, j] += J
    return L


def krylov(n, bonds, d=0):
    L = laplacian(n, bonds)
    v = sp.zeros(n, 1); v[d] = 1
    K = [v]
    for _ in range(n - 1):
        K.append(L * K[-1])
    return L, sp.Matrix.hstack(*K)


def blind_dim(n, bonds, d=0):
    return n - krylov(n, bonds, d)[1].rank()


def connected(n, bonds):
    adj = {i: set() for i in range(n)}
    for (i, j), J in bonds.items():
        if J:
            adj[i].add(j); adj[j].add(i)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x] - seen:
            seen.add(y); st.append(y)
    return len(seen) == n


def on_one_line(fields):
    vs = [v for v in fields.values() if v != (0, 0)]
    return all(a[0] * b[1] - a[1] * b[0] == 0 for a in vs for b in vs)


def random_row(n, rng, p_edge=0.6, cpl=(-3, -2, -1, 1, 2, 3)):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    while True:
        b = {e: rng.choice(cpl) for e in edges if rng.random() < p_edge}
        if connected(n, b):
            return b


def swap_symmetric_row(n, rng):
    """a row symmetric under the swap of sites 1 and 2 (which fixes d = 0): d blind by the antisymmetric mode"""
    while True:
        b = random_row(n, rng)
        sym = {}
        for (i, j), J in b.items():
            for (a, c) in ((i, j), tuple(sorted(({1: 2, 2: 1}.get(i, i), {1: 2, 2: 1}.get(j, j))))):
                sym[(a, c)] = J
        if connected(n, sym):
            return sym


# ---- part 1: the block lemma -----------------------------------------------------------------------------------------
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]])
PZ = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(n, j, M):
    out = np.array([[1]], dtype=complex)
    for k in range(n):
        out = np.kron(out, M if k == j else np.eye(2))
    return out


def bond_matrix(n, bonds):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for (i, j), J in bonds.items():
        for M in (PX, PY, PZ):
            H += J * site_op(n, i, M) @ site_op(n, j, M)
    assert np.all(H.imag == 0) and np.all(H.real == np.round(H.real))
    return sp.Matrix(np.round(H.real).astype(int).tolist())


def sector(n, k):
    # site 0 is the leftmost tensor factor, so bit (n - 1 - j) of a basis index is site j; |1> is the excitation
    return [x for x in range(2 ** n) if bin(x).count('1') == k]


def block_solutions(n, bonds):
    H = bond_matrix(n, bonds)
    Zd = sp.diag(*[(-1 if (x >> (n - 1)) & 1 else 1) for x in range(2 ** n)])
    s1, sN = sector(n, 1), sector(n, n - 1)
    B1, BN = H.extract(s1, s1), H.extract(sN, sN)
    Z1, ZN = Zd.extract(s1, s1), Zd.extract(sN, sN)
    m = len(s1)
    T = sp.Matrix(m, m, sp.symbols(f't0:{m * m}'))
    eqs = list(BN * T - T * B1) + list(ZN * T + T * Z1)
    A, _ = sp.linear_eq_to_matrix(eqs, list(T))
    ns = A.nullspace()
    # Xi = X^(x)N maps |x> to |x xor 11..1>; its block from sector 1 to sector N-1
    full = (1 << n) - 1
    Xi = sp.Matrix(m, m, lambda r, c: 1 if sN[r] == (s1[c] ^ full) else 0)
    phi_ok = (BN * Xi - Xi * B1) == sp.zeros(m, m) and (ZN * Xi + Xi * Z1) == sp.zeros(m, m)
    return len(ns), phi_ok, Xi, ns


def commutant_on_complement(n, bonds, d=0):
    L, K = krylov(n, bonds, d)
    r = K.rank()
    if r == n:
        return 0
    comp = sp.Matrix.hstack(*K.T.nullspace())          # orthogonal complement of the Krylov space (L symmetric)
    Lc = (comp.T * comp).inv() * comp.T * L * comp       # L restricted to that invariant complement
    k = Lc.shape[0]
    X = sp.Matrix(k, k, sp.symbols(f'c0:{k * k}'))
    A, _ = sp.linear_eq_to_matrix(list(Lc * X - X * Lc), list(X))
    return len(A.nullspace())


say("1. the block lemma: blocks T from one excitation to N - 1 with [B, T] = 0 and {Z_d, T} = 0")
rng = random.Random(20261118)
for n, count in ((3, 40), (4, 40), (5, 25)):
    stats = {'nonblind': 0, 'nonblind_ok': 0, 'blind': 0, 'blind_ok': 0, 'phi_ok': 0}
    rows = [random_row(n, rng) for _ in range(count)] + [swap_symmetric_row(n, rng) for _ in range(count // 2)]
    if n == 3:
        rows += [{(0, 1): 3, (0, 2): 6, (1, 2): -2}, {(0, 1): 2, (0, 2): 2, (1, 2): -1}, {(0, 1): 1, (0, 2): -2, (1, 2): 2}]
    for b in rows:
        dim, phi_ok, _, _ = block_solutions(n, b)
        bl = blind_dim(n, b)
        stats['phi_ok'] += phi_ok
        if bl == 0:
            stats['nonblind'] += 1
            stats['nonblind_ok'] += (dim == 1)
        else:
            stats['blind'] += 1
            stats['blind_ok'] += (dim == 1 + commutant_on_complement(n, b) and dim >= 2)
    check(stats['phi_ok'] == len(rows), f"N = {n}: Xi's block solves both conditions at all {len(rows)} rows")
    check(stats['nonblind'] > 0 and stats['nonblind_ok'] == stats['nonblind'],
          f"N = {n}: at the {stats['nonblind']} rows where d is not blind the solution space is one-dimensional (Xi's block)")
    check(stats['blind'] > 0 and stats['blind_ok'] == stats['blind'],
          f"N = {n}: at the {stats['blind']} blind rows it is 1 + dim commutant on the Krylov complement, at least two")

# ---- part 2: even height at N = 3, bonds and jump alone --------------------------------------------------------------
say()
say("2. N = 3, even height, bonds and jump alone (no field)")
rng = random.Random(7)
seen, counts = set(), []
for _ in range(400):
    b = random_row(3, rng, p_edge=0.85, cpl=(-3, -2, -1, 1, 2, 3, 6))
    key = tuple(sorted(b.items()))
    if key in seen or blind_dim(3, b) > 0:
        continue
    seen.add(key)
    counts.append(kernel_dim_mod_p(h_terms(3, b, {}), lit_span(3, 0)))
check(len(counts) > 150 and all(c == 0 for c in counts),
      f"at {len(counts)} distinct rows where d is not blind: no far element of even height (rank mod p, an exact zero)")
grid = range(-4, 5)
swap_rows, recip_rows = [], []
for s_, t_, r_ in itertools.product(grid, grid, grid):
    if [s_, t_, r_].count(0) > 1:
        continue
    e2 = s_ * t_ + s_ * r_ + t_ * r_
    b = {e: J for e, J in (((0, 1), s_), ((0, 2), t_), ((1, 2), r_)) if J}
    if s_ == t_ and e2 != 0:
        swap_rows.append(b)
    elif e2 == 0 and s_ != t_:
        recip_rows.append(b)
FLD = {1: (3, 4), 2: (5, 0)}
sw0 = [kernel_dim_mod_p(h_terms(3, b, {}), lit_span(3, 0)) for b in swap_rows]
swF = [kernel_dim_mod_p(h_terms(3, b, FLD), lit_span(3, 0)) for b in swap_rows]
re0 = [len(exact_far(3, b, {}, 0)[1]) for b in recip_rows]                  # nonzero: read exactly
reF = [kernel_dim_mod_p(h_terms(3, b, FLD), lit_span(3, 0)) for b in recip_rows]
say(f"  read: couplings in {{-4..4}}, at most one zero; swap kind J_du = J_dv, e2 != 0: {len(swap_rows)} rows, "
    f"even-height counts without fields {sorted(set(sw0))}, with f_u = (3, 4), f_v = (5, 0) {sorted(set(swF))}")
say(f"  read: reciprocal kind e2 = 0, J_du != J_dv: {len(recip_rows)} rows, "
    f"even-height counts without fields (exact) {sorted(set(re0))}, with those fields {sorted(set(reF))}")
check(len(swap_rows) == 68 and set(sw0) == {0} and set(swF) == {0} and len(recip_rows) == 8 and set(re0) == {4}
      and set(reF) == {0}, "the readings beside the proof, pinned as the document quotes them: 68 rows 0 and 0, 8 rows 4 and 0")

# ---- part 3: odd height at N = 3, class A ----------------------------------------------------------------------------
say()
say("3. N = 3, odd height, d not blind, fields on one line (class A): T = tau·Xi and tau·omega_j = chi·conj(omega_j), Xi = X^(x)3")
PAULI = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
         'Z': sp.Matrix([[1, 0], [0, -1]])}


def string_matrix(s):
    out = sp.Matrix([[1]])
    for a in s:
        out = sp.kronecker_product(out, PAULI[a])
    return out


s1, s2 = sector(3, 1), sector(3, 2)
full = 7
Xi3 = sp.Matrix(3, 3, lambda r, c: 1 if s2[r] == (s1[c] ^ full) else 0)
A_rows = [({(0, 1): 1, (0, 2): 2, (1, 2): 3}, {1: (3, 4), 2: (-3, -4)}),
          ({(0, 1): 2, (0, 2): -1, (1, 2): 1}, {1: (1, 1), 2: (2, 2)}),
          ({(0, 1): 3, (1, 2): -2}, {1: (0, 1), 2: (0, 0)}),
          ({(0, 1): 1, (0, 2): 3}, {1: (4, 3), 2: (4, 3)})]
for b, f in A_rows:
    assert blind_dim(3, b) == 0 and on_one_line(f)
    span, ker = exact_far(3, b, f, 1)
    ok, nonvacuous = len(ker) >= 1, False
    for vec in ker:
        W = sum((c * string_matrix(s) for s, c in vec.items()), sp.zeros(8))
        T = W.extract(s2, s1)
        i0, j0 = next((i, j) for i in range(3) for j in range(3) if Xi3[i, j])
        tau = T[i0, j0]
        chi = W[7, 0]                                        # <top| W |vac>, top = |111> (index 7), vac = |000>
        om = {j: sp.Integer(fx) + sp.I * fy for j, (fx, fy) in f.items()}
        ok = ok and sp.simplify(T - tau * Xi3) == sp.zeros(3, 3)
        ok = ok and all(sp.simplify(tau * om[j] - chi * sp.conjugate(om[j])) == 0 for j in om)
        nonvacuous = nonvacuous or tau != 0
    check(ok and nonvacuous, f"J = {b}, fields {f}: {len(ker)} far element(s) of odd height, each with T = tau·Xi and tau·omega = chi·conj(omega), tau != 0 for one")
for b, _ in A_rows[:2]:
    f = {1: (3, 4), 2: (5, 0)}
    check(len(exact_far(3, b, f, 1)[1]) == 0, f"control, J = {b} with fields {f} not on one line: no far element of odd height")

# ---- part 4: measured at N = 4 and N = 5 -----------------------------------------------------------------------------
say()
say("4. measured, N = 4 and N = 5: d not blind, fields not on one line, far elements homogeneous under Z^(x)N")
FIELDS = [(1, 0), (0, 1), (-1, 0), (0, -1), (3, 4), (4, 3), (-3, 4), (5, 0), (0, 5), (4, -3), (1, 1), (-1, 1), (2, 0),
          (0, 2), (0, 0)]
for n, count, seed in ((4, 300, 41), (5, 120, 51)):
    rng = random.Random(seed)
    spans = (lit_span(n, 0), lit_span(n, 1))
    cell, hits, line_rows, line_hits, seen = 0, 0, 0, 0, set()
    for _ in range(count):
        b = random_row(n, rng)
        f = {j: rng.choice(FIELDS) for j in range(1, n)}
        key = (tuple(sorted(b.items())), tuple(sorted(f.items())))
        if key in seen or blind_dim(n, b) > 0:
            continue
        seen.add(key)
        found = any(kernel_dim_mod_p(h_terms(n, b, f), sp_) for sp_ in spans)
        if on_one_line(f):
            line_rows += 1
            line_hits += found
        else:
            cell += 1
            hits += found
    check(cell > 0 and hits == 0,
          f"N = {n}: {cell} distinct rows in the cell (of {count} drawn), none with a far element in either parity class")
    check(line_rows > 0 and line_hits == line_rows,
          f"N = {n}: control, the {line_rows} distinct rows with d not blind and the fields on one line all carry one (class A's colouring)")

say()
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED")
Path(__file__).parent.joinpath("results", "blind_seat_block_lemma_gate.txt").write_bytes(("\n".join(OUT) + "\n").encode("utf-8"))
raise SystemExit(1 if FAILS else 0)
