"""Gate for old open questions answered on the letter cube (exact unless a line says "reading").

G5  The Marrakesh f83 fingerprint cells. Z dephasing is diagonal in the Pauli basis and keeps every string; a
    commutator -i[h, P] with an H term adds h's Klein cell (n_XY parity, w_YZ parity) and turns the y-parity (the
    parity of the number of Y letters, real versus imaginary) by y(h) + 1. So a string is reachable only at
    (cell(P) + sum_j n_j cell(h_j), y(P) + sum_j n_j (y(h_j) + 1)), P in the start's support. Exact Pauli-string
    propagation of L^n(rho0) with Fraction coefficients, run until the Krylov span stops growing (so every order is
    covered), on N = 3 with the four f83 Hamiltonians and product starts in the X, Y and Z bases: every string
    reached obeys the rule, and the two-site observables P0 I P2 reached fall in exactly the literal cell sets below
    (X basis: (0,0) and the "X-axis flip"). The pair rule itself is checked on every anticommuting pair of
    three-site strings.
G6  The kernel of L (Z dephasing on every site) split by Ad_{X^N}. Theorem 2 of the complement connection: the
    kernel is spanned by the component projectors of the hopping graph; X^N pairs components, so even : odd =
    (p + s) : p. The gate builds P_C + P_Cbar and P_C - P_Cbar and checks L v = 0 exactly (the lower bound), and
    the nullity of L on each Ad_X parity mod a prime 1 mod 4 (an upper bound); equality certifies the split.
G7  F137 at H != 0, pure amplitude damping, any site rates: the ingredients of the proof, as exact matrix
    identities on XXZ chains and a ring (dyadic rates and couplings, so floating point is exact here): L is block-triangular in
    the joint popcount, its diagonal blocks are -i(H_eff (x) I - I (x) H_eff^*) with H_eff = H - (i/2) sum g_l n_l, and
    X^N H_eff X^N = H_eff^dagger - (i/2) Gamma. The palindrome about -Gamma/2 is printed as a reading, and so is
    the pairing of two Hamiltonians outside the conditions (XX bonds alone, DM bonds XY - YX). For DM the weaker
    X^N step is exact: X^N H X^N = conj(H) and X^N H_eff X^N = conj(H_eff) - (i/2) Gamma.
G8  F49's cross term under amplitude damping: D[sigma^-] = (g/4)(D_X + D_Y) + M, M the move I -> Z;
    {L_H, lights} and {L_H, M} are orthogonal and the move's part is 32 (N-1)^2, 32 (N-1)(2N-3), 32 (N^2-N-1),
    32 (3N^2-6N+2) times g^2 4^(N-2) for XX, XY, ZZ, Heisenberg on the chain; on ring, star and complete graphs
    the XX and XX+YY parts are bond-additive, |E| times 32 (N-1) and 32 (2N-3); the Ising part is
    32 [(N-1)|E| + sum_v C(deg v, 2)], each pair of bonds sharing a site adding 32 (times g^2 4^(N-2)), and Heisenberg the sum.
    Entries are dyadic, so the norms are exact.
G9  n_XY parity (the character of Ad_{Z^N}, (-1)^(p - q) on |a><b|) under amplitude damping: sigma^- rho sigma^+ lowers
    bra and ket together, so p - q stays a symmetry and the parity is its shadow; the blocks of L (components of its
    sparsity graph, an exact count) are the 2N + 1 values of p - q for XXZ, the two parities for XX bonds alone (an H
    that changes the excitation number by two), and one block once a transverse X field is on.
G10 The V-Effect census's 36 two-term chain pairs at every N >= 3 (on one bond all 36 pair). The 22 palindromic pairs
    have a product mirror, X, Y or X-Y per site with period two, commuting with H and anticommuting with every Z_l (F158
    section (e)), the same at N = 3 and 4 (the bond check holds both orientations from N = 3); the 14 hard ones have none.
    Each of the 14 has an F158 odd word at the chain's end (powers of H between jumps Z_l, an odd number of jumps) with a
    nonzero trace, exact and constant from N = 3 or 4 to N = 2P + 6; far bonds enter only as clusters of at least two
    factors that commute with the rest, so from N = 2P + 4 the trace is a polynomial in N of degree at most P // 2 <= 2,
    and three equal values make it constant. The F158 end count at N = 3 is a second route for all 14, with XX+XZ
    as the palindromic control. Swept for G10: experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md (the 36 at N = 3, 4, 5),
    OPERATOR_RIGIDITY_ACROSS_CUSP.md, NON_HEISENBERG_PALINDROME.md (the 22 proved bond by bond), SOFTNESS_IS_N_DEPENDENT.md (k >= 3 and multi-term, not these), V_EFFECT_PALINDROME.md
    (all 36 pair at N = 2), F158 (f5) and PalindromeSoftCertifier; none held the hard side at every N.
G11 F108's orbit under local unitaries. A setting is a pair (P, Q) of distinct letters: jumps P on every site and an H of
    two-site strings commuting with Q^N (Part 1 (Z, X), Part 2 (X, Z), Part 3 (Y, X)). Over all 24^3 site-dependent
    local Cliffords at N = 3, Part 1 reaches all six settings, 64 ways each, always one letter permutation on every site,
    and its stabilizer is the local Pauli group; the quarter turn about X carries Part 1 to Part 3; the three settings
    no Part names have equal F158 end counts (a field along the jump, the control, makes them differ). Swept for G11:
    the F108 Klein-V4 proof, the F108 Parts' Open sections and typed claims, PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4 (its
    quarter turn about X carries Z-dephasing to Y-dephasing, setting (Z, X) to (Y, X)), PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE, THE_THREE_DIAGONALS
    (the dissipator's Clifford orbit), CAUGHT_ERRORS (vi) of 2026-09-05, F103 section 8 (the quarter turn between the two
    flip letters), NON_HEISENBERG_PALINDROME (the P4 family is the setting (Z, Y)), MirrorGroup's letter S3. New: the
    subgroup, its stabilizer, and the exclusion of continuous rotations (argued in the proof, not gated here).
Sweep before writing: docs/ANALYTICAL_FORMULAS.md (F4, F49e, F88a, F88b, F137, F155), docs/proofs/ (the complement
connection's Theorem 2, PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS, MIRROR_SYMMETRY_PROOF's Scope, PROOF_CROSS_TERM_FORMULA,
PROOF_F155), experiments/THERMAL_BREAKING.md, data/ibm_f83_signature_april2026, fw.Confirmations (the f83 entry),
the OpenArcs registry (nothing on these questions), the Pi2 open questions in Core (two of them answered here),
NoJumpGenerator and the C# test T1BreakingInformationalTests (the triangularity). The deleted Pi2 open question,
F88b and the complement connection's own sweep (F4, XOR_SPACE) held the split for the popcount case; the general
component statement (p + s) : p is new here. None held the mechanism, the every-N proof or the move's values.
Adjacency: the triangularity sits beside F155's no-jump generator and the T1 test, now cited from F137.
Run: python simulations/cube_old_questions_gate.py"""
import itertools, sys
from fractions import Fraction as F
import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

FAILS = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

# ---------------------------------------------------------------- Pauli algebra (exact)
MUL = {("I", c): (1, c) for c in "IXYZ"}
MUL.update({(c, "I"): (1, c) for c in "IXYZ"})
MUL.update({(c, c): (1, "I") for c in "XYZ"})
MUL.update({("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"), ("Y", "Z"): (1j, "X"),
            ("Z", "Y"): (-1j, "X"), ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y")})
def smul(a, b):
    ph, out = 1, []
    for x, y in zip(a, b):
        p, c = MUL[(x, y)]; ph *= p; out.append(c)
    return ph, "".join(out)
def anti(a, b):
    return sum(1 for x, y in zip(a, b) if x != "I" and y != "I" and x != y) % 2 == 1
def cell(w):
    return (sum(c in "XY" for c in w) % 2, sum(c in "YZ" for c in w) % 2)
def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def cdiv(a, b):
    d = b[0] * b[0] + b[1] * b[1]; return ((a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d)
def basis_name(x): return x
def cadd(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2)

# ---------------------------------------------------------------- G5
N3 = 3
def pair_rule_violations(n):
    words = ["".join(t) for t in itertools.product("IXYZ", repeat=n)]
    bad = count = 0
    for h in words:
        for P in words:
            if not anti(h, P): continue
            count += 1; ph, W = smul(h, P); c = -2j * ph          # -i[h, P] = c W, c must be real
            yh, yP, yW = (w.count("Y") % 2 for w in (h, P, W))
            if c.imag != 0 or cell(W) != cadd(cell(h), cell(P)) or yW != (yh + yP + 1) % 2: bad += 1
    return count, bad
_cnt, _bad = pair_rule_violations(N3)
check(f"G5 pair rule on all {_cnt} anticommuting pairs of 3-site strings: cell adds, y turns by y(h)+1",
      _cnt == 2016 and _bad == 0, f"{_bad} violations")
CATS = [("truly", [("X", "X"), ("Y", "Y")]), ("pi2_odd_pure", [("X", "Y"), ("Y", "X")]),
        ("pi2_even_nontruly", [("Y", "Z"), ("Z", "Y")]), ("mixed", [("X", "Y"), ("Y", "Z")])]
EXPECT = {  # two-site outer observables reached, by start basis
    "X": {"truly": {(0, 0), (1, 0)}, "pi2_odd_pure": {(0, 0), (1, 1)}, "pi2_even_nontruly": {(0, 0)},
          "mixed": {(0, 0), (1, 1)}},
    "Y": {"truly": {(0, 0), (1, 1)}, "pi2_odd_pure": {(0, 0), (1, 0)}, "pi2_even_nontruly": {(0, 0), (0, 1)},
          "mixed": {(0, 0), (0, 1), (1, 0), (1, 1)}},
    "Z": {"truly": {(0, 0)}, "pi2_odd_pure": {(0, 0)}, "pi2_even_nontruly": {(0, 0), (1, 1)},
          "mixed": {(0, 0), (1, 1)}},
}
GAM = [F(1, 2), F(3, 4), F(1, 4)]
def L_apply(rho, Hterms):
    out = {}
    for P, c in rho.items():
        for h in Hterms:
            if anti(h, P):
                ph, w = smul(h, P)            # -i [h, P] = -2i hP
                re, im = F(round((-2j * ph).real)) , F(round((-2j * ph).imag))
                a, b = c                      # (a + ib)(re + i im)
                nr, ni = a * re - b * im, a * im + b * re
                o = out.get(w, (F(0), F(0))); out[w] = (o[0] + nr, o[1] + ni)
        d = sum(2 * g for l, g in enumerate(GAM) if P[l] in "XY")
        if d:
            o = out.get(P, (F(0), F(0))); out[P] = (o[0] - d * c[0], o[1] - d * c[1])
    return {w: c for w, c in out.items() if c != (F(0), F(0))}
for basis_letter in "XYZ":
    basis = basis_letter
    signs = (1, -1, 1)
    rho0 = {}
    for choice in itertools.product((0, 1), repeat=N3):
        w = "".join(basis if ch else "I" for ch in choice)
        s = 1
        for l, ch in enumerate(choice):
            if ch: s *= signs[l]
        rho0[w] = (F(s, 8), F(0))
    c0 = {cell(w) for w in rho0}
    for name, terms in CATS:
        Hterms = [("".join(a if i == b else (bb if i == b + 1 else "I") for i in range(N3))) for b in range(N3 - 1) for a, bb in terms]
        gen = {cell(h) for h in Hterms}
        group = {(0, 0)}
        for _ in range(3):
            group |= {cadd(x, y) for x in group for y in gen}
        allowed = {cadd(x, y) for x in c0 for y in group}
        ypar = lambda w: w.count("Y") % 2
        moves = {(cell(h), (ypar(h) + 1) % 2) for h in Hterms}
        states = {(cell(w), ypar(w)) for w in rho0}
        for _ in range(4):
            states |= {(cadd(c, m[0]), (y + m[1]) % 2) for (c, y) in states for m in moves}
        reached = dict(rho0); cur = rho0
        kbasis = []                                  # exact Krylov span, echelon form over Q(i)
        def reduce(vec):
            v = dict(vec)
            for piv, b in kbasis:
                if piv in v:
                    a = v[piv]; fac = cdiv(a, b[piv])
                    for w2, c2 in b.items():
                        cur2 = v.get(w2, (F(0), F(0))); d = cmul(fac, c2)
                        v[w2] = (cur2[0] - d[0], cur2[1] - d[1])
                    v = {w2: c2 for w2, c2 in v.items() if c2 != (F(0), F(0))}
            return v
        def add_to_span(vec):
            v = reduce(vec)
            if not v: return False
            kbasis.append((min(v), v)); return True
        add_to_span(rho0)
        saturated = False
        for n in range(64):
            cur = L_apply(cur, Hterms)
            for w in cur: reached[w] = True
            if not add_to_span(cur): saturated = True; break
        check(f"G5 {basis_name(basis_letter)}-basis start, {name}: the Krylov span saturates (order {n + 1})", saturated)
        inside = all(cell(w) in allowed and (cell(w), ypar(w)) in states for w in reached)
        outer = {cell(w) for w in reached if w[1] == "I" and w[0] != "I" and w[2] != "I"}
        check(f"G5 {basis_letter}-basis start, {name}: every reached string obeys the cell and y-parity rule", inside)
        check(f"G5 {basis_letter}-basis start, {name}: outer two-site cells reached = {sorted(EXPECT[basis_letter][name])}",
              outer == EXPECT[basis_letter][name], f"reached {sorted(outer)}")

# ---------------------------------------------------------------- dense helpers
PM = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], complex),
      "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]).astype(complex)}
SM = np.array([[0, 1], [0, 0]], complex)
def op(N, d):
    m = np.array([[1]], complex)
    for i in range(N): m = np.kron(m, d.get(i, np.eye(2)))
    return m
def lr(A, B): return np.kron(A, B.T)       # rho -> A rho B, row-major vec
def dissip(c):
    D = c.shape[0]; I = np.eye(D); cd = c.conj().T
    return lr(c, cd) - 0.5 * lr(cd @ c, I) - 0.5 * lr(I, cd @ c)
def bonds_chain(N): return [(i, i + 1) for i in range(N - 1)]
def ham(N, terms, bonds, coeffs=None):
    coeffs = coeffs or [1] * len(bonds)
    return sum(coeffs[k] * op(N, {a: PM[x], b: PM[y]}) for k, (a, b) in enumerate(bonds) for x, y in terms)

# ---------------------------------------------------------------- G6
PRIME = 998244353  # 1 mod 4, so -1 has a square root
def nullity_mod_p(M):
    A = M.copy() % PRIME; rows, cols = A.shape; r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i, c] != 0), None)
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), PRIME - 2, PRIME)
        A[r] = (A[r] * inv) % PRIME
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r: A[i] = (A[i] - A[i, c] * A[r]) % PRIME
        r += 1
        if r == rows: break
    return cols - r
def kernel_split_mod_p(H, N):
    """Nullity of L = -i[H, .] + sum_l (l+1)(Z_l rho Z_l - rho) on the Ad_X = +1 and -1 vectors, mod PRIME with i a
    square root of -1 there; H must have integer entries."""
    D = 2 ** N
    Hi = np.round(H.real).astype(np.int64); assert np.all(np.abs(H.imag) < 1e-12) and np.all(np.abs(H.real - Hi) < 1e-12)
    pop = [bin(x).count("1") for x in range(D)]
    full = D - 1
    pairs = {}
    for x in range(D):
        for y in range(D):
            key = min((x, y), (full ^ x, full ^ y))
            pairs.setdefault(key, None)
    reps = sorted(pairs)
    idx = {r: k for k, r in enumerate(reps)}
    out = []
    for sign in (1, -1):
        M = np.zeros((len(reps), len(reps)), dtype=np.int64)
        for (x, y), k in idx.items():
            # vector v = e_{xy} + sign e_{x̄ȳ}; i*L v = [H, v] + i*D v; D diagonal: -2 sum_l l [x_l != y_l]
            # i*L has the imaginary dephasing; L itself: -i[H,.] + D. Nullity of L = nullity of (L) over Q(i);
            # write L = -i (C) + Dg with C=[H,.] integer, Dg integer diagonal; use i -> sqrt(-1) mod p.
            for (a, b), s in (((x, y), 1), ((full ^ x, full ^ y), sign)):
                # column of L at e_{ab}: -i (sum_c H[c,a] e_{cb} - sum_c H[b,c] e_{ac}) + Dg[a,b] e_{ab}
                for c_ in range(D):
                    if Hi[c_, a]:
                        t = (c_, b); rep = min(t, (full ^ c_, full ^ b)); sg = 1 if rep == t else sign
                        M[idx[rep], k] += -SQRTM1 * s * sg * Hi[c_, a]
                    if Hi[b, c_]:
                        t = (a, c_); rep = min(t, (full ^ a, full ^ c_)); sg = 1 if rep == t else sign
                        M[idx[rep], k] += SQRTM1 * s * sg * Hi[b, c_]
                dg = -2 * sum((l + 1) for l in range(N) if ((a ^ b) >> (N - 1 - l)) & 1)
                rep = min((a, b), (full ^ a, full ^ b)); sg = 1 if rep == (a, b) else sign
                M[idx[rep], k] += s * sg * dg
        out.append(nullity_mod_p(M % PRIME))
    return tuple(out)
SQRTM1 = next(r for r in (pow(a, (PRIME - 1) // 4, PRIME) for a in range(2, 100)) if r * r % PRIME == PRIME - 1)
def graph_split(H, N):
    D = 2 ** N; parent = list(range(D))
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for x in range(D):
        for y in range(D):
            if x != y and abs(H[x, y]) > 0: parent[f(x)] = f(y)
    comps = {}
    for x in range(D): comps.setdefault(f(x), set()).add(x)
    comps = [frozenset(c) for c in comps.values()]
    p = s = 0; seen = set()
    for c in comps:
        if c in seen: continue
        cb = frozenset((D - 1) ^ x for x in c)
        if cb == c: s += 1
        else: p += 1; seen.add(cb)
        seen.add(c)
    return p + s, p
G6CASES = []
for N in (3, 4, 5):
    G6CASES.append((f"XY chain N={N}", ham(N, [("X", "X"), ("Y", "Y")], bonds_chain(N)), N))
    G6CASES.append((f"Heisenberg chain N={N}", ham(N, [("X", "X"), ("Y", "Y"), ("Z", "Z")], bonds_chain(N)), N))
for N in (4, 5):
    four = sum(op(N, {i: PM["X"], i + 1: PM["X"], i + 2: PM["X"], i + 3: PM["X"]}) for i in range(N - 3))
    G6CASES.append((f"XY chain + XXXX N={N}", ham(N, [("X", "X"), ("Y", "Y")], bonds_chain(N)) + 2 * four, N))
    G6CASES.append((f"XXXX alone N={N}", four, N))
G6CASES.append(("two disjoint XY bonds N=4", ham(4, [("X", "X"), ("Y", "Y")], [(0, 1), (2, 3)]), 4))
G6CASES.append(("XY ring N=4", ham(4, [("X", "X"), ("Y", "Y")], [(0, 1), (1, 2), (2, 3), (3, 0)]), 4))
def components(H, N):
    D = 2 ** N; parent = list(range(D))
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for x in range(D):
        for y in range(D):
            if x != y and abs(H[x, y]) > 0: parent[f(x)] = f(y)
    comps = {}
    for x in range(D): comps.setdefault(f(x), set()).add(x)
    return [frozenset(c) for c in comps.values()]
def projectors_in_kernel(H, N):
    D = 2 ** N; I = np.eye(D)
    L = -1j * (lr(H, I) - lr(I, H))
    for l in range(N):
        Z = op(N, {l: PM["Z"]}); L += (l + 1) * (lr(Z, Z) - np.eye(D * D))
    comps = components(H, N); seen = set(); ok = True; even = odd = 0
    for c in comps:
        if c in seen: continue
        cb = frozenset((D - 1) ^ x for x in c)
        if cb not in comps: return False, 0, 0          # the complement must map components onto components
        seen |= {c, cb}
        Pc = np.diag([1.0 if x in c else 0.0 for x in range(D)]); Pb = np.diag([1.0 if x in cb else 0.0 for x in range(D)])
        for v, kind in (((Pc + Pb).reshape(-1), "even"), ((Pc - Pb).reshape(-1), "odd")):
            if not np.any(v): continue
            ok = ok and np.abs(L @ v).max() == 0.0
            if kind == "even": even += 1
            else: odd += 1
    return ok, even, odd
for name, H, N in G6CASES:
    pred = graph_split(H, N)
    okp, pe, po = projectors_in_kernel(H, N)
    check(f"G6 {name}: the component projectors P_C +- P_Cbar are exact kernel vectors ({pe} even, {po} odd)",
          okp and (pe, po) == pred)
    got = kernel_split_mod_p(H, N)
    check(f"G6 {name}: kernel Ad_X split even:odd = (p+s):p = {pred[0]}:{pred[1]}", got == pred, f"nullity mod p {got[0]}:{got[1]}")

# ---------------------------------------------------------------- G7
def popc(x): return bin(x).count("1")
for N, topo in ((3, "chain"), (4, "chain"), (4, "ring")):
    bonds = bonds_chain(N) if topo == "chain" else [(i, (i + 1) % N) for i in range(N)]
    J = [1, 2, 1, 3][:len(bonds)] if topo == "chain" else [1, 2, 1, 3]
    D = 2 ** N
    H = sum(J[k] * (op(N, {a: PM["X"], b: PM["X"]}) + op(N, {a: PM["Y"], b: PM["Y"]}) + 0.5 * op(N, {a: PM["Z"], b: PM["Z"]}))
            for k, (a, b) in enumerate(bonds))
    gam = [0.25 * (l + 1) for l in range(N)]; Gam = sum(gam)
    I = np.eye(D)
    L = -1j * (lr(H, I) - lr(I, H)) + sum(gam[l] * dissip(op(N, {l: SM})) for l in range(N))
    tri = True
    for r in range(D * D):
        for c in range(D * D):
            if L[r, c] != 0:
                a, b = divmod(r, D); a2, b2 = divmod(c, D)
                d1, d2 = popc(a) - popc(a2), popc(b) - popc(b2)
                if not ((d1, d2) == (0, 0) or (d1, d2) == (-1, -1)): tri = False
    check(f"G7 {topo} N={N}: L has entries only within a joint-popcount block or one step down in both", tri)
    n_op = sum(gam[l] * op(N, {l: np.diag([0, 1]).astype(complex)}) for l in range(N))
    Heff = H - 0.5j * n_op
    Xn = op(N, {l: PM["X"] for l in range(N)})
    ident = Xn @ Heff @ Xn - (Heff.conj().T - 0.5j * Gam * I)
    check(f"G7 {topo} N={N}: X^N H_eff X^N = H_eff^dagger - (i/2) Gamma exactly", np.abs(ident).max() == 0.0,
          f"max {np.abs(ident).max():.1e}")
    blockdiag = L.copy()
    for r in range(D * D):
        for c in range(D * D):
            a, b = divmod(r, D); a2, b2 = divmod(c, D)
            if not (popc(a) == popc(a2) and popc(b) == popc(b2)): blockdiag[r, c] = 0
    eff = -1j * (lr(Heff, I) - lr(I, Heff.conj().T))
    for r in range(D * D):
        for c in range(D * D):
            a, b = divmod(r, D); a2, b2 = divmod(c, D)
            if not (popc(a) == popc(a2) and popc(b) == popc(b2)): eff[r, c] = 0
    check(f"G7 {topo} N={N}: the diagonal blocks are -i(H_eff rho - rho H_eff^dagger) exactly", np.abs(blockdiag - eff).max() == 0.0,
          f"max {np.abs(blockdiag - eff).max():.1e}")
    ev = np.linalg.eigvals(L)
    from scipy.optimize import linear_sum_assignment
    C = np.abs(ev[:, None] - (-Gam - ev)[None, :]); rr, cc = linear_sum_assignment(C)
    print(f"     reading: {topo} N={N} palindrome distance about -Gamma/2 = {C[rr, cc].max():.1e}")
    for oname, oterms in (("XX bonds alone", [("X", "X", 1)]), ("DM bonds XY - YX", [("X", "Y", 1), ("Y", "X", -1)])):
        Ho = sum(J[k] * s * op(N, {a: PM[p], b: PM[q]}) for k, (a, b) in enumerate(bonds) for p, q, s in oterms)
        Lo = -1j * (lr(Ho, I) - lr(I, Ho)) + sum(gam[l] * dissip(op(N, {l: SM})) for l in range(N))
        evo = np.linalg.eigvals(Lo); Co = np.abs(evo[:, None] - (-Gam - evo)[None, :]); ro, co = linear_sum_assignment(Co)
        print(f"     reading: {topo} N={N} {oname}, outside the conditions: palindrome distance = {Co[ro, co].max():.1e}")
        if oname.startswith("DM"):
            Heo = Ho - 0.5j * n_op
            okc = np.abs(Xn @ Ho @ Xn - Ho.conj()).max() == 0.0 and np.abs(Xn @ Ho @ Xn - Ho).max() > 0
            okd = np.abs(Xn @ Heo @ Xn - (Heo.conj() - 0.5j * Gam * I)).max() == 0.0
            check(f"G7 {topo} N={N} DM: X^N H X^N = conj(H) != H and X^N H_eff X^N = conj(H_eff) - (i/2) Gamma exactly", okc and okd)

# ---------------------------------------------------------------- G8
def g8(N, terms, bonds=None):
    D = 2 ** N; Id = np.eye(D * D)
    H = ham(N, terms, bonds if bonds is not None else bonds_chain(N))
    LH = -1j * (lr(H, np.eye(D)) - lr(np.eye(D), H))
    LD = sum(dissip(op(N, {l: SM})) for l in range(N))
    light = sum(0.25 * (lr(op(N, {l: PM["X"]}), op(N, {l: PM["X"]})) + lr(op(N, {l: PM["Y"]}), op(N, {l: PM["Y"]})) - 2 * Id) for l in range(N))
    M = LD - light
    Dl = light + (N / 2) * Id
    A = LH @ Dl + Dl @ LH; B = LH @ M + M @ LH
    nrm = lambda X: float(np.real(np.vdot(X, X)))
    return nrm(A), nrm(B), float(np.real(np.vdot(A, B))), nrm(LH)
FORMS = {"XX": ([("X", "X")], lambda N: 32 * (N - 1) ** 2),
         "XY": ([("X", "X"), ("Y", "Y")], lambda N: 32 * (N - 1) * (2 * N - 3)),
         "ZZ": ([("Z", "Z")], lambda N: 32 * (N * N - N - 1)),
         "Heisenberg": ([("X", "X"), ("Y", "Y"), ("Z", "Z")], lambda N: 32 * (3 * N * N - 6 * N + 2))}
for name, (terms, form) in FORMS.items():
    for N in (2, 3, 4, 5, 6):
        A2, B2, AB, LH2 = g8(N, terms)
        check(f"G8 {name} chain N={N}: lights and move orthogonal, move part = {form(N)}*4^(N-2)",
              AB == 0.0 and B2 == form(N) * 4 ** (N - 2) and A2 == (N - 2) / 2 * LH2,
              f"B2 {B2:.0f}, cross {AB}, A2/|L_H|^2 {A2 / LH2}")

GRAPHS = {"ring": lambda N: [(i, (i + 1) % N) for i in range(N)], "star": lambda N: [(0, i) for i in range(1, N)],
          "complete": lambda N: [(i, j) for i in range(N) for j in range(i + 1, N)]}
PER_BOND = {"XX": ([("X", "X")], lambda N: 32 * (N - 1)), "XY": ([("X", "X"), ("Y", "Y")], lambda N: 32 * (2 * N - 3))}
for name, (terms, per) in PER_BOND.items():
    for gname, gb in GRAPHS.items():
        for N in (3, 4, 5):
            E = gb(N)
            A2, B2, AB, LH2 = g8(N, terms, E)
            check(f"G8 {name} {gname} N={N}: lights and move orthogonal, move part = |E|*{per(N)}*4^(N-2), bond-additive",
                  AB == 0.0 and B2 == len(E) * per(N) * 4 ** (N - 2) and A2 == (N - 2) / 2 * LH2,
                  f"B2 {B2:.0f}, cross {AB}")
def adj_pairs(N, E):                                        # sum over sites of C(deg, 2): bond pairs sharing a site
    deg = [sum(1 for e in E if v in e) for v in range(N)]
    return sum(d * (d - 1) // 2 for d in deg)
CLOSED = {"ZZ": ([("Z", "Z")], lambda N, E: 32 * ((N - 1) * len(E) + adj_pairs(N, E))),
          "Heisenberg": ([("X", "X"), ("Y", "Y"), ("Z", "Z")], lambda N, E: 32 * ((3 * N - 4) * len(E) + adj_pairs(N, E)))}
for name, (terms, form) in CLOSED.items():
    for gname, gb in list(GRAPHS.items()) + [("chain", bonds_chain)]:
        for N in (3, 4, 5):
            E = gb(N)
            A2, B2, AB, LH2 = g8(N, terms, E)
            check(f"G8 {name} {gname} N={N}: move part = {form(N, E)}*4^(N-2) = 32[(..)|E| + sum_v C(deg v, 2)], orthogonal",
                  AB == 0.0 and B2 == form(N, E) * 4 ** (N - 2) and A2 == (N - 2) / 2 * LH2, f"B2 {B2:.0f}")
_, zz_chain, _, _ = g8(4, [("Z", "Z")])                    # so Ising is not bond-additive: per bond, chain != ring
_, zz_ring, _, _ = g8(4, [("Z", "Z")], GRAPHS["ring"](4))
check("G8 ZZ N=4: move part 5632 on the chain (3 bonds), 8192 on the ring (4 bonds), so no single per-bond value",
      zz_chain == 5632 and zz_ring == 8192 and zz_chain * 4 != zz_ring * 3, f"chain {zz_chain:.0f}, ring {zz_ring:.0f}")

# ---------------------------------------------------------------- G9
def blocks(L):                                              # connected components of L's sparsity graph on |a><b|
    n = L.shape[0]; parent = list(range(n))
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    rows, cols = np.nonzero(L)
    for r, c in zip(rows, cols): parent[f(r)] = f(c)
    lab = [f(x) for x in range(n)]
    return lab, len(set(lab))
for N in (3, 4):
    D = 2 ** N; I = np.eye(D); b = bonds_chain(N)
    pc = [popc(x) for x in range(D)]
    diff = [pc[a] - pc[c] for a in range(D) for c in range(D)]
    deph = sum(0.25 * (l + 1) * dissip(op(N, {l: SM})) + 0.125 * (l + 2) * dissip(op(N, {l: PM["Z"]})) for l in range(N))
    xxz = ham(N, [("X", "X"), ("Y", "Y")], b) + 0.5 * ham(N, [("Z", "Z")], b)
    cases = (("XXZ", xxz, 2 * N + 1, "diff"), ("XX bonds alone", ham(N, [("X", "X")], b), 2, "parity"),
             ("XXZ + X field on site 0", xxz + 0.5 * op(N, {0: PM["X"]}), 1, None))
    for name, H, expect, label in cases:
        L = -1j * (lr(H, I) - lr(I, H)) + deph
        lab, nb = blocks(L)
        if label == "diff":
            keyed = all(len({diff[k] for k in range(D * D) if lab[k] == r}) == 1 for r in set(lab))
            labels_ok = keyed and len({diff[k] for k in range(D * D)}) == nb
        elif label == "parity":
            labels_ok = all(len({diff[k] % 2 for k in range(D * D) if lab[k] == r}) == 1 for r in set(lab))
        else:
            labels_ok = True
        check(f"G9 N={N} {name} + amplitude damping + Z dephasing: {nb} blocks (expected {expect})"
              + ({"diff": ", one per value of p - q", "parity": ", the two n_XY parities"}.get(label, "")),
              nb == expect and labels_ok)

# ---------------------------------------------------------------- G10
def pmul(A, B):                                             # operators as {Pauli string: Gaussian rational (re, im)}
    out = {}
    for a, (ar, ai) in A.items():
        for b, (br, bi) in B.items():
            ph, w = smul(a, b)
            cr, ci = ar * br - ai * bi, ar * bi + ai * br
            pr, pi = int(ph.real), int(ph.imag)
            o = out.get(w, (F(0), F(0))); out[w] = (o[0] + cr * pr - ci * pi, o[1] + cr * pi + ci * pr)
    return {w: c for w, c in out.items() if c != (F(0), F(0))}
def odd_word_trace(N, terms, hp, js):
    H = {}
    for b in range(N - 1):
        for t in terms:
            w = ["I"] * N; w[b], w[b + 1] = t[0], t[1]; H["".join(w)] = (F(1), F(0))
    one = {"I" * N: (F(1), F(0))}; pw = [one]
    for _ in range(max(hp)): pw.append(pmul(pw[-1], H))
    W = one
    for k, j in zip(hp, js):
        W = pmul(pmul(W, pw[k]), {"".join("Z" if i == j else "I" for i in range(N)): (F(1), F(0))})
    return W.get("I" * N, (F(0), F(0)))                     # trace / 2^N
ODD_WORDS = {  # pair: (powers of H, jump sites, first N, trace/2^N as Gaussian integer (re, im))
    "XX+XY": ([0, 1, 3], [0, 0, 1], 3, (0, 8)), "XX+YX": ([0, 1, 3], [0, 0, 1], 3, (0, -8)),
    "XY+YY": ([0, 1, 3], [0, 0, 1], 3, (0, 8)), "YX+YY": ([0, 1, 3], [0, 0, 1], 3, (0, -8)),
    "XY+ZX": ([0, 1, 2], [0, 1, 2], 3, (0, 4)), "XZ+YX": ([0, 1, 2], [0, 1, 2], 3, (0, 4)),
    "XY+YZ": ([0, 1, 2], [0, 1, 2], 3, (0, -4)), "YX+ZY": ([0, 1, 2], [0, 1, 2], 3, (0, -4)),
    "YZ+ZX": ([0, 1, 3], [0, 1, 2], 4, (0, 4)), "XZ+ZY": ([0, 1, 3], [0, 1, 2], 4, (0, -4)),
    "XY+XZ": ([0, 0, 1, 1, 1], [0, 1, 0, 1, 2], 3, (0, 4)), "YX+YZ": ([0, 0, 1, 1, 1], [0, 1, 0, 1, 2], 3, (0, -4)),
    "XY+ZY": ([0, 0, 1, 1, 1], [0, 1, 1, 2, 2], 3, (0, 4)), "YX+ZX": ([0, 0, 1, 1, 1], [0, 1, 1, 2, 2], 3, (0, -4))}
for pair, (hp, js, n0, want) in ODD_WORDS.items():
    P = sum(hp); top = 2 * P + 6                            # polynomial of degree <= P//2 <= 2 from N = 2P + 4
    vals = [odd_word_trace(N, pair.split("+"), hp, js) for N in range(n0, top + 1)]
    check(f"G10 {pair}: odd word, H powers {hp} between jumps Z at sites {js}, trace/2^N = {want[0]}+{want[1]}i "
          f"at every N = {n0}..{top} (constant past N = {2 * P + 4} on three points)",
          all(v == (F(want[0]), F(want[1])) for v in vals), f"{sorted(set(vals))}")
N3_WORDS = {"YZ+ZX": ([0, 1, 3], [0, 2, 1], (0, 4)), "XZ+ZY": ([0, 1, 3], [0, 2, 1], (0, -4))}
for pair, (hp, js, want) in N3_WORDS.items():
    v = odd_word_trace(3, pair.split("+"), hp, js)
    check(f"G10 {pair} N=3: odd word, H powers {hp}, jumps at {js}, trace/2^N = {want[1]}i", v == (F(want[0]), F(want[1])))
def gauss_nullity(M):                                       # exact over Q(i); M has Gaussian-integer entries
    A = [[(F(int(round(x.real))), F(int(round(x.imag)))) for x in row] for row in M]
    n, m = len(A), len(A[0]); r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c] != (0, 0)), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        a, b = A[r][c]; d = a * a + b * b; inv = (a / d, -b / d)
        A[r] = [(x * inv[0] - y * inv[1], x * inv[1] + y * inv[0]) for x, y in A[r]]
        for i in range(n):
            if i != r and A[i][c] != (0, 0):
                f0, f1 = A[i][c]
                A[i] = [(x - (f0 * u - f1 * v), y - (f0 * v + f1 * u)) for (x, y), (u, v) in zip(A[i], A[r])]
        r += 1
    return m - r
U_LETTERS = {"X": PM["X"], "Y": PM["Y"], "X-Y": PM["X"] - PM["Y"]}
HARD14 = set(ODD_WORDS)
for a, b in itertools.combinations([p + q for p in "XYZ" for q in "XYZ"], 2):
    pair = a + "+" + b; found = {}
    for N in (3, 4):
        H = ham(N, [tuple(a), tuple(b)], bonds_chain(N))
        found[N] = next(((ue, uo) for ue, uo in itertools.product(U_LETTERS, repeat=2)
                         if np.array_equal(op(N, {l: U_LETTERS[ue if l % 2 == 0 else uo] for l in range(N)}) @ H,
                                           H @ op(N, {l: U_LETTERS[ue if l % 2 == 0 else uo] for l in range(N)}))), None)
    want = pair not in HARD14
    check(f"G10 {pair}: a period-two product mirror (X, Y or X-Y per site) commuting with H "
          + ("exists, the same at N = 3 and 4" if want else "exists at neither N = 3 nor 4"),
          (found[3] is not None and found[3] == found[4]) if want else (found[3] is None and found[4] is None), f"{found}")
END_COUNTS = {p: ((2, 0) if p in ("XX+XY", "XX+YX", "XY+YY", "YX+YY") else (1, 0)) for p in ODD_WORDS}
END_COUNTS["XX+XZ"] = (1, 1)
for pair, expect in END_COUNTS.items():
    N = 3; D = 2 ** N; I = np.eye(D)
    H = ham(N, [tuple(x) for x in pair.split("+")], bonds_chain(N))
    L = -1j * (lr(H, I) - lr(I, H)) + sum(lr(op(N, {l: PM["Z"]}), op(N, {l: PM["Z"]})) - np.eye(D * D) for l in range(N))
    got = (gauss_nullity(L), gauss_nullity(L + 2 * N * np.eye(D * D)))
    check(f"G10 {pair} N=3, unit rates: dim ker L, dim ker(L + 2 sigma) = {expect}"
          + (" (no palindrome, F158)" if expect[0] != expect[1] else " (palindrome; the control)"), got == expect, f"{got}")

# ---------------------------------------------------------------- G11
CLIFF1 = []; CLIFF1S = []                                   # single-site Cliffords: letter maps, and with their signs (det +1)
for perm in itertools.permutations(range(3)):
    par = 1 if sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j]) % 2 == 0 else -1
    for sg in itertools.product((1, -1), repeat=3):
        if sg[0] * sg[1] * sg[2] * par == 1:
            CLIFF1.append({"XYZ"[i]: "XYZ"[perm[i]] for i in range(3)})
            CLIFF1S.append(tuple((sg[i], "XYZ"[perm[i]]) for i in range(3)))
def commutant_set(Q): return {a + b for a in "XYZ" for b in "XYZ" if ((a != Q) + (b != Q)) % 2 == 0}
check("G11 the 24 single-site Cliffords; Part 1's bilinears are the strings commuting with X x X",
      len(CLIFF1) == 24 and commutant_set("X") == {"XX", "YY", "YZ", "ZY", "ZZ"})
N11 = 3; orbit11 = {}; stab11 = 0; uniform_letters = True
for cs in itertools.product(range(24), repeat=N11):
    jl = {CLIFF1[c]["Z"] for c in cs}
    if len(jl) != 1: continue
    Pj = jl.pop()
    imgs = {CLIFF1[cs[l]][w[0]] + CLIFF1[cs[l + 1]][w[1]] for l in range(N11 - 1) for w in commutant_set("X")}
    hit = [Q for Q in "XYZ" if Q != Pj and imgs == commutant_set(Q)]
    if hit:
        orbit11[(Pj, hit[0])] = orbit11.get((Pj, hit[0]), 0) + 1
        uniform_letters &= len({tuple(sorted(CLIFF1[c].items())) for c in cs}) == 1
        stab11 += (Pj, hit[0]) == ("Z", "X")
check("G11 N=3, all 24^3 site-dependent local Cliffords: Part 1 (jumps Z, commutant X) reaches all six settings (P, Q), "
      "64 each, every hit one letter permutation on all sites", len(orbit11) == 6 and set(orbit11.values()) == {64}
      and uniform_letters, f"{orbit11}")
stab_signed = set()
for cs in itertools.product(range(24), repeat=N11):
    if all(CLIFF1[c]["Z"] == "Z" for c in cs):
        imgs = {CLIFF1[cs[l]][w[0]] + CLIFF1[cs[l + 1]][w[1]] for l in range(N11 - 1) for w in commutant_set("X")}
        if imgs == commutant_set("X"): stab_signed.add(tuple(CLIFF1S[c] for c in cs))
PAULI1 = {tuple((s_, l_) for s_, l_ in zip(sg, "XYZ")) for sg in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]}
check("G11 the stabilizer of Part 1 is the local Pauli group: the distinct signed products fixing (Z, X) are exactly "
      "the 4^3 = 64 products of I, X, Y, Z (sign patterns of the letter-fixing maps)",
      len(stab_signed) == 64 and all(all(site in PAULI1 for site in prod) for prod in stab_signed), f"{len(stab_signed)}")
VX = np.eye(2) + 1j * PM["X"]; VXi = (np.eye(2) - 1j * PM["X"]) / 2   # the quarter turn about X, unnormalised: exact entries
imgs_q = {}
for w in commutant_set("X"):
    M = np.kron(VX, VX) @ np.kron(PM[w[0]], PM[w[1]]) @ np.kron(VXi, VXi)
    hits = [(a + b, sg) for a in "XYZ" for b in "XYZ" for sg in (1, -1) if np.array_equal(M, sg * np.kron(PM[a], PM[b]))]
    imgs_q[w] = hits
zimg = VX @ PM["Z"] @ VXi
check("G11 the quarter turn about X sends Z to +-Y and Part 1's five bilinears onto the same five strings (Part 3's set)",
      (np.array_equal(zimg, -PM["Y"]) or np.array_equal(zimg, PM["Y"])) and {h[0][0] for h in imgs_q.values() if h} == commutant_set("X")
      and all(len(h) == 1 for h in imgs_q.values()), f"{imgs_q}")
for (Pj, Qc) in [("Z", "Y"), ("X", "Y"), ("Y", "Z")]:           # the three settings no named Part covers: F158 end counts
    D = 2 ** N11; I = np.eye(D)
    H = sum((k + 1) * op(N11, {b: PM[w[0]], b + 1: PM[w[1]]}) for b in range(N11 - 1)
            for k, w in enumerate(sorted(commutant_set(Qc))))
    def Lof(Hh):
        return -1j * (lr(Hh, I) - lr(I, Hh)) + sum(lr(op(N11, {l: PM[Pj]}), op(N11, {l: PM[Pj]})) - np.eye(D * D) for l in range(N11))
    k0, k1 = gauss_nullity(Lof(H)), gauss_nullity(Lof(H) + 2 * N11 * np.eye(D * D))
    H_bad = H + op(N11, {0: PM[Pj]})                            # a field along the jump on one site: outside the setting
    b0, b1 = gauss_nullity(Lof(H_bad)), gauss_nullity(Lof(H_bad) + 2 * N11 * np.eye(D * D))
    check(f"G11 setting (jumps {Pj}, commutant {Qc}), unit rates, N=3: F158 end counts equal ({k0}, {k1}); "
          f"with a {Pj} field on site 0 they differ ({b0}, {b1})", k0 == k1 and b0 != b1)

print("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
