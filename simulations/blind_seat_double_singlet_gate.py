"""Gate: the double singlet at N = 5 (PROOF_ONE_DEPHASED_SEAT, section on blind dephased seats, bullet
"At N = 5, the double singlet" under the paragraph "The sheet and blindness, at N = 3 and N = 4").

Setting: N = 5, Heisenberg bonds sum J_ij sigma_i.sigma_j (integer couplings, signs free, zero allowed), the jump Z on
seat d = 0, planar integer fields on the sites 1..4, none on 0; bit 1 = excited; s_ij = |1_i 0_j> - |0_i 1_j> the
singlet of a pair (unnormalised here; |s s><s s| below is the projector onto s_ij s_kl); m1 = s_12 s_34, m2 = s_13 s_24; a pairing {i,j},{k,l} of 1..4 is a double
singlet when J_0i = J_0j, J_0k = J_0l and J_ik + J_jl = J_il + J_jk. Exact or certified, no threshold.

1. The lemma's algebraic steps, symbolically in the couplings, every operator built from the bit action: the
   two-excitation states of 1..4 killed by the total lowering are exactly span{m1, m2} (the singlets of 1..4); the
   commutators [S+, B-] and [S+, Bz] of the components of B = sum_r J_0r sigma_r with the raising S+ of 1..4 are
   multiples of Bz and B+; Bz (x m1 + z m2) = 0 is exactly the three equations z (J01 + J02 - J03 - J04) = 0,
   x (J01 - J02) = z (J02 - J04), x (J03 - J04) = -z (J02 - J04), and s_14 s_23 = m2 - m1; for each pairing
   the cross bonds act on s_ij s_kl as (J_ik + J_jl - J_il - J_jk) sigma_i.sigma_k, and sigma_i.sigma_k (s_ij s_kl) is
   orthogonal to s_ij s_kl with three times its squared norm. X^5 commutes with H_b at generic couplings.
2. On the double-singlet locus, symbolically in the couplings and the fields: |0>_0 s_ij s_kl and |1>_0 s_ij s_kl are
   eigenvectors of H_b at -3 J_ij - 3 J_kl; with the fields equal within each pair, X_0 (x) |s s><s s| and
   Y_0 (x) |s s><s s| commute with H, anticommute with Z_0, are odd under Z^5 and independent, the projectors onto
   |0>_0 s s and |1>_0 s s commute with H, and F+ and F- kill |0>_0 s_ij s_kl;
   the field part on s_ij s_kl is (f_i - f_j).sigma_i + (f_k - f_l).sigma_k applied to it, the six vectors sigma_i^x,y,z (s s),
   sigma_k^x,y,z (s s) independent, so an element O_0 (x) |s s><s s| commutes with the fields exactly when they are equal
   within each pair.
3. Rows, by ranks modulo two primes (a kernel modulo p is at least the rational one, so a zero is exact; a full rank
   certifies, so seat 0's Krylov rank 5 certifies that 0 is not blind): drawn rows, planted double singlets (some
   planted at -3 J_ij - 3 J_kl = E0, the sum of the couplings) and near misses (J_0i, J_0k or the cross sum broken by
   one, off both kinds of planted row), every one with 0 not blind, each group's construction checked. Predicted: the d-unexcited eigenvectors of H_b in the two-excitation sector
   (the largest subspace of V0 that H_b keeps) have dimension 1 at a double singlet and 0 elsewhere, and the fence at
   n = 3 (vectors with 0 excited at E0) has dimension 1 at a double singlet with -3 J_ij - 3 J_kl = E0 and 0 elsewhere;
   the ones are bounded below exactly by part 2 and by X^5 applied to |0>_0 s_ij s_kl.
4. Far elements at planted rows: homogeneous under Z^5, counted by parity by a rank modulo p, predicted (odd, even) =
   (3, 0) with the fields equal within the pairs on one line, (2, 0) off one line, and class A's (1, 0) or (0, 0)
   with a pair unequal or opposite; the 3 and 2 bounded below by class A's colouring (checked exactly at the rows on
   the x line) and the two elements of part 2. Condition (I): the cyclic space alg(H_b, Z_0, F+)|vac> whole with the
   fields opposite within the pairs, and read short with them equal (part 2's subspace forces that). The full two-end
   count (F158) at the N = 5 row of the bullet 'At N = 4, a seat that is not blind': near 3 and far 2 with X on (1, 2) and 2Y on (3, 4), 3 and 3 with 3X on
   (3, 4), and 3 and 2 at a planted double singlet off the vacuum's energy, where (II) holds; the near 3 bounded below by 1 and
   the projectors onto |0>_0 s s and |1>_0 s s, the far 2 by part 2.
5. A row off the locus whose two-excitation sector is reducible: J = (4, 0, 0, 6, 0, -3, -6, 4, 6, -6) on the edges
   01, 02, 03, 04, 12, 13, 14, 23, 24, 34; 0 not blind (exactly), no double singlet, the commutant of
   alg(H_b, Z_0) on the two-excitation sector two-dimensional over Q (exact rank), no d-unexcited eigenvector in the
   two-excitation sector, the fence at n = 3 empty (exact), and at drawn field sets (I) and far elements of class A only.

Run: python simulations/blind_seat_double_singlet_gate.py   (writes results/blind_seat_double_singlet_gate.txt)
"""
import itertools
import random
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import pauli_product

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")   # the console's default code page cannot print every character below

OUT, FAILS = [], []
T0 = time.time()


def say(s=""):
    print(s)
    OUT.append(s)


def check(cond, what):
    say(("  PASS  " if cond else "  FAIL  ") + what)
    if not cond:
        FAILS.append(what)


def _is_prime(q):
    return q > 1 and all(q % r for r in range(2, int(q ** 0.5) + 1))


def _prime_1mod4(start):
    q = start
    while not (_is_prime(q) and q % 4 == 1):
        q += 1
    return q


PRIMES = [_prime_1mod4(10 ** 9), _prime_1mod4(10 ** 9 + 100)]


def sqrt_m1(p):
    g = next(g for g in range(2, 200) if pow(g, (p - 1) // 2, p) == p - 1)
    return pow(g, (p - 1) // 4, p)


def rank_mod(M, p):
    M = np.array(np.array(M, dtype=object) % p, dtype=np.int64)
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r


N, D = 5, 32
EDGES = list(itertools.combinations(range(N), 2))
PAIRINGS = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]
S2 = [x for x in range(D) if bin(x).count('1') == 2]
S3 = [x for x in range(D) if bin(x).count('1') == 3]
V0 = [k for k, x in enumerate(S2) if not x & 1]
V1 = [k for k, x in enumerate(S2) if x & 1]


def g(b, i, j):
    return b.get((min(i, j), max(i, j)), 0)


def ds_pairings(b):
    return [((i_, j_), (k_, l_)) for (i_, j_), (k_, l_) in PAIRINGS
            if g(b, 0, i_) == g(b, 0, j_) and g(b, 0, k_) == g(b, 0, l_)
            and g(b, i_, k_) + g(b, j_, l_) == g(b, i_, l_) + g(b, j_, k_)]


def connected(b):
    adj = {i: set() for i in range(N)}
    for (i, j), J in b.items():
        if J:
            adj[i].add(j); adj[j].add(i)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x] - seen:
            seen.add(y); st.append(y)
    return len(seen) == N


def krylov_rank_mod(b, p):
    L = np.zeros((N, N), dtype=np.int64)
    for (i, j), J in b.items():
        L[i, j] -= J; L[j, i] -= J; L[i, i] += J; L[j, j] += J
    v = np.zeros(N, dtype=np.int64); v[0] = 1
    K = [v]
    for _ in range(N - 1):
        K.append(L.dot(K[-1]) % p)
    return rank_mod(np.array(K), p)


def not_blind(b):
    return any(krylov_rank_mod(b, p) == N for p in PRIMES)


def hb_int(b):
    M = np.zeros((D, D), dtype=np.int64)
    for x in range(D):
        for (i, j), J in b.items():
            bi, bj = (x >> i) & 1, (x >> j) & 1
            M[x, x] += J * (1 if bi == bj else -1)
            if bi != bj:
                M[x ^ (1 << i) ^ (1 << j), x] += 2 * J
    return M


def joint_dim(b):
    """dim of the largest subspace of V0 (0 unexcited, two excitations) kept by H_b, modulo each prime; the min."""
    H2 = hb_int(b)[np.ix_(S2, S2)]
    best = None
    for p in PRIMES:
        A = H2[np.ix_(V0, V0)] % p; Bm = H2[np.ix_(V1, V0)] % p
        St, Ak = [], np.eye(6, dtype=np.int64)
        for _ in range(6):
            St.append(Bm.dot(Ak) % p); Ak = Ak.dot(A) % p
        k = 6 - rank_mod(np.vstack(St), p)
        best = k if best is None else min(best, k)
    return best


def fence3(b):
    H = hb_int(b); E0 = sum(b.values())
    cols = [x for x in S3 if x & 1]
    best = None
    for p in PRIMES:
        M = (H[:, cols] - E0 * np.eye(D, dtype=np.int64)[:, cols]) % p
        k = len(cols) - rank_mod(M, p)
        best = k if best is None else min(best, k)
    return best


# ---- symbolic operators on the 32-dimensional space ----------------------------------------------------------------------
def e(*exc):
    v = sp.zeros(D, 1); v[sum(1 << s for s in exc)] = 1
    return v


def ss(i_, j_, k_, l_):
    """s_ij (x) s_kl on the sites 1..4, site 0 unexcited."""
    return e(i_, k_) - e(i_, l_) - e(j_, k_) + e(j_, l_)


def pauli_op(L, k):
    M = sp.zeros(D, D)
    for x in range(D):
        bit = (x >> k) & 1
        if L == 'Z':
            M[x, x] = -1 if bit else 1
        elif L == 'X':
            M[x ^ (1 << k), x] = 1
        else:   # Y|0> = i|1>, Y|1> = -i|0>
            M[x ^ (1 << k), x] = sp.I if bit == 0 else -sp.I
    return M


SIG = {(L, k): pauli_op(L, k) for L in 'XYZ' for k in range(N)}


def dot(a, b):
    return sum((SIG[(L, a)] * SIG[(L, b)] for L in 'XYZ'), sp.zeros(D, D))


def raise_op(k):   # adds an excitation at k
    M = sp.zeros(D, D)
    for x in range(D):
        if not (x >> k) & 1:
            M[x | (1 << k), x] = 1
    return M


def is_zero(M):
    return all(sp.expand(z) == 0 for z in M)


# ---- part 1 -------------------------------------------------------------------------------------------------------------
say("Part 1. The lemma's algebraic steps, symbolically")
x = sp.symbols('J01 J02 J03 J04')
Sp = sum((raise_op(k) for k in range(1, 5)), sp.zeros(D, D))
Sm = Sp.T
R2 = [xx for xx in S2 if not xx & 1]
low = sp.Matrix([[Sm[r, c] for c in R2] for r in range(D)])
ker = low.nullspace()
m_, n_ = ss(1, 2, 3, 4), ss(1, 3, 2, 4)
emb = lambda v: sp.Matrix([v[r] for r in R2])
check(len(ker) == 2 and sp.Matrix.hstack(emb(m_), emb(n_)).rank() == 2 and sp.Matrix.hstack(*ker, emb(m_), emb(n_)).rank() == 2
      and is_zero(Sp * m_) and is_zero(Sp * n_),
      "the two-excitation states of 1..4 killed by the total lowering of 1..4 are exactly span{m1, m2}, m1 = s_12 s_34, "
      "m2 = s_13 s_24 independent, both killed by the total raising of 1..4 too")
Bz = sum((x[r - 1] * SIG[('Z', r)] for r in range(1, 5)), sp.zeros(D, D))
Bp = sum((x[r - 1] * raise_op(r) for r in range(1, 5)), sp.zeros(D, D))
Bm = Bp.T
check(is_zero((Sp * Bm - Bm * Sp) + Bz) and is_zero((Sp * Bz - Bz * Sp) - 2 * Bp),
      "[S+, B-] = -Bz and [S+, Bz] = 2 B+ for B = sum_r J_0r sigma_r on 1..4 (raising and lowering components), so a "
      "singlet with B- of it zero has Bz and B+ of it zero")
al, be = sp.symbols('alpha beta')
v = Bz * (al * m_ + be * n_)
eqs = [sp.expand(c) for c in v if sp.expand(c) != 0]
claimed = [be * (x[0] + x[1] - x[2] - x[3]), al * (x[0] - x[1]) - be * (x[1] - x[3]), al * (x[2] - x[3]) + be * (x[1] - x[3])]
coef = lambda q: [sp.expand(q).coeff(al).coeff(t) for t in x] + [sp.expand(q).coeff(be).coeff(t) for t in x]
Me, Mc = sp.Matrix([coef(q) for q in eqs]), sp.Matrix([coef(q) for q in claimed])
check(Me.rank() == 3 and Mc.rank() == 3 and Me.col_join(Mc).rank() == 3,
      "Bz (x m1 + z m2) = 0 is exactly z (J01 + J02 - J03 - J04) = 0, x (J01 - J02) = z (J02 - J04), "
      "x (J03 - J04) = -z (J02 - J04)")
check(sp.expand(ss(1, 4, 2, 3) - (n_ - m_)) == sp.zeros(D, 1), "s_14 s_23 = m2 - m1")
Jc = {pq: sp.Symbol(f"J{pq[0]}{pq[1]}") for pq in itertools.combinations(range(1, 5), 2)}
Jg = lambda i, j: Jc[(min(i, j), max(i, j))]
ok_cross = True
for (i_, j_), (k_, l_) in PAIRINGS:
    chi = ss(i_, j_, k_, l_)
    K = Jg(i_, k_) * dot(i_, k_) + Jg(j_, l_) * dot(j_, l_) + Jg(i_, l_) * dot(i_, l_) + Jg(j_, k_) * dot(j_, k_)
    w = dot(i_, k_) * chi
    ok_cross &= sp.expand(K * chi - (Jg(i_, k_) + Jg(j_, l_) - Jg(i_, l_) - Jg(j_, k_)) * w) == sp.zeros(D, 1)
    ok_cross &= sp.expand((chi.T * w)[0]) == 0 and sp.expand((w.T * w)[0] - 3 * (chi.T * chi)[0]) == 0
    ok_cross &= sp.expand((dot(i_, j_) + dot(k_, l_)) * chi + 6 * chi) == sp.zeros(D, 1)
check(ok_cross, "for each pairing the cross bonds act on s_ij s_kl as (J_ik + J_jl - J_il - J_jk) sigma_i.sigma_k, "
                "sigma_i.sigma_k (s_ij s_kl) is orthogonal to s_ij s_kl with three times its squared norm, and "
                "sigma_i.sigma_j + sigma_k.sigma_l gives -6 on it")

Jall = sp.symbols('K01 K02 K03 K04 K12 K13 K14 K23 K24 K34')
Hgen = sum((J * dot(i, j) for (i, j), J in zip(EDGES, Jall)), sp.zeros(D, D))
X5 = SIG[('X', 0)] * SIG[('X', 1)] * SIG[('X', 2)] * SIG[('X', 3)] * SIG[('X', 4)]
check(is_zero(X5 * Hgen - Hgen * X5), "X^5 commutes with H_b at generic couplings (so it carries the fence at n = 3 onto the "
      "d-unexcited eigenvectors with two excitations at the same value)")

# ---- part 2 -------------------------------------------------------------------------------------------------------------
say()
say("Part 2. On the double-singlet locus, symbolically in the couplings and the fields")
pp, qq, jab, jcd, jac, jad, jbc = sp.symbols('p q J_ab J_cd J_ac J_ad J_bc')
fa = sp.symbols('f1x f1y f2x f2y')
ok2 = True
for (i_, j_), (k_, l_) in PAIRINGS:
    bonds = {(0, i_): pp, (0, j_): pp, (0, k_): qq, (0, l_): qq, (min(i_, j_), max(i_, j_)): jab, (min(k_, l_), max(k_, l_)): jcd,
             (min(i_, k_), max(i_, k_)): jac, (min(i_, l_), max(i_, l_)): jad, (min(j_, k_), max(j_, k_)): jbc,
             (min(j_, l_), max(j_, l_)): jad + jbc - jac}
    Hb = sum((J * dot(i, j) for (i, j), J in bonds.items()), sp.zeros(D, D))
    chi0 = ss(i_, j_, k_, l_)
    chi1 = SIG[('X', 0)] * chi0
    lam = -3 * jab - 3 * jcd
    ok2 &= sp.expand(Hb * chi0 - lam * chi0) == sp.zeros(D, 1) and sp.expand(Hb * chi1 - lam * chi1) == sp.zeros(D, 1)
    F = (fa[0] * (SIG[('X', i_)] + SIG[('X', j_)]) + fa[1] * (SIG[('Y', i_)] + SIG[('Y', j_)])
         + fa[2] * (SIG[('X', k_)] + SIG[('X', l_)]) + fa[3] * (SIG[('Y', k_)] + SIG[('Y', l_)]))
    H = Hb + F
    P0, P1 = chi0 * chi0.T / 4, chi1 * chi1.T / 4
    P = P0 + P1                                   # 1_0 (x) |s s><s s|, the projector on the sites 1..4
    ok2 &= is_zero(H * P0 - P0 * H) and is_zero(H * P1 - P1 * H)
    Ws = []
    for L in 'XY':
        W = SIG[(L, 0)] * P
        Ws.append(W)
        Z5 = SIG[('Z', 0)] * SIG[('Z', 1)] * SIG[('Z', 2)] * SIG[('Z', 3)] * SIG[('Z', 4)]
        ok2 &= is_zero(H * W - W * H) and is_zero(SIG[('Z', 0)] * W + W * SIG[('Z', 0)]) and is_zero(Z5 * W * Z5 + W)
    Fp = sum((((fa[0] + sp.I * fa[1]) if k in (i_, j_) else (fa[2] + sp.I * fa[3])) * raise_op(k) for k in range(1, 5)),
             sp.zeros(D, D))
    Fm = Fp.H
    ok2 &= sp.Matrix.hstack(*[sp.Matrix(Wm.reshape(D * D, 1)) for Wm in Ws]).rank() == 2
    ok2 &= sp.expand(Fp * chi0) == sp.zeros(D, 1) and sp.expand(Fm * chi0) == sp.zeros(D, 1)
check(ok2, "for each pairing, at J_0i = J_0j = p, J_0k = J_0l = q, J_jl = J_il + J_jk - J_ik: |0>_0 s_ij s_kl and "
           "|1>_0 s_ij s_kl are eigenvectors of H_b at -3 J_ij - 3 J_kl; with the fields equal within each pair, "
           "X_0 (x) |s s><s s| and Y_0 (x) |s s><s s| commute with H, anticommute with Z_0 and are odd under Z^5, independent (rank 2 together), "
           "the projectors onto |0>_0 s s and |1>_0 s s commute with H, and F+ and F- kill "
           "|0>_0 s_ij s_kl")
gs = sp.symbols('g1x g1y g1z g2x g2y g2z g3x g3y g3z g4x g4y g4z')
gv = {1: gs[0:3], 2: gs[3:6], 3: gs[6:9], 4: gs[9:12]}
ok_f = True
for (i_, j_), (k_, l_) in PAIRINGS:
    chi0 = ss(i_, j_, k_, l_)
    Ffull = sum((gv[k][t] * SIG[('XYZ'[t], k)] for k in range(1, 5) for t in range(3)), sp.zeros(D, D))
    red = sum(((gv[i_][t] - gv[j_][t]) * SIG[('XYZ'[t], i_)] + (gv[k_][t] - gv[l_][t]) * SIG[('XYZ'[t], k_)]
               for t in range(3)), sp.zeros(D, D))
    ok_f &= sp.expand(Ffull * chi0 - red * chi0) == sp.zeros(D, 1)
    six = sp.Matrix.hstack(*[SIG[('XYZ'[t], s)] * chi0 for s in (i_, k_) for t in range(3)])
    ok_f &= six.rank() == 6
check(ok_f, "for each pairing a one-site field sum_r f_r.sigma_r on 1..4 acts on s_ij s_kl as (f_i - f_j).sigma_i + "
            "(f_k - f_l).sigma_k, the six vectors sigma_i^x,y,z (s s), sigma_k^x,y,z (s s) independent: it kills s_ij s_kl "
            "exactly when f_i = f_j and f_k = f_l")

# ---- part 3 -------------------------------------------------------------------------------------------------------------
say()
say("Part 3. Rows: d-unexcited eigenvectors in sector 2 and the fence at n = 3, by ranks modulo two primes")
rng = random.Random(20261119)
CPL = (-3, -2, -1, 0, 1, 2, 3)


def plant(rng, at_e0):
    while True:
        (i_, j_), (k_, l_) = rng.choice(PAIRINGS)
        p_, j_ab, j_cd, j_ac, j_ad, j_bc = (rng.choice(CPL) for _ in range(6))
        q_ = -p_ - 2 * j_ab - 2 * j_cd - j_ad - j_bc if at_e0 else rng.choice(CPL)
        b = {(0, i_): p_, (0, j_): p_, (0, k_): q_, (0, l_): q_, tuple(sorted((i_, j_))): j_ab, tuple(sorted((k_, l_))): j_cd,
             tuple(sorted((i_, k_))): j_ac, tuple(sorted((i_, l_))): j_ad, tuple(sorted((j_, k_))): j_bc,
             tuple(sorted((j_, l_))): j_ad + j_bc - j_ac}
        b = {k: v for k, v in b.items() if v}
        if connected(b) and not_blind(b) and (at_e0 or -3 * j_ab - 3 * j_cd != sum(b.values())):
            return b


def near_miss(rng):
    while True:
        b = dict(plant(rng, rng.random() < 0.5))
        (i_, j_), (k_, l_) = ds_pairings(b)[0]
        k = rng.choice([(0, i_), (0, k_), tuple(sorted((j_, l_)))])
        b[k] = b.get(k, 0) + 1
        b = {kk: vv for kk, vv in b.items() if vv}
        if connected(b) and not_blind(b) and not ds_pairings(b):
            return b


def draw(rng):
    while True:
        b = {ed: rng.choice(CPL) for ed in EDGES}
        b = {k: v for k, v in b.items() if v}
        if connected(b) and not_blind(b):
            return b


groups = {'drawn': [draw(rng) for _ in range(400)],
          'planted, -3 J_ij - 3 J_kl != E0': [plant(rng, False) for _ in range(60)],
          'planted at -3 J_ij - 3 J_kl = E0': [plant(rng, True) for _ in range(30)],
          'near misses': [near_miss(rng) for _ in range(60)]}
BUILT = {'planted, -3 J_ij - 3 J_kl != E0': (1, False), 'planted at -3 J_ij - 3 J_kl = E0': (1, True), 'near misses': (0, False)}
for name, rows in groups.items():
    ok_j = ok_f3 = ok_built = True
    tally = {}
    for b in rows:
        parts = ds_pairings(b)
        ds = 1 if parts else 0
        lam_e0 = bool(parts) and -3 * g(b, *parts[0][0]) - 3 * g(b, *parts[0][1]) == sum(b.values())
        if name in BUILT:
            ok_built &= (ds, lam_e0) == BUILT[name]
        jd, f3 = joint_dim(b), fence3(b)
        tally[(jd, f3)] = tally.get((jd, f3), 0) + 1
        ok_j &= jd == ds
        ok_f3 &= f3 == (1 if lam_e0 else 0)
    say(f"  read, {name}: {len(rows)} rows, 0 not blind at every one, (eigenvector dim, fence dim) -> {sorted(tally.items())}")
    if name in BUILT:
        check(ok_built, f"{name}: every row is what its construction says (double singlet {BUILT[name][0]}, at E0 {BUILT[name][1]})")
    check(ok_j, f"{name}: the d-unexcited eigenvectors in sector 2 have dimension 1 exactly at the rows with a double "
                f"singlet and 0 at the others")
    check(ok_f3, f"{name}: the fence at n = 3 has dimension 1 exactly at the double singlets with -3 J_ij - 3 J_kl = E0 "
                 f"and 0 at the others")
row5 = {(0, 1): 1, (0, 2): 1, (0, 3): 2, (0, 4): 2, (1, 3): 1, (2, 4): 2, (2, 3): 3, (1, 2): -1, (3, 4): -2}
check(not_blind(row5) and ds_pairings(row5) == [((1, 2), (3, 4))] and joint_dim(row5) == 1 and fence3(row5) == 1,
      "the N = 5 row of the bullet 'At N = 4, a seat that is not blind' is a double singlet on (1, 2), (3, 4) with -3 J12 - 3 J34 = 9 = E0: "
      "eigenvector dimension 1, fence dimension 1")

# ---- part 4 -------------------------------------------------------------------------------------------------------------
say()
say("Part 4. Far elements, condition (I) and the two-end count at planted rows")
P0 = PRIMES[0]
I_P = sqrt_m1(P0)
IPOW = [1, I_P, P0 - 1, P0 - I_P]


def anticomm(s, t):
    return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1


def far_count(b, fields, parity):
    terms = {}
    for (i, j), J in b.items():
        for L in 'XYZ':
            s = ['I'] * N; s[i] = L; s[j] = L
            terms[''.join(s)] = J
    for j, (fx, fy) in fields.items():
        for L, c in (('X', fx), ('Y', fy)):
            if c:
                s = ['I'] * N; s[j] = L
                terms[''.join(s)] = c
    span = [''.join(w) for w in itertools.product('IXYZ', repeat=N)
            if w[0] in 'XY' and sum(1 for a in w if a in 'XY') % 2 == parity]
    rowidx, ent = {}, {}
    for col, s in enumerate(span):
        for Pp, c in terms.items():
            if not anticomm(Pp, s):
                continue
            ex, u = pauli_product(Pp, s)
            r = rowidx.setdefault(u, len(rowidx))
            ent[(r, col)] = (ent.get((r, col), 0) + 2 * c * IPOW[ex]) % P0
    M = np.zeros((len(rowidx), len(span)), dtype=np.int64)
    for (r, c), vv in ent.items():
        M[r, c] = vv
    return len(span) - rank_mod(M, P0)


def colouring_exact(b, fields):
    """(m.sigma)^5 for the common axis m of the fields: commutes with H and anticommutes with Z_0, exactly."""
    m = next(v for v in fields.values() if v != (0, 0))
    ms = sp.Matrix([[0, m[0] - sp.I * m[1]], [m[0] + sp.I * m[1], 0]])
    W = ms
    for _ in range(N - 1):
        W = sp.kronecker_product(ms, W)
    H = sp.Matrix(hb_int(b).tolist())
    for j, (fx, fy) in fields.items():
        H += fx * SIG[('X', j)] + fy * SIG[('Y', j)]
    return is_zero(H * W - W * H) and is_zero(SIG[('Z', 0)] * W + W * SIG[('Z', 0)])


def h_full_mod(b, fields, p):
    ip = sqrt_m1(p)
    H = hb_int(b) % p
    for j, (fx, fy) in fields.items():
        for xx in range(D):
            y = xx ^ (1 << j)
            H[y, xx] = (H[y, xx] + fx + fy * (ip if not (xx >> j) & 1 else p - ip)) % p
    return H


def cyclic_dim(b, fields):
    best = 0
    for p in PRIMES:
        ip = sqrt_m1(p)
        Fp = np.zeros((D, D), dtype=np.int64)
        for j, (fx, fy) in fields.items():
            for xx in range(D):
                if not (xx >> j) & 1:
                    Fp[xx | (1 << j), xx] = (Fp[xx | (1 << j), xx] + fx + fy * ip) % p
        Z = np.diag([(p - 1 if xx & 1 else 1) for xx in range(D)]).astype(np.int64)
        gens = [hb_int(b) % p, Z, Fp]
        basis, piv, queue = [], [], [np.eye(D, dtype=np.int64)[0]]
        while queue:
            vv = np.array(np.array(queue.pop(), dtype=object) % p, dtype=np.int64)
            r = vv.copy()
            for bv, c in zip(basis, piv):
                if r[c]:
                    r = np.array((np.array(r, dtype=object) - int(r[c]) * np.array(bv, dtype=object)) % p, dtype=np.int64)
            nz = np.nonzero(r)[0]
            if len(nz) == 0:
                continue
            c = nz[0]
            r = np.array((np.array(r, dtype=object) * pow(int(r[c]), p - 2, p)) % p, dtype=np.int64)
            for t in range(len(basis)):
                if basis[t][c]:
                    basis[t] = np.array((np.array(basis[t], dtype=object) - int(basis[t][c]) * np.array(r, dtype=object)) % p,
                                        dtype=np.int64)
            basis.append(r); piv.append(c)
            for G in gens:
                queue.append(np.array((np.array(G, dtype=object).dot(np.array(vv, dtype=object))) % p, dtype=np.int64))
        best = max(best, len(basis))
    return best


def two_end(b, fields, p=P0):
    """(near, far) of F158: operators commuting with H that commute / anticommute with Z_0, kernels modulo p."""
    H = h_full_mod(b, fields, p)
    out = []
    for same in (True, False):
        cells = [(a, c) for a in range(D) for c in range(D) if ((a & 1) == (c & 1)) == same]
        M = np.zeros((D * D, len(cells)), dtype=np.int64)
        for k, (a, c) in enumerate(cells):
            for r in range(D):
                if H[r, a]:
                    M[r * D + c, k] = (M[r * D + c, k] + H[r, a]) % p
                if H[c, r]:
                    M[a * D + r, k] = (M[a * D + r, k] - H[c, r]) % p
        M = M[np.any(M, axis=1)]
        out.append(len(cells) - rank_mod(M, p))
    return tuple(out)


rows4 = [row5] + [plant(rng, False) for _ in range(4)]
ok_far, ok_I, ok_col = True, True, True
for b in rows4:
    (i_, j_), (k_, l_) = ds_pairings(b)[0]
    cases = [({i_: (1, 0), j_: (1, 0), k_: (2, 0), l_: (2, 0)}, (3, 0), 'equal within the pairs, one line'),
             ({i_: (1, 0), j_: (1, 0)}, (3, 0), 'equal on one pair, the other pair zero'),
             ({i_: (1, 0), j_: (1, 0), k_: (0, 1), l_: (0, 1)}, (2, 0), 'equal within the pairs, off one line'),
             ({i_: (3, 4), j_: (3, 4), k_: (1, -2), l_: (1, -2)}, (2, 0), 'equal within the pairs, off one line'),
             ({i_: (1, 0), j_: (2, 0), k_: (0, 1), l_: (0, 1)}, (0, 0), 'one pair unequal, off one line'),
             ({i_: (1, 0), j_: (-1, 0), k_: (0, 1), l_: (0, -1)}, (0, 0), 'opposite within the pairs'),
             ({i_: (1, 0), j_: (2, 0), k_: (3, 0), l_: (3, 0)}, (1, 0), 'one pair unequal, one line')]
    for fl, want, what in cases:
        got = (far_count(b, fl, 1), far_count(b, fl, 0))
        ok_far &= got == want
        if what.startswith('equal'):
            ok_I &= cyclic_dim(b, fl) < D
        if want[0] == 3 or want == (1, 0):
            ok_col &= colouring_exact(b, fl)
        if what.startswith('opposite'):
            ok_I &= cyclic_dim(b, fl) == D
check(ok_far, f"at {len(rows4)} double-singlet rows (the N = 5 row of the bullet 'At N = 4, a seat that is not blind' among them), seven field sets each: far "
              f"elements homogeneous under Z^5 (odd, even) = (3, 0) with the fields equal within the pairs on one line "
              f"(a pair of zeros included), (2, 0) off one line, and class A's (1, 0) or (0, 0) with a pair unequal or "
              f"opposite; the 3 and 2 reached by class A's colouring and part 2's two elements")
check(ok_col, "class A's colouring X^5 commutes with H and anticommutes with Z_0 exactly at every case on the x line "
              "above, the lower bound 1 of the counts 3 and (1, 0)")
check(ok_I, "there the cyclic space alg(H_b, Z_0, F+)|vac> is 32 with the fields opposite within the pairs, which "
            "certifies (I), and reads short of 32 with them equal (a read: part 2's subspace forces it over Q)")
fl3x = {1: (1, 0), 2: (1, 0), 3: (3, 0), 4: (3, 0)}
chiR = ss(1, 2, 3, 4)
PR = chiR * chiR.T / 4 + (SIG[('X', 0)] * chiR) * (SIG[('X', 0)] * chiR).T / 4
X5r = SIG[('X', 0)] * SIG[('X', 1)] * SIG[('X', 2)] * SIG[('X', 3)] * SIG[('X', 4)]
indep3 = sp.Matrix.hstack(*[sp.Matrix(M.reshape(D * D, 1)) for M in (X5r, SIG[('X', 0)] * PR, SIG[('Y', 0)] * PR)]).rank() == 3
check(colouring_exact(row5, fl3x) and indep3,
      "at that row with 3X on (3, 4), class A's colouring X^5 is exact, and X^5, X_0 (x) |s s><s s| and Y_0 (x) |s s><s s| "
      "are independent: the far end there is at least 3")
ne_xy = two_end(row5, {1: (1, 0), 2: (1, 0), 3: (0, 2), 4: (0, 2)})
ne_xx = two_end(row5, {1: (1, 0), 2: (1, 0), 3: (3, 0), 4: (3, 0)})
check(ne_xy == (3, 2) and ne_xx == (3, 3),
      f"the full two-end count at the N = 5 row of the bullet 'At N = 4, a seat that is not blind': near, far = {ne_xy} with X on (1, 2) and 2Y on (3, 4), "
      f"{ne_xx} with X on (1, 2) and 3X on (3, 4); near at least 3 (1 and the projectors onto |0>_0 s s and "
      f"|1>_0 s s), far at least 2 (part 2), so exactly: off one line the row does not pair, on one line it pairs")

bu = rows4[1]
(ui, uj), (uk, ul) = ds_pairings(bu)[0]
ne_u = two_end(bu, {ui: (1, 0), uj: (1, 0), uk: (0, 2), ul: (0, 2)})
check(fence3(bu) == 0 and ne_u == (3, 2),
      f"at a planted double singlet off the vacuum's energy, where (II) holds (fence at n = 3 empty), with X on one pair "
      f"and 2Y on the other: near, far = {ne_u}, the same bounds making it exact")

# ---- part 5 -------------------------------------------------------------------------------------------------------------
say()
say("Part 5. A row off the locus whose two-excitation sector is reducible")
Jr = (4, 0, 0, 6, 0, -3, -6, 4, 6, -6)
br = {ed: J for ed, J in zip(EDGES, Jr) if J}
L = sp.zeros(N, N)
for (i, j), J in br.items():
    L[i, j] -= J; L[j, i] -= J; L[i, i] += J; L[j, j] += J
vv = sp.zeros(N, 1); vv[0] = 1
Kr = sp.Matrix.hstack(*[L ** k * vv for k in range(N)]).rank()
H2 = sp.Matrix(hb_int(br)[np.ix_(S2, S2)].tolist())
Zs = sp.diag(*[(-1 if xx & 1 else 1) for xx in S2])
I10 = sp.eye(10)
rows = sp.Matrix.vstack(sp.kronecker_product(I10, H2) - sp.kronecker_product(H2.T, I10),
                        sp.kronecker_product(I10, Zs) - sp.kronecker_product(Zs.T, I10))
comm = 100 - DomainMatrix.from_Matrix(rows).convert_to(sp.QQ).rank()
Hs = sp.Matrix(hb_int(br).tolist()); E0 = sum(Jr)
cols = [xx for xx in S3 if xx & 1]
f3r = len(cols) - (Hs - E0 * sp.eye(D)).extract(list(range(D)), cols).rank()
check(Kr == N and not ds_pairings(br) and comm == 2 and f3r == 0 and joint_dim(br) == 0,
      f"J = {Jr} on the edges 01..34: 0 not blind (Krylov rank {Kr}, exact), no double singlet, the commutant of "
      f"alg(H_b, Z_0) on the two-excitation sector has dimension {comm} over Q (exact rank), no d-unexcited eigenvector in "
      f"that sector (a near miss for part 3), the fence at n = 3 is empty")
LET = [(0, 0), (1, 0), (-1, 0), (0, 1), (2, 0), (0, 2), (1, 1)]
okr, nsets = True, 0
rr = random.Random(5)
while nsets < 24:
    fl = {k: rr.choice(LET) for k in range(1, 5)}
    fl = {k: v for k, v in fl.items() if v != (0, 0)}
    if not fl:
        continue
    nsets += 1
    vs = list(fl.values())
    line = all(a_[0] * b2[1] - a_[1] * b2[0] == 0 for a_ in vs for b2 in vs)
    okr &= cyclic_dim(br, fl) == D and (far_count(br, fl, 1), far_count(br, fl, 0)) == ((1, 0) if line else (0, 0))
    if line:
        okr &= colouring_exact(br, fl)
check(okr, f"there, at {nsets} drawn field sets: (I) holds (cyclic space 32 by a full rank) and the far elements are "
           f"class A's alone, (1, 0) on one line, its colouring exact there, and (0, 0) off it")

say()
say(f"({time.time() - T0:.0f} s)")
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} FAILED: " + "; ".join(FAILS))
(Path(__file__).parent / "results" / "blind_seat_double_singlet_gate.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
sys.exit(1 if FAILS else 0)
