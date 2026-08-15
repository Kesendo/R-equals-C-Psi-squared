"""THE MIXED-SPACE REFLECTION GATE: the N = 8, Delta = 0 mixed-space law of
F154, derived in docs/proofs/PROOF_MIXED_SPACE_REFLECTION_LAW.md, verified
from below.

The proof's chain: at Delta = 0 the ad_H eigenbasis of a coherence block is
the Slater pair basis (F2b sine modes, e_k = cos(k pi / M), M = N + 1); each
pair state has definite R parity eps(A) eps(B) times a sector constant
(R psi_k = (-1)^(k+1) psi_k); N_l is (one-body) x (one-body), so inside one
eigenspace it connects only pairs solving e_i + e_n = e_m + e_j; and at
M = 9 the PAIR-SUM LEMMA says the only solutions are the trivial and the
chiral ones, both parity-EVEN. Hence comp(N_l) is parity-block-diagonal on
every eigenspace, for every site l separately, which is strictly stronger
than the measured C_l = 0 of endpoint_density_gate.py section (B).

Four sections. (1) FORCED combs, exact: the pair-sum census of the comb
{1..M-1} over GF(p) for two primes p = 1 mod 2M (distinctness mod one prime
proves exact distinctness, the LevelCollisionCensus soundness direction);
gate: at M = 6, 9, 10, 11, 13, 14, 16, 25, 27 the ONLY colliding group is
the chiral-zero group ({x, M-x}, plus the doubled zero mode {M/2, M/2} at
even M), covering the proof's section 7 claims (M = 9 the flagship, odd
prime powers, primes, 2p, 2^a, and M = 6 as F146's recorded divisible
exception). (2) DOOR combs, exact: at M = 12, 15, 18, 21 a non-chiral
colliding group exists, a named parity-ODD instance is verified mod both
primes ({1,9}~{5,6} at 12, {1,5}~{3,4} at 15, {2,14}~{8,9} at 18,
{3,9}~{6,7} at 21), and the general 6|M door family
e_(k-M/3) + e_(k+M/3) = e_k + e_(M/2) is verified for EVERY admissible k at
M = 12, 18, 24, 30, with the parity factor (-1)^(3k+M/2) alternating (at
M = 6 the family is empty, the exception). (3) The sharpened prediction at
N = 8, Delta = 0 on five blocks, the endpoint census trio (1,3), (1,4),
(2,4) plus (2,3) and the population-carrying (4,4): per-l cross-parity
spectral norms of comp(N_l) on every parity-mixed eigenspace, plus a random
combination Sum c_l N_l, gated at 64 * eps*||A||/gap each (the endpoint
gate's error model); contrast: a generic MULTI-body diagonal on the same
parity halves is gated from below (the theorem is a one-body statement).
The eigenvalue grouping gtol = 1e-9 sits far below the difference
spectrum's smallest true gap at N = 8, Delta = 0 (measured 4.5e-2 on these
blocks; grouping cannot merge distinct eigenspaces). (4) The frontier
census in MODE SPACE, exact: bucketing Slater pairs by E_A - E_B over two
primes with the parity eps(A) eps(B), the mixed-bucket counts on (1,3),
(1,4), (2,4) are 0 at M = 6, 7, 8, 10, 11, exactly (20, 21, 50) at M = 9
(reproducing the endpoint census from arithmetic alone), and nonzero at
M = 12: mixed collisions fire exactly at 3|M with M > 6 on the tested
frontier (at M = 6 the only collisions are chiral, hence parity-even, the
same exception F146 records at rung 2).

Writes simulations/results/mixed_space_reflection_gate.txt; last line VERDICT.
"""

import numpy as np
from itertools import combinations, combinations_with_replacement

OUT = "simulations/results/mixed_space_reflection_gate.txt"
lines = []
fails = []
eps = np.finfo(float).eps


def log(s=""):
    lines.append(s)
    print(s)


def popcount(x):
    return bin(x).count("1")


def revbits(a, n):
    r = 0
    for l in range(n):
        if (a >> l) & 1:
            r |= 1 << (n - 1 - l)
    return r


# ------------------------------------------------- exact comb arithmetic
def is_prime(p):
    if p < 2:
        return False
    d = 2
    while d * d <= p:
        if p % d == 0:
            return False
        d += 1
    return True


def two_primes(m, lo=10 ** 6):
    """The two smallest primes p = 1 mod m with p > lo."""
    out, k = [], lo // m + 1
    while len(out) < 2:
        p = k * m + 1
        if is_prime(p):
            out.append(p)
        k += 1
    return tuple(out)


def find_zeta(p, order):
    """An element of exact multiplicative order `order` in GF(p)."""
    small = (2, 3, 5, 7, 11, 13)
    rest = order
    for q in small:
        while rest % q == 0:
            rest //= q
    assert rest == 1, f"order {order} has a prime factor outside {small}"
    for g in range(2, p):
        x, ok = pow(g, (p - 1) // order, p), True
        for q in small:
            if order % q == 0 and pow(x, order // q, p) == 1:
                ok = False
                break
        if ok:
            return x
    raise RuntimeError("no element of the required order")


def comb_mod(M, p):
    """lambda_k = zeta^k + zeta^(-k) mod p, zeta of exact order 2M."""
    z = find_zeta(p, 2 * M)
    return {k: (pow(z, k, p) + pow(z, 2 * M - k, p)) % p for k in range(1, M)}


def pair_sum_groups(M, p):
    """Groups of pair-sum 2-multisets colliding mod p (>= the exact groups;
    distinctness mod p PROVES exact distinctness). None on a mod-p
    single-particle degeneracy (never exact; retry another prime)."""
    lam = comb_mod(M, p)
    if len(set(lam.values())) != M - 1:
        return None
    sums = {}
    for x, y in combinations_with_replacement(range(1, M), 2):
        sums.setdefault((lam[x] + lam[y]) % p, []).append((x, y))
    return sorted(tuple(sorted(v)) for v in sums.values() if len(v) > 1)


def chiral_group(M):
    """The chiral-zero group: {x, M-x}, plus {M/2, M/2} at even M."""
    g = [(x, M - x) for x in range(1, (M - 1) // 2 + 1)]
    if M % 2 == 0:
        g.append((M // 2, M // 2))
    return tuple(sorted(g))


log("=" * 78)
log("THE MIXED-SPACE REFLECTION GATE: forced combs, doors, prediction, frontier")
log("=" * 78)

# --------------------------------------------------- (1) the forced combs
log()
log("(1) FORCED COMBS, EXACT. Pair-sum census over GF(p), two primes")
log("    p = 1 mod 2M each: the only colliding group is the chiral-zero")
log("    group. M = 9 is the theorem's comb; 6 is F146's recorded divisible")
log("    exception; 10, 14 = 2p; 16 = 2^4; 11, 13 prime; 25, 27 the odd")
log("    prime powers beyond 9 (proof section 7).")
for M in (6, 9, 10, 11, 13, 14, 16, 25, 27):
    want = [chiral_group(M)]
    for p in two_primes(2 * M):
        groups = pair_sum_groups(M, p)
        ok = groups == want
        if not ok:
            fails.append(f"forced comb M={M} at p={p}")
        log(f"    M={M:2d} p={p}: colliding groups == chiral zeros: {ok}")

# ----------------------------------------------------- (2) the door combs
log()
log("(2) DOOR COMBS, EXACT. At M = 12, 15, 18, 21 (the rung-2 doors of")
log("    PROOF_SCALAR_COUNT section 7) a non-chiral colliding group exists")
log("    and a named parity-ODD instance holds mod both primes: the forcing")
log("    of the law ends there. The 6|M family")
log("    e_(k-M/3) + e_(k+M/3) = e_k + e_(M/2) is verified for EVERY")
log("    admissible k (M/3 < k < 2M/3, k != M/2), parity (-1)^(3k+M/2); at")
log("    M = 6 the admissible range is empty, the exception's shape.")
door_instances = {12: ((1, 9), (5, 6)), 15: ((1, 5), (3, 4)),
                  18: ((2, 14), (8, 9)), 21: ((3, 9), (6, 7))}
for M, ((a, b), (c, d)) in sorted(door_instances.items()):
    par = (-1) ** (a + b + c + d)
    ok_par = par == -1
    if not ok_par:
        fails.append(f"door instance M={M} not parity-odd")
    for p in two_primes(2 * M):
        lam = comb_mod(M, p)
        groups = pair_sum_groups(M, p)
        ok_eq = (lam[a] + lam[b]) % p == (lam[c] + lam[d]) % p
        ok_door = groups is not None and any(g != chiral_group(M) for g in groups)
        if not ok_eq:
            fails.append(f"door instance M={M} fails at p={p}")
        if not ok_door:
            fails.append(f"no non-chiral group at door M={M}, p={p}")
        log(f"    M={M:2d} p={p}: e_{a}+e_{b} == e_{c}+e_{d}: {ok_eq}, "
            f"non-chiral group present: {ok_door}, parity {par:+d}")
assert not [k for k in range(6 // 3 + 1, 2 * 6 // 3) if k != 3], \
    "M = 6 family unexpectedly nonempty"
log("    M= 6: admissible-k range of the family is empty (the exception)")
for M in (12, 18, 24, 30):
    p1, p2 = two_primes(2 * M)
    n_odd = n_all = 0
    for p in (p1, p2):
        lam = comb_mod(M, p)
        for k in range(M // 3 + 1, 2 * M // 3):
            if k == M // 2:
                continue
            if (lam[k - M // 3] + lam[k + M // 3]) % p != (lam[k] + lam[M // 2]) % p:
                fails.append(f"6|M family fails at M={M}, k={k}, p={p}")
    for k in range(M // 3 + 1, 2 * M // 3):
        if k == M // 2:
            continue
        n_all += 1
        if (-1) ** (3 * k + M // 2) == -1:
            n_odd += 1
    ok = n_odd >= 1
    if not ok:
        fails.append(f"6|M family has no parity-odd k at M={M}")
    log(f"    M={M:2d}: family holds at every admissible k (both primes); "
        f"parity-odd for {n_odd} of {n_all} k  {'ok' if ok else 'FAIL'}")

# ------------------------------------------- (3) the sharpened prediction
log()
log("(3) SHARPENED PREDICTION AT N = 8, Delta = 0: comp(N_l) is parity-")
log("    block-diagonal on EVERY multi-dim eigenspace, per site l, gated at")
log("    64 * eps*||A||/gap; ditto a random combination Sum c_l N_l; the")
log("    generic MULTI-body diagonal contrast is gated from below (> 0.01).")
log("    Blocks: the endpoint census trio (1,3), (1,4), (2,4) plus (2,3)")
log("    and the population-carrying (4,4). Grouping note: the smallest")
log("    true inter-eigenvalue gap on these blocks is printed and gated")
log("    >> gtol = 1e-9 (grouping cannot merge distinct eigenspaces).")


def build_H(n, J, D):
    Jb = np.full(n - 1, float(J)) if np.isscalar(J) else np.asarray(J, float)
    d = 1 << n
    H = np.zeros((d, d))
    for a in range(d):
        for l in range(n - 1):
            za = 1 - 2 * ((a >> l) & 1)
            zb = 1 - 2 * ((a >> (l + 1)) & 1)
            H[a, a] += Jb[l] * D * za * zb
            if ((a >> l) & 1) != ((a >> (l + 1)) & 1):
                H[a, a ^ (1 << l) ^ (1 << (l + 1))] += 2 * Jb[l]
    return H


def block_A(n, p, q, D, J=1.0):
    H = build_H(n, J, D)
    kets = [a for a in range(1 << n) if popcount(a) == p]
    bras = [b for b in range(1 << n) if popcount(b) == q]
    cells = [(a, b) for a in kets for b in bras]
    idx = {c: i for i, c in enumerate(cells)}
    m = len(cells)
    A = np.zeros((m, m), dtype=complex)
    for (a, b), i in idx.items():
        for c in kets:
            if H[a, c] != 0:
                A[i, idx[(c, b)]] += -1j * H[a, c]
        for d_ in bras:
            if H[b, d_] != 0:
                A[i, idx[(a, d_)]] += 1j * H[b, d_]
    return A, cells, idx


def site_indicators(n, cells):
    return [np.array([1.0 if (a ^ b) & (1 << l) else 0.0 for (a, b) in cells])
            for l in range(n)]


def reflection_matrix(n, cells, idx):
    m = len(cells)
    R = np.zeros((m, m))
    for (a, b), i in idx.items():
        R[idx[(revbits(a, n), revbits(b, n))], i] = 1.0
    return R


def eigenspaces(A, gtol=1e-9):
    w, V = np.linalg.eigh(1j * A)
    nrm = np.abs(w).max()
    out = []
    start = 0
    for i in range(1, len(w) + 1):
        if i == len(w) or abs(w[i] - w[start]) > gtol:
            gap = np.inf
            if start > 0:
                gap = min(gap, w[start] - w[start - 1])
            if i < len(w):
                gap = min(gap, w[i] - w[i - 1])
            out.append((w[start], V[:, start:i], gap))
            start = i
    return out, nrm


rng = np.random.default_rng(20260815)
n = 8
trio_mixed = 0
min_gap = np.inf
for (p, q) in [(1, 3), (1, 4), (2, 4), (2, 3), (4, 4)]:
    A, cells, idx = block_A(n, p, q, 0.0)
    Nl = site_indicators(n, cells)
    R = reflection_matrix(n, cells, idx)
    spaces, nrm = eigenspaces(A)
    min_gap = min(min_gap, min(g for (_, _, g) in spaces if np.isfinite(g)))
    c = rng.uniform(-1, 1, n)
    Ncomb = sum(c[l] * Nl[l] for l in range(n))
    worst = worst_comb = contrast = 0.0
    n_multi = n_mixed = 0
    for wv, P, gap in spaces:
        if P.shape[1] < 2:
            continue
        n_multi += 1
        RW = P.conj().T @ R @ P
        inv_res = float(np.max(np.abs(RW @ RW - np.eye(P.shape[1]))))
        if inv_res > 1e-9:
            fails.append(f"R not an involution on an eigenspace of ({p},{q})")
        ev, U = np.linalg.eigh((RW + RW.conj().T) / 2)
        plus, minus = U[:, ev > 0], U[:, ev < 0]
        if plus.shape[1] == 0 or minus.shape[1] == 0:
            continue  # scalar eigenspace: no cross block exists
        n_mixed += 1
        Qp, Qm = P @ plus, P @ minus
        floor = 64 * eps * nrm / gap
        r = max(float(np.linalg.norm(
            Qm.conj().T @ (Nl[l][:, None] * Qp), 2)) for l in range(n))
        worst = max(worst, r / floor)
        rc = float(np.linalg.norm(Qm.conj().T @ (Ncomb[:, None] * Qp), 2))
        worst_comb = max(worst_comb, rc / floor)
        d = rng.uniform(-1, 1, len(cells))
        contrast = max(contrast, float(np.linalg.norm(
            Qm.conj().T @ (d[:, None] * Qp), 2)))
    if (p, q) in [(1, 3), (1, 4), (2, 4)]:
        trio_mixed += n_mixed
    ok = worst <= 1.0 and worst_comb <= 1.0 and contrast > 0.01
    log(f"    N=8 ({p},{q}): {n_mixed} mixed of {n_multi} multi-dim "
        f"eigenspaces; worst per-l cross/floor {worst:.3f}, random-"
        f"combination cross/floor {worst_comb:.3f}, generic multi-body "
        f"contrast {contrast:.3f}  {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"sharpened prediction N=8 ({p},{q})")
log(f"    mixed eigenspaces on the endpoint trio: {trio_mixed} "
    f"(the endpoint gate's 91); smallest true gap {min_gap:.2e} "
    f"(>> gtol = 1e-9)")
if trio_mixed != 91:
    fails.append("mixed-space count changed on the endpoint trio")
if not min_gap > 1e-3:
    fails.append("difference-spectrum gap approaches gtol")

# ------------------------------------- (4) the frontier census, mode space
log()
log("(4) FRONTIER CENSUS IN MODE SPACE, EXACT: Slater pairs bucketed by")
log("    E_A - E_B over two primes, a bucket MIXED when it carries both")
log("    parities eps(A) eps(B). Gate: 0 mixed buckets at M = 6, 7, 8, 10,")
log("    11; exactly (20, 21, 50) on (1,3)/(1,4)/(2,4) at M = 9 (the")
log("    endpoint census, reproduced from arithmetic alone); nonzero at")
log("    M = 12: on this frontier, mixed collisions fire exactly at 3|M")
log("    with M > 6 (M = 6 the recorded exception, chiral-only).")


def mode_census(M, p_exc, q_exc):
    """Mixed-bucket count on block (p_exc, q_exc) of the N = M-1 chain."""
    primes = two_primes(2 * M)
    lams = [comb_mod(M, p) for p in primes]
    modes = list(range(1, M))
    buckets = {}
    for A in combinations(modes, p_exc):
        eA = tuple(sum(lam[k] for k in A) % p for lam, p in zip(lams, primes))
        pA = (-1) ** sum(k + 1 for k in A)
        for B in combinations(modes, q_exc):
            key = tuple((eA[i] - sum(lam[k] for k in B)) % p
                        for i, (lam, p) in enumerate(zip(lams, primes)))
            sgn = pA * (-1) ** sum(k + 1 for k in B)
            e = buckets.setdefault(key, [0, 0])
            e[0 if sgn > 0 else 1] += 1
    return sum(1 for v in buckets.values() if v[0] > 0 and v[1] > 0)


for M in (6, 7, 8, 9, 10, 11, 12):
    counts = tuple(mode_census(M, p_exc, q_exc)
                   for (p_exc, q_exc) in [(1, 3), (1, 4), (2, 4)])
    if M == 9:
        ok = counts == (20, 21, 50)
        note = "== endpoint census (20, 21, 50)"
    elif M == 12:
        ok = all(c > 0 for c in counts)
        note = "nonzero (the door is populated)"
    else:
        ok = counts == (0, 0, 0)
        note = "zero"
    log(f"    M={M:2d} (N={M - 1:2d}): mixed buckets (1,3)/(1,4)/(2,4) = "
        f"{counts}  {note}: {'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append(f"frontier census at M={M}")

log()
log("=" * 78)
if fails:
    log("VERDICT: FAIL: " + "; ".join(fails))
else:
    log("VERDICT: all checks pass")
with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
