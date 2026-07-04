"""Seed-existence reduction: the nullity identity behind the real defective seed on the (1,2) block.

Companion verifier for experiments/F89_SEED_EXISTENCE_REDUCTION.md (the EXISTENCE side of the
codim-1 containment corollary's per-N input; complements the EXCLUSION shell census of
experiments/F89_MULTI_SECTOR_MONODROMY.md and the census table of experiments/F89_PATH_K_DIABOLIC.md).

The (1,2) block pencil is L(q) = A + q*C on coherences |a><b| (a a 2-excitation ket, b a
1-excitation bra) of the XY chain under uniform Z-dephasing:
    A = -2 * hamming(a,b)          the dephasing diagonal, rungs -2 (n_diff=1) and -6 (n_diff=3)
    q*C = -i q (H2 (x) I - I (x) H1)   the coherent XY hop; C is anti-Hermitian, C = i*K.
Normalization pin (q vs Q, the factor 2): H1/H2 are UNIT-hop hopping blocks (spec H1 =
{2 cos(k pi/(N+1))}), so this q coincides numerically with the carrier-clock Q and is TWICE the
octic-book q of the arc's census/seed loci (WeightCoherenceBlock hops -2iq per unit octic q):
L here at knob q = the census builder's block at coupling q/2. Every check below is a nullity or
multiset identity, invariant under jointly rescaling C and q, so no result depends on the choice.
Canonical statement: docs/GLOSSARY.md, "The coupling ratio q and Q".
The number of real eigenvalues at small q>0 is r0 = nullity(P2 C P2) + nullity(P6 C P6) (the
levels with a real, i.e. zero, first-order shift), and at large q it is r_inf = nullity(C)
(the asymptotically-real modes). A finite-q>0 real<->complex transition (hence, by the
discriminant-simple-zero lemma, a real defective seed) is forced whenever r0 > r_inf.

This script proves-by-computation the load-bearing facts (N = 3..13, BOTH parities):
  (F1)  r0 - r_inf = (N - 1) * [N odd]   exactly  (the seed-forcing surplus; zero at even N).
  (N2)  nullity(P2 C P2) = (N - 1) * [N odd]  via a decomposition into N-1 disjoint paths of N vertices.
  (FF)  nullity(C) = #{ (a<b, c) : lambda_a + lambda_b = lambda_c }, the free-fermion FUSION
        resonance count (lambda_k = 2 cos(k pi/(N+1)) the single-magnon energies; the two-magnon
        energies are the Slater sums lambda_a + lambda_b, cf. E_DE = eps_j + eps_k in
        F89_PATH_K_DIABOLIC.md), and n6 = nullity(P6 C P6) equals the SAME count.
  (N1P) the (N1') THEOREM (the -6 rung by ordering sectors; F89_SEED_EXISTENCE_REDUCTION.md Piece 3):
        K6 = -i*(P6 C P6) splits into three no-passing bra-rank sectors (bra left/middle/right of the
        ket pair), each gauge-equivalent via diag((-1)^{z_bra}) to MINUS the 3-magnon hopping H3, so
        spec(C6) = 3 copies of -i*(lambda_a+lambda_b+lambda_c) and n6 = 3*Z3 = the resonance count
        (Z3 = #zero-sum triples; the 3-to-1 chiral bijection has no degenerate D-cases, a theorem).
  (SI)  spectral inheritance as a per-block multiset THEOREM: spec(C2) = (N-1) x {i*lambda_k} and
        spec(C6) = 3 x {-i*(la+lb+lc)}, each an exact WITH-MULTIPLICITY sub-multiset of spec(C).
        (The joint union is NOT a partition of spec(C); the images overlap.)

Tolerance scope: the near-miss law min_{x!=y} |2*lambda_x + lambda_y| ~ 2*pi^3/(N+1)^3 (odd N;
itself a valid fusion near-resonance via c = N+1-x) makes the resonance tol 1e-7 safe only to
N ~ 850 and the SVD nullity floor 1e-9 safe only to N ~ 3959; at this script's N <= 13 the
margin is >= 7 orders.

Run:  python simulations/seed_existence_nullity_check.py     (asserts everything; prints a table)
"""
import numpy as np
from itertools import combinations


# ---------------------------------------------------------------- block builders
def _exc(N, k):
    return [tuple(c) for c in combinations(range(N), k)]


def _hop(N, k):
    """Nearest-neighbour XY hopping on the k-excitation sector (real symmetric)."""
    states = _exc(N, k)
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        occ = set(s)
        for i in range(N - 1):
            j = i + 1
            if (i in occ) ^ (j in occ):
                src = i if i in occ else j
                dst = j if i in occ else i
                H[idx[tuple(sorted((occ - {src}) | {dst}))], idx[s]] += 1.0
    return H


def build(N):
    """Return (A, C): the diagonal dephasing vector A and coherent coefficient C, L(q)=diag(A)+q*C."""
    H2, H1 = _hop(N, 2), _hop(N, 1)
    d1 = H1.shape[0]
    C = -1j * (np.kron(H2, np.eye(d1)) - np.kron(np.eye(H2.shape[0]), H1))
    A = np.array([-2.0 * len(set(a) ^ set(b)) for a in _exc(N, 2) for b in _exc(N, 1)])
    return A, C


def _nullity(M, tol=1e-9):
    if M.size == 0:
        return 0
    return int(np.sum(np.linalg.svd(M, compute_uv=False) < tol))


# ---------------------------------------------------------------- (N2) reduced -2 operator
def _minus2_states(N):
    return [(p, q) for p in range(N) for q in range(N) if p != q]  # p = shared site, q = ket-only


def _K_red(N):
    """The hand-derived -2-level reduced operator K_red (P2 C P2 = i*K_red)."""
    st = _minus2_states(N)
    idx = {s: i for i, s in enumerate(st)}
    K = np.zeros((len(st), len(st)))
    for (p, q) in st:
        for dq in (-1, 1):
            qp = q + dq
            if 0 <= qp < N and qp != p:            # q-hop: ket-only site cannot hop over p
                K[idx[(p, qp)], idx[(p, q)]] += -1.0
        if abs(p - q) == 1:                        # swap: bra hop onto the ket site
            K[idx[(q, p)], idx[(p, q)]] += +1.0
    return K, st


def _path_components(N):
    """Connected components of the K_red graph; returns list of (size, is_path)."""
    st = _minus2_states(N)
    adj = {s: set() for s in st}
    for (p, q) in st:
        for dq in (-1, 1):
            qp = q + dq
            if 0 <= qp < N and qp != p:
                adj[(p, q)].add((p, qp))
        if abs(p - q) == 1:
            adj[(p, q)].add((q, p))
    seen, comps = set(), []
    for s in st:
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        degs = [len(adj[v]) for v in comp]
        is_path = sorted(degs).count(1) == 2 and sorted(degs).count(2) == len(comp) - 2
        comps.append((len(comp), is_path))
    return comps


# ---------------------------------------------------------------- (FF) fusion-resonance count
def _resonance_count(N, tol=1e-7):
    lam = np.array([2 * np.cos(k * np.pi / (N + 1)) for k in range(1, N + 1)])
    two = [lam[a] + lam[b] for a in range(N) for b in range(a + 1, N)]
    return sum(1 for c in range(N) for e in two if abs(e - lam[c]) < tol)


def _spec(M):
    K = M / 1j
    return np.sort(np.linalg.eigvalsh((K + K.conj().T) / 2).real)


def _is_submultiset(sub, sup, tol=1e-9):
    """Exact with-multiplicity containment of two sorted real arrays (two-pointer greedy match)."""
    j = 0
    for x in sub:
        while j < len(sup) and sup[j] < x - tol:
            j += 1
        if j >= len(sup) or abs(sup[j] - x) > tol:
            return False
        j += 1
    return True


# ---------------------------------------------------------------- (N1P) ordering sectors on the -6 rung
def _check_n1prime(N, C, i6, basis, tol=1e-9):
    """The (N1') theorem, gate by gate. Returns Z3 (the zero-sum triple count)."""
    K6 = (C[np.ix_(i6, i6)] / 1j).real                      # C = i*K, so K6 = -i*C6, real symmetric
    b6 = [basis[i] for i in i6]

    def bra_rank(pair, c):                                  # bra left / middle / right of the ket pair
        lo, hi = min(pair), max(pair)
        return 0 if c < lo else (1 if c < hi else 2)

    sec = np.array([bra_rank(a, b[0]) for (a, b) in b6])
    for s in range(3):                                      # no-passing: cross-sector elements exactly 0
        for t in range(3):
            if s != t and K6[np.ix_(sec == s, sec == t)].size:
                assert np.abs(K6[np.ix_(sec == s, sec == t)]).max() == 0.0, \
                    f"(N1P) cross-sector element at N={N}"

    H3 = _hop(N, 3)
    b3 = _exc(N, 3)
    cfg = {b: i for i, b in enumerate(b3)}
    lam = np.array([2 * np.cos(k * np.pi / (N + 1)) for k in range(1, N + 1)])
    triple_sums = np.array([lam[a] + lam[b] + lam[c] for (a, b, c) in combinations(range(N), 3)])
    for s in range(3):                                      # gauge: U K_sec U = -H3 exactly
        ids = np.where(sec == s)[0]
        entries = [b6[i] for i in ids]
        perm = [cfg[tuple(sorted(set(e[0]) | {e[1][0]}))] for e in entries]
        assert len(set(perm)) == len(perm), f"(N1P) config map not injective at N={N}"
        P = np.zeros((len(b3), len(ids)))
        for col, row in enumerate(perm):
            P[row, col] = 1.0
        gauge = np.diag([(-1.0) ** e[1][0] for e in entries])
        assert np.abs(P @ (gauge @ K6[np.ix_(ids, ids)] @ gauge) @ P.T + H3).max() == 0.0, \
            f"(N1P) gauge identity U K_sec U != -H3 at N={N} sector {s}"

    spec6 = _spec(C[np.ix_(i6, i6)])                        # spec(C6) = 3 copies of -(la+lb+lc)
    assert np.abs(spec6 - np.sort(np.concatenate([-triple_sums] * 3))).max() < 1e-10, \
        f"(N1P) spec(C6) != 3 x (-triple sums) at N={N}"

    Z3 = int(np.sum(np.abs(triple_sums) < tol))             # D = 0: no 2*l_x + l_y = 0 with x != y
    D = sum(1 for x in range(N) for y in range(N)
            if x != y and abs(2 * lam[x] + lam[y]) < tol)
    assert D == 0, f"(N1P) degenerate resonance D={D} at N={N} (theorem says impossible)"
    return Z3


# ---------------------------------------------------------------- driver
def check(N):
    A, C = build(N)
    i2 = np.where(np.isclose(A, -2.0))[0]
    i6 = np.where(np.isclose(A, -6.0))[0]
    C2, C6 = C[np.ix_(i2, i2)], C[np.ix_(i6, i6)]
    n2, n6, nC = _nullity(C2), _nullity(C6), _nullity(C)
    r0, r_inf = n2 + n6, nC
    res = _resonance_count(N)

    # (F1) the surplus
    expected = (N - 1) * (N % 2)
    assert r0 - r_inf == expected, f"(F1) r0-r_inf = {r0 - r_inf} != {expected} at N={N}"
    # (N2) reduced operator matches, and the path decomposition
    Kred, st = _K_red(N)
    assert _nullity(Kred) == n2, f"(N2) reduced nullity mismatch at N={N}"
    assert _nullity(Kred) == expected, f"(N2) nullity(P2 C P2) = {_nullity(Kred)} != {expected} at N={N}"
    comps = _path_components(N)
    assert len(comps) == N - 1 and all(sz == N and pth for sz, pth in comps), \
        f"(N2) not (N-1) paths of N vertices at N={N}: {comps}"
    # (FF) fusion-resonance count = nullity(C) = n6
    assert nC == res == n6, f"(FF) nullity(C)={nC}, n6={n6}, resonances={res} disagree at N={N}"
    # (N1P) the ordering-sector theorem on the -6 rung
    basis = [(a, b) for a in _exc(N, 2) for b in _exc(N, 1)]
    Z3 = _check_n1prime(N, C, i6, basis)
    assert n6 == 3 * Z3, f"(N1P) n6={n6} != 3*Z3={3 * Z3} at N={N}"
    # (SI) spectral inheritance, per block, exact multisets (a theorem now; Piece 3 corollary)
    lam = np.array([2 * np.cos(k * np.pi / (N + 1)) for k in range(1, N + 1)])
    spC, spC2, spC6 = _spec(C), _spec(C2), _spec(C6)
    assert np.abs(spC2 - np.sort(np.tile(lam, N - 1))).max() < 1e-10, \
        f"(SI) spec(C2) != (N-1) x {{lambda_k}} at N={N}"
    assert _is_submultiset(spC2, spC) and _is_submultiset(spC6, spC), \
        f"(SI) per-block multiset inheritance fails at N={N}"
    return dict(dim=len(A), n2=n2, n6=n6, nC=nC, r0=r0, r_inf=r_inf, surplus=r0 - r_inf,
                res=res, Z3=Z3)


if __name__ == "__main__":
    print(f"{'N':>3} {'dim':>5} {'n2':>4} {'n6=3*Z3=res':>12} {'nullC':>6} "
          f"{'r0':>4} {'r_inf':>6} {'surplus':>8}")
    for N in range(3, 14):
        r = check(N)
        print(f"{N:>3} {r['dim']:>5} {r['n2']:>4} {r['n6']:>12} {r['nC']:>6} "
              f"{r['r0']:>4} {r['r_inf']:>6} {r['surplus']:>8}")
    print("\nAll assertions passed: (F1) surplus = (N-1)*[N odd], (N2) path decomposition, "
          "(FF) fusion-resonance count, (N1P) ordering-sector theorem, "
          "(SI) per-block multiset inheritance.")
