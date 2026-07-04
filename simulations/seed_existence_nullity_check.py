"""Seed-existence reduction: the nullity identity behind the real defective seed on the (1,2) block.

Companion verifier for experiments/F89_SEED_EXISTENCE_REDUCTION.md (the EXISTENCE side of the
codim-1 containment corollary's per-N input; complements the EXCLUSION shell census of
experiments/F89_MULTI_SECTOR_MONODROMY.md and the census table of experiments/F89_PATH_K_DIABOLIC.md).

The (1,2) block pencil is L(q) = A + q*C on coherences |a><b| (a a 2-excitation ket, b a
1-excitation bra) of the XY chain under uniform Z-dephasing:
    A = -2 * hamming(a,b)          the dephasing diagonal, rungs -2 (n_diff=1) and -6 (n_diff=3)
    q*C = -i q (H2 (x) I - I (x) H1)   the coherent XY hop; C is anti-Hermitian, C = i*K.
The number of real eigenvalues at small q>0 is r0 = nullity(P2 C P2) + nullity(P6 C P6) (the
levels with a real, i.e. zero, first-order shift), and at large q it is r_inf = nullity(C)
(the asymptotically-real modes). A finite-q>0 real<->complex transition (hence, by the
discriminant-simple-zero lemma, a real defective seed) is forced whenever r0 > r_inf.

This script proves-by-computation the four load-bearing facts:
  (F1)  r0 - r_inf = N - 1   exactly, for every odd N  (the seed-forcing surplus).
  (N2)  nullity(P2 C P2) = N - 1 (odd N)  via a decomposition into N-1 disjoint paths of N vertices.
  (FF)  nullity(C) = #{ (a<b, c) : lambda_a + lambda_b = lambda_c }, the free-fermion FUSION
        resonance count (lambda_k = 2 cos(k pi/(N+1)) the single-magnon energies; the two-magnon
        energies are the Slater sums lambda_a + lambda_b, cf. E_DE = eps_j + eps_k in
        F89_PATH_K_DIABOLIC.md), and n6 = nullity(P6 C P6) equals the SAME count.
  (SI)  spectral inheritance: spec(P2 C P2) and spec(P6 C P6) are exact sub-multisets of spec(C).

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
    assert r0 - r_inf == N - 1, f"(F1) r0-r_inf = {r0 - r_inf} != N-1 at N={N}"
    # (N2) reduced operator matches, and the path decomposition
    Kred, st = _K_red(N)
    assert _nullity(Kred) == n2, f"(N2) reduced nullity mismatch at N={N}"
    assert _nullity(Kred) == N - 1, f"(N2) nullity(P2 C P2) = {_nullity(Kred)} != N-1 at N={N}"
    comps = _path_components(N)
    assert len(comps) == N - 1 and all(sz == N and pth for sz, pth in comps), \
        f"(N2) not (N-1) paths of N vertices at N={N}: {comps}"
    # (FF) fusion-resonance count = nullity(C) = n6
    assert nC == res == n6, f"(FF) nullity(C)={nC}, n6={n6}, resonances={res} disagree at N={N}"
    # (SI) spectral inheritance
    spC, spC2, spC6 = _spec(C), _spec(C2), _spec(C6)
    inC = lambda sp: all(any(abs(x - y) < 1e-9 for y in spC) for x in sp)
    assert inC(spC2) and inC(spC6), f"(SI) spectral inheritance fails at N={N}"
    return dict(dim=len(A), n2=n2, n6=n6, nC=nC, r0=r0, r_inf=r_inf, surplus=r0 - r_inf, res=res)


if __name__ == "__main__":
    print(f"{'N':>3} {'dim':>5} {'n2=N-1':>7} {'n6=res':>7} {'nullC':>6} "
          f"{'r0':>4} {'r_inf':>6} {'surplus=N-1':>12}")
    for N in (3, 5, 7, 9, 11, 13):
        r = check(N)
        print(f"{N:>3} {r['dim']:>5} {r['n2']:>7} {r['n6']:>7} {r['nC']:>6} "
              f"{r['r0']:>4} {r['r_inf']:>6} {r['surplus']:>12}")
    print("\nAll assertions passed: (F1) surplus = N-1, (N2) path decomposition, "
          "(FF) fusion-resonance count, (SI) spectral inheritance.")
