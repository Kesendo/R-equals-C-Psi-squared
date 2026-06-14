"""Probe (WIP): is the V-Effect (Pauli-string weight w=N/2 self-pair) the SAME object as the
{0,2}-coherence (n_diff) survivor, or a genuinely DIFFERENT decomposition?

The open question, typed verbatim in two claims (VacuumBlockReductionClaim, SurvivalIncompleteness-
MirrorClaim): the V-Effect (Pauli weight w=N/2) self-pair "co-locates" with the {0,2}-coherence
survivor, but they are DIFFERENT decompositions (Pauli weight vs n_diff bra-ket disagreement count);
their identity is OPEN, not claimed. This probe settles it on the ACTUAL survivor eigenvector.

EXACT single-site algebra (the reason the two gradings are not the same axis):
  a density component |i><j| factorizes per site:
    agreeing  site (i_s = j_s):  |0><0|=(I+Z)/2, |1><1|=(I-Z)/2     -> letters {I, Z}
    disagreeing site (i_s != j_s): |1><0|=(X-iY)/2, |0><1|=(X+iY)/2 -> letters {X, Y}
  => XY-weight (count of X,Y Pauli letters) == n_diff (count of disagreeing sites) EXACTLY.
     total Pauli weight  w == n_diff + (#Z letters on agreeing sites)  == light + Z-shadow.
  So "{0,2}-coherence" (n_diff in {0,2}) is a LOW-LIGHT statement; "V-Effect w=N/2" is a
  TOTAL-WEIGHT statement. They coincide on the survivor ONLY if its Z-shadow lands the total
  weight at N/2. That is an empirical fact about the eigenvector, measured below.

TEST: build the slowest non-kernel RIGHT eigenvector of each candidate low-light density block,
  project it onto the Pauli basis via the tensor-factorized single-site transform, and report
    - the n_diff (== XY-weight) histogram   [the known {0,2} signature]
    - the total Pauli-weight histogram      [does it peak at N/2 -> identity; else -> distinct]
    - the (XY-weight, Z-weight) split        [shows the Z-shadow that separates the two axes]
    - the survivor's squared mass at total weight w = N/2 (the V-Effect "self-pair" shell)
  Also disentangles (2,2) [the junction / ring survivor] vs (3,3)=(N/2,N/2) [the half-filling
  "V-Effect" block, where the chain survivor is pinned by the C# battery] -- same "half" word, two
  blocks. Carbon convention: XY (free fermions, no ZZ), J=1, gamma=1/Q.
"""
import numpy as np
from itertools import combinations

TOL = 1e-7


def popcount(x):
    return bin(x).count("1")


def bonds(N, topo):
    if topo == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topo == "ring":
        return [(i, i + 1) for i in range(N - 1)] + [(N - 1, 0)]
    if topo == "star":
        return [(0, k) for k in range(1, N)]
    raise ValueError(topo)


def basis(N, p):
    return [sum(1 << i for i in c) for c in combinations(range(N), p)]


def H_p(N, p, J, bnds):
    """p-excitation XY hopping Hamiltonian (amplitude J per hop)."""
    states = basis(N, p)
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), complex)
    for a, b in bnds:
        for s in states:
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += J
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += J
    return H, states


def block_L(N, pr, pc, J, g, bnds):
    """The (pr,pc) coherence-block Liouvillian and its row/col computational states."""
    Hr, sr = H_p(N, pr, J, bnds)
    Hc, sc = H_p(N, pc, J, bnds)
    L = -1j * (np.kron(Hr, np.eye(len(sc))) - np.kron(np.eye(len(sr)), Hc.T))
    deph = np.array([-2.0 * g * popcount(sr[a] ^ sc[b])
                     for a in range(len(sr)) for b in range(len(sc))])
    return L + np.diag(deph), sr, sc


def slowest_eigvec(L):
    """(lambda, right eigenvector v normalized) of the slowest non-kernel mode, or None."""
    w, V = np.linalg.eig(L)
    nz = [k for k in range(len(w)) if w[k].real < -TOL]
    if not nz:
        return None
    k = max(nz, key=lambda k: w[k].real)        # least-negative Re = longest-lived
    v = V[:, k]
    return w[k], v / np.linalg.norm(v)


def embed(v, sr, sc, N):
    """Embed the block eigenvector (index a*len(sc)+b -> |sr[a]><sc[b]|) into the full 2^N x 2^N rho."""
    d = 1 << N
    M = np.zeros((d, d), complex)
    nc = len(sc)
    for a in range(len(sr)):
        for b in range(nc):
            M[sr[a], sc[b]] = v[a * nc + b]
    return M


# single-site change of basis |i><j| (row=ket i, col=bra j, order m00,m01,m10,m11) -> (cI,cX,cY,cZ)
#   cI=(m00+m11)/2  cX=(m01+m10)/2  cY=i(m01-m10)/2  cZ=(m00-m11)/2
_B = np.array([[0.5, 0.0, 0.0, 0.5],
               [0.0, 0.5, 0.5, 0.0],
               [0.0, 0.5j, -0.5j, 0.0],
               [0.5, 0.0, 0.0, -0.5]], complex)


def pauli_coeffs(M, N):
    """Full Pauli-coefficient tensor c[a_0,...,a_{N-1}], letters 0=I,1=X,2=Y,3=Z (M = sum c_P P)."""
    T = M.reshape((2,) * N + (2,) * N)
    perm = [x for s in range(N) for x in (s, N + s)]      # interleave ket_s, bra_s
    T = np.transpose(T, perm).reshape((4,) * N)
    for s in range(N):
        T = np.tensordot(_B, T, axes=([1], [s]))
        T = np.moveaxis(T, 0, s)
    return T


def _letter_weights(N):
    """xy[idx], z[idx] arrays over the flattened (4,)*N index: #{X,Y} letters and #Z letters."""
    grids = np.meshgrid(*([np.arange(4)] * N), indexing="ij")
    xy = np.zeros((4,) * N, int)
    z = np.zeros((4,) * N, int)
    for g in grids:
        xy += ((g == 1) | (g == 2)).astype(int)
        z += (g == 3).astype(int)
    return xy.ravel(), z.ravel()


def signatures(M, N):
    """n_diff hist (from M), Pauli total-weight hist, (xy,z) joint, all normalized to sum 1."""
    # n_diff histogram from the |i><j| support (the NDiffHistogram convention)
    d = 1 << N
    ndh = {}
    nz = np.argwhere(np.abs(M) > 1e-12)
    for i, j in nz:
        nd = popcount(int(i) ^ int(j))
        ndh[nd] = ndh.get(nd, 0.0) + abs(M[i, j]) ** 2
    s = sum(ndh.values()) or 1.0
    ndh = {k: v / s for k, v in sorted(ndh.items())}

    c = pauli_coeffs(M, N).ravel()
    p = np.abs(c) ** 2
    xy, z = _letter_weights(N)
    w = xy + z
    tot = p.sum() or 1.0
    wh = {int(k): float(p[w == k].sum() / tot) for k in range(N + 1) if p[w == k].sum() > 1e-12}
    # XY-weight marginal of the Pauli power (must equal n_diff hist -> transform validation)
    xyh = {int(k): float(p[xy == k].sum() / tot) for k in range(N + 1) if p[xy == k].sum() > 1e-12}
    return ndh, wh, xyh


def global_survivor(N, J, g, bnds):
    """Slowest over the (p,p) p=1..N-1 blocks + the (0,1) band edge; returns (sector, lambda, M)."""
    cands = [(p, p) for p in range(1, N)] + [(0, 1)]
    best = None
    for pr, pc in cands:
        L, sr, sc = block_L(N, pr, pc, J, g, bnds)
        r = slowest_eigvec(L)
        if r is None:
            continue
        lam, v = r
        if best is None or lam.real > best[1].real:
            best = ((pr, pc), lam, embed(v, sr, sc, N))
    return best


def analyse_block(N, pr, pc, J, g, bnds, label):
    L, sr, sc = block_L(N, pr, pc, J, g, bnds)
    r = slowest_eigvec(L)
    if r is None:
        print(f"   ({pr},{pc}) {label}: all-kernel")
        return None
    lam, v = r
    M = embed(v, sr, sc, N)
    ndh, wh, xyh = signatures(M, N)
    peak_w = max(wh, key=wh.get)
    half = N // 2
    frac_half = wh.get(half, 0.0)
    # transform self-check: XY-weight marginal == n_diff histogram
    assert all(abs(xyh.get(k, 0.0) - ndh.get(k, 0.0)) < 1e-9 for k in set(xyh) | set(ndh)), \
        f"XY-marginal != n_diff hist in ({pr},{pc})"
    print(f"   ({pr},{pc}) {label}:  -Re(lam)={-lam.real:.4f}")
    print(f"        n_diff (=XY-wt) : {{{', '.join(f'{k}:{v:.3f}' for k, v in ndh.items())}}}")
    print(f"        Pauli total-wt  : {{{', '.join(f'{k}:{v:.3f}' for k, v in wh.items())}}}  peak w={peak_w}")
    print(f"        mass at w=N/2={half}: {frac_half:.3f}   (the V-Effect 'self-pair' shell)")
    return dict(sector=(pr, pc), ndh=ndh, wh=wh, peak_w=peak_w, frac_half=frac_half)


def _validate_sector_method():
    """Sanity: the (p,p)/(0,1) sector survivor reproduces the full 4^N slowest at N=4."""
    def full_slowest(N, J, g, bnds):
        I2 = np.eye(2, dtype=complex)
        X = np.array([[0, 1], [1, 0]], complex)
        Y = np.array([[0, -1j], [1j, 0]])
        Z = np.diag([1, -1]).astype(complex)

        def site(op, l):
            m = np.array([[1.0 + 0j]])
            for k in range(N):
                m = np.kron(m, op if k == l else I2)
            return m
        d = 2 ** N
        H = sum((J / 2) * (site(X, a) @ site(X, b) + site(Y, a) @ site(Y, b)) for a, b in bnds)
        Id = np.eye(d)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        for l in range(N):
            Zl = site(Z, l)
            L += g * (np.kron(Zl, Zl) - np.kron(Id, Id))
        ev = np.linalg.eigvals(L)
        return ev[ev.real < -TOL].real.max()
    for topo in ("chain", "ring", "star"):
        bnds = bonds(4, topo)
        g = 1.0 / 1.5
        sec_re = global_survivor(4, 1.0, g, bnds)[1].real
        full_re = full_slowest(4, 1.0, g, bnds)
        assert abs(sec_re - full_re) < 1e-4, f"{topo}: sector {sec_re} != full {full_re}"
    print("[validate] sector survivor == full 4^N slowest at N=4 (chain/ring/star) - method sound\n")


if __name__ == "__main__":
    _validate_sector_method()
    Q = 1.5
    g = 1.0 / Q

    print(f"Carbon XY, J=1, gamma=1/Q={g:.4f} (Q={Q}, below the handover)\n")
    print("== GLOBAL survivor sector per topology (where the longest-lived mode actually lives) ==")
    for N in (4, 6):
        for topo in ("chain", "ring"):
            sec, lam, M = global_survivor(N, 1.0, g, bonds(N, topo))
            ndh, wh, xyh = signatures(M, N)
            print(f"   N={N:>2} {topo:>5}: survivor sector {sec}, -Re(lam)={-lam.real:.4f}, "
                  f"n_diff peak={max(ndh, key=ndh.get)}, Pauli-wt peak={max(wh, key=wh.get)}, "
                  f"mass@w=N/2={wh.get(N // 2, 0.0):.3f}")

    print("\n== N=6 block-by-block: (2,2) [junction/ring] vs (3,3)=(N/2,N/2) [half-filling 'V-Effect'] ==")
    for topo in ("chain", "ring"):
        print(f" {topo}:")
        analyse_block(6, 2, 2, 1.0, g, bonds(6, topo), "low-light interior (the {0,2} junction)")
        analyse_block(6, 3, 3, 1.0, g, bonds(6, topo), "half-filling centre (the V-Effect block)")

    print("\n== per-sector slowest rate -Re(lam) (full precision, to expose ties / true winner) ==")
    for N in (4, 6):
        for topo in ("chain", "ring", "star"):
            rates = {}
            for p in range(1, N):
                r = slowest_eigvec(block_L(N, p, p, 1.0, g, bonds(N, topo))[0])
                rates[(p, p)] = -r[0].real if r else float("nan")
            r01 = slowest_eigvec(block_L(N, 0, 1, 1.0, g, bonds(N, topo))[0])
            rates[(0, 1)] = -r01[0].real if r01 else float("nan")
            win = min((k for k in rates if rates[k] == rates[k]), key=rates.get)
            print(f"   N={N} {topo:>5}: " + "  ".join(f"{k}:{v:.6f}" for k, v in rates.items())
                  + f"   -> winner {win}")

    print("\n== VERDICT INPUTS ==")
    print("If the survivor's Pauli total-weight mass concentrates at w=N/2 -> identity HOLDS.")
    print("If it spreads / peaks away from N/2 while n_diff stays in {0,2} -> the two are DISTINCT")
    print("decompositions (low-light survivor + Z-shadow, NOT a half-Pauli-weight object).")
