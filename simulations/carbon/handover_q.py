"""The handover Q: where the incompleteness survivor meets the F50 floor (a closed condition,
a known chain solution, a distinct ring one). Resolves the open loose end of the
survival_incompleteness_mirror + clock_hand_ladder arcs and the F2b corollary "Open remainder".

THE QUESTION (arc, verbatim): "whether the handover Q (where the incomplete stops winning and
the odd band edge takes over) has a closed form."

THE ANSWER, in three parts (all self-validated below; XY/free-fermion convention to match
incompleteness_survivor.py + benzene_two_clocks.py: H=(J/2)sum(XX+YY), dephasing g sum_l D[Z_l],
Re(lambda) = -2g <n_XY>, Q = J/g; the handover is a SPECTRAL property, state-independent, and
depends only on Q=J/g, NOT on J and g separately).

1. THE CONDITION (exact, universal, F50-grounded). The diagonal (p,p) "incompleteness survivor"
   decays at -2g<n_XY> with FRACTIONAL darkness <n_XY> < 1, so it out-survives the bare band edge.
   As Q rises <n_XY> brightens; the HANDOVER is where it reaches the F50-pinned OFF-diagonal floor
   <n_XY> = 1 (the (0,1) band edge / Uhr 1, Re = -2g EXACTLY, F50 weight-1 degeneracy). Below: the
   dressed interior mode is darker (the incomplete survives). Above: it crosses and the F50 floor
   takes over. So Q_handover := the Q where <n_XY>_interior(Q) = 1. Monotone => a single crossing.

2. THE CHAIN SOLUTION = the coherence horizon Q*(N). The open XY chain is FILLING-DEGENERATE
   (free-fermion/OBC: every (p,p) block ties; (1,1)-only handover == all-p handover), so the
   crossing is the single-excitation {0,2}-coherence point = Q*(N) (F2b corollary): a square-root
   EP / coalescence. It coincides with Q*(N) EXACTLY at the clean-2x2 N=2,3 (a tangency: darkness
   touches 1 at the EP, no transversal crossing) and sits just BELOW it by the trace dressing
   O((tr-1)^2) at N>=4. So the chain handover inherits an already-characterized ladder
   (closed 2/sqrt(c) at N=2,3, transcendental ~0.59N at N>=4).

3. THE RING SOLUTION = a distinct (2,2) free-fermion level crossing, GROWING ~linearly (NOT
   saturating). The wrap bond breaks filling-degeneracy (Fermi degeneracy at half-filling); the
   darkest interior is the (2,2) two-fermion seam (in pure XY a free-fermion dephasing mode, NOT a
   Hamiltonian bound pair). Its handover is a LEVEL CROSSING (|Im| ~ 1e-15, the frozen (2,2) mode
   meets the floor), a different SECTOR/mechanism than the single-excitation SE-EP (a coalescence).
   It grows ~linearly, Q_h ~ N/sqrt(c_eff) with c_eff ~ 12 FLAT in N (so faster than sqrt(N); not
   saturating) = ~4x the chain's darkness constant, so ~half the chain's handover slope (~0.29N vs
   ~0.59N). The handover and the ring SE-EP are mechanistically distinct but their VALUES CROSS near
   N~=10 (handover slope ~0.24 < SE-EP slope ~0.32): the benzene "split" (2.0 vs 1.609, gap 0.39) is
   a small-N feature, NOT a clean universal separation. N=6 ~ exactly 2 is a hexagon coincidence.
"""
import numpy as np
from itertools import combinations


# ---------- builders (XY, sector-projected; convention of incompleteness_survivor.py) ----------
def bonds(N, topology):
    if topology == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topology == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    raise ValueError(topology)


def H_p(N, p, J, bnds):
    states = [sum(1 << i for i in c) for c in combinations(range(N), p)]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), complex)
    for a, b in bnds:
        for s in states:
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += J
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += J
    return H, states


def block_slowest(N, prow, pcol, J, g, bnds):
    """(<n_XY>, |Im|) of the slowest non-kernel mode of the (prow,pcol) coherence block; None if all-kernel."""
    Hr, sr = H_p(N, prow, J, bnds)
    Hc, sc = H_p(N, pcol, J, bnds)
    L = -1j * (np.kron(Hr, np.eye(len(sc))) - np.kron(np.eye(len(sr)), Hc.T))
    deph = np.array([-2.0 * g * bin(sr[a] ^ sc[b]).count("1")
                     for a in range(len(sr)) for b in range(len(sc))])
    ev = np.linalg.eigvals(L + np.diag(deph))
    nz = ev[ev.real < -1e-9]
    if len(nz) == 0:
        return None
    gap = nz.real.max()
    im = float(np.abs(nz[np.abs(nz.real - gap) < 1e-7].imag).max())
    return -gap / (2.0 * g), im


def interior_nxy(N, J, g, bnds, ps):
    """darkness of the global interior survivor = min over given (p,p) of the block slowest <n_XY>."""
    best = None
    for p in ps:
        r = block_slowest(N, p, p, J, g, bnds)
        if r is None:
            continue
        if best is None or r[0] < best[0]:
            best = (r[0], p, r[1])
    return best


def handover_Q(N, topology, ps, J=1.0, lo=0.3, hi=12.0, tol=1e-6):
    """Q where <n_XY>_interior(Q) = 1 (interior survivor as bright as the F50 band-edge floor)."""
    bnds = bonds(N, topology)
    f = lambda Q: interior_nxy(N, J, J / Q, bnds, ps)[0] - 1.0
    if f(lo) > 0 or f(hi) < 0:
        return None
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if f(m) < 0 else (lo, m)
    return 0.5 * (lo + hi)


def qstar_se(N, J=1.0, ring=False, lo=0.02, hi=12.0):
    """The single-excitation {0,2}-coherence EP Q*(N) (F2b corollary), via benzene_two_clocks' bisection."""
    def L_se(g):
        h = np.zeros((N, N), complex)
        for i in range(N - 1):
            h[i, i + 1] = h[i + 1, i] = J
        if ring:
            h[0, N - 1] = h[N - 1, 0] = J
        I = np.eye(N)
        L = -1j * (np.kron(h, I) - np.kron(I, h.T))
        deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
        return L + np.diag(deph)
    for _ in range(80):
        m = 0.5 * (lo + hi)
        nz = np.linalg.eigvals(L_se(m))
        nz = nz[nz.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        band = nz[np.abs(nz.real - nz.real.max()) < 1e-7]
        lo, hi = (m, hi) if np.abs(band.imag).max() > 1e-7 else (lo, m)
    return J / (0.5 * (lo + hi))


# ============================ 1. THE CONDITION: the F50 floor ============================
def _assert_condition():
    # 1a. the (0,1) band edge sits at darkness EXACTLY 1 (Re=-2g), the F50 weight-1 floor, every N/topology.
    for N in (4, 5, 6):
        for topo in ("chain", "ring"):
            nxy, _ = block_slowest(N, 0, 1, 1.0, 1.0 / 1.7, bonds(N, topo))
            assert abs(nxy - 1.0) < 1e-9, f"{topo} N={N}: (0,1) band edge <n_XY>={nxy} != 1 (F50 floor)"
    # 1b. monotone single crossing: <n_XY>_interior(Q) rises through 1 monotonically (ring N=6, (2,2)).
    bnds = bonds(6, "ring")
    grid = np.linspace(1.0, 2.5, 16)
    vals = [interior_nxy(6, 1.0, 1.0 / Q, bnds, [2])[0] for Q in grid]
    assert all(vals[i + 1] > vals[i] for i in range(len(vals) - 1)), "interior <n_XY>(Q) not monotone"
    assert vals[0] < 1.0 < vals[-1], "interior <n_XY> does not cross 1 on [1.0, 2.5]"
    print("[1] THE CONDITION: the (0,1) band edge sits at the F50 floor <n_XY>=1 exactly (every N, topology);")
    print("    the interior survivor's <n_XY>(Q) rises MONOTONICALLY through it -> a single handover Q.")


# ============================ 2. THE CHAIN SOLUTION = Q*(N) ============================
def _assert_chain():
    print("[2] CHAIN: filling-degenerate => handover = the coherence horizon Q*(N).")
    print(f"      {'N':>2} | {'(1,1)-only':>10} {'all-p':>9} {'Q*(N)':>9} | gap below Q* (trace dressing)")
    prev_gap = -1.0
    for N in (4, 5, 6):
        q11 = handover_Q(N, "chain", [1])
        qall = handover_Q(N, "chain", list(range(1, N)))
        qs = qstar_se(N, ring=False)
        # 2a. free-fermion filling-degeneracy: the (1,1) block alone == every block.
        assert abs(q11 - qall) < 2e-3, f"chain N={N}: (1,1)-only {q11} != all-p {qall} (not filling-degenerate)"
        # 2b. handover sits just BELOW Q*(N) by the trace dressing, the gap growing with N.
        gap = qs - qall
        assert 0 < gap < 0.02, f"chain N={N}: handover {qall} not just below Q*={qs} (gap {gap})"
        assert gap > prev_gap, f"chain N={N}: trace-dressing gap {gap} not growing (prev {prev_gap})"
        prev_gap = gap
        print(f"      {N:>2} | {q11:>10.5f} {qall:>9.5f} {qs:>9.5f} | {gap:>+.5f}")
    print("    (1,1)-only == all-p (free-fermion degeneracy); handover just below Q*(N), gap growing")
    print("    (= the trace dressing; they coincide exactly only at the clean-2x2 N=2,3).")


# ============================ 3. THE RING SOLUTION = (2,2) level crossing, growing ============================
def _assert_ring():
    print("[3] RING: a distinct (2,2) free-fermion level crossing, GROWING ~linearly (not saturating).")
    # 3a. the darkest interior is (2,2), darker than the single-excitation (1,1).
    bnds6 = bonds(6, "ring")
    nxy22 = block_slowest(6, 2, 2, 1.0, 1.0 / 1.5, bnds6)[0]
    nxy11 = block_slowest(6, 1, 1, 1.0, 1.0 / 1.5, bnds6)[0]
    assert nxy22 < nxy11, f"ring N=6 Q=1.5: (2,2) {nxy22} not darker than (1,1) {nxy11}"
    print(f"      darkest interior is (2,2): <n_XY>(2,2)={nxy22:.4f} < (1,1)={nxy11:.4f} (the V-Effect seam)")

    # 3b/3c/3d: handover value, level crossing (|Im|~0), the SE-EP curve it crosses, growth.
    print(f"      {'N':>2} | {'ring Q_h':>9} {'|Im|@Qh':>9} {'Q*ring':>8} {'qh-q*':>7} {'c_eff':>7}")
    Qh, split = {}, {}
    for N in (6, 8, 10):
        qh = handover_Q(N, "ring", [2])
        nxy, _, im = interior_nxy(N, 1.0, 1.0 / qh, bonds(N, "ring"), [2])
        qs = qstar_se(N, ring=True)
        # robust mechanism distinction: the (2,2) seam handover is a FROZEN level crossing (|Im|~0),
        # a different sector than the single-excitation SE-EP. (The numerical gap to Q*ring is NOT a
        # clean separation: the two curves have different slopes and CROSS near N~10 - see below.)
        assert im < 1e-10, f"ring N={N}: handover |Im|={im} not a frozen level crossing"
        Qh[N], split[N] = qh, qh - qs
        print(f"      {N:>2} | {qh:>9.5f} {im:>9.1e} {qs:>8.4f} {qh-qs:>+7.3f} {(N/qh)**2:>7.3f}")
    # 3c. the benzene (N=6) split is robust (the V-Effect seam), but it SHRINKS and the curves cross:
    assert split[6] > 0.3, f"benzene split {split[6]} not the robust V-Effect seam (>0.3)"
    assert split[6] > split[8] > split[10], f"the handover-vs-SE-EP gap does not shrink: {split}"
    assert split[10] < 0.1, f"the curves have not nearly crossed by N=10 (split {split[10]})"
    print(f"    the (2,2) seam handover is a frozen LEVEL CROSSING (|Im|~0), distinct sector from the")
    print(f"    single-excitation SE-EP; their VALUES cross near N~10 (gap {split[6]:+.2f}/{split[8]:+.2f}/{split[10]:+.2f}")
    print(f"    at N=6/8/10) - the benzene split is a small-N feature, not a clean universal separation.")
    # 3e. GROWS (not saturating), and faster than sqrt(N): c_eff = (N/Q_h)^2 stays ~flat (linear),
    #     whereas sqrt(N) growth would force c_eff ∝ N (c_eff(10)/c_eff(8) -> 10/8 = 1.25).
    assert Qh[10] > Qh[8] > Qh[6], f"ring handover not growing: {Qh}"
    ceff_ratio = (10 / Qh[10]) ** 2 / ((8 / Qh[8]) ** 2)
    assert ceff_ratio < 1.15, f"ring c_eff rising like sqrt(N) (ratio {ceff_ratio} ~ 1.25), not linear"
    print(f"    GROWS: Q_h(6,8,10) = {Qh[6]:.3f}, {Qh[8]:.3f}, {Qh[10]:.3f}; c_eff flat "
          f"(ratio 10:8 = {ceff_ratio:.3f} << 1.25=sqrt-N) -> LINEAR ~0.29N, NOT saturating, faster than sqrt(N).")


if __name__ == "__main__":
    _assert_condition()
    _assert_chain()
    _assert_ring()
    print("\nAll asserts passed. The handover Q has a CLOSED CONDITION (the diagonal incompleteness")
    print("survivor reaches the F50 off-diagonal floor <n_XY>=1); the CHAIN solution is the coherence")
    print("horizon Q*(N) (free-fermion filling-degeneracy); the RING solution is a distinct (2,2)")
    print("free-fermion level crossing that grows ~linearly (~0.29N), neither co-located nor saturating.")
