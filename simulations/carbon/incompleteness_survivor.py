"""(b) The dynamic survival probe: is the incomplete (C=0.5, interior-filling) coherence ALWAYS
the longest-lived survivor, over N and topology - and how does its lifetime inherit in N?

Thread (a) established the survival LAW and its algebra: longest-lived = darkest = min <n_XY>
(the Absorption Theorem, a_0 = 2 on the Pi2 dyadic ladder), and the V-Effect/incompleteness
= C=0.5 = a_2 = the Absorption Theorem's INVERSION-MIRROR PARTNER (a_0*a_2 = 1, both forced by
d^2-2d=0; Pi2DyadicLadderClaim). This script tests the dynamic ENACTMENT: WHERE the survivor
lives (does the incomplete win?) and HOW its lifetime scales with N (the separate inheritance).

Carbon convention (XY, free fermions, no ZZ; matches benzene_two_clocks): H = (J/2) sum(XX+YY),
dephasing g*sum_l D[Z_l], Re(lambda) = -2g*<n_XY>, Q = J/g. Survivor lifetime proxy = <n_XY>
(smaller = darker = longer-lived). NOTE this is XY: CHAIN_GAP_SECTOR_DIAGNOSTIC's half-filling
results are HEISENBERG (with ZZ); see _xy_vs_heisenberg_slowmode.py for why the two differ.

Sector-projected: the XY H conserves excitation number and Z-dephasing is diagonal, so L is
block-diagonal in (popcount_ket, popcount_bra). The global slowest lives in a low-light block:
a diagonal (p,p) dressed magnon-admixture (fractional <n_XY>), or the (0,1) odd band edge
(<n_XY>=1). Validated bit-for-bit vs the full 4^N L at N=4 (all three topologies).

FINDINGS (run it; asserts pin the anchors):
  WHERE: the even ring puts the strong-dephasing survivor in the INTERIOR (2 <= p <= N-2, the
    incompleteness region, NOT the w=0/N extremes) and the XY chain and the odd ring tie every filling; above
    a handover Q it switches to the (0,1) odd band edge. In XY the even ring's survivor sits in its even fillings,
    (2,2)/(N-2,N-2) at N = 6, off the dead centre there (the wrapped string's split; the half filling joins them when
    4 divides N), the odd ring ties every filling (read; a pi flux through it changes no block's spectrum,
    PROOF_WEIGHT1_DEGENERACY), and the open chain ties every filling (measured
    below; at the Zeno end proved by the ladder of PROOF_CODIM1_BY_ADDITIVITY with Theorem E (c)), so the
    printed chain sector is the solver's pick from a tie; the
    Heisenberg ZZ picks the half filling (N/2,N/2), the central pair at odd N (CHAIN_GAP; PROOF_WEIGHT1_DEGENERACY Theorem E (d) at
    the Zeno end). The STAR puts its survivor at the popcount BOUNDARY (1,1)/(N-1,N-1): with ZZ a hop
    between the hub and an arm is detuned by the other arms' imbalance (the Zeno end's detuning form,
    the hops' own share read on the stars N = 4..7), and in XY the
    boundary wins only at the order (J/gamma)^6 at the Zeno end, consistent with its exchange statistics
    (free fermions on the same star tie). So "the incomplete survives
    longest" is a statement about the even ring and the ZZ chain; the star is the counterexample.
  HOW LONG: the interior survivor's light content <n_XY> ~ c*Q^2/N^2 (the magnon-admixture inheritance),
    ring ~ 4x chain - a SEPARATE 1/N^2 inheritance from the Pi2 dyadic ladder (which carries the
    constants, not the N-scaling)."""
import sys
import numpy as np
from itertools import combinations

sys.path.insert(0, "simulations/carbon")


def bonds(N, topology):
    if topology == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topology == "ring":
        return [(i, i + 1) for i in range(N - 1)] + [(N - 1, 0)]
    if topology == "star":
        return [(0, k) for k in range(1, N)]
    raise ValueError(f"unknown topology {topology!r}")


def basis(N, p):
    return [sum(1 << i for i in c) for c in combinations(range(N), p)]


def H_p(N, p, J, bnds):
    """p-excitation hopping Hamiltonian (XX+YY, amplitude J per hop) on the bond list."""
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


def slowest_in_block(N, prow, pcol, J, g, bnds):
    """Slowest non-kernel (Re, |Im|) of the (prow,pcol) coherence block, or None if all-kernel."""
    Hr, sr = H_p(N, prow, J, bnds)
    Hc, sc = H_p(N, pcol, J, bnds)
    L = -1j * (np.kron(Hr, np.eye(len(sc))) - np.kron(np.eye(len(sr)), Hc.T))
    deph = np.array([-2.0 * g * bin(sr[a] ^ sc[b]).count("1")
                     for a in range(len(sr)) for b in range(len(sc))])
    ev = np.linalg.eigvals(L + np.diag(deph))
    nz = ev[ev.real < -1e-7]
    if len(nz) == 0:
        return None
    gap = nz.real.max()
    im = float(np.abs(nz[np.abs(nz.real - gap) < 1e-6].imag).max())
    return gap, im


def survivor(N, J, g, bnds):
    """Global slowest over the diagonal (p,p) blocks (p=1..N-1) + the (0,1) odd band edge.
    The slowest of the full L lives among these low-light candidates. Returns
    (Re, |Im|, sector, n_xy=-Re/2g)."""
    cands = [(p, p) for p in range(1, N)] + [(0, 1)]
    best = None
    for prow, pcol in cands:
        r = slowest_in_block(N, prow, pcol, J, g, bnds)
        if r is None:
            continue
        re, im = r
        if best is None or re > best[0]:
            best = (re, im, (prow, pcol))
    re, im, sec = best
    return re, im, sec, -re / (2 * g)


# ---- full 4^N L (generalized to a bond list) for the N=4 validation ----
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
    H = np.zeros((d, d), complex)
    for a, b in bnds:
        H += (J / 2) * (site(X, a) @ site(X, b) + site(Y, a) @ site(Y, b))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l)
        L += g * (np.kron(Zl, Zl) - np.kron(Id, Id))
    ev = np.linalg.eigvals(L)
    re = ev[ev.real < -1e-7].real.max()
    return re, -re / (2 * g)


# ---- 1. VALIDATION: the sector survivor reproduces the full 4^N L slowest (N=4, all topologies) ----
def _assert_sector_validates():
    for topo in ("chain", "ring", "star"):
        bnds = bonds(4, topo)
        for Q in (1.5, 10.0):
            g = 1.0 / Q
            re_s = survivor(4, 1.0, g, bnds)[0]
            re_f = full_slowest(4, 1.0, g, bnds)[0]
            assert abs(re_s - re_f) < 1e-4, f"{topo} N=4 Q={Q}: sector Re={re_s} != full Re={re_f}"
    print("[1] sector survivor == full 4^N L slowest at N=4 (chain, ring, star; Q=1.5, 10) - method sound")


# ---- 2. WHERE: the even ring's survivor is INTERIOR (incompleteness), the XY chain and the odd ring tie every filling, the star is BOUNDARY ----
def _assert_where():
    g = 1.0 / 1.5  # strong dephasing, below the handover
    print("[2] WHERE the strong-dephasing (Q=1.5) survivor lives:")
    for topo in ("chain", "ring", "star"):
        for N in ((4, 5, 6) if topo == "ring" else (4, 6)):
            re, im, sec, nxy = survivor(N, 1.0, g, bonds(N, topo))
            p = sec[0] if sec[0] == sec[1] else None
            tied = topo == "chain" or (topo == "ring" and N % 2 == 1)
            if tied:
                per_p = [slowest_in_block(N, q, q, 1.0, g, bonds(N, topo))[0] for q in range(1, N)]
                spread = (max(per_p) - min(per_p)) / abs(min(per_p))
            kind = (("band-edge (0,1)" + (f", the diagonal blocks every filling tied, relative spread {spread:.0e}"
                                          if tied else "")) if sec == (0, 1)
                    else f"every filling tied, relative spread {spread:.0e} (the printed sector a pick)" if tied
                    else f"INTERIOR (incompleteness)" if p is not None and 2 <= p <= N - 2
                    else f"BOUNDARY" if p in (1, N - 1)
                    else "?")
            print(f"     {topo:>5} N={N}: survivor sector={sec} <n_XY>={nxy:.4f} "
                  f"|Im|={im:.3f} -> {kind}")
            if topo == "ring" and N % 2 == 0:
                assert sec == (0, 1) or (p is not None and 2 <= p <= N - 2), \
                    f"{topo} N={N}: ring survivor should be interior/band-edge, got {sec}"
            # the open XY chain and the odd XY ring tie every filling: the printed sector is a pick from a tie, no assert
            if topo == "star":
                assert p in (1, N - 1), f"star N={N}: survivor should be boundary (1,1)/(N-1,N-1), got {sec}"
    print("     -> even ring: interior (the incomplete survives); XY chain and odd ring: every filling tied; star: boundary")


# ---- 3. HOW LONG: the interior survivor's lifetime <n_XY> ~ c*Q^2/N^2; ring ~ 4x chain ----
def _report_scaling(Ns=(4, 5, 6, 7)):
    Q = 1.5
    g = 1.0 / Q
    print(f"[3] HOW LONG: interior survivor light content <n_XY> vs N (Q={Q}); c = <n_XY>*N^2/Q^2:")
    print(f"     {'N':>3} | {'chain <n_XY>':>13} {'c':>6} | {'ring <n_XY>':>12} {'c':>6} | {'ring/chain':>10}")
    rows = {}
    for N in Ns:
        rc = survivor(N, 1.0, g, bonds(N, "chain"))
        rr = survivor(N, 1.0, g, bonds(N, "ring"))
        nxy_c, nxy_r = rc[3], rr[3]
        cc = nxy_c * N * N / (Q * Q)
        cr = nxy_r * N * N / (Q * Q)
        ratio = nxy_r / nxy_c if nxy_c > 0 else float("nan")
        rows[N] = (nxy_c, nxy_r, ratio)
        flag = "   (ring: the (0,1) band edge, past its handover)" if rr[2] == (0, 1) else ""
        print(f"     {N:>3} | {nxy_c:>13.5f} {cc:>6.3f} | {nxy_r:>12.5f} {cr:>6.3f} | {ratio:>10.3f}{flag}")
    print("     (chain c ~ 0.55, ring c ~ 2.2 in HEISENBERG/CHAIN_GAP; XY values are this model's own")
    print("      magnon-admixture inheritance - a separate 1/N^2 ladder from the Pi2 dyadic constants.)")
    return rows


if __name__ == "__main__":
    _assert_sector_validates()
    _assert_where()
    _report_scaling()
    print("\nAll asserts passed: the incomplete (interior) survives longest on the even ring, the star is the boundary")
    print("counterexample, and the survivor's lifetime is a 1/N^2 magnon-admixture inheritance distinct from the Pi2")
    print("dyadic constant ladder; the XY chain's and the odd ring's ties across the fillings are printed above, not")
    print("asserted.")
