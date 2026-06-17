#!/usr/bin/env python3
"""birth_canal_horizon_junction (Arc A): the NATURE of the interior {0,2} survivor, and whether the
gamma-PROFILE or the TOPOLOGY sets it.

The junction: as Q=J/g falls, the global-slowest Liouvillian mode hands the crown from the odd (0,1)
|1-exc><vac| band edge (F50 floor, rate=2g exactly) to an even diagonal {0,2}-coherence interior
survivor. handover_q.py pinned the interior survivor's nature at UNIFORM gamma:
  - CHAIN: filling-degenerate (free-fermion OBC) => the single-excitation (1,1) {0,2} SE-EP, a
    square-root COALESCENCE.
  - RING: the wrap bond breaks filling-degeneracy => the (2,2) two-excitation seam, a frozen LEVEL
    CROSSING.
birth_canal_n6.py found that on a DEEP-EDGE (non-uniform) chain the low-Q global slowest is the (2,2)
{0,2}-coherence -- but never measured its NATURE. The open question:

  *** Does a non-uniform gamma-profile turn the CHAIN's (2,2) {0,2} mode from part of the oscillating
      SE-EP coalescence into a frozen (2,2) level crossing -- i.e. can the chain host the ring's
      V-Effect seam once the profile breaks free-fermion filling-degeneracy? ***

THE DISCRIMINATOR (the lens-correct one -- NOT |Im|). Below Q* the SE-EP pair is OVERDAMPED REAL
(|Im|=0), so |Im|~0 alone CANNOT tell a coalescence from a frozen level crossing. The sharp test is
the PETERMANN PHASE RIGIDITY r = |vL^H vR| / (||vL|| ||vR||): at an EP the two eigenvectors coalesce
=> r -> 0; at a level crossing the modes stay distinct => r stays O(1). (Same instrument as
CoherenceHorizonWitness.isEp = rigidity < 0.05.) We scan Q and read the MIN rigidity of a sector's
slowest mode: a dip to ~0 = EP/coalescence; staying O(1) = level crossing.

GATE-FIRST hypothesis (Stage 0, falsifiable; from handover_q's chain=coalescence / ring=crossing):
  the chain's interior survivor is an EP COALESCENCE (rigidity dips ~0) REGARDLESS of profile; only
  the ring is the frozen level crossing (rigidity stays O(1)). If the deep-edge chain (2,2) does NOT
  dip (rigidity stays O(1)), the gate FIRES -> the profile changed the nature. Diagnose, do not loosen.

Convention (matches handover_q.py / birth_canal_n6.py): H = sum_bond (XX+YY) hopping (off-diag J=1),
dephasing rate -2*sum_l g_l*[bit differs]; profile shapes sum to N (mean gamma = g), Q = J/g.
"""
import numpy as np
from itertools import combinations
from scipy.linalg import eig


# ---------------- sector-projected builders (handover_q.py generalized to a per-site gamma profile) ----
def bonds(N, topology):
    if topology == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topology == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    if topology == "star":
        return [(0, i) for i in range(1, N)]
    raise ValueError(topology)


def H_p(N, p, bnds, J=1.0):
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


def sector_L(N, prow, pcol, profile, bnds, J=1.0):
    """The (prow,pcol) joint-popcount coherence block of L with a per-site gamma PROFILE."""
    Hr, sr = H_p(N, prow, bnds, J)
    Hc, sc = H_p(N, pcol, bnds, J)
    L = -1j * (np.kron(Hr, np.eye(len(sc))) - np.kron(np.eye(len(sr)), Hc.T))
    deph = np.empty(len(sr) * len(sc))
    k = 0
    for a in range(len(sr)):
        for b in range(len(sc)):
            diffbits = sr[a] ^ sc[b]
            deph[k] = -2.0 * sum(profile[l] for l in range(N) if (diffbits >> l) & 1)
            k += 1
    return L + np.diag(deph)


def slowest(L):
    """(rate, |Im|, rigidity) of the slowest non-kernel mode. rigidity = Petermann r of THAT mode."""
    w, vl, vr = eig(L, left=True, right=True)
    nz = np.where(w.real < -1e-9)[0]
    if len(nz) == 0:
        return None
    i = nz[np.argmax(w[nz].real)]
    r = abs(np.vdot(vl[:, i], vr[:, i])) / (np.linalg.norm(vl[:, i]) * np.linalg.norm(vr[:, i]))
    return -w[i].real, abs(w[i].imag), float(r)


def uniform(N):
    return [1.0] * N


def deep_edge(N, edge=0.25):
    rest = (N - 2 * edge) / (N - 2)
    return [edge] + [rest] * (N - 2) + [edge]


def scan_sector(N, topo, prow, pcol, shape, qlo, qhi, npts=25, J=1.0):
    """Scan Q in [qlo,qhi]; return (min_rigidity, Q@min, max|Im|) of the (prow,pcol) slowest mode.
    min_rigidity -> ~0 signals an EP coalescence somewhere in the window; staying O(1) = level crossing."""
    bnds = bonds(N, topo)
    min_rig, q_at, max_im = 1.0, None, 0.0
    for Q in np.linspace(qlo, qhi, npts):
        g = J / Q
        res = slowest(sector_L(N, prow, pcol, [c * g for c in shape], bnds, J))
        if res is None:
            continue
        _, im, rig = res
        if rig < min_rig:
            min_rig, q_at = rig, Q
        max_im = max(max_im, im)
    return min_rig, q_at, max_im


# ============================ CONTROL: the (0,1) band edge = the F50 floor ============================
def _control_floor():
    for N in (5, 6):
        for topo in ("chain", "ring"):
            g = 1.0 / 1.7
            rate, _, _ = slowest(sector_L(N, 0, 1, [g] * N, bonds(N, topo)))
            assert abs(rate - 2 * g) < 1e-9, f"{topo} N={N}: (0,1) rate {rate} != 2g (F50 floor)"
    print("[control] (0,1) band edge sits at the F50 floor rate=2g exactly (uniform g, every N/topology).")


# ===================== STAGE 0 (THE GATE): frozen level crossing vs oscillating EP ====================
# Calibrated discriminator (from the first run): the CLEAN separator is max|Im| over a Q-scan that
# brackets the EP -- a frozen level crossing never oscillates (max|Im| ~ 0); an EP oscillates above Q*
# (max|Im| = O(1)), even though it is overdamped-REAL below Q* (so a single-Q |Im| is blind). The
# rigidity DIP (Petermann r -> 0) is the secondary EP confirmation (chain (1,1) dips to ~0.09 at Q*).
FROZEN_IM = 1e-6
EP_RIG = 0.12


def _gate_nature():
    print("\n[GATE / Stage 0] interior {0,2} survivor: FROZEN level crossing (max|Im|~0 over the scan)")
    print("  vs oscillating EP coalescence (max|Im|=O(1), rigidity dips ~0 at Q*).\n")
    N = 6
    cases = [
        ("chain (1,1) SE",   "chain", 1, 1, "unif", uniform(N),   "EP anchor (coherence horizon Q*)"),
        ("ring  (2,2)",      "ring",  2, 2, "unif", uniform(N),   "FROZEN anchor (V-Effect seam)"),
        ("chain (2,2)",      "chain", 2, 2, "unif", uniform(N),   "filling-degenerate with (1,1)?"),
        ("chain (2,2)",      "chain", 2, 2, "deep", deep_edge(N), "THE OPEN ONE"),
        ("ring  (2,2)",      "ring",  2, 2, "deep", deep_edge(N), "ring under the profile"),
    ]
    res = {}
    print(f"  {'sector/topo':<18}{'prof':<6}{'max|Im|':>10}{'min rig':>9}{'Q@minrig':>10}   nature")
    for label, topo, pr, pc, prof, shape, note in cases:
        mr, qa, mi = scan_sector(N, topo, pr, pc, shape, qlo=1.0, qhi=4.0, npts=25)
        res[(label, prof)] = (mr, qa, mi)
        frozen = mi < FROZEN_IM
        nat = "FROZEN crossing" if frozen else ("oscillating EP" if mr < EP_RIG else "oscillating")
        print(f"  {label:<18}{prof:<6}{mi:>10.2e}{mr:>9.3f}{(qa or 0):>10.2f}   {nat} | {note}")

    # Anchors: chain (1,1) is the SE-EP (oscillates + rigidity dips at Q*); ring (2,2) uniform is the
    # FROZEN V-Effect level crossing (never oscillates). The clean separator is max|Im|.
    se_im, se_rig = res[("chain (1,1) SE", "unif")][2], res[("chain (1,1) SE", "unif")][0]
    ring_im = res[("ring  (2,2)", "unif")][2]
    assert se_im > 1e-2 and se_rig < EP_RIG, f"anchor: chain (1,1) not an oscillating EP (|Im|={se_im}, rig={se_rig})"
    assert ring_im < FROZEN_IM, f"anchor: ring (2,2) uniform not frozen (|Im|={ring_im})"
    # filling-degeneracy: chain (2,2) uniform == chain (1,1) (free-fermion OBC).
    assert abs(res[("chain (2,2)", "unif")][2] - se_im) < 1e-3, "chain (2,2) unif should equal (1,1) (filling-degenerate)"
    print(f"\n  anchors: chain(1,1) SE-EP oscillates |Im|={se_im:.3f} & rigidity dips {se_rig:.3f} at Q*; "
          f"ring(2,2)-unif FROZEN |Im|={ring_im:.1e}. chain(2,2)-unif equals (1,1) (filling-degenerate).")

    # THE FINDING-GATE: the deep-edge chain (2,2) does NOT inherit the ring's frozen character.
    deep_im = res[("chain (2,2)", "deep")][2]
    assert deep_im > 1e-2, (
        f"GATE FIRED: deep-edge chain (2,2) is FROZEN (|Im|={deep_im:.1e}) -- it DID inherit the ring's "
        f"V-Effect level crossing under the profile. Re-think; bank it.")
    print(f"  FINDING: deep-edge chain (2,2) OSCILLATES (|Im|={deep_im:.3f}) -> its slowest (2,2) survivor")
    print("  does NOT inherit the ring's frozen level crossing. The frozen seam is RING-specific (a")
    print("  topology feature of the wrap bond); Stage 1 shows it is robust to the gamma-profile.")
    return res


# ===================== STAGE 1: the freezing is a ring+uniform symmetry (continuation) ================
def _stage1_unfreeze_continuation():
    """Follow the ONE ring (2,2) seam mode (frozen at uniform gamma) as the profile is turned on
    (uniform -> deep-edge). If its |Im| rises from ~0, the freezing is symmetry-protected (broken by
    the non-uniform profile) -- the same mode, tracked, not a relabel."""
    print("\n[Stage 1] is the ring (2,2) freezing symmetry-protected? Track the SAME seam mode as the")
    print("  profile turns on (uniform -> deep-edge), Q fixed at the uniform handover.")
    N = 6
    Q = 1.62                                      # the ring (2,2) handover Q (from the run, Q@frozen)
    g = 1.0 / Q
    u, d = np.array(uniform(N)), np.array(deep_edge(N))
    bnds = bonds(N, "ring")
    prev_lam, ims = None, []
    for t in np.linspace(0.0, 1.0, 11):
        prof = list(((1 - t) * u + t * d) * g)
        L = sector_L(N, 2, 2, prof, bnds)
        w = np.linalg.eigvals(L)
        nz = w[w.real < -1e-9]
        if prev_lam is None:                      # seed: the slowest (the frozen seam) at uniform
            lam = nz[np.argmax(nz.real)]
        else:                                     # continuation: nearest eigenvalue to the tracked one
            lam = nz[np.argmin(np.abs(nz - prev_lam))]
        prev_lam = lam
        ims.append(abs(lam.imag))
    print(f"  |Im| of the tracked seam mode, t=0(uniform)->1(deep-edge): "
          + ", ".join(f"{x:.2e}" for x in ims[::2]))
    # CORRECTED FINDING (the Stage-0 'ring deep-edge oscillates |Im|=3.4' was a different, oscillating
    # (2,2) mode RELABELED as the slowest -- the mode-tracking trap, hazard 7). Tracking the SAME seam
    # mode shows it stays frozen: the freezing is ROBUST to the gamma-profile, not a uniform-gamma artifact.
    assert ims[0] < FROZEN_IM, f"seed mode not frozen at uniform (|Im|={ims[0]})"
    assert max(ims) < FROZEN_IM, (
        f"GATE: the tracked ring (2,2) seam acquired |Im|={max(ims):.1e} under the profile -- it un-froze. "
        f"Re-think (this would mean the freezing IS a uniform-gamma symmetry).")
    print(f"  the SAME seam mode stays FROZEN: |Im| <= {max(ims):.1e} for all t. The ring (2,2) V-Effect")
    print("  freezing is ROBUST to the gamma-profile (intrinsic to the wrap-bond two-fermion structure);")
    print("  the Stage-0 'ring deep-edge oscillates' was a DIFFERENT (2,2) mode relabeled slowest (hazard 7:")
    print("  the boundary is a MIN OVER A MODE CROSSING, vindicating the arc's own warning).")


# ===================== STAGE 2: THE SEAM -- sterile/birth-canal IS the odd<->even junction ============
# PostEpFlowField (C#) reads the sterile<->birth-canal boundary as the Q-dependence of the GLOBAL
# slowest rate: BirthCanalDeviation = rate(Q=1000) - rate(Q=1.5); IsInSterileZone = |Deviation|~0.
# Arc-A reads the SAME boundary as the slowest-SECTOR identity switch (odd (0,1) <-> even (2,2)).
# The connection (the seam both arcs overlooked): the (0,1) band edge is the -2g Q-FLAT floor
# (Absorption Theorem), the interior (2,2) is -2g<n_XY>(Q) Q-DEPENDENT. So:
#   sterile  <=> the (0,1) floor is the slowest at BOTH Q (rate Q-flat, Deviation~0)
#   birth canal <=> the interior (2,2) has undercut the floor at low Q (rate Q-modulated, Deviation>0)
# GATE: IsInBirthCanal MUST coincide bit-for-bit with "the slowest at low Q is NOT the (0,1) band edge".
def _full_L(N, Q, profile):
    """The full 4^N Liouvillian (birth_canal_n6 convention: Q*H_unit + per-site profile dephasing)."""
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    Zp = np.array([[1, 0], [0, -1]], complex)

    def op_at(s, P):
        o = np.array([[1]], complex)
        for i in range(N):
            o = np.kron(o, P if i == s else I2)
        return o

    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * Q * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = op_at(l, Zp)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def _global_slowest_sector(N, Q, profile):
    """(rate, dominant (popcount_ket, popcount_bra) sector) of the global slowest non-kernel mode."""
    d = 2 ** N
    w, V = np.linalg.eig(_full_L(N, Q, profile))
    idx = np.where(np.abs(w) > 1e-7)[0]
    slow = idx[np.argmax(w[idx].real)]
    rho = V[:, slow].reshape(d, d)
    wts = np.abs(rho) ** 2
    pc = np.array([bin(x).count("1") for x in range(d)])
    blk = {}
    for a in range(d):
        for b in range(d):
            if wts[a, b] > 1e-9:
                blk[(pc[a], pc[b])] = blk.get((pc[a], pc[b]), 0.0) + wts[a, b]
    dom = max(blk.items(), key=lambda kv: kv[1])[0]
    return -float(w[slow].real), dom


def _stage2_seam():
    print("\n[Stage 2 / THE SEAM] how PostEpFlowField's sterile<->birth-canal boundary relates to Arc-A's")
    print("  odd(0,1)<->even(2,2) junction. Both read rate_slow(Q) = min over sectors. CORRECTED (a naive")
    print("  'same boundary' gate fired): the birth canal has TWO mechanisms, and the junction is one of them.\n")
    QLO, QHI, DEV_TOL = 1.5, 1000.0, 1e-4         # PostEpFlowField's probe points + tolerance (full_L convention)
    # the SURVIVOR'S CHANGE-NUMBER dn=|p_col-p_row| is the robust junction discriminator (NOT the sector
    # tuple, which is fooled by F1-conjugate degeneracy: (4,3) and (3,4) are the SAME dn=1 band edge).
    # dn=1 = number-CHANGING band edge (the (0,1)-class, -2g floor); dn=0 = number-CONSERVING interior.
    cases = [
        (6, deep_edge(6),                "deep-edge N=6"),   # junction (the n6 finding: interior wins at lo Q)
        (5, [0.25, 1.5, 1.5, 1.5, 0.25], "canal N=5"),       # odd-drift (edge-protected (0,1) survivor, drifts)
        (5, [1.0, 1.0, 1.0, 1.0, 1.0],   "uniform N=5"),     # sterile (window above Q*(5), band edge Q-flat)
    ]
    print(f"  {'profile':<16}{'Deviation':>10}{'birth?':>7}{'dn@loQ':>7}{'dn@hiQ':>7}  mechanism")
    n_junction = n_odddrift = n_sterile = 0
    for Nc, prof, label in cases:
        r_lo, sec_lo = _global_slowest_sector(Nc, QLO, prof)
        r_hi, sec_hi = _global_slowest_sector(Nc, QHI, prof)
        dn_lo, dn_hi = abs(sec_lo[0] - sec_lo[1]), abs(sec_hi[0] - sec_hi[1])
        deviation = r_hi - r_lo
        birth = abs(deviation) > DEV_TOL
        switched = dn_lo != dn_hi                    # the junction: number-change of the survivor flips
        if switched:
            mech, n_junction = "JUNCTION (interior dn=0 wins at lo Q -> band edge dn=1 at hi Q)", n_junction + 1
        elif birth:
            mech, n_odddrift = "odd-drift (same dn=1 survivor, non-uniform gamma)", n_odddrift + 1
        else:
            mech, n_sterile = "sterile (dn=1 band edge, Q-flat -2g)", n_sterile + 1
        print(f"  {label:<16}{deviation:>10.4f}{str(birth):>7}{dn_lo:>7}{dn_hi:>7}  {mech}")
        # CORRECTED SEAM GATE: a junction (the survivor's dn switching) IMPLIES the birth canal (different-
        # dn slowest modes at the two Q's => rate must change). The converse is FALSE (birth canal can be
        # the dn=1 survivor drifting under non-uniform gamma, no dn switch -- the canal profile).
        assert (not switched) or birth, (
            f"GATE FIRED: '{label}' dn switched ({dn_lo}->{dn_hi}) but Deviation={deviation:.1e} reads sterile "
            f"-- a junction with no rate change is impossible. Diagnose.")
    # all three mechanisms realized: the junction is a STRICT sub-mechanism of the birth canal.
    assert n_junction >= 1, "no junction seen -- expected N=6 deep-edge (interior wins at the low probe Q)"
    assert n_odddrift >= 1, "no birth-canal-without-junction -- expected the edge-protected N=5 canal (odd-drift)"
    assert n_sterile >= 1, "no sterile case -- expected uniform N=5 (probe window above Q*(5))"
    print(f"\n  SEAM (corrected): junction => birth canal (held); converse FALSE. {n_junction} profile(s) enter")
    print("  the birth canal by a JUNCTION (the even (2,2) survivor overtaking = Arc-A / the handover Q*(N)),")
    print(f"  {n_odddrift} by ODD-DRIFT (the same (0,1) edge-survivor drifting under non-uniform gamma; the")
    print("  (0,1) block -iQh-2diag(g) is Q-flat ONLY at uniform g). So PostEpFlowField's birth canal is the")
    print("  UNION; Arc-A's junction is the survivor-identity-change sub-mechanism. The two arcs are one")
    print("  object (rate_slow(Q)), meeting where the (2,2) interior undercuts the (0,1) floor.")


if __name__ == "__main__":
    _control_floor()
    _gate_nature()
    _stage1_unfreeze_continuation()
    _stage2_seam()
    print("\nAll stages passed.")
