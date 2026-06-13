"""Does the V-Effect coherence seam track AROMATICITY (Hueckel 4n+2)? It does NOT.

Tested by comparing the anti-aromatic 4n rings cyclobutadiene C4 and cyclooctatetraene C8 (both with
degenerate zero-energy non-bonding MOs) against the aromatic 4n+2 benzene C6. RESULT (2026-06-13):

  The V-Effect seam (the frozen DOUBLE-excitation mode that overtakes the band-edge beat) is
  RING-UNIVERSAL, not aromaticity-specific. At strong dephasing (low Q) ALL rings (C4, C6, C8) hand
  over to a frozen double-excitation mode: the C=0.5 HALF-FILLING boundary, a sibling of the
  incompleteness V-Effect (docs/HIERARCHY_OF_INCOMPLETENESS.md) on the filling axis.

  Aromaticity (4n+2 vs 4n) is REFUTED as the discriminant: the 4n anti-aromatics C4 and C8 do NOT
  group together. C4 is a small-ring ANOMALY (its double-excitation mode is anomalously slow,
  <n_diff> < 1, so it dominates even at WEAK dephasing); C8 (the next 4n) behaves like aromatic C6
  (the oscillating band-edge survivor leads at weak dephasing). The seam is the C=0.5 half-filling
  boundary made dynamic, not Hueckel's rule.

The full 4^8 Liouvillian (C8) is too big for dense eigendecomposition, so we use a SECTOR-PROJECTED
Liouvillian: the XY ring H conserves excitation number and the Z-dephasing is diagonal, so L is
block-diagonal in (popcount(ket), popcount(bra)); the slowest modes live in the |dp|<=1 blocks
(the odd |vac><psi| Uhr-1 in (0,1), the even {0,2}-coherence seam in the diagonal (p,p)). Largest
C8 block (4,4) = 4900. Cross-checked against the full-L slowest (benzene_two_clocks.slowest_full).

Self-validating: run it, every assert must pass. Wall-clock note: the C8 (4,4)/(3,4) blocks make
this ~5-8 minutes; the C4/C6 sector blocks are tiny (instant)."""
import sys
import numpy as np
from itertools import combinations
sys.path.insert(0, 'simulations/carbon')
from benzene_two_clocks import slowest_full, huckel_mos


def basis(N, p):
    return [sum(1 << i for i in c) for c in combinations(range(N), p)]


def H_p(N, p, J, ring):
    """The p-excitation hopping Hamiltonian (XX+YY ring, amplitude J per hop), C(N,p) x C(N,p)."""
    states = basis(N, p)
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), complex)
    bonds = [(b, b + 1) for b in range(N - 1)] + ([(N - 1, 0)] if ring else [])
    for a, b in bonds:
        for s in states:
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += J
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += J
    return H, states


def slowest_in_block(N, prow, pcol, J, g, ring):
    """Slowest non-kernel mode (Re, |Im|) of the (prow,pcol) coherence block, or None if all-kernel."""
    Hr, sr = H_p(N, prow, J, ring)
    Hc, sc = H_p(N, pcol, J, ring)
    L = -1j * (np.kron(Hr, np.eye(len(sc))) - np.kron(np.eye(len(sr)), Hc.T))
    deph = np.array([-2.0 * g * bin(sr[a] ^ sc[b]).count('1') for a in range(len(sr)) for b in range(len(sc))])
    ev = np.linalg.eigvals(L + np.diag(deph))
    nz = ev[ev.real < -1e-7]
    if len(nz) == 0:
        return None
    gap = nz.real.max()
    return gap, np.abs(nz[np.abs(nz.real - gap) < 1e-6].imag).max()


def global_slowest(N, J, g, ring):
    """Global slowest non-kernel mode over the |dp|<=1 blocks (upper triangle; conjugate Re-equal).
    Returns (Re, |Im|, (prow,pcol))."""
    best = None
    for prow in range(N + 1):
        for pcol in (prow, prow + 1):
            if pcol > N:
                continue
            r = slowest_in_block(N, prow, pcol, J, g, ring)
            if r is not None and (best is None or r[0] > best[0]):
                best = (r[0], r[1], (prow, pcol))
    return best


def _assert_sector_method_validates():
    """The sector global-slowest reproduces the full 4^N L slowest (C4, C6)."""
    for N in (4, 6):
        g = 1.0 / 2.0
        re_s, im_s, sec = global_slowest(N, 1.0, g, ring=True)
        sf, _ = slowest_full(N, 1.0, g, ring=True)
        assert abs(re_s - sf.real.max()) < 1e-4, f"C{N}: sector Re {re_s} != full-L {sf.real.max()}"
        print(f"[1] cross-check C{N} ring Q=2: sector Re={re_s:+.5f} |Im|={im_s:.3f} sec={sec} == full-L Re={sf.real.max():+.5f}")


def _assert_seam_is_ring_universal():
    """At strong dephasing (Q=1.5) EVERY ring hands over to a FROZEN double-excitation mode (a diagonal
    (a,a) sector with 2 <= a <= N-2). The seam is the C=0.5 half-filling boundary, not aromaticity."""
    for N in (4, 6, 8):
        re, im, (a, b) = global_slowest(N, 1.0, 1.0 / 1.5, ring=True)
        assert im < 1e-6, f"C{N} Q=1.5: slowest should be FROZEN, got |Im|={im}"
        assert a == b and 2 <= a <= N - 2, f"C{N} Q=1.5: slowest sector {(a, b)} is not double-excitation (a,a)"
        print(f"[2] ring-universal seam, C{N} Q=1.5: FROZEN double-excitation sector {(a, b)} (Re={re:+.4f})")


def _assert_aromaticity_refuted():
    """At weak dephasing (Q=10) the discriminant: C4 (anomaly) is FROZEN double-excitation, while C6
    AND C8 are the OSCILLATING band-edge survivor (0,1). The 4n anti-aromatics C4 and C8 do NOT
    group -> aromaticity (4n+2 vs 4n) is NOT the discriminant."""
    res = {}
    for N in (4, 6, 8):
        re, im, sec = global_slowest(N, 1.0, 1.0 / 10.0, ring=True)
        res[N] = (re, im, sec)
        kind = "OSCILLATING band-edge" if im > 1e-6 else "FROZEN double-excitation"
        print(f"[3] aromaticity test, C{N} Q=10: {kind} sector={sec} |Im|={im:.3f}")
    assert res[4][1] < 1e-6 and res[4][2] == (2, 2), "C4 should be FROZEN (2,2) at weak dephasing (the anomaly)"
    assert res[6][1] > 1.9 and res[6][2] == (0, 1), "C6 should be the OSCILLATING band-edge (0,1)"
    assert res[8][1] > 1.9 and res[8][2] == (0, 1), "C8 should behave like C6, NOT like C4 (band-edge survivor)"
    print("    -> C4 (4n) != C8 (4n): the anti-aromatics do NOT group; aromaticity is REFUTED.")
    print("       The seam tracks ring + half-filling (C=0.5), not Hueckel's 4n+2. C4 is a small-ring anomaly.")


if __name__ == "__main__":
    for N in (4, 6, 8):
        mos = huckel_mos(N, ring=True)
        print(f"  C{N} ring MOs {np.round(mos, 2)} (zero non-bonding: {int(np.sum(np.abs(mos) < 1e-9))})")
    print()
    _assert_sector_method_validates()
    _assert_seam_is_ring_universal()
    _assert_aromaticity_refuted()
    print("\nAll asserts passed: the V-Effect seam is the ring/half-filling C=0.5 boundary (universal),")
    print("NOT aromaticity (C4 != C8 refutes 4n+2); C4 is a small-ring anomaly.")
