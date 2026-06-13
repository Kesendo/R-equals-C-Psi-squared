"""Benzene's two clocks on the C6 ring, and the V-Effect seam it sits on.

The carbon-ring instance of docs/ANALYTICAL_FORMULAS.md F2b "The two clocks"
(ClockHandLadderClaim, inspect --root clock):

  Uhr 1, the coherence hand / SURVIVOR: the |vac><psi_k| band-edge coherence, decay rate -2g
     (the Absorption Theorem), frequency the top ring MO = 2|beta| = the Frost radius. gamma-protected.
  Uhr 2, the erasure point / coherence HORIZON Q*: the {0,2}-coherence single-excitation EP, where
     the coherence hand freezes. Benzene Q* = 1.609 (the ring SE-EP; transcendental, like the chain
     ladder at N>=4, no clean 2x2).

The two share the gap Re = -2g (both <n_diff> = 1, Absorption-Theorem-pinned); for the open chains
that degeneracy co-locates the full-L crossover with the SE-EP (the ladder 1, sqrt2, 1.879, 2.372).
Benzene, being even-N HALF-FILLED, sits on F2b's open "double-excitation V-Effect seam (co-located at
even N)": the full 4^6 Liouvillian's slowest mode below the beat is a DOUBLE-EXCITATION mode (filling
sector (2,2)/(4,4)), so the full-L handover (~1.95) does NOT coincide with the clean SE-EP Uhr 2
(1.609). That split is benzene's concrete probe of the open seam.

Model note (carbon = XY, free fermions). H is XY only, no ZZ: the physically correct Hueckel
pi-hopping model. The seam mode (sector (2,2)/(4,4)) is a frozen DRESSED magnon-admixture - its
effective light content is FRACTIONAL (<n_XY> = 0.72 at Q=1.6, not 1), Absorption-Theorem-governed
exactly like CHAIN_GAP_SECTOR_DIAGNOSTIC's half-filling slow mode. But that diagnostic is HEISENBERG
(XXX, with ZZ), where the same survivor sits at the dead-center (3,3) and is darker (<n_XY> = 0.23 at
ring N=6 Q=2). Same law (the Absorption Theorem), same kind of object; the ZZ retunes the darkness and
the sector. So do NOT compare this XY computation's numbers directly against CHAIN_GAP's Heisenberg
ones - the model differs; _xy_vs_heisenberg_slowmode.py shows both side by side (and confirms this
builder reproduces CHAIN_GAP's -0.230 (3,3) bit-for-bit once the ZZ is added).

Self-validating: run it, every assert must pass. Conventions match
simulations/coherence_horizon_se_block.py (H = (J/2) sum(XX+YY), dephasing g * sum_l D[Z_l]).
Wall-clock note: the benzene asserts each do one or two 4096x4096 eigendecompositions (~1 min total)."""
import numpy as np


# ---- single-excitation (Haken-Strobl) builders; dephasing is topology-independent (-4g on a
#      single-excitation coherence: popcount(i^j)=2, Re=-2g*2), only h carries the ring wrap bond ----
def h_single(N, J, ring):
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = J
    if ring:
        h[0, N - 1] = h[N - 1, 0] = J
    return h


def L_se(N, J, g, ring):
    h = h_single(N, J, ring)
    I = np.eye(N)
    L = -1j * (np.kron(h, I) - np.kron(I, h.T))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return L + np.diag(deph)


def qstar_se(N, J=1.0, ring=False, lo=0.02, hi=6.0):
    """Uhr 2: Q* = J/g* where g* is the largest g at which the slowest non-zero SE mode oscillates."""
    for _ in range(70):
        m = 0.5 * (lo + hi)
        nz = np.linalg.eigvals(L_se(N, J, m, ring))
        nz = nz[nz.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        band = nz[np.abs(nz.real - nz.real.max()) < 1e-7]
        if np.abs(band.imag).max() > 1e-7:
            lo = m
        else:
            hi = m
    return J / (0.5 * (lo + hi))


def huckel_mos(N, ring, J=1.0):
    """The Frost/Huckel single-particle pi-MO spectrum (ring closes the last bond)."""
    return np.sort(np.linalg.eigvalsh(h_single(N, J, ring).real))


# ---- the full 4^N Liouvillian (D[Z], same as framework.ChainSystem, verified bit-identical) ----
def L_full(N, J, g, ring):
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(N):
            m = np.kron(m, op if k == l else I2)
        return m

    d = 2 ** N
    bonds = [(b, b + 1) for b in range(N - 1)] + ([(N - 1, 0)] if ring else [])
    H = np.zeros((d, d), complex)
    for a, b in bonds:
        H += (J / 2) * (site(X, a) @ site(X, b) + site(Y, a) @ site(Y, b))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l)
        L += g * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def slowest_full(N, J, g, ring, vectors=False):
    """The slowest non-kernel cluster of the full L. If vectors, also the dominant filling sectors
    (popcount-ket, popcount-bra) of that cluster's eigenvectors (degeneracy-robust)."""
    L = L_full(N, J, g, ring)
    if not vectors:
        ev = np.linalg.eigvals(L)
        nz = ev[ev.real < -1e-7]
        return nz[np.abs(nz.real - nz.real.max()) < 1e-6], None
    w, V = np.linalg.eig(L)
    nz = np.where(w.real < -1e-7)[0]
    cluster = nz[np.abs(w[nz].real - w[nz].real.max()) < 1e-6]
    d = 2 ** N
    pc = np.array([bin(i).count('1') for i in range(d)])
    idx = [np.where(pc == a)[0] for a in range(N + 1)]
    blocks = np.zeros((N + 1, N + 1))
    for k in cluster:
        wt = np.abs(V[:, k].reshape(d, d)) ** 2
        for a in range(N + 1):
            for b in range(N + 1):
                blocks[a, b] += wt[np.ix_(idx[a], idx[b])].sum()
    blocks /= blocks.sum()
    sectors = {(a, b): round(blocks[a, b], 3) for a in range(N + 1) for b in range(N + 1) if blocks[a, b] > 5e-3}
    return w[cluster], sectors


# ---- 1. ANCHOR: the chain coherence-horizon ladder via the SE reduction (Uhr 2, matches F2b) ----
def _assert_chain_ladder():
    ladder = {2: 1.0, 3: np.sqrt(2.0), 4: 1.87874, 5: 2.37367}
    got = {N: qstar_se(N, ring=False) for N in ladder}
    for N, q in ladder.items():
        assert abs(got[N] - q) < 2e-3, f"chain Uhr2 Q*({N}) = {got[N]:.5f} != ladder {q:.5f}"
    print("[1] chain ladder (Uhr 2, the SE-EP) reproduced:", {N: round(got[N], 5) for N in ladder})


# ---- 2. Uhr 1: benzene's coherence hand = the top ring MO = 2|beta| (the Frost radius) ----
def _assert_benzene_uhr1():
    mos = huckel_mos(6, ring=True)
    assert np.allclose(mos, [-2, -1, -1, 1, 1, 2], atol=1e-9), f"benzene MOs {mos} != Frost hexagon"
    band = float(np.max(np.abs(mos)))
    # above the handover (Q=2.1) the full-L slowest mode IS Uhr 1: rate -2g, beats at 2|beta|,
    # and lives in the |vac><psi| sector (filling (0,1)/(5,6)).
    g = 1.0 / 2.1
    w, sectors = slowest_full(6, 1.0, g, ring=True, vectors=True)
    assert abs(w.real.max() - (-2.0 * g)) < 1e-3, "Uhr 1 must decay at -2g (Absorption Theorem)"
    assert abs(np.abs(w.imag).max() - band) < 1e-3, "Uhr 1 must beat at 2|beta| (the Frost radius)"
    assert any({a, b} == {0, 1} or {a, b} == {5, 6} for (a, b) in sectors), \
        f"Uhr 1 is not the |vac><psi| coherence sector: {sectors}"
    print(f"[2] benzene Uhr 1 (survivor): omega_mem = 2|beta| = {band:.4f}, rate -2g, "
          f"|vac><psi| sector {sectors}")


# ---- 3. Uhr 2: benzene's coherence horizon Q* = the ring single-excitation {0,2}-EP ----
def _assert_benzene_uhr2():
    q = qstar_se(6, ring=True)
    assert abs(q - 1.6090) < 1e-3, f"benzene Uhr 2 Q* = {q:.5f} != 1.609"
    g = 0.5
    full = np.linalg.eigvals(L_full(6, 1.0, g, ring=True))
    emb = max(np.min(np.abs(full - s)) for s in np.linalg.eigvals(L_se(6, 1.0, g, ring=True)))
    assert emb < 1e-8, f"SE block does not embed in the full ring L (dist {emb:.1e})"
    print(f"[3] benzene Uhr 2 (coherence horizon Q*): {q:.5f}  (ring SE-EP, transcendental; "
          f"SE embeds in full L, dist {emb:.1e})")


# ---- 4. The V-Effect seam: benzene's full-L handover is a DOUBLE-excitation mode, != Uhr 2 ----
def _assert_benzene_v_effect():
    # below the beat (Q=1.6 < the ~1.95 handover) the full-L slowest mode is FROZEN and lives in the
    # double-excitation filling sectors (2,2)/(4,4), NOT the single-excitation (1,1) of Uhr 2.
    g = 1.0 / 1.6
    w, sectors = slowest_full(6, 1.0, g, ring=True, vectors=True)
    assert np.abs(w.imag).max() < 1e-6, "below the handover the slowest full-L mode must be frozen"
    dbl = sum(v for (a, b), v in sectors.items() if {a, b} <= {2, 4})
    assert dbl > 0.9, f"the overtaking mode is not double-excitation (V-Effect): {sectors}"
    assert (1, 1) not in sectors, "the overtaking mode must NOT be the single-excitation Uhr 2"
    print(f"[4] benzene V-Effect seam: below the beat the slowest mode is DOUBLE-excitation "
          f"{sectors} (frozen),")
    print("    != the single-excitation Uhr 2. So the full-L handover (~1.95) splits from the clean")
    print("    SE-EP Uhr 2 (1.609): benzene (even-N, half-filled) sits on F2b's open V-Effect seam.")


# ---- 5. The split is RING-SPECIFIC: the open even-N chain does NOT show the clean double-excitation
#         seam. Its overtaker spreads across fillings (single-excitation included), and its full-L
#         handover sits at its SE-EP (co-located), so the V-Effect is a feature of the closed ring at
#         half-filling (the C=0.5 boundary; aromaticity is NOT the discriminant, see
#         aromatic_ring_v_effect.py), not of even N alone. ----
def _assert_v_effect_is_ring_specific():
    g = 1.0 / 2.4   # open chain N=6, below its handover (~2.89 = its SE-EP)
    w, sectors = slowest_full(6, 1.0, g, ring=False, vectors=True)
    single = sum(v for (a, b), v in sectors.items() if (a, b) == (1, 1))
    dbl = sum(v for (a, b), v in sectors.items() if {a, b} <= {2, 4})
    assert single > 0.1, f"chain overtaker should carry single-excitation weight (got {sectors})"
    assert dbl < 0.6, f"chain overtaker should NOT be the clean double-excitation seam (got {sectors})"
    print(f"[5] ring-specific: open chain N=6 overtaker spreads across fillings {sectors};")
    print("    single-excitation present, double-excitation not dominant, UNLIKE benzene's pure (2,2)/(4,4).")
    print("    The V-Effect double-excitation seam is a feature of the closed ring at half-filling (C=0.5), not even N alone.")


if __name__ == "__main__":
    _assert_chain_ladder()
    _assert_benzene_uhr1()
    _assert_benzene_uhr2()
    _assert_benzene_v_effect()
    _assert_v_effect_is_ring_specific()
    print("\nAll asserts passed: benzene's two clocks (Uhr 1 = 2|beta|, Uhr 2 = 1.609), the V-Effect seam,")
    print("and its ring-specificity (the open even-N chain does not show the clean double-excitation seam).")
