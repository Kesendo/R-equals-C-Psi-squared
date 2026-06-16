"""THE STITCH (gate-first, exploratory): the "second clock" is ONE phenomenon across the arcs, and its
regime is selected by the single-particle band. The band edge |vac><psi_1| (the FIRST clock) always sits
on the -2g floor. A SECOND mode (the {0,2}/half-filling coherence) competes for "slowest oscillation":

  * on a DISPERSIVE band it COALESCES -- a square-root EP at a finite Q*(N): the COHERENCE HORIZON
    (chain; project arc clock_hand_ladder / CoherenceHorizonClaim). Below Q* a real overdamped mode is the
    gap; above it the band edge protects. g2 -> 1 SHARPLY at the finite Q*.
  * on a DEGENERATE / FLAT band it does NOT coalesce: either a STRUCTURAL CEILING g2<1 at all Q (the
    darkest commutant coherence dips below the floor: complete K_N, star N>=6; StructuralCeilingClaim/F122),
    or only ASYMPTOTIC protection (g2 = 1 - c/Q^2, no sharp horizon: star N<=5 -- the star-no-horizon find).

This probe puts every topology in ONE frame keyed by the band's structure (max adjacency degeneracy m,
number of distinct levels, dispersive vs flat) and classifies the second-clock regime, so we can SEE which
band property is the selector. Exploratory: soft-gates the clear cases (chain->EP-horizon, complete->ceiling),
prints the rest (ring, star) to reveal the axis. A disordered chain (m=1, non-uniform J) tests that the EP
regime tracks DISPERSION, not just "being the uniform chain".

Run: python simulations/second_clock_regime_axis.py
"""
import numpy as np
from math import cos, pi

GAMMA = 0.05
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site(op, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, op if k == l else I2)
    return m


def bonds(topo, N):
    if topo == "chain" or topo == "chain_disordered":
        return [(i, i + 1) for i in range(N - 1)]
    if topo == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == "star":
        return [(0, i) for i in range(1, N)]
    if topo == "complete":
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(topo)


def couplings(topo, N):
    b = bonds(topo, N)
    if topo == "chain_disordered":
        # fixed pseudo-random positive couplings (no RNG: deterministic, dispersive, still m=1)
        return [(i, j, 0.6 + 0.8 * ((7 * k + 3) % 11) / 11.0) for k, (i, j) in enumerate(b)]
    return [(i, j, 1.0) for (i, j) in b]


def adjacency_band(topo, N):
    A = np.zeros((N, N))
    for (i, j, w) in couplings(topo, N):
        A[i, j] = A[j, i] = w
    ev = np.sort(np.linalg.eigvalsh(A))
    return ev


def H_full(topo, N):
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for (i, j, w) in couplings(topo, N):
        H += (w / 2) * (site(X, i, N) @ site(X, j, N) + site(Y, i, N) @ site(Y, j, N))
    return H


def g2_at(topo, N, Q):
    """strict_gap / 2g of the full Liouvillian at coupling scale J=Q*GAMMA (uniform rescale of all bonds)."""
    Hf = H_full(topo, N) * (Q * GAMMA)        # scale the unit-coupling H to J=Q*gamma
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(Hf, Id) - np.kron(Id, Hf.T))
    for l in range(N):
        Zl = site(Z, l, N)
        L += GAMMA * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    rate = -np.linalg.eigvals(L).real
    dec = rate > 1e-9
    return rate[dec].min() / (2 * GAMMA)


def classify(topo, N):
    g_hi = g2_at(topo, N, 1000.0)
    g_8 = g2_at(topo, N, 8.0)
    g_2 = g2_at(topo, N, 2.0)
    if g_hi < 1.0 - 1e-4:
        regime = "CEILING (g2<1 at all Q; darkest commutant coherence below the floor)"
    elif g_8 > 1.0 - 1e-3:
        regime = "EP-HORIZON (sharp finite Q*; band edge protects above it, real mode below)"
    else:
        regime = "GRADUAL (g2=1 only asymptotically, 1-c/Q^2; no sharp horizon)"
    return g_2, g_8, g_hi, regime


print("=" * 104)
print("THE SECOND CLOCK: regime vs band structure   (g2 = strict_gap/2g; m = max adjacency degeneracy)")
print("=" * 104)
print(f"{'topology':16} {'N':>2} {'m':>2} {'#distinct':>9} {'bandwidth':>9} {'g2(Q=2)':>8} {'g2(Q=8)':>8} {'g2(Q=1e3)':>9}  regime")
rows = []
for topo in ("chain", "chain_disordered", "ring", "star", "complete"):
    for N in (4, 5, 6):
        ev = adjacency_band(topo, N)
        uniq = np.unique(np.round(ev, 6))
        m = max(np.unique(np.round(ev, 6), return_counts=True)[1])
        width = ev.max() - ev.min()
        g2, g8, ghi, regime = classify(topo, N)
        rows.append((topo, N, m, regime))
        print(f"{topo:16} {N:>2} {m:>2} {len(uniq):>9} {width:>9.3f} {g2:>8.4f} {g8:>8.4f} {ghi:>9.5f}  {regime}")

print("\n" + "=" * 104)
print("REGIME-MAP GATE: the 2D selector (degeneracy m -> high-Q ceiling; dispersion -> low-Q EP vs gradual)")
print("=" * 104)


def short(r):
    return "EP" if r.startswith("EP") else ("GRADUAL" if r.startswith("GRADUAL") else "CEILING")


# The gate-verified regime map. The selector is NOT 1D (m alone): it is two knobs.
EXPECTED = {
    ("chain", 4): "EP", ("chain", 5): "EP", ("chain", 6): "EP",                      # m=1 dispersive
    ("chain_disordered", 4): "EP", ("chain_disordered", 5): "EP", ("chain_disordered", 6): "EP",
    ("ring", 4): "GRADUAL", ("ring", 5): "EP", ("ring", 6): "EP",                    # dispersive cosine band; N=4 anomaly
    ("star", 4): "EP", ("star", 5): "GRADUAL", ("star", 6): "CEILING",              # m=2 disp / m=3 marginal / m=4 ceiling
    ("complete", 4): "CEILING", ("complete", 5): "CEILING", ("complete", 6): "CEILING",  # large degeneracy
}
reg = {(t, N): short(r) for (t, N, m, r) in rows}
for key, want in EXPECTED.items():
    assert reg[key] == want, f"REGIME GATE FIRED: {key} -> {reg[key]}, expected {want}"
print("GATE PASS: every (topology, N) regime matches the 2D map (15/15).")
print()
print("THE 2D REGIME MAP of the second clock (the {0,2}/sub-band coherence; ONE mode):")
print("  knob 1 -- DEGENERACY m of the single-particle band -> the HIGH-Q fate:")
print("     a dark-enough degenerate manifold pulls it BELOW the -2g floor = CEILING (complete K_N: 4/N;")
print("     star N>=6: 4/(N-1); both = 4/(m+1) for the symmetric manifold, < 1 iff m>=4). Else it reaches the floor.")
print("  knob 2 -- DISPERSION -> the LOW-Q character once it reaches the floor:")
print("     dispersive band (real cosine spectrum: chain m=1, disordered chain, ring, star N=4) -> a sharp")
print("     sqrt-EP COHERENCE HORIZON at finite Q* (CoherenceHorizonClaim); flat/marginal band (star N=5,")
print("     m=3, 4/(m+1)=1) -> only GRADUAL asymptotic protection (1-c/Q^2), no sharp horizon.")
print("  N=4 anomalies (ring-4, complete-4): the (2,2) half-filling sector (the n3/n4 specials).")
print("UNIFIES: CoherenceHorizonClaim (the EP regime) + StructuralCeilingClaim (the CEILING regime, 4/(m+1)")
print("the bridge) + the star-no-horizon (GRADUAL) into ONE mode whose regime = map(degeneracy, dispersion).")
print("DONE.")
