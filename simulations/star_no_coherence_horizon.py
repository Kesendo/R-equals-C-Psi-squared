"""Gate-first: the STAR has NO chain-like coherence horizon (topology_band_edge arc, the "star Q*(N) for
N<=5" follow-on -- RESOLVED as a null + unification, not a closed form).

The chain coherence horizon Q*(N) is a DISPERSION effect: the {0,2}-coherence pair of a dispersive band
coalesces into a square-root EP (Q*(N) ~ 2N/pi; CoherenceHorizonClaim). It reduces to the single-excitation
(Haken-Strobl) Liouvillian -- the chain's SE-EP reproduces the full-L horizon exactly. The star's
single-particle band is FLAT (adjacency spectrum +-sqrt(N-1) and 0 with multiplicity N-2): no dispersion,
no coalescence, no EP. The chain's SE-EP mechanism does NOT port.

Instead the star's band-edge protection at EVERY Q is governed by the SAME (1,1)-commutant value 4/(N-1)
that sets the high-Q structural ceiling (StructuralCeilingClaim / F122):
   4/(N-1) > 1  (N=4):           the sub-band mode is above the 2g floor -> the band edge protects down to a
                                 low-Q real-mode CROSSING (a real mode crosses the floor; NOT an EP);
   4/(N-1) = 1  (N=5, marginal): the sub-band mode sits exactly on the floor and is pulled marginally below
                                 at finite Q as g2 = 1 - 1/Q^2 -> NO horizon, no oscillating memory, only
                                 asymptotic protection (the apparent "Q*~316" is a tolerance artifact);
   4/(N-1) < 1  (N>=6):          the structural ceiling g2 = 4/(N-1).
   N=3 = the path P_3 = a chain (the lone exception: it IS a chain, the genuine {0,2}-EP at sqrt2).

So the low-Q question is SUBSUMED by the high-Q ceiling; there is no separate star Q*(N) closed form.

  STAGE 0  port fidelity: the chain SE-EP reproduces Q*(3,4,5) = sqrt2, 1.8787, 2.3737.
  STAGE 1  the gate: the star SE-EP is SPURIOUS -- it predicts a horizon (N=4: Q~261) that the full L does
           not show (full-L star N=4 is already protected far below it). The chain mechanism does not port.
  STAGE 2  governance by 4/(N-1): full-L star high-Q g2 = min(1, 4/(N-1)); N=5 is marginal ((1-g2)*Q^2 -> 1,
           no horizon); N>=6 is the ceiling 4/(N-1).
  STAGE 3  the dichotomy: the star sub-band gap mode is REAL (overdamped) when g2<1 -- no coalescence, a
           crossing, not the chain's oscillating-pair EP.

Run: python simulations/star_no_coherence_horizon.py
"""
import numpy as np
from math import sqrt, cos, pi

GAMMA = 0.05


# ----------------------------------------------------------------- single-particle Hamiltonians
def h_chain(N, J=1.0):
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = J
    return h


def h_star(N, J=1.0):
    h = np.zeros((N, N), complex)
    for i in range(1, N):
        h[0, i] = h[i, 0] = J
    return h


def rho_adj(topo, N):
    return sqrt(N - 1) if topo == "star" else 2.0 * cos(pi / (N + 1))


# ----------------------------------------------------------------- SE Haken-Strobl Liouvillian + EP
def L_se(h, g):
    N = h.shape[0]
    I = np.eye(N)
    L = -1j * (np.kron(h, I) - np.kron(I, h.T))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return L + np.diag(deph)


def qstar_se(hfunc, N, J=1.0, lo=1e-4, hi=12.0):
    """Q* = J/g*; g* = largest g where the slowest nonzero SE mode still oscillates (the EP)."""
    for _ in range(90):
        m = 0.5 * (lo + hi)
        ev = np.linalg.eigvals(L_se(hfunc(N, J), m))
        nz = ev[ev.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        gap = nz.real.max()
        band = nz[np.abs(nz.real - gap) < 1e-7]
        if np.abs(band.imag).max() > 1e-7:
            lo = m
        else:
            hi = m
    return J / (0.5 * (lo + hi))


# ----------------------------------------------------------------- full Liouvillian (ground truth)
_X = np.array([[0, 1], [1, 0]], complex)
_Y = np.array([[0, -1j], [1j, 0]])
_Z = np.diag([1, -1]).astype(complex)
_I2 = np.eye(2)


def _site(op, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, op if k == l else _I2)
    return m


def L_full(topo, N, J, g):
    bonds = ([(0, i) for i in range(1, N)] if topo == "star"
             else [(i, i + 1) for i in range(N - 1)])
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for (a, b) in bonds:
        H += (J / 2) * (_site(_X, a, N) @ _site(_X, b, N) + _site(_Y, a, N) @ _site(_Y, b, N))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = _site(_Z, l, N)
        L += g * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def full_clock(topo, N, Q):
    """(g2, omega_mem/band, gap_is_real) from the full L at coupling J=Q*GAMMA."""
    J = Q * GAMMA
    ev = np.linalg.eigvals(L_full(topo, N, J, GAMMA))
    rate = -ev.real
    dec = rate > 1e-9
    gap = rate[dec].min()
    at_gap = dec & (np.abs(rate - gap) < 1e-6)
    om = np.abs(ev.imag[at_gap]).max()
    band = J * rho_adj(topo, N)
    return gap / (2 * GAMMA), (om / band if band > 0 else 0.0), om < 1e-4 * band


# ================================================================= STAGE 0 -- port fidelity (chain)
print("=" * 96)
print("STAGE 0 -- PORT FIDELITY: the SE-EP reproduces the chain horizon ladder")
print("=" * 96)
ladder = {2: 1.0, 3: sqrt(2.0), 4: 1.87874, 5: 2.37367}
for N, q in ladder.items():
    got = qstar_se(h_chain, N)
    assert abs(got - q) < 3e-3, f"STAGE 0 FAILED: chain Q*({N})={got:.5f} != {q:.5f}"
print("  chain Q* =", {N: round(qstar_se(h_chain, N), 5) for N in ladder}, "(the {0,2} dispersion EP)")
print("STAGE 0 PASS: the SE-EP harness reproduces the chain coherence horizon.")

# ================================================================= STAGE 1 -- the SE-EP is spurious for the star
print("\n" + "=" * 96)
print("STAGE 1 -- THE GATE: the star SE-EP is SPURIOUS (the chain dispersion mechanism does not port)")
print("=" * 96)
se4 = qstar_se(h_star, 4)
g2_star4_at20, _, _ = full_clock("star", 4, 20.0)
print(f"  star N=4: SE-EP predicts a horizon at Q~{se4:.1f}, but the full L is already PROTECTED at Q=20 "
      f"(g2={g2_star4_at20:.5f}).")
assert se4 > 50.0, f"STAGE 1: expected a large spurious star SE-EP, got {se4:.2f}"
assert abs(g2_star4_at20 - 1.0) < 1e-4, f"STAGE 1: star N=4 should be protected at Q=20, g2={g2_star4_at20:.5f}"
print("STAGE 1 PASS: the SE-EP horizon is an artifact of the flat band; the chain mechanism does not port.")

# ================================================================= STAGE 2 -- governance by 4/(N-1)
print("\n" + "=" * 96)
print("STAGE 2 -- GOVERNANCE BY 4/(N-1): full-L star high-Q g2 = min(1, 4/(N-1)); N=5 marginal, N>=6 ceiling")
print("=" * 96)
print(f"  {'N':>2} {'4/(N-1)':>8} {'min(1,.)':>9} {'g2(Q=1000)':>11} {'(1-g2)*Q^2':>11}  regime")
for N in (4, 5, 6):
    g2, ob, _ = full_clock("star", N, 1000.0)
    pred = min(1.0, 4.0 / (N - 1))
    reg = ("protected (band edge)" if 4.0 / (N - 1) > 1 + 1e-9
           else "marginal: g2=1-1/Q^2, NO horizon" if abs(4.0 / (N - 1) - 1) < 1e-9
           else f"structural ceiling 4/{N-1}")
    print(f"  {N:>2} {4.0/(N-1):>8.4f} {pred:>9.4f} {g2:>11.6f} {(1-g2)*1000*1000:>11.4f}  {reg}")
    assert abs(g2 - pred) < 2e-3, f"STAGE 2: star N={N} g2={g2:.5f} != min(1,4/(N-1))={pred:.5f}"
# N=5 marginal: (1-g2)*Q^2 -> 1 (the 1/Q^2 approach, NOT a horizon)
for Q in (100.0, 1000.0):
    g2, _, _ = full_clock("star", 5, Q)
    c = (1 - g2) * Q * Q
    assert abs(c - 1.0) < 0.05, f"STAGE 2: star N=5 (1-g2)*Q^2={c:.4f} should -> 1 (marginal 1/Q^2)"
print("STAGE 2 PASS: 4/(N-1) governs the star at all Q. N=5 is marginal (1-g2=1/Q^2, no horizon).")

# ================================================================= STAGE 3 -- the dichotomy (real vs EP)
print("\n" + "=" * 96)
print("STAGE 3 -- THE DICHOTOMY: the star sub-band gap mode is REAL (overdamped) when g2<1 -- a crossing,")
print("           not the chain's oscillating-pair EP. (N=4 has a low-Q crossing; N=5 has none.)")
print("=" * 96)
# N=4 crossing exists: broken (g2<1) at Q=1.5, protected (g2=1) at Q=3
g2_4lo, _, real_4lo = full_clock("star", 4, 1.5)
g2_4hi, _, _ = full_clock("star", 4, 3.0)
print(f"  star N=4: g2(Q=1.5)={g2_4lo:.5f} (gap mode real={real_4lo}), g2(Q=3)={g2_4hi:.5f} -> a low-Q crossing in (1.5,3)")
assert g2_4lo < 0.99 and real_4lo, "STAGE 3: star N=4 should be broken with a REAL gap mode at Q=1.5"
assert abs(g2_4hi - 1.0) < 1e-4, "STAGE 3: star N=4 should protect by Q=3"
# N=5 no crossing: g2<1 with a real gap mode at every finite Q (asymptotes, never crosses)
for Q in (3.0, 20.0, 100.0):
    g2, _, real = full_clock("star", 5, Q)
    assert g2 < 1.0 and real, f"STAGE 3: star N=5 gap mode should be real and g2<1 at Q={Q} (g2={g2:.5f})"
print("  star N=5: gap mode REAL and g2<1 at Q=3,20,100 -- no oscillating coalescence, no crossing.")
print("STAGE 3 PASS: dispersion EP (chain) vs flat-band real-mode crossing/asymptote (star).")

print("\nSUMMARY: the star has no chain-like coherence horizon. The chain's is a dispersion EP; the star's")
print("flat band has none. One number, the (1,1)-commutant 4/(N-1), governs the star at every Q: >1 protects")
print("(N=4, down to a low-Q crossing), =1 marginal (N=5, g2=1-1/Q^2, no horizon), <1 ceilings (N>=6). The")
print("low-Q question is subsumed by the high-Q ceiling -- no separate Q*(N) closed form. (N=3 = path = chain.)")
print("DONE.")
