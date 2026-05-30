"""
The Frost circle as the face of the open-system clock (closes docs/carbon README Q5).

A chemist draws a conjugated ring's pi-MO energies as the vertices of an inscribed
N-gon on a circle of radius 2|beta| (Frost-Musulin, 1953). The R=CPsi^2 clock
(e^{lambda t} = e^{-alpha t} * e^{i omega t}) is also a circle, but a running, open-system
one: under Holstein phonon dephasing each pi-coherence winds inward (radius -> lifetime
tau = 1/(2 gamma)) while it turns (angle -> the Bohr frequency). The static Frost circle
is the chemist's snapshot; the clock is the same circle in motion, and it adds two things
the static picture does not carry:

  1. WHICH coherence outlives the rest, and how long. The longest-lived pi-coherence is
     the band-edge one, beating at omega = 2|beta|*cos(pi/(N+1)) for an N-carbon polyene
     chain (2|beta| for the benzene ring), with lifetime tau = 1/(2 gamma_phonon).
  2. A crossover Q* = J/gamma (coherent <-> incoherent). Below Q* the band-edge coherence
     stops beating (the slowest Liouvillian mode becomes pure decay); above it, it beats.
     Q* grows with chain length: longer conjugated systems tolerate less dephasing before
     their slowest pi-coherence freezes.

Carbon map (docs/carbon README, resolved 2026-05-22): each carbon pi-site is a qubit
(occupied/empty), Huckel hopping beta = the XX+YY coupling J, Holstein phonon = the
Z-dephasing rate gamma. Energies are in units of |beta| = J.
"""
import sys
import math
import warnings
sys.path.insert(0, 'simulations')
import numpy as np
import framework as fw

warnings.filterwarnings('ignore')  # ChainSystem(N=2) structural-degeneracy notice
GAMMA = 0.05


def huckel_mos(N, topology, J=1.0):
    """The single-particle pi-MO energies (the Frost-circle / Huckel spectrum): the N x N
    tight-binding matrix, nearest-neighbour hopping J, ring closes the last bond."""
    A = np.zeros((N, N))
    for i in range(N - 1):
        A[i, i + 1] = A[i + 1, i] = J
    if topology == 'ring':
        A[0, N - 1] = A[N - 1, 0] = J
    return np.sort(np.linalg.eigvalsh(A))


def clock(N, J, topology, gamma=GAMMA):
    """(gap, omega_mem) of the open-system clock: gap = slowest decay rate (= 2 gamma,
    tau = 1/gap), omega_mem = the band-edge pi-coherence frequency the memory hand reads."""
    cs = fw.ChainSystem(N=N, gamma_0=gamma, J=J, topology=topology, H_type='xy')
    ev = np.linalg.eigvals(cs.L)
    rate = -ev.real
    om = np.abs(ev.imag)
    nz = rate > 1e-9
    gap = float(rate[nz].min())
    omega = float(om[np.abs(rate - gap) <= 1e-6].max())
    return gap, omega


# --- 1. Benzene ring: the longest-lived coherence sits at the Frost-circle radius -------
print("1. Benzene C6 ring (the Frost circle, energies in |beta|):")
mos = huckel_mos(6, 'ring')
gap, om = clock(6, 1.0, 'ring')
print(f"   pi-MO energies: {np.round(mos, 3)}   (the inscribed hexagon on radius 2|beta|)")
print(f"   longest-lived pi-coherence: omega_mem = {om:.4f} = 2|beta| = the Frost radius")
print(f"   lifetime tau = 1/(2 gamma) = {1.0 / gap:.1f}   (the radial hand the static circle lacks)")
print()

# --- 2. Polyene chains: omega_mem = 2|beta| cos(pi/(N+1)), the top pi-MO (band edge) ----
print("2. Polyene chains (open):  the band-edge coherence omega_mem = 2|beta| cos(pi/(N+1))")
print("   N    omega_mem      2 cos(pi/(N+1))   top |MO|     |diff|")
for N in (4, 5, 6):
    gap, om = clock(N, 1.0, 'chain')
    band = 2.0 * math.cos(math.pi / (N + 1))
    top_mo = float(np.max(np.abs(huckel_mos(N, 'chain'))))
    print(f"   {N}   {om:.6f}      {band:.6f}      {top_mo:.6f}    {abs(om - band):.1e}")
print()

# --- 3. The crossover Q* (coherent <-> incoherent), grows with chain length -------------
print("3. Crossover Q* = J/gamma  (below: band-edge coherence frozen; above: beating):")


def omega_mem(N, Q):
    _, om = clock(N, Q * GAMMA, 'chain')
    return om


def q_star(N):
    lo, hi = 0.3, 5.0
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if omega_mem(N, mid) > 1e-6:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


for N in (3, 4, 5):
    print(f"   N={N}:  Q* = {q_star(N):.3f}")
print("   Q* grows with N (sqrt(2) at N=3): the longer the conjugated chain, the less")
print("   phonon dephasing it tolerates before its slowest pi-coherence stops beating.")
print("   Real conjugated systems sit at Q ~ 100 (room T), deep in the beating regime;")
print("   Q* is the strong-dephasing edge, and it tightens with conjugation length.")
