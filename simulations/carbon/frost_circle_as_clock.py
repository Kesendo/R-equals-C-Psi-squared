"""
The Frost circle as the face of the open-system clock (closes docs/carbon README Q5).

Frost's circle is a static Hückel-spectrum construction.  This script separately
computes selected open XX+YY models: their band-edge frequency, dephasing lifetime,
and crossover Q* = J/gamma.  The two calculations share a hopping-style algebraic
comparison; they do not make beta = J or identify a carbon degree of freedom, bath,
gamma, T2, or Q.

Within the selected model, the band-edge frequency is 2J*cos(pi/(N+1)) for an open
chain and 2J for the C6 ring.  The Q* ladder is likewise a selected-model result.
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
    tau = 1/gap), omega_mem = the band-edge pi-coherence frequency the memory hand reads.

    This is a small-N selected-model dense full-Liouvillian check, not a scalable or
    independent validation of the C6 handover or the general clock law.
    """
    cs = fw.ChainSystem(N=N, gamma_0=gamma, J=J, topology=topology, H_type='xy')
    ev = np.linalg.eigvals(cs.L)
    rate = -ev.real
    om = np.abs(ev.imag)
    nz = rate > 1e-9
    gap = float(rate[nz].min())
    omega = float(om[np.abs(rate - gap) <= 1e-6].max())
    return gap, omega


# --- 1. C6 ring: structural Frost comparison and selected-model clock -------------------
print("1. C6 ring: Frost spectrum and selected XX+YY-clock comparison:")
mos = huckel_mos(6, 'ring')
gap, om = clock(6, 1.0, 'ring')
print(f"   Hückel levels / |beta|: {np.round(mos, 3)}   (the inscribed hexagon)")
print(f"   selected XX+YY band-edge clock: omega_mem = {om:.4f} = 2J")
print(f"   selected-model lifetime tau = 1/(2 gamma) = {1.0 / gap:.1f}")
print()

# --- 2. Selected open XY chains ----------------------------------------------------------
print("2. Selected open XY chains: omega_mem = 2J cos(pi/(N+1))")
print("   N    omega_mem      2 cos(pi/(N+1))   top |MO|     |diff|")
for N in (4, 5, 6):
    gap, om = clock(N, 1.0, 'chain')
    band = 2.0 * math.cos(math.pi / (N + 1))
    top_mo = float(np.max(np.abs(huckel_mos(N, 'chain'))))
    print(f"   {N}   {om:.6f}      {band:.6f}      {top_mo:.6f}    {abs(om - band):.1e}")
print()

# --- 3. Selected-model crossover Q* (coherent <-> incoherent) ---------------------------
print("3. Selected-model crossover Q* = J/gamma  (below: a silent relaxation outlives the protected beat; above: beating leads):")


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
    print(f"   N={N}:  selected-model Q* = {q_star(N):.3f}")
print("   Q* grows with N (sqrt(2) at N=3): in this selected model, a longer chain")
print("   tolerates less chosen dephasing before a silent relaxation overtakes the (protected)")
print("   band-edge beat as the longest-lived mode (the longest-lived ROLE swaps; the beat")
print("   itself never stops, mechanism sharpened 2026-06-13).")

# --- 4. Selected C6 ring crossover -------------------------------------------------------
print()
print("4. Selected C6 ring crossover (the half-filled even-N ring sits on the V-Effect seam):")
print("   the selected-model single-excitation erasure point (Uhr 2) is Q* = 1.609; the selected full-L")
print("   band-edge handover is a DOUBLE-excitation mode near Q ~ 1.95. The two SPLIT here, though they coincide for the")
print("   open chains. Computed and asserted in simulations/carbon/benzene_two_clocks.py.")

# --- self-check: selected-model ladder and C6 band edge ----------------------------------
_LADDER = {3: 2.0 ** 0.5, 4: 1.87874, 5: 2.37367}
for _n, _q in _LADDER.items():
    _got = q_star(_n)
    assert abs(_got - _q) < 5e-3, f"chain Q*({_n})={_got:.4f} != canonical {_q:.4f}"
assert abs(clock(6, 1.0, 'ring')[1] - 2.0) < 1e-6, "C6 selected-ring hand omega_mem != 2J"
print()
print("   [self-check] selected-chain ladder Q*(3,4,5) and C6 omega_mem = 2J reproduced.")
