"""
The Frost circle as the face of the open-system clock (closes docs/carbon README Q5).

Frost's circle is a static Hückel-spectrum construction.  This script separately
computes selected open XX+YY models: their band-edge frequency, dephasing lifetime,
and crossover Q_h = J/gamma.  The two calculations share a hopping-style algebraic
comparison; they do not make beta = J or identify a carbon degree of freedom, bath,
gamma, T2, or Q.

Within the selected model, the band-edge frequency is 2J*cos(pi/(N+1)) for an open
chain and 2J for the C6 ring.  The Q_h ladder is likewise a selected-model result; it is the
HANDOVER, not the single-excitation EP Q*, which sits just above it from N=4 (see the self-check).
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

# --- 3. Selected-model crossover Q_h (coherent <-> incoherent) --------------------------
print("3. Selected-model crossover Q_h = J/gamma  (below: a silent relaxation outlives the protected beat; above: beating leads):")


def omega_mem(N, Q):
    _, om = clock(N, Q * GAMMA, 'chain')
    return om


def q_handover(N):
    lo, hi = 0.3, 5.0
    if omega_mem(N, hi) <= 1e-6:                  # Q_h grows ~0.59 N, so the bracket runs out
        raise ValueError(f"q_handover({N}): the handover is above the bracket end {hi}")
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if omega_mem(N, mid) > 1e-6:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


for N in (3, 4, 5):
    print(f"   N={N}:  selected-model handover Q_h = {q_handover(N):.5f}")
print("   Q_h grows with N (sqrt(2) at N=3): in this selected model, a longer chain")
print("   tolerates less chosen dephasing before a silent relaxation overtakes the (protected)")
print("   band-edge beat as the longest-lived mode (the longest-lived ROLE swaps; the beat")
print("   itself never stops, mechanism sharpened 2026-06-13).")

# --- 4. Selected C6 ring crossover -------------------------------------------------------
print()
print("4. Selected C6 ring crossover (the ring sits on the V-Effect seam through its 2-excitation")
print("   (2,2)/(N-2,N-2) doublet, not through half filling):")
print("   the selected-model single-excitation erasure point (Uhr 2) is Q* = 1.609; the selected full-L")
print("   band-edge handover is a TWO-excitation (2,2) mode at Q_h = 2.0000000. The ring SPLITS them")
print("   widely; on the open chains they coincide at N=2,3 and separate from N=4 by 2e-4 / 1.5e-3.")
print("   The 2.0000000 is a bisection over the joint-popcount sectors; benzene_two_clocks.py")
print("   asserts the split itself (frozen (2,2)/(4,4) below, Uhr 1 above).")

# --- self-check: selected-model ladder and C6 band edge ----------------------------------
# The bisection above measures the HANDOVER Q_h (where the FULL Liouvillian's slowest mode stops
# oscillating), not the single-excitation EP Q*. They are equal at N=2,3 and separate from N=4
# (Q* = 1.87874 / 2.37367), by 2e-4 / 1.5e-3.
#
# The error law of this estimator is NOT the bisection's own 4.7/2**24 = 2.8e-7. The binding term
# is the rate-matching window in clock(), |rate - gap| <= 1e-6: the gap approaches the 2*gamma
# floor LINEARLY in Q (slope ~2.5 at N=4), so the discriminator flips early by a fraction of that
# window, and the measured residual against a window-free solve is 1.7e-7 at N=4 and 8.4e-7 at
# N=5, i.e. of order the window and three times the bisection term. Widening the window to 1e-4
# moves Q_h(4) by -2.8e-5 and turns the gate below red, which is the mutation this gate is for.
# 1e-5 therefore bounds the residual with about a factor 10 to spare while still rejecting the EP
# at both rungs (2.0e-4 at N=4, 1.5e-3 at N=5, 20x and 150x the window).
_LADDER = {3: 2.0 ** 0.5, 4: 1.878541, 5: 2.372175}
_EP = {4: 1.87874, 5: 2.37367}
_MEASURED = {}
for _n, _q in _LADDER.items():
    _MEASURED[_n] = q_handover(_n)
    _msg = f"chain Q_h({_n})={_MEASURED[_n]:.6f} != {_q:.6f}"
    assert abs(_MEASURED[_n] - _q) < 1e-5, _msg
# Anti-vacuity, on the MEASUREMENT: the same comparison must REJECT the EP at both rungs, so a
# run that returned the EP instead of the handover cannot pass the loop above.
for _n in _EP:
    assert abs(_MEASURED[_n] - _EP[_n]) > 1e-5, f"the N={_n} rung no longer rejects the EP {_EP[_n]}"
assert abs(clock(6, 1.0, 'ring')[1] - 2.0) < 1e-6, "C6 selected-ring hand omega_mem != 2J"
print()
print("   [self-check] selected-chain handover Q_h(3,4,5) reproduced to 1e-5, the EP rejected at")
print("                N=4 and N=5, and C6 omega_mem = 2J.")
