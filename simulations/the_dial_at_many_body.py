"""
The dial at many-body: what the clock's memory hand actually does
(a closed form, a crossover, and a skipped band).

We asked whether the clock's slowest-mode angle walks through the element markings as
Q = J/gamma sweeps (DEPTH_3: theta -> alpha = sin^2(theta)/2 -> valences; 60 deg = Boron
3/8, 45 = Beryllium 1/4, 30 = Lithium 1/8, at the canonical Q = sqrt(3), 1, 1/sqrt(3)).
Reading the Heisenberg chain under Z-dephasing, the honest answer is more interesting
than a smooth sweep, and the clock corrected the clean story itself:

1. CLOSED FORM (real, bit-exact). The rate-2gamma rotating coherence turns at
       omega = 8 J cos^2(pi/2N)   [Pauli H = J*(XX+YY+ZZ), this framework]
             = 2 J cos^2(pi/2N)   [spin  H = (J/4)*(XX+YY+ZZ), the C# clock]
   with gap = 2 gamma, so theta = arctan(cos^2(pi/2N) * Q): a cos^2-geared arctan(Q),
   gear g_N = cos^2(pi/2N). Bulk N->inf: g->1, theta->arctan(Q) = the F95 2-level reading.

2. A CROSSOVER (the clean story breaks). The rotating coherence is the SLOWEST (memory)
   mode only above a threshold Q*. Below Q* an F3 sub-2gamma hybrid mode (pure decay,
   omega = 0) is slowest, so the memory hand is FROZEN at theta = 0. At N=3, Q* ~ 3.2.
   (F3 was noted at N>=4; the clock surfaces it at N=3 too, as a Q-threshold.)

3. A SKIPPED BAND (the point). At the crossover the memory hand jumps from 0 deg straight
   to arctan(g_N * Q*) ~ 69 deg (N=3) and climbs toward 90. The element markings 60/45/30
   sit in the gap the memory hand skips. The 60 deg coherence DOES exist (at Q = sqrt(3)/g_N,
   on the rate-2gamma band), but there it is an EXCITED mode, not the memory mode.

So the face (the DEPTH_3 markings) is real, but the memory (slowest) hand of a real chain
does not trace that arc, it skips it. The element-angle coherences live on the rate-2gamma
band as excited modes, below the crossover, not as the memory mode.

Convention bridge: the C# clock builds H = (J/4)(XX+YY+ZZ) (spin S = sigma/2); this Python
framework builds H = J(XX+YY+ZZ) (Pauli). Energies, hence omega, differ by exactly 4; the
gap (2gamma) does not. So CLI J = framework J * 4. The sanity check below reproduces the
CLI memory reading (N=3, CLI J=1, Q=20: omega 1.5 spin, 86.2 deg).
"""
import sys
import math
import warnings
sys.path.insert(0, 'simulations')
import numpy as np
import framework as fw

warnings.filterwarnings('ignore')  # ChainSystem(N=2) emits a structural-degeneracy notice
GAMMA = 0.05


def memory_mode(N, J_fw, gamma=GAMMA, gap_tol=1e-6):
    """(gap, |omega|) of the slowest nonstationary Liouvillian mode (the clock's memory
    hand). J_fw is the Pauli J of this framework (= CLI/spin J divided by 4)."""
    cs = fw.ChainSystem(N=N, gamma_0=gamma, J=J_fw, topology='chain', H_type='heisenberg')
    ev = np.linalg.eigvals(cs.L)
    rate = -ev.real
    om = np.abs(ev.imag)
    gap = float(rate[rate > 1e-9].min())
    omega = float(om[np.abs(rate - gap) <= gap_tol].max())
    return gap, omega


def cos2(N):
    return math.cos(math.pi / (2 * N)) ** 2


# --- 1. the closed form for the rate-2gamma rotating coherence, bit-exact ---------------
print("1. CLOSED FORM  omega = 8 J cos^2(pi/2N)  and  gap = 2*gamma   (Pauli H, full L,")
print("   read at J=1 / Q=80 where the rotating coherence is firmly the slowest mode)")
print("   N    gap      2*gamma    omega (J=1)        8*cos^2(pi/2N)     |diff|")
for N in range(2, 7):
    gap, omega = memory_mode(N, 1.0)
    closed = 8.0 * cos2(N)
    print(f"   {N}   {gap:.4f}   {2*GAMMA:.4f}    {omega:.12f}   {closed:.12f}   {abs(omega-closed):.1e}")
gap3, om3 = memory_mode(3, 0.25)  # CLI J=1 (Q=20) -> framework J=0.25
print(f"   sanity: N=3 at CLI J=1 (Q=20) -> omega {om3:.4f}, "
      f"theta {math.degrees(math.atan2(om3, gap3)):.1f} deg  [C# CLI: 1.5, 86.2]")
print("   gear g_N = cos^2(pi/2N) (spin); bulk N->inf: g->1, theta -> arctan(Q) = F95.")
print()

# --- 2. the crossover: below Q* the memory hand is frozen (F3 hybrid is slowest) ---------
print("2. CROSSOVER  Q* (CLI Q = J/gamma): below it the F3 sub-2gamma hybrid is slowest")
print("   (pure decay, theta=0); above it the rate-2gamma rotating coherence is slowest.")
print("   N    Q*       theta_mem just above Q*  (the bottom of the rotating band)")
crossover = {}
for N in (3, 4, 5):
    qstar, th_above = None, None
    Q = 0.5
    while Q <= 10.0:
        gap, om = memory_mode(N, (Q * GAMMA) / 4.0)   # CLI Q -> framework J = Q*gamma/4
        if om > 1e-6:                                  # the rotating coherence is now slowest
            qstar, th_above = Q, math.degrees(math.atan2(om, gap))
            break
        Q += 0.05
    crossover[N] = (qstar, th_above, cos2(N))
    print(f"   {N}    {qstar:.2f}     {th_above:5.1f} deg")
print()

# --- 3. the skipped band: every element-Q falls below Q*, in the frozen regime ----------
print("3. SKIPPED BAND  element-Q = tan(theta_canonical)/g_N (where the rotating coherence")
print("   would sit at the element angle). If element-Q < Q*, the memory hand is still")
print("   FROZEN there, so it never points at that element as its slowest mode.")
anchors = [("Boron 60", 60.0), ("Beryllium 45", 45.0), ("Lithium 30", 30.0)]
print("   N    Q*      Boron-Q  Be-Q   Li-Q     verdict")
for N in (3, 4, 5):
    qstar, _, g = crossover[N]
    qs = [math.tan(math.radians(d)) / g for _, d in anchors]
    verdict = "memory hand skips all 3" if all(q < qstar for q in qs) else "memory hand LANDS on one"
    print(f"   {N}    {qstar:.2f}    {qs[0]:5.2f}   {qs[1]:5.2f}  {qs[2]:5.2f}    {verdict}")
print()
print("The face (DEPTH_3 markings) is real; the memory (slowest) hand of a real chain does")
print("not trace that arc, it skips it. The element-angle coherences live on the rate-2gamma")
print("band as excited modes, below the crossover, not as the memory mode.")
