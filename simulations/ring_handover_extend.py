"""Extension of ring_handover_qh: one more Q_h point (N=14) for the slope trend, and the (2,2) seam
HIGH-Q darkness vs the parallel-session commutant closed form 2(N-2)/N (b191df3) -- a cross-session
convergence test (does my low-Q seam's high-Q limit meet Tom's high-Q commutant?)."""
import sys
import numpy as np
sys.path.insert(0, "simulations")
from ring_handover_qh import darkness, qh_sector

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

print("=" * 90)
print("(a) Q_h(k=2) extended to N=14 (shift-invert), for the asymptotic-slope trend")
print("=" * 90)
# N=6/8/10 by a root solve on the same darkness=1 criterion (N=6 is EXACTLY 2 = N/3);
# N=12/14 from qh_sector, which its own docstring caps at +-0.01 (it breaks at hi-lo<0.02).
qh = {6: 2.00000000, 8: 2.35038584, 10: 2.83519426, 12: 3.35977}
COARSE = {12, 14}                     # the rungs carrying the +-0.01, marked in the table
qh[14] = qh_sector(14, 2)
print(f"{'N':>4} {'Q_h':>10} {'Q_h/N':>9} {'local slope (dQ_h/dN)':>22}")
Ns = sorted(qh)
for i, N in enumerate(Ns):
    ls = "" if i == 0 else f"{(qh[N]-qh[Ns[i-1]])/(N-Ns[i-1]):>22.5f}"
    mark = " +-0.01" if N in COARSE else "       "
    print(f"{N:>4} {qh[N]:>10.5f}{mark} {qh[N]/N:>9.5f} {ls}")

print("\n" + "=" * 90)
print("(b) the (2,2) seam HIGH-Q darkness (Q=200) vs the commutant forms (cross-session: b191df3)")
print("    Tom's (1,1)/(N/2,N/2) commutant = 2(N-2)/N (even). Does my (2,2) seam meet it, or differ?")
print("=" * 90)
print(f"{'N':>4} {'(2,2) darkness(Q=200)':>22} {'2(N-2)/N':>10} {'2(N-1)/N':>10} {'match 2(N-2)/N?':>16}")
for N in (6, 8, 10, 12):
    dk = darkness(N, 2, 200.0)
    a = 2 * (N - 2) / N
    b = 2 * (N - 1) / N
    print(f"{N:>4} {dk:>22.5f} {a:>10.5f} {b:>10.5f} {str(abs(dk-a)<1e-2):>16}")
print("\n  (if the (2,2) seam high-Q darkness != 2(N-2)/N, the seam is a DIFFERENT object than Tom's")
print("   commutant -- the two sessions touch the same ring but track different modes; the handover")
print("   Q_h is the low-Q crossover, his g2=1 the high-Q ceiling, both consistent.)")
