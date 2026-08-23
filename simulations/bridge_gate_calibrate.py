"""Calibration gate: does the sector propagator reproduce the committed
March numbers of simulations/results/mediator_bridge_scale.txt exactly?

Targets, read off the committed file:
  Test 1 (N=11, dt_meas=0.5 grid, gamma=0.05 uniform), MI(BridgeA:BridgeB):
      t=2  0.400630   t=4  0.685503   t=6  0.576939   t=8  0.524207
      t=10 0.510783   t=12 0.468091   t=14 0.456110   t=20 0.457382
      peak over the 41-point grid: 0.777324
  Test 2(b) (N=11, coarse 8-point grid {0,2,4,6,8,10,15,20}), gamma_M sweep
  on qubit 5, MI(BridgeA:BridgeB):
      0.000 -> 0.699531   0.010 -> 0.695988   0.050 -> 0.685503
      0.100 -> 0.677099   0.200 -> 0.666886   0.500 -> 0.646889
"""

import sys, time
sys.path.insert(0, __file__.rsplit("\\", 1)[0])

import numpy as np
from bridge_sector import (sector_basis, SectorPropagator, bell_on_vacuum,
                            mutual_information, mediator_bridge)

N = 11
GAMMA = 0.05
BRIDGE_A = [0, 1, 2, 3, 4]
BRIDGE_B = [6, 7, 8, 9, 10]
PAIR_A = [0, 1]
PAIR_D = [9, 10]

states, index = sector_basis(N)
print(f"sector dimension at N={N}: {len(states)}  (dense would be {2**N})")

bonds = mediator_bridge(3)
print("bonds:", bonds)

rho0 = bell_on_vacuum(N, 0, 1, states, index)

# ---------------------------------------------------------------- Test 1
t0 = time.time()
gammas = [GAMMA] * N
prop = SectorPropagator(N, bonds, gammas, states, index)
t_meas = [i * 0.5 for i in range(41)]
rows = {}
peak = [0.0, 0.0]


def cb(t, rho):
    mi = mutual_information(rho, N, states, BRIDGE_A, BRIDGE_B)
    rows[round(t, 3)] = mi
    if mi > peak[0]:
        peak[0], peak[1] = mi, t


prop.propagate(rho0, 20.0, 0.05, t_meas, cb)
print(f"\nTest 1 reproduction ({time.time() - t0:.1f}s)")
targets = {2.0: 0.400630, 4.0: 0.685503, 6.0: 0.576939, 8.0: 0.524207,
           10.0: 0.510783, 12.0: 0.468091, 14.0: 0.456110, 16.0: 0.456791,
           18.0: 0.457758, 20.0: 0.457382}
worst = 0.0
for t, want in targets.items():
    got = rows[round(t, 3)]
    d = abs(got - want)
    worst = max(worst, d)
    print(f"  t={t:5.1f}  got {got:.6f}  march {want:.6f}  diff {d:.2e}")
print(f"  peak: got {peak[0]:.6f} at t={peak[1]:.2f}   march 0.777324")
print(f"  worst curve diff: {worst:.2e}")

# ---------------------------------------------------------------- Test 2b
print("\nTest 2(b) reproduction, coarse 8-point grid")
coarse = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]
march = {0.0: 0.699531, 0.01: 0.695988, 0.05: 0.685503,
         0.1: 0.677099, 0.2: 0.666886, 0.5: 0.646889}
worst2 = 0.0
for gm in [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]:
    g = [GAMMA] * N
    g[5] = gm
    p = SectorPropagator(N, bonds, g, states, index)
    best = [0.0]

    def cb2(t, rho, best=best):
        mi = mutual_information(rho, N, states, BRIDGE_A, BRIDGE_B)
        if mi > best[0]:
            best[0] = mi

    p.propagate(rho0, 20.0, 0.05, coarse, cb2)
    d = abs(best[0] - march[gm])
    worst2 = max(worst2, d)
    print(f"  gM={gm:5.3f}  got {best[0]:.6f}  march {march[gm]:.6f}  diff {d:.2e}")
print(f"  worst sweep diff: {worst2:.2e}")

print("\nGATE:", "PASS" if max(worst, worst2) < 1e-5 else "FAIL",
      f"(worst overall {max(worst, worst2):.2e})")

# ---------------------------------------------------------------- Test 2a
print("\nTest 2(a) reproduction, J_meta sweep on the coarse 8-point grid")
march_j = {0.25: 0.438644, 0.5: 0.634487, 0.75: 0.739250, 1.0: 0.685503,
           1.5: 0.800290, 2.0: 0.847763, 3.0: 0.891614}
worst3 = 0.0
for jm in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    p = SectorPropagator(N, mediator_bridge(3, j_meta=jm), [GAMMA] * N,
                         states, index)
    best = [0.0]

    def cb3(t, rho, best=best):
        mi = mutual_information(rho, N, states, BRIDGE_A, BRIDGE_B)
        if mi > best[0]:
            best[0] = mi

    p.propagate(rho0, 20.0, 0.05, coarse, cb3)
    d = abs(best[0] - march_j[jm])
    worst3 = max(worst3, d)
    print(f"  J_meta={jm:5.2f}  got {best[0]:.6f}  march {march_j[jm]:.6f}"
          f"  diff {d:.2e}")
print(f"  worst J_meta diff: {worst3:.2e}")

print("\nFULL GATE:", "PASS" if max(worst, worst2, worst3) < 1e-5 else "FAIL",
      f"(worst overall {max(worst, worst2, worst3):.2e})")


# ---------------------------------------------------------------- Test 1, MI(A:D)
# The page states this column reproduces too; pin it rather than assert it.
print()
print("Test 1 reproduction, the MI(PairA : PairD) column")
march_ad = {2.0: 0.032949, 4.0: 0.071576, 6.0: 0.040135, 8.0: 0.026251,
            10.0: 0.035291, 12.0: 0.014361, 14.0: 0.009307, 16.0: 0.009462,
            18.0: 0.009035, 20.0: 0.008562}
p_ad = SectorPropagator(N, bonds, [GAMMA] * N, states, index)
rows_ad = {}
p_ad.propagate(rho0, 20.0, 0.05, t_meas,
               lambda t, rho: rows_ad.__setitem__(
                   round(t, 3), mutual_information(rho, N, states, PAIR_A, PAIR_D)))
worst4 = 0.0
for t, want in march_ad.items():
    got = rows_ad[round(t, 3)]
    worst4 = max(worst4, abs(got - want))
    print(f"  t={t:5.1f}  got {got:.6f}  march {want:.6f}  diff {abs(got - want):.2e}")
print(f"  worst MI(A:D) diff: {worst4:.2e}")
print()
print("GATE INCLUDING THE SECOND COLUMN:",
      "PASS" if max(worst, worst2, worst3, worst4) < 1e-5 else "FAIL",
      f"(worst overall {max(worst, worst2, worst3, worst4):.2e})")
