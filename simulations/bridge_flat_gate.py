"""NextStep (1) of the arc `the_gate_that_does_not_gate`: read the flat gate.

The question: at Level 3 (N=11) the March sweep of gamma_M (the meta-mediator's
own dephasing) moved peak MI(BridgeA:BridgeB) by only 7.5% over 0 -> 0.5, while
the SAME parameter at Level 2 (N=5) is reported to move an MI from 0.44 to 0.12.
Is the Level-3 flatness real insensitivity, or an artefact of the coarse
8-point measurement grid (whose own grid error, 13%, is LARGER than the effect)?

Engine: simulations/bridge_sector.py, calibrated against every committed March
number to 4.4e-7 by bridge_gate_calibrate.py.

Four readings:
  A. grid convergence of the Level-3 gamma_M sweep (coarse / 0.5 / 0.1 / 0.05)
  B. the same sweep pushed far past 0.5, to see whether the gate ever bites
  C. the background control: is the flatness about the LEVEL, or about the
     other ten sites already carrying gamma=0.05 while gamma_M is varied?
  D. the Level-2 (N=5) sweep on the same fine grid and the same observable
     shape (the two flanking pairs), in both backgrounds
"""

import sys, time
sys.path.insert(0, __file__.rsplit("\\", 1)[0])

import numpy as np
from bridge_sector import (sector_basis, SectorPropagator, bell_on_vacuum,
                            mutual_information, mediator_bridge)

TMAX = 20.0
DT = 0.05


def peak_mi(n, bonds, gammas, keep_a, keep_b, t_meas, states, index, rho0):
    p = SectorPropagator(n, bonds, gammas, states, index)
    best = [0.0, 0.0]

    def cb(t, rho):
        mi = mutual_information(rho, n, states, keep_a, keep_b)
        if mi > best[0]:
            best[0], best[1] = mi, t

    p.propagate(rho0, TMAX, DT, t_meas, cb)
    return best[0], best[1]


def grid(dt_meas):
    k = int(round(TMAX / dt_meas))
    return [i * dt_meas for i in range(k + 1)]


COARSE = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]

# ------------------------------------------------------------------ setup
N3 = 11
states3, index3 = sector_basis(N3)
rho3 = bell_on_vacuum(N3, 0, 1, states3, index3)
bonds3 = mediator_bridge(3)
A3, B3 = [0, 1, 2, 3, 4], [6, 7, 8, 9, 10]

N2 = 5
states2, index2 = sector_basis(N2)
rho2 = bell_on_vacuum(N2, 0, 1, states2, index2)
bonds2 = mediator_bridge(2)
A2, B2 = [0, 1], [3, 4]

GM = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]


def sweep(n, bonds, med, a, b, states, index, rho0, gms, t_meas, bg):
    out = []
    for gm in gms:
        g = [bg] * n
        g[med] = gm
        mi, tp = peak_mi(n, bonds, g, a, b, t_meas, states, index, rho0)
        out.append((gm, mi, tp))
    return out


def sweep_fine(n, bonds, med, a, b, states, index, rho0, gms, bg):
    """Every-RK4-step measurement, with dt shrunk where a large gamma_M would
    otherwise make RK4 unstable.  Returns (gamma_M, peak MI, t*, dt used)."""
    out = []
    for gm in gms:
        g = [bg] * n
        g[med] = gm
        p = SectorPropagator(n, bonds, g, states, index)
        dt = p.stable_dt(DT)
        best = [0.0, 0.0]

        def cb(t, rho, best=best):
            mi = mutual_information(rho, n, states, a, b)
            if mi > best[0]:
                best[0], best[1] = mi, t

        p.propagate_every_step(rho0, TMAX, dt, cb)
        out.append((gm, best[0], best[1], dt))
    return out


def show_fine(title, rows):
    print(f"\n{title}")
    print(f"  {'gamma_M':>8}  {'peak MI':>10}  {'t*':>6}  {'dt':>8}")
    print("  " + "-" * 38)
    for gm, mi, tp, dt in rows:
        print(f"  {gm:8.3f}  {mi:10.6f}  {tp:6.2f}  {dt:8.5f}")
    lo, hi = min(r[1] for r in rows), max(r[1] for r in rows)
    print(f"  span (max-min)/max = {(hi - lo) / hi * 100:.2f}%"
          f"   [{lo:.6f} .. {hi:.6f}]")


def show(title, rows):
    print(f"\n{title}")
    print(f"  {'gamma_M':>8}  {'peak MI':>10}  {'t*':>6}")
    print("  " + "-" * 28)
    for gm, mi, tp in rows:
        print(f"  {gm:8.3f}  {mi:10.6f}  {tp:6.2f}")
    lo, hi = min(r[1] for r in rows), max(r[1] for r in rows)
    print(f"  span (max-min)/max = {(hi - lo) / hi * 100:.2f}%"
          f"   [{lo:.6f} .. {hi:.6f}]")
    return (hi - lo) / hi * 100


# ------------------------------------------------------- A. grid convergence
print("=" * 70)
print("A. LEVEL 3 (N=11), gamma_M sweep, background gamma=0.05, four t-grids")
print("=" * 70)
t0 = time.time()
spans = {}
for label, tm in [("March coarse 8-point", COARSE), ("dt_meas=0.5", grid(0.5)),
                  ("dt_meas=0.1", grid(0.1)), ("dt_meas=0.05 (every RK4 step)",
                                               grid(0.05))]:
    rows = sweep(N3, bonds3, 5, A3, B3, states3, index3, rho3, GM, tm, 0.05)
    spans[label] = show(label, rows)
print(f"\n  spans by grid: " +
      ", ".join(f"{k}: {v:.2f}%" for k, v in spans.items()))
print(f"  [{time.time() - t0:.1f}s]")

# ------------------------------------------------------- B. push the gate
print("\n" + "=" * 70)
print("B. LEVEL 3, gamma_M pushed far past 0.5 (every-step measurement)")
print("=" * 70)
wide = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0, 1000.0]
show_fine("gamma_M far sweep, background 0.05",
          sweep_fine(N3, bonds3, 5, A3, B3, states3, index3, rho3, wide, 0.05))

# ------------------------------------------------------- C. background control
print("\n" + "=" * 70)
print("C. LEVEL 3, SAME sweep with the other ten sites QUIET (background 0)")
print("=" * 70)
show_fine("gamma_M sweep, background 0",
          sweep_fine(N3, bonds3, 5, A3, B3, states3, index3, rho3, GM, 0.0))
show_fine("gamma_M far sweep, background 0",
          sweep_fine(N3, bonds3, 5, A3, B3, states3, index3, rho3, wide, 0.0))

# ------------------------------------------------------- D. level 2
print("\n" + "=" * 70)
print("D. LEVEL 2 (N=5), MI({0,1}:{3,4}), mediator = site 2, every-step grid")
print("=" * 70)
show_fine("gamma_M sweep, background 0.05",
          sweep_fine(N2, bonds2, 2, A2, B2, states2, index2, rho2, GM, 0.05))
show_fine("gamma_M far sweep, background 0.05",
          sweep_fine(N2, bonds2, 2, A2, B2, states2, index2, rho2, wide, 0.05))
show_fine("gamma_M sweep, background 0",
          sweep_fine(N2, bonds2, 2, A2, B2, states2, index2, rho2, GM, 0.0))
show_fine("gamma_M far sweep, background 0",
          sweep_fine(N2, bonds2, 2, A2, B2, states2, index2, rho2, wide, 0.0))

# ------------------------------------------------------- E. the other arm
# The J_meta sweep it is compared against was read on the SAME coarse grid,
# so it needs the same correction before the two can be set side by side.
print(SEP if False else "\n" + "=" * 70)
print("E. LEVEL 3, J_meta sweep: coarse grid vs converged grid")
print("=" * 70)
print(f"  {'J_meta':>7}  {'coarse max':>11}  {'fine peak':>10}  {'t*':>6}")
print("  " + "-" * 40)
for jm in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]:
    bonds = mediator_bridge(3, j_meta=jm)
    pc = SectorPropagator(11, bonds, [0.05] * 11, states3, index3)
    bc = [0.0]
    pc.propagate(rho3, TMAX, DT, COARSE,
                 lambda t, r, b=bc: b.__setitem__(
                     0, max(b[0], mutual_information(r, 11, states3, A3, B3))))
    pf = SectorPropagator(11, bonds, [0.05] * 11, states3, index3)
    bf = [0.0, 0.0]

    def cbe(t, r, b=bf):
        mi = mutual_information(r, 11, states3, A3, B3)
        if mi > b[0]:
            b[0], b[1] = mi, t

    pf.propagate_every_step(rho3, TMAX, DT, cbe)
    print(f"  {jm:7.2f}  {bc[0]:11.6f}  {bf[0]:10.6f}  {bf[1]:6.2f}")
