"""Level 2 vs Level 3: the gamma_M contrast, decomposed one factor at a time.

The contrast under test:
  Level 2 (N=5),  hypotheses/THE_OTHER_SIDE.md section 22 Test 3:  0.44 -> 0.12
  Level 3 (N=11), simulations/results/mediator_bridge_scale.txt:   0.6995 -> 0.6469

They are not the same measurement.  Three things differ besides the level:
  (i)   FUNCTIONAL: Level 2 reads MI at the fixed time t = 5.0
        (simulations/mediator_bridge.py:565 t_measure, :621 the sweep);
        Level 3 takes a max over the coarse grid {0,2,4,6,8,10,15,20}
        (Propagate/Program.cs:879).
  (ii)  BLOCKS: Level 2 uses the 2-qubit end pairs {0,1}:{3,4}; Level 3 uses
        the 5-qubit halves {0..4}:{6..10}.
  (iii) INITIAL STATE: Level 2 uses 'bell_A_0M_pp_B' = Bell(0,1) x |0> x |+>|+>
        (mediator_bridge.py:223-226); Level 3 uses Bell(0,1) x |0>^9.

This script changes one at a time.  N=5 runs dense (32x32) through the committed
primitives of simulations/mediator_bridge.py, so its conventions are the
repo's own; N=11 runs through the sector engine bridge_sector.py, calibrated
to 4.4e-7 against every committed March number.
"""

import sys
import numpy as np
from scipy.linalg import expm

HERE = __file__.rsplit("\\", 1)[0]
sys.path.insert(0, HERE)

import mediator_bridge as mb
from bridge_sector import (sector_basis, SectorPropagator, bell_on_vacuum,
                            mutual_information, mediator_bridge)

SEP = chr(10) + '=' * 74
NL = chr(10)
GM = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
GAMMA = 0.05
TMAX = 20.0
DT = 0.05


# ----------------------------------------------------------------- N=5 dense
def n5_state(name):
    return mb.make_initial_state(name)


def n5_bell_vacuum():
    up = np.array([1, 0], dtype=complex)
    psi_a = (np.kron(up, up) + np.kron([0, 1], [0, 1])) / np.sqrt(2)
    psi = np.kron(np.kron(psi_a, up), np.kron(up, up))
    return np.outer(psi, psi.conj())


def n5_run(gm, rho0, functional, keep_a=(0, 1), keep_b=(3, 4)):
    gammas = [GAMMA, GAMMA, gm, GAMMA, GAMMA]
    L, _, _ = mb.build_mediator_system(gammas=gammas)
    if functional == "fixed5":
        rho = mb.evolve_rho(L, rho0, 5.0)
        return mb.mutual_info(rho, 5, list(keep_a), list(keep_b)), 5.0
    # fine-grid peak: step the propagator on the same 0.05 grid
    best, tbest = 0.0, 0.0
    step = expm(L * DT)
    vec = rho0.flatten()
    for k in range(int(TMAX / DT) + 1):
        rho = vec.reshape(32, 32)
        rho = (rho + rho.conj().T) / 2
        mi = mb.mutual_info(rho, 5, list(keep_a), list(keep_b))
        if mi > best:
            best, tbest = mi, k * DT
        vec = step @ vec
    return best, tbest


# ---------------------------------------------------------------- N=11 sector
S3, I3 = sector_basis(11)
RHO3 = bell_on_vacuum(11, 0, 1, S3, I3)
BONDS3 = mediator_bridge(3)


def n11_run(gm, functional, keep_a, keep_b):
    g = [GAMMA] * 11
    g[5] = gm
    p = SectorPropagator(11, BONDS3, g, S3, I3)
    if functional == "fixed5":
        best = [0.0, 5.0]

        def cb(t, rho):
            if abs(t - 5.0) < 1e-9:
                best[0] = mutual_information(rho, 11, S3, keep_a, keep_b)

        p.propagate(RHO3, TMAX, DT, [5.0], cb)
        return best[0], 5.0
    best = [0.0, 0.0]

    def cb2(t, rho):
        mi = mutual_information(rho, 11, S3, keep_a, keep_b)
        if mi > best[0]:
            best[0], best[1] = mi, t

    p.propagate_every_step(RHO3, TMAX, DT, cb2)
    return best[0], best[1]


def report(title, runner):
    rows = [(gm,) + runner(gm) for gm in GM]
    print(f"\n{title}")
    print(f"  {'gamma_M':>8}  {'MI':>10}  {'t':>6}")
    print("  " + "-" * 28)
    for gm, mi, t in rows:
        print(f"  {gm:8.3f}  {mi:10.6f}  {t:6.2f}")
    lo, hi = rows[-1][1], rows[0][1]
    print(f"  drop over 0 -> 0.5: {(hi - lo) / hi * 100:.2f}%"
          f"   [{hi:.6f} -> {lo:.6f}]")
    return (hi - lo) / hi * 100


print("=" * 74)
print("LEVEL 2 (N=5), dense, through the committed mediator_bridge.py primitives")
print("=" * 74)
rho_pp = n5_state('bell_A_0M_pp_B')
rho_vac = n5_bell_vacuum()

s = {}
s["L2 march (fixed t=5, pairs, |++> on B)"] = report(
    "1. AS PUBLISHED: fixed t=5.0, pairs {0,1}:{3,4}, Bell x |0> x |++>",
    lambda gm: n5_run(gm, rho_pp, "fixed5"))
s["L2 peak instead of fixed t"] = report(
    "2. ONE CHANGE: fine-grid peak instead of the fixed t=5.0 snapshot",
    lambda gm: n5_run(gm, rho_pp, "peak"))
s["L2 vacuum instead of |++>"] = report(
    "3. ONE CHANGE: Bell x |0>^3 (vacuum) instead of |++> on B, fixed t=5.0",
    lambda gm: n5_run(gm, rho_vac, "fixed5"))
s["L2 both changes"] = report(
    "4. BOTH: vacuum initial state AND fine-grid peak",
    lambda gm: n5_run(gm, rho_vac, "peak"))

print(SEP)
print("LEVEL 3 (N=11), sector engine, Bell x |0>^9 throughout")
print("=" * 74)

COARSE = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]


def n11_coarse(gm, keep_a, keep_b):
    g = [GAMMA] * 11
    g[5] = gm
    p = SectorPropagator(11, BONDS3, g, S3, I3)
    best = [0.0, 0.0]

    def cb(t, rho):
        mi = mutual_information(rho, 11, S3, keep_a, keep_b)
        if mi > best[0]:
            best[0], best[1] = mi, t

    p.propagate(RHO3, TMAX, DT, COARSE, cb)
    return best[0], best[1]


HALVES_A, HALVES_B = [0, 1, 2, 3, 4], [6, 7, 8, 9, 10]
ENDS_A, ENDS_B = [0, 1], [9, 10]

s["L3 march (coarse max, 5-blocks)"] = report(
    "5. AS PUBLISHED: max over the coarse 8-point grid, 5-qubit halves",
    lambda gm: n11_coarse(gm, HALVES_A, HALVES_B))
s["L3 fine peak, 5-blocks"] = report(
    "6. ONE CHANGE: fine-grid peak, still the 5-qubit halves",
    lambda gm: n11_run(gm, "peak", HALVES_A, HALVES_B))
s["L3 fixed t=5, 5-blocks"] = report(
    "7. ONE CHANGE: fixed t=5.0 snapshot, still the 5-qubit halves",
    lambda gm: n11_run(gm, "fixed5", HALVES_A, HALVES_B))
s["L3 fixed t=5, end pairs"] = report(
    "8. TWO CHANGES: fixed t=5.0 AND the 2-qubit end pairs {0,1}:{9,10}"
    "  (= the Level-2 functional, carried to Level 3)",
    lambda gm: n11_run(gm, "fixed5", ENDS_A, ENDS_B))
s["L3 fine peak, end pairs"] = report(
    "9. fine-grid peak on the 2-qubit end pairs",
    lambda gm: n11_run(gm, "peak", ENDS_A, ENDS_B))

print(SEP)
print("SUMMARY: drop in MI over gamma_M = 0 -> 0.5")
print("=" * 74)
for k, v in s.items():
    print(f"  {k:<44} {v:6.2f}%")
