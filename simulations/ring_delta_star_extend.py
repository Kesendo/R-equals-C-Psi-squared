"""Gate-first: extend the ring XXZ Delta*(N) descent PAST N=14 -- does the non-monotone hump return
toward Delta=1 (like the chain) or settle above 1? (the xxz_axis_handover arc's open ring question)

RESULT (2026-06-17): Delta*(15) = 1.27413 (odd parity), via the CERTIFIED gamma->0 reduction (reused from
ring_xxz_delta_star_descent.py: Delta* = the Neel-side DOWNWARD crossing of gap(R)=2, half-filling sector,
factor-4 Z<->n; never the 4^N Liouvillian). The odd-parity descent 9->15:
    N:   9      11      13      15
    D*:  1.33095 1.31596 1.29499 1.27413
    step:    -0.0150  -0.0210  -0.0209
The steps do NOT shrink (they grow then hold ~0.021) -- a descent SETTLING to a plateau above 1 would
DECELERATE (shrinking steps toward the asymptote). The steady, non-decelerating descent further refutes the
Round-2 "plateau trap" (N<=11 looked settled at ~1.31) one N deeper, and DISFAVORS "settles above 1"; it is
consistent with a slow descent toward Delta=1 (the chain's limit), though the N->inf limit stays FORMALLY
OPEN: at N=15 Delta* is still 1.274 (far from 1), and resolving 1-vs-settles needs N>>16 (infeasible -- dense
eigh past C(N,p)~25000). A power-law fit degenerates even on the descending tail (no fit-based limit).

How Delta*(N>14) is obtained: a coarse scan of [1.15,1.45] for the Neel-side downward crossing of gap(R)=2,
then the committed _bisect_cross -- a handful of gapR calls (each = eigh of the C(N,p) sector H + eigvalsh of
R). N=15 (C(15,8)=6435) ~16 min; N=16 (C(16,8)=12870) ~1 h (even parity, not run -- the limit is unreachable
at feasible N either way). REGIME (load-bearing): gamma->0 reduction NOT finite-gamma Q=20; half-filling
p=(N+1)//2; Neel-side downward crossing; non-monotone-aware; N=4 excluded (no handover).

  Default run (FAST): port-check dstar_fast vs the committed RING_DSTAR_SEQUENCE at small N + assert the baked
    N=15 continues the odd descent with non-shrinking steps. No slow recompute.
  Regeneration (SLOW): `python ring_delta_star_extend.py 15 [16 ...]` recomputes Delta*(N) for those N.

Run: python simulations/ring_delta_star_extend.py
"""
import sys
import numpy as np

sys.path.insert(0, 'simulations')
from ring_xxz_delta_star_descent import gapR, _bisect_cross, RING_DSTAR_SEQUENCE  # certified machinery

# Delta*(N) for N>14 (gamma->0 reduction; slow -- regenerate with the CLI args below). Odd parity, N=15.
RING_DSTAR_HIGHN = {15: 1.27413}


def dstar_fast(N, lo=1.15, hi=1.45, step=0.05):
    """Delta* via a coarse bracket of the Neel-side downward crossing of gap(R)=2, then _bisect_cross."""
    grid = np.arange(lo, hi + 1e-9, step)
    g = [gapR(N, float(D), ring=True) for D in grid]
    for i in range(len(grid) - 1):
        if g[i] - 2 > 0 >= g[i + 1] - 2:
            return _bisect_cross(N, ring=True, lo=float(grid[i]), hi=float(grid[i + 1]))
    return None


# ---- SLOW regeneration path (opt-in via CLI args) ----
if len(sys.argv) > 1:
    for N in (int(a) for a in sys.argv[1:]):
        ds = dstar_fast(N)
        committed = RING_DSTAR_SEQUENCE.get(N) or RING_DSTAR_HIGHN.get(N)
        tag = f"(committed {committed:.5f}, |d|={abs(ds-committed):.1e})" if committed else "(NEW)"
        print(f"  ring Delta*({N}) = {ds:.5f}  {tag}")
    sys.exit(0)

# ---- FAST default: port-check + the descent-continues gate ----
print("=" * 92)
print("RING Delta*(N) extension -- the descent continues past N=14 (Delta*(15) = 1.27413)")
print("=" * 92)

# STAGE 0 -- port fidelity at small (fast) N: dstar_fast reproduces the committed grid-scan values
print("STAGE 0 -- port fidelity (coarse-bracket+bisect == committed delta_star_ring):")
worst = 0.0
for N in (8, 10):
    got = dstar_fast(N)
    ref = RING_DSTAR_SEQUENCE[N]
    worst = max(worst, abs(got - ref))
    print(f"  N={N}: dstar_fast={got:.5f}  committed={ref:.5f}  |d|={abs(got-ref):.1e}")
    assert abs(got - ref) < 3e-3, f"STAGE 0 GATE FIRED: dstar_fast({N})={got} != committed {ref}"
print(f"STAGE 0 PASS: the targeted bracketing reproduces the committed method (worst |d|={worst:.1e}).")

# STAGE 1 -- the descent continues, and the steps do NOT shrink (no plateau deceleration)
print("\nSTAGE 1 -- the odd-parity descent 9->15 and its steps:")
odd = {9: RING_DSTAR_SEQUENCE[9], 11: RING_DSTAR_SEQUENCE[11], 13: RING_DSTAR_SEQUENCE[13],
       15: RING_DSTAR_HIGHN[15]}
ns = sorted(odd)
steps = [odd[ns[i + 1]] - odd[ns[i]] for i in range(len(ns) - 1)]
for i, n in enumerate(ns):
    s = f"   step {steps[i-1]:+.4f}" if i else ""
    print(f"  N={n}: Delta* = {odd[n]:.5f}{s}")
assert odd[15] < odd[13] < odd[11] < odd[9], "STAGE 1 GATE FIRED: odd parity not monotone-descending 9->15"
# a plateau-above-1 would DECELERATE (each |step| smaller); here the later steps are NOT smaller than the first
assert abs(steps[-1]) >= abs(steps[0]) - 1e-9, \
    f"STAGE 1 GATE FIRED: descent is DECELERATING (steps {steps}) -- would indicate settling to a plateau"
print(f"STAGE 1 PASS: descent continues to N=15 (Delta*=1.27413); the latest step is NOT smaller than the")
print(f"  first ({[f'{s:+.4f}' for s in steps]}) -> NO overall deceleration toward a plateau; 'settles above 1' disfavored.")

# STAGE 2 -- the honest limit read
print("\nSTAGE 2 -- honest limit read:")
print("  At N=15 Delta* = 1.274 is still far above 1; the steady (non-decelerating) descent is consistent")
print("  with a slow approach to Delta=1 (the chain's BKT limit) and DISFAVORS settling above 1, but the")
print("  N->inf limit stays FORMALLY OPEN -- resolving it needs N>>16 (infeasible: dense eigh past")
print("  C(N,p)~25000), and a power-law fit degenerates even on the descending tail. New evidence shifts")
print("  the weight (no plateau), it does not close the question.")
print("\nDONE.")
