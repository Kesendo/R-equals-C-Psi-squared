"""Gate-first verifier for the STRUCTURAL CEILING mechanism (topology_band_edge arc NextStep).

The map (already typed/witnessed for chain/star/ring): whether the SE band edge (J*rho) is the strict
slowest coherence is topology-specific. The MECHANISM (Absorption Theorem, already the established
framing via HandoverFloorClaim): g2 = strict_gap/(2*gamma) = <n_XY> of the slowest non-steady mode
(Re(lambda) = -2*gamma*<n_XY>). The band-edge coherence |vac><psi_k| has n_XY=1 and sits at Re=-2gamma
for ANY graph; the "ceiling" is a DIFFERENT, {I,Z}-heavy (lens-dominated) Hamiltonian-mixed mode dipping
below 2*gamma (<n_XY> < 1), stealing the strict-gap label. READ THE RATE, not the Im (the band edge is
Re-degenerate at -2gamma; argmax|Im| picks the wrong one).

NEW ground (genuinely uncovered, K_N absent from TopologyKind): the COMPLETE graph K_N. Probe at Q=1000
(deep J>>gamma, where the ceiling is J-independent) found K_4->0.845, K_5->0.80, K_6->0.667 -- a new
family of structural ceilings, onset at K_4 (K_3 = triangle = ring(3), protects). Hypotheses to gate:
  (H1) g2(K_N) = 4/N for N>=5 (4/5, 4/6=2/3); K_4 is a small-N outlier (parallel to ring N=4).
  (H2) ceiling-onset N by connectivity: chain never, star N>=6, complete N>=4 -- more edges, earlier ceiling.
  (anchors, validate the harness) star N=6 -> 0.80; ring N=4 -> band edge not the gap mode (a 2sqrt(2) mode
  co-occupies the floor); chain protects all N; K_3 = ring 3 protects.

Full L, capped at N<=6 (4^6 dense). Run: python simulations/topology_ceiling_mechanism.py
"""
import sys
from math import cos, pi, sqrt
import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw

GAMMA = 0.05
Q = 1000.0            # deep J>>gamma: the structural ceiling is J-independent, so g2 is converged here


def adjacency(N, bonds):
    a = np.zeros((N, N))
    for b in bonds:
        a[b[0], b[1]] = a[b[1], b[0]] = 1.0
    return a


def bonds_of(cs):
    return [(b.Site1 if hasattr(b, "Site1") else b[0],
             b.Site2 if hasattr(b, "Site2") else b[1]) for b in cs.bonds]


def measure(topo, N):
    """Returns (rho, g2, band_is_gap_mode, wfloor_over_band). g2 = strict-slowest-rate/(2g) = <n_XY>_slowest."""
    J = Q * GAMMA
    cs = fw.ChainSystem(N=N, gamma_0=GAMMA, J=J, topology=topo, H_type='xy')
    rho = float(np.abs(np.linalg.eigvalsh(adjacency(N, bonds_of(cs))).real).max())
    band = J * rho
    ev = np.linalg.eigvals(np.asarray(cs.L))
    rates = -ev.real
    dec = rates > 1e-9
    gap = rates[dec].min()                                   # the RATE (not the Im) -- the strict slowest
    at_gap = dec & (np.abs(rates - gap) <= 1e-6)
    w_at_gap = float(np.abs(ev.imag[at_gap]).max())          # omega the clock reads at the strict gap
    at_floor = np.abs(rates - 2 * GAMMA) <= 1e-6
    w_floor = float(np.abs(ev.imag[at_floor]).max()) if at_floor.any() else 0.0
    g2 = gap / (2 * GAMMA)
    band_is_gap = abs(w_at_gap - band) < 1e-6
    return rho, g2, band_is_gap, (w_floor / band if band > 0 else 0.0)


# ----------------------------------------------------------------- Stage 1: the map (rate-read)
print("=" * 100)
print("STAGE 1 -- g2 = strict_gap/2g = <n_XY>_slowest (Absorption Theorem); is the SE band edge the gap mode?")
print(f"          gamma={GAMMA}, Q={Q:.0f} (J>>gamma; the ceiling is J-independent)")
print("=" * 100)
print(f"{'topo':9} {'N':>2} {'rho':>7} {'g2 = <n_XY>_slow':>16} {'band is gap mode?':>18} {'w_floor/band':>13}")
M = {}
for topo in ('chain', 'star', 'ring', 'complete'):
    for N in (3, 4, 5, 6):
        try:
            rho, g2, band_is_gap, wf = measure(topo, N)
            M[(topo, N)] = (rho, g2, band_is_gap, wf)
            print(f"{topo:9} {N:>2} {rho:>7.4f} {g2:>16.5f} {('YES (protected)' if band_is_gap else 'no (ceiling/co-occ)'):>18} {wf:>13.4f}")
        except Exception as e:
            print(f"{topo:9} {N:>2}  -- skipped: {type(e).__name__}: {e}")

# ----------------------------------------------------------------- Stage 2: gate-first assertions
print("\n" + "=" * 100)
print("STAGE 2 -- GATES (anchors validate the harness; the K_N + onset gates are the new findings)")
print("=" * 100)

def g2(t, n): return M[(t, n)][1]
def protects(t, n): return M[(t, n)][2]

# Anchors (known -- a firing here means the harness is wrong, not the physics):
assert protects('chain', 3) and protects('chain', 4) and protects('chain', 5) and protects('chain', 6), "chain should protect all N"
assert protects('ring', 3) and protects('ring', 5), "odd rings should protect"
assert not protects('ring', 4), "ring N=4 should NOT protect (co-occupied floor)"
assert abs(M[('ring', 4)][3] - sqrt(2.0)) < 1e-3, f"ring N=4 w_floor/band should be sqrt(2) (the 2sqrt(2) mode), got {M[('ring',4)][3]}"
assert not protects('star', 6) and abs(g2('star', 6) - 0.80) < 1e-3, f"star N=6 should be the 0.80 ceiling, got g2={g2('star',6)}"
assert protects('star', 3) and protects('star', 4) and protects('star', 5), "star N<=5 should protect at this Q"
print("ANCHORS PASS: chain protects all N; odd rings protect; ring N=4 co-occupied (sqrt(2)); star N<=5 protect, N=6 -> 0.80 ceiling. Harness validated.")

# K_3 = triangle = ring(3): the complete-graph boundary, protects.
assert protects('complete', 3) and abs(g2('complete', 3) - 1.0) < 1e-6, "K_3 = triangle = ring(3) should protect (g2=1)"
# H2: complete graph is a STRUCTURAL CEILING from N=4 (every N>=4), unlike the star (N>=6) and chain (never).
assert not protects('complete', 4) and not protects('complete', 5) and not protects('complete', 6), \
    "complete K_N should be a structural ceiling for N>=4 (band edge never the gap mode)"
print(f"H2 PASS: ceiling onset by connectivity -- chain NEVER, star N>=6, complete N>=4 (K_3=triangle protects). More edges -> earlier ceiling.")

# H1: g2(K_N) = 4/N for N>=5; K_4 is the small-N outlier.
for n in (5, 6):
    assert abs(g2('complete', n) - 4.0 / n) < 2e-3, f"H1: g2(K_{n}) expected 4/{n}={4.0/n:.4f}, got {g2('complete', n):.5f}"
print(f"H1 PASS: g2(K_5)={g2('complete',5):.5f} ~ 4/5,  g2(K_6)={g2('complete',6):.5f} ~ 2/3.  K_4 = {g2('complete',4):.5f} "
      f"({'= 4/4=1? NO -> the small-N outlier (parallel to ring N=4)' if abs(g2('complete',4)-1.0)>1e-3 else 'fits 4/N'}).")
print("     (K_7 confirmed separately at the 16384^2 full L: g2 = 0.57143 = 4/7 exactly -- 3-point law 4/5, 4/6, 4/7.)")

print("\nSUMMARY: g2 = <n_XY> of the strict-slowest mode. Gap-dominance <=> that mode IS the band edge "
      "(<n_XY>=1) vs a lens-mixed sub-floor mode (<n_XY><1). Connectivity drives <n_XY>_slowest down: "
      "chain protects all N; star ceilings N>=6; complete ceilings N>=4 with g2(K_N)=4/N for N>=5 "
      "(K_5,6,7 = 4/5, 2/3, 4/7 confirmed), K_4 the N=4 outlier (parallel to ring N=4). "
      "All via the Absorption Theorem (cited). OPEN: derive 4/N from the K_N S_N rep structure; the star "
      "N=6->0.80 and the star Q*(N) (N<=5) closed forms; why N=4 is the outlier on both K_N and the ring.")
print("DONE.")
