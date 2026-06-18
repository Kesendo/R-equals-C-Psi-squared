#!/usr/bin/env python3
"""Gate-first: the STAR survivor NEVER un-freezes -- its own seam, the third member of the trichotomy.

A single-Q snapshot does not distinguish the topologies: below its coherence horizon EVERY topology's
slowest mode is overdamped (real, |Im|=0). incompleteness_survivor.py at Q=1.5 finds chain (2,2), ring
(2,2), star (1,1) ALL frozen. The distinction is the |Im|(Q) CURVE -- does the survivor UN-FREEZE
(acquire a frequency) as the dephasing weakens (Q grows)?

  chain : un-freezes at Q*(N) (the coherence horizon; the {0,2}-EP pair oscillates above it).
  ring  : un-freezes at its handover (the oscillating band edge overtakes the frozen (2,2) seam).
  star  : frozen at EVERY Q -- iff g2 = 4/(N-1) <= 1, i.e. N >= 5 (N=4 is the known outlier).

So the star has its OWN frozen seam, distinct from the ring's (2,2) level-crossing one: for N >= 5 the
star's (1,1) BOUNDARY survivor is frozen at every Q. And the THRESHOLD is exactly the structural ceiling:
the survivor is the darkest [H,A]=0 commutant (1,1) coherence only when it is darker than the -2g floor,
i.e. g2 = 4/(N-1) <= 1 (N >= 5); at N=4 (g2 = 4/3 > 1) the commutant mode is not dark enough, an
oscillating band-edge mode wins, and the star un-freezes. So the frozen seam IS PROOF_STRUCTURAL_CEILING's
g2 = 4/(N-1) <= 1, read DYNAMICALLY (sec.7's "no horizon" at the survivor level). It is the third case of
the chain(oscillates) / ring(frozen crossing) / star(frozen iff g2<=1) trichotomy the typed layer named
only for chain and ring.

Reuses simulations/carbon/incompleteness_survivor.py: survivor(N,J,g,bnds) -> (Re,|Im|,sector,<n_XY>),
the sector-projected global slowest, validated bit-for-bit vs the full 4^N L at N=4 (all topologies).
"""
import sys
import numpy as np

sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import survivor, bonds   # noqa: E402

FROZEN = 1e-6        # |Im| below this = frozen (overdamped / non-oscillating)
QS = [0.8, 1.5, 3.0, 6.0, 12.0, 25.0, 50.0]


def im_of(N, topo, Q, J=1.0):
    _, im, sec, _ = survivor(N, J, J / Q, bonds(N, topo))
    return im, sec


def _pin_at_Q15():
    """Stage 0 (pin to source): at Q=1.5 ALL THREE survivors are frozen -- a single Q does not separate
    them (this is why the |Im|(Q) sweep is needed)."""
    print("[Stage 0 / pin] at Q=1.5 every topology's survivor is frozen (overdamped below its horizon):")
    for topo in ("chain", "ring", "star"):
        for N in (4, 6):
            im, sec = im_of(N, topo, 1.5)
            print(f"     {topo:>5} N={N}: sector={sec} |Im|={im:.3e}")
            assert im < FROZEN, f"{topo} N={N} Q=1.5 should be frozen, got |Im|={im}"
    print("     => |Im|=0 at one Q is NOT a topology signature. The signature is whether it un-freezes.")


def _sweep_and_gate(N=6):
    """Stage 1 (the finding): sweep Q. The star survivor stays frozen at EVERY Q (never un-freezes);
    the chain and ring un-freeze (acquire |Im|>0) as Q grows past their horizon/handover."""
    print(f"\n[Stage 1 / the |Im|(Q) sweep] N={N}: does the survivor un-freeze as Q grows?")
    print(f"  {'Q':>6}" + "".join(f"{t:>22}" for t in ("chain |Im| (sector)", "ring |Im| (sector)",
                                                        "star |Im| (sector)")))
    star_max = 0.0
    chain_max = ring_max = 0.0
    for Q in QS:
        cells = {}
        for topo in ("chain", "ring", "star"):
            im, sec = im_of(N, topo, Q)
            cells[topo] = f"{im:.2e} {sec}"
            if topo == "star":
                star_max = max(star_max, im)
            elif topo == "chain":
                chain_max = max(chain_max, im)
            else:
                ring_max = max(ring_max, im)
        print(f"  {Q:>6.1f}" + "".join(f"{cells[t]:>22}" for t in ("chain", "ring", "star")))

    # the finding: the star NEVER un-freezes; the chain and ring DO.
    assert star_max < FROZEN, (
        f"GATE FIRED: the star survivor un-freezes somewhere (max |Im|={star_max:.2e} > {FROZEN}). The star "
        f"would then have a horizon/oscillating survivor -- contradicts PROOF_STRUCTURAL_CEILING sec.7. Diagnose.")
    assert chain_max > 1e-2, (
        f"GATE FIRED: the chain survivor never un-freezes in the sweep (max |Im|={chain_max:.2e}); it should "
        f"oscillate above Q*(N). Diagnose (widen Q, or the sector survivor missed the band edge).")
    assert ring_max > 1e-2, (
        f"GATE FIRED: the ring survivor never un-freezes (max |Im|={ring_max:.2e}); it should oscillate above "
        f"its handover. Diagnose.")
    print(f"\n  star max |Im| over the sweep = {star_max:.2e} (frozen at EVERY Q); "
          f"chain {chain_max:.2f}, ring {ring_max:.2f} (un-freeze).")
    print("  => the STAR survivor never acquires a frequency: its own frozen seam, the third member of")
    print("  chain(un-freezes at Q*) / ring(un-freezes at the handover) / star(frozen at all Q). This is")
    print("  the survivor-level form of PROOF_STRUCTURAL_CEILING sec.7's no-horizon.")


def _cross_N():
    """Stage 2 (the threshold): the star is frozen at all Q exactly when g2 = 4/(N-1) <= 1, i.e. N >= 5.
    N=4 (g2 = 4/3 > 1) is the known star outlier: the commutant mode is not dark enough, so an oscillating
    band-edge mode wins and the star un-freezes. The frozen seam IS the structural ceiling, dynamically."""
    print("\n[Stage 2 / the threshold] star frozen at all Q  <=>  g2 = 4/(N-1) <= 1:")
    res = {}
    for N in (4, 5, 6, 7, 8):
        sm = max(im_of(N, "star", Q)[0] for Q in QS)
        res[N] = sm
        g2 = 4.0 / (N - 1)
        print(f"     N={N}: g2=4/(N-1)={g2:.3f}  star max|Im| over Q = {sm:.2e}  -> "
              f"{'UN-FREEZES' if sm > 1e-2 else 'frozen at all Q'}")
    # N=4 (g2=4/3>1) un-freezes; N>=5 (g2<=1) frozen at all Q.
    assert res[4] > 1e-2, f"N=4 star should un-freeze (g2=4/3>1, the (2,2)/K_4 outlier), got {res[4]:.2e}"
    for N in (5, 6, 7, 8):
        assert res[N] < FROZEN, f"N={N} (g2={4/(N-1):.3f}<=1) star should be frozen at all Q, got {res[N]:.2e}"
    print("     => frozen seam <=> g2 = 4/(N-1) <= 1 (N>=5). The dynamic reading of the structural ceiling:")
    print("     the survivor is the darkest commutant (1,1) coherence only when it undercuts the -2g floor.")


if __name__ == "__main__":
    _pin_at_Q15()
    _sweep_and_gate()
    _cross_N()
    print("\nAll stages passed: the star has its own frozen seam -- its survivor never un-freezes for")
    print("N>=5, and that threshold IS the structural ceiling g2 = 4/(N-1) <= 1 read dynamically (N=4,")
    print("g2=4/3>1, is the known outlier that un-freezes).")
