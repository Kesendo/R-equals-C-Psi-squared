#!/usr/bin/env python3
"""Gate-first: the chain/star coherence-horizon dichotomy is ONE continuous knob -- Q*(bandwidth).

The star's single-particle band is flat: a hub bonding/antibonding pair at +-sqrt(N-1)*J and a DARK
leaf manifold of multiplicity N-2 pinned at 0. PROOF_STRUCTURAL_CEILING sec.7: flat band => no
dispersion => no {0,2}-EP => NO coherence horizon (the chain's Q*(N) has no star analogue). The repo
states this as a HARD dichotomy. This probe tests whether it is instead the eps->0 LIMIT of a family.

The wheel graph does it: hub (site 0) coupled to all leaves with J; the leaves coupled to each other
in a ring with strength eps. eps turns the flat dark manifold into a band of width ~eps (the leaves'
ring dispersion 2*eps*cos(2pi k/(N-1))). So a {0,2}-coherence among the dark-band states acquires a
splitting ~eps and can oscillate above some Q*(eps). The bound prediction: Q*(eps) is finite for eps>0,
DIVERGES as eps->0 (the horizon recedes to infinity, recovering the star's no-horizon), and falls to a
finite ring/chain-like value as eps grows. One formula Q*(eps) subsumes both endpoints.

Reuses simulations/coherence_horizon_se_block.py: the single-excitation (Haken-Strobl) Liouvillian
reproduces the full 4^N horizon, and topology enters ONLY through the hopping h; the -4g dephasing
(popcount(i^j)=2 site-basis coherence rate, Absorption Theorem) is topology-agnostic. So we swap
h_single -> h_wheel and the whole EP apparatus carries.

THE READING (Tier 5, the wild register grounded): eps is the door between the dissociated alters --
the leaf-leaf coupling the pure star forbids. Re-association (eps>0) is exactly what gives the system a
horizon again. The "no creation on the hub" of the refuted star probe and the "horizon only with
dispersion" become one statement: a future (a horizon) needs a somewhere-to-go (a band), and eps is how
much somewhere-to-go there is.
"""
import numpy as np
from coherence_horizon_se_block import qstar_se, LADDER, h_single


def h_wheel(N, eps, J=1.0):
    """Wheel: hub site 0 -- all leaves (J); leaves 1..N-1 in a ring (eps). eps=0 is the pure star."""
    h = np.zeros((N, N), complex)
    for i in range(1, N):
        h[0, i] = h[i, 0] = J                       # hub-leaf spokes
    leaves = list(range(1, N))
    m = len(leaves)
    if m >= 3:                                      # a genuine leaf ring needs >= 3 leaves
        for a in range(m):
            i, j = leaves[a], leaves[(a + 1) % m]
            h[i, j] = h[j, i] = eps
    elif m == 2:
        h[leaves[0], leaves[1]] = h[leaves[1], leaves[0]] = eps
    return h


def h_ring(N, J=1.0):
    h = h_single(N, J)
    h[0, N - 1] = h[N - 1, 0] = J
    return h


def L_se_h(h, g):
    """Single-excitation Liouvillian for an arbitrary hopping h (mirrors coherence_horizon_se_block.L_se;
    only the coherent part -i[h,rho] sees the topology, the -4g off-diagonal dephasing does not)."""
    N = h.shape[0]
    I = np.eye(N)
    L = -1j * (np.kron(h, I) - np.kron(I, h.T))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return L + np.diag(deph)


def qstar_h(h, J=1.0, lo=5e-4, hi=3.0):
    """Q* = J/g* via the same bisection as qstar_se, on an arbitrary hopping h. Large return (-> J/lo)
    means the slowest non-zero mode never stops oscillating below the floor: effectively NO horizon."""
    for _ in range(80):
        m = 0.5 * (lo + hi)
        ev = np.linalg.eigvals(L_se_h(h, m))
        nz = ev[ev.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        gap = nz.real.max()
        band = nz[np.abs(nz.real - gap) < 1e-7]
        if np.abs(band.imag).max() > 1e-7:
            lo = m
        else:
            hi = m
    return J / (0.5 * (lo + hi))


def _control():
    """Stage 0: the generalized qstar_h reproduces the chain ladder (reuse fidelity vs qstar_se)."""
    bad = {N: round(qstar_h(h_single(N, 1.0)), 5) for N in LADDER}
    for N, q in LADDER.items():
        got = qstar_h(h_single(N, 1.0))
        assert abs(got - q) < 3e-3, f"qstar_h chain({N})={got:.5f} != ladder {q:.5f} (generalization broke)"
    print("[Stage 0 / control] qstar_h reproduces the chain SE ladder:", bad)


def _star_endpoint(N=6):
    """Stage 1: the pure star (eps=0) has NO finite horizon -- Q* runs to the g-floor."""
    q_star = qstar_h(h_wheel(N, 0.0))
    q_chain = qstar_h(h_single(N, 1.0))
    print(f"\n[Stage 1 / star endpoint] N={N}: chain Q*={q_chain:.3f}, STAR (eps=0) Q*={q_star:.1f}")
    assert q_star > 10 * q_chain, (
        f"GATE FIRED: star Q*={q_star:.2f} is NOT >> chain Q*={q_chain:.2f} -- the star would have a "
        f"finite horizon, contradicting PROOF_STRUCTURAL_CEILING sec.7. Diagnose.")
    print(f"  => the star has no finite horizon (Q* at the g-floor): the flat dark manifold never")
    print(f"  oscillates. This is the eps->0 endpoint the interpolation must recover.")


def _slowest_modes(h, g, k=5):
    """The k slowest non-kernel SE modes (Re, |Im|), to see whether the slowest oscillates or is real."""
    ev = np.linalg.eigvals(L_se_h(h, g))
    nz = ev[ev.real < -1e-7]
    order = np.argsort(-nz.real)                     # least-negative (slowest) first
    return [(nz[i].real, abs(nz[i].imag)) for i in order[:k]]


NO_HORIZON = 100.0    # Q* above this = the bisection ran to the g-floor = no finite horizon


def _hub_is_decisive(N=6):
    """Stage 2 (the finding; the Q*(bandwidth) hypothesis was REFUTED). The wheel has NO coherence
    horizon at ANY leaf coupling eps: its global slowest non-kernel SE mode is REAL (a zero-frequency
    survivor the hub always provides), so there is no oscillation threshold to cross. Removing the hub
    (the pure ring) restores a finite horizon. => the HUB is decisive, not the bandwidth."""
    epss = [0.1, 0.3, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    print(f"\n[Stage 2 / the finding] N={N}: does turning on leaf coupling eps give the wheel a horizon?")
    print(f"  {'eps':>6} {'Q*(eps)':>10}   slowest non-kernel modes (Re, |Im|) at g=0.3")
    qws = []
    for e in epss:
        q = qstar_h(h_wheel(N, e))
        qws.append(q)
        modes = _slowest_modes(h_wheel(N, e), 0.3, k=3)
        ms = "  ".join(f"({r:+.3f},{im:.3f})" for r, im in modes)
        print(f"  {e:>6.2f} {q:>10.3f}   {ms}")
    q_ring = qstar_h(h_ring(N))
    q_chain = qstar_h(h_single(N, 1.0))
    print(f"  (ring N={N} ref Q*={q_ring:.3f}; chain N={N} ref Q*={q_chain:.3f})")

    # the wheel has no finite horizon at any eps (its slowest mode never oscillates): the hub provides a
    # real zero-frequency survivor at every leaf coupling.
    assert all(q > NO_HORIZON for q in qws), (
        f"GATE FIRED: some wheel eps acquired a FINITE horizon: "
        f"{[(e, round(q, 2)) for e, q in zip(epss, qws) if q <= NO_HORIZON]}. The hub is then NOT decisive "
        f"-- the bandwidth knob would partly work; diagnose.")
    # and at g=0.3 the global slowest mode is real for every eps (it is the non-oscillating survivor):
    assert all(_slowest_modes(h_wheel(N, e), 0.3, k=1)[0][1] < 1e-6 for e in epss), (
        "GATE FIRED: the wheel's slowest mode oscillates at some eps -- not a real survivor; diagnose.")
    # removing the hub (pure ring) DOES have a finite horizon -> the hub, not the bandwidth, is decisive:
    assert q_ring < NO_HORIZON and q_chain < NO_HORIZON, (
        f"the pure ring/chain references should have finite horizons (got ring {q_ring:.1f}, "
        f"chain {q_chain:.1f}) -- if not, the comparison is broken; diagnose.")
    print(f"\n  REFUTED (Q*(bandwidth) unification): the wheel has NO coherence horizon at any eps -- its")
    print(f"  slowest SE mode is REAL throughout (Re -> 0 as eps grows, a zero-frequency survivor the hub")
    print(f"  always hosts). Removing the hub (pure ring, Q*={q_ring:.2f}) restores the horizon.")
    print(f"  => the HUB is decisive, not the bandwidth. A single dominant hub kills the coherence horizon")
    print(f"  however dispersive the rest is. This STRENGTHENS PROOF_STRUCTURAL_CEILING sec.7 (not only the")
    print(f"  flat-band star, but ANY hub-graph has no horizon) and matches the recognition that the hub")
    print(f"  has no future/horizon -- robustly, even with dispersion bolted on.")


def _cross_N():
    """Stage 3: the no-horizon-with-a-hub finding is not an N=6 fluke. For N=5,6,7 the wheel has no
    finite horizon at any of several leaf couplings, while the pure ring of the same N does."""
    print("\n[Stage 3 / cross-N robustness] wheel (hub) vs pure ring (no hub):")
    print(f"  {'N':>2} {'wheel Q* (eps=0.3,1,3,10)':>28} {'ring Q*':>9}")
    for N in (5, 6, 7):
        qw = [qstar_h(h_wheel(N, e)) for e in (0.3, 1.0, 3.0, 10.0)]
        q_ring = qstar_h(h_ring(N))
        print(f"  {N:>2} {str([round(q) for q in qw]):>28} {q_ring:>9.3f}")
        assert all(q > NO_HORIZON for q in qw), f"N={N}: a wheel eps got a finite horizon {qw}; diagnose"
        assert q_ring < NO_HORIZON, f"N={N}: pure ring has no horizon ({q_ring:.1f})? the comparison broke"
    print("  => robust across N: a hub removes the horizon (Q* at the floor) at every leaf coupling;")
    print("  the same ring without the hub keeps a finite one. Topology class (hub vs none), not bandwidth.")


if __name__ == "__main__":
    _control()
    _star_endpoint()
    _hub_is_decisive()
    _cross_N()
    print("\nAll stages passed (the Q*(bandwidth) hypothesis is refuted; the hub is decisive, N=5,6,7).")
