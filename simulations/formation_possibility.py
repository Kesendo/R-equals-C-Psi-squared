#!/usr/bin/env python3
"""A formation-possibility scanner for bound complexes.

Reframed from a triple-alpha mirror into the tool it actually is. Given a chain (the substrate),
scan the binding strength and read WHERE a k-body bound complex can form, and where the marginal,
near-threshold window sits: the regime where a resonance (like carbon's Hoyle state, barely above
the Be-8 + alpha threshold) would live. It does not reproduce one fine-tuned point; it maps the
landscape of possibility, so you can see what is generic (the window itself) and where a fine-tuned
resonance would have to land.

The order parameter is the CLUSTERING of the energetically-bound complex: the weight, on the
fully-clustered (k-adjacent) configurations, of the actual bound eigenstate (the band-edge state,
which for Delta > 0 is the highest-energy one, since adjacency is energetically high here). An
earlier version took the argmax over all eigenstates of this weight; that proxy was gauge-junk, it
picked a tightly-localized mid-spectrum state rather than the bound complex, and falsely reported
the pair as deeply bound at every Delta. Reading the band-edge state directly fixes it.
  clustering -> 0    the complex does not form (extended, dissolves into the free continuum)
  clustering ~ 0.5   MARGINAL: barely bound, large extent, near the threshold (the resonance window)
  clustering -> 1    deeply bound, tightly clustered (a ground-state-like complex)

Carbon's Hoyle resonance sits in the marginal window. The scan finds that window in parameter space.

  python simulations/formation_possibility.py [N] [k] [j2]
"""
import sys
from itertools import combinations

import numpy as np


def Hn(N, k, delta, j2=0.0, ring=False):
    """k-excitation sector of the XXZ chain (or ring): nearest-neighbour XY hopping (amplitude 2) +
    Delta*ZZ (diagonal) + j2 next-nearest-neighbour XY hopping (amplitude 2*j2). The j2 term
    breaks Bethe integrability. ring=True closes the loop (the wrap bond N-1 <-> 0 and its NNN
    partners), the aromatic substrate (benzene C6, etc.)."""
    basis = list(combinations(range(N), k))
    index = {b: i for i, b in enumerate(basis)}
    M = len(basis)
    H = np.zeros((M, M))

    nn = [(l, l + 1) for l in range(N - 1)]
    nnn = [(l, l + 2) for l in range(N - 2)]
    if ring and N >= 3:
        nn.append((N - 1, 0))                       # the wrap bond closes the ring
        nnn += [(N - 2, 0), (N - 1, 1)]             # next-nearest across the seam

    def hop(S, l, m, amp, col):
        if (l in S) == (m in S):
            return
        nS = set(S)
        if l in S:
            nS.discard(l); nS.add(m)
        else:
            nS.discard(m); nS.add(l)
        H[index[tuple(sorted(nS))], col] += amp

    for bi, occ in enumerate(basis):
        S = set(occ)
        e = 0.0
        for (l, m) in nn:
            sl = -1 if l in S else 1
            sm = -1 if m in S else 1
            e += sl * sm
        H[bi, bi] = delta * e
        for (l, m) in nn:
            hop(S, l, m, 2.0, bi)
        if j2 != 0.0:
            for (l, m) in nnn:
                hop(S, l, m, 2.0 * j2, bi)
    return H, basis


def adjacency_links(basis, N, ring=False):
    """Per config, the number of adjacent links among its excitations; the fully clustered
    k-string has the maximum, k-1. On a ring, sites 0 and N-1 count as adjacent (the seam), so a
    cluster wrapping the seam is still recognized as clustered."""
    out = []
    for occ in basis:
        links = sum(1 for a, b in zip(occ, occ[1:]) if b == a + 1)
        if ring and (0 in occ) and (N - 1 in occ):
            links += 1
        out.append(links)
    return np.array(out)


def formation_order(N, k, delta, j2=0.0, ring=False):
    """The formation order parameter: the weight on the fully-clustered (k-string) configs carried
    by the actual bound complex, the band-edge eigenstate. For Delta > 0 adjacency is energetically
    high, so the bound complex is the HIGHEST-energy eigenstate; reading it directly avoids the
    argmax proxy that mistook a localized mid-spectrum state for the complex."""
    H, basis = Hn(N, k, delta, j2, ring)
    full = adjacency_links(basis, N, ring) == (k - 1)   # the k-adjacent (k-string) configurations
    w, V = np.linalg.eigh(H)
    p = np.abs(V[:, -1]) ** 2                        # band-edge (highest-E) state = the bound complex
    return float(p[full].sum())


def regime(clustering):
    if clustering < 0.35:
        return "unformed (extended, dissolves into the continuum)"
    if clustering < 0.65:
        return "MARGINAL near-threshold (resonance window)"
    return "deeply bound (tightly clustered complex)"


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    j2 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    ring = len(sys.argv) > 4 and sys.argv[4].lower() == "ring"
    topo = "ring" if ring else "chain"
    print(f"Formation-possibility scan: {k}-body bound complex on an N={N} {topo}"
          f"{f', j2={j2}' if j2 else ''}, vs binding Delta\n")
    print(f"  {'Delta':>6}  {'cluster':>7}   regime")
    deltas = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0]
    marginal = []
    for d in deltas:
        c = formation_order(N, k, d, j2, ring)
        if 0.35 <= c < 0.65:
            marginal.append(d)
        print(f"  {d:>6.2f}  {c:>7.3f}   {regime(c)}")
    if marginal:
        print(f"\n  near-threshold (resonance) window: Delta ~ {min(marginal):.2f} .. {max(marginal):.2f}")
    else:
        print("\n  no marginal window on this grid (refine Delta to find the crossover)")
    print("  Below the window the complex does not form; above it, it is deeply bound and needs no")
    print("  resonance. A near-threshold resonance (carbon's Hoyle state) lives in the marginal window.")


if __name__ == "__main__":
    main()
