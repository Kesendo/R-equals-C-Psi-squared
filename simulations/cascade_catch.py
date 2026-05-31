#!/usr/bin/env python3
"""The cascade step: does catching form the next element?

A slow bound pair sits on the left; a fast free excitation comes in from the right and catches
it. The thesis: the slowness of the new element is what lets the next one catch it, and the
catching IS the next birth. So we ask, in the three-excitation sector of an XXZ chain: when the
free excitation reaches the bound pair, does a three-excitation complex (the next, slower
element) form and travel on, or do they just pass through?

We track, over time, how tightly the three excitations cluster:
  - mean span  <max_site - min_site>  of the three (small = clustered, large = spread)
  - P(tight)   probability all three sit within 3 sites (a 3-complex)
The test is the LATE-TIME value after the encounter: free (Delta=0) should disperse again;
interacting (Delta>0) should leave a lasting tight lump, a complex that formed and moved on.
"""
import sys

import numpy as np

BLOCKS = "▁▂▃▄▅▆▇█"


def three_exc_basis(N):
    return [(a, b, c) for a in range(N) for b in range(a + 1, N) for c in range(b + 1, N)]


def H3(N, delta, j2=0.0):
    """Hamiltonian in the 3-excitation sector: nearest-neighbour XY hopping (amplitude 2) +
    Delta*ZZ (diagonal) + j2 next-nearest-neighbour XY hopping (amplitude 2*j2). The j2 term
    breaks Bethe integrability: the magnon scattering becomes inelastic, so the conserved string
    content is broken and a 3-complex can form from a collision."""
    basis = three_exc_basis(N)
    index = {b: i for i, b in enumerate(basis)}
    M = len(basis)
    H = np.zeros((M, M))

    def hop(S, l, m, amp, col):
        """Move the single occupied site between l and m (if exactly one is occupied)."""
        if (l in S) == (m in S):
            return
        newS = set(S)
        if l in S:
            newS.discard(l); newS.add(m)
        else:
            newS.discard(m); newS.add(l)
        H[index[tuple(sorted(newS))], col] += amp

    for bi, occ in enumerate(basis):
        S = set(occ)
        e = 0.0
        for l in range(N - 1):
            sl = -1 if l in S else 1
            sr = -1 if (l + 1) in S else 1
            e += sl * sr
        H[bi, bi] = delta * e
        for l in range(N - 1):
            hop(S, l, l + 1, 2.0, bi)               # nearest-neighbour hop
        if j2 != 0.0:
            for l in range(N - 2):
                hop(S, l, l + 2, 2.0 * j2, bi)       # next-nearest-neighbour hop (breaks integrability)
    return H, basis


def initial_state(N, basis, cp, cf, wp, wf, k_free):
    """Bound pair (two adjacent excitations) near cp, plus a free excitation near cf with
    leftward momentum k_free. Only pair+separated-free configs get amplitude at t=0."""
    psi = np.zeros(len(basis), complex)
    for i, (a, b, c) in enumerate(basis):
        pair_site, free_site = None, None
        if b == a + 1 and c >= b + 2:       # pair (a,a+1), free c on the right
            pair_site, free_site = a, c
        elif c == b + 1 and a <= b - 2:     # pair (b,b+1), free a on the left
            pair_site, free_site = b, a
        if pair_site is None:
            continue
        amp_pair = np.exp(-((pair_site - cp) ** 2) / (2 * wp ** 2))
        amp_free = np.exp(-((free_site - cf) ** 2) / (2 * wf ** 2)) * np.exp(-1j * k_free * free_site)
        psi[i] = amp_pair * amp_free
    psi /= np.linalg.norm(psi)
    return psi


def sparkline(ys):
    ys = np.asarray(ys, float)
    lo, hi = ys.min(), ys.max()
    if hi - lo < 1e-12:
        return BLOCKS[0] * len(ys)
    idx = np.clip(((ys - lo) / (hi - lo) * 7).round().astype(int), 0, 7)
    return "".join(BLOCKS[i] for i in idx)


def run(N, delta, j2, ts):
    H, basis = H3(N, delta, j2)
    spans = np.array([c - a for (a, b, c) in basis], float)
    tight_mask = spans <= 2.0  # all three within 3 sites
    psi0 = initial_state(N, basis, cp=2.0, cf=N - 3.0, wp=1.0, wf=1.4, k_free=np.pi / 2)
    w, V = np.linalg.eigh(H)
    coeff = V.T.conj() @ psi0
    mean_span, p_tight = [], []
    for t in ts:
        p = np.abs(V @ (np.exp(-1j * w * t) * coeff)) ** 2
        mean_span.append(float((p * spans).sum()))
        p_tight.append(float(p[tight_mask].sum()))
    return np.array(mean_span), np.array(p_tight)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    ts = np.linspace(0.0, 6.0, 120)
    late = slice(len(ts) * 2 // 3, None)  # last third = after the encounter
    print(f"Cascade catch with integrability breaking, N={N}, 3-excitation sector")
    print(f"  slow bound pair + fast free excitation collide; j2 = next-nearest hop breaks Bethe")
    print(f"  integrability. Does a 3-complex form (P_tight late rising above the integrable floor)?\n")
    # the free, integrable floor
    _, pt0 = run(N, 0.0, 0.0, ts)
    floor = pt0[late].mean()
    print(f"  free integrable floor (Delta=0, j2=0): P_tight late = {floor:.3f}\n")
    print(f"  {'Delta':>6}  {'j2':>5}  {'P_tight late':>12}  {'span late':>9}  {'vs floor':>8}  P(tight)(t)")
    for delta in [1.0, 2.0]:
        for j2 in [0.0, 0.3, 0.6, 1.0]:
            ms, pt = run(N, delta, j2, ts)
            pt_late, ms_late = pt[late].mean(), ms[late].mean()
            print(f"  {delta:>6.2f}  {j2:>5.2f}  {pt_late:>12.3f}  {ms_late:>9.2f}  {pt_late-floor:>+8.3f}  {sparkline(pt)}")
        print()
    print(f"  reading: at j2=0 (integrable) the string content is conserved, the collision is elastic,")
    print(f"  no 3-complex forms (vs-floor near 0 or negative). If breaking integrability (j2>0) lifts")
    print(f"  P_tight clearly above the floor, the next element is born from the catching where order breaks.")


if __name__ == "__main__":
    main()
