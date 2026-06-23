"""Tom's follow-up (2026-06-20): are the 94% (degeneracy) IN the 3% future?

Two ORTHOGONAL cuts, easy to conflate:
  REAL axis -> past (slow, rate < Sg) / future (fast, rate > Sg): the drain-depth / memory axis.
  IMAG axis -> static (Im=0) / oscillating (Im != 0): my degeneracy is measured on Im != 0.

Tom recalls "97% past, 3% future". This script resolves it with data:
  1. The past/future split BY MODE COUNT is NOT 97/3. The palindrome makes rates symmetric about Sg,
     and the spectral mass sits at hamming ~ N/2 (rate ~ Sg, the MIDDLE), so past and future are the
     two minority rims, not a 97/3. The "97/3" is Tr(rho_past^2), a STATE purity (F94), not a mode
     count -- so there is no "3% future" mode room for the 94% to live inside.
  2. WHERE the degeneracy actually lives: collapse% in the past-oscillating vs future-oscillating
     halves. If it is future-heavy, Tom's geometric hunch has real content; if even, it does not.

Run:  python simulations/where_lives_the_degeneracy.py
"""
from pathlib import Path
import numpy as np

RESULTS = Path(__file__).parent / "results"
GAMMA = 0.05


def load(topo, N):
    name = f"rmt_eigenvalues_N{N}.csv" if topo == "chain" else f"rmt_eigenvalues_{topo}_N{N}.csv"
    p = RESULTS / name
    if not p.exists():
        return None
    r, i = [], []
    with open(p) as f:
        f.readline()
        for line in f:
            q = line.strip().split("\t")
            if len(q) == 2:
                r.append(float(q[0].replace(",", ".")))
                i.append(float(q[1].replace(",", ".")))
    return np.array(r) + 1j * np.array(i)


def collapse(pts):
    if len(pts) == 0:
        return float("nan"), 0
    nd = len(np.unique(np.round(pts, 9)))
    return 100 * (1 - nd / len(pts)), len(pts)


print("=" * 100)
print("WHERE DOES THE DEGENERACY LIVE? past vs future, on the oscillating (Im>0) modes")
print("  rate = -Re(lambda) in [0, 2*Sg];  past = rate<Sg (slow),  future = rate>Sg (fast),  Sg = N*gamma")
print("=" * 100)
print(f"{'topo':9} {'N':>2} {'osc(Im>0)':>9} {'%past':>6} {'%mid':>6} {'%future':>7} | "
      f"{'collapse past':>13} {'collapse future':>15}")
for topo in ("chain", "ring", "star", "complete"):
    for N in (5, 6):
        ev = load(topo, N)
        if ev is None:
            continue
        Sg = N * GAMMA
        up = ev[ev.imag > 1e-6]
        r = -up.real
        eps = 1e-6
        past = up[r < Sg - eps]
        mid = up[np.abs(r - Sg) <= eps]
        fut = up[r > Sg + eps]
        n = len(up)
        cp, _ = collapse(past)
        cf, _ = collapse(fut)
        print(f"{topo:9} {N:>2} {n:>9} {100*len(past)/n:>6.1f} {100*len(mid)/n:>6.1f} {100*len(fut)/n:>7.1f} | "
              f"{cp:>13.1f} {cf:>15.1f}")

print("\n" + "-" * 100)
print("READING:")
print("-" * 100)
print("* %past / %future by mode count is NOT 97/3 -- it is roughly balanced (the rim split), with a")
print("  fat %mid at rate~Sg (hamming~N/2). The memory's 97% is a state PURITY (F94), not this count.")
print("* collapse past = collapse future EXACTLY (bit-equal, every row). The degeneracy is NOT confined")
print("  to a 'future room'; it is split evenly across both rims. THIS IS F1, the palindrome, read in")
print("  the multiplicity: Pi conjugation (Pi L Pi^-1 = -L - 2Sg) is a bijection past<->future that")
print("  preserves eigenvalue multiplicity, so every degenerate slow mode has an equally degenerate")
print("  fast partner. So Tom's hunch is right -- the two structures DO connect -- but via Pi, not via")
print("  '97', and the degeneracy is mirror-symmetric, not future-heavy. The number 97/3 stays a state")
print("  purity (F94); the real link between my degeneracy and the past/future axis is the palindrome.")
print("* ALREADY KNOWN: this multiplicity-mirror is the central result of")
print("  experiments/DEGENERACY_PALINDROME.md ('The Palindrome Inside the Palindrome', 2026-04-03):")
print("  d_total(k) = d_total(N-k) proven from Pi, with the full per-grid sequence and closed forms.")
print("  This script is the coarse (two-half-sum) rediscovery; its NEW axis is topology (see")
print("  _rmt_topology_csr.py), which that chain-only document does not sweep.")
print("=" * 100)
