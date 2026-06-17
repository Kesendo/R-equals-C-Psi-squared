"""THE N=4 POLYGON CONFLUENCE IS UNIQUE (gate-first): why N=4 keeps surfacing as the special case.

The surd census (surd_census.py) found that at N=4 three sectors light three distinct QUADRATIC polygons
(pentagon φ / hexagon √3 / square √2). This probe asks the (b) question — does N=6 (the next even, half-
fillable N) reproduce it? — and answers NO. The clean confluence is an N=4-only event, and both of its
ingredients trace to the SAME small-N cyclotomic accident:

  (1) BAND-EDGE DEGREE. The chain single-excitation band edge 2cos(π/(N+1)) is a quadratic surd only for
      N+1 ∈ {4,5,6} (Niven / niven_rationality_root.py). Among EVEN N, that is ONLY N=4 (N+1=5, the pentagon,
      φ — the unique even-N golden). N=6 → N+1=7, the heptagon, is CUBIC (degree 3): no quadratic surd.

  (2) THE Re-SIDE SURD LEAK. At N=4 the (1,1) commutant ceiling ladder 4/N hits 1 and vacates the sub-floor,
      so the (2,2) half-filling sector BECOMES the ceiling (2−2/√3, a hexagon surd via a rep principal angle).
      At N=6 the (1,1) ceiling is already 4/6 = 2/3 < 1 (rational) and IS the floor; the (2,2)/(3,3) surds sit
      ABOVE it, so nothing leaks onto the Re ceiling. The leak needs the ladder to vacate — an N=4 accident.

So N=4 is the unique even, half-fillable N whose band edge is a quadratic golden AND whose ceiling ladder
vacates the sub-floor. The pentagon is the hinge.

Run: python simulations/n4_polygon_uniqueness.py
"""
import sys
import numpy as np
import sympy as sp
from itertools import combinations

sys.stdout.reconfigure(encoding="utf-8")
NULLTOL, NZTOL = 1e-7, 1e-7
t = sp.Symbol("t")


def bonds(topo, n):
    if topo == "chain":
        return [(i, i + 1) for i in range(n - 1)]
    if topo == "ring":
        return [(i, (i + 1) % n) for i in range(n)]
    if topo == "complete":
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    raise ValueError(topo)


def sector_states(n, p):
    return [sum(1 << i for i in c) for c in combinations(range(n), p)]


def sector_H(topo, n, states):
    idx = {s: a for a, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for a, s in enumerate(states):
        for (i, j) in bonds(topo, n):
            if ((s >> i) & 1) != ((s >> j) & 1):
                H[idx[s ^ (1 << i) ^ (1 << j)], a] += 1.0
    return H


def commutant_darkest(topo, n, p, q):
    """Smallest nonzero N_XY eigenvalue on the ad_H kernel of the (p,q) sector (the high-Q g2 contribution)."""
    A, B = sector_states(n, p), sector_states(n, q)
    Hp, Hq = sector_H(topo, n, A), sector_H(topo, n, B)
    na, nb = len(A), len(B)
    adH = np.kron(Hp, np.eye(nb)) - np.kron(np.eye(na), Hq.T)
    w, V = np.linalg.eigh(adH)
    diag = np.array([bin(A[a] ^ B[b]).count("1") for a in range(na) for b in range(nb)], float)
    cols = [k for k in range(len(w)) if abs(w[k]) < NULLTOL]
    if not cols:
        return None
    U = V[:, cols]
    ev = np.linalg.eigvalsh(U.T @ (diag[:, None] * U))
    pos = ev[ev > NZTOL]
    return float(pos.min()) if len(pos) else None


def band_edge_degree(N):
    return int(sp.degree(sp.minimal_polynomial(2 * sp.cos(sp.pi / (N + 1)), t), t))


def is_rational(x):
    r = sp.nsimplify(x, tolerance=1e-9)
    return r.is_rational and abs(float(r) - x) < 1e-7


def darkest_sector(topo, n, fillings):
    vals = {(p, q): commutant_darkest(topo, n, p, q) for (p, q) in fillings}
    floor = min(v for v in vals.values() if v is not None)
    where = [k for k, v in vals.items() if v is not None and abs(v - floor) < 1e-7]
    return vals, floor, where


def adjacency(topo, n):
    A = np.zeros((n, n))
    for (i, j) in bonds(topo, n):
        A[i, j] = A[j, i] = 1.0
    return A


print("=" * 92)
print("THE N=4 POLYGON CONFLUENCE IS UNIQUE — does N=6 reproduce it? (no)")
print("=" * 92)

print("\n(1) chain band-edge algebraic degree (quadratic surd only for N+1 in {4,5,6}):")
even_quadratic = []
for N in (2, 4, 6, 8):
    d = band_edge_degree(N)
    if d == 2:
        even_quadratic.append(N)
    print(f"  N={N} (N+1={N+1}): degree {d}  {'<- quadratic surd' if d == 2 else ('rational' if d == 1 else 'higher-degree')}")

print("\n(2) commutant ceiling (darkest sector) — does a surd leak onto the Re floor?")
c4_vals, c4_floor, c4_where = darkest_sector("complete", 4, [(1, 1), (2, 2)])
c6_vals, c6_floor, c6_where = darkest_sector("complete", 6, [(1, 1), (2, 2), (3, 3)])
print(f"  complete-4: ceiling = {c4_floor:.6f} at {c4_where}  ({'rational' if is_rational(c4_floor) else 'SURD'})  "
      f"<- the (2,2) leak (2-2/sqrt3)")
print(f"  complete-6: ceiling = {c6_floor:.6f} at {c6_where}  ({'rational' if is_rational(c6_floor) else 'SURD'})  "
      f"<- rational (1,1)=2/3 floor; (2,2)/(3,3) surds sit ABOVE it (no leak)")

print("\n(3) ring half-filling — surd (N=4) vs fully rational (N=6):")
r6_vals, r6_floor, r6_where = darkest_sector("ring", 6, [(1, 1), (2, 2), (3, 3)])
ring6_band = np.sort(np.linalg.eigvalsh(adjacency("ring", 6)))
print(f"  ring-6 commutant: all sectors {r6_floor:.6f} ({'rational' if is_rational(r6_floor) else 'surd'}); "
      f"band = {np.round(ring6_band, 3)} (integer hexagon)")

print("\n" + "=" * 92)
print("GATES (a firing gate is a finding: diagnose, don't loosen)")
print("=" * 92)

# G1: among even N, ONLY N=4 has a quadratic-surd band edge (the unique even-N golden, the pentagon)
assert even_quadratic == [4], f"G1 FIRED: even N with quadratic band edge = {even_quadratic}, expected [4]"
assert band_edge_degree(4) == 2 and band_edge_degree(6) == 3, "G1 FIRED: N=4 not quadratic or N=6 not cubic"
print("G1 PASS: among even N, ONLY N=4 has a quadratic-surd band edge (pentagon φ); N=6 is the cubic heptagon.")

# G2: the Re-side surd leak happens at N=4 (the (2,2) ceiling is a surd) but NOT at N=6 (rational (1,1) floor)
assert c4_where == [(2, 2)] and not is_rational(c4_floor), "G2 FIRED: complete-4 ceiling is not the (2,2) surd"
assert c6_where == [(1, 1)] and is_rational(c6_floor), "G2 FIRED: complete-6 ceiling is not the rational (1,1)"
print("G2 PASS: a surd leaks onto the Re ceiling at N=4 ((2,2) = 2-2/√3) but NOT at N=6 (rational (1,1) = 2/3 floor).")

# G3: the N=6 ring is fully rational (the integer hexagon band) — no frequency surd, unlike ring-4's 2√2
assert is_rational(r6_floor), "G3 FIRED: ring-6 commutant not rational"
assert np.allclose(ring6_band, np.round(ring6_band)), "G3 FIRED: ring-6 band not integer"
print("G3 PASS: ring-6 is the integer hexagon (band {±2,±1}, commutant 4/3) — fully rational, no surd (cf. ring-4's 2√2).")

print("\nROOT: N=4 is the UNIQUE polygon confluence. Both ingredients are the same cyclotomic accident: the")
print("pentagon (N+1=5) is the only even-N quadratic golden band edge, AND only at N=4 does the (1,1) ceiling")
print("ladder vacate the sub-floor so a half-filling surd becomes the ceiling. N=6 (heptagon, rational floors,")
print("integer hexagon ring) reproduces neither. The pentagon is the hinge. DONE.")
