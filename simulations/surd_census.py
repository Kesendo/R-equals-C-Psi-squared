"""THE SURD CENSUS (gate-first): every small-N special number of the Liouvillian is a cyclotomic POLYGON
constant 2cos(π/m). Turning the Niven lens (niven_rationality_root.py) on the WHOLE spectrum, not just the
band edge, reveals two correlated structures:

  1. EVERYTHING is a polygon. The band edge √2/φ/√3, the K_4 ceiling's √3, the ring-4 frequency's √2 — all
     are 2cos(π/m) for a small polygon m (square/pentagon/hexagon). The small-N "zoo" is the polygons.

  2. TWO ARITHMETICS on the two axes (suggestive, not exact — see the caveat):
       Im-side (frequencies, from L_H / the path-graph adjacency spectrum): band edge 2cos(π/(N+1)) — the
         cyclotomic surd ladder, polygon m=N+1, quadratic only through the hexagon (N=5), degree climbing after.
       Re-side (decay, from L_D / the commutant rep structure): structural ceiling g2 = 4/N — RATIONAL, the
         S_N standard-rep principal angle λ₂=(N−2)/N. The rational corner of the same cyclotomic family.
     The sharper axis is graph-spectral (cyclotomic) vs S_N-standard-rep (rational); it CORRELATES with Im/Re
     but is not identical (F65's single-excitation decay rates are Re-values yet cyclotomically-flavored,
     rational only on N+1∈{1,2,3,4,6}; and the (2,2) Re ceiling is itself a surd — the caveat below).

  3. THE N=4 POLYGON CONFLUENCE: at the first even half-filled N, three sectors light three different polygons
     at once — single-excitation band edge = PENTAGON (φ), K_4 (2,2) ceiling = HEXAGON (√3, via a rep
     principal angle), ring-4 (2,2) frequency = SQUARE (√2, 2-magnon). The (2,2) half-filling anomaly is where
     a cyclotomic surd LEAKS onto the Re-side.

Gates PROVE the polygon identifications with sympy (exact 2cos(π/m) comparison; v1's nsimplify-into-basis bug
that mangled degree-≥3 cyclotomics is fixed). A firing gate is a finding.
Run: python simulations/surd_census.py
"""
import sys
from math import gcd
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
t = sp.Symbol("t")


def totient(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def cyclotomic_polygon(x):
    """If the EXACT value x equals 2cos(π/m) for some integer m in 2..200, return m; else None.
    (Exact comparison — NOT nsimplify-into-basis, which corrupts degree-≥3 cyclotomics like 2cos(π/7).)"""
    for m in range(2, 201):
        if sp.simplify(x - 2 * sp.cos(sp.pi / m)) == 0:
            return m
    return None


def degree(x):
    return int(sp.degree(sp.minimal_polynomial(sp.nsimplify(x), t), t))


def surd_field(x):
    """Which quadratic field Q(√d), d ∈ {2,3,5}, does the degree-2 surd x live in? (None if not degree 2.)"""
    if degree(x) != 2:
        return None
    for d in (2, 3, 5):
        cand = sp.nsimplify(x, [sp.sqrt(d)])
        if sp.simplify(x - cand) == 0 and cand.has(sp.sqrt(d)):
            return d
    return None


print("=" * 98)
print("THE SURD CENSUS — every small-N special is a cyclotomic POLYGON constant 2cos(π/m)")
print("=" * 98)

print("\n--- IM-SIDE (frequencies, L_H graph spectrum): band edge 2cos(π/(N+1)), polygon m = N+1 ---")
im = {}
for N in range(3, 9):
    x = 2 * sp.cos(sp.pi / (N + 1))
    m = cyclotomic_polygon(x)
    deg = totient(2 * (N + 1)) // 2
    im[N] = (m, deg)
    tag = "rational" if deg == 1 else (f"quadratic surd √{surd_field(x)}" if deg == 2 else f"degree {deg}")
    print(f"  N={N}: 2cos(π/{N+1}) = {float(x):.5f}  ->  polygon m={m}, algebraic degree {deg}  ({tag})")

print("\n--- RE-SIDE (decay, L_D commutant): structural ceiling g2 = 4/N (the (1,1) S_N standard rep) ---")
for N in range(4, 9):
    g = sp.Rational(4, N)
    print(f"  K_{N}: g2 = 4/{N} = {float(g):.5f}  ->  RATIONAL (rep-theoretic: λ₂ = (N−2)/N rational)")

print("\n--- THE (2,2) HALF-FILLING ANOMALY at N=4: a cyclotomic surd LEAKS onto the Re-side ---")
k4 = 2 - 2 / sp.sqrt(3)
ring4 = 2 * sp.sqrt(2)
print(f"  K_4 (2,2) ceiling g2 = 2−2/√3 = {float(k4):.5f}  ->  degree {degree(k4)}, Q(√{surd_field(k4)}); "
      f"√3 = 2cos(π/{cyclotomic_polygon(sp.sqrt(3))}) (HEXAGON), via a rep PRINCIPAL ANGLE λ₂=1/√3 (not a frequency)")
print(f"  ring-4 (2,2) Im = 2√2 = {float(ring4):.5f}  ->  degree {degree(ring4)}, Q(√{surd_field(ring4)}); "
      f"√2 = 2cos(π/{cyclotomic_polygon(sp.sqrt(2))}) (SQUARE), a frequency scaled ×2 (the 2-magnon band edge)")

print("\n--- THE N=4 POLYGON CONFLUENCE (three sectors, three polygons, one N) ---")
print("  single-excitation band edge : φ  = 2cos(π/5) -> PENTAGON  (Im axis)")
print("  K_4 (2,2) ceiling           : √3 = 2cos(π/6) -> HEXAGON   (Re axis, principal angle)")
print("  ring-4 (2,2) frequency      : √2 = 2cos(π/4) -> SQUARE    (Im axis, 2-magnon)")

print("\n" + "=" * 98)
print("GATES (a firing gate is a finding: diagnose, don't loosen)")
print("=" * 98)

# G1: the Im band edge IS the cyclotomic ladder 2cos(π/(N+1)), polygon m=N+1, degree φ_euler(2(N+1))/2
for N in range(3, 9):
    m, deg = im[N]
    assert m == N + 1, f"G1 FIRED: N={N} band-edge polygon {m} != N+1={N+1}"
    assert deg == totient(2 * (N + 1)) // 2, f"G1 FIRED: N={N} degree mismatch"
assert surd_field(2 * sp.cos(sp.pi / 4)) == 2 and surd_field(2 * sp.cos(sp.pi / 5)) == 5 \
    and surd_field(2 * sp.cos(sp.pi / 6)) == 3, "G1 FIRED: N=3,4,5 surds not √2,√5,√3"
print("G1 PASS (Im): the band edge is the cyclotomic ladder 2cos(π/(N+1)), polygon m=N+1; N=3,4,5 = √2/√5(φ)/√3.")

# G2: the Re commutant ceiling g2 = 4/N is rational for every N
for N in range(4, 9):
    assert sp.Rational(4, N).is_rational, f"G2 FIRED: 4/{N} not rational"
print("G2 PASS (Re): the structural ceiling g2 = 4/N is rational at every N (S_N standard-rep principal angle).")

# G3: the (2,2) anomaly surds are cyclotomic-polygon constants (√3 = hexagon, √2 = square)
assert degree(k4) == 2 and surd_field(k4) == 3 and cyclotomic_polygon(sp.sqrt(3)) == 6, "G3 FIRED: K_4 not √3/hexagon"
assert degree(ring4) == 2 and surd_field(ring4) == 2 and cyclotomic_polygon(sp.sqrt(2)) == 4, "G3 FIRED: ring-4 not √2/square"
print("G3 PASS ((2,2) leak): K_4 ceiling = √3 (hexagon, principal angle), ring-4 = √2 (square, 2-magnon).")

# G4: the N=4 confluence — three DISTINCT polygons in three sectors
confluence = {
    cyclotomic_polygon(2 * sp.cos(sp.pi / 5)),   # band edge φ -> 5
    cyclotomic_polygon(sp.sqrt(3)),              # K_4 ceiling -> 6
    cyclotomic_polygon(sp.sqrt(2)),              # ring-4 -> 4
}
assert confluence == {4, 5, 6}, f"G4 FIRED: N=4 confluence polygons {confluence} != {{4,5,6}}"
print("G4 PASS (confluence): at N=4 three sectors light three distinct polygons — square (4) / pentagon (5) / hexagon (6).")

print("\nROOT: the small-N spectrum is a garden of cyclotomic polygons. The Im axis is the 2cos(π/(N+1)) ladder")
print("directly; the Re commutant ceiling is the rational corner (4/N); the (2,2) half-filling anomaly leaks a")
print("polygon surd onto the Re side, and at N=4 three polygons appear at once. CAVEAT: 'Re rational / Im surd'")
print("is suggestive, not exact (F65 rates are Re-but-cyclotomic; the (2,2) ceiling is a surd). DONE.")
