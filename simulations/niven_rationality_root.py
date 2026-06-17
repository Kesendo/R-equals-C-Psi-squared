"""THE NIVEN ROOT (gate-first): the number-theoretic root of the small-N specials, with TWO faces.

The single-excitation chain is cyclotomic: its rates and frequencies are trig values at angles k·π/(N+1).
Niven's theorem (the only rational cosines of rational-π angles are 0, ±1/2, ±1) puts a number-theoretic
ceiling on how clean these can be, and THAT is the shared root of "why small N is special":

  RE-FACE (the dissipator decay rates, F65: α_k/γ₀ = −2·sin²(kπ/(N+1)) = cos(2kπ/(N+1)) − 1).
    α_k/γ₀ are ALL rational iff N+1 ∈ {1,2,3,4,6} (i.e. N ∈ {0,1,2,3,5}) — the crystallographic set, since
    cos(2π/m) ∈ ℚ ⟺ m ∈ {1,2,3,4,6}. So N=3 is the LAST rational before the first irrational at N=4
    (golden family, the rates carry √5); N=5 is a rational ISLAND (m=6), not a return for good.
    [This is the documented F65 "Niven rationality" paragraph, ANALYTICAL_FORMULAS.md.]

  IM-FACE (the band edge / coherence-hand frequency, F2b / TopologyBandEdge: ω/J = 2·cos(π/(N+1))).
    2cos(π/m) is RATIONAL iff m ≤ 3 (N ≤ 2; last rational N=2 = 1), a single QUADRATIC SURD (degree 2,
    an a±√b form) iff m ∈ {4,5,6} (N=3 √2, N=4 φ, N=5 √3), and degree ≥ 3 (first a CUBIC, no a±√b) from
    m=7 (N=6). The exact degree is [Q(2cos(π/m)):Q] = φ_euler(2m)/2. [NOT yet documented anywhere.]

  V-FACE (the V-Effect gain, docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md: V(N) = 1 + cos(π/N)).
    Niven-rational iff N ∈ {2,3}; for N ≥ 4 it lands on the named constant of the NEXT ring polygon —
    SILVER (√2) at N=4 (square), GOLDEN at N=5 (pentagon, cos π/5), √3 at N=6 (hexagon). Its golden is
    at N=5, NOT N=4, because the angle is π/N, not π/(N+1). [Already documented; referenced here.]

The three faces use the SAME Niven theorem on cyclotomic angles, differing only by the angle convention
(2kπ/m for the rate, π/m for the band edge, π/N for the V-Effect), so they give DIFFERENT cutoffs and
DIFFERENT golden-N. "First golden" is angle-convention-dependent: N=4 on the two SE faces (rate + band
edge, angle π/(N+1)), N=5 on the V-Effect face (angle π/N). The real content is that the golden ratio is
forced by the cyclotomic geometry, not chosen.

This is the number-theoretic root of the SE band-edge / dissipator-rate / V-Effect family of small-N
specials. It is NOT the single root of every special: the (n,n)/{0,2} filling-maximality modes (the
n3_special_cases Re-side entries) are a SEPARATE combinatorial root. Two roots, both real.

Gates PROVE the rationality/degree verdicts with sympy (exact, not numeric). A firing gate is a finding.
Run: python simulations/niven_rationality_root.py
"""
import sys
from math import cos, pi, gcd
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")   # Windows console defaults to cp1252; the output is Unicode


def euler_totient(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def rates_all_rational(m):
    """Are all dissipator rates α_k/γ₀ = −2 sin²(kπ/m), k=1..m−1, rational? (sympy-exact)"""
    for k in range(1, m):
        val = -2 * sp.sin(sp.pi * k / m) ** 2
        if not sp.nsimplify(val).is_rational:
            return False
    return True


def band_edge_degree_formula(m):
    """[Q(2cos(π/m)) : Q] = φ(2m)/2 for m >= 2."""
    return euler_totient(2 * m) // 2


def band_edge_degree_exact(m):
    """sympy ground truth: degree of the minimal polynomial of 2cos(π/m) over Q."""
    x = sp.Symbol("x")
    return int(sp.degree(sp.minimal_polynomial(2 * sp.cos(sp.pi / m), x), x))


SURDS = {sp.sqrt(2): "√2", (1 + sp.sqrt(5)) / 2: "φ = (1+√5)/2", sp.sqrt(3): "√3"}


def band_edge_identity(m):
    val = sp.simplify(sp.nsimplify(2 * sp.cos(sp.pi / m), [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)]))
    if val.is_rational:
        return f"{val} (rational)"
    for surd, label in SURDS.items():
        if sp.simplify(val - surd) == 0:
            return label
    return f"degree {band_edge_degree_exact(m)}"


print("=" * 100)
print("THE NIVEN ROOT: two faces of Niven's theorem on the single-excitation cyclotomic angle π/(N+1)")
print("=" * 100)
print(f"{'N':>2} {'m=N+1':>5}   RE: rates −2sin²(kπ/m)          IM: band edge 2cos(π/m)")
print(f"{'':>2} {'':>5}   rational? (set {{1,2,3,4,6}})       value / algebraic degree")
rows = []
for N in range(1, 9):
    m = N + 1
    re_rat = rates_all_rational(m)
    im_deg = band_edge_degree_exact(m)
    im_id = band_edge_identity(m)
    rows.append((N, m, re_rat, im_deg, im_id))
    print(f"{N:>2} {m:>5}   {'rational' if re_rat else 'IRRATIONAL':<10}                    "
          f"{im_id:<18} (deg {im_deg})")

print("\n" + "=" * 100)
print("GATES (a firing gate is a finding: diagnose, don't loosen)")
print("=" * 100)

re_rational_Ns = [N for (N, m, re_rat, _, _) in rows if re_rat]
deg = {N: im_deg for (N, m, _, im_deg, _) in rows}

# G0: the band-edge degree formula φ(2m)/2 matches sympy's exact minimal-polynomial degree
for (N, m, _, im_deg, _) in rows:
    assert band_edge_degree_formula(m) == im_deg, \
        f"G0 FIRED: N={N} φ(2m)/2={band_edge_degree_formula(m)} != sympy {im_deg}"
print("G0 PASS: [Q(2cos(π/m)):Q] = φ_euler(2m)/2 matches sympy minimal_polynomial for N=1..8.")

# G1 (RE-FACE, the documented F65 root): rates all rational iff N ∈ {1,2,3,5} (N+1 ∈ {2,3,4,6} of {1,2,3,4,6})
assert re_rational_Ns == [1, 2, 3, 5], f"G1 FIRED: rate-rational at {re_rational_Ns}, expected [1,2,3,5]"
print("G1 PASS (RE-face, F65): the dissipator rates are all rational iff N+1 ∈ {1,2,3,4,6} — N=3 the LAST "
      "rational before the gap, N=5 a rational island (m=6).")

# G2 (the SE hinge): N=4 is the FIRST irrational on BOTH SE faces, and golden on both
assert 4 not in re_rational_Ns, "G2 FIRED: N=4 rates unexpectedly rational"
assert min(N for N in range(1, 9) if N not in re_rational_Ns) == 4, "G2 FIRED: first rate-irrational != N=4"
# the N=4 rates carry √5 (golden family): sin²(π/5) = (5−√5)/8
assert sp.simplify(sp.sin(sp.pi / 5) ** 2 - (5 - sp.sqrt(5)) / 8) == 0, "G2 FIRED: sin²(π/5) != (5−√5)/8"
# the N=4 band edge IS φ exactly
assert sp.simplify(2 * sp.cos(sp.pi / 5) - (1 + sp.sqrt(5)) / 2) == 0, "G2 FIRED: 2cos(π/5) != φ"
print("G2 PASS (the SE hinge): N=4 is the FIRST GOLDEN on the two SE faces — rates carry √5 "
      "(sin²(π/5)=(5−√5)/8) AND the band edge is φ = 2cos(π/5) = (1+√5)/2 exactly. The golden ratio enters "
      "both SE faces at once (angle π/(N+1)).")

# G2b (the V-Effect face, OFF_NIVEN_AS_WAVE_BREAKING.md): V(N)=1+cos(π/N), Niven-rational iff N∈{2,3},
# golden at N=5 (NOT N=4 — angle π/N), silver at N=4, √3 at N=6. Guards the "first golden" overclaim.
def v_gain(N):
    return 1 + sp.cos(sp.pi / N)
v_rational_Ns = [N for N in range(2, 9) if sp.nsimplify(v_gain(N)).is_rational]
assert v_rational_Ns == [2, 3], f"G2b FIRED: V-Effect rational at {v_rational_Ns}, expected [2,3]"
assert sp.simplify(v_gain(4) - (1 + sp.sqrt(2) / 2)) == 0, "G2b FIRED: V(4) != 1+√2/2 (silver)"
assert sp.simplify(v_gain(5) - (5 + sp.sqrt(5)) / 4) == 0, "G2b FIRED: V(5) != (5+√5)/4 (golden)"
print("G2b PASS (V-Effect face, OFF_NIVEN_AS_WAVE_BREAKING.md): V(N)=1+cos(π/N) is Niven-rational iff "
      "N∈{2,3}; its golden is at N=5 (pentagon), silver at N=4, √3 at N=6 — golden SHIFTED off N=4 because "
      "the angle is π/N, not π/(N+1). 'First golden' is angle-convention-dependent.")

# G3 (IM-FACE, new): band edge rational iff N ≤ 2 (Niven), single quadratic surd iff N ≤ 5, first cubic N=6
im_rational_Ns = [N for N in deg if deg[N] == 1]
im_quadratic_Ns = sorted(N for N in deg if deg[N] <= 2)
assert im_rational_Ns == [1, 2], f"G3 FIRED: band edge rational at {im_rational_Ns}, expected [1,2]"
assert im_quadratic_Ns == [1, 2, 3, 4, 5], f"G3 FIRED: band edge deg<=2 at {im_quadratic_Ns}, expected [1..5]"
assert deg[6] == 3, f"G3 FIRED: N=6 band-edge degree {deg[6]} != 3 (first cubic)"
print("G3 PASS (IM-face, new): the band edge is rational only at N≤2 (last N=2 = 1), a single quadratic surd "
      "(a±√b) for N=1..5 (√2, φ, √3 at N=3,4,5), and degree ≥3 (first CUBIC) from N=6.")

# G4: the two faces have DIFFERENT cutoffs (guard the classic conflation) but the SAME first-golden hinge
assert set(re_rational_Ns) != set(im_quadratic_Ns), "G4 FIRED: the two faces' N-sets coincide (they must not)"
print("G4 PASS: the RE-face cutoff {1,2,3,4,6} (cos(2π/m), crystallographic) and the IM-face cutoff N+1≤6 "
      "(cos(π/m), quadratic surd) are DISTINCT — same Niven theorem, different double-angle — meeting only at "
      "the N=4 first-golden hinge.")

print("\nROOT: the small-N specials of the band edge / dissipator rates / V-Effect share ONE number-theoretic")
print("root — Niven's theorem on cyclotomic angles. Three faces (angles π/(N+1) for the SE band edge & rates,")
print("π/N for the V-Effect); golden enters at N=4 (SE) or N=5 (V-Effect), convention-dependent. A SEPARATE")
print("combinatorial root (filling-maximality, the (n,n)/{0,2} modes) explains the other n3 specials. DONE.")
