"""THE NIVEN ROOT (gate-first): the number-theoretic root of the small-N specials, with THREE faces.

The single-excitation chain is cyclotomic. For the uniform open XX chain with one dephased endpoint, F65's
H-eigenmode rate comb has exact first-order-in-γ₀/J coefficients a_k = (4/(N+1))·sin²(kπ/(N+1)), with
α_k^full = γ₀·a_k + O(γ₀³/J²) and α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality belongs to the first-order coefficient comb. At finite γ₀/J the
relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so
no exact finite-γ₀/J full-L rationality is claimed. Niven's theorem
(the only rational cosines of rational-π angles are 0, ±1/2, ±1) puts a number-theoretic ceiling on the
comb's exact coefficient arithmetic, and THAT is the shared root of "why small N is special":

  RE-FACE (F65's first-order endpoint comb: a_k = (4/(N+1))·sin²(kπ/(N+1))).
    These coefficients are ALL rational iff N+1 ∈ {1,2,3,4,6} (i.e. N ∈ {0,1,2,3,5}) — the crystallographic set, since
    cos(2π/m) ∈ ℚ ⟺ m ∈ {1,2,3,4,6}. So N=3 is the LAST rational before the first irrational at N=4
    (golden family, the first-order comb carries √5); N=5 is a rational ISLAND (m=6), not a return for good.
    [This is the documented F65 "Niven rationality" paragraph, ANALYTICAL_FORMULAS.md.]

  IM-FACE (the band edge / coherence-hand frequency, F2b / TopologyBandEdge: ω/J = 2·cos(π/(N+1))).
    2cos(π/m) is RATIONAL iff m ≤ 3 (N ≤ 2; last rational N=2 = 1), a single QUADRATIC SURD (degree 2,
    an a±√b form) iff m ∈ {4,5,6} (N=3 √2, N=4 φ, N=5 √3), and degree ≥ 3 (first a CUBIC, no a±√b) from
    m=7 (N=6). The exact degree is [Q(2cos(π/m)):Q] = φ_euler(2m)/2.
    [Documented in docs/ANALYTICAL_FORMULAS.md and the typed Niven claim/witness.]

  F6-FACE (the F6 Q-edge gain, historical alias "V-Effect gain",
    docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md: V(N) = 1 + cos(π/N)).
    Niven-rational iff N ∈ {2,3}; for N ≥ 4 it lands on the named constant of the NEXT ring polygon —
    SILVER (√2) at N=4 (square), GOLDEN at N=5 (pentagon, cos π/5), √3 at N=6 (hexagon). Its golden is
    at N=5, NOT N=4, because the angle is π/N, not π/(N+1). [Already documented; referenced here.]

The three faces use the SAME Niven theorem on cyclotomic angles, differing only by the angle convention
(2kπ/m for the first-order comb, π/m for the band edge, π/N for the F6 Q-edge gain), so they give DIFFERENT cutoffs and
DIFFERENT golden-N. "First golden" is angle-convention-dependent: N=4 on the two SE faces (comb + band
edge, angle π/(N+1)), N=5 on the F6 face (angle π/N). The real content is that the golden ratio is
forced by the cyclotomic geometry, not chosen.

This is the number-theoretic root of the SE band-edge / first-order endpoint-rate-comb / F6-Q-edge-gain family of small-N
specials. It is NOT the single root of every special: the (n,n)/{0,2} filling-maximality modes (the
n3_special_cases Re-side entries) are a SEPARATE combinatorial root. Two roots, both real.

The all-N IM cutoff comes from the complete low-totient classification φ(r)≤4; G0 is a bounded SymPy
cross-check through N=10, not the all-N proof. All arithmetic is exact. A firing gate is a finding.
Run: python simulations/niven_rationality_root.py
"""
import sys
from math import gcd
import sympy as sp


def euler_totient(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def f65_first_order_coefficient(n, k):
    """Return the exact F65 first-order endpoint-comb coefficient a_k."""
    if n < 1:
        raise ValueError("F65 rate ratio requires n >= 1")
    if not 1 <= k <= n:
        raise ValueError("F65 rate ratio requires 1 <= k <= n")
    m = n + 1
    return sp.simplify(sp.Rational(4, m) * sp.sin(sp.pi * k / m) ** 2)


def first_order_comb_all_rational(m):
    """Are all F65 first-order coefficients (4/m) sin²(kπ/m), k=1..m−1, rational?"""
    for k in range(1, m):
        val = f65_first_order_coefficient(m - 1, k)
        if not sp.nsimplify(val).is_rational:
            return False
    return True


def band_edge_degree_formula(m):
    """[Q(2cos(π/m)) : Q] = φ(2m)/2 for m >= 2."""
    return euler_totient(2 * m) // 2


def low_totient_moduli_at_most_four():
    """Complete classification of r with φ(r) <= 4, independent of SymPy.

    If a prime p divides r, then p-1 divides φ(r).  Thus φ(r) <= 4 permits
    only p in {2,3,5}.  The same divisibility applied to the prime powers
    bounds their exponents by 2^3, 3^1 and 5^1.  The finite candidate set
    below is therefore exhaustive, rather than a cutoff sweep.
    """
    candidates = {
        (2 ** a) * (3 ** b) * (5 ** c)
        for a in range(4)
        for b in range(2)
        for c in range(2)
    }
    return tuple(
        r for r in sorted(candidates) if euler_totient(r) <= 4
    )


def band_edge_low_degree_partition():
    """All-N degree-1/degree-2 partition from d(N)=φ(2(N+1))/2."""
    rational_n = []
    quadratic_n = []
    for twice_m in low_totient_moduli_at_most_four():
        if twice_m < 4 or twice_m % 2:
            continue
        n = twice_m // 2 - 1
        degree = euler_totient(twice_m) // 2
        if degree == 1:
            rational_n.append(n)
        elif degree == 2:
            quadratic_n.append(n)
    low_degree_n = rational_n + quadratic_n
    return {
        "rational_n": tuple(rational_n),
        "quadratic_n": tuple(quadratic_n),
        "higher_degree_from_n": max(low_degree_n) + 1,
    }


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


def v_gain(N):
    return 1 + sp.cos(sp.pi / N)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    print("=" * 100)
    print("THE NIVEN ROOT: three faces of Niven's theorem on cyclotomic angles")
    print("=" * 100)
    print(f"{'N':>2} {'m=N+1':>5}   RE: first-order comb (4/m)sin²(kπ/m) IM: band edge 2cos(π/m)")
    print(f"{'':>2} {'':>5}   all rational? (set {{1,2,3,4,6}})   value / algebraic degree")
    rows = []
    for n in range(1, 11):
        m = n + 1
        re_rat = first_order_comb_all_rational(m)
        im_deg = band_edge_degree_exact(m)
        im_id = band_edge_identity(m)
        rows.append((n, m, re_rat, im_deg, im_id))
        print(f"{n:>2} {m:>5}   {'all rational' if re_rat else 'contains irrational coefficients':<31}  "
              f"{im_id:<18} (deg {im_deg})")

    print("\n" + "=" * 100)
    print("GATES (a firing gate is a finding: diagnose, don't loosen)")
    print("=" * 100)

    re_rational_ns = [n for (n, _, re_rat, _, _) in rows if re_rat]
    deg = {n: im_deg for (n, _, _, im_deg, _) in rows}

    # G0: bounded SymPy cross-check of the exact degree formula, including N=9 beyond the old N<=8 sweep.
    for n, m, _, im_deg, _ in rows:
        assert band_edge_degree_formula(m) == im_deg, \
            f"G0 FIRED: N={n} φ(2m)/2={band_edge_degree_formula(m)} != sympy {im_deg}"
    print("G0 PASS (bounded cross-check): [Q(2cos(π/m)):Q] = φ_euler(2m)/2 matches "
          "SymPy minimal_polynomial for N=1..10; the all-N cutoff is proved separately in G3.")

    # G1: the documented F65 first-order endpoint comb is rational exactly on the crystallographic set.
    assert re_rational_ns == [1, 2, 3, 5], \
        f"G1 FIRED: first-order-comb rational at {re_rational_ns}, expected [1,2,3,5]"
    print("G1 PASS (RE-face, F65): for the uniform open XX chain with one dephased endpoint, the first-order "
          "H-eigenmode rate-comb coefficients are all rational iff N+1 ∈ {1,2,3,4,6} — N=3 the LAST rational "
          "before the gap, N=5 a rational island (m=6); a_k = (4/(N+1))·sin²(kπ/(N+1)), α_k^full = γ₀·a_k + O(γ₀³/J²), "
          "and α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality belongs to the first-order coefficient comb. At finite γ₀/J the "
          "relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so "
          "no exact finite-γ₀/J full-L rationality is claimed.")

    # G1a: feed the old transcription through the same N=3, k=1 angle.
    exact_n3_k1 = f65_first_order_coefficient(3, 1)
    assert exact_n3_k1 == sp.Rational(1, 2), \
        f"G1a FIRED: F65 N=3,k=1 = {exact_n3_k1}, expected 1/2"
    old_misnormalization_n3_k1 = sp.simplify(-2 * sp.sin(sp.pi / 4) ** 2)
    assert old_misnormalization_n3_k1 != sp.Rational(1, 2), \
        "G1a CONTROL FIRED: the old -2·sin² transcription unexpectedly passed the positive F65 pin"
    print("G1a PASS (F65 normalization): a_1(3) = 1/2; the old -2·sin² transcription gives -1 and fails.")

    # G2: N=4 is the first irrational on both SE faces, and golden on both.
    assert 4 not in re_rational_ns, "G2 FIRED: N=4 rates unexpectedly rational"
    assert min(n for n in range(1, 9) if n not in re_rational_ns) == 4, \
        "G2 FIRED: first rate-irrational != N=4"
    assert sp.simplify(sp.sin(sp.pi / 5) ** 2 - (5 - sp.sqrt(5)) / 8) == 0, \
        "G2 FIRED: sin²(π/5) != (5−√5)/8"
    assert sp.simplify(2 * sp.cos(sp.pi / 5) - (1 + sp.sqrt(5)) / 2) == 0, \
        "G2 FIRED: 2cos(π/5) != φ"
    print("G2 PASS (the SE hinge): N=4 is the FIRST GOLDEN on the two SE faces — the first-order comb carries √5 "
          "(sin²(π/5)=(5−√5)/8) AND the band edge is φ = 2cos(π/5) = (1+√5)/2 exactly. The golden ratio enters "
          "both SE faces at once (angle π/(N+1)).")

    # G2b: the F6 Q-edge-gain face has its golden value at N=5, not N=4.
    v_rational_ns = [n for n in range(2, 9) if sp.nsimplify(v_gain(n)).is_rational]
    assert v_rational_ns == [2, 3], \
        f"G2b FIRED: F6 Q-edge gain rational at {v_rational_ns}, expected [2,3]"
    assert sp.simplify(v_gain(4) - (1 + sp.sqrt(2) / 2)) == 0, \
        "G2b FIRED: V(4) != 1+√2/2 (silver)"
    assert sp.simplify(v_gain(5) - (5 + sp.sqrt(5)) / 4) == 0, \
        "G2b FIRED: V(5) != (5+√5)/4 (golden)"
    print("G2b PASS (F6 Q-edge-gain face; historical V-Effect-gain alias): V(N)=1+cos(π/N) is Niven-rational iff "
          "N∈{2,3}; its golden is at N=5 (pentagon), silver at N=4, √3 at N=6 — golden SHIFTED off N=4 because "
          "the angle is π/N, not π/(N+1). 'First golden' is angle-convention-dependent.")

    # G3: complete all-N low-degree classification, plus the bounded SymPy rows as independent anchors.
    low_totient_moduli = low_totient_moduli_at_most_four()
    assert low_totient_moduli == (1, 2, 3, 4, 5, 6, 8, 10, 12), \
        f"G3 FIRED: complete φ(r)<=4 classification changed to {low_totient_moduli}"
    all_n_partition = band_edge_low_degree_partition()
    assert all_n_partition == {
        "rational_n": (1, 2),
        "quadratic_n": (3, 4, 5),
        "higher_degree_from_n": 6,
    }, f"G3 FIRED: all-N band-edge degree partition changed to {all_n_partition}"

    im_rational_ns = [n for n in deg if deg[n] == 1]
    im_degree_at_most_two_ns = sorted(n for n in deg if deg[n] <= 2)
    im_quadratic_surd_ns = sorted(n for n in deg if deg[n] == 2)
    assert im_rational_ns == [1, 2], \
        f"G3 FIRED: band edge rational at {im_rational_ns}, expected [1,2]"
    assert deg[2] == 1, f"G3 FIRED: N=2 degree {deg[2]} != 1 (rational anchor)"
    assert im_degree_at_most_two_ns == [1, 2, 3, 4, 5], \
        f"G3 FIRED: band edge deg<=2 at {im_degree_at_most_two_ns}, expected [1..5]"
    assert im_quadratic_surd_ns == [3, 4, 5], \
        f"G3 FIRED: quadratic-surd band edge at {im_quadratic_surd_ns}, expected [3,4,5]"
    assert deg[6] == 3, f"G3 FIRED: N=6 band-edge degree {deg[6]} != 3 (first cubic)"
    print("G3 PASS (IM-face, all N): the complete low-totient classification φ(r)≤4 gives degree at most 2 "
          "iff N≤5; it is rational only at N≤2 "
          "(last N=2 = 1), a single quadratic surd (a±√b) exactly at N=3,4,5 (√2, φ, √3), and "
          "degree ≥3 (first CUBIC) from N=6.")

    # G4: the two faces have distinct cutoffs even though they meet at the N=4 golden hinge.
    assert set(re_rational_ns) != set(im_degree_at_most_two_ns), \
        "G4 FIRED: the two faces' N-sets coincide (they must not)"
    print("G4 PASS: the RE-face cutoff {1,2,3,4,6} (cos(2π/m), crystallographic) and the IM-face cutoff N+1≤6 "
          "(cos(π/m), degree at most 2) are DISTINCT — same Niven theorem, different double-angle — meeting only at "
          "the N=4 first-golden hinge.")

    print("\nROOT: the small-N specials of the band edge / first-order endpoint comb / F6 Q-edge gain share ONE number-theoretic root.")
    print("For the uniform open XX chain with one dephased endpoint, its RE face is F65's first-order H-eigenmode rate comb; Niven's theorem fixes that comb's coefficient arithmetic.")
    print("Canonical RE scope: a_k = (4/(N+1))·sin²(kπ/(N+1)); α_k^full = γ₀·a_k + O(γ₀³/J²); α_k^full/γ₀ = a_k + O((γ₀/J)²); exact Niven rationality belongs to the first-order coefficient comb.")
    print("At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so no exact finite-γ₀/J full-L rationality is claimed.")
    print("Three faces remain: angles π/(N+1) for the SE band edge and comb, π/N for the F6 face; golden enters at N=4 (SE) or N=5 (F6), convention-dependent. A SEPARATE")
    print("combinatorial root (filling-maximality, the (n,n)/{0,2} modes) explains the other n3 specials. DONE.")


if __name__ == "__main__":
    main()
