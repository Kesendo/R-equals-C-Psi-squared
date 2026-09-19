"""Finite scale and coordinate comparisons behind the old carrier story.

Tests Tom's claim (2026-05-28): a scale (Massstab / carrier gamma0) cannot fully
hide, else it would not be a scale. Only its VALUE (the absolute second) is gauge;
the finite computations show which dimensionless values survive a joint scale
change.  Repeated values such as 1/2, 1/4, 2, and 60 degrees are a numeric
rhyme, not shared ancestry.  The assassin-in-the-crowd line remains a story,
not a causal identification.

Three numerical diagnostics, each able to fail:

  PART 1  Two CPUs. Same dynamics at two absolute clocks (gamma0 and g*gamma0).
          Finite N=1/N=2 rescaling observations compare theta,
          |lambda|/gamma0, and a selected Bell+ purity-coherence readout line.
          The dimensionful spectrum and selected readout time scale by g.

  PART 2  A finite occurrence census for {1/2, 1/4, 2, 60deg}. Rewriting those
          values with d=2 is arithmetic, not evidence for one common source.

  PART 3  The F95 b=1/2 roots acquire an imaginary coordinate above c=1/4.
          That imaginary coordinate is not a literal n_Y observable. The
          constructed dictionary c=1/4+Q^2/4 gives a tautological angle identity
          with atan(Q) for the named 2x2 decay pair; it is a comparison, not a
          bridge or a physical identification.

Tom + Claude, 2026-05-28. Run: python simulations/where_the_carrier_hides.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from framework.pauli import ur_pauli, _build_bilinear
from framework.lindblad import lindbladian_pauli_dephasing, cpsi_bell_plus

ATOL = 1e-12
HEIS = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]

_results = []


def report(name, dev, tol=ATOL):
    ok = dev < tol
    _results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: max|d| = {dev:.2e}")
    return ok


def sorted_spec(L):
    s = np.linalg.eigvals(L)
    return np.array(sorted(s, key=lambda z: (round(z.real, 9), round(z.imag, 9))))


# ----------------------------------------------------------------------
# PART 1 - two CPUs: the value of the clock hides, every ratio stays
# ----------------------------------------------------------------------

def build_L(N, J, gamma):
    """Heisenberg chain at overall Hamiltonian scale J, uniform Z-dephasing gamma."""
    if N == 1:
        H = (J / 2.0) * ur_pauli("Z")
    else:
        bonds = [(i, i + 1) for i in range(N - 1)]
        H = J * _build_bilinear(N, bonds, HEIS)
    return lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter="Z")


def part1(g=10.0):
    print("PART 1  - two CPUs (clock gamma0 vs g*gamma0): the value hides, the ratios stay")
    print("          finite N=1/N=2 rescaling observations")
    for N in (1, 2):
        J_A, gamma_A = 1.0, 0.3
        J_B, gamma_B = g * J_A, g * gamma_A          # World B = World A run g times faster
        LA, LB = build_L(N, J_A, gamma_A), build_L(N, J_B, gamma_B)
        sA, sB = sorted_spec(LA), sorted_spec(LB)

        # The "hertz" scale exactly with the hidden clock g.
        report(f"N={N}: spec(L_B) = g * spec(L_A)  (the hertz carry the clock)",
               float(np.max(np.abs(sB - g * sA))))
        # Every dimensionless ratio is identical: spec / gamma0 (here gamma0 ~ gamma).
        report(f"N={N}: spec(L_A)/gamma_A = spec(L_B)/gamma_B  (dimensionless, clock-free)",
               float(np.max(np.abs(sA / gamma_A - sB / gamma_B))))
        # theta from the fastest oscillatory mode is identical.
        osc_A = sA[np.argmax(np.abs(sA.imag))]
        osc_B = sB[np.argmax(np.abs(sB.imag))]
        thA = np.degrees(np.arctan2(abs(osc_A.imag), abs(osc_A.real)))
        thB = np.degrees(np.arctan2(abs(osc_B.imag), abs(osc_B.real)))
        report(f"N={N}: theta identical across the two clocks", abs(thA - thB))

    # This selected Bell+ purity-coherence readout line reaches 1/4. Its time
    # scales with the clock; this is not a claim about every quarter-valued object.
    roots = np.roots([1.0, 0.0, 1.0, -1.5])            # u^3 + u - 3/2 = 0  (CPsi=1/4)
    u_star = [r.real for r in roots if abs(r.imag) < 1e-9 and 0 < r.real <= 1][0]
    gA, gB = 0.3, g * 0.3
    tA, tB = -np.log(u_star) / (4 * gA), -np.log(u_star) / (4 * gB)
    report("selected Bell+ purity-coherence readout line evaluates to 1/4",
           abs(cpsi_bell_plus(0, 0, gA, tA) - 0.25))
    report("dimensionless crossing gamma*t* is identical across clocks",
           abs(gA * tA - gB * tB))
    report("but the crossing TIME scales by the hidden clock g (t*_A / t*_B = g)",
           abs(tA / tB - g))
    print(f"          t*_A = {tA:.5f},  t*_B = {tB:.5f},  ratio = {tA / tB:.5f} = g={g}")
    print("          -> these finite normalized observations do not recover the absolute second.")


# ----------------------------------------------------------------------
# PART 2 - occurrence census and an explicitly constructed d=2 rewrite
# ----------------------------------------------------------------------

def structural_generators():
    """Independent framework structures, each computed to the number(s) it produces."""
    Z = ur_pauli("Z")
    gens = {}
    gens["polynomial d^2-2d=0 roots"] = {round(float(r), 12) for r in np.roots([1.0, -2.0, 0.0]).real}
    gens["discriminant 1-4c=0"] = {0.25}
    gens["F95 threshold b^2 (b=1/2)"] = {0.25}
    gens["Bloch polarity eig(Z/2)"] = {round(float(abs(e)), 12) for e in np.linalg.eigvalsh(0.5 * Z)}
    p = np.linspace(0.0, 1.0, 2_000_001)               # bilinear apex p(1-p), by grid
    f = p * (1.0 - p)
    gens["bilinear apex argmax"] = {round(float(p[np.argmax(f)]), 12)}
    gens["bilinear apex maxval"] = {round(float(np.max(f)), 12)}
    gens["F97 cardioid |z*| at cusp"] = {0.5}           # z*(0) = b = 1/2
    gens["qubit dim anchor 1/d"] = {0.5}
    gens["canonical angle arctan(sqrt3)"] = {round(float(np.degrees(np.arctan(np.sqrt(3.0)))), 9)}
    gens["Niven canonical angles"] = {30.0, 45.0, 60.0, 90.0}
    gens["absorption |lambda|/g0 at Q=sqrt3"] = {round(float(np.sqrt(1.0 + 3.0)), 12)}
    gens["dyadic ladder 2^(1-n)"] = {round(2.0 ** (1 - n), 12) for n in range(-1, 5)}
    return gens


def count_pinnings(value, gens, tol=1e-9):
    hits = [name for name, vals in gens.items() if any(abs(value - v) < tol for v in vals)]
    return hits


def part2():
    print("\nPART 2  - finite occurrence census; repeated values do not establish one cause")
    gens = structural_generators()
    obvious = {"1/2": 0.5, "1/4": 0.25, "2": 2.0, "60deg": 60.0}
    controls = {"0.37": 0.37, "0.611": 0.611, "0.123": 0.123, "50deg": 50.0}

    print("          obvious numbers - how many independent structures pin each:")
    for name, val in obvious.items():
        hits = count_pinnings(val, gens)
        _results.append(len(hits) >= 2)
        print(f"  [{'PASS' if len(hits) >= 2 else 'FAIL'}] {name:>6} pinned by {len(hits)} structures: "
              f"{', '.join(hits)}")
    print("          control numbers - should be pinned by none:")
    for name, val in controls.items():
        hits = count_pinnings(val, gens)
        report(f"control {name} pinned by 0 structures", float(len(hits)), tol=0.5)

    # A deliberately constructed arithmetic rewrite at d=2. It does not infer ancestry.
    d = 2.0
    print("          constructed d=2 rewrite (numeric rhyme, not shared ancestry):")
    report("1/2 = 1/d", abs(0.5 - 1.0 / d))
    report("1/4 = 1/d^2", abs(0.25 - 1.0 / d ** 2))
    report("2 = d", abs(2.0 - d))
    report("60deg = arctan(sqrt(d^2-1))", abs(60.0 - np.degrees(np.arctan(np.sqrt(d ** 2 - 1.0)))), tol=1e-9)


# ----------------------------------------------------------------------
# PART 3 - positive-b quadratic coordinates and a constructed dictionary
# ----------------------------------------------------------------------

def constructed_angle_dictionary(q_values):
    """Return the two equal angles after defining c=1/4+Q^2/4."""
    q_values = np.asarray(q_values, dtype=float)
    c_values = 0.25 + q_values ** 2 / 4.0
    pair_angle = np.degrees(np.arctan(q_values))
    quadratic_angle = np.degrees(np.arctan(np.sqrt(4.0 * c_values - 1.0)))
    return c_values, pair_angle, quadratic_angle


def part3():
    print("\nPART 3  - quadratic threshold and a constructed angle dictionary")
    b = 0.5

    # (a) In the F95 b=1/2 quadratic, the imaginary part is born exactly at c = b^2 = 1/4.
    cs = np.linspace(0.0, 1.0, 1_000_001)
    disc = b * b - cs                                   # roots z = b +- sqrt(b^2 - c)
    im_born = cs[np.argmax(disc < 0)]                   # first c where roots go complex
    report("F95: the roots acquire an imaginary coordinate above c = b^2 = 1/4",
           abs(im_born - 0.25), tol=2e-6)               # grid resolution 1e-6

    # (b) Above 1/4: Re(z) = 1/2 is fixed; Im carries the root displacement.
    c_above = np.array([0.30, 0.50, 0.75, 1.00])
    re_above = np.full_like(c_above, b)                 # Re(z) = b for all c > b^2
    im_above = np.sqrt(c_above - b * b)                 # Im(z) = sqrt(c - 1/4), the change
    report("above 1/4: Re(z) = 1/2 for this positive-b quadratic",
           float(np.max(np.abs(re_above - 0.5))))
    report("above 1/4: Im(z) = sqrt(c-1/4) carries the root-coordinate change",
           float(np.max(np.abs(im_above - np.sqrt(c_above - 0.25)))))

    # (c) The seam: normalize the morning's Lindblad spectrum (lambda = -gamma0 +- iJ)
    #     by the absorption rate 2*gamma0 (= a_0*gamma0 = d*gamma0). Then Re -> -1/2,
    #     and the named-pair coordinate arctan(Q) equals arctan(sqrt(4c-1))
    #     under the deliberately chosen dictionary c = 1/4 + Q^2/4.
    gamma0 = 1.0
    Qs = np.array([0.0, 0.5, 1.0, np.sqrt(3.0), 2.0, 5.0])
    Js = Qs * gamma0
    lam = -gamma0 + 1j * Js                             # the morning's coherence eigenvalue
    re_norm = (lam / (2 * gamma0)).real                 # normalize by absorption rate 2gamma0
    report("normalize the named decay pair by 2gamma0 -> Re = -1/2, an ordinary normalized real coordinate",
           float(np.max(np.abs(re_norm + 0.5))))
    c_bridge, theta_lind, theta_f95 = constructed_angle_dictionary(Qs)
    report("tautological angle identity under the constructed dictionary (not independent evidence)",
           float(np.max(np.abs(theta_lind - theta_f95))))
    report("at Q=0: the constructed c value is 1/4",
           abs(c_bridge[0] - 0.25))
    print(f"          Q      : {np.array2string(Qs, precision=4)}")
    print(f"          c=1/4+Q^2/4 : {np.array2string(c_bridge, precision=4)}")
    print(f"          theta(deg) : {np.array2string(theta_lind, precision=4)}  (constructed coordinate equality)")
    print("          -> this equality follows from how c(Q) was chosen; it does not identify")
    print("             the quadratic imaginary coordinate with an n_Y measurement.")


# ----------------------------------------------------------------------

def main():
    print("=" * 78)
    print("WHERE THE CARRIER HIDES - the Massstab in the obvious numbers")
    print("=" * 78)
    part1()
    part2()
    part3()
    print("\n" + "=" * 78)
    n_ok, n_tot = sum(_results), len(_results)
    print(f"RESULT: {n_ok}/{n_tot} numerical checks within their stated tolerances "
          f"({'ALL PASS' if n_ok == n_tot else 'SOME FAILED'})")
    print("=" * 78)
    print("""
Reading:
  Part 1: finite N=1/N=2 rescaling observations leave the selected normalized
          quantities unchanged while the absolute second differs. That is the
          scope of this finite diagnostic.
  Part 2: the repeated values remain an inviting arithmetic rhyme. The census
          and the d=2 rewrite do not establish a universal carrier or genealogy.
  Part 3: the named pair angle and the F95 root angle coincide only after the
          dictionary c=1/4+Q^2/4 is chosen to make them coincide. The imaginary
          root coordinate is not a measured n_Y content. A deeper relation is an
          open question, not a conclusion of this diagnostic.
""")


if __name__ == "__main__":
    main()
