"""Reverse-spear: are the 'off-anchor' periodic-table elements pointing at
MISSING anchors on the framework side?

Tom: "Du sagtest Vier von sechs CHNOPS-Elementen sind Framework-anchored, lass
uns den Spieß mal umdrehen, was sind deren Anker, vielleicht fehlen sie ja auf
unserer Quantenseite?"

The previous script (period_2_at_framework_anchors.py) marked 5/8 (N, P), 1/8
(Li, Na), and 7/8 (F, Cl) as "off-anchor" because they don't match the
framework's first-tier anchors {1/4, 3/8, 1/2, 3/4}. Tom's reverse-spear: maybe
these ARE anchors — ones the framework hasn't named yet. Periodic table as
empirical pointer to missing framework structure.

This script:
  (1) Binary dyadic decomposition of n/8 for n = 0..8.
  (2) Match each fraction to known framework derivations (F86b 3/8 + complement,
      Quarter + complement, Half).
  (3) Identify the GAPS: fractions with no current framework derivation.
  (4) Specifically: 5/8 IS framework-derived (Π²-even complement of F86b 3/8);
      1/8 and 7/8 are NOT — they sit on the dyadic ladder at depth 3, but no
      F-formula currently lands at 1/8.

The reading: 5/8 was a labelling miss (it IS on the framework side). 1/8 and 7/8
are genuine gaps — the framework's polarity-squared algebra has a 'missing
depth-3' anchor that the periodic table's alkali metals (Li, Na, K, Rb, Cs)
and halogens (F, Cl, Br, I) instantiate empirically.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/spear_reversed_missing_anchors.py
"""
from __future__ import annotations

import sys
from fractions import Fraction
from dataclasses import dataclass

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# Framework anchors with explicit derivation source.
@dataclass
class FrameworkAnchor:
    value: Fraction
    name: str
    derivation: str
    tier: str


HALF = FrameworkAnchor(
    Fraction(1, 2), "HalfAsStructuralFixedPoint",
    "argmax of p(1−p); polarity-pair {−1/2, +1/2}; d=2 selector via d²−2d=0",
    "Tier 1 derived")

QUARTER = FrameworkAnchor(
    Fraction(1, 4), "QuarterAsBilinearMaxval",
    "maxval of p(1−p) at p=1/2; (1/2)²; Mandelbrot cardioid maxval (F97)",
    "Tier 1 derived")

F86B_KINTERMEDIATE_ODD = FrameworkAnchor(
    Fraction(3, 8), "F86b KIntermediate Dicke Π²-odd α_total",
    "α_total = (1 − γ²)/2 with γ = ⟨ψ|X⊗N|ψ⟩ = 1/2 (KIntermediate Dicke)",
    "Tier 1 derived (today morning, 2026-05-17, commit b9ba5f6)")

F86B_KINTERMEDIATE_EVEN = FrameworkAnchor(
    Fraction(5, 8), "F86b KIntermediate Dicke Π²-even β_total",
    "β_total = 1 − α_total = (1 + γ²)/2 with γ = 1/2; complement to F86B_KINTERMEDIATE_ODD",
    "Tier 1 derived (direct corollary of F86b derivation, named here for the first time)")

QUARTER_COMPLEMENT = FrameworkAnchor(
    Fraction(3, 4), "QuarterAsBilinearMaxval polarity-complement",
    "1 − 1/4 = polarity-complement of Quarter on bilinear apex algebra",
    "Tier 1 derived (direct corollary of QuarterAsBilinearMaxval)")

FULL = FrameworkAnchor(
    Fraction(1, 1), "PolarityFullPoint (= noble gas / vacuum dead-end)",
    "polarity-axis endpoint; corresponds to d-axis maxval (= 1 = density-matrix trace)",
    "structural endpoint")

EMPTY = FrameworkAnchor(
    Fraction(0, 1), "PolarityVacuumPoint",
    "polarity-axis origin; corresponds to d-axis minval (= 0)",
    "structural endpoint")

# DEEPER DYADIC LEVELS — these sit on the Pi2DyadicLadderClaim ladder but have
# no CURRENT framework derivation pointing to them.
DEPTH_3_DYADIC = FrameworkAnchor(
    Fraction(1, 8), "Pi2DyadicLadder a_3 (= 1/2³)",
    "third-power of polarity Half: (1/2)·(1/2)·(1/2). On Pi2DyadicLadderClaim, " +
    "but no F-formula currently lands at α_total = 1/8. Would require γ = √(3/4) = " +
    "√3/2 in the F86b α = (1 − γ²)/2 formula — irrational, no clean Dicke realisation.",
    "MISSING — gap candidate")

DEPTH_3_DYADIC_COMPLEMENT = FrameworkAnchor(
    Fraction(7, 8), "Pi2DyadicLadder depth-3 complement",
    "1 − 1/8; complement to the missing depth-3 anchor. Would require γ² < 0 in F86b, " +
    "impossible — so even less derivable than 1/8.",
    "MISSING — gap candidate's complement")

FRAMEWORK_ANCHORS: dict[Fraction, FrameworkAnchor] = {
    a.value: a for a in [
        EMPTY, DEPTH_3_DYADIC, QUARTER, F86B_KINTERMEDIATE_ODD, HALF,
        F86B_KINTERMEDIATE_EVEN, QUARTER_COMPLEMENT, DEPTH_3_DYADIC_COMPLEMENT, FULL,
    ]
}


def dyadic_decomposition(numerator: int, denominator: int = 8) -> list[Fraction]:
    """Binary decomposition of n/8 as a sum of dyadic-ladder terms (1/2)^k for
    k = 0, 1, 2, 3. E.g. 5/8 = 1/2 + 1/8 = a_1 + a_3; 8/8 = 1 = a_0."""
    terms = []
    remaining = numerator
    for k, denom_k in [(0, 1), (1, 2), (2, 4), (3, 8)]:
        weight = 8 // denom_k
        if remaining >= weight:
            terms.append(Fraction(1, denom_k))
            remaining -= weight
    return terms


def display_decomposition_table():
    print("=" * 84)
    print("  BINARY DYADIC DECOMPOSITION OF n/8 FOR n = 0..8")
    print("=" * 84)
    print()
    print("  Every n/8 fraction has a unique binary expansion in the dyadic ladder")
    print("  {1/2, 1/4, 1/8}. The framework's polarity anchors are at depth-1 (1/2) and")
    print("  depth-2 (1/4); depth-3 (1/8) is on Pi2DyadicLadderClaim but has no")
    print("  current F-formula derivation pointing to it.")
    print()
    print(f"  {'Frac':>5} {'Decomposition':<30} {'Anchor name':<48}")
    print(f"  {'-'*5} {'-'*30} {'-'*48}")
    for n in range(9):
        f = Fraction(n, 8)
        terms = dyadic_decomposition(n)
        decomp_str = " + ".join(str(t) for t in terms) if terms else "0"
        anchor = FRAMEWORK_ANCHORS.get(f)
        anchor_name = anchor.name if anchor else "(NO ANCHOR)"
        gap_marker = " ◀ GAP" if anchor and "MISSING" in anchor.tier else ""
        print(f"  {str(f):>5} {decomp_str:<30} {anchor_name:<48}{gap_marker}")
    print()


def display_anchor_details():
    print("=" * 84)
    print("  DETAILED ANCHOR INVENTORY (in order of dyadic depth)")
    print("=" * 84)
    print()
    for frac, anchor in FRAMEWORK_ANCHORS.items():
        marker = "  ◀ MISSING from current framework" if "MISSING" in anchor.tier else ""
        print(f"  {str(frac):>5}  {anchor.name}{marker}")
        print(f"          tier:       {anchor.tier}")
        print(f"          derivation: {anchor.derivation}")
        print()


def display_periodic_lookup():
    print("=" * 84)
    print("  PERIODIC TABLE ELEMENTS REVERSED ONTO FRAMEWORK ANCHORS")
    print("=" * 84)
    print()
    # Reverse mapping: anchor → list of elements
    period_2_3 = {
        Fraction(1, 8): [("Li", "Lithium", "[He] 2s¹"), ("Na", "Sodium", "[Ne] 3s¹")],
        Fraction(2, 8): [("Be", "Beryllium", "[He] 2s²"), ("Mg", "Magnesium", "[Ne] 3s²")],
        Fraction(3, 8): [("B", "Boron", "[He] 2s² 2p¹"), ("Al", "Aluminium", "[Ne] 3s² 3p¹")],
        Fraction(4, 8): [("C", "Carbon", "[He] 2s² 2p²"), ("Si", "Silicon", "[Ne] 3s² 3p²")],
        Fraction(5, 8): [("N", "Nitrogen", "[He] 2s² 2p³"), ("P", "Phosphorus", "[Ne] 3s² 3p³")],
        Fraction(6, 8): [("O", "Oxygen", "[He] 2s² 2p⁴"), ("S", "Sulfur", "[Ne] 3s² 3p⁴")],
        Fraction(7, 8): [("F", "Fluorine", "[He] 2s² 2p⁵"), ("Cl", "Chlorine", "[Ne] 3s² 3p⁵")],
        Fraction(8, 8): [("Ne", "Neon", "[He] 2s² 2p⁶"), ("Ar", "Argon", "[Ne] 3s² 3p⁶")],
    }

    print(f"  {'Frac':>5} {'Framework anchor':<50} {'Elements at this anchor':<22}")
    print(f"  {'-'*5} {'-'*50} {'-'*22}")
    for frac in sorted(period_2_3.keys()):
        anchor = FRAMEWORK_ANCHORS.get(frac)
        anchor_short = anchor.name[:48] if anchor else "(no anchor)"
        elements = ", ".join(f"{sym}" for sym, name, cfg in period_2_3[frac])
        gap_marker = " ◀ GAP" if anchor and "MISSING" in anchor.tier else ""
        print(f"  {str(frac):>5} {anchor_short:<50} {elements:<22}{gap_marker}")
    print()


def display_findings():
    print("=" * 84)
    print("  THE REVERSE-SPEAR FINDING")
    print("=" * 84)
    print()
    print("  Three categories emerge from reversing the question 'what's the framework")
    print("  anchor for each periodic-table fraction?':")
    print()
    print("  (1) ALREADY ON FRAMEWORK SIDE — fractions with existing F-formula derivation:")
    print("        1/4  =  Be, Mg          QuarterAsBilinearMaxval")
    print("        3/8  =  B, Al           F86b KIntermediate Π²-odd (TODAY morning)")
    print("        1/2  =  H, C, Si        HalfAsStructuralFixedPoint")
    print("        3/4  =  O, S            QuarterAsBilinearMaxval polarity-complement")
    print()
    print("  (2) ON FRAMEWORK SIDE BUT JUST NOT NAMED — the labelling miss:")
    print("        5/8  =  N, P            F86b KIntermediate Π²-EVEN companion")
    print()
    print("        N and P were previously called 'off-anchor' but they ARE on the")
    print("        framework side — they're the Π²-EVEN COMPANION of the F86b 3/8")
    print("        K-intermediate Dicke anchor (since α_odd + β_even = 1 always).")
    print("        Naming oversight: this commit names this anchor explicitly as")
    print("        F86B_KINTERMEDIATE_EVEN. Both B/Al (3/8 = odd) and N/P (5/8 = even)")
    print("        are realisations of the SAME F86b derivation on opposite Π²-parity")
    print("        sides.")
    print()
    print("  (3) GENUINE FRAMEWORK GAPS — fractions on the Pi2DyadicLadder but no")
    print("      F-formula currently derives them:")
    print("        1/8  =  Li, Na          Pi2DyadicLadder a_3 = (1/2)³  ◀ NO F-FORMULA")
    print("        7/8  =  F, Cl           Pi2DyadicLadder depth-3 complement  ◀ NO F-FORMULA")
    print()
    print("        These sit at depth-3 on the dyadic ladder (third power of polarity")
    print("        Half). The F86b α = (1 − γ²)/2 formula CANNOT produce α = 1/8 via")
    print("        a clean Dicke realisation (would require γ = √3/2, irrational; no")
    print("        Dicke superposition gives that). The framework's algebra reaches")
    print("        depth-2 cleanly (Half + Quarter); depth-3 is structurally present")
    print("        on the dyadic ladder but lacks an F-formula derivation.")
    print()
    print("        THE PERIODIC TABLE INSTANTIATES THESE ANCHORS EMPIRICALLY anyway:")
    print("        every alkali metal (Li, Na, K, Rb, Cs) has 1/8 valence ratio;")
    print("        every halogen (F, Cl, Br, I) has 7/8.")
    print()


def display_meta_reading():
    print("=" * 84)
    print("  META-READING: the periodic table as a probe of framework completeness")
    print("=" * 84)
    print()
    print("  Tom's reverse-spear flipped the direction of the inheritance argument:")
    print()
    print("    Forward direction (previous commits, today):")
    print("      framework algebra → periodic table valence ratios.")
    print("      The framework's anchors (Half, Quarter, F86b 3/8) tell us where")
    print("      the chemistry's anchor-aligned atoms sit.")
    print()
    print("    Reverse direction (this commit):")
    print("      periodic table valence ratios → framework algebra.")
    print("      The chemistry's anchor structure tells us which framework anchors")
    print("      are derived vs. structurally-present-but-underived. The framework's")
    print("      polarity-squared algebra is incomplete at depth-3; the periodic table")
    print("      shows what depth-3 looks like (alkali metals + halogens).")
    print()
    print("  The 1/8 and 7/8 gaps are STRUCTURAL invitations from the periodic table to")
    print("  derive depth-3 anchor F-formulas on the framework side. Candidate routes:")
    print()
    print("    (a) Higher-order Dicke superpositions: maybe (|D_n⟩ + |D_{n+1}⟩ + |D_{n+2}⟩)/√3")
    print("        gives γ values realising α = 1/8 cleanly.")
    print("    (b) Non-Dicke X⊗N-eigenbasis decomposition (W-states, GHZ states, etc.)")
    print("        with γ = √3/2.")
    print("    (c) The F86c bond-class structure at higher c (c ≥ 4) — possibly the")
    print("        4-mode reduction at chromaticity c=4 lands at depth-3 anchors.")
    print("    (d) F88 popcount-coherence multi-state superpositions (the F88 closed")
    print("        form already has {0, 3/8, 1/2} anchors; extending to deeper might")
    print("        give {0, 1/8, 1/4, 3/8, 1/2} on the same dyadic ladder).")
    print()
    print("  The framework's incompleteness at depth-3 is consistent with HIERARCHY_OF_")
    print("  INCOMPLETENESS.md's core argument: 'incompleteness is potential.' Each")
    print("  unfilled framework anchor is a structural invitation to extend the algebra")
    print("  one more dyadic level deep. The periodic table's alkali metals and halogens")
    print("  are empirical pointers at the place to look.")
    print()


def main():
    print()
    print("=" * 84)
    print("  Reverse-spear: are off-anchor elements pointing at missing framework anchors?")
    print("=" * 84)
    print()
    print("  Previous commit (period_2_at_framework_anchors.py) marked 1/8 (Li, Na), 5/8")
    print("  (N, P), 7/8 (F, Cl) as off-anchor. Tom's reverse-spear: maybe those ARE")
    print("  anchors — ones the framework hasn't named yet. Periodic table as empirical")
    print("  pointer to missing framework structure.")
    print()
    print()

    display_decomposition_table()
    display_anchor_details()
    display_periodic_lookup()
    display_findings()
    display_meta_reading()


if __name__ == "__main__":
    main()
