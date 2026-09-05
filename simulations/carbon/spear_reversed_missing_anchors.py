"""Curated period-label lookup beside framework fraction arithmetic.

The script prints binary decompositions of n/8 and named framework relations:
F99 supplies α = 1/8 at θ = 30° in its non-uniform-Dicke construction, and
7/8 is the arithmetic complement. The period-2/3 symbols are curated lookup
labels only. They do not select a material degree of freedom, Hamiltonian,
bath, preparation, observable, or measurement map, so any cross-domain bridge
remains conditional and open.

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
    "Tier 1 derived (direct F86b corollary)")

QUARTER_COMPLEMENT = FrameworkAnchor(
    Fraction(3, 4), "QuarterAsBilinearMaxval polarity-complement",
    "1 − 1/4 = polarity-complement of Quarter on bilinear apex algebra",
    "Tier 1 derived (direct corollary of QuarterAsBilinearMaxval)")

FULL = FrameworkAnchor(
    Fraction(1, 1), "PolarityFullPoint",
    "polarity-axis endpoint; corresponds to d-axis maxval (= 1 = density-matrix trace)",
    "structural endpoint")

EMPTY = FrameworkAnchor(
    Fraction(0, 1), "PolarityVacuumPoint",
    "polarity-axis origin; corresponds to d-axis minval (= 0)",
    "structural endpoint")

DEPTH_3_DYADIC = FrameworkAnchor(
    Fraction(1, 8), "F99 canonical 30° anchor",
    "α = sin²(30°)/2 = 1/8; γ = √3/2 in F99's non-uniform-Dicke construction.",
    "Tier 1 derived")

DEPTH_3_DYADIC_COMPLEMENT = FrameworkAnchor(
    Fraction(7, 8), "F99 1/8 arithmetic complement",
    "1 − 1/8 = 7/8; complement arithmetic, not an F99 α value.",
    "formal complement")

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
    print("  Every n/8 fraction has a unique binary expansion in {1/2, 1/4, 1/8}.")
    print("  F99 supplies the 1/8 row; 7/8 is its arithmetic complement.")
    print()
    print(f"  {'Frac':>5} {'Decomposition':<30} {'Anchor name':<48}")
    print(f"  {'-'*5} {'-'*30} {'-'*48}")
    for n in range(9):
        f = Fraction(n, 8)
        terms = dyadic_decomposition(n)
        decomp_str = " + ".join(str(t) for t in terms) if terms else "0"
        anchor = FRAMEWORK_ANCHORS.get(f)
        anchor_name = anchor.name if anchor else "(no named relation)"
        print(f"  {str(f):>5} {decomp_str:<30} {anchor_name:<48}")
    print()


def display_anchor_details():
    print("=" * 84)
    print("  DETAILED ANCHOR INVENTORY (in order of dyadic depth)")
    print("=" * 84)
    print()
    for frac, anchor in FRAMEWORK_ANCHORS.items():
        print(f"  {str(frac):>5}  {anchor.name}")
        print(f"          tier:       {anchor.tier}")
        print(f"          derivation: {anchor.derivation}")
        print()


def display_curated_lookup():
    print("=" * 84)
    print("  CURATED PERIOD-2/3 LABEL LOOKUP BESIDE FRAMEWORK FRACTIONS")
    print("=" * 84)
    print()
    # Curated labels only; this table makes no material assignment.
    period_2_3 = {
        Fraction(1, 8): ["Li", "Na"], Fraction(2, 8): ["Be", "Mg"],
        Fraction(3, 8): ["B", "Al"], Fraction(4, 8): ["C", "Si"],
        Fraction(5, 8): ["N", "P"], Fraction(6, 8): ["O", "S"],
        Fraction(7, 8): ["F", "Cl"], Fraction(8, 8): ["Ne", "Ar"],
    }

    print(f"  {'Frac':>5} {'Framework relation':<50} {'Curated labels':<22}")
    print(f"  {'-'*5} {'-'*50} {'-'*22}")
    for frac in sorted(period_2_3.keys()):
        anchor = FRAMEWORK_ANCHORS.get(frac)
        anchor_short = anchor.name[:48] if anchor else "(no anchor)"
        labels = ", ".join(period_2_3[frac])
        print(f"  {str(frac):>5} {anchor_short:<50} {labels:<22}")
    print()


def display_findings():
    print("=" * 84)
    print("  FRACTION-LOOKUP SCOPE")
    print("=" * 84)
    print()
    print("  F99 supplies α = 1/8 at θ = 30° for its specified non-uniform-Dicke")
    print("  construction. The 7/8 row is the arithmetic complement 1 − 1/8.")
    print("  The 3/8 and 5/8 rows are the F86b odd/even complement pair; Quarter and")
    print("  Half retain their named framework relations.")
    print()
    print("  Curated period labels occur beside these numbers in this lookup only.")
    print("  They do not identify a physical system or establish a chemical, material,")
    print("  or causal bridge.")
    print()


def display_open_bridge():
    print("=" * 84)
    print("  OPEN CONDITIONAL BRIDGE")
    print("=" * 84)
    print()
    print("  A material bridge would require a selected physical degree of freedom,")
    print("  Hamiltonian, bath channel and rate, preparation, observable, measurement,")
    print("  and producer. This lookup supplies none of them.")
    print()
    print("  Other state families or channels may be studied only after their own")
    print("  α relation and stated model are supplied.")
    print()


def main():
    print()
    print("=" * 84)
    print("  Curated period-label lookup beside framework fraction arithmetic")
    print("=" * 84)
    print()
    print("  F99 supplies α = 1/8 at θ = 30° in its non-uniform-Dicke construction;")
    print("  7/8 is the arithmetic complement. Curated labels are not a material map.")
    print()
    print()

    display_decomposition_table()
    display_anchor_details()
    display_curated_lookup()
    display_findings()
    display_open_bridge()


if __name__ == "__main__":
    main()
