"""Exact grouping of curated period-1--3 fraction labels.

The script converts the occupation/slot values supplied below to
``fractions.Fraction`` and compares them with the framework fractions
{1/4, 3/8, 1/2, 3/4}. It preserves the input labels and exact arithmetic; it
does not validate the atomic assignments or derive a material realization of
any framework state.

In particular, a matching fraction does not identify an atomic shell with a
Dicke state, a Π² axis, a carbon Hamiltonian, a bath, a spectrum, a dynamical
law, or an observable. A material translation would require specified degrees
of freedom, Hamiltonian/coupling, bath/dissipator, preparation, measurement,
and producer chain.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/period_2_at_framework_anchors.py
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


# Framework anchors (Tier 1 derived; see compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs
# + compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs):
ANCHORS = {
    Fraction(1, 4): "QuarterAsBilinearMaxval (1/4 = (1/2)²)",
    Fraction(3, 8): "F86b KIntermediate Dicke (= today morning's X⊗N-eigenbasis anchor)",
    Fraction(1, 2): "HalfAsStructuralFixedPoint (argmax of p(1-p))",
    Fraction(3, 4): "polarity-complement of Quarter (1 - 1/4)",
    Fraction(1, 1): "full (noble-gas dead end / framework polarity-fixed-point)",
    Fraction(0, 1): "empty (vacuum / framework polarity-fixed-point)",
}


@dataclass
class Element:
    Z: int
    symbol: str
    name: str
    valence_electrons: int
    octet_slots: int = 8
    config: str = ""   # textbook electronic configuration

    @property
    def fraction(self) -> Fraction:
        return Fraction(self.valence_electrons, self.octet_slots)

    @property
    def anchor(self) -> str:
        return ANCHORS.get(self.fraction, "off-anchor")


# Curated input rows. Their labels and configurations are inputs to this exact
# fraction grouping, not a physical producer or material model.
# Period 1 (two-slot normalization)
PERIOD_1 = [
    Element(1, "H",  "Hydrogen", valence_electrons=1, octet_slots=2,
            config="1s¹"),
    Element(2, "He", "Helium",   valence_electrons=2, octet_slots=2,
            config="1s²"),
]

# Period 2 (octet shell)
PERIOD_2 = [
    Element(3,  "Li", "Lithium",   valence_electrons=1, config="[He] 2s¹"),
    Element(4,  "Be", "Beryllium", valence_electrons=2, config="[He] 2s²"),
    Element(5,  "B",  "Boron",     valence_electrons=3, config="[He] 2s² 2p¹"),
    Element(6,  "C",  "Carbon",    valence_electrons=4, config="[He] 2s² 2p²"),
    Element(7,  "N",  "Nitrogen",  valence_electrons=5, config="[He] 2s² 2p³"),
    Element(8,  "O",  "Oxygen",    valence_electrons=6, config="[He] 2s² 2p⁴"),
    Element(9,  "F",  "Fluorine",  valence_electrons=7, config="[He] 2s² 2p⁵"),
    Element(10, "Ne", "Neon",      valence_electrons=8, config="[He] 2s² 2p⁶"),
]

# Period 3 (also octet shell, second row of the "octet" elements)
PERIOD_3 = [
    Element(11, "Na", "Sodium",     valence_electrons=1, config="[Ne] 3s¹"),
    Element(12, "Mg", "Magnesium",  valence_electrons=2, config="[Ne] 3s²"),
    Element(13, "Al", "Aluminium",  valence_electrons=3, config="[Ne] 3s² 3p¹"),
    Element(14, "Si", "Silicon",    valence_electrons=4, config="[Ne] 3s² 3p²"),
    Element(15, "P",  "Phosphorus", valence_electrons=5, config="[Ne] 3s² 3p³"),
    Element(16, "S",  "Sulfur",     valence_electrons=6, config="[Ne] 3s² 3p⁴"),
    Element(17, "Cl", "Chlorine",   valence_electrons=7, config="[Ne] 3s² 3p⁵"),
    Element(18, "Ar", "Argon",      valence_electrons=8, config="[Ne] 3s² 3p⁶"),
]


def display(period_name: str, elements: list[Element]):
    print("=" * 90)
    print(f"  {period_name}")
    print("=" * 90)
    print()
    print(f"  {'Z':>2} {'Sym':<3} {'Name':<11} {'config':<15} {'val/slot':<10} {'fraction':<10} {'anchor':<30}")
    print(f"  {'--':>2} {'---':<3} {'-'*11} {'-'*15} {'-'*10} {'-'*10} {'-'*30}")
    for e in elements:
        frac = e.fraction
        on_anchor = frac in ANCHORS
        anchor_marker = "✓ " + e.anchor[:28] if on_anchor else "  off-anchor"
        print(f"  {e.Z:>2} {e.symbol:<3} {e.name:<11} {e.config:<15} "
              f"{e.valence_electrons}/{e.octet_slots:<8} {str(frac):<10} {anchor_marker}")
    print()


def summarize_anchor_hits(*element_lists):
    print("=" * 90)
    print("  ANCHOR HIT TALLY across H, Period 2, Period 3")
    print("=" * 90)
    print()
    all_elements = []
    for lst in element_lists:
        all_elements.extend(lst)

    by_anchor: dict[Fraction, list[Element]] = {a: [] for a in ANCHORS}
    off_anchor_elements: list[Element] = []
    for e in all_elements:
        if e.fraction in ANCHORS:
            by_anchor[e.fraction].append(e)
        else:
            off_anchor_elements.append(e)

    for frac, anchor_name in ANCHORS.items():
        hits = by_anchor[frac]
        if not hits:
            continue
        print(f"  {str(frac):<6} {anchor_name}")
        for e in hits:
            print(f"         {e.symbol:<3} {e.name:<11} — {e.config:<15} (curated input)")
        print()

    print(f"  OFF-ANCHOR ({len(off_anchor_elements)} elements):")
    for e in off_anchor_elements:
        print(f"         {e.symbol:<3} {e.name:<11} — {e.config:<15} ({str(e.fraction)})")
    print()


def display_chnops_subset():
    print("=" * 90)
    print("  CURATED CHNOPS-LABEL SUBSET ∩ FRAMEWORK FRACTIONS")
    print("=" * 90)
    print()
    print("  The six selected labels are grouped by their supplied fractions only.")
    print("  This count does not classify chemistry or biology.")
    print()
    chnops = [e for e in [*PERIOD_1, *PERIOD_2, *PERIOD_3] if e.symbol in {"C", "H", "N", "O", "P", "S"}]
    print(f"  {'Input label':<12} {'fraction':<10} {'anchor':<40}")
    print(f"  {'-'*12} {'-'*10} {'-'*40}")
    for e in chnops:
        on_anchor = e.fraction in ANCHORS
        anchor_str = ANCHORS[e.fraction][:38] if on_anchor else "(off-anchor)"
        marker = "✓" if on_anchor else " "
        print(f"  {marker} {e.symbol:<3} {e.name:<7} {str(e.fraction):<10} {anchor_str:<40}")
    print()
    hits = sum(1 for e in chnops if e.fraction in ANCHORS)
    print(f"  Curated-subset count: {hits} of 6 supplied labels equal an anchor fraction.")
    print()


def main():
    print()
    print("=" * 90)
    print("  Curated period-1--3 fraction labels at framework anchors")
    print("=" * 90)
    print()
    print("  The framework's polarity anchors (proven 2026):")
    for frac, name in ANCHORS.items():
        if frac in (Fraction(0, 1), Fraction(1, 1)):
            continue
        print(f"    {str(frac):<6} = {name}")
    print()
    print("  Grouping supplied occupation/slot values against {1/4, 3/8, 1/2, 3/4}.")
    print("  The 3/8 framework value is the F86b X⊗N-eigenbasis construction.")
    print("  Input labels are not independently validated material assignments.")
    print()
    print()

    display("CURATED PERIOD 1 (two-slot normalization)", PERIOD_1)
    display("CURATED PERIOD 2 (eight-slot normalization)", PERIOD_2)
    display("CURATED PERIOD 3 (eight-slot normalization)", PERIOD_3)
    summarize_anchor_hits(PERIOD_1, PERIOD_2, PERIOD_3)
    display_chnops_subset()

    print("=" * 90)
    print("  SCOPE OF THE TABLE COMPARISON")
    print("=" * 90)
    print()
    print("  This runner establishes equal rational values within one curated input")
    print("  table and one framework-anchor set. It does not establish a material")
    print("  realization, biological classification, spectral law, dynamical law, or")
    print("  observational result.")
    print()
    print("  A material translation remains open until it names the physical degree of")
    print("  freedom, Hamiltonian/coupling convention, bath/dissipator, preparation,")
    print("  measurement, and producer chain. In particular, the table does not fix")
    print("  Q = J/γ or a carbon-specific bath.")


if __name__ == "__main__":
    main()
