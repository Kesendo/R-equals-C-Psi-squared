using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>SymmetryFamilyInventory (Tier 1 derived; 2026-05-12): aggregator-Claim
/// for the discrete symmetries of the chain XY+Z-dephasing Liouvillian, classified by
/// axis (operator, antilinear, parameter) and effect (splitting, pairing, spectral
/// reflection, criterion). The typed members wired through this Claim are the five
/// splitting axes, the two sector-pairing primitives and the one negative result;
/// the page <c>docs/SYMMETRY_FAMILY_INVENTORY.md</c> lists the whole family, the later
/// members (the mirror group D4, the antilinear triangle, the Klein group on the
/// self-paired sector, the gauge, the fold lattice, the order sorting, the two
/// criteria) anchored to their own Claims.</summary>
public sealed class SymmetryFamilyInventory : Claim
{
    public SymmetryFamilyInventory()
        : base("SymmetryFamilyInventory: typed family of discrete symmetries of chain XY+Z-deph L, classified by axis (operator / antilinear / parameter) and effect (splitting / pairing / spectral reflection / criterion); this Claim wires 5 axes + 2 sector-pairing primitives + 1 negative result, the page docs/SYMMETRY_FAMILY_INVENTORY.md lists the whole family.",
               Tier.Tier1Derived,
               "docs/SYMMETRY_FAMILY_INVENTORY.md + reflections/ON_THE_SYMMETRY_FAMILY.md")
    { }

    public override string DisplayName =>
        "SymmetryFamilyInventory: typed family of discrete L-symmetries on operator, antilinear and parameter axes";

    public override string Summary =>
        $"family inventory aggregating the substrate U(1)×U(1), the F71 splitting, the three parameter-side invariances F91 + F92 + F93, the two sector-pairing primitives Z⊗N + X⊗N, and the F71_col×F71_row negative result; the whole family (Π = R·D and D₄, the antilinear triangle, K with the transpose, the dead-set antiunitaries, the Klein group on the self-paired sector, the fold lattice, the bipartite gauge, the order sorting, the gamma fold, the frozen divisor, the two-end count, the complement connection) is listed on docs/SYMMETRY_FAMILY_INVENTORY.md ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("anchor markdown",
                summary: "docs/SYMMETRY_FAMILY_INVENTORY.md (the tables: operator, antilinear, parameter, criteria)");
            yield return new InspectableNode("anchor reflection",
                summary: "reflections/ON_THE_SYMMETRY_FAMILY.md (the synthesis)");
            yield return new InspectableNode("operator-axis members wired here",
                summary: "F71, Z⊗N, X⊗N, the Π orbit of sectors");
            yield return new InspectableNode("parameter-axis members wired here",
                summary: "F91 (γ), F92 (J), F93 (h): the mirror-pair reshuffle with the site reflection, a Klein four-group on each axis");
            yield return new InspectableNode("members anchored to their own Claims",
                summary: "MirrorGroupD4Claim (F118), CommutatorDConjugationSign (F114), AntilinearTriangleClaim (F119), ChiralKClaim (K with the transpose), DeadSetLawClaim (F132), KleinFourGroupSelfPairedRefinement, FoldCheckerboardParityClaim (F150), BipartiteGaugeCriterionClaim (F151), MirrorOrderSortingClaim (F131), FrozenDivisorClaim (F140), PalindromeTwoEndCountClaim (F158), PalindromeComplementConnectionClaim (F165)");
            yield return new InspectableNode("substrate",
                summary: "U(1) × U(1) per-side popcount (continuous; JointPopcountSectors)");
        }
    }
}
