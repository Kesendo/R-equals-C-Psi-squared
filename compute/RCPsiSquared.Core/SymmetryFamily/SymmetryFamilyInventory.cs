using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>SymmetryFamilyInventory (Tier 1 derived; 2026-05-12): aggregator-Claim
/// listing all discrete symmetries of the chain XY+Z-dephasing Liouvillian, classified
/// by axis (operator vs parameter) and effect (sector-splitting vs sector-pairing).
/// See <c>docs/SYMMETRY_FAMILY_INVENTORY.md</c> for the full table.</summary>
public sealed class SymmetryFamilyInventory : Claim
{
    public SymmetryFamilyInventory()
        : base("SymmetryFamilyInventory: typed family of discrete symmetries (5 axes + 2 sector-pairing primitives + 1 negative result) of chain XY+Z-deph L, classified by axis (operator/parameter) and effect (splitting/pairing).",
               Tier.Tier1Derived,
               "docs/SYMMETRY_FAMILY_INVENTORY.md + reflections/ON_THE_SYMMETRY_FAMILY.md")
    { }

    public override string DisplayName =>
        "SymmetryFamilyInventory: typed family of discrete L-symmetries on operator + parameter axes";

    public override string Summary =>
        $"family inventory aggregating F71 + F91 + F92 + F93 + Z⊗N + X⊗N + F71_col×F71_row negative result; see docs/SYMMETRY_FAMILY_INVENTORY.md ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("anchor markdown",
                summary: "docs/SYMMETRY_FAMILY_INVENTORY.md (the table)");
            yield return new InspectableNode("anchor reflection",
                summary: "reflections/ON_THE_SYMMETRY_FAMILY.md (the synthesis)");
            yield return new InspectableNode("operator-axis members",
                summary: "F71, Z⊗N, X⊗N");
            yield return new InspectableNode("parameter-axis members (F-chain twins)",
                summary: "F91 (γ), F92 (J), F93 (h)");
            yield return new InspectableNode("substrate",
                summary: "U(1) × U(1) per-side popcount (continuous; JointPopcountSectors)");
        }
    }
}
