using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.SymmetryFamily;

/// <summary>ZGlobalMirrorRefinement (Tier 1 derived; 2026-05-12): the global Z-string
/// operator Z⊗N = ⊗_l Z_l acts on a Liouville-space basis pair |a⟩⟨b| as the
/// scalar (−1)^{popcount(a) + popcount(b)} = (−1)^{p_c + p_r}. This is exactly the
/// parity of the joint-popcount label, so Z⊗N is trivially redundant with the
/// joint-popcount-parity classification. Typed for inventory completeness only; does
/// not produce a new sector splitting beyond <see cref="BlockSpectrum.JointPopcountSectors"/>.
/// See <c>compute/RCPsiSquared.Core/Symmetry/ZGlobalMirror.cs</c> for the raw operator.</summary>
public sealed class ZGlobalMirrorRefinement : Claim
{
    private readonly SymmetryFamilyInventory _inventory;

    public ZGlobalMirrorRefinement(SymmetryFamilyInventory inventory)
        : base("ZGlobalMirrorRefinement: Z⊗N eigenvalue (-1)^{p_c+p_r} on |a⟩⟨b| equals joint-popcount parity; trivially redundant with JointPopcountSectors classification.",
               Tier.Tier1Derived,
               "Z⊗N |a⟩⟨b| = (-1)^{popcount(a)+popcount(b)} |a⟩⟨b| algebraic identity")
    {
        _inventory = inventory ?? throw new ArgumentNullException(nameof(inventory));
    }

    /// <summary>Joint-popcount parity (p_c + p_r) mod 2 = Z⊗N eigenvalue index for the
    /// Liouville-space basis pair labeled by (p_c, p_r). Zero corresponds to +1 eigenvalue,
    /// one corresponds to −1.</summary>
    public static int JointPopcountParity(int pCol, int pRow) => (pCol + pRow) % 2;

    public override string DisplayName =>
        "ZGlobalMirrorRefinement: Z⊗N is redundant with joint-popcount parity";

    public override string Summary =>
        $"Z⊗N classification = joint-popcount parity (p_c+p_r) mod 2; redundant, no new sector split ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parity formula",
                summary: "Z⊗N |a⟩⟨b| eigenvalue = (-1)^{popcount(a)+popcount(b)} = (-1)^{p_c+p_r}");
            yield return new InspectableNode("inventory role",
                summary: "typed for completeness; does NOT refine joint-popcount blocks");
            yield return new InspectableNode("raw operator",
                summary: "compute/RCPsiSquared.Core/Symmetry/ZGlobalMirror.cs");
        }
    }
}
