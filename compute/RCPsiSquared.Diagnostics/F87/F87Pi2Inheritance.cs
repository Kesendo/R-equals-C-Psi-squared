using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F87's trichotomy {truly, soft, hard} and its Π²-class 4-way refinement
/// {truly, pi2_odd_pure, pi2_even_nontruly, mixed} both inherit from the
/// Pi2-Foundation. Two structural inheritance edges, both Tier1Derived in their
/// own right:
///
/// <list type="bullet">
///   <item><b>F1 residual M as discriminator</b> (number-level): F87 classifies via
///         the F1 palindrome residual <c>M = Π·L·Π⁻¹ + L + 2σ·I</c>. The "2"
///         coefficient in F1's M is already typed via <see cref="F1Pi2Inheritance"/>
///         as <c>a_0</c> on the Pi2 dyadic ladder. F87 inherits this anchor
///         transitively.</item>
///   <item><b>4-way Π²-class refinement aligns with KleinFour</b> (structural-level):
///         F87's Π²-class refinement {truly / pi2_odd_pure / pi2_even_nontruly / mixed}
///         is structurally the same partition as <see cref="KleinFourCellClaim"/>'s
///         four cells (Pp / Pm / Mp / Mm = truly / Π²-even non-truly / Π²-odd subgroup
///         A / Π²-odd subgroup B; with "mixed" capturing F87 Hamiltonians spanning
///         multiple Klein cells). The 4-cell count is a Klein structural fact.</item>
/// </list>
///
/// <para>Tier1Derived: pure inheritance composition. F87 itself is Tier1Derived
/// (<see cref="F87TrichotomyClassification"/>); F1's M residual is Tier1Derived
/// in <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c>; the KleinFour decomposition is
/// Tier1Derived in F88. This claim makes the inheritance explicit through the
/// Object Manager.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F87 +
/// <c>experiments/V_EFFECT_FINE_STRUCTURE.md</c> +
/// <c>compute/RCPsiSquared.Diagnostics/F87/F87TrichotomyClassification.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (KleinFourCellClaim).</para></summary>
public sealed class F87Pi2Inheritance : Claim
{
    private readonly F1Pi2Inheritance _f1Inheritance;

    /// <summary>The number of Π²-class cells in the F87 4-way refinement, structurally
    /// equal to the KleinFour decomposition's 4 cells (Pp / Pm / Mp / Mm).</summary>
    public const int KleinClassCount = 4;

    /// <summary>The F87 4-way Π²-class names, in the order matching KleinFour cells
    /// (Pp / Pm / Mp / Mm).</summary>
    public IReadOnlyList<string> KleinAlignedClassNames => new[]
    {
        "truly",                  // Pp = (+, +)
        "pi2_even_nontruly",      // Pm = (+, −)
        "pi2_odd_subgroup_A",     // Mp = (−, +)
        "pi2_odd_subgroup_B",     // Mm = (−, −)
    };

    /// <summary>The "2" coefficient transitively inherited from F1's residual via
    /// <see cref="F1Pi2Inheritance.TwoFactor"/>. F87's discriminator <c>M</c> uses this
    /// "2" inside the F1 closed form; F87 inherits its number-level constant from
    /// the Pi2 dyadic ladder via F1.</summary>
    public double TransitivelyInheritedTwoFactor => _f1Inheritance.TwoFactor;

    public F87Pi2Inheritance(F1Pi2Inheritance f1Inheritance)
        : base("F87 trichotomy inherits from Pi2-Foundation: discriminator via F1 (2 = a_0); 4-way Π²-classes align with KleinFour",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F87 + " +
               "experiments/V_EFFECT_FINE_STRUCTURE.md + " +
               "compute/RCPsiSquared.Diagnostics/F87/F87TrichotomyClassification.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    {
        _f1Inheritance = f1Inheritance ?? throw new ArgumentNullException(nameof(f1Inheritance));
    }

    public override string DisplayName =>
        "F87 ←(F1 residual + KleinFour 4-cell)← Pi2-Foundation";

    public override string Summary =>
        $"F87 trichotomy {{truly, soft, hard}} discriminated via F1 residual M (2 = a_0 = {TransitivelyInheritedTwoFactor}); " +
        $"4-way Π²-class refinement = KleinFour {KleinClassCount} cells ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F87 closed form",
                summary: "trichotomy via F1 residual: truly iff ‖M‖ < ε; soft iff M ≠ 0 but spectrum pairs; hard iff no pairing");
            yield return new InspectableNode("number-level inheritance via F1",
                summary: $"F87's discriminator M uses F1's 2σ·I shift; the 2 = a_0 = {TransitivelyInheritedTwoFactor} comes from Pi2 dyadic ladder");
            yield return new InspectableNode("structural inheritance via KleinFour",
                summary: $"F87's 4-way Π²-class refinement aligns with KleinFour cells; {KleinClassCount} classes total");
            yield return new InspectableNode("KleinAlignedClassNames",
                summary: string.Join(" / ", KleinAlignedClassNames) + " (matching Pp / Pm / Mp / Mm)");
            yield return InspectableNode.RealScalar("TransitivelyInheritedTwoFactor (= a_0 via F1)", TransitivelyInheritedTwoFactor);
            yield return new InspectableNode("F87 hardware confirmation",
                summary: "Marrakesh d7mjnjjaq2pc73a1pk4g (2026-04-26) Δ(soft − truly) = −0.722; Kingston regime-uniformity (2026-05-05)");
            yield return new InspectableNode("V-Effect counts",
                summary: "36-enum N=3 → 14/19/3 (truly/soft/hard); 120-enum N=3..5 → 15/46/59 N-stable (combinatorial proof commit 81caf67)");
        }
    }
}
